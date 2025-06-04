"""Barchart organized in subplots (facets)."""

import base64
import pathlib
from dataclasses import dataclass, field
from math import copysign
from pathlib import Path

import folium
import geopandas as gpd
import pandas as pd
from folium import GeoJson, plugins

from evals.constants import ALIAS_COUNTRY, ALIAS_REGION, DataModel
from evals.data.icons import RIGHT_TO_BRACKET_SOLID
from evals.utils import filter_by, prettify_number


@dataclass
class GridMapConfig:
    """Transmission grip map configuration."""

    # This layer will be visible by default
    show_year: str = "2030"

    crs: int = 4326

    map_center: list = field(default_factory=lambda: [41.9, 15])

    zoom_start: int = 5
    zoom_min: int = 4
    zoom_max: int = 8

    bounds: dict = field(default_factory=lambda: {"N": 65, "E": 30, "S": 15, "W": -5})
    max_bounds: bool = True  # cannot move map away from bounds

    tile_provider: str = (
        "https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png"
    )

    # required for copyright and licensing reasons
    attribution: str = (
        '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> '
        'contributors, &copy; <a href="https://cartodb.com/attributions">CartoDB</a>'
    )

    # grid lines
    line_capacity_threshold: float = 0.1  # GWh
    line_weight_min: float = 2.0  # px
    line_weight_max: float = 20.0  # px
    # ToDo: Add colors
    carrier_style: dict = field(
        default_factory=lambda: {
            "AC": {  # AC Lines
                "color": "#D90429",
                "nice_name": "AC",
                "offset": -10,
            },
            "DC": {
                "color": "#E19990",
                "nice_name": "DC",
                "offset": 10,  # px
            },
            "gas pipeline": {
                "color": "#63A452",
                "nice_name": "Methane (brownfield)",
            },
            "gas pipeline new": {
                "color": "#8BA352",
                "nice_name": "Methane (new)",
            },
            "H2 pipeline": {
                "color": "#258994",
                "nice_name": "H2",
                "offset": -10,
            },
            "H2 pipeline retrofitted": {
                "color": "#255194",
                "nice_name": "H2 (retrofitted)",
                "dash_array": "10",  # px, equal gaps
                "offset": 10,  # px
            },
            "H2 pipeline (Kernnetz)": {
                "color": "#259468",
                "nice_name": "H2 (Kernnetz)",
                "dash_array": "10",  # px, equal gaps
                "offset": 10,  # px
            },
        }
    )


class TransmissionGridMap:
    """
    Creates a map with transmission capacities.

    Parameters
    ----------
    grid
        A data frame with the transmission grid capacities.
    import_energy
        A data frame with global import energy import amounts and
        capacity.
    import_capacity
        A data frame with global import capacities.
    buses
        The buses data frame, used to determine the bus coordinates.
    config
        The GridMapConfig object with configuration options.
    """

    def __init__(
        self,
        grid: pd.DataFrame,
        import_energy: pd.DataFrame,
        import_capacity: pd.DataFrame,
        buses: pd.DataFrame,
        config: GridMapConfig,
    ) -> None:
        self.grid = grid
        self.import_energy = import_energy
        self.import_capacity = import_capacity
        self.buses = buses  # todo: do not uses buses from one year
        self.cfg = config
        self.fmap = folium.Map(
            tiles=None,
            location=self.cfg.map_center,
            zoom_start=self.cfg.zoom_start,
            max_bounds=self.cfg.max_bounds,
            max_lat=self.cfg.bounds["N"],
            max_lon=self.cfg.bounds["E"],
            min_lat=self.cfg.bounds["S"],
            min_lon=self.cfg.bounds["W"],
        )

        # feature groups are layers in the map and can be shown or hid
        self.feature_groups = {}
        for year in sorted(grid.index.unique(DataModel.YEAR)):
            fg = folium.FeatureGroup(name=year, show=True)
            self.feature_groups[year] = fg
            self.fmap.add_child(fg)  # register the feature group

    def save(
        self, output_path: pathlib.Path, file_name: str, subdir: str = "HTML"
    ) -> None:
        """
        Write the map to a html file.

        We want to store the HTML inside the JSON folder by default,
        because Folium does not support the import of JSON files.
        Therefore, we dump HTML files and include them as iFrames in
        the web UI instead of importing JSONs via the plotly library.

        Parameters
        ----------
        output_path
            The path to save the map in.
        file_name
            The name of the file to export the map to.
        subdir
            An optional subdirectory to store files at. Leave emtpy
            to skip, or change to html.
        """
        output_path = self.make_evaluation_result_directories(output_path, subdir)
        self.fmap.save(output_path / f"{file_name}.html")

    def make_evaluation_result_directories(
        self, result_path: Path, subdir: Path | str
    ) -> Path:
        """
        Create all directories needed to store evaluations results.

        Parameters
        ----------
        result_path
            The path of the result folder.
        subdir
            A relative path inside the result folder.

        Returns
        -------
        :
            The joined path: result_dir / subdir.
        """
        output_path = self.make_directory(result_path, subdir)
        output_path = self.make_directory(output_path, "HTML")

        return output_path

    def draw_grid_by_carrier_groups_myopic(self) -> None:
        """Plot carrier groups for all years to one map."""
        self.add_basemap_layers()

        plot_grid = self._calculate_line_weights(self.grid)

        _groups = [DataModel.YEAR, "bus0", "bus1"]
        year_edge = plot_grid.groupby(_groups, group_keys=False)

        plot_grid = year_edge.apply(self._calculate_line_offset)

        grid_line = plot_grid.groupby(plot_grid.index.names, group_keys=False)
        plot_grid = grid_line.apply(self._calculate_line_center)

        plot_grid.groupby([DataModel.YEAR, DataModel.CARRIER]).apply(
            self._draw_grid_polyline_with_circle_marker
        )

        self.draw_country_markers()
        # self.draw_import_locations()
        self.add_control_widgets()

    def add_control_widgets(self) -> None:
        """Add UI elements to the map."""
        plugins.GroupedLayerControl(
            groups={"Year": list(self.feature_groups.values())},
            collapsed=False,
            position="topleft",
        ).add_to(self.fmap)

        plugins.Fullscreen(
            position="topright",
            title="Full Screen",
            title_cancel="Exit Full Screen",
            force_separate_button=True,
        ).add_to(self.fmap)

    def draw_import_locations(self) -> None:
        """
        Add import location icons and lines to the map.

        Notes
        -----
        Available icons: https://fontawesome.com/icons/categories
        """
        icon_locations = {
            # node: [y, x]  Lat, Lon
            "BE0 0": [51.21868, 2.86993],
            "DE0 4": [53.92445, 8.67684],
            "EE6 0": [58.78505, 23.15726],
            "ES0 0": [43.41430, -4.27864],
            "FI3 0": [60.89107, 22.68793],
            "FR0 0": [48.18790, -3.68987],
            "GB5 0": [54.63442, -0.70133],
            "GR0 0": [38.67136, 26.65004],
            "HU0 0": [48.20228, 22.60233],
            "IT0 0": [37.16396, 13.49807],
            "LT6 0": [55.71138, 21.07711],
            "LV6 0": [56.99505, 27.72035],
            "NL0 0": [53.03334, 4.96787],
            "NO3 0": [60.05473, 5.00377],
            "PL0 0": [51.99998, 22.13991],
            "PT0 0": [37.98446, -8.88731],
            "RO0 0": [44.51848, 28.89059],
            "SK0 0": [48.78314, 22.35254],
        }

        _idx = self.import_capacity.index.names
        row_slices = self.import_capacity.to_frame().groupby(_idx, group_keys=False)
        import_capacity = row_slices.apply(self._calculate_line_weights)

        for (year, node), capas in import_capacity.groupby(
            [DataModel.YEAR, DataModel.LOCATION]
        ):
            fg = self.feature_groups[year]
            # need a new instance for every icon
            icon = self._get_icon(RIGHT_TO_BRACKET_SOLID)
            popup_table = filter_by(self.import_energy, year=year, location=node)
            popup_table = popup_table.droplevel(
                [DataModel.YEAR, DataModel.CARRIER, DataModel.BUS_CARRIER]
            )
            bootstrap5_classes = (
                "table table-striped table-hover table-condensed table-responsive"
            )
            popup_html = popup_table.to_frame().to_html(classes=bootstrap5_classes)
            folium.Marker(
                location=icon_locations[node],
                icon=folium.CustomIcon(icon, icon_size=(10, 10)),
                popup=folium.Popup(popup_html),
                tooltip="Global Import",
            ).add_to(fg)

            # draw line from import icon location to node location
            node_y = self.buses.loc[node, "y"]
            node_x = self.buses.loc[node, "x"]
            capacity = capas[self.import_capacity.name].iloc[0]
            label = f"{capas.attrs['name']}: {capacity:.2f} {capas.attrs['unit']}"
            folium.PolyLine(
                locations=[icon_locations[node], [node_y, node_x]],
                color="black",
                weight=capas["line_weight"].iloc[0],
                tooltip=label,
                popup=label,
            ).add_to(fg)

    def add_basemap_layers(self) -> None:
        """Add common background layer to the map."""
        self._add_wms_tiles()
        self._load_geojson(
            "regions_onshore_base_s_adm.geojson",
            style={
                "weight": 1,
                "color": "grey",
                "fillColor": "white",
                "opacity": 0.5,
            },
        )

        # self._load_geojson(
        #     "neighbors.geojson",
        #     style={
        #         "weight": 0.5,
        #         "color": "black",
        #         "fillColor": "black",
        #         "opacity": 0.2,
        #     },
        # )

    def draw_country_markers(self) -> None:
        """
        Draw markers for countries on the map.

        Retrieves bus information from networks, iterates over unique
        bus locations, creates CircleMarker and Marker objects for
        each location with corresponding short and nice names, and
        adds them to the map.
        """
        fg_labels = folium.FeatureGroup(
            name="Country Marker", overlays=True, interactive=False
        )

        buses0 = self.grid.index.unique("bus0")
        buses1 = self.grid.index.unique("bus1")
        icon_css = "margin-top:1.5px; font-size:10px; font-family:sans-serif"

        for bus in buses0.union(buses1):
            # keep region ID for AT and DE, else just the country code
            short_name = bus[:2] + bus[-1] if bus.startswith(("AT", "DE")) else bus[:2]
            nice_name = ALIAS_REGION.get(bus, ALIAS_COUNTRY[bus[:2]])
            location = self.buses.loc[bus, ["y", "x"]].to_numpy()

            icon = plugins.BeautifyIcon(
                icon_shape="circle",
                border_width=2,
                border_color="black",
                # background_color="white",
                text_color="black",
                inner_icon_style=icon_css,
                number=short_name,
            )
            marker = folium.Marker(location=location, popup=nice_name, icon=icon)
            marker.add_to(fg_labels)

        fg_labels.add_to(self.fmap)

    @staticmethod
    def _get_icon(icon: str) -> str:
        """
        Encode a raw SVG string to bytes.

        Parameters
        ----------
        icon
            The utf-8 encoded HTML representation of an SVG icon.

        Returns
        -------
        :
            The base64 encoded SVG icon as a string.
        """
        data = base64.b64encode(icon.strip().encode("utf-8")).decode("utf-8")
        return f"data:image/svg+xml;base64,{data}"

    def _calculate_line_weights(self, df_slice: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate the line weights for a grid.

        Parameters
        ----------
        df_slice
            The grids that will be plotted to the map.

        Returns
        -------
        :
            The grids with an additional column for the line weight in px.
        """
        # prevent assignment to copies of a data view
        df_slice = df_slice.copy()

        col = f"{df_slice.attrs['name']} ({df_slice.attrs['unit']})"

        _min, _max = df_slice[col].min(), df_slice[col].max()
        _max_width = self.cfg.line_weight_max  # 20.0
        _min_width = self.cfg.line_weight_min  # 2.0

        def linear_scale(ser: pd.Series) -> pd.Series:
            """
            Scale values between lower and upper line weight values.

            Returns the linear equation k * x + d. Where x is the ratio,
            k is the min-max range and d the lower bound constant.

            Parameters
            ----------
            ser
                The values to be scaled between 2 and 5.

            Returns
            -------
            The scaled value used as line width in pixel.
            """
            min_max_ratio = (ser - _min) / (_max - _min)
            return min_max_ratio * (_max_width - _min_width) + _min_width

        if _min == _max:
            df_slice.loc[:, "line_weight"] = self.cfg.line_weight_min
        else:
            df_slice.loc[:, "line_weight"] = df_slice[col].apply(linear_scale)

        return df_slice

    @staticmethod
    def _calculate_line_offset(df_slice: pd.DataFrame, gap: float = 1) -> pd.DataFrame:
        """
        Add a column with the offset values for a grid.

        The offset is used to prevent overplotting lines and labels.
        It is only required if multiple edges exist between the same
        nodes.

        Parameters
        ----------
        df_slice
            A data frame slice for every unique node connection
            and for all displayed carriers in a map.
        gap : float, optional
            Number of pixels to insert between the edges of adjacent lines
            (default is 1 px), preventing strokes from touching.

        Returns
        -------
        :
            The input data slice with the offset in pixel in an
            additional column.
        """
        # if df_slice.shape[0] == 1:
        #     df_slice["offset"] = 0
        # elif df_slice.shape[0] == 2:
        #     # move lines (down and up) by half their combined line
        #     # weights plus 1 px for a visible gap
        #     half_weight = df_slice["line_weight"].sum() / 2 + 1
        #     df_slice["offset"] = [-0.5 * half_weight, 0.5 * half_weight]
        # else:
        #     raise NotImplementedError(f"Number of rows: {df_slice} not supported.")
        #
        # return df_slice
        weights = df_slice["line_weight"].astype(float).tolist()
        n = len(weights)

        # Total envelope width = sum of all stroke widths + (n-1) gaps
        total_width = sum(weights) + (n - 1) * gap

        # Start at left edge of that envelope
        current = -total_width / 2
        offsets = []

        # For each stroke:
        #   • move by half its width  → centerline of this band
        #   • record that as the offset
        #   • then advance by (half its width + gap) to get to the next band’s start
        for w in weights:
            current += w / 2
            offsets.append(current)
            current += w / 2 + gap

        # Attach offsets and return
        df_slice["offset"] = offsets
        return df_slice

    @staticmethod
    def _calculate_line_center(df_slice: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate the line center for all lines.

        In case the line has an offset, the center is moved by 10%
        of the line length and along the line in the direction of the
        offset.

        Parameters
        ----------
        df_slice
            The data frame with a "line" column that contains
            coordinate pairs and an "offset" column that contains
            a positive or negative float for the line offset.

        Returns
        -------
        :
            The input data frame with additional column containing the
            line center.
        """

        def compute_center(row):
            offset = row["offset"]
            line = row["line"]
            if offset != 0:
                x0, x1 = line[0][0], line[1][0]
                y0, y1 = line[0][1], line[1][1]
                # Move center by +-10% of line length depending on offset sign.
                ratio = 0.5 + copysign(0.1, offset)
                x = x0 + ratio * (x1 - x0)
                y = y0 + ratio * (y1 - y0)
                return [x, y]
            else:
                # Compute the simple midpoint
                return [(line[0][i] + line[1][i]) / 2 for i in range(len(line[0]))]

        df_slice["line_center"] = df_slice.apply(compute_center, axis=1)

        return df_slice

    def _draw_grid_polyline_with_circle_marker(self, grid: pd.DataFrame) -> None:
        """
        Draw grid lines on the map for a specific carrier.

        Retrieves the nice name, color, and dash array configuration
        for the carrier. Iterates over grid data, creates PolyLine
        objects for capacities larger than the threshold in the config,
        sets color, weight, dash array, tooltip, and popup based on the
        capacity, and adds them to the FeatureGroup.

        Parameters
        ----------
        grid
            The metric dataframe with the capacities and the lines.
        """
        year = grid.index.unique(DataModel.YEAR)[0]
        carrier = grid.index.unique(DataModel.CARRIER)[0]

        fg = self.feature_groups[year]
        style = self.cfg.carrier_style[carrier]
        nice_name = style["nice_name"]
        color = style["color"]
        unit = grid.attrs["unit"]

        col = f"{grid.attrs['name']} ({grid.attrs['unit']})"
        significant_edges = grid[grid[col] >= self.cfg.line_capacity_threshold]

        for _, row in significant_edges.iterrows():
            capacity = row[col]  # / 1000  # GW
            tooltip = f"{nice_name}: {capacity:.2f} {unit}"

            plugins.PolyLineOffset(
                locations=row["line"],
                offset=row["offset"],
                color=color,
                weight=row["line_weight"],
                dash_array=f"{row['line_weight']}" if style.get("dash_array") else None,
                # https://www.w3schools.com/graphics/svg_stroking.asp
                line_cap="butt",  # or "round"
                tooltip=tooltip,
                popup=tooltip,
            ).add_to(fg)

            # https://github.com/masajid390/BeautifyMarker
            icon_css = "margin-top:2.5px; font-size:10px; font-family:sans-serif"
            icon = plugins.BeautifyIcon(
                icon_shape="circle",
                border_width=1,
                border_color=color,
                background_color="white",
                text_color=color,
                number=prettify_number(capacity),
                inner_icon_style=icon_css,
            )

            folium.Marker(
                location=row["line_center"],
                popup=folium.Popup(nice_name),
                icon=icon,
                tooltip=tooltip,
            ).add_to(fg)

    def _add_wms_tiles(self) -> None:
        """Add a web map tile service layer to the map."""
        folium.TileLayer(
            name="WMS tiles",
            control=False,  # no layer controls
            tiles=self.cfg.tile_provider,
            attr=self.cfg.attribution,
            min_zoom=self.cfg.zoom_min,
            max_zoom=self.cfg.zoom_max,
            # additional leaflet.js args: https://leafletjs.com/reference.html#tilelayer
        ).add_to(self.fmap)

    def _load_geojson(self, file_name: str, style: dict = None) -> None:
        """
        Add the geojson layer.

        Parameters
        ----------
        file_name
            The name of the geojson file under esmtools/data.
        style
            The style dictionary to pass to the geojson layer.
        """
        # res = resources.files("evals") / "data"
        # gdf = gpd.read_file(res / file_name).to_crs(crs=f"EPSG:{self.cfg.crs}")
        gdf = gpd.read_file(Path("resources") / file_name).to_crs(
            crs=f"EPSG:{self.cfg.crs}"
        )
        if style:  # applies the same style to all features
            gdf["style"] = [style] * gdf.shape[0]

        gj = GeoJson(gdf, control=False, overlay=True)
        gj.add_to(self.fmap)

    @staticmethod
    def make_directory(base: Path, subdir: Path | str) -> Path:
        """
        Create a directory and return its path.

        Parameters
        ----------
        base
            The path to base of the new folder.
        subdir
            A relative path inside the base folder.

        Returns
        -------
        :
            The joined path: result_dir / subdir / now.
        """
        base = Path(base).resolve()
        assert base.is_dir(), f"Base path does not exist: {base}."
        directory_path = base / subdir
        directory_path.mkdir(parents=True, exist_ok=True)

        return directory_path
