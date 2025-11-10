# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Creates European CO2 pipeline network from project collection KML file.
https://www.google.com/maps/d/u/0/viewer?mid=1prz_ns6tdj_1kacbrcm47q_299-3QxA

"""

import logging
from itertools import chain

import fiona
import geopandas as gpd
import pandas as pd
import pypsa
import shapely
from pypsa.geo import haversine_pts
from shapely import segmentize, unary_union
from shapely.algorithms.polylabel import polylabel
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point
from shapely.ops import nearest_points, linemerge

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)


logger = logging.getLogger(__name__)

CLUSTER_TOL = 25000 # in meters
DISTANCE_CRS = "EPSG:3035"
GEO_CRS = "EPSG:4326"
OFFSHORE_BUS_RADIUS = 10000 # in meters
COASTAL_DISTANCE = 50000 # in meters
PIPELINE_LABEL = "Infrastruktur"
REGIONS_ONSHORE_BUFFER = 30000
PIPELINE_LABELS = {
    "Aramis Projekt - PCI Projekt": "Aramis",
    "Arcon": "Acorn",
    "Acorn": "Acorn",
    "CO2TransPorts - PCI Liste": "CO2TransPorts",
    "CO2TransPorts  - PCI Liste": "CO2TransPorts",
    "CarbonConnect - PCI Liste": "CarbonConnect",
    "Delta Rhine Corridor - Pipeline - PCI Liste": "Delta Rhine Corridor",
    "EU2NSEA - PCI-Projekt": "EU2NSEA",
    "N-Lites - PMI liste": "N-Lites",
    "Nordsee CO₂-Korridor Deutschland-Belgien": "Belgium backbone",
    "OGE - Cluster Elbmündung": "OGE",
    "OGE - North Sea CO₂ Corridor - Cluster rheinisches Revier": "OGE",
    "OGE - WHV CO2 Corridor ": "OGE",
    "OGE - WHV CO2 Corridor": "OGE",
    "OGE CO₂-Transportnetz - Pipelinenetz - PCI Liste": "OGE",
    "PYCASSO - PCI Liste": "PYCASSO",
}
PIPELINE_COLS = ["Name", "description", "geometry"]
STORE_LABEL = "Speicherstätten"
STORE_LABELS = {
    "Acorn - East May Storage Site": "East Mey",
    "Acorn - South Storage Site": "Acorn",
    "Aramis Projekt": "Aramis",
    "CarbonConnect - PCI Liste": "CarbonConnect",
    "Erdgasfeld P18-A": "CO2TransPorts",
    "EU2NSEA - Luna storage site - Wintershall Dea": "EU2NSEA",
    "EU2NSEA - Smeaheia storage site - Equinor": "EU2NSEA",
    "N-LiTES - Aquifer Aurora - PMI Liste": "N-Lites",
}
STORE_COLS = ["Name", "description", "geometry"]
MAX_STORE_DISTANCE = 10000 # in meters


def clean_text(s):
    if isinstance(s, str):
        return s.replace("\xa0", " ").replace("\n", "").strip()
    return s


def calculate_haversine_distance(buses, lines, line_length_factor):
    coords = buses[["x", "y"]]

    lines.loc[:, "length"] = (
        haversine_pts(coords.loc[lines["bus0"]], coords.loc[lines["bus1"]])
        * line_length_factor
    ).round(1)

    return lines


def create_new_buses(
    gdf: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    scope: gpd.GeoDataFrame,
    carrier: str,
    regions_onshore_buffer: int = REGIONS_ONSHORE_BUFFER,
    tol: int = CLUSTER_TOL,
    offset: int = 0,
):
    source_crs = gdf.crs
    
    buffered_regions = (
        regions_onshore
        .buffer(regions_onshore_buffer)
        .union_all()  # Coastal buffer
    )

    # filter all rows in gdf where at least one of the geometry linestring endings is outside unary_union(regions_onshore)
    # create a list of Points of all linestring endings in gdf
    list_points = list(
        chain(*gdf.geometry.apply(lambda x: [Point(x.coords[0]), Point(x.coords[-1])]))
    )
    # create multipoint geometry of all points
    list_points = MultiPoint(list_points)
    list_points = list_points.intersection(scope.union_all())
    # Drop all points that are within unary_union(regions_onshore) and a buffer of 5000 meters
    list_points = list_points.difference(buffered_regions)

    gdf_points = gpd.GeoDataFrame(
        geometry=[geom for geom in list_points.geoms],
        crs = source_crs,
    )

    gdf_points["geometry"] = gdf_points.buffer(tol)

    # Aggregate rows with touching polygons
    gdf_points = gdf_points.dissolve()
    # split into separate polygons
    gdf_points = gdf_points.explode().reset_index(drop=True)

    gdf_points["poi"] = (
        gdf_points["geometry"]
        .apply(lambda polygon: polylabel(polygon, tolerance=tol / 2))
    )

    # Extract x and y coordinates into separate columns
    gdf_points_geo = gdf_points[["poi"]].copy()
    gdf_points_geo = gpd.GeoDataFrame(gdf_points_geo, geometry="poi", crs=gdf_points.crs)
    gdf_points_geo.to_crs(GEO_CRS, inplace=True)
    gdf_points["x"] = gdf_points_geo["poi"].x
    gdf_points["y"] = gdf_points_geo["poi"].y
    gdf_points["name"] = gdf_points.apply(
        lambda x: f"OFFSHORE {int(x.name)+1+offset}", axis=1
    )
    gdf_points.set_index("name", inplace=True)
    gdf_points["carrier"] = carrier

    return gdf_points[["x", "y", "carrier", "geometry"]]


def drop_z_dim(geom):
    if geom is None or geom.is_empty:
        return geom
    return shapely.force_2d(geom)


def explode_linestrings_to_segments(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    seg_records = []
    for group_id, row in enumerate(gdf.itertuples(index=False), start=1):
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "LineString":
            lines = [geom]
        elif geom.geom_type == "MultiLineString":
            lines = geom.geoms
        else:
            continue
        for line in lines:
            coords = list(line.coords)
            for a, b in zip(coords[:-1], coords[1:]):
                rec = {**row._asdict(), "geometry": LineString([a, b]), "group_id": group_id}
                seg_records.append(rec)
    return gpd.GeoDataFrame(seg_records, crs=gdf.crs)


def split_multilinestring(row):
    """
    Splits rows containing a MultiLineString geometry into multiple rows,
    converting them to a single LineString. New rows inherit all other
    attributes from the original row. Non-MultiLineString rows are returned as-
    is.

    Parameters:
        row (pd.Series): A pandas Series containing a 'geometry' column with a MultiLineString or LineString.

    Returns:
        row (pd.Series): A row containing a LineString geometry including their original attributes.
    """
    geom = row["geometry"]
    if isinstance(geom, MultiLineString):
        # Convert MultiLineString into a list of LineStrings
        lines = [line for line in geom.geoms]
        # Create a DataFrame with the new rows, including all other columns
        return pd.DataFrame(
            {
                "geometry": lines,
                **{
                    col: [row[col]] * len(lines)
                    for col in row.index
                    if col != "geometry"
                },
            }
        )
    else:
        # Return the original row as a DataFrame, including all columns
        return pd.DataFrame([row])


def find_points_on_line_overpassing_region(
    link, regions,
):

    overlap = gpd.overlay(link, regions)

    # All rows with multilinestrings, split them into their individual linestrings and fill the rows with the same data
    overlap = pd.concat(
        overlap.apply(split_multilinestring, axis=1).tolist(), ignore_index=True
    )

    overlap["center_point"] = overlap["geometry"].apply(
        lambda l: l.interpolate(l.length / 2)
    )

    overlap["on_point"] = overlap.apply(
        lambda row: nearest_points(row["center_point"], row["geometry"])[1], axis=1
    )

    return overlap[["on_point"]].rename(columns={"on_point": "geometry"})


def count_intersections(line, polygons):
    return sum(line.intersects(polygon) for polygon in polygons)


def split_to_overpassing_segments(
    gdf: gpd.GeoDataFrame,
    regions: gpd.GeoDataFrame,
    distance_crs: str = DISTANCE_CRS,
):
    logger.info("Splitting linestrings into segments that connect overpassing regions.")
    buffer_radius = 1 # m

    ## Delete later
    gdf_split = gdf.copy().to_crs(distance_crs)
    regions_dist = regions.to_crs(distance_crs)

    # Increase resolution of both geometries
    gdf_split["geometry"] = gdf_split["geometry"].apply(lambda x: segmentize(x, 200))
    regions_dist["geometry"] = regions_dist["geometry"].apply(
        lambda x: segmentize(x, 300)
    )

    # Do the following splitting operation only for lines that overpass multiple regions
    crosses_multiple = gdf_split.geometry.apply(lambda line: count_intersections(line, regions_dist.geometry)) > 2

    if crosses_multiple.any():
        gdf_points = find_points_on_line_overpassing_region(gdf_split.loc[crosses_multiple], regions_dist)
        gdf_points = gpd.GeoDataFrame(gdf_points, crs=distance_crs)

        gdf_points["buffer"] = gdf_points["geometry"].buffer(buffer_radius)

        # Split linestrings of gdf by union of points[buffer]
        gdf_split["geometry"] = gdf_split["geometry"].apply(
            lambda x: x.difference(gdf_points["buffer"].union_all())
        )

    # Drop empty geometries
    gdf_split = gdf_split[~gdf_split["geometry"].is_empty]

    gdf_split.reset_index(inplace=True)
    # All rows with multilinestrings, split them into their individual linestrings and fill the rows with the same data
    gdf_split = pd.concat(
        gdf_split.apply(split_multilinestring, axis=1).tolist(), ignore_index=True
    )

    gdf_split = gpd.GeoDataFrame(gdf_split, geometry="geometry", crs=distance_crs)

    # Drop empty geometries
    gdf_split = gdf_split[~gdf_split["geometry"].is_empty]

    # Recalculate lengths
    gdf_split["length"] = (
        gdf_split["geometry"].length.div(1e3).round(1)
    )  # Calculate in km, round to 1 decimal

    return gdf_split.reset_index(drop=True)


def map_endpoints_to_closest_region(
    gdf,
    regions,
    max_distance=OFFSHORE_BUS_RADIUS,
    coords=0,
    lines=True,
):
    """
    Maps endpoints in a GeoDataFrame to their closest regions within a specified maximum distance.
    
    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame containing geometries (line geometries if lines=True, point geometries otherwise).
    regions : GeoDataFrame
        GeoDataFrame containing region geometries with a 'name' column.
    max_distance : float, optional
        Maximum allowed distance between points and regions. Points farther than this 
        will have their region set to None. Default is OFFSHORE_BUS_RADIUS.
    coords : int, optional
        Index of the coordinate to extract from line geometries when lines=True. Default is 0.
    lines : bool, optional
        Whether gdf contains line geometries. If True, points are extracted from line
        geometries using coords. If False, gdf geometries are treated as points. Default is True.
    
    Returns
    -------
    pandas.Series
        Series containing the name of the closest region for each endpoint, or None if
        the closest region is farther than max_distance.
    """
    if lines:
        gdf_points = gdf.geometry.apply(lambda x: Point(x.coords[coords]))
    else:
        gdf_points = gdf.geometry

    gdf_points = gpd.GeoDataFrame(geometry=gdf_points)
    # Spatial join nearest with regions

    # Find nearest region index
    regions = regions.to_crs(DISTANCE_CRS)
    gdf_points = gdf_points.to_crs(DISTANCE_CRS)

    gdf_points = gpd.sjoin_nearest(gdf_points, regions, how="left")
    gdf_points = gdf_points.join(
        regions, on="name", lsuffix="_point", rsuffix="_region"
    )
    gdf_points["distance"] = gdf_points.apply(
        lambda x: x.geometry_point.distance(x.geometry_region), axis=1
    )

    bool_too_far = gdf_points["distance"] > max_distance
    gdf_points.loc[bool_too_far, "name"] = None

    return gdf_points["name"]


def map_to_closest_region(
    gdf, regions, max_distance=OFFSHORE_BUS_RADIUS, add_suffix=None
):
    # add Suffix to regions index
    regions = regions.copy()
    if add_suffix:
        regions.index = regions.index + " " + add_suffix

    gdf = gdf.copy()
    # if columns bus0 and bus1 dont exist, create them
    if "bus0" not in gdf.columns:
        gdf["bus0"] = None
    if "bus1" not in gdf.columns:
        gdf["bus1"] = None

    # Apply mapping to rows where 'bus0' is None
    gdf.loc[gdf["bus0"].isna(), "bus0"] = map_endpoints_to_closest_region(
        gdf[gdf["bus0"].isna()], regions, max_distance, coords=0
    )

    # Apply mapping to rows where 'bus1' is None
    gdf.loc[gdf["bus1"].isna(), "bus1"] = map_endpoints_to_closest_region(
        gdf[gdf["bus1"].isna()], regions, max_distance, coords=-1
    )

    return gdf


def safe_linemerge(geoms):
    merged = unary_union(geoms)
    # linemerge requires a MultiLineString or list of lines
    if isinstance(merged, shapely.LineString):
        return merged
    return linemerge(merged)


def set_underwater_fraction(links, regions_offshore):
    links = links.copy()
    links.loc[:, "underwater_fraction"] = (
        links.intersection(regions_offshore.union_all()).to_crs(DISTANCE_CRS).length
        / links.to_crs(DISTANCE_CRS).length
    ).round(2)

    return links


def create_unique_ids(df):
    """
    Create unique IDs for each project, starting with the index and adding a
    two-digit numerical suffix "-01", "-02", etc. only if there are multiple
    geometries for the same project.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        df (pd.DataFrame): An indexed DataFrame with unique IDs for each project.
    """
    df = df.copy().reset_index()

    # Count the occurrences of each 'pci_code'
    counts = df["id"].value_counts()

    # Generate cumulative counts within each group
    df["count"] = df.groupby("id").cumcount() + 1  # Start counting from 1, not 0

    # Add a two-digit suffix if the id appears more than once
    df["suffix"] = df.apply(
        lambda row: (
            f"-{str(row['count']).zfill(2)}"
            if counts[row["id"]] > 1
            else ""
        ),
        axis=1,
    )

    # Create the 'id' 
    df["id"] = df["id"] + df["suffix"]

    # Clean up by dropping the helper columns
    df = df.drop(columns=["count", "suffix"])

    df.set_index("id", inplace=True)

    return df


def aggregate_links(gdf):

    crs = gdf.crs

    gdf['bus0'], gdf['bus1'] = zip(*gdf[['bus0', 'bus1']].apply(lambda x: sorted(x), axis=1))

    gdf = gdf.groupby(["label", "bus0", "bus1"]).agg(
        {
            "length": "max",
            "p_nom": "sum",
            "underwater_fraction": "mean",
            "geometry": "first",
        }
    ).reset_index()

    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs=crs)

    return gdf


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_european_co2_pipelines",
            clusters="adm",
            opts="",
            run="KN2045_Mix",
            configfiles=["config/config.nrw.yaml"],        
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    settings = snakemake.params.european_co2_pipelines
    kml_path = snakemake.input.kml
    length_factor = settings.get("length_factor", 1.25)

    # Import PyPSA-Eur regular data
    n = pypsa.Network(snakemake.input.network)
    buses_coords = n.buses.loc[n.buses.carrier=="AC", ["x", "y"]].copy()
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name").to_crs(DISTANCE_CRS)
    regions_offshore = gpd.read_file(snakemake.input.regions_offshore).set_index("name").to_crs(DISTANCE_CRS)
    scope = gpd.read_file(snakemake.input.scope).to_crs(DISTANCE_CRS)

    # Import KML file by layers
    kml_layers = fiona.listlayers(kml_path)
    gdfs = []

    for layer in kml_layers:
        df = gpd.read_file(kml_path, driver="KML", layer=layer)
        df["layer"] = layer
        gdfs.append(df)

    gdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True)).to_crs(DISTANCE_CRS)

    # Clean text fields
    gdf["Name"] = gdf["Name"].apply(clean_text)

    # Only keep the rows within scope
    scope_union = scope.union_all()
    gdf = gdf[gdf.geometry.apply(lambda x: x.intersects(scope_union))]

    # Pipeline processing 
    pipelines = gdf[
        (gdf["layer"] == PIPELINE_LABEL)
        & (gdf.geometry.type.isin(["LineString", "MultiLineString"]))
    ][PIPELINE_COLS].copy()

    # Remove z-dimension if present
    pipelines["geometry"] = pipelines["geometry"].apply(drop_z_dim)

    # Only keep lines that are not completely within onshore regions
    # regions_offshore_union = regions_offshore.union_all().buffer(3000) # buffer to 3000 m to avoid precision issues
    # pipelines = pipelines[
    #     pipelines.geometry.apply(
    #         lambda x: not x.within(regions_offshore_union)
    #     )
    # ]

    pipelines = explode_linestrings_to_segments(pipelines)
    pipelines = split_to_overpassing_segments(pipelines, regions_onshore)
    pipelines = pipelines.drop_duplicates(subset=["Name", "geometry"])

    # Create a gdf of all endpoints of pipelines
    offshore_endpoints = pipelines.geometry.apply(lambda x: MultiPoint([Point(x.coords[0]), Point(x.coords[-1])]))
    offshore_endpoints = gpd.GeoDataFrame(geometry=offshore_endpoints.explode().reset_index(drop=True), crs=pipelines.crs)
    
    # Keep only endpoints that are outside onshore regions
    regions_onshore_union_buffer = regions_onshore.union_all().buffer(5000) # buffer to 5000 m to avoid precision issues

    offshore_endpoints = offshore_endpoints[
        offshore_endpoints.geometry.apply(
            lambda x: not x.within(regions_onshore_union_buffer)
        )
    ]

    # Buffer by 10 meters
    offshore_endpoints = offshore_endpoints.buffer(10)

    # Keep all that have an intersection of exactly two pipelines
    sindex = pipelines.sindex
    offshore_endpoints = offshore_endpoints[
        offshore_endpoints.geometry.apply(
            lambda pt: sum(pipelines.iloc[list(sindex.query(pt))].intersects(pt)) == 2
        )
    ]
    # Remove duplicate geometries and reset index
    offshore_endpoints = offshore_endpoints.drop_duplicates().reset_index(drop=True)
    
    # add column "intersects_offshore_endpoints" that containts the index of offshore_endpoints that each pipeline intersects
    pipelines["intersects_offshore_endpoints"] = pipelines.geometry.apply(
        lambda line: offshore_endpoints.index[
            offshore_endpoints.geometry.apply(lambda pt: line.intersects(pt))
        ].tolist()
    )
    pipelines["intersects_offshore_endpoints"] = (
        pipelines["intersects_offshore_endpoints"]
        .apply(lambda x: int(x[0]) if len(x) > 0 else None)
    )

    b_merge_candidates = pipelines["intersects_offshore_endpoints"].notna()
    merge_candidates = pipelines.loc[b_merge_candidates].copy()

    # Drop b_merge_candidates from pipelines
    pipelines = pipelines.loc[~b_merge_candidates]

    # Group by group_id, intersects_offshore_endpoints and merge linestrings in geometry columns. Rest of columns take the first value
    merged = merge_candidates.groupby(["group_id", "intersects_offshore_endpoints", "Name"]).agg({
        "description": "first",
        "length": "sum",
        "geometry": safe_linemerge,
    }).reset_index()
    merged = gpd.GeoDataFrame(merged, geometry="geometry", crs=pipelines.crs)

    # Append merged back to pipelines
    pipelines = pd.concat([pipelines, merged], ignore_index=True)

    buses_co2_offshore = create_new_buses(
        pipelines,
        regions_onshore,
        scope,
        "AC",
    )
    # Append to existing buses
    buses_coords = pd.concat([buses_coords, buses_co2_offshore[["x", "y"]]])

    pipelines = map_to_closest_region(pipelines, buses_co2_offshore, max_distance=OFFSHORE_BUS_RADIUS)
    pipelines = map_to_closest_region(pipelines, regions_onshore, max_distance=COASTAL_DISTANCE)

    # Calculate haversine distances
    pipelines = calculate_haversine_distance(buses_coords, pipelines, length_factor)

    # Set underwater fraction
    pipelines = set_underwater_fraction(
        pipelines,
        regions_offshore,
    )
    
    ### Clean up
    # Drop rows that are na in bus0 or bus1
    pipelines = pipelines.dropna(subset=["bus0", "bus1"])

    # Drop lines that connect the same bus
    pipelines = pipelines[pipelines["bus0"] != pipelines["bus1"]]

    # Drop duplicates with same name, bus0, bus1 keeping the longest
    pipelines = pipelines.sort_values("length", ascending=False)
    pipelines = pipelines.drop_duplicates(subset=["Name", "bus0", "bus1"], keep="first")

    ### Capacities
    #  Map model names
    pipelines["label"] = pipelines["Name"].map(PIPELINE_LABELS)
    pipelines["mtpa"] = pipelines["label"].map(
        settings["transport_volume_mtpa"]
    )
    pipelines["p_nom"] = (
        pipelines["mtpa"] * 1e6 / 8760 / settings.get("mtpa_utilisation_factor", 1)
    ).round(0) # result is tonnes per hour

    pipelines = pipelines.dropna(subset=["bus0", "bus1", "p_nom"])
    pipelines = pipelines.reset_index(drop=True)

    pipelines = aggregate_links(pipelines)

    # Create unique IDs
    pipelines = create_unique_ids(pipelines.rename(columns={"label": "id"}))

    # Add missing columns
    pipelines["carrier"] = "CO2 pipeline"
    pipelines["underground"] = "t"


    ### STORES
    # Store processing
    stores = gdf[
        (gdf["layer"] == STORE_LABEL)
        & (gdf.geometry.type.isin(["Point", "MultiPoint"]))
    ][STORE_COLS].copy()

    # Remove z-dimension if present
    stores["geometry"] = stores["geometry"].apply(drop_z_dim)

    # Only keep stores that have a key in STORE_LABELS
    stores = stores[stores["Name"].isin(STORE_LABELS.keys())]

    # Map model names
    stores["label"] = stores["Name"].map(STORE_LABELS)
    
    # Group by label
    stores = stores.groupby("label").agg({
        "geometry": "first",
    }).reset_index()
    stores = gpd.GeoDataFrame(stores, geometry="geometry", crs=gdf.crs)

    # Map to closest offshore region
    stores = map_to_closest_region(
        stores, buses_co2_offshore, max_distance=MAX_STORE_DISTANCE
    )
    stores = map_to_closest_region(
        stores, regions_offshore, max_distance=COASTAL_DISTANCE,
    )
    # Rename bus0 to bus
    stores = stores.rename(columns={"bus0": "bus"})
    stores = stores.drop(columns=["bus1"])
    stores["mtpa"] = stores["label"].map(
        settings["co2_sequestration_potential"]
    )
    stores["e_nom"] = stores["mtpa"] * 1e6  # in tonnes

    # Add missing columns
    stores["carrier"] = "co2 sequestered"
    
    # index, rename to id
    stores = stores.rename(columns={"label": "id"})
    stores.set_index("id", inplace=True)

    # Export
    buses_co2_offshore.to_csv(snakemake.output.buses_offshore, index=True)
    pipelines.to_csv(snakemake.output.links_co2_pipeline, index=True)
    stores.to_csv(snakemake.output.stores_co2, index=True)


    # # Debugging
    # map = regions_onshore.explore()
    # map = regions_offshore.explore(m=map, color = "lightblue")
    # map = pipelines.explore(m=map, color="purple")
    # map = buses_co2_offshore.explore(m=map, color="red")
    # map

    # # Create another gdf that has the x and y mapped from buses_coords
    # pipelines_ptp = pipelines.copy().to_crs(GEO_CRS)
    # pipelines_ptp["point0"] = pipelines_ptp.bus0.map(buses_coords.apply(lambda row: Point(row["x"], row["y"]), axis=1))
    # pipelines_ptp["point1"] = pipelines_ptp.bus1.map(buses_coords.apply(lambda row: Point(row["x"], row["y"]), axis=1))

    # # Map the remaining where na from buses_co2_offshore
    # missing_bus0 = pipelines_ptp["point0"].isna()
    # pipelines_ptp.loc[missing_bus0, "point0"] = pipelines_ptp.loc[missing_bus0].bus0.map(buses_co2_offshore.apply(lambda row: Point(row["x"], row["y"]), axis=1))
    # missing_bus1 = pipelines_ptp["point1"].isna()
    # pipelines_ptp.loc[missing_bus1, "point1"] = pipelines_ptp.loc[missing_bus1].bus1.map(buses_co2_offshore.apply(lambda row: Point(row["x"], row["y"]), axis=1))   
    # pipelines_ptp["geometry"] = pipelines_ptp.apply(
    #     lambda row: LineString([row["point0"], row["point1"]]), axis=1
    # )