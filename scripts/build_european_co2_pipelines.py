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
from shapely import segmentize, unary_union
from shapely.algorithms.polylabel import polylabel
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point
from shapely.ops import nearest_points

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)


logger = logging.getLogger(__name__)

CLUSTER_TOL = 25000 # in meters
DISTANCE_CRS = "EPSG:3035"
GEO_CRS = "EPSG:4326"
OFFSHORE_BUS_RADIUS = 5000 # in meters
PIPELINE_LABEL = "Infrastruktur"


def create_new_buses(
    gdf: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    scope: gpd.GeoDataFrame,
    carrier: str,
    regions_onshore_buffer: int = 1000,
    distance_crs: str = DISTANCE_CRS,
    geo_crs: str = GEO_CRS,
    tol: int = CLUSTER_TOL,
    offset: int = 0,
):
    buffered_regions = (
        regions_onshore.to_crs(distance_crs)
        .buffer(regions_onshore_buffer)
        .to_crs(geo_crs)
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
        geometry=[geom for geom in list_points.geoms], crs=geo_crs
    )

    gdf_points["geometry"] = gdf_points.to_crs(distance_crs).buffer(tol).to_crs(geo_crs)

    # Aggregate rows with touching polygons
    gdf_points = gdf_points.dissolve()
    # split into separate polygons
    gdf_points = gdf_points.explode().reset_index(drop=True)

    gdf_points["poi"] = (
        gdf_points["geometry"]
        .to_crs(distance_crs)
        .apply(lambda polygon: polylabel(polygon, tolerance=tol / 2))
        .to_crs(geo_crs)
    )

    # Extract x and y coordinates into separate columns
    gdf_points["x"] = gdf_points["poi"].x
    gdf_points["y"] = gdf_points["poi"].y
    gdf_points["name"] = gdf_points.apply(
        lambda x: f"OFFSHORE {int(x.name)+1+offset}", axis=1
    )
    gdf_points.set_index("name", inplace=True)
    gdf_points["carrier"] = carrier

    return gdf_points[["x", "y", "carrier", "geometry"]]


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
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    regions_offshore = gpd.read_file(snakemake.input.regions_offshore).set_index("name")
    scope = gpd.read_file(snakemake.input.scope)

    # Import KML file by layers
    kml_layers = fiona.listlayers(kml_path)
    gdfs = []

    for layer in kml_layers:
        df = gpd.read_file(kml_path, driver="KML", layer=layer)
        df["layer"] = layer
        gdfs.append(df)

    gdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True))

    # Pipeline processing 
    pipelines = gdf[
        (gdf["layer"] == PIPELINE_LABEL)
        & (gdf.geometry.type.isin(["LineString", "MultiLineString"]))
    ].copy()

    buses_co2_offshore = create_new_buses(
        pipelines,
        regions_onshore,
        scope,
        "AC",
    )

    # Debugging
    map = regions_onshore.explore()
    map = regions_offshore.explore(m=map, color = "lightblue")
    map = pipelines.explore(m=map, color="purple")
    map = buses_co2_offshore.explore(m=map, color="red")
    map
