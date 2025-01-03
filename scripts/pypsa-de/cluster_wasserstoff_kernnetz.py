# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Cluster Wasserstoff Kernnetz to clustered model regions.
"""

import logging

logger = logging.getLogger(__name__)

import os
import sys

import geopandas as gpd
import pandas as pd
import pyproj
from pypsa.geo import haversine_pts
from shapely import wkt
from shapely.geometry import LineString, Point
from shapely.ops import transform

from scripts._helpers import configure_logging, mock_snakemake
from scripts.cluster_gas_network import load_bus_regions, reindex_pipes

# Define a function for projecting points to meters
project_to_meters = pyproj.Transformer.from_proj(
    pyproj.Proj("epsg:4326"),  # assuming WGS84
    pyproj.Proj(proj="utm", zone=33, ellps="WGS84"),  # adjust the projection as needed
    always_xy=True,
).transform

# Define a function for projecting points back to decimal degrees
project_to_degrees = pyproj.Transformer.from_proj(
    pyproj.Proj(proj="utm", zone=33, ellps="WGS84"),  # adjust the projection as needed
    pyproj.Proj("epsg:4326"),
    always_xy=True,
).transform


def split_line_by_length(line, segment_length_km):
    """
    Split a Shapely LineString into segments of a specified length.

    Parameters:
    - line (LineString): The original LineString to be split.
    - segment_length_km (float): The desired length of each resulting segment in kilometers.

    Returns:
    list: A list of Shapely LineString objects representing the segments.
    """

    # Convert segment length from kilometers to meters
    segment_length_meters = segment_length_km * 1000

    # Project the LineString to a suitable metric projection
    projected_line = transform(project_to_meters, line)

    total_length = projected_line.length
    num_segments = int(total_length / segment_length_meters)

    # Return early if no segmentation required
    if num_segments <= 1:
        return [line]

    segments = []
    for i in range(1, num_segments + 1):
        start_point = projected_line.interpolate((i - 1) * segment_length_meters)
        end_point = projected_line.interpolate(i * segment_length_meters)

        # Extract x and y coordinates from the tuples
        start_point_coords = (start_point.x, start_point.y)
        end_point_coords = (end_point.x, end_point.y)

        # Create Shapely Point objects
        start_point_degrees = Point(start_point_coords)
        end_point_degrees = Point(end_point_coords)

        # Project the points back to decimal degrees
        start_point_degrees = transform(project_to_degrees, start_point_degrees)
        end_point_degrees = transform(project_to_degrees, end_point_degrees)

        # last point without interpolation
        if i == num_segments:
            end_point_degrees = Point(line.coords[-1])

        segment = LineString([start_point_degrees, end_point_degrees])
        segments.append(segment)

    return segments


def divide_pipes(df, segment_length=10):
    """
    Divide a GeoPandas DataFrame of LineString geometries into segments of a
    specified length.

    Parameters:
    - df (GeoDataFrame): The input DataFrame containing LineString geometries.
    - segment_length (float): The desired length of each resulting segment in kilometers.

    Returns:
    GeoDataFrame: A new GeoDataFrame with additional rows representing the segmented pipes.
    """

    result = pd.DataFrame(columns=df.columns)

    for index, pipe in df.iterrows():
        segments = split_line_by_length(pipe.geometry, segment_length)

        for i, segment in enumerate(segments):
            res_row = pipe.copy()
            res_row.geometry = segment
            res_row.point0 = Point(segment.coords[0])
            res_row.point1 = Point(segment.coords[1])
            res_row.length_haversine = segment_length
            result.loc[f"{index}-{i}"] = res_row

    return result


def build_clustered_h2_network(
    df, bus_regions, recalculate_length=True, length_factor=1.25
):
    for i in [0, 1]:
        gdf = gpd.GeoDataFrame(geometry=df[f"point{i}"], crs="EPSG:4326")

        bus_mapping = gpd.sjoin(gdf, bus_regions, how="left", predicate="within")[
            "name"
        ]
        bus_mapping = bus_mapping.groupby(bus_mapping.index).first()

        df[f"bus{i}"] = bus_mapping

        df[f"point{i}"] = df[f"bus{i}"].map(
            bus_regions.to_crs(3035).centroid.to_crs(4326)
        )

    # drop pipes where not both buses are inside regions
    df = df.loc[~df.bus0.isna() & ~df.bus1.isna()]

    # drop pipes within the same region
    df = df.loc[df.bus1 != df.bus0]

    if df.empty:
        return df

    if recalculate_length:
        logger.info("Recalculating pipe lengths as center to center * length factor")
        # recalculate lengths as center to center * length factor
        df["length"] = df.apply(
            lambda p: length_factor
            * haversine_pts([p.point0.x, p.point0.y], [p.point1.x, p.point1.y]),
            axis=1,
        )

    # tidy and create new numbered index
    df.drop(["point0", "point1"], axis=1, inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def aggregate_parallel_pipes(df, aggregate_build_years="mean"):
    strategies = {
        "bus0": "first",
        "bus1": "first",
        "p_nom": "sum",
        "p_nom_diameter": "sum",
        "max_pressure_bar": "mean",
        "build_year": aggregate_build_years,
        "diameter_mm": "mean",
        "length": "mean",
        "name": " ".join,
        "p_min_pu": "min",
        "investment_costs (Mio. Euro)": "sum",
        "removed_gas_cap": "sum",
        "ipcei": " ".join,
        "pci": " ".join,
        "retrofitted": lambda x: (x.sum() / len(x))
        > 0.6,  # consider as retrofit if more than 60% of pipes are retrofitted (relevant for costs)
    }
    return df.groupby(df.index).agg(strategies)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "cluster_wasserstoff_kernnetz",
            simpl="",
            clusters=27,
            run="KN2045_Bal_v4",
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2020",
        )

    configure_logging(snakemake)

    fn = snakemake.input.cleaned_h2_network
    df = pd.read_csv(fn, index_col=0)
    for col in ["point0", "point1", "geometry"]:
        df[col] = df[col].apply(wkt.loads)

    bus_regions = load_bus_regions(
        snakemake.input.regions_onshore, snakemake.input.regions_offshore
    )
    logger.info(f"Clustering Wasserstoff Kernnetz for {list(bus_regions.index)}")
    kernnetz_cf = snakemake.params.kernnetz

    if kernnetz_cf["divide_pipes"]:
        segment_length = kernnetz_cf["pipes_segment_length"]
        df = divide_pipes(df, segment_length=segment_length)

    wasserstoff_kernnetz = build_clustered_h2_network(
        df,
        bus_regions,
        recalculate_length=kernnetz_cf["recalculate_length"],
        length_factor=1.25,
    )

    if kernnetz_cf["divide_pipes"] & (not kernnetz_cf["aggregate_parallel_pipes"]):
        # Set length to 0 for duplicates from the 2nd occurrence onwards and make name unique
        logger.info(
            f"Setting length to 0 for split pipes as Kernnetz pipes are segmented (divide pipes: {kernnetz_cf["divide_pipes"]}) and parallel pipes not aggregated (aggregate_parallel_pipes: {kernnetz_cf["aggregate_parallel_pipes"]})."
        )
        wasserstoff_kernnetz["occurrence"] = (
            wasserstoff_kernnetz.groupby("name").cumcount() + 1
        )
        wasserstoff_kernnetz.loc[wasserstoff_kernnetz["occurrence"] > 1, "length"] = 0
        wasserstoff_kernnetz["name"] = wasserstoff_kernnetz.apply(
            lambda row: (
                f"{row['name']}-split{row['occurrence']}"
                if row["occurrence"] > 1
                else row["name"]
            ),
            axis=1,
        )
        wasserstoff_kernnetz = wasserstoff_kernnetz.drop(columns="occurrence")

    if not wasserstoff_kernnetz.empty:
        wasserstoff_kernnetz[["bus0", "bus1"]] = (
            wasserstoff_kernnetz[["bus0", "bus1"]]
            .apply(sorted, axis=1)
            .apply(pd.Series)
        )

        wasserstoff_kernnetz["p_min_pu"] = 0
        wasserstoff_kernnetz["p_nom_diameter"] = 0

        if kernnetz_cf["aggregate_parallel_pipes"]:
            reindex_pipes(wasserstoff_kernnetz, prefix="H2 pipeline")
            wasserstoff_kernnetz = aggregate_parallel_pipes(
                wasserstoff_kernnetz, kernnetz_cf["aggregate_build_years"]
            )

        else:
            wasserstoff_kernnetz.index = wasserstoff_kernnetz.name.astype(str)

    wasserstoff_kernnetz.to_csv(snakemake.output.clustered_h2_network)
