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

import pandas as pd
import pyproj
from shapely import wkt
from shapely.geometry import LineString, Point
from shapely.ops import transform

paths = ["workflow/submodules/pypsa-eur/scripts", "../submodules/pypsa-eur/scripts"]
for path in paths:
    sys.path.insert(0, os.path.abspath(path))
from cluster_gas_network import (
    build_clustered_gas_network,
    load_bus_regions,
    reindex_pipes,
)

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


def aggregate_parallel_pipes(df):
    strategies = {
        "bus0": "first",
        "bus1": "first",
        "p_nom": "sum",
        "p_nom_diameter": "sum",
        "max_pressure_bar": "mean",
        "build_year": "mean",
        "diameter_mm": "mean",
        "length": "mean",
        "name": " ".join,
        "p_min_pu": "min",
        "removed_gas_cap": "sum",
    }
    return df.groupby(df.index).agg(strategies)


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "cluster_wasserstoff_kernnetz",
            simpl="",
            clusters=22,
        )

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    fn = snakemake.input.cleaned_h2_network
    df = pd.read_csv(fn, index_col=0)
    for col in ["point0", "point1", "geometry"]:
        df[col] = df[col].apply(wkt.loads)

    bus_regions = load_bus_regions(
        snakemake.input.regions_onshore, snakemake.input.regions_offshore
    )

    kernnetz_cf = snakemake.params.kernnetz
    if kernnetz_cf["divide_pipes"]:
        segment_length = kernnetz_cf["pipes_segment_length"]
        df = divide_pipes(df, segment_length=segment_length)

    wasserstoff_kernnetz = build_clustered_gas_network(df, bus_regions)

    if not wasserstoff_kernnetz.empty:
        wasserstoff_kernnetz[["bus0", "bus1"]] = (
            wasserstoff_kernnetz[["bus0", "bus1"]]
            .apply(sorted, axis=1)
            .apply(pd.Series)
        )

        reindex_pipes(wasserstoff_kernnetz, prefix="H2 pipeline")

        wasserstoff_kernnetz["p_min_pu"] = 0
        wasserstoff_kernnetz["p_nom_diameter"] = 0
        wasserstoff_kernnetz = aggregate_parallel_pipes(wasserstoff_kernnetz)

    wasserstoff_kernnetz.to_csv(snakemake.output.clustered_h2_network)
