# -*- coding: utf-8 -*-
import logging

logger = logging.getLogger(__name__)

import geopandas as gpd
import pandas as pd
import numpy as np
from typing import Union
import xarray as xr
import shapely
import rasterio
from atlite.gis import ExclusionContainer
from atlite.gis import shape_availability
from tqdm import tqdm
import sys
import os

path = "dev/pypsa-de"
sys.path.insert(0, os.path.abspath(path))

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)


# Function to encode city names in UTF-8
def encode_utf8(city_name: str) -> bytes:
    """
    Encode a city name as a UTF-8 byte string.

    Parameters
    ----------
    city_name : str
        The name of the city to be encoded.

    Returns
    -------
    bytes
        The UTF-8 encoded byte string of the city name.
    """
    return city_name.encode("utf-8")


def process_eligible_points(
    x_coords: list[float],
    y_coords: list[float],
    values: np.ndarray,
    crs: str,
    lau_shapes: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Process points to create eligible area geometries for PTES potential assessment.

    Parameters
    ----------
    x_coords : list[float]
        X-coordinates of eligible points.
    y_coords : list[float]
        Y-coordinates of eligible points.
    values : np.ndarray
        Values associated with each point (eligibility flags).
    crs : str
        Coordinate reference system of the input points.
    lau_shapes : gpd.GeoDataFrame
        GeoDataFrame containing LAU (Local Administrative Unit) shapes to intersect with eligible areas.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing the intersection of eligible areas with LAU shapes.
    """
    eligible_areas = pd.DataFrame({"x": x_coords, "y": y_coords, "eligible": values})
    eligible_areas = gpd.GeoDataFrame(
        eligible_areas,
        geometry=gpd.points_from_xy(eligible_areas.x, eligible_areas.y),
        crs=crs,
    )

    chunk_size = 100000  # Adjust based on your memory constraints
    buffers = []

    for i in range(0, len(eligible_areas), chunk_size):
        chunk = eligible_areas.iloc[i : i + chunk_size].copy()
        chunk["geometry"] = chunk.geometry.buffer(5, cap_style="square")
        buffers.append(chunk)

    eligible_areas = pd.concat(buffers, ignore_index=True)

    # Use spatial indexing for more efficient overlay
    merged_data = eligible_areas.union_all()
    result = gpd.GeoDataFrame(geometry=[merged_data], crs=eligible_areas.crs)
    result = result.explode(index_parts=False).reset_index(drop=True)

    # Overlay with dh_systems using spatial indexing
    return gpd.overlay(result, lau_shapes, how="intersection")


def prepare_subnodes(
    subnodes: pd.DataFrame,
    cities: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    lau: gpd.GeoDataFrame,
    heat_techs: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Prepare subnodes by filtering district heating systems data for largest systems and assigning the corresponding LAU and onshore region shapes.
    Parameters
    ----------
    subnodes : pd.DataFrame
        DataFrame containing information about district heating systems.
    cities : gpd.GeoDataFrame
        GeoDataFrame containing city coordinates with columns 'Stadt' and 'geometry'.
    regions_onshore : gpd.GeoDataFrame
        GeoDataFrame containing onshore region geometries of clustered network.
    lau : gpd.GeoDataFrame
        GeoDataFrame containing LAU (Local Administrative Units) geometries and IDs.
    heat_techs : gpd.GeoDataFrame
        GeoDataFrame containing NUTS3 region geometries of heat technologies and data from eGo^N project.
    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with processed subnodes, including geometries, clusters, LAU IDs, and NUTS3 shapes.
    """

    subnodes["Stadt"] = subnodes["Stadt"].str.split("_").str[0]

    # Drop duplicates if Gelsenkirchen, Kiel, or Flensburg is included and keep the one with higher Wärmeeinspeisung in GWh/a
    subnodes = subnodes.drop_duplicates(subset="Stadt", keep="first")

    subnodes["yearly_heat_demand_MWh"] = subnodes["Wärmeeinspeisung in GWh/a"] * 1e3

    logger.info(
        f"The selected district heating networks have an overall yearly heat demand of {subnodes['yearly_heat_demand_MWh'].sum()} MWh/a. "
    )

    subnodes["geometry"] = subnodes["Stadt"].apply(
        lambda s: cities.loc[cities["Stadt"] == s, "geometry"].values[0]
    )

    subnodes = subnodes.dropna(subset=["geometry"])
    # Convert the DataFrame to a GeoDataFrame
    subnodes = gpd.GeoDataFrame(subnodes, crs="EPSG:4326")

    # Assign cluster to subnodes according to onshore regions
    subnodes["cluster"] = subnodes.apply(
        lambda x: regions_onshore.geometry.contains(x.geometry).idxmax(), axis=1
    )
    # For cities that are assigned to onshore regions outside Germany assign closest German cluster
    subnodes.loc[~subnodes.cluster.str.contains("DE"), "cluster"] = subnodes.loc[
        ~subnodes.cluster.str.contains("DE")
    ].apply(
        lambda x: (
            regions_onshore.filter(like="DE", axis=0)
            .geometry.distance(x.geometry)
            .idxmin()
        ),
        axis=1,
    )
    subnodes["lau"] = subnodes.apply(
        lambda x: lau.loc[lau.geometry.contains(x.geometry).idxmax(), "LAU_ID"], axis=1
    )
    subnodes["lau_shape"] = subnodes.apply(
        lambda x: lau.loc[lau.geometry.contains(x.geometry).idxmax(), "geometry"].wkt,
        axis=1,
    )
    subnodes["nuts3"] = subnodes.apply(
        lambda x: heat_techs.geometry.contains(x.geometry).idxmax(),
        axis=1,
    )
    subnodes["nuts3_shape"] = subnodes.apply(
        lambda x: heat_techs.loc[
            heat_techs.geometry.contains(x.geometry).idxmax(), "geometry"
        ].wkt,
        axis=1,
    )

    return subnodes


def add_ptes_limit(
    subnodes: gpd.GeoDataFrame,
    osm_land_cover: rasterio.io.DatasetReader,
    natura: rasterio.io.DatasetReader,
    groundwater: xr.Dataset,
    codes: list,
    max_groundwater_depth: float,
    ptes_potential_scalar: float,
) -> gpd.GeoDataFrame:
    """
    Add PTES limit to subnodes according to land availability within city regions.

    Parameters
    ----------
    subnodes : gpd.GeoDataFrame
        GeoDataFrame containing information about district heating subnodes.
    osm_land_cover : rasterio.io.DatasetReader
        OSM land cover raster dataset.
    natura : rasterio.io.DatasetReader
        NATURA 2000 protected areas raster dataset.
    groundwater : xr.Dataset
        Groundwater depth dataset.
    codes : list
        List of CORINE land cover codes to include.
    max_groundwater_depth : float
        Maximum allowable groundwater depth for PTES installation.
    ptes_potential_scalar : float
        Scalar to adjust PTES potential.

    Returns
    -------
    gpd.GeoDataFrame
        Updated GeoDataFrame with PTES potential added.
    """
    dh_systems = subnodes.copy()
    dh_systems["lau_shape"] = dh_systems["lau_shape"].apply(shapely.wkt.loads)
    dh_systems = dh_systems.set_geometry("lau_shape")
    dh_systems.crs = "EPSG:4326"
    dh_systems = dh_systems.to_crs(3035)

    # Process in batches by region/cluster
    batch_results = []
    batch_size = 1  # Define the batch size

    for i in tqdm(
        range(0, len(dh_systems), batch_size), desc="Processing LAU shapes", unit="area"
    ):
        batch = dh_systems.iloc[i : i + batch_size]

        # Create exclusion container for this batch
        excluder = ExclusionContainer(crs=3035, res=10)
        excluder.add_raster(osm_land_cover, codes=codes, invert=True, crs=3035)
        excluder.add_raster(natura, codes=[1], invert=True, crs=3035)

        # Process just this batch
        batch_shapes = batch["lau_shape"]
        band, transform = shape_availability(batch_shapes, excluder)
        masked_data = band
        row_indices, col_indices = np.where(masked_data != osm_land_cover.nodata)
        values = masked_data[row_indices, col_indices]

        x_coords, y_coords = rasterio.transform.xy(transform, row_indices, col_indices)

        # Create and process eligible areas for this batch
        if len(x_coords) > 0:  # Only process if there are eligible areas
            eligible_areas = process_eligible_points(
                x_coords, y_coords, values, osm_land_cover.crs, batch
            )
            batch_results.append(eligible_areas)

        # Clear memory
        del band, masked_data, excluder
        import gc

        gc.collect()

    # Combine results from all batches
    if batch_results:
        eligible_areas = pd.concat(batch_results, ignore_index=True)

    eligible_areas = gpd.sjoin(
        eligible_areas, dh_systems.drop("Stadt", axis=1), how="left", rsuffix=""
    )[["Stadt", "geometry"]].set_geometry("geometry")

    # filter for eligible areas that are larger than 10000 m^2
    eligible_areas = eligible_areas[eligible_areas.area > 10000]

    # Find closest value in groundwater dataset and kick out areas with groundwater level > threshold
    eligible_areas["groundwater_level"] = eligible_areas.to_crs("EPSG:4326").apply(
        lambda a: groundwater.sel(
            lon=a.geometry.centroid.x, lat=a.geometry.centroid.y, method="nearest"
        )["WTD"].values[0],
        axis=1,
    )
    eligible_areas = eligible_areas[
        eligible_areas.groundwater_level < max_groundwater_depth
    ]

    # Combine eligible areas by city
    eligible_areas = eligible_areas.dissolve("Stadt")

    # Calculate PTES potential according to Dronninglund and DEA parameters
    eligible_areas["area_m2"] = eligible_areas.area
    eligible_areas["nstorages_pot"] = eligible_areas.area_m2 / 10000
    eligible_areas["storage_pot_mwh"] = eligible_areas["nstorages_pot"] * 4500

    subnodes.set_index("Stadt", inplace=True)
    subnodes["ptes_pot_mwh"] = (
        eligible_areas.loc[subnodes.index.intersection(eligible_areas.index)][
            "storage_pot_mwh"
        ]
        * ptes_potential_scalar
    )
    subnodes["ptes_pot_mwh"] = subnodes["ptes_pot_mwh"].fillna(0)
    subnodes.reset_index(inplace=True)

    return subnodes


def extend_regions_onshore(
    regions_onshore: gpd.GeoDataFrame,
    subnodes_all: gpd.GeoDataFrame,
    head: Union[int, bool] = 40,
) -> Union[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Extend onshore regions to include city LAU regions and restrict geometries.

    Parameters
    ----------
    regions_onshore : geopandas.GeoDataFrame
        GeoDataFrame containing the onshore regions with geometries.
    subnodes_all : pandas.DataFrame
        DataFrame containing information about subnodes, including city names,
        clusters, and LAU shapes in WKT format.
    head : int or bool, optional
        Number of top cities to include based on "Wärmeeinspeisung in GWh/a".
        If set to True, defaults to 40. Default is 40.

    Returns
    -------
    regions_onshore_extended : geopandas.GeoDataFrame
        GeoDataFrame with extended onshore regions including city LAU regions.
    regions_onshore_restricted : geopandas.GeoDataFrame
        GeoDataFrame with restricted onshore regions, replacing geometries
        with the remaining city areas.
    """
    subnodes_all["lau_shape"] = subnodes_all["lau_shape"].apply(shapely.wkt.loads)
    if isinstance(head, bool):
        head = 40
    # Extend regions_onshore to include the cities' lau regions
    subnodes = subnodes_all.sort_values(
        by="Wärmeeinspeisung in GWh/a", ascending=False
    ).head(head)[["Stadt", "cluster", "lau_shape"]]

    subnodes = gpd.GeoDataFrame(subnodes, crs="EPSG:4326", geometry="lau_shape")
    # Create column name that comprises cluster and Stadt
    subnodes["name"] = subnodes["cluster"] + " " + subnodes["Stadt"]
    # Crop city regions from onshore regions
    regions_onshore["geometry"] = regions_onshore.geometry.difference(
        subnodes.union_all()
    )

    # Rename lau_shape to geometry
    subnodes = subnodes.rename(columns={"lau_shape": "geometry"}).drop(
        columns=["Stadt", "cluster"]
    )
    # Concat regions_onshore and subnodes
    regions_onshore_extended = pd.concat([regions_onshore, subnodes.set_index("name")])

    # Restrict regions_onshore geometries to only consist of the remaining city areas
    subnodes_rest = subnodes_all[~subnodes_all["Stadt"].isin(subnodes["name"])]

    subnodes_rest_dissolved = subnodes_rest.set_geometry("lau_shape").dissolve(
        "cluster"
    )
    # regions_onshore_restricted should replace geometries of regions_onshore with the geometries of subnodes_rest
    regions_onshore_restricted = regions_onshore_extended.copy()
    regions_onshore_restricted.loc[subnodes_rest_dissolved.index, "geometry"] = (
        subnodes_rest_dissolved["lau_shape"]
    )

    return regions_onshore_extended, regions_onshore_restricted


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        path = "../../"
        sys.path.insert(0, os.path.abspath(path))
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_district_heating_subnodes",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2045",
            run="LowGroundWaterDepth",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    logger.info("Adding SysGF-specific functionality")

    heat_techs = gpd.read_file(snakemake.input.heating_technologies_nuts3).set_index(
        "index"
    )
    lau = gpd.read_file(
        # "/home/cpschau/Code/dev/pypsa-ariadne/.snakemake/storage/http/gisco-services.ec.europa.eu/distribution/v2/lau/download/ref-lau-2021-01m.geojson/LAU_RG_01M_2021_3035.geojson",
        f"{snakemake.input.lau}!LAU_RG_01M_2021_3035.geojson",
        crs="EPSG:3035",
    ).to_crs("EPSG:4326")

    fernwaermeatlas = pd.read_excel(
        snakemake.input.fernwaermeatlas,
        sheet_name="Fernwärmeatlas_öffentlich",
    )
    cities = gpd.read_file(snakemake.input.cities)
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    # Assign onshore region to heat techs based on geometry
    heat_techs["cluster"] = heat_techs.apply(
        lambda x: regions_onshore.geometry.contains(x.geometry).idxmax(),
        axis=1,
    )

    subnodes = prepare_subnodes(
        fernwaermeatlas,
        cities,
        regions_onshore,
        lau,
        heat_techs,
    )

    # Add PTES limit to subnodes according to land availability within city regions
    osm_land_cover = rasterio.open(snakemake.input.osm_land_cover)
    natura = rasterio.open(snakemake.input.natura)
    groundwater = xr.open_dataset(snakemake.input.groundwater_depth).sel(
        lon=slice(subnodes["geometry"].x.min(), subnodes["geometry"].x.max()),
        lat=slice(subnodes["geometry"].y.min(), subnodes["geometry"].y.max()),
    )
    subnodes = add_ptes_limit(
        subnodes,
        osm_land_cover,
        natura,
        groundwater,
        snakemake.params.district_heating["osm_landcover_codes"],
        snakemake.params.district_heating["max_groundwater_depth"],
        snakemake.params.district_heating["ptes_potential_scalar"],
    )

    subnodes.to_file(snakemake.output.district_heating_subnodes, driver="GeoJSON")

    regions_onshore_extended, regions_onshore_restricted = extend_regions_onshore(
        regions_onshore,
        subnodes,
        head=snakemake.params.district_heating["add_subnodes"],
    )

    regions_onshore_extended.to_file(
        snakemake.output.regions_onshore_extended, driver="GeoJSON"
    )

    regions_onshore_restricted.to_file(
        snakemake.output.regions_onshore_restricted, driver="GeoJSON"
    )
