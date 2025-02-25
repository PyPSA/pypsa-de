# -*- coding: utf-8 -*-
import logging

logger = logging.getLogger(__name__)

import geopandas as gpd
import pandas as pd
from typing import Union

import shapely


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
        subnodes.unary_union
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

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

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
