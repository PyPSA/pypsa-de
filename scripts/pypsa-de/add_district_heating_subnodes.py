# -*- coding: utf-8 -*-


import logging

logger = logging.getLogger(__name__)
from random import randint
from time import sleep

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from typing import Union


def prepare_subnodes(
    subnodes: pd.DataFrame,
    cities: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    lau: gpd.GeoDataFrame,
    head: Union[int, bool] = 40,
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
    head : Union[int, bool], optional
        Number of largest district heating networks to keep. Defaults to 40. If set to True, it will be set to 40.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with processed subnodes, including geometries, clusters, LAU IDs, and NUTS3 shapes.
    """
    # If head is boolean set it to 40 for default behavior
    if isinstance(head, bool):
        head = 40

    subnodes["Stadt"] = subnodes["Stadt"].str.split("_").str[0]

    # Drop duplicates if Gelsenkirchen, Kiel, or Flensburg is included and keep the one with higher Wärmeeinspeisung in GWh/a
    subnodes = subnodes.drop_duplicates(subset="Stadt", keep="first")

    # Keep only n largest district heating networks according to head parameter
    subnodes = subnodes.sort_values(
        by="Wärmeeinspeisung in GWh/a", ascending=False
    ).head(head)

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
    subnodes["lau"] = subnodes.apply(
        lambda x: lau.loc[lau.geometry.contains(x.geometry).idxmax(), "LAU_ID"], axis=1
    )
    subnodes["lau_shape"] = subnodes.apply(
        lambda x: lau.loc[lau.geometry.contains(x.geometry).idxmax(), "geometry"].wkt,
        axis=1,
    )

    return subnodes


def add_subnodes(n: pypsa.Network, subnodes: gpd.GeoDataFrame) -> None:
    """
    Add largest district heating systems subnodes to the network.

    They are initialized with:
     - the total annual heat demand taken from the mother node, that is assigned to urban central heat and low-temperature heat for industry,
     - the heat demand profiles taken from the mother node,
     - the district heating investment options (stores, links) from the mother node,
     - and heat vents as generator components.
    The district heating loads in the mother nodes are reduced accordingly.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which subnodes will be added.
    subnodes : gpd.GeoDataFrame
        GeoDataFrame containing information about district heating subnodes.

    Returns
    -------
    None
    """

    # Add subnodes to network
    for _, subnode in subnodes.iterrows():
        name = f'{subnode["cluster"]} {subnode["Stadt"]} urban central heat'

        # Add buses
        n.madd(
            "Bus",
            [name],
            y=subnode.geometry.y,
            x=subnode.geometry.x,
            country="DE",
            location=f"{subnode['cluster']} {subnode['Stadt']}",
            carrier="urban central heat",
            unit="MWh_th",
        )

        # Get heat loads for urban central heat and low-temperature heat for industry
        uch_load_cluster = (
            n.snapshot_weightings.generators
            @ n.loads_t.p_set[f"{subnode['cluster']} urban central heat"]
        )
        lti_load_cluster = (
            n.loads.loc[
                f"{subnode['cluster']} low-temperature heat for industry", "p_set"
            ]
            * 8760
        )

        # Calculate share of low-temperature heat for industry in total district heating load of cluster
        dh_load_cluster = uch_load_cluster + lti_load_cluster
        lti_share = lti_load_cluster / dh_load_cluster

        # Calculate demand ratio between load of subnode according to Fernwärmeatlas and remaining load of assigned cluster
        demand_ratio = min(
            1,
            (subnode["yearly_heat_demand_MWh"] / dh_load_cluster),
        )

        lost_load = subnode["yearly_heat_demand_MWh"] - dh_load_cluster

        # District heating demand exceeding the original cluster load is disregarded
        if demand_ratio == 1:
            logger.info(
                f"District heating load of {subnode['Stadt']} exceeds load of its assigned cluster {subnode['cluster']}. {lost_load} MWh/a are disregarded."
            )

        # Add load components to subnode preserving the share of low-temperature heat for industry of the cluster
        uch_load = (
            demand_ratio
            * (1 - lti_share)
            * n.loads_t.p_set.filter(
                regex=f"{subnode['cluster']} urban central heat"
            ).rename(
                {
                    f"{subnode['cluster']} urban central heat": f"{subnode['cluster']} {subnode['Stadt']} urban central heat"
                },
                axis=1,
            )
        )
        n.madd(
            "Load",
            [name],
            bus=name,
            p_set=uch_load,
            carrier="urban central heat",
            location=f"{subnode['cluster']} {subnode['Stadt']}",
        )

        lti_load = (
            demand_ratio
            * lti_share
            * n.loads.filter(
                regex=f"{subnode['cluster']} low-temperature heat for industry", axis=0
            ).p_set.rename(
                {
                    f"{subnode['cluster']} low-temperature heat for industry": f"{subnode['cluster']} {subnode['Stadt']} low-temperature heat for industry"
                },
                axis=0,
            )
        )
        n.madd(
            "Load",
            [
                f"{subnode['cluster']} {subnode['Stadt']} low-temperature heat for industry"
            ],
            bus=name,
            p_set=lti_load,
            carrier="low-temperature heat for industry",
            location=f"{subnode['cluster']} {subnode['Stadt']}",
        )

        # Adjust loads of cluster buses
        n.loads_t.p_set.loc[
            :, f'{subnode["cluster"]} urban central heat'
        ] *= 1 - demand_ratio * (1 - lti_share)
        n.loads.loc[
            f'{subnode["cluster"]} low-temperature heat for industry', "p_set"
        ] *= (1 - demand_ratio * lti_share)

        # Replicate district heating stores and links of mother node for subnodes
        n.madd(
            "Bus",
            [f"{subnode['cluster']} {subnode['Stadt']} urban central water tanks"],
            location=f"{subnode['cluster']} {subnode['Stadt']}",
            carrier="urban central water tanks",
            unit="MWh_th",
        )

        stores = (
            n.stores.filter(like=f"{subnode['cluster']} urban central", axis=0)
            .reset_index()
            .replace(
                {
                    f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                },
                regex=True,
            )
            .set_index("Store")
        )
        n.madd("Store", stores.index, **stores)

        links = (
            n.links.loc[~n.links.carrier.str.contains("heat pump")]
            .filter(like=f"{subnode['cluster']} urban central", axis=0)
            .reset_index()
            .replace(
                {
                    f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                },
                regex=True,
            )
            .set_index("Link")
        )
        n.madd("Link", links.index, **links)

        # Add heat pumps to subnode
        heat_pumps = (
            n.links.filter(
                regex=f"{subnode['cluster']} urban central.*heat pump", axis=0
            )
            .reset_index()
            .replace(
                {
                    f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                },
                regex=True,
            )
            .set_index("Link")
        ).drop("efficiency", axis=1)
        heat_pumps_t = n.links_t.efficiency.filter(
            regex=f"{subnode['cluster']} urban central.*heat pump"
        )
        heat_pumps_t.columns = heat_pumps_t.columns.str.replace(
            f"{subnode['cluster']} urban central",
            f"{subnode['cluster']} {subnode['Stadt']} urban central",
        )
        n.madd("Link", heat_pumps.index, efficiency=heat_pumps_t, **heat_pumps)

        # Add heat vent to subnode
        n.madd(
            "Generator",
            [f"{name} heat vent"],
            bus=name,
            location=f"{subnode['cluster']} {subnode['Stadt']}",
            carrier="urban central heat vent",
            p_nom_extendable=True,
            p_min_pu=-1,
            p_max_pu=0,
            unit="MWh_th",
        )

    return


def extend_cops(cops: xr.DataArray, subnodes: gpd.GeoDataFrame) -> xr.DataArray:
    """
    Extend COPs (Coefficient of Performance) by subnodes mirroring the timeseries of the corresponding
    mother node.

    Parameters
    ----------
    cops : xr.DataArray
        DataArray containing COP timeseries data.
    subnodes : gpd.GeoDataFrame
        GeoDataFrame containing information about district heating subnodes.

    Returns
    -------
    xr.DataArray
        Extended DataArray with COP timeseries for subnodes.
    """
    cops_extended = cops.copy()

    # Iterate over the DataFrame rows
    for _, subnode in subnodes.iterrows():
        cluster_name = subnode["cluster"]
        city_name = subnode["Stadt"]

        # Select the xarray entry where name matches the cluster
        selected_entry = cops.sel(name=cluster_name)

        # Rename the selected entry
        renamed_entry = selected_entry.assign_coords(name=f"{cluster_name} {city_name}")

        # Combine the renamed entry with the extended dataset
        cops_extended = xr.concat([cops_extended, renamed_entry], dim="name")

    # Change dtype of the name dimension to string
    cops_extended.coords["name"] = cops_extended.coords["name"].astype(str)

    return cops_extended


def extend_heating_distribution(
    existing_heating_distribution: pd.DataFrame, subnodes: gpd.GeoDataFrame
) -> pd.DataFrame:
    """
    Extend heating distribution by subnodes mirroring the distribution of the
    corresponding mother node.

    Parameters
    ----------
    existing_heating_distribution : pd.DataFrame
        DataFrame containing the existing heating distribution.
    subnodes : gpd.GeoDataFrame
        GeoDataFrame containing information about district heating subnodes.

    Returns
    -------
    pd.DataFrame
        Extended DataFrame with heating distribution for subnodes.
    """
    # Merge the existing heating distribution with subnodes on the cluster name
    mother_nodes = (
        existing_heating_distribution.loc[subnodes.cluster.unique()]
        .unstack(-1)
        .to_frame()
    )
    cities_within_cluster = subnodes.groupby("cluster")["Stadt"].apply(list)
    mother_nodes["cities"] = mother_nodes.apply(
        lambda i: cities_within_cluster[i.name[2]], axis=1
    )
    # Explode the list of cities
    mother_nodes = mother_nodes.explode("cities")

    # Reset index to temporarily flatten it
    mother_nodes_reset = mother_nodes.reset_index()

    # Append city name to the third level of the index
    mother_nodes_reset["name"] = (
        mother_nodes_reset["name"] + " " + mother_nodes_reset["cities"]
    )

    # Set the index back
    mother_nodes = mother_nodes_reset.set_index(["heat name", "technology", "name"])

    # Drop the temporary 'cities' column
    mother_nodes.drop("cities", axis=1, inplace=True)

    # Reformat to match the existing heating distribution
    mother_nodes = mother_nodes.squeeze().unstack(-1).T

    # Combine the exploded data with the existing heating distribution
    existing_heating_distribution_extended = pd.concat(
        [existing_heating_distribution, mother_nodes]
    )
    return existing_heating_distribution_extended


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_district_heating_subnodes",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2020",
            run="KN2045_Bal_v4",
        )

    logger.info("Adding SysGF-specific functionality")

    n = pypsa.Network(snakemake.input.network)

    lau = gpd.read_file(
        f"{snakemake.input.lau}!LAU_RG_01M_2021_3035.geojson",
        crs="EPSG:3035",
    ).to_crs("EPSG:4326")

    fernwaermeatlas = pd.read_excel(
        snakemake.input.fernwaermeatlas,
        sheet_name="Fernwärmeatlas_öffentlich",
    )
    cities = gpd.read_file(snakemake.input.cities)
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    subnodes = prepare_subnodes(
        fernwaermeatlas,
        cities,
        regions_onshore,
        lau,
        head=snakemake.params.district_heating["add_subnodes"],
    )
    subnodes.to_file(snakemake.output.district_heating_subnodes, driver="GeoJSON")

    add_subnodes(n, subnodes)

    if snakemake.config["foresight"] == "myopic":
        cops = xr.open_dataarray(snakemake.input.cop_profiles)
        cops_extended = extend_cops(cops, subnodes)
        cops_extended.to_netcdf(snakemake.output.cop_profiles_extended)

    if snakemake.wildcards.planning_horizons == str(snakemake.params["baseyear"]):
        existing_heating_distribution = pd.read_csv(
            snakemake.input.existing_heating_distribution,
            header=[0, 1],
            index_col=0,
        )
        existing_heating_distribution_extended = extend_heating_distribution(
            existing_heating_distribution, subnodes
        )
        existing_heating_distribution_extended.to_csv(
            snakemake.output.existing_heating_distribution_extended
        )
    else:
        # write empty file to output
        with open(snakemake.output.existing_heating_distribution_extended, "w") as f:
            pass
    n.export_to_netcdf(snakemake.output.network)
