# -*- coding: utf-8 -*-
import logging

logger = logging.getLogger(__name__)

import geopandas as gpd
import pandas as pd
import pypsa
import xarray as xr
from typing import Union

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from scripts.prepare_network import maybe_adjust_costs_and_potentials
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)


def add_subnodes(
    n: pypsa.Network,
    subnodes: gpd.GeoDataFrame,
    cop: xr.DataArray,
    direct_heat_source_utilisation_profile: xr.DataArray,
    head: Union[int, bool] = 40,
) -> None:
    """
    Add largest district heating systems subnodes to the network.
    They are initialized with:
     - the total annual heat demand taken from the mother node, that is assigned to urban central heat and low-temperature heat for industry,
     - the heat demand profiles taken from the mother node,
     - the district heating investment options (stores, storage units, links, generators) from the mother node,
    The district heating loads in the mother nodes are reduced accordingly.
    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which subnodes will be added.
    subnodes : gpd.GeoDataFrame
        GeoDataFrame containing information about district heating subnodes.
    cop : xr.DataArray
        COPs for heat pumps.
    direct_heat_source_utilisation_profile : xr.DataArray
        Direct heat source utilisation profiles.
    Returns
    -------
    None
    """

    # If head is boolean set it to 40 for default behavior
    if isinstance(head, bool):
        head = 40

    # Keep only n largest district heating networks according to head parameter
    subnodes_head = subnodes.sort_values(
        by="Wärmeeinspeisung in GWh/a", ascending=False
    ).head(head)
    subnodes_head.to_file(snakemake.output.district_heating_subnodes, driver="GeoJSON")

    subnodes_rest = subnodes[~subnodes.index.isin(subnodes_head.index)]

    # Add subnodes to network
    for _, subnode in subnodes_head.iterrows():
        name = f'{subnode["cluster"]} {subnode["Stadt"]} urban central'
        location = f"{subnode['cluster']} {subnode['Stadt']}"

        # Add buses
        buses = (
            n.buses.filter(like=f"{subnode['cluster']} urban central", axis=0)
            .reset_index()
            .replace(
                {
                    f"{subnode['cluster']} urban central": name,
                    f"{subnode['cluster']}$": location,
                },
                regex=True,
            )
            .set_index("Bus")
        )
        n.add("Bus", buses.index, **buses)

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

        # Calculate demand ratio between load of subnode according to Fernwärmeatlas and remaining load of assigned cluster
        demand_ratio = min(
            1,
            (subnode["yearly_heat_demand_MWh"] / dh_load_cluster),
        )

        dh_load_cluster_subnodes = subnodes_head.loc[
            subnodes_head.cluster == subnode["cluster"], "yearly_heat_demand_MWh"
        ].sum()
        lost_load = dh_load_cluster_subnodes - dh_load_cluster

        # District heating demand from Fernwärmeatlas exceeding the original cluster load is disregarded. The shares of the subsystems are set according to Fernwärmeatlas, while the aggregate load of cluster is preserved.
        if dh_load_cluster_subnodes > dh_load_cluster:
            logger.warning(
                f"Aggregated district heating load of systems within {subnode['cluster']} exceeds load of cluster. {lost_load} MWh/a are disregarded."
            )
            demand_ratio *= subnode["yearly_heat_demand_MWh"] / dh_load_cluster_subnodes

        # Add load components to subnode preserving the share of low-temperature heat for industry of the cluster
        uch_load = demand_ratio * n.loads_t.p_set[
            f"{subnode['cluster']} urban central heat"
        ].rename(f"{subnode['cluster']} {subnode['Stadt']} urban central heat")

        n.add(
            "Load",
            f"{name} heat",
            bus=f"{name} heat",
            p_set=uch_load,
            carrier="urban central heat",
            location=f"{subnode['cluster']} {subnode['Stadt']}",
        )

        lti_load = (
            demand_ratio
            * n.loads.loc[
                f"{subnode['cluster']} low-temperature heat for industry", "p_set"
            ]
        )
        n.add(
            "Load",
            f"{subnode['cluster']} {subnode['Stadt']} low-temperature heat for industry",
            bus=f"{name} heat",
            p_set=lti_load,
            carrier="low-temperature heat for industry",
            location=location,
        )

        # Adjust loads of cluster buses
        n.loads_t.p_set.loc[:, f'{subnode["cluster"]} urban central heat'] *= (
            1 - demand_ratio
        )

        n.loads.loc[
            f'{subnode["cluster"]} low-temperature heat for industry', "p_set"
        ] *= (1 - demand_ratio)

        # Replicate district heating stores of mother node for subnodes
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

        # Restrict PTES capacity in subnodes if modeled as store
        if stores.carrier.str.contains("pits$").any():
            stores.loc[stores.carrier.str.contains("pits$").index, "e_nom_max"] = (
                subnode["ptes_pot_mwh"]
            )
        n.add("Store", stores.index, **stores)

        # Replicate district heating storage units of mother node for subnodes
        storage_units = (
            n.storage_units.filter(like=f"{subnode['cluster']} urban central", axis=0)
            .reset_index()
            .replace(
                {
                    f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                },
                regex=True,
            )
            .set_index("StorageUnit")
        )

        # Restrict PTES capacity in subnodes if modeled as storage unit
        if storage_units.carrier.str.contains("pits$").any():
            storage_units.loc[
                storage_units.carrier.str.contains("pits$"), "p_nom_max"
            ] = (subnode["ptes_pot_mwh"] / storage_units["max_hours"])
        n.add("StorageUnit", storage_units.index, **storage_units)

        # restrict PTES capacity in mother nodes
        mother_nodes_ptes_pot = subnodes_rest.groupby("cluster").ptes_pot_mwh.sum()

        mother_nodes_ptes_pot.index = (
            mother_nodes_ptes_pot.index + " urban central water pits"
        )

        if not n.storage_units.filter(like="urban central water pits", axis=0).empty:
            n.storage_units.loc[mother_nodes_ptes_pot.index, "p_nom_max"] = (
                mother_nodes_ptes_pot
                / n.storage_units.loc[mother_nodes_ptes_pot.index, "max_hours"]
            )
        elif not n.stores.filter(like="urban central water pits", axis=0).empty:
            n.stores.loc[mother_nodes_ptes_pot.index, "e_nom_max"] = (
                mother_nodes_ptes_pot
            )

        # Replicate district heating generators of mother node for subnodes
        generators = (
            n.generators.filter(like=f"{subnode['cluster']} urban central", axis=0)
            .reset_index()
            .replace(
                {
                    f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                },
                regex=True,
            )
            .set_index("Generator")
        )
        n.add("Generator", generators.index, **generators)

        # Replicate district heating links of mother node for subnodes with separate treatment for links with dynamic efficiencies
        links = (
            n.links.loc[~n.links.carrier.str.contains("heat pump|direct", regex=True)]
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
        n.add("Link", links.index, **links)

        # Add heat pumps and direct heat source utilization to subnode
        for heat_source in snakemake.params.heat_pump_sources:
            cop_heat_pump = (
                cop.sel(
                    heat_system="urban central",
                    heat_source=heat_source,
                    name=f"{subnode['cluster']} {subnode['Stadt']}",
                )
                .to_pandas()
                .to_frame(name=f"{name} {heat_source} heat pump")
                .reindex(index=n.snapshots)
                if snakemake.params.sector["time_dep_hp_cop"]
                else n.links.filter(like=heat_source, axis=0).efficiency.mode()
            )

            heat_pump = (
                n.links.filter(
                    regex=f"{subnode['cluster']} urban central.*{heat_source}.*heat pump",
                    axis=0,
                )
                .reset_index()
                .replace(
                    {
                        f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                    },
                    regex=True,
                )
                .drop(["efficiency", "efficiency2"], axis=1)
                .set_index("Link")
            )
            if heat_pump["bus2"].str.match("$").any():
                n.add("Link", heat_pump.index, efficiency=cop_heat_pump, **heat_pump)
            else:
                n.add(
                    "Link",
                    heat_pump.index,
                    efficiency=-(cop_heat_pump - 1),
                    efficiency2=cop_heat_pump,
                    **heat_pump,
                )

            if heat_source in snakemake.params.direct_utilisation_heat_sources:
                # Add direct heat source utilization to subnode
                efficiency_direct_utilisation = (
                    direct_heat_source_utilisation_profile.sel(
                        heat_source=heat_source,
                        name=f"{subnode['cluster']} {subnode['Stadt']}",
                    )
                    .to_pandas()
                    .to_frame(name=f"{name} {heat_source} heat direct utilisation")
                    .reindex(index=n.snapshots)
                )

                direct_utilization = (
                    n.links.filter(
                        regex=f"{subnode['cluster']} urban central.*{heat_source}.*direct",
                        axis=0,
                    )
                    .reset_index()
                    .replace(
                        {
                            f"{subnode['cluster']} urban central": f"{subnode['cluster']} {subnode['Stadt']} urban central"
                        },
                        regex=True,
                    )
                    .set_index("Link")
                    .drop("efficiency", axis=1)
                )

                n.add(
                    "Link",
                    direct_utilization.index,
                    efficiency=efficiency_direct_utilisation,
                    **direct_utilization,
                )

            # Restrict heat source potential in subnodes
            if heat_source in snakemake.params.district_heating["limited_heat_sources"]:
                # get potential
                p_max_source = pd.read_csv(
                    snakemake.input[heat_source],
                    index_col=0,
                ).squeeze()[f"{subnode['cluster']} {subnode['Stadt']}"]
                # add potential to generator
                n.generators.loc[
                    f"{subnode['cluster']} {subnode['Stadt']} urban central {heat_source} heat",
                    "p_nom_max",
                ] = p_max_source

    return


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

        from scripts._helpers import mock_snakemake

        # Change directory to this script directory
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        snakemake = mock_snakemake(
            "add_district_heating_subnodes",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2045",
            run="Baseline",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    logger.info("Adding SysGF-specific functionality")

    n = pypsa.Network(snakemake.input.network)

    lau = gpd.read_file(
        f"{snakemake.input.lau_regions}!LAU_RG_01M_2019_3035.geojson",
        crs="EPSG:3035",
    ).to_crs("EPSG:4326")

    fernwaermeatlas = pd.read_excel(
        snakemake.input.fernwaermeatlas,
        sheet_name="Fernwärmeatlas_öffentlich",
    )
    cities = gpd.read_file(snakemake.input.cities)
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    subnodes = gpd.read_file(snakemake.input.subnodes)

    add_subnodes(
        n,
        subnodes,
        cop=xr.open_dataarray(snakemake.input.cop_profiles),
        direct_heat_source_utilisation_profile=xr.open_dataarray(
            snakemake.input.direct_heat_source_utilisation_profiles
        ),
        head=snakemake.params.district_heating["add_subnodes"],
    )

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

    maybe_adjust_costs_and_potentials(
        n, snakemake.params["adjustments"], snakemake.wildcards.planning_horizons
    )
    n.export_to_netcdf(snakemake.output.network)
