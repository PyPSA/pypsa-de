import logging

logger = logging.getLogger(__name__)

import os
import sys

import geopandas as gpd
import pandas as pd
import pypsa
import xarray as xr
from typing import Dict, List

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.prepare_network import maybe_adjust_costs_and_potentials


def add_buses(n: pypsa.Network, subnode: pd.Series, name: str) -> None:
    """
    Add buses for a district heating subnode.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which buses will be added.
    subnode : pd.Series
        Series containing information about the district heating subnode.
    name : str
        Name prefix for the buses.

    Returns
    -------
    None
    """
    buses = (
        n.buses.filter(like=f"{subnode['cluster']} urban central", axis=0)
        .reset_index()
        .replace(
            {
                f"{subnode['cluster']} urban central": name,
                f"{subnode['cluster']}$": f"{subnode['cluster']} {subnode['Stadt']}",
            },
            regex=True,
        )
        .set_index("Bus")
    )
    n.add("Bus", buses.index, **buses)


def add_loads(
    n: pypsa.Network,
    n_copy: pypsa.Network,
    subnode: pd.Series,
    name: str,
    subnodes_head: gpd.GeoDataFrame,
) -> None:
    """
    Add loads for a district heating subnode and adjust mother node loads.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which loads will be added.
    n_copy : pypsa.Network
        Copy of the original network.
    subnode : pd.Series
        Series containing data about the district heating subnode to be added.
    name : str
        Name prefix for the loads.
    subnodes_head : gpd.GeoDataFrame
        GeoDataFrame containing data about largest district heating systems.

    Returns
    -------
    None
    """
    # Get heat loads for urban central heat and low-temperature heat for industry
    urban_central_heat_load_cluster = (
        n_copy.snapshot_weightings.generators
        @ n_copy.loads_t.p_set[f"{subnode['cluster']} urban central heat"]
    )
    low_temperature_heat_for_industry_load_cluster = (
        n_copy.loads.loc[
            f"{subnode['cluster']} low-temperature heat for industry", "p_set"
        ]
        * 8760
    )

    # Calculate share of low-temperature heat for industry in total district heating load of cluster
    dh_load_cluster = (
        urban_central_heat_load_cluster + low_temperature_heat_for_industry_load_cluster
    )

    dh_load_cluster_subnodes = subnodes_head.loc[
        subnodes_head.cluster == subnode["cluster"], "yearly_heat_demand_MWh"
    ].sum()
    lost_load = dh_load_cluster_subnodes - dh_load_cluster

    # District heating demand from Fernwärmeatlas exceeding the original cluster load is disregarded. The shares of the subsystems are set according to Fernwärmeatlas, while the aggregate load of cluster is preserved.
    if lost_load > 0:
        logger.warning(
            f"Aggregated district heating load of systems within {subnode['cluster']} exceeds load of cluster."
        )
        demand_ratio = subnode["yearly_heat_demand_MWh"] / dh_load_cluster_subnodes

        urban_central_heat_load = demand_ratio * n_copy.loads_t.p_set.filter(
            regex=f"{subnode['cluster']}.*urban central heat"
        ).sum(1).rename(f"{subnode['cluster']} {subnode['Stadt']} urban central heat")

        low_temperature_heat_for_industry_load = (
            demand_ratio
            * n_copy.loads.filter(
                regex=f"{subnode['cluster']}.*low-temperature heat for industry",
                axis=0,
            )["p_set"].sum()
        )
    else:
        # Calculate demand ratio between load of subnode according to Fernwärmeatlas and remaining load of assigned cluster
        demand_ratio = subnode["yearly_heat_demand_MWh"] / dh_load_cluster

        urban_central_heat_load = demand_ratio * n_copy.loads_t.p_set[
            f"{subnode['cluster']} urban central heat"
        ].rename(f"{subnode['cluster']} {subnode['Stadt']} urban central heat")

        low_temperature_heat_for_industry_load = (
            demand_ratio
            * n_copy.loads.loc[
                f"{subnode['cluster']} low-temperature heat for industry", "p_set"
            ]
        )

    # Add load components to subnode preserving the share of low-temperature heat for industry of the cluster
    n.add(
        "Load",
        f"{name} heat",
        bus=f"{name} heat",
        p_set=urban_central_heat_load,
        carrier="urban central heat",
    )

    n.add(
        "Load",
        f"{subnode['cluster']} {subnode['Stadt']} low-temperature heat for industry",
        bus=f"{name} heat",
        p_set=low_temperature_heat_for_industry_load,
        carrier="low-temperature heat for industry",
    )

    # Adjust loads of cluster buses
    n.loads_t.p_set.loc[
        :, f'{subnode["cluster"]} urban central heat'
    ] -= urban_central_heat_load

    n.loads.loc[
        f'{subnode["cluster"]} low-temperature heat for industry', "p_set"
    ] -= low_temperature_heat_for_industry_load

    if lost_load > 0:
        lost_load_subnode = subnode["yearly_heat_demand_MWh"] - (
            n.snapshot_weightings.generators @ urban_central_heat_load
            + low_temperature_heat_for_industry_load * 8760
        )
        logger.warning(
            f"District heating load of {subnode['cluster']} {subnode['Stadt']} is reduced by {lost_load_subnode} MWh/a."
        )


def add_stores(
    n: pypsa.Network,
    subnode: pd.Series,
    name: str,
    subnodes_rest: gpd.GeoDataFrame,
    dynamic_ptes_capacity: bool = False,
    limit_ptes_potential_mother_nodes: bool = False,
) -> None:
    """
    Add stores for a district heating subnode.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which stores will be added.
    subnode : pd.Series
        Series containing information about the district heating subnode.
    name : str
        Name prefix for the stores.
    subnodes_rest : gpd.GeoDataFrame
        GeoDataFrame containing information about remaining district heating subnodes.
    dynamic_ptes_capacity : bool, optional
        Whether to use dynamic PTES capacity, by default False
    limit_ptes_potential_mother_nodes : bool, optional
        Whether to limit PTES potential in mother nodes, by default False

    Returns
    -------
    None
    """
    # Replicate district heating stores of mother node for subnodes
    stores = (
        n.stores.filter(like=f"{subnode['cluster']} urban central", axis=0)
        .reset_index()
        .replace(
            {f"{subnode['cluster']} urban central": name},
            regex=True,
        )
        .set_index("Store")
    )

    # Restrict PTES capacity in subnodes
    stores.loc[stores.carrier.str.contains("pits$").index, "e_nom_max"] = subnode[
        "ptes_pot_mwh"
    ]

    if dynamic_ptes_capacity:
        e_max_pu_static = stores.e_max_pu
        e_max_pu = (
            n.stores_t.e_max_pu[f"{subnode['cluster']} urban central water pits"]
            .rename(f"{name} water pits")
            .to_frame()
            .reindex(columns=stores.index)
            .fillna(e_max_pu_static)
        )
        n.add(
            "Store",
            stores.index,
            e_max_pu=e_max_pu,
            **stores.drop("e_max_pu", axis=1),
        )
    else:
        n.add("Store", stores.index, **stores)

    # Limit storage potential in mother nodes
    if limit_ptes_potential_mother_nodes:
        mother_nodes_ptes_pot = subnodes_rest.groupby("cluster").ptes_pot_mwh.sum()

        mother_nodes_ptes_pot.index = (
            mother_nodes_ptes_pot.index + " urban central water pits"
        )
        n.stores.loc[mother_nodes_ptes_pot.index, "e_nom_max"] = mother_nodes_ptes_pot


def add_storage_units(n: pypsa.Network, subnode: pd.Series, name: str) -> None:
    """
    Add storage units for a district heating subnode.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which storage units will be added.
    subnode : pd.Series
        Series containing information about the district heating subnode.
    name : str
        Name prefix for the storage units.

    Returns
    -------
    None
    """
    # Replicate district heating storage units of mother node for subnodes
    storage_units = (
        n.storage_units.filter(like=f"{subnode['cluster']} urban central", axis=0)
        .reset_index()
        .replace(
            {f"{subnode['cluster']} urban central": name},
            regex=True,
        )
        .set_index("StorageUnit")
    )

    n.add("StorageUnit", storage_units.index, **storage_units)


def add_generators(n: pypsa.Network, subnode: pd.Series, name: str) -> None:
    """
    Add generators for a district heating subnode.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which generators will be added.
    subnode : pd.Series
        Series containing information about the district heating subnode.
    name : str
        Name prefix for the generators.

    Returns
    -------
    None
    """
    # Replicate district heating generators of mother node for subnodes
    generators = (
        n.generators.filter(like=f"{subnode['cluster']} urban central", axis=0)
        .reset_index()
        .replace(
            {f"{subnode['cluster']} urban central": name},
            regex=True,
        )
        .set_index("Generator")
    )
    n.add("Generator", generators.index, **generators)


def add_links(
    n: pypsa.Network,
    subnode: pd.Series,
    name: str,
    cop: xr.DataArray,
    direct_heat_source_utilisation_profile: xr.DataArray,
    heat_pump_sources: List[str],
    direct_utilisation_heat_sources: List[str],
    time_dep_hp_cop: bool,
    limited_heat_sources: List[str],
    heat_source_potentials: Dict[str, str],
) -> None:
    """
    Add links for a district heating subnode.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which links will be added.
    subnode : pd.Series
        Series containing information about the district heating subnode.
    name : str
        Name prefix for the links.
    cop : xr.DataArray
        COPs for heat pumps.
    direct_heat_source_utilisation_profile : xr.DataArray
        Direct heat source utilisation profiles.
    heat_pump_sources : List[str]
        List of heat pump sources.
    direct_utilisation_heat_sources : List[str]
        List of heat sources that can be directly utilized.
    time_dep_hp_cop : bool
        Whether to use time-dependent COPs for heat pumps.
    limited_heat_sources : List[str]
        List of heat sources with limited potential.
    heat_source_potentials : Dict[str, str]
        Dictionary mapping heat sources to paths with potential data.

    Returns
    -------
    None
    """
    # Replicate district heating links of mother node for subnodes with separate treatment for links with dynamic efficiencies
    links = (
        n.links.loc[~n.links.carrier.str.contains("heat pump|direct", regex=True)]
        .filter(like=f"{subnode['cluster']} urban central", axis=0)
        .reset_index()
        .replace(
            {f"{subnode['cluster']} urban central": name},
            regex=True,
        )
        .set_index("Link")
    )
    n.add("Link", links.index, **links)

    # Add heat pumps and direct heat source utilization to subnode
    for heat_source in heat_pump_sources:
        cop_heat_pump = (
            cop.sel(
                heat_system="urban central",
                heat_source=heat_source,
                name=f"{subnode['cluster']} {subnode['Stadt']}",
            )
            .to_pandas()
            .to_frame(name=f"{name} {heat_source} heat pump")
            .reindex(index=n.snapshots)
            if time_dep_hp_cop
            else n.links.filter(like=heat_source, axis=0).efficiency.mode()
        )

        heat_pump = (
            n.links.filter(
                regex=f"{subnode['cluster']} urban central.*{heat_source}.*heat pump",
                axis=0,
            )
            .reset_index()
            .replace(
                {f"{subnode['cluster']} urban central": name},
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

        if heat_source in direct_utilisation_heat_sources:
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
                    {f"{subnode['cluster']} urban central": name},
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
        if heat_source in limited_heat_sources:
            # get potential
            p_max_source = pd.read_csv(
                heat_source_potentials[heat_source],
                index_col=0,
            ).squeeze()[f"{subnode['cluster']} {subnode['Stadt']}"]
            # add potential to generator
            n.generators.loc[
                f"{subnode['cluster']} {subnode['Stadt']} urban central {heat_source} heat",
                "p_nom_max",
            ] = p_max_source


def add_subnodes(
    n: pypsa.Network,
    subnodes: gpd.GeoDataFrame,
    cop: xr.DataArray,
    direct_heat_source_utilisation_profile: xr.DataArray,
    head: int = 40,
    dynamic_ptes_capacity: bool = False,
    limit_ptes_potential_mother_nodes: bool = True,
    heat_pump_sources: List[str] = None,
    direct_utilisation_heat_sources: List[str] = None,
    time_dep_hp_cop: bool = False,
    limited_heat_sources: List[str] = None,
    heat_source_potentials: Dict[str, str] = None,
    output_path: str = None,
) -> None:
    """
    Add the largest district heating systems as subnodes to the network based on
    their heat demand. For each subnode, create individual district heating components
    from the corresponding mother node template, including heat sources, storage options,
    and heat pumps. Adjust loads of original mother nodes to maintain the overall
    energy balance.

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
    head : int
        Number of largest district heating systems to be added as subnodes.
    dynamic_ptes_capacity : bool
        Whether to use dynamic PTES capacity.
    limit_ptes_potential_mother_nodes : bool
        Whether to limit PTES potential in mother nodes.
    heat_pump_sources : List[str]
        List of heat pump sources.
    direct_utilisation_heat_sources : List[str]
        List of heat sources that can be directly utilized.
    time_dep_hp_cop : bool
        Whether to use time-dependent COPs for heat pumps.
    limited_heat_sources : List[str]
        List of heat sources with limited potential.
    heat_source_potentials : Dict[str, str]
        Dictionary mapping heat sources to paths with potential data.
    output_path : str
        Path to save the subnodes_head GeoDataFrame.
    Returns
    -------
    None
    """

    # Keep only n largest district heating networks according to head parameter
    subnodes_head = subnodes.sort_values(
        by="Wärmeeinspeisung in GWh/a", ascending=False
    ).head(head)

    if output_path:
        subnodes_head.to_file(output_path, driver="GeoJSON")

    subnodes_rest = subnodes[~subnodes.index.isin(subnodes_head.index)]

    n_copy = n.copy()

    # Add subnodes to network
    for _, subnode in subnodes_head.iterrows():
        name = f'{subnode["cluster"]} {subnode["Stadt"]} urban central'

        # Add different component types
        add_buses(n, subnode, name)
        add_loads(n, n_copy, subnode, name, subnodes_head)
        add_stores(
            n,
            subnode,
            name,
            subnodes_rest,
            dynamic_ptes_capacity,
            limit_ptes_potential_mother_nodes,
        )
        add_storage_units(n, subnode, name)
        add_generators(n, subnode, name)
        add_links(
            n,
            subnode,
            name,
            cop,
            direct_heat_source_utilisation_profile,
            heat_pump_sources,
            direct_utilisation_heat_sources,
            time_dep_hp_cop,
            limited_heat_sources,
            heat_source_potentials,
        )


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

    # Create a dictionary of heat source potentials for the limited heat sources
    heat_source_potentials = {}
    for source in snakemake.params.district_heating["limited_heat_sources"]:
        heat_source_potentials[source] = snakemake.input[source]

    add_subnodes(
        n,
        subnodes,
        cop=xr.open_dataarray(snakemake.input.cop_profiles),
        direct_heat_source_utilisation_profile=xr.open_dataarray(
            snakemake.input.direct_heat_source_utilisation_profiles
        ),
        head=snakemake.params.district_heating["subnodes"]["nlargest"],
        dynamic_ptes_capacity=snakemake.params.district_heating["ptes"][
            "dynamic_capacity"
        ],
        limit_ptes_potential_mother_nodes=snakemake.params.district_heating["subnodes"][
            "limit_ptes_potential"
        ]["limit_mother_nodes"],
        heat_pump_sources=snakemake.params.heat_pump_sources,
        direct_utilisation_heat_sources=snakemake.params.district_heating[
            "direct_utilisation_heat_sources"
        ],
        time_dep_hp_cop=snakemake.params.sector["time_dep_hp_cop"],
        limited_heat_sources=snakemake.params.district_heating["limited_heat_sources"],
        heat_source_potentials=heat_source_potentials,
        output_path=snakemake.output.district_heating_subnodes,
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
