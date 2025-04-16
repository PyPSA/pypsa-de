# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves optimal operation and capacity for a network with the option to
iteratively optimize while updating line reactances.

This script is used for optimizing the electrical network as well as the
sector coupled network.

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.

The optimization is based on the :func:`network.optimize` function.
Additionally, some extra constraints specified in :mod:`solve_network` are added.

.. note::

    The rules ``solve_elec_networks`` and ``solve_sector_networks`` run
    the workflow for all scenarios in the configuration file (``scenario:``)
    based on the rule :mod:`solve_network`.
"""

import importlib
import logging
import os
import pathlib
import re
import sys
from functools import partial
from typing import Any

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
import yaml
from _benchmark import memory_logger
from pypsa.descriptors import get_activity_mask
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.prepare_sector_network import get

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


class ObjectiveValueError(Exception):
    pass


def add_land_use_constraint_perfect(n: pypsa.Network) -> None:
    """
    Add global constraints for tech capacity limit.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance

    Returns
    -------
    pypsa.Network
        Network with added land use constraints
    """
    logger.info("Add land-use constraint for perfect foresight")

    def compress_series(s):
        def process_group(group):
            if group.nunique() == 1:
                return pd.Series(group.iloc[0], index=[None])
            else:
                return group

        return s.groupby(level=[0, 1]).apply(process_group)

    def new_index_name(t):
        # Convert all elements to string and filter out None values
        parts = [str(x) for x in t if x is not None]
        # Join with space, but use a dash for the last item if not None
        return " ".join(parts[:2]) + (f"-{parts[-1]}" if len(parts) > 2 else "")

    def check_p_min_p_max(p_nom_max):
        p_nom_min = n.generators[ext_i].groupby(grouper).sum().p_nom_min
        p_nom_min = p_nom_min.reindex(p_nom_max.index)
        check = (
            p_nom_min.groupby(level=[0, 1]).sum()
            > p_nom_max.groupby(level=[0, 1]).min()
        )
        if check.sum():
            logger.warning(
                f"summed p_min_pu values at node larger than technical potential {check[check].index}"
            )

    grouper = [n.generators.carrier, n.generators.bus, n.generators.build_year]
    ext_i = n.generators.p_nom_extendable
    # get technical limit per node and investment period
    p_nom_max = n.generators[ext_i].groupby(grouper).min().p_nom_max
    # drop carriers without tech limit
    p_nom_max = p_nom_max[~p_nom_max.isin([np.inf, np.nan])]
    # carrier
    carriers = p_nom_max.index.get_level_values(0).unique()
    gen_i = n.generators[(n.generators.carrier.isin(carriers)) & (ext_i)].index
    n.generators.loc[gen_i, "p_nom_min"] = 0
    # check minimum capacities
    check_p_min_p_max(p_nom_max)
    # drop multi entries in case p_nom_max stays constant in different periods
    # p_nom_max = compress_series(p_nom_max)
    # adjust name to fit syntax of nominal constraint per bus
    df = p_nom_max.reset_index()
    df["name"] = df.apply(
        lambda row: f"nom_max_{row['carrier']}"
        + (f"_{row['build_year']}" if row["build_year"] is not None else ""),
        axis=1,
    )

    for name in df.name.unique():
        df_carrier = df[df.name == name]
        bus = df_carrier.bus
        n.buses.loc[bus, name] = df_carrier.p_nom_max.values


def add_land_use_constraint(n: pypsa.Network, planning_horizons: str) -> None:
    """
    Add land use constraints for renewable energy potential.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    planning_horizons : str
        The planning horizon year as string

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with constraints added
    """
    # warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in [
        "solar",
        "solar rooftop",
        "solar-hsat",
        "onwind",
        "offwind-ac",
        "offwind-dc",
        "offwind-float",
    ]:
        ext_i = (n.generators.carrier == carrier) & ~n.generators.p_nom_extendable
        existing = (
            n.generators.loc[ext_i, "p_nom"]
            .groupby(n.generators.bus.map(n.buses.location))
            .sum()
        )
        existing.index += f" {carrier}-{planning_horizons}"
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    # check if existing capacities are larger than technical potential
    existing_large = n.generators[
        n.generators["p_nom_min"] > n.generators["p_nom_max"]
    ].index
    if len(existing_large):
        logger.warning(
            f"Existing capacities larger than technical potential for {existing_large},\
                        adjust technical potential to existing capacities"
        )
        n.generators.loc[existing_large, "p_nom_max"] = n.generators.loc[
            existing_large, "p_nom_min"
        ]

    n.generators["p_nom_max"] = n.generators["p_nom_max"].clip(lower=0)


def add_solar_potential_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add constraint to make sure the sum capacity of all solar technologies (fixed, tracking, ets. ) is below the region potential.

    Example:
    ES1 0: total solar potential is 10 GW, meaning:
           solar potential : 10 GW
           solar-hsat potential : 8 GW (solar with single axis tracking is assumed to have higher land use)
    The constraint ensures that:
           solar_p_nom + solar_hsat_p_nom * 1.13 <= 10 GW
    """
    land_use_factors = {
        "solar-hsat": config["renewable"]["solar"]["capacity_per_sqkm"]
        / config["renewable"]["solar-hsat"]["capacity_per_sqkm"],
    }
    rename = {"Generator-ext": "Generator"}

    solar_carriers = ["solar", "solar-hsat"]
    solar = n.generators[
        n.generators.carrier.isin(solar_carriers) & n.generators.p_nom_extendable
    ].index

    solar_today = n.generators[
        (n.generators.carrier == "solar") & (n.generators.p_nom_extendable)
    ].index
    solar_hsat = n.generators[(n.generators.carrier == "solar-hsat")].index

    if solar.empty:
        return

    land_use = pd.DataFrame(1, index=solar, columns=["land_use_factor"])
    for carrier, factor in land_use_factors.items():
        land_use = land_use.apply(
            lambda x: (x * factor) if carrier in x.name else x, axis=1
        )

    location = pd.Series(n.buses.index, index=n.buses.index)
    ggrouper = n.generators.loc[solar].bus
    rhs = (
        n.generators.loc[solar_today, "p_nom_max"]
        .groupby(n.generators.loc[solar_today].bus.map(location))
        .sum()
        - n.generators.loc[solar_hsat, "p_nom"]
        .groupby(n.generators.loc[solar_hsat].bus.map(location))
        .sum()
        * land_use_factors["solar-hsat"]
    ).clip(lower=0)

    lhs = (
        (n.model["Generator-p_nom"].rename(rename).loc[solar] * land_use.squeeze())
        .groupby(ggrouper)
        .sum()
    )

    logger.info("Adding solar potential constraint.")
    n.model.add_constraints(lhs <= rhs, name="solar_potential")


def add_co2_sequestration_limit(
    n: pypsa.Network,
    limit_dict: dict[str, float],
    planning_horizons: str | None,
) -> None:
    """
    Add a global constraint on the amount of Mt CO2 that can be sequestered.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    limit_dict : dict[str, float]
        CO2 sequestration potential limit constraints by year.
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    """

    if not n.investment_periods.empty:
        periods = n.investment_periods
        limit = pd.Series(
            {
                f"co2_sequestration_limit-{period}": limit_dict.get(period, 200)
                for period in periods
            }
        )
        names = limit.index
    else:
        limit = get(limit_dict, int(planning_horizons))
        periods = [np.nan]
        names = pd.Index(["co2_sequestration_limit"])

    n.add(
        "GlobalConstraint",
        names,
        sense=">=",
        constant=-limit * 1e6,
        type="operational_limit",
        carrier_attribute="co2 sequestered",
        investment_period=periods,
    )


def add_carbon_constraint(n: pypsa.Network, snapshots: pd.DatetimeIndex) -> None:
    glcs = n.global_constraints.query('type == "co2_atmosphere"')
    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last = n.snapshot_weightings.reset_index().groupby("period").last()
            last_i = last.set_index([last.index, last.timestep]).index
            final_e = n.model["Store-e"].loc[last_i, stores.index]
            time_valid = int(glc.loc["investment_period"])
            time_i = pd.IndexSlice[time_valid, :]
            lhs = final_e.loc[time_i, :] - final_e.shift(snapshot=1).loc[time_i, :]

            rhs = glc.constant
            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def add_carbon_budget_constraint(n: pypsa.Network, snapshots: pd.DatetimeIndex) -> None:
    glcs = n.global_constraints.query('type == "Co2Budget"')
    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last = n.snapshot_weightings.reset_index().groupby("period").last()
            last_i = last.set_index([last.index, last.timestep]).index
            final_e = n.model["Store-e"].loc[last_i, stores.index]
            time_valid = int(glc.loc["investment_period"])
            time_i = pd.IndexSlice[time_valid, :]
            weighting = n.investment_period_weightings.loc[time_valid, "years"]
            lhs = final_e.loc[time_i, :] * weighting

            rhs = glc.constant
            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def add_max_growth(n: pypsa.Network, opts: dict) -> None:
    """
    Add maximum growth rates for different carriers.
    """

    # take maximum yearly difference between investment periods since historic growth is per year
    factor = n.investment_period_weightings.years.max() * opts["factor"]
    for carrier in opts["max_growth"].keys():
        max_per_period = opts["max_growth"][carrier] * factor
        logger.info(
            f"set maximum growth rate per investment period of {carrier} to {max_per_period} GW."
        )
        n.carriers.loc[carrier, "max_growth"] = max_per_period * 1e3

    for carrier in opts["max_relative_growth"].keys():
        max_r_per_period = opts["max_relative_growth"][carrier]
        logger.info(
            f"set maximum relative growth per investment period of {carrier} to {max_r_per_period}."
        )
        n.carriers.loc[carrier, "max_relative_growth"] = max_r_per_period


def add_retrofit_gas_boiler_constraint(
    n: pypsa.Network, snapshots: pd.DatetimeIndex
) -> None:
    """
    Allow retrofitting of existing gas boilers to H2 boilers and impose load-following must-run condition on existing gas boilers.
    Modifies the network in place, no return value.

    n : pypsa.Network
        The PyPSA network to be modified
    snapshots : pd.DatetimeIndex
        The snapshots of the network
    """
    c = "Link"
    logger.info("Add constraint for retrofitting gas boilers to H2 boilers.")
    # existing gas boilers
    mask = n.links.carrier.str.contains("gas boiler") & ~n.links.p_nom_extendable
    gas_i = n.links[mask].index
    mask = n.links.carrier.str.contains("retrofitted H2 boiler")
    h2_i = n.links[mask].index

    n.links.loc[gas_i, "p_nom_extendable"] = True
    p_nom = n.links.loc[gas_i, "p_nom"]
    n.links.loc[gas_i, "p_nom"] = 0

    # heat profile
    cols = n.loads_t.p_set.columns[
        n.loads_t.p_set.columns.str.contains("heat")
        & ~n.loads_t.p_set.columns.str.contains("industry")
        & ~n.loads_t.p_set.columns.str.contains("agriculture")
    ]
    profile = n.loads_t.p_set[cols].div(
        n.loads_t.p_set[cols].groupby(level=0).max(), level=0
    )
    # to deal if max value is zero
    profile.fillna(0, inplace=True)
    profile.rename(columns=n.loads.bus.to_dict(), inplace=True)
    profile = profile.reindex(columns=n.links.loc[gas_i, "bus1"])
    profile.columns = gas_i

    rhs = profile.mul(p_nom)

    dispatch = n.model["Link-p"]
    active = get_activity_mask(n, c, snapshots, gas_i)
    rhs = rhs[active]
    p_gas = dispatch.sel(Link=gas_i)
    p_h2 = dispatch.sel(Link=h2_i)

    lhs = p_gas + p_h2

    n.model.add_constraints(lhs == rhs, name="gas_retrofit")


def prepare_network(
    n: pypsa.Network,
    solve_opts: dict,
    foresight: str,
    planning_horizons: str | None,
    co2_sequestration_potential: dict[str, float],
    limit_max_growth: dict[str, Any] | None = None,
) -> None:
    """
    Prepare network with various constraints and modifications.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    solve_opts : Dict
        Dictionary of solving options containing clip_p_max_pu, load_shedding etc.
    foresight : str
        Planning foresight type ('myopic' or 'perfect')
    planning_horizons : str or None
        The current planning horizon year or None for perfect foresight
    co2_sequestration_potential : Dict[str, float]
        CO2 sequestration potential constraints by year

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with added constraints
    """
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.links_t.p_max_pu,
            n.links_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    if load_shedding := solve_opts.get("load_shedding"):
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        # TODO: retrieve color and nice name from config
        n.add("Carrier", "load", color="#dd2e23", nice_name="Load shedding")
        buses_i = n.buses.index
        if not np.isscalar(load_shedding):
            # TODO: do not scale via sign attribute (use Eur/MWh instead of Eur/kWh)
            load_shedding = 1e2  # Eur/kWh

        n.add(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,  # Eur/kWh
            p_nom=1e9,  # kW
        )

    if solve_opts.get("curtailment_mode"):
        n.add("Carrier", "curtailment", color="#fedfed", nice_name="Curtailment")
        n.generators_t.p_min_pu = n.generators_t.p_max_pu
        buses_i = n.buses.query("carrier == 'AC'").index
        n.add(
            "Generator",
            buses_i,
            suffix=" curtailment",
            bus=buses_i,
            p_min_pu=-1,
            p_max_pu=0,
            marginal_cost=-0.1,
            carrier="curtailment",
            p_nom=1e6,
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if foresight == "myopic":
        add_land_use_constraint(n, planning_horizons)

    if foresight == "perfect":
        add_land_use_constraint_perfect(n)
        if limit_max_growth is not None and limit_max_growth["enable"]:
            add_max_growth(n, limit_max_growth)

    if n.stores.carrier.eq("co2 sequestered").any():
        limit_dict = co2_sequestration_potential
        add_co2_sequestration_limit(
            n, limit_dict=limit_dict, planning_horizons=planning_horizons
        )


def add_CCL_constraints(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> None:
    """
    Add CCL (country & carrier limit) constraint to the network.

    Add minimum and maximum levels of generator nominal capacity per carrier
    for individual countries. Opts and path for agg_p_nom_minmax.csv must be defined
    in config.yaml. Default file is available at data/agg_p_nom_minmax.csv.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : dict
        Configuration dictionary
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Example
    -------
    scenario:
        opts: [Co2L-CCL-24h]
    electricity:
        agg_p_nom_limits: data/agg_p_nom_minmax.csv
    """

    assert planning_horizons is not None, (
        "add_CCL_constraints are not implemented for perfect foresight, yet"
    )

    agg_p_nom_minmax = pd.read_csv(
        config["solving"]["agg_p_nom_limits"]["file"], index_col=[0, 1], header=[0, 1]
    )[planning_horizons]
    logger.info("Adding generation capacity constraints per carrier and country")
    p_nom = n.model["Generator-p_nom"]

    gens = n.generators.query("p_nom_extendable").rename_axis(index="Generator-ext")
    if config["solving"]["agg_p_nom_limits"]["agg_offwind"]:
        rename_offwind = {
            "offwind-ac": "offwind-all",
            "offwind-dc": "offwind-all",
            "offwind": "offwind-all",
        }
        gens = gens.replace(rename_offwind)
    grouper = pd.concat([gens.bus.map(n.buses.country), gens.carrier], axis=1)
    lhs = p_nom.groupby(grouper).sum().rename(bus="country")

    if config["solving"]["agg_p_nom_limits"]["include_existing"]:
        gens_cst = n.generators.query("~p_nom_extendable").rename_axis(
            index="Generator-cst"
        )
        gens_cst = gens_cst[
            (gens_cst["build_year"] + gens_cst["lifetime"]) >= int(planning_horizons)
        ]
        if config["solving"]["agg_p_nom_limits"]["agg_offwind"]:
            gens_cst = gens_cst.replace(rename_offwind)
        rhs_cst = (
            pd.concat(
                [gens_cst.bus.map(n.buses.country), gens_cst[["carrier", "p_nom"]]],
                axis=1,
            )
            .groupby(["bus", "carrier"])
            .sum()
        )
        rhs_cst.index = rhs_cst.index.rename({"bus": "country"})
        rhs_min = agg_p_nom_minmax["min"].dropna()
        idx_min = rhs_min.index.join(rhs_cst.index, how="left")
        rhs_min = rhs_min.reindex(idx_min).fillna(0)
        rhs = (rhs_min - rhs_cst.reindex(idx_min).fillna(0).p_nom).dropna()
        rhs[rhs < 0] = 0
        minimum = xr.DataArray(rhs).rename(dim_0="group")
    else:
        minimum = xr.DataArray(agg_p_nom_minmax["min"].dropna()).rename(dim_0="group")

    index = minimum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) >= minimum.loc[index], name="agg_p_nom_min"
        )

    if config["solving"]["agg_p_nom_limits"]["include_existing"]:
        rhs_max = agg_p_nom_minmax["max"].dropna()
        idx_max = rhs_max.index.join(rhs_cst.index, how="left")
        rhs_max = rhs_max.reindex(idx_max).fillna(0)
        rhs = (rhs_max - rhs_cst.reindex(idx_max).fillna(0).p_nom).dropna()
        rhs[rhs < 0] = 0
        maximum = xr.DataArray(rhs).rename(dim_0="group")
    else:
        maximum = xr.DataArray(agg_p_nom_minmax["max"].dropna()).rename(dim_0="group")

    index = maximum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) <= maximum.loc[index], name="agg_p_nom_max"
        )


def add_EQ_constraints(n, o, scaling=1e-1):
    """
    Add equity constraints to the network.

    Currently this is only implemented for the electricity sector only.

    Opts must be specified in the config.yaml.

    Parameters
    ----------
    n : pypsa.Network
    o : str

    Example
    -------
    scenario:
        opts: [Co2L-EQ0.7-24h]

    Require each country or node to on average produce a minimal share
    of its total electricity consumption itself. Example: EQ0.7c demands each country
    to produce on average at least 70% of its consumption; EQ0.7 demands
    each node to produce on average at least 70% of its consumption.
    """
    # TODO: Generalize to cover myopic and other sectors?
    float_regex = r"[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])
    if o[-1] == "c":
        ggrouper = n.generators.bus.map(n.buses.country)
        lgrouper = n.loads.bus.map(n.buses.country)
        sgrouper = n.storage_units.bus.map(n.buses.country)
    else:
        ggrouper = n.generators.bus
        lgrouper = n.loads.bus
        sgrouper = n.storage_units.bus
    load = (
        n.snapshot_weightings.generators
        @ n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    )
    inflow = (
        n.snapshot_weightings.stores
        @ n.storage_units_t.inflow.groupby(sgrouper, axis=1).sum()
    )
    inflow = inflow.reindex(load.index).fillna(0.0)
    rhs = scaling * (level * load - inflow)
    p = n.model["Generator-p"]
    lhs_gen = (
        (p * (n.snapshot_weightings.generators * scaling))
        .groupby(ggrouper.to_xarray())
        .sum()
        .sum("snapshot")
    )
    # TODO: double check that this is really needed, why do have to subtract the spillage
    if not n.storage_units_t.inflow.empty:
        spillage = n.model["StorageUnit-spill"]
        lhs_spill = (
            (spillage * (-n.snapshot_weightings.stores * scaling))
            .groupby(sgrouper.to_xarray())
            .sum()
            .sum("snapshot")
        )
        lhs = lhs_gen + lhs_spill
    else:
        lhs = lhs_gen
    n.model.add_constraints(lhs >= rhs, name="equity_min")


def add_BAU_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add business-as-usual (BAU) constraints for minimum capacities.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network instance
    config : dict
        Configuration dictionary containing BAU minimum capacities
    """
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    p_nom = n.model["Generator-p_nom"]
    ext_i = n.generators.query("p_nom_extendable")
    ext_carrier_i = xr.DataArray(ext_i.carrier.rename_axis("Generator-ext"))
    lhs = p_nom.groupby(ext_carrier_i).sum()
    rhs = mincaps[lhs.indexes["carrier"]].rename_axis("carrier")
    n.model.add_constraints(lhs >= rhs, name="bau_mincaps")


# TODO: think about removing or make per country
def add_SAFE_constraints(n, config):
    """
    Add a capacity reserve margin of a certain fraction above the peak demand.
    Renewable generators and storage do not contribute. Ignores network.

    Parameters
    ----------
        n : pypsa.Network
        config : dict

    Example
    -------
    config.yaml requires to specify opts:

    scenario:
        opts: [Co2L-SAFE-24h]
    electricity:
        SAFE_reservemargin: 0.1
    Which sets a reserve margin of 10% above the peak demand.
    """
    peakdemand = n.loads_t.p_set.sum(axis=1).max()
    margin = 1.0 + config["electricity"]["SAFE_reservemargin"]
    reserve_margin = peakdemand * margin
    conventional_carriers = config["electricity"]["conventional_carriers"]  # noqa: F841
    ext_gens_i = n.generators.query(
        "carrier in @conventional_carriers & p_nom_extendable"
    ).index
    p_nom = n.model["Generator-p_nom"].loc[ext_gens_i]
    lhs = p_nom.sum()
    exist_conv_caps = n.generators.query(
        "~p_nom_extendable & carrier in @conventional_carriers"
    ).p_nom.sum()
    rhs = reserve_margin - exist_conv_caps
    n.model.add_constraints(lhs >= rhs, name="safe_mintotalcap")


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.

    Parameters
    ----------
        n : pypsa.Network
        sns: pd.DatetimeIndex
        config : dict

    Example:
    --------
    config.yaml requires to specify operational_reserve:
    operational_reserve: # like https://genxproject.github.io/GenX/dev/core/#Reserves
        activate: true
        epsilon_load: 0.02 # percentage of load at each snapshot
        epsilon_vres: 0.02 # percentage of VRES at each snapshot
        contingency: 400000 # MW
    """
    reserve_config = config["electricity"]["operational_reserve"]
    EPSILON_LOAD = reserve_config["epsilon_load"]
    EPSILON_VRES = reserve_config["epsilon_vres"]
    CONTINGENCY = reserve_config["contingency"]

    # Reserve Variables
    n.model.add_variables(
        0, np.inf, coords=[sns, n.generators.index], name="Generator-r"
    )
    reserve = n.model["Generator-r"]
    summed_reserve = reserve.sum("Generator")

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        p_nom_vres = (
            n.model["Generator-p_nom"]
            .loc[vres_i.intersection(ext_i)]
            .rename({"Generator-ext": "Generator"})
        )
        lhs = summed_reserve + (
            p_nom_vres * (-EPSILON_VRES * xr.DataArray(capacity_factor))
        ).sum("Generator")

        # Total demand per t
        demand = get_as_dense(n, "Load", "p_set").sum(axis=1)

        # VRES potential of non extendable generators
        capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
        renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
        potential = (capacity_factor * renewable_capacity).sum(axis=1)

        # Right-hand-side
        rhs = EPSILON_LOAD * demand + EPSILON_VRES * potential + CONTINGENCY

        n.model.add_constraints(lhs >= rhs, name="reserve_margin")

    # additional constraint that capacity is not exceeded
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = n.model["Generator-p"]
    reserve = n.model["Generator-r"]

    capacity_variable = n.model["Generator-p_nom"].rename(
        {"Generator-ext": "Generator"}
    )
    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    lhs = dispatch + reserve - capacity_variable * xr.DataArray(p_max_pu[ext_i])

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)

    n.model.add_constraints(lhs <= rhs, name="Generator-p-reserve-upper")


def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.index.str.contains("battery discharger")
    charger_bool = n.links.index.str.contains("battery charger")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


def add_lossy_bidirectional_link_constraints(n):
    if not n.links.p_nom_extendable.any() or not any(n.links.get("reversed", [])):
        return

    carriers = n.links.loc[n.links.reversed, "carrier"].unique()  # noqa: F841
    backwards = n.links.query(
        "carrier in @carriers and p_nom_extendable and reversed"
    ).index
    forwards = backwards.str.replace("-reversed", "")
    lhs = n.model["Link-p_nom"].loc[backwards]
    rhs = n.model["Link-p_nom"].loc[forwards]
    n.model.add_constraints(lhs == rhs, name="Link-bidirectional_sync")


def add_chp_constraints(n):
    electric = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("electric")
    )
    heat = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("heat")
    )

    electric_ext = n.links[electric].query("p_nom_extendable").index
    heat_ext = n.links[heat].query("p_nom_extendable").index

    electric_fix = n.links[electric].query("~p_nom_extendable").index
    heat_fix = n.links[heat].query("~p_nom_extendable").index

    p = n.model["Link-p"]  # dimension: [time, link]

    # output ratio between heat and electricity and top_iso_fuel_line for extendable
    if not electric_ext.empty:
        p_nom = n.model["Link-p_nom"]

        lhs = (
            p_nom.loc[electric_ext]
            * (n.links.p_nom_ratio * n.links.efficiency)[electric_ext].values
            - p_nom.loc[heat_ext] * n.links.efficiency[heat_ext].values
        )
        n.model.add_constraints(lhs == 0, name="chplink-fix_p_nom_ratio")

        rename = {"Link-ext": "Link"}
        lhs = (
            p.loc[:, electric_ext]
            + p.loc[:, heat_ext]
            - p_nom.rename(rename).loc[electric_ext]
        )
        n.model.add_constraints(lhs <= 0, name="chplink-top_iso_fuel_line_ext")

    # top_iso_fuel_line for fixed
    if not electric_fix.empty:
        lhs = p.loc[:, electric_fix] + p.loc[:, heat_fix]
        rhs = n.links.p_nom[electric_fix]
        n.model.add_constraints(lhs <= rhs, name="chplink-top_iso_fuel_line_fix")

    # back-pressure
    if not electric.empty:
        lhs = (
            p.loc[:, heat] * (n.links.efficiency[heat] * n.links.c_b[electric].values)
            - p.loc[:, electric] * n.links.efficiency[electric]
        )
        n.model.add_constraints(lhs <= rhs, name="chplink-backpressure")


def add_pipe_retrofit_constraint(n):
    """
    Add constraint for retrofitting existing CH4 pipelines to H2 pipelines.
    """
    if "reversed" not in n.links.columns:
        n.links["reversed"] = False
    gas_pipes_i = n.links.query(
        "carrier == 'gas pipeline' and p_nom_extendable and ~reversed"
    ).index
    h2_retrofitted_i = n.links.query(
        "carrier == 'H2 pipeline retrofitted' and p_nom_extendable and ~reversed"
    ).index

    if h2_retrofitted_i.empty or gas_pipes_i.empty:
        return

    p_nom = n.model["Link-p_nom"]

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    lhs = p_nom.loc[gas_pipes_i] + CH4_per_H2 * p_nom.loc[h2_retrofitted_i]
    rhs = n.links.p_nom[gas_pipes_i].rename_axis("Link-ext")

    n.model.add_constraints(lhs == rhs, name="Link-pipe_retrofit")


def add_flexible_egs_constraint(n):
    """
    Upper bounds the charging capacity of the geothermal reservoir according to
    the well capacity.
    """
    well_index = n.links.loc[n.links.carrier == "geothermal heat"].index
    storage_index = n.storage_units.loc[
        n.storage_units.carrier == "geothermal heat"
    ].index

    p_nom_rhs = n.model["Link-p_nom"].loc[well_index]
    p_nom_lhs = n.model["StorageUnit-p_nom"].loc[storage_index]

    n.model.add_constraints(
        p_nom_lhs <= p_nom_rhs,
        name="upper_bound_charging_capacity_of_geothermal_reservoir",
    )


def add_co2_atmosphere_constraint(n, snapshots):
    glcs = n.global_constraints[n.global_constraints.type == "co2_atmosphere"]

    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last_i = snapshots[-1]
            lhs = n.model["Store-e"].loc[last_i, stores.index]
            rhs = glc.constant

            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")



def determine_planning_year(n):
    """
    Determine the planning year from the network using various methods.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network

    Returns
    -------
    int
        The identified planning year
    """
    import logging

    # Approach 1: Try to get year from wildcards
    planning_year = None
    try:
        planning_year = int(snakemake.wildcards.planning_horizons)
        logging.info(f"Using planning year {planning_year} from wildcards")
        return planning_year
    except (NameError, AttributeError, ValueError):
        logging.info("Could not get planning year from wildcards")

    # Approach 2: Check if network has investment_period attribute
    planning_year = getattr(n, 'investment_period', None)
    if planning_year is not None:
        logging.info(f"Using planning year {planning_year} from network investment_period attribute")
        return planning_year

    # Approach 3: Try to extract from snapshots
    if not n.snapshots.empty:
        try:
            first_snapshot = n.snapshots[0]
            if hasattr(first_snapshot, 'year'):
                planning_year = first_snapshot.year
                logging.info(f"Using planning year {planning_year} from first snapshot year attribute")
                return planning_year
            else:
                # Try to parse year from snapshot string
                try:
                    planning_year = int(str(first_snapshot)[:4])
                    logging.info(f"Using planning year {planning_year} from first snapshot string")
                    return planning_year
                except:
                    pass
        except:
            pass

    # Approach 4: Try to extract from filename if network has a 'filename' attribute
    filename = getattr(n, 'filename', '')
    if filename:
        import re
        year_match = re.search(r'_(\d{4})\.', filename)
        if year_match:
            planning_year = int(year_match.group(1))
            logging.info(f"Using planning year {planning_year} from filename")
            return planning_year

    # Approach 5: Default to 2030 if all else fails
    planning_year = 2030
    logging.warning(f"Could not determine planning year, using default: {planning_year}")
    return planning_year


def phase_out_conventional_generators(n, config=None):
    """
    Phase out conventional power plants in Germany based on technology-specific schedules.
    Each technology type can have its own phase-out timeline with custom capacity limits.
    All components are treated as links in the PyPSA network.
    """
    import logging

    logging.info("=== Applying technology-specific phase-out constraints ===")

    # Use provided config or get from network
    if config is None:
        config = getattr(n, 'config', {})
        if not config:
            try:
                config = snakemake.config
            except (NameError, AttributeError):
                logging.warning("Cannot find configuration, using default phase-out years")
                config = {}

    # Get current planning year
    planning_year = None
    try:
        planning_year = int(snakemake.wildcards.planning_horizons)
        logging.info(f"Using planning year {planning_year} from snakemake wildcards")
    except (NameError, AttributeError):
        # Try to extract year from component names
        for link in n.links.index:
            if '-20' in link:  # Look for patterns like "-2030" or "-2045"
                try:
                    year = int(link.split('-')[-1])
                    if 2000 <= year <= 2100:  # Sanity check
                        planning_year = year
                        logging.info(f"Extracted planning year {planning_year} from component names")
                        break
                except (ValueError, IndexError):
                    continue

    if not planning_year:
        # Try to get year from config
        planning_year = config.get("year", 2045)
        logging.info(f"Using planning year {planning_year} from config")

    # Get phase-out settings from config
    phaseout_settings = config.get("solving", {}).get("constraints", {}).get("generator_phaseout", {})

    # Check if the constraint is enabled
    if not phaseout_settings.get("apply_constraints", True):
        logging.info("Generator phase-out constraints are disabled in config")
        return n

    # Define technology types and their identifiers - all as links
    technology_patterns = {
        "CCGT": {
            "patterns": ["CCGT", "ccgt", "combined cycle gas"],
            "exclude": ["H2", "h2", "hydrogen"]  # Explicitly exclude H2 CCGT
        },
        "OCGT": {
            "patterns": ["OCGT", "ocgt", "open cycle gas"],
            "exclude": ["H2", "h2", "hydrogen"]  # Explicitly exclude H2 OCGT
        },
        "urban_central_gas_CHP": {
            "patterns": ["urban central gas", "gas CHP"],
            "exclude": ["CC", "cc", "H2", "h2", "hydrogen"]
        },
        "urban_central_gas_CHP_CC": {
            "patterns": ["urban central gas", "gas CHP"],
            "require": ["CC", "cc"],
            "exclude": ["H2", "h2", "hydrogen"]  # Also exclude hydrogen from CC
        },
        "urban_central_oil_CHP": {
            "patterns": ["urban central oil", "oil CHP"]
        },
        "waste_CHP": {
            "patterns": ["waste CHP", "waste chp"],
            "exclude": ["CC", "cc"]
        },
        "waste_CHP_CC": {
            "patterns": ["waste CHP", "waste chp"],
            "require": ["CC", "cc"]
        },
        "urban_central_H2_CHP": {
            "patterns": ["urban central H2", "hydrogen CHP", "H2 CHP"]
        },
        "H2_OCGT": {
            "patterns": ["H2 OCGT", "hydrogen OCGT"],
            "exclude": ["retrofit"]
        },
        "H2_retrofit_OCGT": {
            "patterns": ["H2 retrofit", "hydrogen retrofit"]
        },
        "urban_central_solid_biomass_CHP": {
            "patterns": ["urban central solid biomass", "solid biomass CHP", "biomass CHP"],
            "exclude": ["CC", "cc"]
        },
        "urban_central_solid_biomass_CHP_CC": {
            "patterns": ["urban central solid biomass", "solid biomass CHP", "biomass CHP"],
            "require": ["CC", "cc"]
        }
    }

    # Process each technology type
    for tech_name, tech_info in technology_patterns.items():
        # Skip if tech doesn't have phase-out settings
        tech_settings = phaseout_settings.get(tech_name, {})
        if not tech_settings or tech_settings == {}:
            logging.info(f"No phase-out settings found for {tech_name}, skipping")
            continue

        # Get capacity limits and phase-out year for this technology
        capacity_limits = tech_settings.get("capacity_limits", {})
        phase_out_year = tech_settings.get("phase_out_year", 2045)

        # Determine current capacity limit
        current_limit = 1.0  # Default: no reduction

        # If the year is at or past complete phase-out, set limit to 0
        if planning_year >= phase_out_year:
            current_limit = 0.0
        # Otherwise, check for specific year limits
        else:
            for year_str, limit in sorted(capacity_limits.items()):
                year = int(year_str)
                if planning_year >= year:
                    current_limit = limit  # Use the most recent limit that applies

        # If no reduction needed, skip this technology
        if current_limit >= 1.0:
            logging.info(f"No capacity reduction needed for {tech_name} in year {planning_year}")
            continue

        logging.info(f"Applying {tech_name} phase-out for year {planning_year} with capacity limit {current_limit:.1%}")

        # Identify components for this technology
        patterns = tech_info.get("patterns", [])
        require_patterns = tech_info.get("require", [])
        exclude_patterns = tech_info.get("exclude", [])

        # Filter components based on patterns
        component_indices = []

        for idx in n.links.index:
            # Must start with DE (Germany)
            if not idx.startswith('DE'):
                continue

            # Must match one of the patterns
            if not any(pattern.lower() in idx.lower() or
                      (pattern.lower() in str(n.links.at[idx, 'carrier']).lower()
                       if 'carrier' in n.links.columns else False)
                      for pattern in patterns):
                continue

            # Must include all required patterns
            if require_patterns and not all(
                    pattern.lower() in idx.lower() or
                    (pattern.lower() in str(n.links.at[idx, 'carrier']).lower()
                     if 'carrier' in n.links.columns else False)
                    for pattern in require_patterns):
                continue

            # Must exclude all exclude patterns
            if exclude_patterns and any(
                    pattern.lower() in idx.lower() or
                    (pattern.lower() in str(n.links.at[idx, 'carrier']).lower()
                     if 'carrier' in n.links.columns else False)
                    for pattern in exclude_patterns):
                continue

            component_indices.append(idx)

        if not component_indices:
            logging.info(f"No {tech_name} components found in Germany")
            continue

        logging.info(f"Found {len(component_indices)} {tech_name} components to adjust")

        # Apply phase-out by setting capacity to zero or reduced capacity
        for comp_idx in component_indices:
            original_p_nom = n.links.at[comp_idx, 'p_nom']

            # Calculate new capacity based on limit
            new_p_nom = original_p_nom * current_limit
            n.links.at[comp_idx, 'p_nom'] = new_p_nom

            # Also adjust p_nom_max if it exists
            if 'p_nom_max' in n.links.columns:
                n.links.at[comp_idx, 'p_nom_max'] = new_p_nom

            # If extendable and complete phase-out, set p_nom_min to zero
            if current_limit == 0 and 'p_nom_min' in n.links.columns and 'p_nom_extendable' in n.links.columns and n.links.at[
                comp_idx, 'p_nom_extendable']:
                n.links.at[comp_idx, 'p_nom_min'] = 0

            # If complete phase-out, set to not extendable to prevent the optimizer from adding capacity
            if current_limit == 0 and 'p_nom_extendable' in n.links.columns:
                n.links.at[comp_idx, 'p_nom_extendable'] = False

            logging.info(
                f"Reduced {comp_idx} capacity from {original_p_nom} to {new_p_nom} ({current_limit:.1%} of original)")

    logging.info("Successfully applied technology-specific phase-out constraints")
    return n


def add_resistive_heater_constraints(n):
    """
    Imposes max capacity constraints on resistive heaters in Germany.
    Uses the same approach as the working heat_pump_constraints function.
    """
    import logging
    import pandas as pd

    logging.info("[Resistive heater constraints] Starting resistive heater constraint application...")

    # Get config and check if constraints are enabled
    if not hasattr(n, 'config'):
        try:
            n.config = snakemake.config
            logging.info("[Resistive heater constraints] Got config from snakemake")
        except (NameError, AttributeError):
            logging.warning("[Resistive heater constraints] No config found, skipping resistive heater constraints")
            return

    # Try multiple possible paths for resistive heater config
    resistive_config = None
    possible_paths = [
        # Path 1: Under solving.constraints
        n.config.get('solving', {}).get('constraints', {}).get('resistive_heaters', {}),
        # Path 2: Directly under constraints
        n.config.get('constraints', {}).get('resistive_heaters', {}),
        # Path 3: At top level
        n.config.get('resistive_heaters', {})
    ]

    # Log the config structure for debugging
    logging.info("[Resistive heater constraints] Config contains these top-level keys: " +
                 str(list(n.config.keys())))
    if 'solving' in n.config:
        logging.info("[Resistive heater constraints] 'solving' section contains: " +
                     str(list(n.config['solving'].keys())))
        if 'constraints' in n.config['solving']:
            logging.info("[Resistive heater constraints] 'solving.constraints' section contains: " +
                         str(list(n.config['solving']['constraints'].keys())))

    # Try each path and use the first one that has apply_constraints=True
    for i, path_config in enumerate(possible_paths):
        if path_config and path_config.get('apply_constraints', False):
            resistive_config = path_config
            logging.info(f"[Resistive heater constraints] Found enabled resistive heater config at path {i + 1}")
            break

    # If none of the paths had apply_constraints=True
    if resistive_config is None or not resistive_config.get('apply_constraints', False):
        logging.info("[Resistive heater constraints] Resistive heater constraints disabled in config, skipping")
        return

    # Get planning year with multiple fallbacks
    try:
        planning_year = int(snakemake.wildcards.planning_horizons)
        logging.info(f"[Resistive heater constraints] Got planning year {planning_year} from snakemake wildcards")
    except:
        # Try other methods to determine planning year
        planning_year = 2045  # Default fallback
        logging.info(f"[Resistive heater constraints] Using default planning year {planning_year}")

    year_str = str(planning_year)
    logging.info(f"[Resistive heater constraints] Final planning year: {year_str}")

    # Get active scenario
    active_scenario = resistive_config.get('active_scenario', 'restrictive')
    logging.info(f"[Resistive heater constraints] Using {active_scenario} resistive heater scenario from config")

    # Get countries to apply constraints to
    countries_to_process = resistive_config.get('countries', [])
    if not countries_to_process:
        # If empty, apply to all countries in the network
        countries_to_process = n.buses.country.unique().tolist()
        logging.info(f"[Resistive heater constraints] Applying to all countries in network: {countries_to_process}")
    else:
        logging.info(f"[Resistive heater constraints] Applying to specified countries: {countries_to_process}")

    # Get scenario config
    scenario_config = resistive_config.get('scenarios', {}).get(active_scenario, {})

    # Debug log the scenario config
    logging.info(
        f"[Resistive heater constraints] Available years in {active_scenario} scenario: {list(scenario_config.keys())}")
    logging.info(
        f"[Resistive heater constraints] Types of year keys: {[type(k).__name__ for k in scenario_config.keys()]}")

    # Find all resistive heater carriers in the network
    all_carriers = n.links.carrier.unique()
    resistive_heater_carriers = [c for c in all_carriers if 'resistive' in c.lower()]
    logging.info(f"[Resistive heater constraints] Detected resistive heater carriers: {resistive_heater_carriers}")

    # Ensure all carriers are in the carriers component
    if "carriers" not in n.components.keys():
        n.add("Carrier", resistive_heater_carriers)
        logging.info(f"[Resistive heater constraints] Added resistive heater carriers to network")
    elif not all(c in n.carriers.index for c in resistive_heater_carriers):
        # Add any missing carriers
        missing_carriers = [c for c in resistive_heater_carriers if c not in n.carriers.index]
        for c in missing_carriers:
            n.add("Carrier", c)
        logging.info(f"[Resistive heater constraints] Added missing carriers: {missing_carriers}")

    # Log the initial capacities
    logging.info("[Resistive heater constraints] Initial capacities before constraints:")
    for rh_carrier in resistive_heater_carriers:
        total_cap = n.links.loc[n.links.carrier == rh_carrier, 'p_nom'].sum()
        num_links = len(n.links[n.links.carrier == rh_carrier])
        logging.info(f"  {rh_carrier}: {total_cap:.2f} MW across {num_links} links")

    # Process each country
    countries_processed = 0
    total_links_constrained = 0

    for country in countries_to_process:
        # Get year config
        year_dict = None
        for key in scenario_config:
            if str(key) == year_str:
                year_dict = scenario_config[key]
                logging.info(f"[Resistive heater constraints] Found year {key} matching {year_str}")
                break

        if year_dict is None:
            logging.warning(f"[Resistive heater constraints] No data for year {year_str} in {active_scenario} scenario")
            continue

        # Get buses for this country
        country_buses = n.buses.index[n.buses.country == country]
        if len(country_buses) == 0:
            logging.warning(f"[Resistive heater constraints] No buses found for {country}")
            continue

        logging.info(f"[Resistive heater constraints] Processing country: {country} with {len(country_buses)} buses")

        # Track country totals
        country_max = 0
        country_links = 0

        # Find all resistive heater links in this country
        for rh_carrier, caps in year_dict.items():
            if rh_carrier in resistive_heater_carriers:
                # Get the maximum capacity constraint
                max_cap = caps.get("max", float('inf'))

                # Look for links with bus0 (input) in this country
                relevant_links = n.links.index[
                    (n.links.carrier == rh_carrier) &
                    (n.links.bus0.isin(country_buses))
                    ]

                if len(relevant_links) == 0:
                    logging.warning(f"[Resistive heater constraints] No {rh_carrier} links found in {country}")
                    continue

                # Calculate current total capacity
                current_capacity = n.links.loc[relevant_links, 'p_nom'].sum()

                # Log what we found
                logging.info(
                    f"[Resistive heater constraints] Found {len(relevant_links)} {rh_carrier} links in {country} "
                    f"with initial capacity {current_capacity:.2f} MW, target: {max_cap} MW")

                # Set capacity equally across all links
                capacity_per_link = max_cap / len(relevant_links)

                # Apply the constraint to each link
                for link in relevant_links:
                    # CRITICAL: Following the same pattern as the heat_pump_constraints function
                    # Set fixed capacity and make non-extendable
                    original_capacity = n.links.at[link, 'p_nom']
                    n.links.at[link, 'p_nom'] = capacity_per_link
                    n.links.at[link, 'p_nom_extendable'] = False

                    logging.info(
                        f"[Resistive heater constraints] RESET: {link} capacity from {original_capacity:.2f} MW to {capacity_per_link:.2f} MW")

                logging.info(f"[Resistive heater constraints] Limited {rh_carrier} capacity to {max_cap} MW "
                             f"({capacity_per_link:.2f} MW per link across {len(relevant_links)} links)")

                country_max += max_cap
                country_links += len(relevant_links)
                total_links_constrained += len(relevant_links)

            elif rh_carrier not in resistive_heater_carriers and rh_carrier != 'total':
                logging.warning(f"[Resistive heater constraints] Carrier {rh_carrier} not found in network")

        # Process 'total' constraint if present
        if 'total' in year_dict and resistive_heater_carriers:
            total_max = year_dict['total'].get('max', float('inf'))

            # Get all resistive heater links in this country
            all_resistive_links = []
            for carrier in resistive_heater_carriers:
                carrier_links = n.links.index[
                    (n.links.carrier == carrier) &
                    (n.links.bus0.isin(country_buses))
                    ]
                all_resistive_links.extend(carrier_links)

            # Remove duplicates
            all_resistive_links = list(set(all_resistive_links))

            if all_resistive_links:
                # Calculate current total capacity
                current_total = n.links.loc[all_resistive_links, 'p_nom'].sum()

                if current_total > total_max:
                    # Scale down all links proportionally
                    scaling_factor = total_max / current_total if current_total > 0 else 0

                    for link in all_resistive_links:
                        original_capacity = n.links.at[link, 'p_nom']
                        new_capacity = original_capacity * scaling_factor
                        n.links.at[link, 'p_nom'] = new_capacity
                        n.links.at[link, 'p_nom_extendable'] = False

                        logging.info(
                            f"[Resistive heater constraints] TOTAL ADJUST: {link} capacity from {original_capacity:.2f} MW to {new_capacity:.2f} MW")

                logging.info(
                    f"[Resistive heater constraints] Applied total resistive heater constraint: {total_max} MW")

        # Log country summary
        if country_links > 0:
            logging.info(
                f"[Resistive heater constraints] {country}: Applied capacity constraints to {country_links} links, max total: {country_max} MW")
            countries_processed += 1

    # Log overall summary
    if total_links_constrained > 0:
        logging.info(
            f"[Resistive heater constraints] Successfully applied constraints to {total_links_constrained} resistive heater links across {countries_processed} countries")
    else:
        logging.warning(f"[Resistive heater constraints] No resistive heater constraints were applied")

    # Add extra log to show the post-constraint capacities
    try:
        resistive_capacities = []
        for idx, row in n.links.iterrows():
            if 'resistive' in row['carrier'].lower():
                if row['bus0'].startswith(tuple(countries_to_process)):
                    resistive_capacities.append({
                        'Link': idx,
                        'Carrier': row['carrier'],
                        'Country': row['bus0'][:2],
                        'Max Capacity': row.get('p_nom_max', 'Not set'),
                        'Current Capacity': row.get('p_nom', 0),
                        'Extendable': row.get('p_nom_extendable', False)
                    })

        if resistive_capacities:
            capacities_df = pd.DataFrame(resistive_capacities)
            # Group by carrier and country
            capacities_summary = capacities_df.groupby(['Country', 'Carrier']).agg({
                'Current Capacity': 'sum',
                'Link': 'count',
                'Max Capacity': lambda x: sum([y for y in x if y != 'Not set'])
            }).reset_index()
            capacities_summary.columns = ['Country', 'Carrier', 'Total Capacity (MW)', 'Number of Links',
                                          'Max Capacity (MW)']
            logging.info(
                f"[Resistive heater constraints] Current resistive heater capacities by carrier:\n{capacities_summary.to_string()}")

            # Add overall summary by country
            country_summary = capacities_df.groupby('Country').agg({
                'Current Capacity': 'sum',
                'Link': 'count',
                'Max Capacity': lambda x: sum([y for y in x if y != 'Not set'])
            }).reset_index()
            country_summary.columns = ['Country', 'Total Capacity (MW)', 'Number of Links', 'Max Capacity (MW)']
            logging.info(
                f"[Resistive heater constraints] Current resistive heater capacities by country:\n{country_summary.to_string()}")

            # Verify final capacities
            total_capacity = capacities_df['Current Capacity'].sum()
            if total_capacity > country_max:
                logging.error(f"[Resistive heater constraints] WARNING: Final capacity ({total_capacity:.2f} MW) "
                              f"still exceeds limit ({country_max} MW)!")
            else:
                logging.info(f"[Resistive heater constraints] SUCCESS: Final capacity ({total_capacity:.2f} MW) "
                             f"now respects limit ({country_max} MW)")
        else:
            logging.info("[Resistive heater constraints] No resistive heater links with installed capacity found")
    except Exception as e:
        logging.error(f"[Resistive heater constraints] Error generating capacity summary: {e}")

    logging.info("[Resistive heater constraints] Completed resistive heater constraint application")


def add_heat_pump_constraints(n):
    """
    Imposes min/max capacity constraints on heat pumps in Germany.
    Follows the exact same pattern as the working resistive_heater_constraints function.
    """
    import logging
    import pandas as pd
    import sys
    import traceback

    try:
        # Add a very visible message that will appear even in minimal logging
        print("=============== HEAT PUMP CONSTRAINTS RUNNING ===============", file=sys.stderr)
        logging.info("=============== HEAT PUMP CONSTRAINTS RUNNING ===============")

        # Get config and check if constraints are enabled
        if not hasattr(n, 'config'):
            try:
                n.config = snakemake.config
                logging.info("[HP constraints] Got config from snakemake")
            except (NameError, AttributeError):
                logging.warning("[HP constraints] No config found, skipping heat pump constraints")
                print("[HP constraints] No config found, skipping", file=sys.stderr)
                return

        # Try multiple possible paths for heat pump config
        hp_config = None
        possible_paths = [
            # Path 1: Under solving.constraints
            n.config.get('solving', {}).get('constraints', {}).get('heat_pumps', {}),
            # Path 2: Directly under constraints
            n.config.get('constraints', {}).get('heat_pumps', {}),
            # Path 3: At top level
            n.config.get('heat_pumps', {})
        ]

        # Log the config structure for debugging
        logging.info("[HP constraints] Config contains these top-level keys: " +
                     str(list(n.config.keys())))
        if 'solving' in n.config:
            logging.info("[HP constraints] 'solving' section contains: " +
                         str(list(n.config['solving'].keys())))
            if 'constraints' in n.config['solving']:
                logging.info("[HP constraints] 'solving.constraints' section contains: " +
                             str(list(n.config['solving']['constraints'].keys())))

        # Try each path and use the first one that has apply_constraints=True
        for i, path_config in enumerate(possible_paths):
            if path_config and path_config.get('apply_constraints', False):
                hp_config = path_config
                logging.info(f"[HP constraints] Found enabled heat pump config at path {i + 1}")
                break

        # If none of the paths had apply_constraints=True
        if hp_config is None or not hp_config.get('apply_constraints', False):
            logging.info("[HP constraints] Heat pump constraints disabled in config, skipping")
            print("[HP constraints] Disabled in config", file=sys.stderr)
            return

        # Configuration verification
        logging.info("[HP constraints] Heat pump constraints are ENABLED in config")
        logging.info(f"[HP constraints] Full HP config: {hp_config}")

        # Verify scenario data exists
        scenario_name = hp_config.get('active_scenario', 'conservative')
        scenarios = hp_config.get('scenarios', {})

        if scenario_name not in scenarios:
            logging.warning(f"[HP constraints] ERROR: Scenario '{scenario_name}' not found in config!")
            for s in scenarios:
                logging.warning(f"  - Available scenario: {s}")
            print(f"[HP constraints] Scenario '{scenario_name}' not found!", file=sys.stderr)
            return

        scenario_data = scenarios[scenario_name]
        logging.info(f"[HP constraints] Using scenario: {scenario_name} with years: {list(scenario_data.keys())}")

        # Get planning year with multiple fallbacks
        planning_year = None

        # Try from snakemake wildcards
        try:
            planning_year = int(snakemake.wildcards.planning_horizons)
            logging.info(f"[HP constraints] Got planning year {planning_year} from snakemake wildcards")
        except Exception as e:
            logging.info(f"[HP constraints] Error getting year from wildcards: {e}")

            # Try from network attributes
            try:
                if hasattr(n, 'investment_periods') and len(n.investment_periods) > 0:
                    planning_year = int(max(n.investment_periods))
                    logging.info(f"[HP constraints] Got year {planning_year} from investment periods")
                elif hasattr(n, 'snapshot_weightings'):
                    # Try to extract year from snapshot names
                    years = [int(str(s).split('-')[0]) for s in n.snapshots[:5] if '-' in str(s)]
                    if years:
                        planning_year = max(years)
                        logging.info(f"[HP constraints] Got year {planning_year} from snapshots")
            except Exception as e2:
                logging.warning(f"[HP constraints] Error getting year from network: {e2}")

            # Last resort: try to extract from filename or use default
            if planning_year is None:
                try:
                    # Check if any attribute has year information
                    for attr in dir(n):
                        if 'year' in attr.lower() and hasattr(n, attr):
                            val = getattr(n, attr)
                            if isinstance(val, (int, float)) and 1900 < val < 2100:
                                planning_year = int(val)
                                logging.info(f"[HP constraints] Got year {planning_year} from {attr}")
                                break
                except:
                    pass

                # Ultimate fallback
                if planning_year is None:
                    planning_year = 2045
                    logging.warning(f"[HP constraints] Using default year {planning_year}")

        year_str = str(planning_year)
        logging.info(f"[HP constraints] Final planning year: {year_str}")
        print(f"[HP constraints] Using planning year: {year_str}", file=sys.stderr)

        # Get active scenario
        active_scenario = hp_config.get('active_scenario', 'conservative')
        logging.info(f"[HP constraints] Using {active_scenario} heat pump scenario from config")

        # Get countries to apply constraints to
        countries_to_process = hp_config.get('countries', [])
        if not countries_to_process:
            # If empty, apply to all countries in the network
            countries_to_process = n.buses.country.unique().tolist()
            logging.info(f"[HP constraints] Applying to all countries in network: {countries_to_process}")
        else:
            logging.info(f"[HP constraints] Applying to specified countries: {countries_to_process}")

        # Get scenario config
        scenario_config = hp_config.get('scenarios', {}).get(active_scenario, {})

        # Debug log the scenario config
        logging.info(
            f"[HP constraints] Available years in {active_scenario} scenario: {list(scenario_config.keys())}")
        logging.info(
            f"[HP constraints] Types of year keys: {[type(k).__name__ for k in scenario_config.keys()]}")

        # Find all heat pump carriers in the network
        all_carriers = n.links.carrier.unique()

        # DETAILED CARRIER ANALYSIS
        logging.info("[HP constraints] DETAILED CARRIER ANALYSIS:")
        for carrier in all_carriers:
            count = len(n.links[n.links.carrier == carrier])
            bus_examples = n.links[n.links.carrier == carrier].bus1.iloc[:3].tolist() if count > 0 else []
            is_heat_related = any(term in str(carrier).lower() for term in ['heat', 'thermal', 'hp', 'pump'])

            if is_heat_related or count > 0:
                logging.info(f"  Carrier '{carrier}': {count} links, example buses: {bus_examples}")
                if count > 0:
                    # Check a sample link
                    sample_link = n.links[n.links.carrier == carrier].index[0]
                    logging.info(f"    Sample link '{sample_link}': p_nom={n.links.at[sample_link, 'p_nom']}, "
                                 f"p_nom_extendable={n.links.at[sample_link, 'p_nom_extendable']}")

        # Find heat pump carriers using multiple methods
        heat_pump_carriers = []

        # Method 1: Exact "heat pump" match
        exact_matches = [c for c in all_carriers if 'heat pump' in str(c).lower()]
        if exact_matches:
            heat_pump_carriers.extend(exact_matches)
            logging.info(f"[HP constraints] Found exact 'heat pump' matches: {exact_matches}")

        # Method 2: Contains both "heat" and "pump"
        if not heat_pump_carriers:
            word_matches = [c for c in all_carriers if
                            ('heat' in str(c).lower() and 'pump' in str(c).lower())]
            if word_matches:
                heat_pump_carriers.extend(word_matches)
                logging.info(f"[HP constraints] Found 'heat'+'pump' matches: {word_matches}")

        # Method 3: Contains "hp" abbreviation
        if not heat_pump_carriers:
            abbrev_matches = [c for c in all_carriers if 'hp' in str(c).lower().split()]
            if abbrev_matches:
                heat_pump_carriers.extend(abbrev_matches)
                logging.info(f"[HP constraints] Found 'hp' abbreviation matches: {abbrev_matches}")

        # Remove duplicates
        heat_pump_carriers = list(set(heat_pump_carriers))

        if not heat_pump_carriers:
            print("[HP constraints] NO HEAT PUMP CARRIERS FOUND IN NETWORK!", file=sys.stderr)
            logging.warning("[HP constraints] ⚠️ NO HEAT PUMP CARRIERS FOUND IN NETWORK!")

            # Try one last approach - use the keys from our config as potential carriers
            year_dict = None
            for key in scenario_config:
                if str(key) == year_str:
                    year_dict = scenario_config[key]
                    break

            if year_dict:
                potential_carriers = list(year_dict.keys())
                logging.info(f"[HP constraints] Using carriers from config as fallback: {potential_carriers}")

                # Check if any of these exist in the network with slight variations
                for config_carrier in potential_carriers:
                    if config_carrier == 'total':
                        continue

                    # Try different variations
                    variations = [
                        config_carrier,
                        config_carrier.lower(),
                        config_carrier.replace(' ', '_'),
                        config_carrier.replace(' ', '-'),
                        config_carrier.replace('heat pump', 'heat_pump'),
                        config_carrier.replace('heat pump', 'HP')
                    ]

                    for var in variations:
                        matches = [c for c in all_carriers if var in c]
                        if matches:
                            logging.info(f"[HP constraints] Found carrier '{config_carrier}' as '{matches}'")
                            heat_pump_carriers.extend(matches)
                            break

                # Remove duplicates again
                heat_pump_carriers = list(set(heat_pump_carriers))

                if heat_pump_carriers:
                    logging.info(f"[HP constraints] Using these carriers: {heat_pump_carriers}")
                else:
                    logging.warning("[HP constraints] Could not find any matching carriers - cannot apply constraints")
                    print("[HP constraints] No matching carriers found", file=sys.stderr)
                    return

        # Log the initial capacities
        logging.info("[HP constraints] Initial capacities before constraints:")
        for hp_carrier in heat_pump_carriers:
            total_cap = n.links.loc[n.links.carrier == hp_carrier, 'p_nom'].sum()
            num_links = len(n.links[n.links.carrier == hp_carrier])
            logging.info(f"  {hp_carrier}: {total_cap:.2f} MW across {num_links} links")

        # Process each country
        countries_processed = 0
        total_links_constrained = 0

        for country in countries_to_process:
            # Get year config
            year_dict = None
            for key in scenario_config:
                if str(key) == year_str:
                    year_dict = scenario_config[key]
                    logging.info(f"[HP constraints] Found year {key} matching {year_str}")
                    break

            if year_dict is None:
                logging.warning(f"[HP constraints] No data for year {year_str} in {active_scenario} scenario")
                continue

            # Get buses for this country
            country_buses = n.buses.index[n.buses.country == country]
            if len(country_buses) == 0:
                logging.warning(f"[HP constraints] No buses found for {country}")
                continue

            logging.info(f"[HP constraints] Processing country: {country} with {len(country_buses)} buses")

            # Track country totals
            country_min = 0
            country_max = 0
            country_links = 0

            # Find all heat pump links in this country
            for hp_carrier_key, caps in year_dict.items():
                if hp_carrier_key == 'total':
                    continue

                # Check if this carrier exists in our network
                matching_carriers = [c for c in heat_pump_carriers if hp_carrier_key.lower() in c.lower()]

                if not matching_carriers:
                    logging.warning(
                        f"[HP constraints] Warning: Constraint carrier '{hp_carrier_key}' not found in network")
                    # Try partial match
                    matching_carriers = [c for c in heat_pump_carriers if
                                         any(part in c.lower() for part in hp_carrier_key.lower().split())]
                    if matching_carriers:
                        logging.info(f"[HP constraints] Using partial matches: {matching_carriers}")
                    else:
                        logging.warning(
                            f"[HP constraints] No matches for '{hp_carrier_key}' - skipping this constraint")
                        continue

                for hp_carrier in matching_carriers:
                    # Get the minimum and maximum capacity constraints
                    min_cap = caps.get("min", 0)
                    max_cap = caps.get("max", float('inf'))
                    fixed_capacity = min_cap == max_cap

                    # Try both bus0 and bus1 to find links
                    links_by_bus1 = n.links.index[
                        (n.links.carrier == hp_carrier) &
                        (n.links.bus1.isin(country_buses))
                        ]

                    links_by_bus0 = n.links.index[
                        (n.links.carrier == hp_carrier) &
                        (n.links.bus0.isin(country_buses))
                        ]

                    # Decide which set to use
                    if len(links_by_bus1) > 0:
                        relevant_links = links_by_bus1
                        bus_type = "bus1 (output)"
                    elif len(links_by_bus0) > 0:
                        relevant_links = links_by_bus0
                        bus_type = "bus0 (input)"
                    else:
                        logging.warning(f"[HP constraints] No {hp_carrier} links found in {country}")
                        continue

                    # Calculate current total capacity
                    current_capacity = n.links.loc[relevant_links, 'p_nom'].sum()

                    # Log what we found
                    logging.info(
                        f"[HP constraints] Found {len(relevant_links)} {hp_carrier} links in {country} using {bus_type} "
                        f"with initial capacity {current_capacity:.2f} MW, target: {min_cap}-{max_cap} MW")

                    # Distribute capacity equally across all links
                    capacity_per_link = min_cap / len(relevant_links)
                    max_per_link = max_cap / len(relevant_links)

                    # Apply the constraint to each link
                    for link in relevant_links:
                        original_cap = n.links.at[link, 'p_nom']
                        original_ext = n.links.at[link, 'p_nom_extendable']

                        if fixed_capacity:
                            # Set fixed capacity and make non-extendable
                            n.links.at[link, 'p_nom'] = capacity_per_link
                            n.links.at[link, 'p_nom_extendable'] = False
                            logging.info(
                                f"[HP constraints] FIXED {link}: {original_cap:.2f} → {capacity_per_link:.2f} MW (fixed)")
                        else:
                            # Set min/max bounds and keep extendable
                            n.links.at[link, 'p_nom_extendable'] = True
                            n.links.at[link, 'p_nom_min'] = capacity_per_link
                            n.links.at[link, 'p_nom_max'] = max_per_link
                            logging.info(
                                f"[HP constraints] BOUNDED {link}: set range {capacity_per_link:.2f}-{max_per_link:.2f} MW")

                    if fixed_capacity:
                        logging.info(f"[HP constraints] Fixed {hp_carrier} capacity to exactly {min_cap} MW "
                                     f"({capacity_per_link:.2f} MW per link across {len(relevant_links)} links)")
                    else:
                        logging.info(f"[HP constraints] Set {hp_carrier} capacity range to {min_cap}-{max_cap} MW "
                                     f"({capacity_per_link:.2f}-{max_per_link:.2f} MW per link)")

                    country_min += min_cap
                    country_max += max_cap
                    country_links += len(relevant_links)
                    total_links_constrained += len(relevant_links)

            # Process 'total' constraint if present
            if 'total' in year_dict and heat_pump_carriers:
                total_caps = year_dict['total']
                total_min = total_caps.get('min', 0)
                total_max = total_caps.get('max', float('inf'))

                logging.info(f"[HP constraints] Processing 'total' constraint: min={total_min}, max={total_max}")

                # Find all heat pump links in this country
                all_hp_links = []
                for carrier in heat_pump_carriers:
                    # Try with both bus0 and bus1
                    links_by_bus1 = n.links.index[
                        (n.links.carrier == carrier) &
                        (n.links.bus1.isin(country_buses))
                        ]

                    links_by_bus0 = n.links.index[
                        (n.links.carrier == carrier) &
                        (n.links.bus0.isin(country_buses))
                        ]

                    all_hp_links.extend(links_by_bus1)
                    all_hp_links.extend(links_by_bus0)

                # Remove duplicates
                all_hp_links = list(set(all_hp_links))

                if all_hp_links:
                    # Calculate current total capacity
                    current_total = n.links.loc[all_hp_links, 'p_nom'].sum()
                    logging.info(
                        f"[HP constraints] Total heat pump capacity: {current_total:.2f} MW across {len(all_hp_links)} links")

                    # Check if we need to adjust
                    if current_total < total_min:
                        # Need to increase capacity
                        increase_factor = total_min / current_total if current_total > 0 else 0
                        logging.info(f"[HP constraints] Need to increase capacity by factor {increase_factor:.2f}")

                        for link in all_hp_links:
                            original_capacity = n.links.at[link, 'p_nom']
                            new_capacity = original_capacity * increase_factor if increase_factor > 0 else total_min / len(
                                all_hp_links)
                            n.links.at[link, 'p_nom'] = new_capacity
                            n.links.at[link, 'p_nom_extendable'] = False
                            logging.info(
                                f"[HP constraints] INCREASED {link}: {original_capacity:.2f} → {new_capacity:.2f} MW")

                    elif current_total > total_max:
                        # Need to decrease capacity
                        decrease_factor = total_max / current_total
                        logging.info(f"[HP constraints] Need to decrease capacity by factor {decrease_factor:.2f}")

                        for link in all_hp_links:
                            original_capacity = n.links.at[link, 'p_nom']
                            new_capacity = original_capacity * decrease_factor
                            n.links.at[link, 'p_nom'] = new_capacity
                            n.links.at[link, 'p_nom_extendable'] = False
                            logging.info(
                                f"[HP constraints] DECREASED {link}: {original_capacity:.2f} → {new_capacity:.2f} MW")

                    logging.info(f"[HP constraints] Applied total heat pump constraint: {total_min}-{total_max} MW")

            # Log country summary
            if country_links > 0:
                if country_min == country_max:
                    logging.info(
                        f"[HP constraints] {country}: Applied fixed capacity to {country_links} links, total: {country_min} MW")
                else:
                    logging.info(
                        f"[HP constraints] {country}: Applied constraints to {country_links} links, range: {country_min}-{country_max} MW")
                countries_processed += 1

        # Log overall summary
        if total_links_constrained > 0:
            logging.info(
                f"[HP constraints] Successfully applied constraints to {total_links_constrained} heat pump links across {countries_processed} countries")
            print(f"[HP constraints] Applied constraints to {total_links_constrained} links", file=sys.stderr)
        else:
            logging.warning(f"[HP constraints] ⚠️ NO HEAT PUMP CONSTRAINTS WERE APPLIED!")
            print("[HP constraints] NO CONSTRAINTS WERE APPLIED!", file=sys.stderr)

            # Extra debug information about what's in the network
            logging.warning("[HP constraints] DEBUG: Let's look at all carrier types in the network:")
            for carrier in n.links.carrier.unique():
                count = len(n.links[n.links.carrier == carrier])
                logging.warning(f"[HP constraints] Carrier '{carrier}': {count} links")

        # Add a final success message that will be visible even with minimal logging
        print("=============== HEAT PUMP CONSTRAINTS COMPLETED ===============", file=sys.stderr)
        logging.info("=============== HEAT PUMP CONSTRAINTS COMPLETED ===============")

        # Final verification
        try:
            # Check final capacities
            logging.info("[HP constraints] VERIFICATION - Final heat pump capacities:")
            for hp_carrier in heat_pump_carriers:
                country_capacities = {}
                for country in countries_to_process:
                    country_buses = n.buses.index[n.buses.country == country]

                    # Try both bus0 and bus1
                    links_by_bus1 = n.links.index[
                        (n.links.carrier == hp_carrier) &
                        (n.links.bus1.isin(country_buses))
                        ]

                    links_by_bus0 = n.links.index[
                        (n.links.carrier == hp_carrier) &
                        (n.links.bus0.isin(country_buses))
                        ]

                    all_links = list(set(list(links_by_bus1) + list(links_by_bus0)))

                    if all_links:
                        capacity = n.links.loc[all_links, 'p_nom'].sum()
                        country_capacities[country] = capacity

                if country_capacities:
                    logging.info(f"  {hp_carrier}: {country_capacities}")
                else:
                    logging.info(f"  {hp_carrier}: No links found")
        except Exception as e:
            logging.error(f"[HP constraints] Error in final verification: {e}")

        return total_links_constrained > 0

    except Exception as e:
        error_msg = f"ERROR IN HEAT PUMP CONSTRAINTS: {e}"
        print(error_msg, file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        logging.error(error_msg)
        logging.error(traceback.format_exc())
        return False


# Add heat pumps with enhanced flexibility logging
def add_heat_pump_flexibility_constraints(n):
    """
    Add constraints to limit heat pump operation during peak hours.
    This is compatible with older PyPSA versions that don't have pypsa.linopt.
    """
    import logging
    import pandas as pd
    import numpy as np

    # Get config settings
    config = getattr(n, 'config', {})
    sector_config = config.get('sector', {})
    flexibility_enabled = sector_config.get("heat_pump_flexibility", False)

    if not flexibility_enabled:
        logging.info("[HP-FLEX] Heat pump flexibility is disabled in config")
        return

    hp_avail_max = sector_config.get("heat_pump_avail_max", 1.0)
    hp_avail_min = sector_config.get("heat_pump_avail_min", 0.3)
    peak_hours = sector_config.get("heat_pump_peak_hours", [7, 8, 17, 18, 19])

    logging.info(
        f"[HP-FLEX] Applying heat pump flexibility: max={hp_avail_max}, min={hp_avail_min}, peak_hours={peak_hours}")

    # Find all heat pumps
    heat_pumps = n.links.index[n.links.carrier.str.contains("heat pump", case=False)]

    if len(heat_pumps) == 0:
        logging.info("[HP-FLEX] No heat pumps found")
        return

    logging.info(f"[HP-FLEX] Found {len(heat_pumps)} heat pumps to apply flexibility to")

    # Create peak hour masks for weekdays
    snapshots = n.snapshots
    peak_hour_mask = pd.Series(False, index=snapshots)
    for t in snapshots:
        if t.weekday() < 5 and t.hour in peak_hours:  # Weekday and peak hour
            peak_hour_mask.loc[t] = True

    peak_hours_count = peak_hour_mask.sum()
    logging.info(f"[HP-FLEX] Identified {peak_hours_count} peak hours out of {len(snapshots)} total snapshots")

    # For older PyPSA versions, we use p_max_pu
    # Ensure p_max_pu exists in links_t
    if not hasattr(n, 'links_t'):
        n.links_t = pd.DataFrame()

    if 'p_max_pu' not in n.links_t:
        n.links_t['p_max_pu'] = pd.DataFrame(index=n.snapshots)

    # Apply the flexibility settings through p_max_pu
    for hp in heat_pumps:
        if hp not in n.links_t.p_max_pu.columns:
            # Initialize with max availability
            n.links_t.p_max_pu[hp] = hp_avail_max

        # Apply peak hour restrictions
        n.links_t.p_max_pu.loc[peak_hour_mask, hp] = hp_avail_min

    logging.info(f"[HP-FLEX] Applied p_max_pu constraints to {len(heat_pumps)} heat pumps")

    # Make sure PyPSA knows to use p_max_pu
    # This might be needed for some versions
    try:
        n.consistency_check()
    except:
        logging.warning("[HP-FLEX] Consistency check failed, but continuing anyway")


def verify_heat_pump_flexibility_settings(n):
    """
    Check that heat pump flexibility settings are properly applied.
    Version for older PyPSA that uses p_max_pu instead of direct constraints.
    """
    import logging
    import sys

    print("=============== VERIFYING HEAT PUMP FLEXIBILITY SETTINGS ===============", file=sys.stderr)
    logging.warning("=============== VERIFYING HEAT PUMP FLEXIBILITY SETTINGS ===============")

    # Get all heat pumps
    heat_pumps = n.links.index[n.links.carrier.str.contains("heat pump", case=False)]
    if len(heat_pumps) == 0:
        logging.warning("[HP-VERIFY] No heat pumps found in network")
        return

    logging.warning(f"[HP-VERIFY] Found {len(heat_pumps)} heat pumps in network")

    # Check config
    config = getattr(n, 'config', {})
    sector_config = config.get('sector', {})
    flexibility_enabled = sector_config.get("heat_pump_flexibility", False)
    hp_avail_max = sector_config.get("heat_pump_avail_max", 1.0)
    hp_avail_min = sector_config.get("heat_pump_avail_min", 0.3)
    peak_hours = sector_config.get("heat_pump_peak_hours", [7, 8, 17, 18, 19])

    logging.warning(f"[HP-VERIFY] Flexibility enabled: {flexibility_enabled}")
    logging.warning(f"[HP-VERIFY] Max availability: {hp_avail_max}")
    logging.warning(f"[HP-VERIFY] Min availability (during peaks): {hp_avail_min}")
    logging.warning(f"[HP-VERIFY] Peak hours: {peak_hours}")

    # Check if p_max_pu exists and has the right values
    if hasattr(n, 'links_t') and 'p_max_pu' in n.links_t:
        # Check a sample heat pump
        if heat_pumps.size > 0 and heat_pumps[0] in n.links_t.p_max_pu.columns:
            sample_hp = heat_pumps[0]
            p_max_values = n.links_t.p_max_pu[sample_hp]

            # Create peak hour masks
            snapshots = n.snapshots
            peak_hour_mask = pd.Series(False, index=snapshots)
            for t in snapshots:
                if t.weekday() < 5 and t.hour in peak_hours:
                    peak_hour_mask.loc[t] = True

            # Check values
            peak_values = p_max_values[peak_hour_mask]
            non_peak_values = p_max_values[~peak_hour_mask]

            logging.warning(f"[HP-VERIFY] Sample heat pump: {sample_hp}")
            logging.warning(f"[HP-VERIFY] Peak hours p_max_pu: {peak_values.unique()}")
            logging.warning(f"[HP-VERIFY] Non-peak hours p_max_pu: {non_peak_values.unique()}")

            # Verify values match the config
            if all(abs(v - hp_avail_min) < 0.01 for v in peak_values.unique()):
                logging.warning("[HP-VERIFY] Peak hour values match configuration ✓")
            else:
                logging.warning("[HP-VERIFY] WARNING: Peak hour values do not match configuration!")

            if all(abs(v - hp_avail_max) < 0.01 for v in non_peak_values.unique()):
                logging.warning("[HP-VERIFY] Non-peak hour values match configuration ✓")
            else:
                logging.warning("[HP-VERIFY] WARNING: Non-peak hour values do not match configuration!")
        else:
            logging.warning("[HP-VERIFY] Could not find heat pumps in p_max_pu")
    else:
        logging.warning("[HP-VERIFY] p_max_pu not found in network! Flexibility will not be applied.")

    print("=============== HEAT PUMP FLEXIBILITY VERIFICATION COMPLETE ===============", file=sys.stderr)
    logging.warning("=============== HEAT PUMP FLEXIBILITY VERIFICATION COMPLETE ===============")

def confirm_p_max_pu_usage(n):
    """Check if p_max_pu is actually used in the optimization."""
    import logging

    # Check if there are any constraints related to p_max_pu
    if hasattr(n, 'model'):
        # Try to find constraints that use p_max_pu
        p_max_constraints = [c for c in n.model.constraints
                             if str(c).find('p_max_pu') >= 0]
        logging.warning(f"Found {len(p_max_constraints)} constraints using p_max_pu")
        if p_max_constraints:
            # Log a sample constraint
            logging.warning(f"Sample constraint: {p_max_constraints[0]}")
    else:
        logging.warning("No model attribute found - can't verify constraint usage")


def analyze_heat_pump_operation(n):
    """Analyze how heat pumps are being used after optimization."""
    import logging
    import numpy as np
    import pandas as pd

    logging.warning("Analyzing heat pump operation patterns...")

    # Get all heat pumps
    heat_pumps = n.links.index[n.links.carrier.str.contains("heat pump", case=False)]
    if not heat_pumps.empty:
        # Get operation data
        if 'p' in n.links_t:
            # Extract operation for peak and non-peak hours
            # Creating masks for peak hours
            config = getattr(n, 'config', {})
            sector_config = config.get('sector', {})
            peak_hours = sector_config.get("heat_pump_peak_hours", [7, 8, 17, 18, 19])

            weekday_mask = np.array([t.weekday() < 5 for t in n.snapshots])
            hour_mask = np.array([t.hour in peak_hours for t in n.snapshots])
            peak_mask = weekday_mask & hour_mask

            peak_hours_df = n.links_t.p.loc[peak_mask, heat_pumps]
            non_peak_hours_df = n.links_t.p.loc[~peak_mask, heat_pumps]

            # Calculate average utilization
            peak_avg = peak_hours_df.mean().mean()
            non_peak_avg = non_peak_hours_df.mean().mean()

            logging.warning(f"Heat pump average usage during peak hours: {peak_avg:.2f} MW")
            logging.warning(f"Heat pump average usage during non-peak hours: {non_peak_avg:.2f} MW")
            logging.warning(f"Usage ratio (peak/non-peak): {peak_avg / non_peak_avg if non_peak_avg else 0:.2f}")

            # Check capacity factors
            for hp in heat_pumps[:5]:  # Check first 5 pumps
                cap = n.links.at[hp, 'p_nom']
                if cap > 0:
                    peak_cf = peak_hours_df[hp].mean() / cap
                    non_peak_cf = non_peak_hours_df[hp].mean() / cap
                    logging.warning(f"{hp}: Peak CF = {peak_cf:.2f}, Non-peak CF = {non_peak_cf:.2f}")
        else:
            logging.warning("No operation data (links_t.p) found in results")
    else:
        logging.warning("No heat pumps found")


def add_storage_capacity_limits(n):
    """
    Adds constraints to enforce minimum/maximum battery storage capacity for Germany (DE),
    with support for different battery types.

    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network to which constraints will be applied
    """
    import logging
    import json

    logging.info("[Storage capacity limits] Starting storage capacity limit constraint application...")

    # Get config
    config = n.config
    storage_config = config.get('solving', {}).get('constraints', {}).get('storage_capacity_limits', {})

    if not storage_config:
        # Try to find storage_capacity_limits at root level
        storage_config = config.get('storage_capacity_limits', {})

    if not storage_config:
        # Try to find it in constraints at root level
        storage_config = config.get('constraints', {}).get('storage_capacity_limits', {})

    # Log the full config structure for debugging
    try:
        logging.info(f"[Storage capacity limits] Full storage config: {json.dumps(storage_config, indent=2)}")
    except:
        logging.info(f"[Storage capacity limits] Storage config: {storage_config}")

    if not storage_config.get('apply_constraints', False):
        logging.info("[Storage capacity limits] Storage capacity limit constraints disabled in config")
        return

    # Get planning year
    try:
        planning_year = int(snakemake.wildcards.planning_horizons)
        logging.info(f"[Storage capacity limits] Using planning year {planning_year} from wildcards")
    except (NameError, AttributeError, ValueError):
        # Fallback to default
        planning_year = 2045
        logging.warning(
            f"[Storage capacity limits] Could not get planning year from wildcards, using default: {planning_year}")

    year_str = str(planning_year)

    # Focus on Germany (DE)
    country = 'DE'

    # Check if we're using the new battery_types structure
    battery_types = storage_config.get('battery_types', {})

    if battery_types:
        # Using new structure with multiple battery types
        for battery_type, type_config in battery_types.items():
            logging.info(f"[Storage capacity limits] Processing battery type: {battery_type}")

            # Get carriers for this battery type
            battery_carriers = type_config.get('carriers', [])
            if not battery_carriers:
                logging.warning(f"[Storage capacity limits] No carriers defined for battery type: {battery_type}")
                continue

            # Get limits for this year and battery type
            year_limits = type_config.get('limits', {}).get(year_str, {})

            # If not found as string, try with integer
            if not year_limits and planning_year in type_config.get('limits', {}):
                year_limits = type_config.get('limits', {}).get(planning_year, {})
                logging.info(f"[Storage capacity limits] Found limits using integer key {planning_year}")

            # If still not found, try nearest year
            if not year_limits:
                available_years = [int(y) for y in type_config.get('limits', {}).keys() if str(y).isdigit()]
                if available_years:
                    nearest_year = min(available_years, key=lambda x: abs(x - planning_year))
                    year_limits = type_config.get('limits', {}).get(str(nearest_year), {})
                    logging.info(
                        f"[Storage capacity limits] Using nearest year {nearest_year} instead of {planning_year}")

            if not year_limits:
                logging.warning(
                    f"[Storage capacity limits] No limits defined for year {planning_year} and battery type {battery_type}")
                continue

            if country not in year_limits:
                logging.warning(
                    f"[Storage capacity limits] No limits defined for {country} in battery type {battery_type}")
                continue

            country_limits = year_limits[country]
            min_capacity = country_limits.get('min', None)  # In GW
            max_capacity = country_limits.get('max', None)  # In GW

            if min_capacity is not None:
                min_capacity = min_capacity * 1000  # Convert GW to MW
            if max_capacity is not None:
                max_capacity = max_capacity * 1000  # Convert GW to MW

            logging.info(
                f"[Storage capacity limits] Limits for {country}, {battery_type}: min={min_capacity} MW, max={max_capacity} MW")

            if min_capacity is None and max_capacity is None:
                logging.warning(f"[Storage capacity limits] No min/max capacity defined for {country}, {battery_type}")
                continue

            # Apply constraints for this battery type
            apply_battery_type_constraints(n, country, battery_carriers, min_capacity, max_capacity, battery_type)

    else:
        # Fallback to old structure with single battery type
        logging.info("[Storage capacity limits] Using legacy config structure with single battery type")

        # Get limits for this year
        year_limits = storage_config.get('limits', {}).get(year_str, {})

        # If not found as string, try with integer
        if not year_limits and planning_year in storage_config.get('limits', {}):
            year_limits = storage_config.get('limits', {}).get(planning_year, {})
            logging.info(f"[Storage capacity limits] Found limits using integer key {planning_year}")

        # If still not found, try nearest year
        if not year_limits:
            available_years = [int(y) for y in storage_config.get('limits', {}).keys() if str(y).isdigit()]
            if available_years:
                nearest_year = min(available_years, key=lambda x: abs(x - planning_year))
                year_limits = storage_config.get('limits', {}).get(str(nearest_year), {})
                logging.info(f"[Storage capacity limits] Using nearest year {nearest_year} instead of {planning_year}")

        if not year_limits:
            logging.warning(f"[Storage capacity limits] No limits defined for year {planning_year}")
            return

        if country not in year_limits:
            logging.warning(f"[Storage capacity limits] No limits defined for {country}")
            return

        country_limits = year_limits[country]
        min_capacity = country_limits.get('min', None)  # In GW
        max_capacity = country_limits.get('max', None)  # In GW

        if min_capacity is not None:
            min_capacity = min_capacity * 1000  # Convert GW to MW
        if max_capacity is not None:
            max_capacity = max_capacity * 1000  # Convert GW to MW

        logging.info(f"[Storage capacity limits] Limits for {country}: min={min_capacity} MW, max={max_capacity} MW")

        if min_capacity is None and max_capacity is None:
            logging.warning(f"[Storage capacity limits] No min/max capacity defined for {country}")
            return

        # Get list of battery carriers from config or use defaults
        default_battery_carriers = ['battery', 'home battery']
        battery_carriers = storage_config.get('battery_carriers', default_battery_carriers)

        # Apply constraints for the legacy single battery type
        apply_battery_type_constraints(n, country, battery_carriers, min_capacity, max_capacity, "all_batteries")

    logging.info("[Storage capacity limits] Completed storage capacity limit constraint application")


def apply_battery_type_constraints(n, country, battery_carriers, min_capacity, max_capacity, battery_type_name):
    """
    Helper function to apply constraints for a specific battery type
    """
    import logging

    logging.info(f"[Storage capacity limits] Applying constraints for battery type: {battery_type_name}")
    logging.info(f"[Storage capacity limits] Battery carriers: {battery_carriers}")

    # Expand to include variations (like 'battery charger', 'battery discharger', etc.)
    expanded_battery_carriers = []
    for carrier in battery_carriers:
        expanded_battery_carriers.append(carrier)
        expanded_battery_carriers.append(f"{carrier} charger")
        expanded_battery_carriers.append(f"{carrier} discharger")

    logging.info(f"[Storage capacity limits] Battery carriers (expanded): {expanded_battery_carriers}")

    try:
        # Find all links in Germany that are battery-related
        # First find buses in Germany
        country_buses = n.buses.index[n.buses.country == country]
        if len(country_buses) == 0:
            logging.warning(f"[Storage capacity limits] No buses found for {country}, skipping")
            return

        # Find battery links connected to these buses
        battery_links = []
        for link in n.links.index:
            link_carrier = n.links.at[link, 'carrier']
            bus0 = n.links.at[link, 'bus0']
            bus1 = n.links.at[link, 'bus1']

            # Check if it's a battery link (more flexible matching)
            is_battery = any(battery in link_carrier.lower() for battery in battery_carriers)
            is_in_germany = (bus0 in country_buses) or (bus1 in country_buses)

            if is_battery and is_in_germany:
                battery_links.append(link)

        # Log all found battery links for debugging
        logging.info(f"[Storage capacity limits] Found {len(battery_links)} {battery_type_name} links in {country}")
        for i, link in enumerate(battery_links):
            if i < 10:  # Limit to first the 10 to avoid too much logging
                logging.info(
                    f"[Storage capacity limits] Battery link {i + 1}: {link}, carrier: {n.links.at[link, 'carrier']}")

        if len(battery_links) == 0:
            logging.warning(f"[Storage capacity limits] No {battery_type_name} links found for {country}, skipping")
            return

        # Group battery links by type
        chargers = []
        dischargers = []
        for link in battery_links:
            carrier = n.links.at[link, 'carrier'].lower()
            if 'charger' in carrier:
                chargers.append(link)
            elif 'discharger' in carrier:
                dischargers.append(link)
            else:
                # For other battery links without explicit charger/discharger in name
                # Check efficiency or bus connection pattern
                efficiency = n.links.at[link, 'efficiency']
                bus0 = n.links.at[link, 'bus0']
                bus1 = n.links.at[link, 'bus1']

                if 'battery' in bus1.lower():
                    # If the destination bus has 'battery' in the name, it's likely a charger
                    chargers.append(link)
                elif 'battery' in bus0.lower():
                    # If the source bus has 'battery' in the name, it's likely a discharger
                    dischargers.append(link)
                elif efficiency > 0:
                    # Fallback - based on efficiency
                    chargers.append(link)
                else:
                    dischargers.append(link)

        logging.info(
            f"[Storage capacity limits] Found {len(chargers)} chargers and {len(dischargers)} dischargers for {battery_type_name}")

        # Check if we have any battery links
        if len(chargers) + len(dischargers) == 0:
            logging.warning(
                f"[Storage capacity limits] No {battery_type_name} chargers/dischargers found for {country}, skipping")
            return

        # Since we're constraining capacity, we need to be careful about double-counting
        # We'll focus on dischargers since they represent the usable capacity
        target_links = dischargers if dischargers else chargers

        # To avoid double counting when both chargers and dischargers exist, we'll use just one type
        if len(dischargers) > 0 and len(chargers) > 0:
            target_links = dischargers
            logging.info(f"[Storage capacity limits] Using only dischargers to avoid double-counting")

        # For extendable links, we need to constrain the p_nom_opt variable
        extendable_links = [link for link in target_links if n.links.at[link, 'p_nom_extendable']]

        # Also check for already fixed capacity in non-extendable links
        fixed_links = [link for link in target_links if not n.links.at[link, 'p_nom_extendable']]
        fixed_capacity = sum(n.links.at[link, 'p_nom'] for link in fixed_links)

        logging.info(
            f"[Storage capacity limits] Found {len(extendable_links)} extendable and {len(fixed_links)} fixed {battery_type_name} links")
        logging.info(f"[Storage capacity limits] Fixed capacity: {fixed_capacity} MW")

        if not extendable_links and fixed_capacity == 0:
            logging.warning(
                f"[Storage capacity limits] No extendable {battery_type_name} links found for {country} and no fixed capacity, skipping")
            return

        # Now add constraints
        try:
            # Adjust limits to account for fixed capacity
            adjusted_min = min_capacity - fixed_capacity if min_capacity is not None else None
            adjusted_max = max_capacity - fixed_capacity if max_capacity is not None else None

            logging.info(
                f"[Storage capacity limits] Adjusted limits for {battery_type_name}: min={adjusted_min} MW, max={adjusted_max} MW")

            # Skip if no extendable links but fixed capacity satisfies constraints
            if not extendable_links:
                if (adjusted_min is None or adjusted_min <= 0) and (
                        adjusted_max is None or fixed_capacity <= max_capacity):
                    logging.info(
                        f"[Storage capacity limits] Fixed capacity ({fixed_capacity} MW) satisfies constraints for {battery_type_name}")
                    return
                else:
                    logging.warning(
                        f"[Storage capacity limits] Fixed capacity ({fixed_capacity} MW) doesn't satisfy constraints for {battery_type_name}, but no extendable links found")
                    return

            # If min constraint is already satisfied by fixed capacity, set to 0
            if adjusted_min is not None and adjusted_min <= 0:
                logging.info(
                    f"[Storage capacity limits] Minimum constraint already satisfied by fixed capacity for {battery_type_name}")
                adjusted_min = 0

            # If max constraint is smaller than fixed capacity, issue warning
            if adjusted_max is not None and adjusted_max < 0:
                logging.warning(
                    f"[Storage capacity limits] Maximum constraint ({max_capacity} MW) is smaller than fixed capacity ({fixed_capacity} MW) for {battery_type_name}")
                # We'll still apply the constraint, but it might be infeasible

            # For minimum capacity: sum of p_nom_opt >= adjusted_min
            if adjusted_min is not None and adjusted_min > 0:
                min_cname = f"storage-capacity-min-{country}-{battery_type_name}"
                n.model.add_constraints(
                    n.model["Link-p_nom"].loc[extendable_links].sum() >= adjusted_min,
                    name=min_cname
                )
                logging.info(
                    f"[Storage capacity limits] Added minimum capacity constraint for {battery_type_name}: {adjusted_min} MW")

            # For maximum capacity: sum of p_nom_opt <= adjusted_max
            if adjusted_max is not None and adjusted_max >= 0:
                max_cname = f"storage-capacity-max-{country}-{battery_type_name}"
                n.model.add_constraints(
                    n.model["Link-p_nom"].loc[extendable_links].sum() <= adjusted_max,
                    name=max_cname
                )
                logging.info(
                    f"[Storage capacity limits] Added maximum capacity constraint for {battery_type_name}: {adjusted_max} MW")

            logging.info(
                f"[Storage capacity limits] Successfully added capacity constraints for {country}, {battery_type_name}")

            # Add a callback to check the final values after solving
            def log_final_battery_capacity():
                if n.links.p_nom_opt.empty:
                    logging.info(f"[Storage capacity limits] No p_nom_opt attribute in results")
                    return

                # Calculate total battery capacity
                total_extendable = 0
                for link in extendable_links:
                    if link in n.links.index and hasattr(n.links, 'p_nom_opt'):
                        total_extendable += n.links.p_nom_opt.get(link, 0)

                total_capacity = total_extendable + fixed_capacity
                logging.info(
                    f"[Storage capacity limits] Final {battery_type_name} capacity: {total_capacity} MW (extendable: {total_extendable} MW, fixed: {fixed_capacity} MW)")

                if min_capacity is not None and total_capacity < min_capacity:
                    logging.warning(
                        f"[Storage capacity limits] Final {battery_type_name} capacity ({total_capacity} MW) is less than minimum ({min_capacity} MW)")

                if max_capacity is not None and total_capacity > max_capacity:
                    logging.warning(
                        f"[Storage capacity limits] Final {battery_type_name} capacity ({total_capacity} MW) exceeds maximum ({max_capacity} MW)")

            # Store the callback function to be called after solving
            if not hasattr(n, 'post_solve_callbacks'):
                n.post_solve_callbacks = []
            n.post_solve_callbacks.append(log_final_battery_capacity)

        except Exception as e:
            logging.error(f"[Storage capacity limits] Error adding constraints for {battery_type_name}: {e}")
            import traceback
            logging.error(traceback.format_exc())

            # Try alternative approach using individual variables
            try:
                logging.info(f"[Storage capacity limits] Trying alternative approach for {battery_type_name}...")

                # Get p_nom variables for extendable links
                p_nom_vars = n.model["Link-p_nom"]

                # Create sum expression manually
                sum_expr = sum(p_nom_vars.loc[link] for link in extendable_links)

                # Add constraints
                if adjusted_min is not None and adjusted_min > 0:
                    n.model.add_constraints(sum_expr >= adjusted_min,
                                            name=f"storage-capacity-min-{country}-{battery_type_name}-alt")
                    logging.info(
                        f"[Storage capacity limits] Added minimum capacity constraint (alt) for {battery_type_name}: {adjusted_min} MW")

                if adjusted_max is not None and adjusted_max >= 0:
                    n.model.add_constraints(sum_expr <= adjusted_max,
                                            name=f"storage-capacity-max-{country}-{battery_type_name}-alt")
                    logging.info(
                        f"[Storage capacity limits] Added maximum capacity constraint (alt) for {battery_type_name}: {adjusted_max} MW")

                logging.info(
                    f"[Storage capacity limits] Successfully added capacity constraints using alternative approach for {battery_type_name}")

            except Exception as e2:
                logging.error(
                    f"[Storage capacity limits] Alternative approach also failed for {battery_type_name}: {e2}")
                import traceback
                logging.error(traceback.format_exc())

    except Exception as e:
        logging.error(f"[Storage capacity limits] Error processing {battery_type_name}: {e}")
        import traceback
        logging.error(traceback.format_exc())


def extra_functionality(
    n: pypsa.Network, snapshots: pd.DatetimeIndex, planning_horizons: str | None = None
) -> None:
    """
    Add custom constraints and functionality.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance with config and params attributes
    snapshots : pd.DatetimeIndex
        Simulation timesteps
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments ``opts`` and
    ``snakemake.config`` are expected to be attached to the network.
    """
    config = n.config
    constraints = config["solving"].get("constraints", {})
    if constraints["BAU"] and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if constraints["SAFE"] and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if constraints["CCL"] and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config, planning_horizons)

        # Apply conventional generators phase-out
    if constraints.get("generator_phaseout", {}).get("apply_constraints", True):
        logging.info("Applying conventional generators phase-out...")
        phase_out_conventional_generators(n, config)
    else:
        logging.info("CHP phase-out disabled in config, skipping")

        # Add resistive heater capacity constraints
    rh_enabled = False
    possible_paths = [
        config.get('solving', {}).get('constraints', {}).get('resistive_heaters', {}),
        config.get('constraints', {}).get('resistive_heaters', {}),
        config.get('resistive_heaters', {})
    ]

    for path_config in possible_paths:
        if path_config and path_config.get('apply_constraints', False):
            rh_enabled = True
            break

    if rh_enabled:
        logging.info("[Setup] Adding resistive heater constraints")
        add_resistive_heater_constraints(n)
    else:
        logging.info("[Setup] Resistive heater constraints disabled in config, skipping")

    # Add heat pump capacity constraints
    hp_enabled = False
    possible_paths = [
        config.get('solving', {}).get('constraints', {}).get('heat_pumps', {}),
        config.get('constraints', {}).get('heat_pumps', {}),
        config.get('heat_pumps', {})
    ]

    for path_config in possible_paths:
        if path_config and path_config.get('apply_constraints', False):
            hp_enabled = True
            break

    if hp_enabled:
        logging.info("[Setup] Adding heat pump constraints")
        add_heat_pump_constraints(n)
    else:
        logging.info("[Setup] Heat pump constraints disabled in config, skipping")

    # Heat pump flexibility constraints (separate from heat pump capacity constraints)
    hp_flex_enabled = config.get('sector', {}).get('heat_pump_flexibility', False)

    if hp_flex_enabled:
        logging.info("[Setup] Adding heat pump flexibility constraints")
        add_heat_pump_flexibility_constraints(n)
    else:
        logging.info("[Setup] Heat pump flexibility disabled in config, skipping")

    # Add storage capacity limits
    add_storage_capacity_limits(n)

    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)

    if EQ_o := constraints["EQ"]:
        add_EQ_constraints(n, EQ_o.replace("EQ", ""))

    if {"solar-hsat", "solar"}.issubset(
        config["electricity"]["renewable_carriers"]
    ) and {"solar-hsat", "solar"}.issubset(
        config["electricity"]["extendable_carriers"]["Generator"]
    ):
        add_solar_potential_constraints(n, config)

    add_battery_constraints(n)
    add_lossy_bidirectional_link_constraints(n)
    add_pipe_retrofit_constraint(n)
    if n._multi_invest:
        add_carbon_constraint(n, snapshots)
        add_carbon_budget_constraint(n, snapshots)
        add_retrofit_gas_boiler_constraint(n, snapshots)
    else:
        add_co2_atmosphere_constraint(n, snapshots)

    if config["sector"]["enhanced_geothermal"]["enable"]:
        add_flexible_egs_constraint(n)

    if n.params.custom_extra_functionality:
        source_path = pathlib.Path(n.params.custom_extra_functionality).resolve()
        assert source_path.exists(), f"{source_path} does not exist"
        sys.path.append(os.path.dirname(source_path))
        module_name = os.path.splitext(os.path.basename(source_path))[0]
        module = importlib.import_module(module_name)
        custom_extra_functionality = getattr(module, module_name)
        custom_extra_functionality(n, snapshots, snakemake)  # pylint: disable=E0601


def check_objective_value(n: pypsa.Network, solving: dict) -> None:
    """
    Check if objective value matches expected value within tolerance.

    Parameters
    ----------
    n : pypsa.Network
        Network with solved objective
    solving : Dict
        Dictionary containing objective checking parameters

    Raises
    ------
    ObjectiveValueError
        If objective value differs from expected value beyond tolerance
    """
    check_objective = solving["check_objective"]
    if check_objective["enable"]:
        atol = check_objective["atol"]
        rtol = check_objective["rtol"]
        expected_value = check_objective["expected_value"]
        if not np.isclose(n.objective, expected_value, atol=atol, rtol=rtol):
            raise ObjectiveValueError(
                f"Objective value {n.objective} differs from expected value "
                f"{expected_value} by more than {atol}."
            )


def solve_network(
    n: pypsa.Network,
    config: dict,
    params: dict,
    solving: dict,
    rule_name: str | None = None,
    planning_horizons: str | None = None,
    **kwargs,
) -> None:
    """
    Solve network optimization problem.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : Dict
        Configuration dictionary containing solver settings
    params : Dict
        Dictionary of solving parameters
    solving : Dict
        Dictionary of solving options and configuration
    rule_name : str, optional
        Name of the snakemake rule being executed
    planning_horizons : str, optional
            The current planning horizon year or None in perfect foresight
    **kwargs
        Additional keyword arguments passed to the solver

    Returns
    -------
    n : pypsa.Network
        Solved network instance
    status : str
        Solution status
    condition : str
        Termination condition

    Raises
    ------
    RuntimeError
        If solving status is infeasible
    ObjectiveValueError
        If objective value differs from expected value
    """
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    kwargs["multi_investment_periods"] = config["foresight"] == "perfect"
    kwargs["solver_options"] = (
        solving["solver_options"][set_of_options] if set_of_options else {}
    )
    kwargs["solver_name"] = solving["solver"]["name"]
    kwargs["extra_functionality"] = partial(
        extra_functionality, planning_horizons=planning_horizons
    )
    kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment", False
    )
    kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)
    kwargs["io_api"] = cf_solving.get("io_api", None)

    if config['run'].get('quick_test', False):
        start_date = '2019-01-1 00:00:00'
        end_date = '2019-01-30 00:00:00'
        selected_snapshots = n.snapshots[(n.snapshots >= start_date) & (n.snapshots <= end_date)]
        kwargs['snapshots'] = selected_snapshots
        n.snapshots = selected_snapshots

    if kwargs["solver_name"] == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    rolling_horizon = cf_solving.pop("rolling_horizon", False)
    skip_iterations = cf_solving.pop("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    # add to network for extra_functionality
    n.config = config
    n.params = params

    if rolling_horizon and rule_name == "solve_operations_network":
        kwargs["horizon"] = cf_solving.get("horizon", 365)
        kwargs["overlap"] = cf_solving.get("overlap", 0)
        n.optimize.optimize_with_rolling_horizon(**kwargs)
        status, condition = "", ""
    elif skip_iterations:
        status, condition = n.optimize(**kwargs)
    else:
        kwargs["track_iterations"] = cf_solving["track_iterations"]
        kwargs["min_iterations"] = cf_solving["min_iterations"]
        kwargs["max_iterations"] = cf_solving["max_iterations"]
        if cf_solving["post_discretization"].pop("enable"):
            logger.info("Add post-discretization parameters.")
            kwargs.update(cf_solving["post_discretization"])
        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
            **kwargs
        )

    if not rolling_horizon:
        if status != "ok":
            logger.warning(
                f"Solving status '{status}' with termination condition '{condition}'"
            )
        check_objective_value(n, solving)

    if "infeasible" in condition:
        labels = n.model.compute_infeasibilities()
        logger.info(f"Labels:\n{labels}")
        n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'")

    if status == "warning":
        raise RuntimeError(
            "Solving status 'warning'. Results may not be reliable. Aborting."
        )

    return n


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_sector_network_perfect",
            configfiles="../config/test/config.perfect.yaml",
            opts="",
            clusters="5",
            ll="v1.0",
            sector_opts="",
            # planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
    )

    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        solve_network(
            n,
            config=snakemake.config,
            params=snakemake.params,
            solving=snakemake.params.solving,
            planning_horizons=planning_horizons,
            rule_name=snakemake.rule,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    with open(snakemake.output.config, "w") as file:
        yaml.dump(
            n.meta,
            file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )