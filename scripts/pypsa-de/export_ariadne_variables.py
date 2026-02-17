import ast
import logging
import math
import os
import re
import sys
from functools import reduce

import numpy as np
import pandas as pd
import pypsa
from numpy import isclose
from pypsa.statistics import get_transmission_carriers

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + "/../.."))

from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import calculate_annuity, load_costs

pypsa.options.params.statistics.round = 10
logger = logging.getLogger(__name__)

# Defining global variables

TWh2MWh = 1e6
MWh2TWh = 1e-6
MW2GW = 1e-3
MW2TW = 1e-6
t2Mt = 1e-6
toe_to_MWh = 11.630  # GWh/ktoe OR MWh/toe


EUR20TOEUR23 = 1.1076


def domestic_length_factor(n, carriers, region="DE"):
    """
    Calculate the length factor for specified carriers within a PyPSA network.

    Parameters
    ----------
    n (pypsa.Network): The PyPSA network object.
    carriers (list or str): List of carrier types to filter, or a single carrier as a string.
    region (str): The region code to match in the buses (e.g., "DE").

    Returns
    -------
    float or dict: A single length factor if one carrier is provided; otherwise, a dictionary
                   of length factors for each carrier and component type.
    """
    # If a single carrier is provided as a string, wrap it in a list
    if isinstance(carriers, str):
        carriers = [carriers]

    length_factors = {}

    # Check if carriers exist in network components
    for carrier in carriers:
        if carrier not in (
            n.links.carrier.unique().tolist() + n.lines.carrier.unique().tolist()
        ):
            logger.info(f"Carrier '{carrier}' is neither in lines nor links.")
            continue  # Skip this carrier if not found in both links and lines

        # Loop through relevant components
        for c in n.components:
            if c.name in ["Link", "Line"] and carrier in c.static["carrier"].unique():
                # Filter based on carrier and region, excluding reversed links
                all_i = c.static[
                    (c.static["carrier"] == carrier)
                    & (c.static.bus0 + c.static.bus1).str.contains(region)
                    & ~c.static.index.str.contains("reversed")
                ].index

                # Separate domestic and cross-border links
                domestic_i = all_i[
                    c.static.loc[all_i, "bus0"].str.contains(region)
                    & c.static.loc[all_i, "bus1"].str.contains(region)
                ]
                cross_border_i = all_i.difference(domestic_i)

                # Ensure indices match expected totals
                assert len(all_i) == len(domestic_i) + len(cross_border_i)

                # Calculate length factor if both sets are non-empty
                if len(domestic_i) > 0 and len(cross_border_i) > 0:
                    length_factor = (
                        c.static.loc[domestic_i, "length"].mean()
                        / c.static.loc[cross_border_i, "length"].mean()
                    )
                    length_factors[(carrier, c.name)] = length_factor
                else:
                    logger.info(
                        f"No domestic or cross-border links found for {carrier} in {c.name}."
                    )

    # Return single length factor if only one carrier was provided and has a length factor
    if len(carriers) == 1 and len(length_factors) == 1:
        return next(iter(length_factors.values()))

    return length_factors


def _get_fuel_fractions(n, region, fuel):
    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "at_port": True,
        "nice_names": False,
    }

    renewable_fuel_supply = (
        n.statistics.supply(bus_carrier=f"renewable {fuel}", **kwargs)
        .groupby(["bus", "carrier"])
        .sum()
    ).round(3)  # rounding for numerical stability

    total_fuel_supply = (
        n.statistics.supply(bus_carrier=f"{fuel}", **kwargs)
        .groupby(["name", "carrier"])
        .sum()
    ).round(3)

    if fuel == "gas":
        fuel_refining = "gas compressing"
    elif fuel == "oil":
        fuel_refining = "oil refining"
    else:
        raise ValueError(f"Fuel {fuel} not supported")

    if "DE" in region:
        domestic_fuel_supply = (
            total_fuel_supply.reindex(
                [
                    (f"DE {fuel_refining}", f"{fuel_refining}"),
                    (f"DE renewable {fuel} -> DE {fuel}", f"renewable {fuel}"),
                    (f"EU renewable {fuel} -> DE {fuel}", f"renewable {fuel}"),
                ]
            )
            .dropna()
            .groupby("carrier")
            .sum()
        )  # If links are not used they are dropped here
        total_imported_renewable_fuel = total_fuel_supply.get(
            (f"EU renewable {fuel} -> DE {fuel}", f"renewable {fuel}"), 0
        )
        total_exported_renewable_fuel = total_fuel_supply.get(
            (f"DE renewable {fuel} -> EU {fuel}", f"renewable {fuel}"), 0
        )
        domestic_renewable_fuel = renewable_fuel_supply.get(f"DE renewable {fuel}")
        foreign_renewable_fuel = renewable_fuel_supply.get(f"EU renewable {fuel}")
    else:
        domestic_fuel_supply = (
            total_fuel_supply.reindex(
                [
                    (f"EU {fuel_refining}", f"{fuel_refining}"),
                    (f"DE renewable {fuel} -> EU {fuel}", f"renewable {fuel}"),
                    (f"EU renewable {fuel} -> EU {fuel}", f"renewable {fuel}"),
                ]
            )
            .dropna()
            .groupby("carrier")
            .sum()
        )
        total_imported_renewable_fuel = total_fuel_supply.get(
            (f"DE renewable {fuel} -> EU {fuel}", f"renewable {fuel}"), 0
        )
        total_exported_renewable_fuel = total_fuel_supply.get(
            (f"EU renewable {fuel} -> DE {fuel}", f"renewable {fuel}"), 0
        )
        domestic_renewable_fuel = renewable_fuel_supply.get(f"EU renewable {fuel}", 0)
        foreign_renewable_fuel = renewable_fuel_supply.get(f"DE renewable {fuel}", 0)

    # Legacy code for dropping gas pipelines (now deactivated)
    drops = ["gas pipeline", "gas pipeline new"]
    for d in drops:
        if d in total_fuel_supply.index:
            total_fuel_supply.drop(d, inplace=True)

    if total_imported_renewable_fuel == 0:
        imported_renewable_fuel = pd.Series(
            0, index=renewable_fuel_supply.index.get_level_values("carrier").unique()
        )
    else:
        imported_renewable_fuel = (
            foreign_renewable_fuel
            / foreign_renewable_fuel.sum()
            * total_imported_renewable_fuel
        )

    if total_exported_renewable_fuel == 0:
        exported_renewable_fuel = pd.Series(
            0, index=renewable_fuel_supply.index.get_level_values("carrier").unique()
        )
    else:
        exported_renewable_fuel = (
            domestic_renewable_fuel
            / domestic_renewable_fuel.sum()
            * total_exported_renewable_fuel
        )

    renewable_fuel_balance = imported_renewable_fuel.add(
        domestic_renewable_fuel, fill_value=0
    ).subtract(exported_renewable_fuel, fill_value=0)
    # These numbers may deviate because the input data from n.statistics
    # can deviate. I don't know why exactly, but at least we can check that
    # the difference stays roughly the same after the calculation.
    assert isclose(
        domestic_fuel_supply.get(f"renewable {fuel}", 0) - renewable_fuel_balance.sum(),
        total_fuel_supply.get([f"DE renewable {fuel} -> DE {fuel}"], pd.Series(0)).sum()
        + total_fuel_supply.get(
            [f"DE renewable {fuel} -> EU {fuel}"], pd.Series(0)
        ).sum()
        - renewable_fuel_supply.get(f"DE renewable {fuel}", pd.Series(0)).sum(),
        rtol=1e-3,
        atol=1e-5,
    )

    fuel_fractions = pd.Series()

    fuel_fractions["Biomass"] = renewable_fuel_balance.filter(like="bio").sum()

    if fuel == "gas":
        fuel_fractions["Natural Gas"] = domestic_fuel_supply.get("gas compressing", 0)
        fuel_fractions["Efuel"] = renewable_fuel_balance.get("Sabatier", 0)
    elif fuel == "oil":
        fuel_fractions["Fossil"] = domestic_fuel_supply.get("oil refining", 0)
        fuel_fractions["Efuel"] = renewable_fuel_balance.get("Fischer-Tropsch", 0)
    else:
        raise ValueError(f"Fuel {fuel} not supported")

    fuel_fractions = fuel_fractions.divide(domestic_fuel_supply.sum()).round(9)

    assert isclose(fuel_fractions.sum(), 1)

    fuel_fractions /= fuel_fractions.sum()

    return fuel_fractions


def _get_h2_fossil_fraction(n):
    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "at_port": True,
        "nice_names": False,
    }
    total_h2_supply = (
        n.statistics.supply(bus_carrier="H2", **kwargs)
        .drop("Store", errors="ignore")
        .groupby("carrier")
        .sum()
    )

    h2_fossil_fraction = (
        total_h2_supply.filter(like="SMR").sum() / total_h2_supply.sum()
    )

    return h2_fossil_fraction


def _get_t_sum(df, df_t, carrier, region, snapshot_weightings, port):
    if isinstance(carrier, list):
        return sum(
            [
                _get_t_sum(df, df_t, car, region, snapshot_weightings, port)
                for car in carrier
            ]
        )
    idx = df[df.carrier == carrier].filter(like=region, axis=0).index

    return (
        df_t[port][idx]
        .multiply(
            snapshot_weightings,
            axis=0,
        )
        .values.sum()
    )


def sum_load(n, carrier, region):
    return _get_t_sum(
        n.loads,
        n.loads_t,
        carrier,
        region,
        n.snapshot_weightings.generators,
        "p",
    )


def sum_co2(n, carrier, region):
    if isinstance(carrier, list):
        return sum([sum_co2(n, car, region) for car in carrier])
    try:
        port = (
            n.links.groupby("carrier")
            .first()
            .loc[carrier]
            .filter(like="bus")
            .tolist()
            .index("co2 atmosphere")
        )
    except KeyError:
        logger.info(
            "Warning: carrier `",
            carrier,
            "` not found in network.links.carrier!",
            sep="",
        )
        return 0

    return (
        -1
        * t2Mt
        * _get_t_sum(
            n.links,
            n.links_t,
            carrier,
            region,
            n.snapshot_weightings.generators,
            f"p{port}",
        )
    )


def get_total_co2(n, region):
    # including international bunker fuels and negative emissions
    df = n.links.filter(like=region, axis=0)
    co2 = 0
    for port in [col[3:] for col in df if col.startswith("bus")]:
        links = df.index[df[f"bus{port}"] == "co2 atmosphere"]
        co2 -= (
            n.links_t[f"p{port}"][links]
            .multiply(
                n.snapshot_weightings.generators,
                axis=0,
            )
            .values.sum()
        )
    return t2Mt * co2


def get_capacities(n, region):
    return _get_capacities(n, region, n.statistics.optimal_capacity)


def add_system_cost_rows(n):
    def fill_if_lifetime_inf(n, carrier, lifetime, component="links"):
        df = getattr(n, component)
        if df.loc[df.carrier == carrier, "lifetime"].sum() == np.inf:
            df.loc[df.carrier == carrier, "lifetime"] = lifetime
        else:
            logger.error(f"Mean lifetime of {carrier} is not infinite!")
            raise ValueError()

    logger.info("Overwriting lifetime of components to compute annuities")

    # Lines
    n.lines.lifetime = 40

    # Generators
    n.generators.loc[n.generators.lifetime == np.inf, "lifetime"] = 0
    n.generators.loc[n.generators.carrier == "ror", "lifetime"] = 60  # hydro lifetime

    # Links
    fill_if_lifetime_inf(n, "DC", 40)
    n.links.loc[n.links.lifetime == np.inf, "lifetime"] = 0
    n.links.loc[n.links.lifetime == 1, "lifetime"] = 0
    n.links.loc[n.links.carrier == "coal", "lifetime"] = 40

    # Stores
    n.stores.loc[
        (n.stores.lifetime == np.inf) & (n.stores.carrier == "H2 Store"), "lifetime"
    ] = 30  # hydrogen storage tank type 1 including compressor
    fill_if_lifetime_inf(n, "co2 stored", 25, "stores")
    n.stores.loc[n.stores.lifetime == np.inf, "lifetime"] = 0
    n.stores.loc[n.stores.carrier == "gas", "lifetime"] = 100
    n.stores.loc[n.stores.carrier == "oil", "lifetime"] = 30
    n.stores.loc[n.stores.carrier == "methanol", "lifetime"] = 30

    # Storage Units
    n.storage_units.lifetime = 60  # hydro lifetime

    for component in ["lines", "links", "generators", "stores", "storage_units"]:
        df = getattr(n, component)
        if df.empty:
            continue
        decentral_idx = df.index[df.index.str.contains("decentral|rural|rooftop")]
        not_decentral_idx = df.index[~df.index.str.contains("decentral|rural|rooftop")]

        for idx, discount_rate in zip([decentral_idx, not_decentral_idx], [0.04, 0.07]):
            df.loc[idx, "annuity"] = (
                calculate_annuity(df.loc[idx, "lifetime"], discount_rate)
                * df.loc[idx, "overnight_cost"]
            )

        df["FOM"] = df["capital_cost"] - df["annuity"]
        # Special case offwind, because it includes the connection
        if component == "generators":
            df.loc[
                df.carrier.str.contains("offwind"),
                "FOM",
            ] = (
                0.023185
                * df.loc[
                    df.carrier.str.contains("offwind"),
                    "overnight_cost",
                ]
            )
        if df["FOM"].min() < 0:
            logger.info(df["FOM"].min())
            logger.error(f"Capital cost is smaller than annuity for {component}")
            # n.links.query("carrier=='DC' and index.str.startswith('DC')")[["carrier","annuity","capital_cost","lifetime","FOM","build_year"]].sort_values("FOM")


"""
    get_system_cost(n, region)

    Calculate total investment, CAPEX, and OPEX in the given region.
"""


def get_system_cost(n, region):
    add_system_cost_rows(n)

    invest = _get_capacities(
        n,
        region,
        lambda **kwargs: n.statistics.expanded_capex(
            **kwargs, cost_attribute="overnight_cost"
        ),
        cap_string="Investment|Energy Supply|",
    )

    # Corresponds to the expanded_capex function
    grid_invest = get_grid_investments(n, region, scope="expanded")

    # For capex we subtract the capex of existing assets before 2020
    capex = _get_capacities(
        n,
        region,
        lambda **kwargs: n.statistics.capex(**kwargs, cost_attribute="annuity"),
        cap_string="System Cost|CAPEX|",
    )

    capex2020 = _get_capacities(
        networks[0],
        region,
        lambda **kwargs: networks[0].statistics.installed_capex(
            **kwargs, cost_attribute="annuity"
        ),
        cap_string="System Cost|CAPEX|",
    )

    # Subtracting all capex of assets built before 2020
    capex -= capex2020

    baseyear_grid_invest = get_grid_investments(n, region, scope="baseyear")

    # Assuming 40 years lifetime, 7% discount rate
    grid_capex = pd.Series(
        data=calculate_annuity(40, 0.07)
        * 5
        * baseyear_grid_invest.values,  # yearly invest for 5 years
        index=baseyear_grid_invest.index.str.replace(
            "Investment|Energy Supply|",
            "System Cost|CAPEX|",
        ),
    )

    # For FOM all existing capacities are considered
    fom = _get_capacities(
        n,
        region,
        lambda **kwargs: n.statistics.capex(**kwargs, cost_attribute="FOM"),
        cap_string="System Cost|FOM|",
    )

    vom = _get_capacities(
        n,
        region,
        n.statistics.opex,
        cap_string="System Cost|OPEX|",
    )

    opex = pd.Series(
        data=vom.values + fom.values,
        index=vom.index,
    )

    all_grid_invest = get_grid_investments(n, region, scope="all")

    # Assuming VOM=0, FOM=2% of Investment
    grid_opex = pd.Series(
        data=0.02 * 5 * all_grid_invest.values,  # yearly invest for 5 years
        index=all_grid_invest.index.str.replace(
            "Investment|Energy Supply|",
            "System Cost|OPEX|",
        ),
    )

    grid_fom = pd.Series(
        data=grid_opex.values,
        index=grid_opex.index.str.replace(
            "System Cost|OPEX|",
            "System Cost|FOM|",
        ),
    )

    for var, grid_var, var_name in zip(
        [invest, capex, opex, fom],
        [grid_invest, grid_capex, grid_opex, grid_fom],
        ["Investment|Energy Supply|", "System Cost|CAPEX|", "System Cost|OPEX|"],
    ):
        var[var_name + "Electricity"] += grid_var[
            var_name + "Electricity|Transmission and Distribution"
        ]
        var[var_name + "Hydrogen"] += grid_var[var_name + "Hydrogen|Transmission"]
        if var_name + "Gas|Transmission" in grid_var.keys():
            var[var_name + "Gas"] += grid_var[var_name + "Gas|Transmission"]

    return pd.concat(
        [invest, grid_invest, capex, fom, grid_capex, opex, grid_opex, grid_fom]
    )


def get_installed_capacities(n, region):
    return _get_capacities(
        n, region, n.statistics.installed_capacity, cap_string="Installed Capacity|"
    )


def get_capacity_additions_simple(n, region):
    caps = get_capacities(n, region)
    incaps = get_installed_capacities(n, region)
    return pd.Series(
        data=caps.values - incaps.values,
        index="Capacity Additions" + caps.index.str[8:],
    )


def get_capacity_additions(n, region):
    def _f(**kwargs):
        return n.statistics.optimal_capacity(**kwargs).sub(
            n.statistics.installed_capacity(**kwargs), fill_value=0
        )

    return _get_capacities(n, region, _f, cap_string="Capacity Additions|")


def get_capacity_additions_nstat(n, region):
    def _f(*args, **kwargs):
        kwargs.pop("storage", None)
        return n.statistics.expanded_capacity(*args, **kwargs)

    return _get_capacities(n, region, _f, cap_string="Capacity Additions Nstat|")


def _get_capacities(n, region, cap_func, cap_string="Capacity|"):
    kwargs = {
        "groupby": ["bus", "carrier"],
        "at_port": True,
        "nice_names": False,
    }

    var = pd.Series()

    capacities_electricity = (
        cap_func(
            bus_carrier=["AC", "low voltage"],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .drop(
            # transmission capacities
            ["AC", "DC", "electricity distribution grid"],
            errors="ignore",
        )
        .multiply(MW2GW)
    )

    var[cap_string + "Electricity|Biomass|w/ CCS"] = capacities_electricity.get(
        "urban central solid biomass CHP CC ", 0
    )

    var[cap_string + "Electricity|Biomass|w/o CCS"] = capacities_electricity.reindex(
        ["urban central solid biomass CHP", "solid biomass", "biogas"]
    ).sum()

    var[cap_string + "Electricity|Biomass|Solids"] = capacities_electricity.filter(
        like="solid biomass"
    ).sum()

    var[cap_string + "Electricity|Biomass|Gases and Liquids"] = (
        capacities_electricity.get("biogas", 0)
    )

    var[cap_string + "Electricity|Biomass"] = (
        var[cap_string + "Electricity|Biomass|Solids"]
        + var[cap_string + "Electricity|Biomass|Gases and Liquids"]
    )

    var[cap_string + "Electricity|Non-Renewable Waste"] = capacities_electricity.filter(
        like="waste CHP"
    ).sum()

    var[cap_string + "Electricity|Coal|Hard Coal"] = capacities_electricity.filter(
        like="coal"
    ).sum()

    var[cap_string + "Electricity|Coal|Lignite"] = capacities_electricity.filter(
        like="lignite"
    ).sum()

    # var[cap_string + "Electricity|Coal|Hard Coal|w/ CCS"] =
    # var[cap_string + "Electricity|Coal|Hard Coal|w/o CCS"] =
    # var[cap_string + "Electricity|Coal|Lignite|w/ CCS"] =
    # var[cap_string + "Electricity|Coal|Lignite|w/o CCS"] =
    # var[cap_string + "Electricity|Coal|w/ CCS"] =
    # var[cap_string + "Electricity|Coal|w/o CCS"] =
    # Q: CCS for coal Implemented, but not activated, should we use it?
    # !: No, because of Kohleausstieg
    # > config: coal_cc

    var[cap_string + "Electricity|Coal"] = var[
        [
            cap_string + "Electricity|Coal|Lignite",
            cap_string + "Electricity|Coal|Hard Coal",
        ]
    ].sum()

    # var[cap_string + "Electricity|Gas|CC|w/ CCS"] =
    # var[cap_string + "Electricity|Gas|CC|w/o CCS"] =
    # ! Not implemented, rarely used

    var[cap_string + "Electricity|Gas|CC"] = capacities_electricity.get("CCGT", 0)

    var[cap_string + "Electricity|Gas|OC"] = capacities_electricity.get("OCGT", 0)

    var[cap_string + "Electricity|Gas|w/ CCS"] = capacities_electricity.get(
        "urban central gas CHP CC", 0
    )

    var[cap_string + "Electricity|Gas|w/o CCS"] = (
        capacities_electricity.get("urban central gas CHP", 0)
        + var[
            [
                cap_string + "Electricity|Gas|CC",
                cap_string + "Electricity|Gas|OC",
            ]
        ].sum()
    )

    var[cap_string + "Electricity|Gas"] = var[
        [
            cap_string + "Electricity|Gas|w/ CCS",
            cap_string + "Electricity|Gas|w/o CCS",
        ]
    ].sum()

    # var[cap_string + "Electricity|Geothermal"] =
    # ! Not implemented

    var[cap_string + "Electricity|Hydro"] = capacities_electricity.reindex(
        ["ror", "hydro"]
    ).sum()

    # Q!: Not counting PHS here, because it is a true storage,
    # as opposed to hydro

    # var[cap_string + "Electricity|Hydrogen|CC"] =
    # ! Not implemented
    # var[cap_string + "Electricity|Hydrogen|OC"] =
    # Q: "H2-turbine"
    # Q: What about retrofitted gas power plants? -> Lisa
    var[cap_string + "Electricity|Hydrogen|CC"] = capacities_electricity.reindex(
        [
            "H2 CCGT",
            "urban central H2 CHP",
            "H2 retrofit CCGT",
            "urban central H2 retrofit CHP",
        ]
    ).sum()

    var[cap_string + "Electricity|Hydrogen|OC"] = capacities_electricity.reindex(
        [
            "H2 OCGT",
            "H2 retrofit OCGT",
        ]
    ).sum()

    var[cap_string + "Electricity|Hydrogen|FC"] = capacities_electricity.get(
        "H2 Fuel Cell", 0
    )

    var[cap_string + "Electricity|Hydrogen"] = (
        var[cap_string + "Electricity|Hydrogen|CC"]
        + var[cap_string + "Electricity|Hydrogen|OC"]
        + var[cap_string + "Electricity|Hydrogen|FC"]
    )

    var[cap_string + "Electricity|Nuclear"] = capacities_electricity.get("nuclear", 0)

    # var[cap_string + "Electricity|Ocean"] =
    # ! Not implemented

    # var[cap_string + "Electricity|Oil|w/ CCS"] =
    # var[cap_string + "Electricity|Oil|w/o CCS"] =
    # ! Not implemented

    var[cap_string + "Electricity|Oil"] = capacities_electricity.filter(
        like="oil"
    ).sum()

    var[cap_string + "Electricity|Solar|PV|Rooftop"] = capacities_electricity.get(
        "solar rooftop", 0
    )

    var[cap_string + "Electricity|Solar|PV|Open Field"] = (
        capacities_electricity.reindex(["solar", "solar-hsat"]).sum()
    )

    var[cap_string + "Electricity|Solar|PV"] = var[
        [
            cap_string + "Electricity|Solar|PV|Open Field",
            cap_string + "Electricity|Solar|PV|Rooftop",
        ]
    ].sum()

    # var[cap_string + "Electricity|Solar|CSP"] =
    # ! not implemented

    var[cap_string + "Electricity|Solar"] = var[cap_string + "Electricity|Solar|PV"]

    var[cap_string + "Electricity|Wind|Offshore"] = capacities_electricity.reindex(
        ["offwind", "offwind-ac", "offwind-dc", "offwind-float"]
    ).sum()

    var[cap_string + "Electricity|Wind|Onshore"] = capacities_electricity.get("onwind")

    var[cap_string + "Electricity|Wind"] = capacities_electricity.filter(
        like="wind"
    ).sum()

    assert (
        var[cap_string + "Electricity|Wind"]
        == var[
            [
                cap_string + "Electricity|Wind|Offshore",
                cap_string + "Electricity|Wind|Onshore",
            ]
        ].sum()
    )

    # var[cap_string + "Electricity|Storage Converter|CAES"] =
    # ! Not implemented

    var[cap_string + "Electricity|Storage Converter|Hydro Dam Reservoir"] = (
        capacities_electricity.get("hydro", 0)
    )

    var[cap_string + "Electricity|Storage Converter|Pump Hydro"] = (
        capacities_electricity.get("PHS", 0)
    )

    var[cap_string + "Electricity|Storage Converter|Stationary Batteries|Utility"] = (
        capacities_electricity.get("battery discharger", 0)
    )

    var[cap_string + "Electricity|Storage Converter|Stationary Batteries|Home"] = (
        capacities_electricity.get("home battery discharger", 0)
    )

    var[cap_string + "Electricity|Storage Converter|Stationary Batteries"] = var[
        [
            cap_string + "Electricity|Storage Converter|Stationary Batteries|Utility",
            cap_string + "Electricity|Storage Converter|Stationary Batteries|Home",
        ]
    ].sum()

    var[cap_string + "Electricity|Storage Converter|Vehicles"] = (
        capacities_electricity.get("V2G", 0)
    )

    var[cap_string + "Electricity|Storage Converter"] = var[
        [
            cap_string + "Electricity|Storage Converter|Hydro Dam Reservoir",
            cap_string + "Electricity|Storage Converter|Pump Hydro",
            cap_string + "Electricity|Storage Converter|Stationary Batteries",
            cap_string + "Electricity|Storage Converter|Vehicles",
        ]
    ].sum()

    if cap_string.startswith("Investment") or cap_string.startswith("System Cost"):
        storage_capacities = (
            cap_func(
                **kwargs,
            )
            .filter(like=region)
            .groupby("carrier")
            .sum()
            .multiply(MW2GW)
        )
    else:
        storage_capacities = (
            cap_func(
                storage=True,
                **kwargs,
            )
            .filter(like=region)
            .groupby("carrier")
            .sum()
            .multiply(MW2GW)
        )
    # var[cap_string + "Electricity|Storage Reservoir|CAES"] =
    # ! Not implemented

    var[cap_string + "Electricity|Storage Reservoir|Hydro Dam Reservoir"] = (
        storage_capacities.get("hydro")
    )

    var[cap_string + "Electricity|Storage Reservoir|Pump Hydro"] = (
        storage_capacities.get("PHS")
    )

    var[cap_string + "Electricity|Storage Reservoir|Stationary Batteries|Utility"] = (
        storage_capacities.get("battery", 0)
    )

    var[cap_string + "Electricity|Storage Reservoir|Stationary Batteries|Home"] = (
        storage_capacities.get("home battery", 0)
    )

    var[cap_string + "Electricity|Storage Reservoir|Stationary Batteries"] = var[
        [
            cap_string + "Electricity|Storage Reservoir|Stationary Batteries|Utility",
            cap_string + "Electricity|Storage Reservoir|Stationary Batteries|Home",
        ]
    ].sum()

    var[cap_string + "Electricity|Storage Reservoir|Vehicles"] = storage_capacities.get(
        "EV battery", 0
    )

    var[cap_string + "Electricity|Storage Reservoir"] = var[
        [
            cap_string + "Electricity|Storage Reservoir|Hydro Dam Reservoir",
            cap_string + "Electricity|Storage Reservoir|Pump Hydro",
            cap_string + "Electricity|Storage Reservoir|Stationary Batteries",
            cap_string + "Electricity|Storage Reservoir|Vehicles",
        ]
    ].sum()

    var[cap_string + "Electricity"] = var[
        [
            cap_string + "Electricity|Wind",
            cap_string + "Electricity|Solar",
            cap_string + "Electricity|Oil",
            cap_string + "Electricity|Coal",
            cap_string + "Electricity|Gas",
            cap_string + "Electricity|Biomass",
            cap_string + "Electricity|Hydro",
            cap_string + "Electricity|Hydrogen",
            cap_string + "Electricity|Nuclear",
            cap_string + "Electricity|Non-Renewable Waste",
        ]
    ].sum()

    capacities_central_heat = (
        cap_func(
            bus_carrier=[
                "urban central heat",
            ],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .drop(
            ["urban central heat vent"],
            errors="ignore",  # drop existing labels or do nothing
        )
        .multiply(MW2GW)
    )
    heat_pump_idxs = capacities_central_heat.filter(like="heat pump").index
    capacities_central_heat[heat_pump_idxs] = abs(
        capacities_central_heat[heat_pump_idxs]
    )

    if cap_string.startswith("Investment") or cap_string.startswith("System Cost"):
        secondary_heat_techs = [
            "DAC",
            "Fischer-Tropsch",
            "H2 Electrolysis",
            "H2 Fuel Cell",
            "methanolisation",
            "Sabatier",
            "CHP",  # We follow REMIND convention and account all CHPs only in electricity
        ]
        secondary_heat_idxs = [
            idx
            for idx in capacities_central_heat.index
            if any([tech in idx for tech in secondary_heat_techs])
        ]
        capacities_central_heat[secondary_heat_idxs] = 0

    var[cap_string + "Heat|Solar thermal"] = capacities_central_heat.filter(
        like="solar thermal"
    ).sum()

    # !!! Missing in the Ariadne database
    #  We could be much more detailed for the heat sector (as for electricity)
    # if desired by Ariadne
    #
    var[cap_string + "Heat|Biomass|w/ CCS"] = capacities_central_heat.get(
        "urban central solid biomass CHP CC", 0
    )
    var[cap_string + "Heat|Biomass|w/o CCS"] = (
        capacities_central_heat.get("urban central solid biomass CHP", 0)
        + capacities_central_heat.filter(like="biomass boiler").sum()
    )

    var[cap_string + "Heat|Biomass"] = (
        var[cap_string + "Heat|Biomass|w/ CCS"]
        + var[cap_string + "Heat|Biomass|w/o CCS"]
    )

    assert isclose(
        var[cap_string + "Heat|Biomass"],
        capacities_central_heat.filter(like="biomass").sum(),
    )

    var[cap_string + "Heat|Non-Renewable Waste"] = capacities_central_heat.filter(
        like="waste CHP"
    ).sum()

    var[cap_string + "Heat|Resistive heater"] = capacities_central_heat.filter(
        like="resistive heater"
    ).sum()

    var[cap_string + "Heat|Processes"] = pd.Series(
        {
            c: capacities_central_heat.get(c)
            for c in [
                "Fischer-Tropsch",
                "H2 Electrolysis",
                "H2 Fuel Cell",
                "Sabatier",
                "methanolisation",
            ]
        }
    ).sum()

    var[cap_string + "Heat|Hydrogen"] = capacities_central_heat.reindex(
        ["urban central H2 CHP", "urban central H2 retrofit CHP"]
    ).sum()

    # !!! Missing in the Ariadne database

    var[cap_string + "Heat|Gas"] = (
        capacities_central_heat.filter(like="gas boiler").sum()
        + capacities_central_heat.filter(like="gas CHP").sum()
    )

    # var[cap_string + "Heat|Geothermal"] =
    # ! Not implemented

    var[cap_string + "Heat|Heat pump"] = capacities_central_heat.filter(
        like="heat pump"
    ).sum()

    var[cap_string + "Heat|Oil"] = capacities_central_heat.filter(
        like="oil boiler"
    ).sum()

    var[cap_string + "Heat|Storage Converter"] = capacities_central_heat[
        capacities_central_heat.index.str.contains("water (?:tanks|pits) discharger")
    ].sum()

    var[cap_string + "Heat|Storage Reservoir"] = storage_capacities[
        storage_capacities.index.str.contains("water (?:tanks|pits)")
    ].sum()

    var[cap_string + "Heat"] = (
        var[cap_string + "Heat|Solar thermal"]
        + var[cap_string + "Heat|Resistive heater"]
        + var[cap_string + "Heat|Biomass"]
        + var[cap_string + "Heat|Oil"]
        + var[cap_string + "Heat|Gas"]
        + var[cap_string + "Heat|Processes"]
        + var[cap_string + "Heat|Hydrogen"]
        + var[cap_string + "Heat|Heat pump"]
        + var[cap_string + "Heat|Non-Renewable Waste"]
    )

    var[cap_string + "Heat|Renewable"] = (
        var[cap_string + "Heat|Solar thermal"]
        + var[cap_string + "Heat|Biomass"]
        + var[cap_string + "Heat|Hydrogen"]
        + var[cap_string + "Heat|Heat pump"]
        + var[cap_string + "Heat|Resistive heater"]
    )

    capacities_decentral_heat = (
        cap_func(
            bus_carrier=[
                "urban decentral heat",
                "rural heat",
            ],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .drop(
            ["DAC"],
            errors="ignore",  # drop existing labels or do nothing
        )
        .multiply(MW2GW)
    )
    heat_pump_idxs = capacities_decentral_heat.filter(like="heat pump").index
    capacities_decentral_heat[heat_pump_idxs] = abs(
        capacities_decentral_heat[heat_pump_idxs]
    )

    var[cap_string + "Decentral Heat|Solar thermal"] = capacities_decentral_heat.filter(
        like="solar thermal"
    ).sum()

    var[cap_string + "Decentral Heat|Heat Pump"] = capacities_decentral_heat.filter(
        like="heat pump"
    ).sum()

    var[cap_string + "Decentral Heat|Resistive Heater"] = (
        capacities_decentral_heat.filter(like="resistive heater").sum()
    )

    var[cap_string + "Decentral Heat|Biomass"] = capacities_decentral_heat.filter(
        like="biomass"
    ).sum()

    var[cap_string + "Decentral Heat|Gas"] = capacities_decentral_heat.filter(
        like="gas"
    ).sum()

    var[cap_string + "Decentral Heat|Oil"] = capacities_decentral_heat.filter(
        like="oil boiler"
    ).sum()

    var[cap_string + "Decentral Heat|Storage Converter"] = capacities_decentral_heat[
        capacities_decentral_heat.index.str.contains("water tanks discharger")
    ].sum()

    var[cap_string + "Decentral Heat|Storage Reservoir"] = storage_capacities[
        storage_capacities.index.str.contains("(?:decentral|rural) water tanks")
    ].sum()

    capacities_h2 = (
        cap_func(
            bus_carrier="H2",
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(MW2GW)
    )

    var[cap_string + "Hydrogen|Gas|w/ CCS"] = capacities_h2.get("SMR CC", 0)

    var[cap_string + "Hydrogen|Gas|w/o CCS"] = capacities_h2.get("SMR", 0)

    var[cap_string + "Hydrogen|Gas"] = capacities_h2.filter(like="SMR").sum()

    assert (
        var[cap_string + "Hydrogen|Gas"]
        == var[cap_string + "Hydrogen|Gas|w/ CCS"]
        + var[cap_string + "Hydrogen|Gas|w/o CCS"]
    )

    var[cap_string + "Hydrogen|Electricity"] = abs(
        capacities_electricity.get("H2 Electrolysis", 0)
    )

    var[cap_string + "Hydrogen"] = (
        capacities_h2.get("H2 Electrolysis", 0) + var[cap_string + "Hydrogen|Gas"]
    )

    var[cap_string + "Hydrogen|Storage Reservoir"] = var[
        cap_string + "Hydrogen|Reservoir"
    ] = storage_capacities.get("H2 Store", 0)

    capacities_gas = (
        cap_func(
            bus_carrier="renewable gas",
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(MW2GW)
    )

    var[cap_string + "Gases|Hydrogen"] = capacities_gas.get("Sabatier", 0)

    var[cap_string + "Gases|Biomass"] = capacities_gas.reindex(
        [
            "biogas to gas",
            "biogas to gas CC",
        ]
    ).sum()

    var[cap_string + "Gases"] = (
        var[cap_string + "Gases|Hydrogen"] + var[cap_string + "Gases|Biomass"]
    )

    # This check requires further changes to n.statistics
    #
    # assert isclose(
    #     var[cap_string + "Gases"],
    #     capacities_gas.sum(),
    # )

    capacities_liquids = (
        cap_func(
            bus_carrier="renewable oil",
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(MW2GW)
    )
    #
    var[cap_string + "Liquids|Hydrogen"] = capacities_liquids.get("Fischer-Tropsch", 0)

    var[cap_string + "Liquids|Biomass"] = capacities_liquids.filter(
        like="biomass to liquid"
    ).sum()

    try:
        capacities_methanol = (
            cap_func(
                bus_carrier="methanol",
                **kwargs,
            )
            .filter(like=region)
            .groupby("carrier")
            .sum()
            .multiply(MW2GW)
        )
        #
        var[cap_string + "Methanol"] = capacities_methanol.get("methanolisation", 0)
    except KeyError:
        logger.info(
            "Warning: carrier `methanol` not found in network.links.carrier! Assuming 0 capacities."
        )
        var[cap_string + "Methanol"] = 0

    var[cap_string + "Liquids|Hydrogen"] += var[cap_string + "Methanol"]

    var[cap_string + "Liquids"] = (
        var[cap_string + "Liquids|Hydrogen"] + var[cap_string + "Liquids|Biomass"]
    )

    if cap_string.startswith("Investment"):
        var = var.div(MW2GW).mul(1e-9).div(5).round(3)  # in bn € / year
    elif cap_string.startswith("System Cost"):
        var = var.div(MW2GW).mul(1e-9).round(3)

    return var


def get_CHP_E_and_H_usage(n, bus_carrier, region, fossil_fraction=1):
    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "nice_names": False,
    }

    usage = (
        n.statistics.withdrawal(
            bus_carrier=bus_carrier,
            **kwargs,
        )
        .filter(like=region)
        .filter(like="CHP")
        .multiply(fossil_fraction)
    )

    E_to_H = (
        n.links.loc[usage.index.get_level_values("name")].efficiency
        / n.links.loc[usage.index.get_level_values("name")].efficiency2
    )

    E_fraction = E_to_H * (1 / (E_to_H + 1))

    E_usage = usage.multiply(E_fraction).sum()
    H_usage = usage.multiply(1 - E_fraction).sum()

    return E_usage, H_usage


def get_primary_energy(n, region):
    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "nice_names": False,
    }

    var = pd.Series()

    oil_fossil_fraction = _get_fuel_fractions(n, region, "oil")["Fossil"]
    primary_oil_factor = (
        n.links.query("carrier=='oil refining'").efficiency.unique().item()
    )
    oil_usage = (
        n.statistics.withdrawal(bus_carrier="oil", **kwargs)
        .filter(like=region)
        .drop("Store", errors="ignore")
        .groupby("carrier")
        .sum()
        .multiply(oil_fossil_fraction)
    )

    oil_CHP_E_usage, oil_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "oil", region, fossil_fraction=oil_fossil_fraction
    )

    ## Primary Energy

    var["Primary Energy|Oil|Heat"] = (
        oil_usage.filter(like="urban central oil boiler").sum() + oil_CHP_H_usage
    ) / primary_oil_factor

    var["Primary Energy|Oil|Electricity"] = (
        oil_usage.get("oil", 0) + oil_CHP_E_usage
    ) / primary_oil_factor

    var["Primary Energy|Oil"] = oil_usage.sum() / primary_oil_factor

    # At the moment, everything that is not electricity or heat is counted as liquid fuel
    var["Primary Energy|Oil|Liquids"] = (
        var["Primary Energy|Oil"]
        - var["Primary Energy|Oil|Electricity"]
        - var["Primary Energy|Oil|Heat"]
    )

    assert isclose(
        var["Primary Energy|Oil"],
        n.statistics.withdrawal(bus_carrier="oil primary", **kwargs)
        .get(("Link", "DE oil refining"), pd.Series(0))
        .item(),
        atol=1,  # MW
    )

    gas_fractions = _get_fuel_fractions(n, region, "gas")

    gas_usage = (
        n.statistics.withdrawal(
            bus_carrier="gas",
            **kwargs,
        )
        .filter(like=region)
        .drop(
            "Store",
            errors="ignore",
        )
        .groupby("carrier")
        .sum()
        .multiply(gas_fractions["Natural Gas"])
    )
    primary_gas_factor = (
        n.links.query("carrier=='gas compressing'").efficiency.unique().item()
    )

    gas_CHP_E_usage, gas_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "gas", region, fossil_fraction=gas_fractions["Natural Gas"]
    )

    var["Primary Energy|Gas|Heat"] = (
        gas_usage.filter(like="urban central gas boiler").sum() + gas_CHP_H_usage
    ) / primary_gas_factor

    var["Primary Energy|Gas|Electricity"] = (
        gas_usage.reindex(
            [
                "CCGT",
                "OCGT",
            ],
        ).sum()
        + gas_CHP_E_usage
    ) / primary_gas_factor

    var["Primary Energy|Gas|Hydrogen"] = (
        gas_usage.filter(like="SMR").sum() / primary_gas_factor
    )

    var["Primary Energy|Gas"] = gas_usage.sum() / primary_gas_factor

    _gas_primary = (
        n.statistics.withdrawal(
            bus_carrier="gas primary",
            **kwargs,
        )
        .get(("Link", "DE gas compressing"), pd.Series(0))
        .item()
    )

    assert isclose(
        var["Primary Energy|Gas"],
        _gas_primary,
        atol=1e-1,
    )

    var["Primary Energy|Gas|Gases"] = (
        var["Primary Energy|Gas"]
        - var["Primary Energy|Gas|Electricity"]
        - var["Primary Energy|Gas|Heat"]
        - var["Primary Energy|Gas|Hydrogen"]
    )

    waste_CHP_E_usage, waste_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "non-sequestered HVC", region
    )

    var["Primary Energy|Waste|Electricity"] = waste_CHP_E_usage
    var["Primary Energy|Waste|Heat"] = waste_CHP_H_usage
    var["Primary Energy|Waste"] = (
        var["Primary Energy|Waste|Electricity"] + var["Primary Energy|Waste|Heat"]
    )

    coal_usage = (
        n.statistics.withdrawal(
            bus_carrier=["lignite", "coal"],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
    )

    coal_CHP_E_usage, coal_CHP_H_usage = get_CHP_E_and_H_usage(n, "coal", region)
    lignite_CHP_E_usage, lignite_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "lignite", region
    )

    var["Primary Energy|Coal|Hard Coal"] = (
        coal_usage.get("coal", 0) + coal_CHP_E_usage + coal_CHP_H_usage
    )

    var["Primary Energy|Coal|Lignite"] = (
        coal_usage.get("lignite", 0) + lignite_CHP_E_usage + lignite_CHP_H_usage
    )

    var["Primary Energy|Coal|Electricity"] = (
        var["Primary Energy|Coal|Hard Coal"]
        - coal_CHP_H_usage
        + var["Primary Energy|Coal|Lignite"]
        - lignite_CHP_H_usage
    )

    var["Primary Energy|Coal|Heat"] = coal_CHP_H_usage + lignite_CHP_H_usage

    var["Primary Energy|Coal"] = (
        var["Primary Energy|Coal|Heat"]
        + var["Primary Energy|Coal|Electricity"]
        + coal_usage.get("coal for industry", 0)
    )

    assert isclose(var["Primary Energy|Coal"], coal_usage.sum())

    var["Primary Energy|Fossil"] = (
        var["Primary Energy|Coal"]
        + var["Primary Energy|Gas"]
        + var["Primary Energy|Oil"]
    )

    biomass_usage = (
        n.statistics.withdrawal(
            bus_carrier=["solid biomass", "biogas"],
            **kwargs,
        )
        .filter(like=region)
        .drop("Store", errors="ignore")
        .groupby("carrier")
        .sum()
    )

    biomass_CHP_E_usage, biomass_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "solid biomass", region
    )

    var["Primary Energy|Biomass|Gases"] = biomass_usage.filter(
        like="biogas to gas"
    ).sum()

    # btl_efficiency = n.links.query(
    #         "carrier == 'biomass to liquid' and (build_year == 2020)"
    #     ).efficiency.unique().item()

    unsus_btl_secondary = (
        n.statistics.supply(
            bus_carrier=["renewable oil"],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .get("unsustainable bioliquids", 0)
    )

    var["Primary Energy|Biomass|Liquids"] = (
        biomass_usage.filter(like="biomass to liquid").sum()
        + unsus_btl_secondary  # divide by the BtL efficiency 2020?
    )

    var["Primary Energy|Biomass|w/ CCS"] = biomass_usage[
        biomass_usage.index.str.contains("CC")
    ].sum()

    var["Primary Energy|Biomass|w/o CCS"] = biomass_usage[
        ~biomass_usage.index.str.contains("CC")
    ].sum()

    var["Primary Energy|Biomass|Electricity"] = (
        biomass_CHP_E_usage + biomass_usage.reindex(["solid biomass", "biogas"]).sum()
    )
    var["Primary Energy|Biomass|Heat"] = biomass_CHP_H_usage + biomass_usage.get(
        "urban central solid biomass boiler", 0
    )

    var["Primary Energy|Biomass|Solids"] = (
        biomass_usage.filter(like="solid biomass for industry").sum()
        + biomass_usage.get("rural biomass boiler", 0)
        + biomass_usage.get("urban decentral biomass boiler", 0)
    )

    var["Primary Energy|Biomass"] = (
        var["Primary Energy|Biomass|Electricity"]
        + var["Primary Energy|Biomass|Heat"]
        + var["Primary Energy|Biomass|Liquids"]
        + var["Primary Energy|Biomass|Solids"]
        + var["Primary Energy|Biomass|Gases"]
    )

    assert isclose(
        var["Primary Energy|Biomass"],
        biomass_usage.sum() + unsus_btl_secondary,
    )

    var["Primary Energy|Nuclear"] = (
        n.statistics.withdrawal(
            bus_carrier=["uranium"],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .get("nuclear", 0)
    )

    # ! This should basically be equivalent to secondary energy
    renewable_electricity = (
        n.statistics.supply(
            bus_carrier=["AC", "low voltage"],
            **kwargs,
        )
        .drop(
            [
                # Assuming renewables are only generators and StorageUnits
                "Link",
                "Line",
            ],
            errors="ignore",
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
    )

    solar_thermal_heat = (
        n.statistics.supply(
            bus_carrier=[
                "urban decentral heat",
                "urban central heat",
                "rural heat",
            ],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .filter(like="solar thermal")
        .sum()
    )

    var["Primary Energy|Hydro"] = renewable_electricity.reindex(
        [
            "ror",
            "PHS",
            "hydro",
        ]
    ).sum()

    var["Primary Energy|Solar"] = (
        renewable_electricity.filter(like="solar").sum() + solar_thermal_heat
    )

    var["Primary Energy|Wind"] = renewable_electricity.filter(like="wind").sum()

    assert isclose(
        renewable_electricity.sum() + solar_thermal_heat,
        (
            var["Primary Energy|Hydro"]
            + var["Primary Energy|Solar"]
            + var["Primary Energy|Wind"]
        ),
    )
    # Primary Energy|Other
    # Not implemented

    var["Primary Energy"] = (
        var["Primary Energy|Fossil"]
        + var["Primary Energy|Biomass"]
        + var["Primary Energy|Hydro"]
        + var["Primary Energy|Solar"]
        + var["Primary Energy|Wind"]
        + var["Primary Energy|Nuclear"]
        + var["Primary Energy|Waste"]
    )

    return var * MWh2TWh


def get_secondary_energy(n, region, _industry_demand):
    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "nice_names": False,
    }
    var = pd.Series()

    electricity_supply = (
        n.statistics.supply(bus_carrier=["low voltage", "AC"], **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .drop(
            ["AC", "DC", "electricity distribution grid"],
            errors="ignore",
        )
    )

    var["Secondary Energy|Electricity|Coal|Hard Coal"] = electricity_supply.filter(
        like="coal"
    ).sum()

    var["Secondary Energy|Electricity|Coal|Lignite"] = electricity_supply.filter(
        like="lignite"
    ).sum()

    var["Secondary Energy|Electricity|Coal"] = (
        var["Secondary Energy|Electricity|Coal|Hard Coal"]
        + var["Secondary Energy|Electricity|Coal|Lignite"]
    )

    var["Secondary Energy|Electricity|Oil"] = electricity_supply.filter(
        like="oil"
    ).sum()

    var["Secondary Energy|Electricity|Gas"] = electricity_supply.reindex(
        [
            "CCGT",
            "OCGT",
            "urban central gas CHP",
            "urban central gas CHP CC",
        ],
    ).sum()

    gas_fractions = _get_fuel_fractions(n, region, "gas")

    var["Secondary Energy|Electricity|Gas|Fossil"] = (
        var["Secondary Energy|Electricity|Gas"] * gas_fractions["Natural Gas"]
    )
    var["Secondary Energy|Electricity|Gas|Biomass"] = (
        var["Secondary Energy|Electricity|Gas"] * gas_fractions["Biomass"]
    )
    var["Secondary Energy|Electricity|Gas|Hydrogen"] = (
        var["Secondary Energy|Electricity|Gas"] * gas_fractions["Efuel"]
    )

    var["Secondary Energy|Electricity|Fossil"] = (
        var["Secondary Energy|Electricity|Gas|Fossil"]
        + var["Secondary Energy|Electricity|Oil"]
        + var["Secondary Energy|Electricity|Coal"]
    )

    var["Secondary Energy|Electricity|Biomass|w/o CCS"] = electricity_supply.reindex(
        ["urban central solid biomass CHP", "solid biomass", "biogas"]
    ).sum()
    var["Secondary Energy|Electricity|Biomass|w/ CCS"] = electricity_supply.get(
        "urban central solid biomass CHP CC", 0
    )
    var["Secondary Energy|Electricity|Biomass|Solid"] = electricity_supply.filter(
        like="solid biomass"
    ).sum()
    var["Secondary Energy|Electricity|Biomass|Gaseous and Liquid"] = (
        electricity_supply.get("biogas", 0)
    )
    var["Secondary Energy|Electricity|Biomass"] = (
        var["Secondary Energy|Electricity|Biomass|w/o CCS"]
        + var["Secondary Energy|Electricity|Biomass|w/ CCS"]
    )

    var["Secondary Energy|Electricity|Hydro"] = electricity_supply.reindex(
        [
            "hydro",
            "ror",
        ]
    ).sum()
    # ! Neglecting PHS here because it is storage infrastructure

    var["Secondary Energy|Electricity|Nuclear"] = electricity_supply.filter(
        like="nuclear"
    ).sum()
    var["Secondary Energy|Electricity|Solar"] = electricity_supply.filter(
        like="solar"
    ).sum()
    var["Secondary Energy|Electricity|Wind|Offshore"] = electricity_supply.filter(
        like="offwind"
    ).sum()
    var["Secondary Energy|Electricity|Wind|Onshore"] = electricity_supply.get("onwind")
    var["Secondary Energy|Electricity|Wind"] = (
        var["Secondary Energy|Electricity|Wind|Offshore"]
        + var["Secondary Energy|Electricity|Wind|Onshore"]
    )
    var["Secondary Energy|Electricity|Non-Biomass Renewables"] = (
        var["Secondary Energy|Electricity|Hydro"]
        + var["Secondary Energy|Electricity|Solar"]
        + var["Secondary Energy|Electricity|Wind"]
    )

    var["Secondary Energy|Electricity|Hydrogen"] = electricity_supply.reindex(
        [
            "H2 Fuel Cell",
            "H2 OCGT",
            "H2 CCGT",
            "urban central H2 CHP",
            "H2 retrofit OCGT",
            "H2 retrofit CCGT",
            "urban central H2 retrofit CHP",
        ]
    ).sum()
    # ! Add H2 Turbines if they get implemented

    var["Secondary Energy|Electricity|Waste"] = electricity_supply.filter(
        like="waste CHP"
    ).sum()

    var["Secondary Energy|Electricity|Curtailment"] = (
        n.statistics.curtailment(bus_carrier=["AC", "low voltage"], **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .reindex(
            [
                "offwind-ac",
                "offwind-dc",
                "onwind",
                "solar",
                "solar rooftop",
                "solar-hsat",
            ]
        )
        .sum()
    )

    electricity_balance = (
        n.statistics.energy_balance(bus_carrier=["AC", "low voltage"], **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    if "V2G" in electricity_balance.index:
        logger.error(
            "The exporter requires changes to correctly account vehicle to grid technology."
        )
    var["Secondary Energy|Electricity|Storage Losses"] = (
        -1
        * electricity_balance.reindex(
            [
                "battery charger",
                "battery discharger",
                "home battery charger",
                "home battery discharger",
                "PHS",
            ]
        ).sum()
    )

    idx = n.lines.query(
        f"bus0.str.contains('{region}') and bus1.str.contains('{region}')"
    ).index
    acl = (
        n.lines_t.p0[idx].sum(axis=1).multiply(n.snapshot_weightings.generators).sum()
        + n.lines_t.p1[idx].sum(axis=1).multiply(n.snapshot_weightings.generators).sum()
    )
    idx = n.links.query(
        f"bus0.str.contains('{region}') and bus1.str.contains('{region}') and carrier == 'DC'"
    ).index
    dcl = (
        n.links_t.p0[idx].sum(axis=1).multiply(n.snapshot_weightings.generators).sum()
        + n.links_t.p1[idx].sum(axis=1).multiply(n.snapshot_weightings.generators).sum()
    )
    idx = n.links.query(
        f"bus0.str.contains('{region}') and bus1.str.contains('{region}') and carrier == 'electricity distribution grid'"
    ).index
    distril = (
        n.links_t.p0[idx].sum(axis=1).multiply(n.snapshot_weightings.generators).sum()
        + n.links_t.p1[idx].sum(axis=1).multiply(n.snapshot_weightings.generators).sum()
    )

    var["Secondary Energy|Electricity|Transmission Losses"] = acl + dcl + distril

    var["Secondary Energy|Electricity"] = (
        var[
            "Secondary Energy|Electricity|Gas"
        ]  # fossil gas, biomethane, sabatier methane
        + var["Secondary Energy|Electricity|Oil"]
        + var["Secondary Energy|Electricity|Coal"]
        + var["Secondary Energy|Electricity|Biomass"]
        + var["Secondary Energy|Electricity|Non-Biomass Renewables"]
        + var["Secondary Energy|Electricity|Nuclear"]
        + var["Secondary Energy|Electricity|Hydrogen"]
        + var["Secondary Energy|Electricity|Waste"]
    )

    assert isclose(
        electricity_supply[
            ~electricity_supply.index.str.contains(
                "PHS|battery discharger|home battery discharger|V2G"
            )
        ].sum(),
        var["Secondary Energy|Electricity"],
    )

    heat_supply = (
        n.statistics.supply(
            bus_carrier=[
                "urban central heat",
                # rural and urban decentral heat do not produce secondary energy
            ],
            **kwargs,
        )
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    var["Secondary Energy|Heat|Gas"] = heat_supply.filter(like="gas").sum()

    var["Secondary Energy|Heat|Biomass"] = heat_supply.filter(like="biomass").sum()

    var["Secondary Energy|Heat|Coal"] = (
        heat_supply.filter(like="coal").sum() + heat_supply.filter(like="lignite").sum()
    )
    # var["Secondary Energy|Heat|Geothermal"] = \
    # var["Secondary Energy|Heat|Nuclear"] = \
    # var["Secondary Energy|Heat|Other"] = \
    # ! Not implemented
    var["Secondary Energy|Heat|Hydrogen"] = heat_supply.filter(
        like="urban central H2"
    ).sum()

    var["Secondary Energy|Heat|Oil"] = heat_supply.filter(
        like="urban central oil"
    ).sum()

    var["Secondary Energy|Heat|Solar"] = heat_supply.filter(like="solar thermal").sum()

    var["Secondary Energy|Heat|Electricity|Heat Pumps"] = heat_supply.filter(
        like="heat pump"
    ).sum()
    var["Secondary Energy|Heat|Electricity|Resistive"] = heat_supply.filter(
        like="resistive heater"
    ).sum()
    var["Secondary Energy|Heat|Electricity"] = (
        var["Secondary Energy|Heat|Electricity|Heat Pumps"]
        + var["Secondary Energy|Heat|Electricity|Resistive"]
    )
    var["Secondary Energy|Heat|Waste"] = heat_supply.filter(like="waste CHP").sum()
    var["Secondary Energy|Heat|Other"] = heat_supply.reindex(
        [
            "Fischer-Tropsch",
            "H2 Fuel Cell",
            "H2 Electrolysis",
            "Sabatier",
            "methanolisation",
        ]
    ).sum()
    # TODO remember to specify in comments

    var["Secondary Energy|Heat"] = (
        var["Secondary Energy|Heat|Gas"]
        + var["Secondary Energy|Heat|Biomass"]
        + var["Secondary Energy|Heat|Oil"]
        + var["Secondary Energy|Heat|Solar"]
        + var["Secondary Energy|Heat|Electricity"]
        + var["Secondary Energy|Heat|Other"]
        + var["Secondary Energy|Heat|Coal"]
        + var["Secondary Energy|Heat|Waste"]
        + var["Secondary Energy|Heat|Hydrogen"]
    )
    assert isclose(
        var["Secondary Energy|Heat"],
        heat_supply[~heat_supply.index.str.contains("discharger")].sum(),
    )

    hydrogen_production = (
        n.statistics.supply(bus_carrier="H2", **kwargs)
        .filter(like=region)
        .drop("Store", errors="ignore")
        .groupby(["carrier"])
        .sum()
    )

    var["Secondary Energy|Hydrogen|Electricity"] = hydrogen_production.get(
        "H2 Electrolysis", 0
    )

    var["Secondary Energy|Hydrogen|Gas"] = hydrogen_production.reindex(
        ["SMR", "SMR CC"]
    ).sum()

    var["Secondary Energy|Hydrogen|Other"] = hydrogen_production.get(
        "H2 for industry", 0
    )

    var["Secondary Energy|Hydrogen"] = (
        var["Secondary Energy|Hydrogen|Electricity"]
        + var["Secondary Energy|Hydrogen|Gas"]
        + var["Secondary Energy|Hydrogen|Other"]
    )

    assert isclose(
        var["Secondary Energy|Hydrogen"],
        hydrogen_production[
            ~hydrogen_production.index.str.startswith("H2 pipeline")
        ].sum(),
        rtol=0.01,
        atol=1e-3,
    )

    # Liquids
    liquids_production = (
        n.statistics.supply(bus_carrier=["oil", "renewable oil", "methanol"], **kwargs)
        .filter(like=region)
        .drop(
            [("Store", "DE oil Store"), ("Store", "DE methanol Store")], errors="ignore"
        )
        .groupby(["carrier"])
        .sum()
        .drop(["renewable oil", "methanol"], errors="ignore")  # Drop trade links
    )
    var["Secondary Energy|Liquids|Fossil"] = var["Secondary Energy|Liquids|Oil"] = (
        liquids_production.get("oil refining", 0)
    )
    var["Secondary Energy|Methanol"] = liquids_production.get("methanolisation", 0)
    var["Secondary Energy|Liquids|Hydrogen"] = liquids_production.get(
        "Fischer-Tropsch", 0
    ) + liquids_production.get("methanolisation", 0)

    var["Secondary Energy|Liquids|Biomass"] = liquids_production.filter(
        like="bio"
    ).sum()

    var["Secondary Energy|Liquids"] = (
        var["Secondary Energy|Liquids|Oil"]
        + var["Secondary Energy|Liquids|Hydrogen"]
        + var["Secondary Energy|Liquids|Biomass"]
    )
    assert isclose(
        var["Secondary Energy|Liquids"],
        liquids_production.sum(),
        rtol=0.01,
        atol=1e-3,
    )

    gas_supply = (
        n.statistics.supply(bus_carrier=["gas", "renewable gas"], **kwargs)
        .filter(like=region)
        .drop(("Store", "DE gas Store"), errors="ignore")
        .groupby(["carrier"])
        .sum()
        .drop(["renewable gas"], errors="ignore")
    )

    # Fraction supplied by Hydrogen conversion
    var["Secondary Energy|Gases|Hydrogen"] = gas_supply.get("Sabatier", 0)

    var["Secondary Energy|Gases|Biomass"] = gas_supply.filter(like="bio").sum()

    var["Secondary Energy|Gases|Natural Gas"] = gas_supply.get("gas compressing", 0)

    var["Secondary Energy|Gases"] = (
        var["Secondary Energy|Gases|Hydrogen"]
        + var["Secondary Energy|Gases|Biomass"]
        + var["Secondary Energy|Gases|Natural Gas"]
    )

    assert isclose(var["Secondary Energy|Gases"], gas_supply.sum())

    industry_demand = _industry_demand.filter(
        like=region,
        axis=0,
    ).sum()
    mwh_coal_per_mwh_coke = 1.366

    coke_fraction = (
        industry_demand.get("coke")
        * mwh_coal_per_mwh_coke
        / (
            industry_demand.get("coke") * mwh_coal_per_mwh_coke
            + industry_demand.get("coal")
        )
    )

    # Coke is added as a coal demand, so we need to convert back to units of coke for secondary energy
    var["Secondary Energy|Solids|Coal"] = var["Secondary Energy|Solids"] = (
        sum_load(n, "coal for industry", region) * coke_fraction / mwh_coal_per_mwh_coke
    )

    biomass_usage = (
        n.statistics.withdrawal(bus_carrier="solid biomass", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    var["Secondary Energy|Solids|Biomass"] = (
        biomass_usage.get("urban decentral biomass boiler", 0)
        + biomass_usage.get("rural biomass boiler", 0)
        + biomass_usage.filter(like="solid biomass for industry").sum()
    )

    electricity_withdrawal = (
        n.statistics.withdrawal(bus_carrier=["low voltage", "AC"], **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    var["Secondary Energy Input|Electricity|Hydrogen"] = electricity_withdrawal.get(
        "H2 Electrolysis", 0
    )

    var["Secondary Energy Input|Electricity|Heat"] = electricity_withdrawal.filter(
        like="urban central"
    ).sum()

    var["Secondary Energy Input|Electricity|Liquids"] = electricity_withdrawal.get(
        "methanolisation", 0
    )

    hydrogen_withdrawal = (
        n.statistics.withdrawal(bus_carrier="H2", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    H2_CHP_E_usage, H2_CHP_H_usage = get_CHP_E_and_H_usage(n, "H2", region)

    var["Secondary Energy Input|Hydrogen|Electricity"] = (
        hydrogen_withdrawal.reindex(
            [
                "H2 Fuel Cell",
                "H2 OCGT",
                "H2 CCGT",
                "H2 retrofit OCGT",
                "H2 retrofit CCGT",
            ]
        ).sum()
        + H2_CHP_E_usage
    )

    var["Secondary Energy Input|Hydrogen|Heat"] = H2_CHP_H_usage

    var["Secondary Energy Input|Hydrogen|Gases"] = hydrogen_withdrawal.get(
        "Sabatier", 0
    )

    var["Secondary Energy Input|Hydrogen|Liquids"] = hydrogen_withdrawal.reindex(
        ["Fischer-Tropsch", "methanolisation"]
    ).sum()

    var["Secondary Energy"] = (
        var["Secondary Energy|Electricity"]
        + var["Secondary Energy|Heat"]
        + var["Secondary Energy|Hydrogen"]
        + var["Secondary Energy|Gases"]
        + var["Secondary Energy|Liquids"]
        + var["Secondary Energy|Solids"]
    )

    return var * MWh2TWh


def get_final_energy(
    n,
    region,
    _industry_demand,
    _energy_totals,
    _sector_ratios,
    _industry_production,
    config,
    config_industry,
):
    var = pd.Series()

    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "nice_names": False,
    }
    h2_fossil_fraction = _get_h2_fossil_fraction(n)
    oil_fractions = _get_fuel_fractions(n, region, "oil")

    if config["sector"]["ammonia"]:
        # MWh/a
        Haber_Bosch_NH3 = (
            n.statistics.supply(bus_carrier="NH3", **kwargs)
            .groupby("carrier")
            .sum()["Haber-Bosch"]
        )

        CH4_for_NH3 = (
            Haber_Bosch_NH3
            * h2_fossil_fraction
            * config_industry["MWh_CH4_per_tNH3_SMR"]
            / config_industry["MWh_NH3_per_tNH3"]
        )
        H2_for_NH3 = (
            Haber_Bosch_NH3
            * (1 - h2_fossil_fraction)
            / config_industry["MWh_H2_per_tNH3_electrolysis"]
        )
        subcategories = ["HVC", "Methanol", "Chlorine"]

    else:
        CH4_for_NH3 = 0
        H2_for_NH3 = 0
        subcategories = ["HVC", "Methanol", "Chlorine", "Ammonia"]

    carrier = ["hydrogen", "methane", "naphtha"]

    ip = _industry_production.loc[region, subcategories]  # kt/a
    sr = _sector_ratios["DE"].loc[carrier, subcategories]  # MWh/tMaterial
    non_energy = sr.multiply(ip).sum(axis=1) * 1e3

    # write var
    var["Final Energy|Non-Energy Use|Gases"] = non_energy.methane + CH4_for_NH3

    gas_fractions = _get_fuel_fractions(n, region, "gas")
    for gas_type in gas_fractions.index:
        var[f"Final Energy|Non-Energy Use|Gases|{gas_type}"] = (
            var["Final Energy|Non-Energy Use|Gases"] * gas_fractions[gas_type]
        )

    var["Final Energy|Non-Energy Use|Liquids"] = non_energy.naphtha

    var["Final Energy|Non-Energy Use|Liquids|Petroleum"] = (
        non_energy.naphtha * oil_fractions["Fossil"]
    )
    var["Final Energy|Non-Energy Use|Liquids|Efuel"] = (
        non_energy.naphtha * oil_fractions["Efuel"]
    )
    var["Final Energy|Non-Energy Use|Liquids|Biomass"] = (
        non_energy.naphtha * oil_fractions["Biomass"]
    )

    var["Final Energy|Non-Energy Use|Solids"] = 0
    var["Final Energy|Non-Energy Use|Solids|Coal"] = 0
    var["Final Energy|Non-Energy Use|Solids|Biomass"] = 0

    var["Final Energy|Non-Energy Use|Hydrogen"] = non_energy.hydrogen + H2_for_NH3

    var["Final Energy|Non-Energy Use"] = non_energy.sum() + CH4_for_NH3 + H2_for_NH3

    assert isclose(
        var["Final Energy|Non-Energy Use"],
        var["Final Energy|Non-Energy Use|Gases"]
        + var["Final Energy|Non-Energy Use|Liquids"]
        + var["Final Energy|Non-Energy Use|Solids"]
        + var["Final Energy|Non-Energy Use|Hydrogen"],
    )

    energy_totals = _energy_totals.loc[region[0:2]]

    industry_demand = _industry_demand.filter(
        like=region,
        axis=0,
    ).sum()

    # !: Pypsa-eur does not strictly distinguish between energy and
    # non-energy use

    var["Final Energy|Industry|Electricity"] = sum_load(
        n, "industry electricity", region
    )
    # electricity is not used for non-energy purposes
    var["Final Energy|Industry excl Non-Energy Use|Electricity"] = var[
        "Final Energy|Industry|Electricity"
    ]

    var["Final Energy|Industry|Heat"] = sum_load(
        n, "low-temperature heat for industry", region
    )
    # heat is not used for non-energy purposes
    var["Final Energy|Industry excl Non-Energy Use|Heat"] = var[
        "Final Energy|Industry|Heat"
    ]

    # var["Final Energy|Industry|Solar"] = \
    # !: Included in |Heat

    # var["Final Energy|Industry|Geothermal"] = \
    # Not implemented

    var["Final Energy|Industry|Gases"] = sum_load(n, "gas for industry", region)

    for gas_type in gas_fractions.index:
        var[f"Final Energy|Industry|Gases|{gas_type}"] = (
            var["Final Energy|Industry|Gases"] * gas_fractions[gas_type]
        )

    # "gas for industry" is now regionally resolved and could be used here
    # subtract non-energy used gases from total gas demand
    var["Final Energy|Industry excl Non-Energy Use|Gases"] = (
        var["Final Energy|Industry|Gases"] - var["Final Energy|Non-Energy Use|Gases"]
    )

    for gas_type in gas_fractions.index:
        var[f"Final Energy|Industry excl Non-Energy Use|Gases|{gas_type}"] = (
            var["Final Energy|Industry excl Non-Energy Use|Gases"]
            * gas_fractions[gas_type]
        )

    # var["Final Energy|Industry|Power2Heat"] = \
    # Q: misleading description

    var["Final Energy|Industry|Hydrogen"] = sum_load(n, "H2 for industry", region)
    # subtract non-energy used hydrogen from total hydrogen demand
    var["Final Energy|Industry excl Non-Energy Use|Hydrogen"] = (
        var["Final Energy|Industry|Hydrogen"]
        - var["Final Energy|Non-Energy Use|Hydrogen"]
    )

    var["Final Energy|Industry|Liquids|Petroleum"] = (
        sum_load(n, "naphtha for industry", region) * oil_fractions["Fossil"]
    )

    # subtract non-energy used petroleum from total petroleum demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids|Petroleum"] = (
        var["Final Energy|Industry|Liquids|Petroleum"]
        - var["Final Energy|Non-Energy Use|Liquids|Petroleum"]
    )

    var["Final Energy|Industry|Liquids|Efuel"] = (
        sum_load(n, "naphtha for industry", region) * oil_fractions["Efuel"]
    )
    # subtract non-energy used efuels from total efuels demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids|Efuel"] = (
        var["Final Energy|Industry|Liquids|Efuel"]
        - var["Final Energy|Non-Energy Use|Liquids|Efuel"]
    )

    var["Final Energy|Industry|Liquids|Biomass"] = (
        sum_load(n, "naphtha for industry", region) * oil_fractions["Biomass"]
    )
    # subtract non-energy used biomass from total biomass demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids|Biomass"] = (
        var["Final Energy|Industry|Liquids|Biomass"]
        - var["Final Energy|Non-Energy Use|Liquids|Biomass"]
    )

    var["Final Energy|Industry|Liquids"] = sum_load(n, "naphtha for industry", region)
    # subtract non-energy used liquids from total liquid demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids"] = (
        var["Final Energy|Industry|Liquids"]
        - var["Final Energy|Non-Energy Use|Liquids"]
    )

    # var["Final Energy|Industry|Other"] = \

    var["Final Energy|Industry|Solids|Biomass"] = sum_load(
        n, "solid biomass for industry", region
    )
    var["Final Energy|Industry excl Non-Energy Use|Solids|Biomass"] = var[
        "Final Energy|Industry|Solids|Biomass"
    ]

    mwh_coal_per_mwh_coke = 1.366

    coke_fraction = (
        industry_demand.get("coke")
        * mwh_coal_per_mwh_coke
        / (
            industry_demand.get("coke") * mwh_coal_per_mwh_coke
            + industry_demand.get("coal")
        )
    )
    # Contains coke demand, which is a coal product
    # Here coke is considered a secondary energy source
    var["Final Energy|Industry|Solids|Coal"] = (
        sum_load(n, "coal for industry", region) * (1 - coke_fraction)
        + sum_load(n, "coal for industry", region)
        * coke_fraction
        / mwh_coal_per_mwh_coke
    )
    var["Final Energy|Industry excl Non-Energy Use|Solids|Coal"] = var[
        "Final Energy|Industry|Solids|Coal"
    ]

    var["Final Energy|Industry|Solids"] = (
        var["Final Energy|Industry|Solids|Biomass"]
        + var["Final Energy|Industry|Solids|Coal"]
    )
    # no solids used for non-energy purposes
    var["Final Energy|Industry excl Non-Energy Use|Solids"] = var[
        "Final Energy|Industry|Solids"
    ]

    # Why is AMMONIA zero?

    # var["Final Energy|Industry excl Non-Energy Use|Non-Metallic Minerals"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Chemicals"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Steel"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Steel|Primary"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Steel|Secondary"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Pulp and Paper"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Food and Tobacco"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Non-Ferrous Metals"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Engineering"] = \
    # var["Final Energy|Industry excl Non-Energy Use|Vehicle Construction"] = \
    # Q: Most of these could be found somewhere, but are model inputs!

    var["Final Energy|Industry"] = var.get(
        [
            "Final Energy|Industry|Electricity",
            "Final Energy|Industry|Heat",
            "Final Energy|Industry|Gases",
            "Final Energy|Industry|Hydrogen",
            "Final Energy|Industry|Liquids",
            "Final Energy|Industry|Solids",
        ]
    ).sum()
    var["Final Energy|Industry excl Non-Energy Use"] = var.get(
        [
            "Final Energy|Industry excl Non-Energy Use|Electricity",
            "Final Energy|Industry excl Non-Energy Use|Heat",
            "Final Energy|Industry excl Non-Energy Use|Gases",
            "Final Energy|Industry excl Non-Energy Use|Hydrogen",
            "Final Energy|Industry excl Non-Energy Use|Liquids",
            "Final Energy|Industry excl Non-Energy Use|Solids",
        ]
    ).sum()
    assert isclose(
        var["Final Energy|Industry"] - var["Final Energy|Non-Energy Use"],
        var["Final Energy|Industry excl Non-Energy Use"],
    )
    # Final energy is delivered to the consumers
    low_voltage_electricity = (
        n.statistics.withdrawal(
            bus_carrier="low voltage",
            **kwargs,
        )
        .filter(
            like=region,
        )
        .groupby("carrier")
        .sum()
    )

    rescom_electricity = low_voltage_electricity[
        # carrier does not contain one of the following substrings
        ~low_voltage_electricity.index.str.contains(
            "urban central|industry|agriculture|charger|distribution"
            # Excluding chargers (battery and EV)
        )
    ]
    var["Final Energy|Residential and Commercial|Electricity|Heat Pumps"] = (
        rescom_electricity.filter(like="heat pump").sum()
    )
    var["Final Energy|Residential and Commercial|Electricity"] = (
        rescom_electricity.sum()
    )
    # urban decentral heat and rural heat are delivered as different forms of energy
    # (gas, oil, biomass, ...)
    decentral_heat_withdrawal = (
        n.statistics.withdrawal(
            bus_carrier=["rural heat", "urban decentral heat"],
            **kwargs,
        )
        .filter(
            like=region,
        )
        .groupby("carrier")
        .sum()
        .drop(
            [  # chargers affect all sectors equally
                "urban decentral water tanks charger",
                "urban decentral water pits charger",
                "rural water tanks charger",
                "rural water pits charger",
            ],
            errors="ignore",
        )
    )

    decentral_heat_residential_and_commercial_fraction = (
        sum_load(n, ["urban decentral heat", "rural heat"], region)
        / decentral_heat_withdrawal.sum()
    )

    decentral_heat_supply_rescom = (
        n.statistics.supply(
            bus_carrier=["rural heat", "urban decentral heat"],
            **kwargs,
        )
        .filter(
            like=region,
        )
        .groupby("carrier")
        .sum()
        .multiply(decentral_heat_residential_and_commercial_fraction)
    )
    # Dischargers probably should not be considered, to avoid double counting

    var["Final Energy|Residential and Commercial|Heat"] = (
        sum_load(
            n, "urban central heat", region
        )  # For urban central Final Energy is delivered as Heat
        + decentral_heat_supply_rescom.filter(like="solar thermal").sum()
    )  # Assuming for solar thermal secondary energy == Final energy

    gas_usage = (
        n.statistics.withdrawal(bus_carrier="gas", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    # !!! Here the final is delivered as gas, not as heat
    var["Final Energy|Residential and Commercial|Gases"] = gas_usage.get(
        "urban decentral gas boiler", 0
    ) + gas_usage.get("rural gas boiler", 0)

    for gas_type in gas_fractions.index:
        var[f"Final Energy|Residential and Commercial|Gases|{gas_type}"] = (
            var["Final Energy|Residential and Commercial|Gases"]
            * gas_fractions[gas_type]
        )

    # var["Final Energy|Residential and Commercial|Hydrogen"] = \
    # ! Not implemented
    oil_usage = (
        n.statistics.withdrawal(bus_carrier="oil", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    var["Final Energy|Residential and Commercial|Liquids"] = oil_usage.get(
        "urban decentral oil boiler", 0
    ) + oil_usage.get("rural oil boiler", 0)

    var["Final Energy|Residential and Commercial|Liquids|Petroleum"] = (
        var["Final Energy|Residential and Commercial|Liquids"] * oil_fractions["Fossil"]
    )
    var["Final Energy|Residential and Commercial|Liquids|Efuel"] = (
        var["Final Energy|Residential and Commercial|Liquids"] * oil_fractions["Efuel"]
    )
    var["Final Energy|Residential and Commercial|Liquids|Biomass"] = (
        var["Final Energy|Residential and Commercial|Liquids"]
        * oil_fractions["Biomass"]
    )

    # var["Final Energy|Residential and Commercial|Other"] = \
    # var["Final Energy|Residential and Commercial|Solids|Coal"] = \
    # ! Not implemented

    biomass_usage = (
        n.statistics.withdrawal(bus_carrier="solid biomass", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
    )

    var["Final Energy|Residential and Commercial|Solids"] = var[
        "Final Energy|Residential and Commercial|Solids|Biomass"
    ] = biomass_usage.get("urban decentral biomass boiler", 0) + biomass_usage.get(
        "rural biomass boiler", 0
    )

    # Q: Everything else seems to be not implemented

    var["Final Energy|Residential and Commercial"] = (
        var["Final Energy|Residential and Commercial|Electricity"]
        + var["Final Energy|Residential and Commercial|Heat"]
        + var["Final Energy|Residential and Commercial|Gases"]
        + var["Final Energy|Residential and Commercial|Liquids"]
        + var["Final Energy|Residential and Commercial|Solids"]
    )

    # TODO double check prices for usage of correct FE carrier

    var["Final Energy|Residential and Commercial|Space and Water Heating"] = (
        # district heating
        var["Final Energy|Residential and Commercial|Heat"]
        # decentral boilers
        + var["Final Energy|Residential and Commercial|Gases"]
        + var["Final Energy|Residential and Commercial|Liquids"]
        + var["Final Energy|Residential and Commercial|Solids"]
        # resistive heaters and heat pumps
        + low_voltage_electricity.filter(like="rural").sum()
        + low_voltage_electricity.filter(like="urban decentral").sum()
    )

    # var["Final Energy|Transportation|Other"] = \

    var["Final Energy|Transportation|Electricity"] = low_voltage_electricity.get(
        "BEV charger", 0
    )

    # var["Final Energy|Transportation|Gases"] = \
    # var["Final Energy|Transportation|Gases|Natural Gas"] = \
    # var["Final Energy|Transportation|Gases|Biomass"] = \
    # var["Final Energy|Transportation|Gases|Efuel"] = \
    # var["Final Energy|Transportation|Gases|Synthetic Fossil"] = \
    # ! Not implemented

    var["Final Energy|Transportation|Hydrogen"] = sum_load(
        n, "land transport fuel cell", region
    )
    # ?? H2 for shipping

    international_aviation_fraction = energy_totals["total international aviation"] / (
        energy_totals["total domestic aviation"]
        + energy_totals["total international aviation"]
    )
    international_navigation_fraction = energy_totals[
        "total international navigation"
    ] / (
        energy_totals["total domestic navigation"]
        + energy_totals["total international navigation"]
    )
    var["Final Energy|Transportation|Domestic Aviation"] = var[
        "Final Energy|Transportation|Domestic Aviation|Liquids"
    ] = sum_load(n, "kerosene for aviation", region) * (
        1 - international_aviation_fraction
    )
    var["Final Energy|Transportation|Domestic Navigation"] = var[
        "Final Energy|Transportation|Domestic Navigation|Liquids"
    ] = sum_load(n, "shipping oil", region) * (1 - international_navigation_fraction)

    var["Final Energy|Transportation|Methanol"] = sum_load(
        n, "shipping methanol", region
    ) * (1 - international_navigation_fraction)

    var["Final Energy|Transportation|Liquids"] = (
        sum_load(n, "land transport oil", region)
        + var["Final Energy|Transportation|Domestic Aviation|Liquids"]
        + var["Final Energy|Transportation|Domestic Navigation|Liquids"]
        + var["Final Energy|Transportation|Methanol"]
    )

    # var["Final Energy|Transportation|Liquids|Biomass"] = \
    # var["Final Energy|Transportation|Liquids|Synthetic Fossil"] = \
    var["Final Energy|Transportation|Liquids|Petroleum"] = (
        var["Final Energy|Transportation|Liquids"] * oil_fractions["Fossil"]
    )

    var["Final Energy|Transportation|Liquids|Efuel"] = (
        var["Final Energy|Transportation|Liquids"] * oil_fractions["Efuel"]
        + var["Final Energy|Transportation|Methanol"]
    )

    var["Final Energy|Transportation|Liquids|Biomass"] = (
        var["Final Energy|Transportation|Liquids"] * oil_fractions["Biomass"]
    )

    var["Final Energy|Bunkers|Aviation"] = var[
        "Final Energy|Bunkers|Aviation|Liquids"
    ] = sum_load(n, "kerosene for aviation", region) * international_aviation_fraction

    for var_key, fraction_key in zip(
        ["Biomass", "Petroleum", "Efuel"], oil_fractions.index
    ):
        var[f"Final Energy|Bunkers|Aviation|Liquids|{var_key}"] = (
            var["Final Energy|Bunkers|Aviation|Liquids"] * oil_fractions[fraction_key]
        )
    # TODO Navigation hydrogen
    var["Final Energy|Bunkers|Navigation|Methanol"] = (
        sum_load(n, "shipping methanol", region) * international_navigation_fraction
    )
    var["Final Energy|Bunkers|Navigation|Liquids"] = (
        sum_load(n, "shipping oil", region) * international_navigation_fraction
    )

    for var_key, fraction_key in zip(
        ["Biomass", "Petroleum", "Efuel"], oil_fractions.index
    ):
        var[f"Final Energy|Bunkers|Navigation|Liquids|{var_key}"] = (
            var["Final Energy|Bunkers|Navigation|Liquids"] * oil_fractions[fraction_key]
        )

    var["Final Energy|Bunkers|Navigation|Liquids"] += var[
        "Final Energy|Bunkers|Navigation|Methanol"
    ]
    var["Final Energy|Bunkers|Navigation|Liquids|Efuel"] += var[
        "Final Energy|Bunkers|Navigation|Methanol"
    ]

    var["Final Energy|Bunkers|Navigation"] = var[
        "Final Energy|Bunkers|Navigation|Liquids"
    ]

    # var["Final Energy|Bunkers|Navigation|Gases"] = \
    # ! Not implemented
    # var["Final Energy|Bunkers|Navigation|Hydrogen"] = \
    # ! Not used

    var["Final Energy|Bunkers|Liquids"] = (
        var["Final Energy|Bunkers|Navigation|Liquids"]
        + var["Final Energy|Bunkers|Aviation|Liquids"]
    )
    var["Final Energy|Bunkers"] = (
        var["Final Energy|Bunkers|Navigation"] + var["Final Energy|Bunkers|Aviation"]
    )

    var["Final Energy|Transportation"] = (
        var["Final Energy|Transportation|Electricity"]
        + var["Final Energy|Transportation|Liquids"]
        + var["Final Energy|Transportation|Hydrogen"]
    )

    var["Final Energy|Agriculture|Electricity"] = sum_load(
        n, "agriculture electricity", region
    ) + sum_load(n, "agriculture machinery electric", region)
    var["Final Energy|Agriculture|Heat"] = sum_load(n, "agriculture heat", region)
    var["Final Energy|Agriculture|Liquids"] = sum_load(
        n, "agriculture machinery oil", region
    )
    # var["Final Energy|Agriculture|Gases"] = \
    var["Final Energy|Agriculture"] = (
        var["Final Energy|Agriculture|Electricity"]
        + var["Final Energy|Agriculture|Heat"]
        + var["Final Energy|Agriculture|Liquids"]
    )

    # assert isclose(
    #     var["Final Energy|Agriculture"],
    #     energy_totals.get("total agriculture")
    # )
    # It's nice to do these double checks, but it's less
    # straightforward for the other categories
    # !!! TODO this assert is temporarily disabled because of https://github.com/PyPSA/pypsa-eur/issues/985

    central_heat_withdrawal = (
        n.statistics.withdrawal(
            bus_carrier=["urban central heat"],
            **kwargs,
        )
        .filter(
            like=region,
        )
        .groupby("carrier")
        .sum()
    )

    var["Final Energy|Carbon Dioxide Removal|Heat"] = decentral_heat_withdrawal.get(
        "DAC", 0
    ) + central_heat_withdrawal.get("DAC", 0)

    electricity = (
        n.statistics.withdrawal(
            bus_carrier="AC",
            **kwargs,
        )
        .filter(
            like=region,
        )
        .groupby("carrier")
        .sum()
    )

    var["Final Energy|Carbon Dioxide Removal|Electricity"] = electricity.get("DAC", 0)

    var["Final Energy|Carbon Dioxide Removal"] = (
        var["Final Energy|Carbon Dioxide Removal|Electricity"]
        + var["Final Energy|Carbon Dioxide Removal|Heat"]
    )

    var["Final Energy|Electricity"] = (
        var["Final Energy|Agriculture|Electricity"]
        + var["Final Energy|Residential and Commercial|Electricity"]
        + var["Final Energy|Transportation|Electricity"]
        + var["Final Energy|Industry excl Non-Energy Use|Electricity"]
        + var["Final Energy|Carbon Dioxide Removal|Electricity"]
    )

    var["Final Energy|Solids"] = (
        # var["Final Energy|Agriculture|Solids"]
        var["Final Energy|Residential and Commercial|Solids"]
        + var["Final Energy|Industry excl Non-Energy Use|Solids"]
    )

    var["Final Energy|Solids|Biomass"] = (
        var["Final Energy|Residential and Commercial|Solids|Biomass"]
        + var["Final Energy|Industry excl Non-Energy Use|Solids|Biomass"]
    )

    var["Final Energy|Solids|Coal"] = var[
        "Final Energy|Industry excl Non-Energy Use|Solids|Coal"
    ]

    var["Final Energy|Gases"] = (
        var["Final Energy|Residential and Commercial|Gases"]
        + var["Final Energy|Industry excl Non-Energy Use|Gases"]
    )
    for gas_type in gas_fractions.index:
        var[f"Final Energy|Gases|{gas_type}"] = (
            var["Final Energy|Gases"] * gas_fractions[gas_type]
        )

    var["Final Energy|Liquids"] = (
        var["Final Energy|Agriculture|Liquids"]
        + var["Final Energy|Residential and Commercial|Liquids"]
        + var["Final Energy|Transportation|Liquids"]
        + var["Final Energy|Industry excl Non-Energy Use|Liquids"]
    )

    var["Final Energy|Liquids|Petroleum"] = (
        var["Final Energy|Liquids"] * oil_fractions["Fossil"]
    )

    var["Final Energy|Liquids|Efuel"] = (
        var["Final Energy|Liquids"] * oil_fractions["Efuel"]
    )

    var["Final Energy|Liquids|Biomass"] = (
        var["Final Energy|Liquids"] * oil_fractions["Biomass"]
    )

    var["Final Energy|Heat"] = (
        var["Final Energy|Agriculture|Heat"]
        + var["Final Energy|Residential and Commercial|Heat"]
        + var["Final Energy|Industry excl Non-Energy Use|Heat"]
        + var["Final Energy|Carbon Dioxide Removal|Heat"]
    )
    # var["Final Energy|Solar"] = \
    var["Final Energy|Hydrogen"] = (
        var["Final Energy|Transportation|Hydrogen"]
        + var["Final Energy|Industry excl Non-Energy Use|Hydrogen"]
    )

    # var["Final Energy|Geothermal"] = \
    # ! Not implemented

    waste_withdrawal = (
        n.statistics.withdrawal(
            bus_carrier=["non-sequestered HVC"],
            **kwargs,
        )
        .filter(
            like=region,
        )
        .groupby("carrier")
        .sum()
    )

    var["Final Energy|Waste"] = waste_withdrawal.get("HVC to air", 0)

    var["Final Energy incl Non-Energy Use incl Bunkers"] = (
        var["Final Energy|Industry"]
        + var["Final Energy|Residential and Commercial"]
        + var["Final Energy|Agriculture"]
        + var["Final Energy|Transportation"]
        + var["Final Energy|Bunkers"]
        + var["Final Energy|Waste"]
        + var["Final Energy|Carbon Dioxide Removal"]
    )

    var["Final Energy"] = (
        var["Final Energy|Industry excl Non-Energy Use"]
        + var["Final Energy|Residential and Commercial"]
        + var["Final Energy|Agriculture"]
        + var["Final Energy|Transportation"]
        + var["Final Energy|Waste"]
        + var["Final Energy|Carbon Dioxide Removal"]
    )

    return var * MWh2TWh


def get_emissions(n, region, _energy_totals, _industry_demand):
    energy_totals = _energy_totals.loc[region[0:2]]

    industry_demand = _industry_demand.filter(
        like=region,
        axis=0,
    ).sum()

    kwargs = {
        "groupby": ["name", "bus", "carrier"],
        "nice_names": False,
    }

    var = pd.Series()

    co2_emissions = (
        n.statistics.supply(bus_carrier="co2", **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(t2Mt)
    )

    co2_atmosphere_withdrawal = (
        n.statistics.withdrawal(bus_carrier="co2", **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(t2Mt)
    )

    var["Emissions|CO2|Model"] = co2_emissions.sum() - co2_atmosphere_withdrawal.sum()

    co2_storage = (
        n.statistics.supply(bus_carrier="co2 stored", **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(t2Mt)
    )

    # Assert negligible numerical errors / leakage in stored CO2
    assert co2_storage.get("co2 stored", 0) < 1.0
    co2_storage.drop("co2 stored", inplace=True, errors="ignore")

    try:
        total_ccs = (
            n.statistics.supply(bus_carrier="co2 sequestered", **kwargs)
            .filter(like=region)
            .get("Link")
            .groupby("carrier")
            .sum()
            .multiply(t2Mt)
            .sum()
        )
    except AttributeError:  # no sequestration in 2020 -> NoneType
        total_ccs = 0.0

    # Correcting the PyPSA emission balance
    # Supply side
    # 1. DACCS and BECCS is regarded as negative emissions
    # 2a. DACCU and BECCU are regarded as neutral
    # 2b. BE and fossilCCS are regarded as neutral
    # 3. fossilCCU is regarded as emissions
    # Demand side
    # 4. CCU (e-fuels) is neutral

    ccs_fraction = total_ccs / co2_storage.sum()
    ccu_fraction = 1 - ccs_fraction

    # Correcting for fossil CCU (3.)
    # Step 1: Add the stored CO2 back to the emissions
    # Step 2 (below): CCU goes to e-fuels -> subtract from e-fuel users
    common_index_emitters = co2_emissions.index.intersection(co2_storage.index)

    co2_emissions.loc[common_index_emitters] += co2_storage.loc[
        common_index_emitters
    ].multiply(ccu_fraction)

    var["Emissions|CO2|Model|CCU"] = (
        co2_storage.loc[common_index_emitters].multiply(ccu_fraction).sum()
    )

    # Correcting for BECCS and DACCS (1.)
    # Only stored and sequestered emissions are accounted as negative
    common_index_withdrawals = co2_atmosphere_withdrawal.index.intersection(
        co2_storage.index
    )
    co2_negative_emissions = co2_storage.loc[common_index_withdrawals].multiply(
        ccs_fraction
    )

    # Now repeat the same for the CHP emissions

    CHP_emissions = (
        (
            n.statistics.supply(bus_carrier="co2", **kwargs)
            .filter(like=region)
            .filter(like="CHP")
            .multiply(t2Mt)
        )
        .groupby(["name", "carrier"])
        .sum()
    )

    CHP_atmosphere_withdrawal = (
        (
            n.statistics.withdrawal(bus_carrier="co2", **kwargs)
            .filter(like=region)
            .filter(like="CHP")
            .multiply(t2Mt)
        )
        .groupby(["name", "carrier"])
        .sum()
    )

    CHP_storage = (
        (
            n.statistics.supply(bus_carrier="co2 stored", **kwargs)
            .filter(like=region)
            .filter(like="CHP")
            .multiply(t2Mt)
        )
        .groupby(["name", "carrier"])
        .sum()
    )

    # Addresses 3.

    common_index_emitters = CHP_emissions.index.intersection(CHP_storage.index)

    CHP_emissions.loc[common_index_emitters] += CHP_storage.loc[
        common_index_emitters
    ].multiply(ccu_fraction)

    # Addresses 1.

    common_index_withdrawals = CHP_atmosphere_withdrawal.index.intersection(
        CHP_storage.index
    )

    CHP_negative_emissions = CHP_storage.loc[common_index_withdrawals].multiply(
        ccs_fraction
    )

    ## E-fuels and Biofuels are assumed to be carbon neutral

    oil_techs = [
        "HVC to air",
        "agriculture machinery oil",
        "kerosene for aviation",
        "land transport oil",
        "oil",
        "shipping oil",
        "waste CHP",
        "waste CHP CC",
        "urban decentral oil boiler",
        "rural oil boiler",
        "urban central oil CHP",
    ]

    # multiply by fossil fraction to disregard e- and biofuel emissions

    oil_fossil_fraction = _get_fuel_fractions(n, region, "oil")["Fossil"]

    # This variable is not in the database, but it might be useful for double checking the totals
    var["Emissions|CO2|Efuels|Liquids"] = co2_emissions.loc[
        co2_emissions.index.isin(oil_techs)
    ].sum() * (1 - oil_fossil_fraction)

    co2_emissions.loc[co2_emissions.index.isin(oil_techs)] *= oil_fossil_fraction

    CHP_emissions.loc[
        CHP_emissions.index.get_level_values("carrier").isin(oil_techs)
    ] *= oil_fossil_fraction

    gas_techs = [
        "CCGT",
        "OCGT",
        "SMR",
        "SMR CC",
        "gas for industry",
        "gas for industry CC",
        "rural gas boiler",
        "urban decentral gas boiler",
        "urban central gas CHP",
        "urban central gas CHP CC",
    ]

    gas_fractions = _get_fuel_fractions(n, region, "gas")

    var["Emissions|CO2|Efuels|Gases"] = co2_emissions.loc[
        co2_emissions.index.isin(gas_techs)
    ].sum() * (1 - gas_fractions["Natural Gas"])

    co2_emissions.loc[co2_emissions.index.isin(gas_techs)] *= gas_fractions[
        "Natural Gas"
    ]
    CHP_emissions.loc[
        CHP_emissions.index.get_level_values("carrier").isin(gas_techs)
    ] *= gas_fractions["Natural Gas"]

    # All methanol is e-methanol and hence considered carbon neutral

    methanol_techs = [
        "industry methanol",
        "shipping methanol",
    ]

    var["Emissions|CO2|Efuels|Methanol"] = co2_emissions.loc[
        co2_emissions.index.isin(methanol_techs)
    ].sum()

    co2_emissions.loc[co2_emissions.index.isin(methanol_techs)] = 0

    # Emissions in DE are:

    var["Emissions|CO2"] = co2_emissions.sum() - co2_negative_emissions.sum()

    # I am not sure any longer how to exactly capture the difference between Model emissions and the reported emissions
    # assert isclose(
    #     var["Emissions|CO2"],
    #     var["Emissions|CO2|Model"]
    #     - var["Emissions|CO2|Efuels|Liquids"]
    #     - var["Emissions|CO2|Efuels|Gases"]
    #     - var["Emissions|CO2|Efuels|Methanol"]
    #     - var["Emissions|CO2|Model|CCU"],
    #     rtol=1e-2, atol=1e-2,
    # )

    # Split CHP emissions between electricity and heat sectors

    CHP_E_to_H = (
        n.links.loc[CHP_emissions.index.get_level_values("name")].efficiency
        / n.links.loc[CHP_emissions.index.get_level_values("name")].efficiency2
    )

    CHP_E_fraction = CHP_E_to_H * (1 / (CHP_E_to_H + 1))

    negative_CHP_E_to_H = (
        n.links.loc[CHP_atmosphere_withdrawal.index.get_level_values("name")].efficiency
        / n.links.loc[
            CHP_atmosphere_withdrawal.index.get_level_values("name")
        ].efficiency2
    )

    negative_CHP_E_fraction = negative_CHP_E_to_H * (1 / (negative_CHP_E_to_H + 1))

    # It would be interesting to relate the Emissions|CO2|Model to Emissions|CO2 reported to the DB by considering imports of carbon, e.g., (exports_oil_renew - imports_oil_renew) * 0.2571 * t2Mt + (exports_gas_renew - imports_gas_renew) * 0.2571 * t2Mt + (exports_meoh - imports_meoh) / 4.0321 * t2Mt
    # Then it would be necessary to consider negative carbon from solid biomass imports as well
    # Actually we might have to include solid biomass imports in the co2 constraints as well

    assert isclose(
        co2_emissions.filter(like="CHP").sum(),
        CHP_emissions.sum(),
    )
    assert isclose(
        co2_atmosphere_withdrawal.filter(like="CHP").sum(),
        CHP_atmosphere_withdrawal.sum(),
    )

    process_emissions = (
        n.statistics.supply(bus_carrier="process emissions", **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
    )

    pe_fossil_fraction = (
        process_emissions.get("process emissions", 0)
        + process_emissions.get("naptha for industry", 0) * oil_fossil_fraction
    ) / process_emissions.sum()

    var["Carbon Sequestration|DACCS"] = co2_negative_emissions.get("DAC", 0)

    var["Carbon Sequestration|BECCS"] = co2_negative_emissions.filter(like="bio").sum()

    # E and Biofuels with CC
    var["Carbon Sequestration|Other"] = co2_storage.mul(ccs_fraction)[
        ~co2_storage.index.str.contains("bio|process")
    ].sum() + co2_storage.mul(ccs_fraction).get("process emissions CC", 0) * (
        1 - pe_fossil_fraction
    )

    var["Carbon Sequestration"] = (
        var["Carbon Sequestration|DACCS"]
        + var["Carbon Sequestration|BECCS"]
        + var["Carbon Sequestration|Other"]
    )

    # assert isclose(
    #     var["Carbon Sequestration"],
    #     co2_storage.mul(ccs_fraction)[~co2_storage.index.str.contains("process")].sum(),
    # )

    # ! LULUCF should also be subtracted (or added??), we get from REMIND,
    # TODO how to consider it here?

    # Make sure these values are about right
    var["Emissions|CO2|Industrial Processes"] = co2_emissions.reindex(
        [
            "process emissions",
            "process emissions CC",
        ]
    ).sum() * pe_fossil_fraction + co2_emissions.get(
        "industry methanol", 0
    )  # considered 0 anyways

    mwh_coal_per_mwh_coke = 1.366  # from eurostat energy balance

    coke_fraction = (
        industry_demand.get("coke")
        * mwh_coal_per_mwh_coke
        / (
            industry_demand.get("coke") * mwh_coal_per_mwh_coke
            + industry_demand.get("coal")
        )
    )
    # 0.3361 t_CO2/MWh
    coking_emissions = (
        sum_load(n, "coal for industry", region)
        * coke_fraction
        * (mwh_coal_per_mwh_coke - 1)
        * 0.3361
        * t2Mt
    )
    var["Emissions|Gross Fossil CO2|Energy|Demand|Industry"] = (
        co2_emissions.reindex(
            [
                "gas for industry",
                "gas for industry CC",
                "coal for industry",
            ]
        ).sum()
    ) - coking_emissions
    var["Emissions|Gross Fossil CO2|Energy|Supply|Solids"] = var[
        "Emissions|CO2|Energy|Supply|Solids"
    ] = coking_emissions

    var["Emissions|CO2|Energy|Demand|Industry"] = var[
        "Emissions|Gross Fossil CO2|Energy|Demand|Industry"
    ] - co2_negative_emissions.get("solid biomass for industry CC", 0)

    var["Emissions|CO2|Industry"] = (
        var["Emissions|CO2|Energy|Demand|Industry"]
        + var["Emissions|CO2|Industrial Processes"]
    )

    var["Emissions|CO2|Energy|Demand|Residential and Commercial"] = (
        co2_emissions.filter(like="urban decentral").sum()
        + co2_emissions.filter(like="rural").sum()
    )

    international_aviation_fraction = energy_totals["total international aviation"] / (
        energy_totals["total domestic aviation"]
        + energy_totals["total international aviation"]
    )
    international_navigation_fraction = energy_totals[
        "total international navigation"
    ] / (
        energy_totals["total domestic navigation"]
        + energy_totals["total international navigation"]
    )

    var["Emissions|CO2|Energy|Demand|Transportation|Domestic Aviation"] = (
        co2_emissions.get("kerosene for aviation")
        * (1 - international_aviation_fraction)
    )

    var["Emissions|CO2|Energy|Demand|Transportation|Domestic Navigation"] = (
        co2_emissions.filter(like="shipping").sum()
        * (1 - international_navigation_fraction)
    )

    var["Emissions|CO2|Energy|Demand|Transportation"] = (
        co2_emissions.get("land transport oil", 0)
        + var["Emissions|CO2|Energy|Demand|Transportation|Domestic Aviation"]
        + var["Emissions|CO2|Energy|Demand|Transportation|Domestic Navigation"]
    )

    var["Emissions|CO2|Energy|Demand|Bunkers|Aviation"] = (
        co2_emissions.get("kerosene for aviation") * international_aviation_fraction
    )

    var["Emissions|CO2|Energy|Demand|Bunkers|Navigation"] = (
        co2_emissions.filter(like="shipping").sum() * international_navigation_fraction
    )

    var["Emissions|CO2|Energy|Demand|Bunkers"] = (
        var["Emissions|CO2|Energy|Demand|Bunkers|Aviation"]
        + var["Emissions|CO2|Energy|Demand|Bunkers|Navigation"]
    )

    var["Emissions|CO2|Energy|Demand|Other Sector"] = co2_emissions.get(
        "agriculture machinery oil"
    )

    var["Emissions|CO2|Energy|Demand"] = var.get(
        [
            "Emissions|CO2|Energy|Demand|Industry",
            "Emissions|CO2|Energy|Demand|Transportation",
            "Emissions|CO2|Energy|Demand|Residential and Commercial",
            "Emissions|CO2|Energy|Demand|Other Sector",
        ]
    ).sum()
    var["Emissions|CO2|Energy incl Bunkers|Demand"] = (
        var["Emissions|CO2|Energy|Demand"] + var["Emissions|CO2|Energy|Demand|Bunkers"]
    )

    CHP_emissions_E = CHP_emissions.multiply(CHP_E_fraction).groupby("carrier").sum()
    CHP_emissions_H = (
        CHP_emissions.multiply(1 - CHP_E_fraction).groupby("carrier").sum()
    )
    CHP_negative_emissions_E = (
        CHP_negative_emissions.multiply(negative_CHP_E_fraction)
        .groupby("carrier")
        .sum()
    )
    CHP_negative_emissions_H = (
        CHP_negative_emissions.multiply(1 - negative_CHP_E_fraction)
        .groupby("carrier")
        .sum()
    )

    var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Gas"] = (
        co2_emissions.reindex(
            [
                "OCGT",
                "CCGT",
            ],
        ).sum()
        + CHP_emissions_E.filter(like="gas").sum()
    )
    var["Emissions|CO2|Energy|Supply|Electricity|Gas"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Gas"]
        - CHP_negative_emissions_E.filter(like="gas").sum()
    )

    var["Emissions|CO2|Energy|Supply|Electricity|Coal"] = var[
        "Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Coal"
    ] = (
        co2_emissions.reindex(
            [
                "coal",
                "lignite",
            ],
        ).sum()
        + CHP_emissions_E.filter(regex="coal|lignite").sum()
    )
    var["Emissions|CO2|Energy|Supply|Electricity|Oil"] = var[
        "Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Oil"
    ] = (
        co2_emissions.reindex(
            [
                "oil",
            ],
        ).sum()
        + CHP_emissions_E.filter(like="oil").sum()
    )

    var["Emissions|CO2|Energy|Supply|Electricity|Biomass"] = (
        CHP_negative_emissions_E.filter(like="bio").sum()
    )

    var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Waste"] = (
        CHP_emissions_E.filter(like="waste CHP").sum()
    )
    var["Emissions|CO2|Energy|Supply|Electricity|Waste"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Waste"]
        - CHP_negative_emissions_E.filter(like="waste").sum()
    )

    var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Gas"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Coal"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Oil"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity|Waste"]
    )

    assert isclose(
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"],
        co2_emissions.reindex(
            [
                "OCGT",
                "CCGT",
                "coal",
                "lignite",
                "oil",
            ],
        ).sum()
        + CHP_emissions_E.sum(),
    )

    var["Emissions|CO2|Energy|Supply|Electricity"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"]
        - CHP_negative_emissions_E.sum()
    )

    var["Emissions|Gross Fossil CO2|Energy|Supply|Heat"] = (
        co2_emissions.filter(like="urban central").filter(like="boiler").sum()
        + CHP_emissions_H.sum()
    )

    var["Emissions|CO2|Energy|Supply|Heat"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Heat"]
        - CHP_negative_emissions_H.sum()
    )

    var["Emissions|CO2|Energy|Supply|Electricity and Heat"] = (
        var["Emissions|CO2|Energy|Supply|Heat"]
        + var["Emissions|CO2|Energy|Supply|Electricity"]
    )

    var["Emissions|CO2|Energy|Supply|Hydrogen"] = var[
        "Emissions|Gross Fossil CO2|Energy|Supply|Hydrogen"
    ] = co2_emissions.filter(like="SMR").sum()

    var["Emissions|Gross Fossil CO2|Energy|Supply|Gases"] = co2_emissions.get(
        "gas compressing", 0
    )

    var["Emissions|CO2|Energy|Supply|Gases"] = (-1) * co2_negative_emissions.get(
        "biogas to gas CC", 0
    ) + var["Emissions|Gross Fossil CO2|Energy|Supply|Gases"]

    var["Emissions|CO2|Supply|Non-Renewable Waste"] = co2_emissions.get("HVC to air", 0)

    var["Emissions|Gross Fossil CO2|Energy|Supply|Liquids"] = co2_emissions.get(
        "oil refining", 0
    )

    var["Emissions|CO2|Energy|Supply|Liquids"] = var[
        "Emissions|Gross Fossil CO2|Energy|Supply|Liquids"
    ] - co2_negative_emissions.get("biomass to liquid CC", 0)

    var["Emissions|CO2|Energy|Supply|Liquids and Gases"] = (
        var["Emissions|CO2|Energy|Supply|Liquids"]
        + var["Emissions|CO2|Energy|Supply|Gases"]
    )

    var["Emissions|Gross Fossil CO2|Energy|Supply"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Heat"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Hydrogen"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Liquids"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Solids"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Gases"]
    )

    var["Emissions|Gross Fossil CO2|Energy"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply"]
        + var["Emissions|Gross Fossil CO2|Energy|Demand|Industry"]
    )

    var["Emissions|CO2|Energy|Supply"] = (
        var["Emissions|CO2|Energy|Supply|Gases"]
        + var["Emissions|CO2|Energy|Supply|Hydrogen"]
        + var["Emissions|CO2|Energy|Supply|Electricity and Heat"]
        + var["Emissions|CO2|Energy|Supply|Liquids"]
        + var["Emissions|CO2|Energy|Supply|Solids"]
    )

    # var["Emissions|CO2|Energy|Supply|Other Sector"] = \

    var["Emissions|CO2|Energy"] = (
        var["Emissions|CO2|Energy|Demand"] + var["Emissions|CO2|Energy|Supply"]
    )

    var["Emissions|CO2|Energy incl Bunkers"] = (
        var["Emissions|CO2|Energy incl Bunkers|Demand"]
        + var["Emissions|CO2|Energy|Supply"]
    )

    var["Emissions|CO2|Energy and Industrial Processes"] = (
        var["Emissions|CO2|Energy"] + var["Emissions|CO2|Industrial Processes"]
    )

    emission_difference = var["Emissions|CO2"] - (
        var["Emissions|CO2|Energy and Industrial Processes"]
        + var["Emissions|CO2|Energy|Demand|Bunkers"]
        + var["Emissions|CO2|Supply|Non-Renewable Waste"]
        - co2_negative_emissions.get("DAC", 0)
    )

    logger.info(
        "Differences in accounting for CO2 emissions: %s",
        emission_difference,
    )

    # assert abs(emission_difference) < 1e-5

    return var


# functions for prices
def get_nodal_flows(n, bus_carrier, region, query="index == index or index != index"):
    """
    Get the nodal flows for a given bus carrier and region.

    Parameters
    ----------
        n (pypsa.Network): The PyPSA network object.
        bus_carrier (str): The bus carrier for which to retrieve the nodal flows.
        region (str): The region for which to retrieve the nodal flows.
        query (str, optional): A query string to filter the nodal flows. Defaults to 'index == index or index != index'.

    Returns
    -------
        pandas.DataFrame: The nodal flows for the specified bus carrier and region.
    """

    groupby = ["name", "bus", "carrier"]

    result = (
        n.statistics.withdrawal(
            bus_carrier=bus_carrier,
            groupby=groupby,
            groupby_time=False,
        )
        .query(query)
        .groupby("bus")
        .sum()
        .T.filter(
            like=region,
            axis=1,
        )
    )

    return result


def get_nodal_supply(n, bus_carrier, query="index == index or index != index"):
    """
    Get the nodal flows for a given bus carrier and region.

    Parameters
    ----------
        n (pypsa.Network): The PyPSA network object.
        bus_carrier (str): The bus carrier for which to retrieve the nodal flows.
        region (str): The region for which to retrieve the nodal flows.
        query (str, optional): A query string to filter the nodal flows. Defaults to 'index == index or index != index'.

    Returns
    -------
        pandas.DataFrame: The nodal flows for the specified bus carrier and region.
    """

    groupby = ["name", "bus", "carrier"]

    result = (
        n.statistics.supply(
            bus_carrier=bus_carrier,
            groupby=groupby,
            groupby_time=False,
        )
        .query(query)
        .groupby("bus")
        .sum()
        .T
    )

    return result


def price_load(n, load_carrier, region):
    """
    Calculate the average price of a specific load carrier in a given region.

    Parameters
    ----------
    - n (pandas.DataFrame): The network model.
    - load_carrier (str): The load carrier to calculate the price for.
    - region (str): The region to calculate the price in.

    Returns
    -------
    - tuple: A tuple containing the average price and the total load of the specified load carrier in the region.
    """

    load = n.loads[
        (n.loads.carrier == load_carrier) & (n.loads.bus.str.contains(region))
    ]
    if n.loads_t.p[load.index].values.sum() < 1:
        return np.nan, 0
    result = (
        n.loads_t.p[load.index] * n.buses_t.marginal_price[load.bus].values
    ).values.sum()
    result /= n.loads_t.p[load.index].values.sum()
    return result, n.loads_t.p[load.index].values.sum()


def costs_gen_generators(n, region, carrier):
    """
    Calculate the cost per unit of generated energy of a generators in a given
    region.

    Parameters
    ----------
    - n (pandas.DataFrame): The network model.
    - region (str): The region to consider.
    - carrier (str): The carrier of the generators.

    Returns
    -------
    - tuple: A tuple containing cost and total generation of the generators.
    """

    gens = n.generators[
        (n.generators.carrier == carrier) & (n.generators.bus.str.contains(region))
    ]
    gen = (
        n.generators_t.p[gens.index]
        .multiply(n.snapshot_weightings.generators, axis="index")
        .sum()
    )
    if gen.empty or gen.sum() < 1:
        return np.nan, 0

    # CAPEX
    capex = (gens.p_nom_opt * gens.capital_cost).sum()

    # OPEX
    opex = (gen * gens.marginal_cost).sum()

    result = (capex + opex) / gen.sum()
    return result, gen.sum()


def costs_gen_links(n, region, carrier, gen_bus="p1"):
    """
    Calculate the cost per unit of generated energy from a specific link.

    Parameters
    ----------
        n (pypsa.Network): The PyPSA network object.
        region (str): The region to consider for the links.
        carrier (str): The carrier of the links.
        gen_bus (str, optional): The bus where the main generation of the link takes place. Defaults to "p1".

    Returns
    -------
        tuple: A tuple containing the costs per unit of generetad energy and the total generation of the specified generator bus.
    """

    links = n.links[(n.links.carrier == carrier) & (n.links.index.str.contains(region))]
    gen = abs(
        n.links_t[gen_bus][links.index].multiply(
            n.snapshot_weightings.generators, axis="index"
        )
    ).sum()
    if gen.empty or gen.sum() < 1:
        return np.nan, 0

    # CAPEX
    capex = (links.p_nom_opt * links.capital_cost).sum()

    # OPEX
    input = abs(
        n.links_t["p0"][links.index].multiply(
            n.snapshot_weightings.generators, axis="index"
        )
    ).sum()
    opex = (input * links.marginal_cost).sum()

    # input costs and output revenues other than main generation @ gen_bus
    sum = 0
    for i in range(0, 5):
        if f"p{i}" == gen_bus:
            continue
        elif links.empty:
            break
        elif n.links.loc[links.index][f"bus{i}"].iloc[0] == "":
            break
        else:
            update_cost = (
                (
                    n.links_t[f"p{i}"][links.index]
                    * n.buses_t.marginal_price[links[f"bus{i}"]].values
                )
                .multiply(n.snapshot_weightings.generators, axis="index")
                .values.sum()
            )
            sum = sum + update_cost

    result = (capex + opex + sum) / gen.sum()
    return result, gen.sum()


def get_weighted_costs_links(carriers, n, region):
    numerator = 0
    denominator = 0

    for c in carriers:
        cost_gen = costs_gen_links(n, region, c)
        if not math.isnan(cost_gen[0]):
            numerator += cost_gen[0] * cost_gen[1]
            denominator += cost_gen[1]

    if denominator == 0:
        return np.nan
    result = numerator / denominator
    result = np.where(result < 1e4, result, np.nan)
    return result


def get_weighted_costs(costs, flows):
    cleaned_costs = []
    cleaned_flows = []

    for cost, flow in zip(costs, flows):
        if not math.isnan(cost) and not math.isnan(flow) and flow != 0:
            cleaned_costs.append(cost)
            cleaned_flows.append(flow)

    if not cleaned_costs or not cleaned_flows:
        return np.nan

    df_cleaned = pd.DataFrame({"costs": cleaned_costs, "flows": cleaned_flows})
    result = (df_cleaned["costs"] * df_cleaned["flows"]).sum() / df_cleaned[
        "flows"
    ].sum()
    return result


def get_prices(n, region):
    """
    Calculate the prices of various energy sources in the Ariadne model.

    Parameters
    ----------
    - n (PyPSa network): The Ariadne model scenario output.
    - region (str): The region for which the prices are calculated.

    Returns
    -------
    - var (pandas.Series): A series containing the calculated prices.

    This function calculates the prices of different energy sources in the Ariadne model
    based on the nodal flows and marginal prices of the model. The calculated prices are
    stored in a pandas Series object and returned.
    """

    # last Update according to template from 2024-08-13
    # in Ariadne database: 200
    # currently reported: 32

    var = pd.Series()

    try:
        co2_limit_de = n.global_constraints.loc["co2_limit-DE", "mu"]
    except KeyError:
        co2_limit_de = 0

    # co2 additions
    co2_price = -n.global_constraints.loc["CO2Limit", "mu"] - co2_limit_de
    # specific emissions in tons CO2/MWh according to n.links[n.links.carrier =="your_carrier].efficiency2.unique().item()
    specific_emissions = {
        "oil": 0.2571,
        "gas": 0.198,  # OCGT
        "hard coal": 0.3361,
        "lignite": 0.4069,
    }

    # vars: Tier 1, Category: energy(price)

    # Price|Primary Energy
    # reported: 4/4

    nodal_flows_bm = get_nodal_flows(n, "solid biomass", region)
    nodal_prices_bm = n.buses_t.marginal_price[nodal_flows_bm.columns]

    # Price|Primary Energy|Biomass
    var["Price|Primary Energy|Biomass"] = (
        nodal_flows_bm.mul(nodal_prices_bm).values.sum() / nodal_flows_bm.values.sum()
    )

    # Price|Primary Energy|Coal
    nf_coal = get_nodal_flows(n, "coal", "EU")
    nodal_prices_coal = n.buses_t.marginal_price[nf_coal.columns]
    coal_price = (
        nf_coal.mul(nodal_prices_coal).values.sum() / nf_coal.values.sum()
        if nf_coal.values.sum() > 0
        else np.nan
    )

    nf_lignite = get_nodal_flows(n, "lignite", "EU")
    nodal_prices_lignite = n.buses_t.marginal_price[nf_lignite.columns]
    lignite_price = (
        nf_lignite.mul(nodal_prices_lignite).values.sum() / nf_lignite.values.sum()
        if nf_lignite.values.sum() > 0
        else np.nan
    )

    coal_fraction = nf_coal.values.sum() / (
        nf_coal.values.sum() + nf_lignite.values.sum()
    )
    lignite_fraction = nf_lignite.values.sum() / (
        nf_coal.values.sum() + nf_lignite.values.sum()
    )
    co2_add_coal = (
        coal_fraction * specific_emissions["hard coal"] * co2_price
        + lignite_fraction * specific_emissions["lignite"] * co2_price
    )

    var["Price|Primary Energy|Coal"] = (
        get_weighted_costs(
            [coal_price, lignite_price], [nf_coal.values.sum(), nf_lignite.values.sum()]
        )
        + co2_add_coal
    )

    # Price|Primary Energy|Gas
    nodal_flows_gas = get_nodal_flows(n, "gas", region)
    nodal_prices_gas = n.buses_t.marginal_price[nodal_flows_gas.columns]

    # co2 part
    co2_cost_gas = specific_emissions["gas"] * co2_price

    var["Price|Primary Energy|Gas"] = (
        nodal_flows_gas.mul(nodal_prices_gas).values.sum()
        / nodal_flows_gas.values.sum()
        + co2_cost_gas
    )

    # Price|Primary Energy|Oil
    nodal_flows_oil = get_nodal_flows(n, "oil primary", region)
    nodal_prices_oil = n.buses_t.marginal_price[nodal_flows_oil.columns]

    # co2 part
    co2_cost_oil = specific_emissions["oil"] * co2_price

    var["Price|Primary Energy|Oil"] = (
        nodal_flows_oil.mul(nodal_prices_oil).values.sum()
        / nodal_flows_oil.values.sum()
        + co2_cost_oil
    )

    # Price|Secondary Energy
    # reported: 8/14

    # Price|Secondary Energy|Electricity
    nodal_flows_ac = get_nodal_flows(
        n, "AC", region, query="not carrier.str.contains('gas')"
    )
    nodal_prices_ac = n.buses_t.marginal_price[nodal_flows_ac.columns]

    var["Price|Secondary Energy|Electricity"] = (
        nodal_flows_ac.mul(nodal_prices_ac).values.sum() / nodal_flows_ac.values.sum()
    )

    # Price|Secondary Energy|Gases|Natural Gas
    var["Price|Secondary Energy|Gases|Natural Gas"] = (
        costs_gen_generators(n, region, "gas")[0] + co2_cost_gas
    )

    # Price|Secondary Energy|Gases|Hydrogen
    var["Price|Secondary Energy|Gases|Hydrogen"] = costs_gen_links(
        n, region, "Sabatier"
    )[0]

    # Price|Secondary Energy|Gases|Biomass
    var["Price|Secondary Energy|Gases|Biomass"] = get_weighted_costs_links(
        ["biogas to gas", "biogas to gas CC"], n, region
    )

    # Price|Secondary Energy|Gases|Efuel
    # Price for gaseous Efuels at the secondary level, i.e. for large scale consumers. Prices should include the effect of carbon prices.
    # what are gaseous Efuels?

    # Price|Secondary Energy|Gases
    nodal_flows_gas = get_nodal_flows(
        n,
        ["gas", "renewable gas"],
        region,
        query="not carrier.str.contains('pipeline')"
        "& not carrier == 'gas'"
        "& not carrier.str.contains('rural')"
        "& not carrier.str.contains('urban decentral')",
    )
    nodal_prices_gas = n.buses_t.marginal_price[nodal_flows_gas.columns]
    nodal_prices_gas.loc[:, "DE gas"] = nodal_prices_gas["DE gas"] + co2_cost_gas

    var["Price|Secondary Energy|Gases"] = (
        nodal_flows_gas.mul(nodal_prices_gas).values.sum()
        / nodal_flows_gas.values.sum()
    )

    # Price|Secondary Energy|Hydrogen
    # (carbon costs not yet included)
    # nodal_flows_h2 = get_nodal_flows(n, "H2", region)
    # nodal_prices_h2 = n.buses_t.marginal_price[nodal_flows_h2.columns]

    # var["Price|Secondary Energy|Hydrogen"] = (
    #     nodal_flows_h2.mul(nodal_prices_h2).values.sum() / nodal_flows_h2.values.sum()
    # )

    carriers = ["SMR", "SMR CC", "H2 Electrolysis"]
    var["Price|Secondary Energy|Hydrogen"] = get_weighted_costs_links(
        carriers, n, region
    )

    # Price|Secondary Energy|Hydrogen|Green
    # incorporate CAPEX and electricity price as electrolysis is forced in in some years
    var["Price|Secondary Energy|Hydrogen|Green"] = costs_gen_links(
        n, region, "H2 Electrolysis"
    )[0]

    # Price|Secondary Energy|Liquids
    nodal_flows_oil = get_nodal_flows(
        n,
        ["oil primary", "renewable oil"],
        "DE",
        query="not carrier.str.contains('rural')"
        "& not carrier.str.contains('urban decentral')",
    )
    nodal_prices_oil = n.buses_t.marginal_price[nodal_flows_oil.columns]
    nodal_prices_oil.loc[:, "DE oil primary"] = (
        nodal_prices_oil["DE oil primary"] + co2_cost_oil
    )

    var["Price|Secondary Energy|Liquids"] = (
        nodal_flows_oil.mul(nodal_prices_oil).values.sum()
        / nodal_flows_oil.values.sum()
    )

    # Price|Secondary Energy|Liquids|Biomass
    # Price|Secondary Energy|Liquids|Oil
    # Price|Secondary Energy|Liquids|Hydrogen
    carriers = ["Fischer-Tropsch", "methanol"]
    var["Price|Secondary Energy|Liquids|Hydrogen"] = get_weighted_costs_links(
        carriers, n, region
    )

    # Price|Secondary Energy|Liquids|Efuel
    # Price|Secondary Energy|Solids|Biomass
    # Price|Secondary Energy|Solids|Coal

    # "Price|Final Energy|
    # reported: 23/182
    # warning: these prices do not incorporate co2 prices yet

    ### Price|Final Energy|Transportation|

    # Price|Final Energy|Transportation|Freight|Electricity x
    # Price|Final Energy|Transportation|Freight|Gases
    # Price|Final Energy|Transportation|Freight|Hydrogen
    var["Price|Final Energy|Transportation|Freight|Hydrogen"] = price_load(
        n, "land transport fuel cell", region
    )[0]
    # Price|Final Energy|Transportation|Freight|Liquids
    carriers = [
        "kerosene for aviation",
        "shipping methanol",
        "shipping oil",
        "land transport oil",
    ]
    df = pd.DataFrame({c: price_load(n, c, region) for c in carriers})
    var["Price|Final Energy|Transportation|Freight|Liquids"] = (
        df.iloc[0] * df.iloc[1]
    ).sum() / df.iloc[1].sum()
    # Price|Final Energy|Transportation|Freight|Solids x

    # Price|Final Energy|Transportation|Passenger|Electricity
    var["Price|Final Energy|Transportation|Passenger|Electricity"] = price_load(
        n, "land transport EV", region
    )[0]
    # Price|Final Energy|Transportation|Passenger|Gases
    # Price|Final Energy|Transportation|Passenger|Hydrogen
    var["Price|Final Energy|Transportation|Passenger|Hydrogen"] = price_load(
        n, "land transport fuel cell", region
    )[0]
    # Price|Final Energy|Transportation|Passenger|Liquids
    carriers = ["kerosene for aviation", "land transport oil"]
    df = pd.DataFrame({c: price_load(n, c, region) for c in carriers})
    var["Price|Final Energy|Transportation|Passenger|Liquids"] = (
        df.iloc[0] * df.iloc[1]
    ).sum() / df.iloc[1].sum()
    # Price|Final Energy|Transportation|Passenger|Solids x

    # Price|Final Energy|Transportation|Liquids|Petroleum
    # Price|Final Energy|Transportation|Liquids|Petroleum|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Petroleum|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Petroleum|Carbon Price Component
    # Price|Final Energy|Transportation|Liquids|Petroleum|Other Taxes
    # Price|Final Energy|Transportation|Liquids|Diesel
    # Price|Final Energy|Transportation|Liquids|Diesel|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Diesel|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Diesel|Carbon Price Component
    # Price|Final Energy|Transportation|Liquids|Diesel|Other Taxes

    # Price|Final Energy|Transportation|Gases|Natural Gas
    # Price|Final Energy|Transportation|Gases|Natural Gas|Sales Margin
    # Price|Final Energy|Transportation|Gases|Natural Gas|Transport and Distribution
    # Price|Final Energy|Transportation|Gases|Natural Gas|Carbon Price Component
    # Price|Final Energy|Transportation|Gases|Natural Gas|Other Taxes

    # Price|Final Energy|Transportation|Liquids|Biomass
    # Price|Final Energy|Transportation|Liquids|Biomass|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Biomass|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Biomass|Other Taxes

    # Price|Final Energy|Transportation|Liquids|Efuel
    # Price|Final Energy|Transportation|Liquids|Efuel|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Efuel|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Efuel|Other Taxes

    # Price|Final Energy|Transportation|Gases|Efuel
    # Price|Final Energy|Transportation|Gases|Efuel|Sales Margin
    # Price|Final Energy|Transportation|Gases|Efuel|Transport and Distribution
    # Price|Final Energy|Transportation|Gases|Efuel|Other Taxes

    # Price|Final Energy|Transportation|Gases|SNG
    # Price|Final Energy|Transportation|Gases|SNG|Sales Margin
    # Price|Final Energy|Transportation|Gases|SNG|Transport and Distribution
    # Price|Final Energy|Transportation|Gases|SNG|Other Taxes
    # Price|Final Energy|Transportation|Gases|SNG|Carbon Price Component

    # Price|Final Energy|Transportation|Hydrogen
    # Price|Final Energy|Transportation|Hydrogen|Sales Margin
    # Price|Final Energy|Transportation|Hydrogen|Transport and Distribution
    # Price|Final Energy|Transportation|Hydrogen|Other Taxes

    # Price|Final Energy|Transportation|Electricity
    var["Price|Final Energy|Transportation|Electricity"] = price_load(
        n, "land transport EV", region
    )[0]
    # Price|Final Energy|Transportation|Electricity|Sales Margin
    # Price|Final Energy|Transportation|Electricity|Transport and Distribution
    # Price|Final Energy|Transportation|Electricity|Other Taxes

    # Price|Final Energy|Transportation|Liquids|Kerosene
    var["Price|Final Energy|Transportation|Electricity"] = price_load(
        n, "kerosene for aviation", region
    )[0]
    # Price|Final Energy|Transportation|Liquids|Kerosene|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Kerosene|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Kerosene|Carbon Price Component
    # Price|Final Energy|Transportation|Liquids|Kerosene|Other Taxes

    # Price|Final Energy|Transportation|Electricity|Carbon Price Component ?
    # Price|Final Energy|Transportation|Gases|Carbon Price Component
    # Price|Final Energy|Transportation|Hydrogen|Carbon Price Component
    # Price|Final Energy|Transportation|Liquids|Carbon Price Component

    ### Price|Final Energy|Residential|

    # Price|Final Energy|Residential|Electricity
    # Price|Final Energy|Residential|Hydrogen
    # Price|Final Energy|Residential|Gases
    # Price|Final Energy|Residential|Gases|Natural Gas
    # Price|Final Energy|Residential|Liquids
    # Price|Final Energy|Residential|Liquids|Biomass
    # Price|Final Energy|Residential|Liquids|Oil
    # Price|Final Energy|Residential|Solids
    # Price|Final Energy|Residential|Solids|Biomass
    # Price|Final Energy|Residential|Solids|Coal

    ### Price|Final Energy|Residential and Commercial

    # Price|Final Energy|Residential and Commercial|Liquids|Oil
    carriers = ["rural oil boiler", "urban decentral oil boiler"]
    var["Price|Final Energy|Residential and Commercial|Liquids|Oil"] = (
        get_weighted_costs_links(carriers, n, region)
    )
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Sales Margin
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Carbon Price Component
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Other Taxes
    # Price|Final Energy|Residential and Commercial|Liquids
    var["Price|Final Energy|Residential and Commercial|Liquids"] = var[
        "Price|Final Energy|Residential and Commercial|Liquids|Oil"
    ]

    # Price|Final Energy|Residential and Commercial|Gases|
    carriers = [
        "urban central gas boiler",
        "urban central gas CHP",
        "urban central gas CHP CC",
        "rural gas boiler",
    ]
    var["Price|Final Energy|Residential and Commercial|Gases"] = (
        get_weighted_costs_links(carriers, n, region)
    )
    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas
    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas|Sales Margin
    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas|Carbon Price Component
    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas|Other Taxes

    # Price|Final Energy|Residential and Commercial|Heat
    nf_rc_heat = get_nodal_flows(
        n,
        ["urban central heat", "rural heat", "urban decentral heat"],
        region,
        query="not carrier.str.contains('agriculture')"
        "& not carrier.str.contains('industry')"
        "& not carrier.str.contains('DAC')",
    )
    np_rc_heat = n.buses_t.marginal_price[nf_rc_heat.columns]
    var["Price|Final Energy|Residential and Commercial|Heat"] = (
        nf_rc_heat.mul(np_rc_heat).values.sum() / nf_rc_heat.values.sum()
    )
    # Price|Final Energy|Residential and Commercial|Heat|Sales Margin
    # Price|Final Energy|Residential and Commercial|Heat|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Heat|Other Taxes

    # Price|Final Energy|Residential and Commercial|Liquids|Biomass x
    # Price|Final Energy|Residential and Commercial|Liquids|Biomass|Sales Margin
    # Price|Final Energy|Residential and Commercial|Liquids|Biomass|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Liquids|Biomass|Other Taxes

    # Price|Final Energy|Residential and Commercial|Solids
    carriers = "rural biomass boiler", "urban decentral biomass boiler"
    var["Price|Final Energy|Residential and Commercial|Solids|Biomass"] = (
        get_weighted_costs_links(carriers, n, region)
    )
    var["Price|Final Energy|Residential and Commercial|Solids"] = var[
        "Price|Final Energy|Residential and Commercial|Solids|Biomass"
    ]
    # Price|Final Energy|Residential and Commercial|Solids|Biomass|Sales Margin
    # Price|Final Energy|Residential and Commercial|Solids|Biomass|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Solids|Biomass|Other Taxes

    # Price|Final Energy|Residential and Commercial|Gases|Biomass x
    # Price|Final Energy|Residential and Commercial|Gases|Biomass|Sales Margin
    # Price|Final Energy|Residential and Commercial|Gases|Biomass|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Gases|Biomass|Other Taxes

    # Price|Final Energy|Residential and Commercial|Liquids|Efuel x
    # Price|Final Energy|Residential and Commercial|Liquids|Efuel|Sales Margin
    # Price|Final Energy|Residential and Commercial|Liquids|Efuel|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Liquids|Efuel|Other Taxes

    # Price|Final Energy|Residential and Commercial|Gases|Efuel x
    # Price|Final Energy|Residential and Commercial|Gases|Efuel|Sales Margin
    # Price|Final Energy|Residential and Commercial|Gases|Efuel|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Gases|Efuel|Other Taxes

    # Price|Final Energy|Residential and Commercial|Hydrogen x
    # Price|Final Energy|Residential and Commercial|Hydrogen|Sales Margin
    # Price|Final Energy|Residential and Commercial|Hydrogen|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Hydrogen|Other Taxes

    # Price|Final Energy|Residential and Commercial|Electricity
    nodal_flows_lv = get_nodal_flows(
        n,
        "low voltage",
        region,
        query="not carrier.str.contains('agriculture')"
        "& not carrier.str.contains('industry')"
        "& not carrier.str.contains('urban central')",
    )
    nodal_prices_lv = n.buses_t.marginal_price[nodal_flows_lv.columns]
    var["Price|Final Energy|Residential and Commercial|Electricity"] = (
        nodal_flows_lv.mul(nodal_prices_lv).values.sum() / nodal_flows_lv.values.sum()
    )
    # Price|Final Energy|Residential and Commercial|Electricity|Sales Margin x
    # Price|Final Energy|Residential and Commercial|Electricity|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Electricity|Other Taxes

    # Price|Final Energy|Residential and Commercial|Electricity Heating
    # Price|Final Energy|Residential and Commercial|Electricity Heating|Sales Margin
    # Price|Final Energy|Residential and Commercial|Electricity Heating|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Electricity Heating|Other Taxes

    # Price|Final Energy|Residential and Commercial|Gases|SNG
    # Price|Final Energy|Residential and Commercial|Gases|SNG|Sales Margin
    # Price|Final Energy|Residential and Commercial|Gases|SNG|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Gases|SNG|Other Taxes
    # Price|Final Energy|Residential and Commercial|Gases|SNG|Carbon Price Component

    ### Price|Final Energy|Industry|

    # Price|Final Energy|Industry|Liquids|Oil
    var["Price|Final Energy|Industry|Liquids|Oil"] = get_weighted_costs_links(
        ["naphtha for industry"], n, region
    )
    # Price|Final Energy|Industry|Liquids|Oil|Sales Margin
    # Price|Final Energy|Industry|Liquids|Oil|Transport and Distribution
    # Price|Final Energy|Industry|Liquids|Oil|Carbon Price Component
    # Price|Final Energy|Industry|Liquids|Oil|Other Taxes

    var["Price|Final Energy|Industry|Solids|Coal"] = price_load(
        n, "coal for industry", region
    )[0]
    # Price|Final Energy|Industry|Solids|Coal|Sales Margin x
    # Price|Final Energy|Industry|Solids|Coal|Transport and Distribution
    # Price|Final Energy|Industry|Solids|Coal|Carbon Price Component
    # Price|Final Energy|Industry|Solids|Coal|Other Taxes

    # var["Price|Final Energy|Industry|Gases|Natural Gas"] ?
    # Price|Final Energy|Industry|Gases|Natural Gas|Sales Margin x
    # Price|Final Energy|Industry|Gases|Natural Gas|Transport and Distribution
    # Price|Final Energy|Industry|Gases|Natural Gas|Carbon Price Component
    # Price|Final Energy|Industry|Gases|Natural Gas|Other Taxes

    # Price|Final Energy|Industry|Gases|SNG
    # Price|Final Energy|Industry|Gases|SNG|Sales Margin
    # Price|Final Energy|Industry|Gases|SNG|Transport and Distribution
    # Price|Final Energy|Industry|Gases|SNG|Carbon Price Component
    # Price|Final Energy|Industry|Gases|SNG|Other Taxes

    # Price|Final Energy|Industry|Heat|Sales Margin x
    # Price|Final Energy|Industry|Heat|Transport and Distribution
    # Price|Final Energy|Industry|Heat|Other Taxes

    # Price|Final Energy|Industry|Liquids|Biomass x
    # Price|Final Energy|Industry|Liquids|Biomass|Sales Margin
    # Price|Final Energy|Industry|Liquids|Biomass|Transport and Distribution
    # Price|Final Energy|Industry|Liquids|Biomass|Other Taxes

    var["Price|Final Energy|Industry|Solids|Biomass"] = price_load(
        n, "solid biomass for industry", region
    )[0]
    # Price|Final Energy|Industry|Solids|Biomass|Sales Margin x
    # Price|Final Energy|Industry|Solids|Biomass|Transport and Distribution
    # Price|Final Energy|Industry|Solids|Biomass|Other Taxes

    # Price|Final Energy|Industry|Gases|Biomass ?
    # Price|Final Energy|Industry|Gases|Biomass|Sales Margin
    # Price|Final Energy|Industry|Gases|Biomass|Transport and Distribution
    # Price|Final Energy|Industry|Gases|Biomass|Other Taxes

    var["Price|Final Energy|Industry|Liquids|Efuel"] = price_load(
        n, "naphtha for industry", region
    )[0]

    # Price|Final Energy|Industry|Liquids|Efuel|Sales Margin x
    # Price|Final Energy|Industry|Liquids|Efuel|Transport and Distribution
    # Price|Final Energy|Industry|Liquids|Efuel|Other Taxes

    # Price|Final Energy|Industry|Gases|Efuel x
    # Price|Final Energy|Industry|Gases|Efuel|Sales Margin
    # Price|Final Energy|Industry|Gases|Efuel|Transport and Distribution
    # Price|Final Energy|Industry|Gases|Efuel|Other Taxes

    # Price|Final Energy|Industry|Hydrogen|Sales Margin x
    # Price|Final Energy|Industry|Hydrogen|Transport and Distribution
    # Price|Final Energy|Industry|Hydrogen|Other Taxes

    # Price|Final Energy|Industry|Electricity|Sales Margin x
    # Price|Final Energy|Industry|Electricity|Transport and Distribution
    # Price|Final Energy|Industry|Electricity|Other Taxes

    # Price|Final Energy|Industry|Electricity Heating
    # Price|Final Energy|Industry|Electricity Heating|Sales Margin
    # Price|Final Energy|Industry|Electricity Heating|Transport and Distribution
    # Price|Final Energy|Industry|Electricity Heating|Other Taxes

    # Price|Final Energy|Industry|Electricity
    var["Price|Final Energy|Industry|Electricity"] = price_load(
        n, "industry electricity", region
    )[0]
    # Price|Final Energy|Industry|Electricity|Sales Margin
    # Price|Final Energy|Industry|Electricity|Transport and Distribution
    # Price|Final Energy|Industry|Electricity|Other Taxes

    # Price|Final Energy|Industry|Gases
    var["Price|Final Energy|Industry|Gases"] = get_weighted_costs_links(
        ["gas for industry", "gas for industry CC"], n, region
    )
    # Price|Final Energy|Industry|Heat
    var["Price|Final Energy|Industry|Heat"] = price_load(
        n, "low-temperature heat for industry", region
    )[0]
    # Price|Final Energy|Industry|Liquids
    carriers = ["naphtha for industry", "industry methanol"]
    df = pd.DataFrame({c: price_load(n, c, region) for c in carriers})
    var["Price|Final Energy|Industry|Liquids"] = (
        df.iloc[0] * df.iloc[1]
    ).sum() / df.iloc[1].sum()
    # Price|Final Energy|Industry|Hydrogen
    var["Price|Final Energy|Industry|Hydrogen"] = price_load(
        n, "H2 for industry", region
    )[0]
    # Price|Final Energy|Industry|Solids
    carriers = [
        "solid biomass for industry",
        "solid biomass for industry CC",
        "coal for industry",
    ]
    var["Price|Final Energy|Industry|Solids"] = get_weighted_costs_links(
        carriers, n, region
    )

    return var


def get_discretized_value(value, disc_int, build_threshold=0.3):
    if value == 0.0:
        return value

    remainder = value % disc_int
    base = value - remainder
    discrete = disc_int if remainder > build_threshold * disc_int else 0.0

    return base + discrete


def get_grid_investments(
    n,
    region,
    scope="all",  # all, baseyear, expanded
    var_name="Investment|Energy Supply|Electricity|Transmission|",
):
    assert scope in ["all", "baseyear", "expanded"]
    # TODO gap between years should be read from config
    var = pd.Series()

    offwind = n.generators.filter(like="offwind", axis=0).filter(like="DE", axis=0)

    offwind_capacity = offwind.p_nom_opt
    if scope == "expanded":
        offwind_capacity -= offwind.p_nom
    elif scope == "baseyear":
        # WARNING USING GLOBAL VARIABLE `networks[0]`!!!
        # Subtracting 2020 capacity
        offwind2020 = (
            networks[0]
            .generators.filter(like="offwind", axis=0)
            .filter(like="DE", axis=0)
            .p_nom
        )
        common_index = offwind.index.intersection(offwind2020.index)
        offwind_capacity.loc[common_index] -= offwind2020.loc[common_index]
    offwind_connection_overnight_cost = (
        offwind_capacity * offwind.connection_overnight_cost
    ) * 1e-9
    offwind_connection_ac = offwind_connection_overnight_cost.filter(like="ac")
    offwind_connection_dc = offwind_connection_overnight_cost.filter(regex="dc|float")

    var[var_name + "AC|Offshore"] = offwind_connection_ac.sum() / 5
    var[var_name + "AC|NEP|Offshore"] = (
        offwind_connection_ac.filter(regex="25|30").sum() / 5
    )

    var[var_name + "DC|Offshore"] = offwind_connection_dc.sum() / 5
    var[var_name + "DC|NEP|Offshore"] = (
        offwind_connection_dc.filter(regex="25|30").sum() / 5
    )

    dc_links = n.links[
        (n.links.carrier == "DC")
        & (n.links.bus0 + n.links.bus1).str.contains(region)
        & ~n.links.reversed
    ]
    current_year = n.generators.build_year.max()  # noqa
    nep_dc = dc_links.query(
        "(index.str.startswith('DC') or index.str.startswith('TYNDP')) and build_year > 2025 and (@current_year - 5 < build_year <= @current_year)"
    ).index

    dc_capacity = dc_links.p_nom_opt
    if scope == "expanded":
        dc_capacity -= dc_links.p_nom_min
    elif scope == "baseyear":
        # WARNING USING GLOBAL VARIABLE!!!
        # Subtracting 2020 capacity
        dc_capacity -= networks[0].links.loc[dc_links.index].p_nom_min

    dc_investments = dc_capacity * dc_links.overnight_cost * 1e-9
    # International dc_projects are only accounted with half the costs
    dc_investments[
        ~(dc_links.bus0.str.contains(region) & dc_links.bus1.str.contains(region))
    ] *= 0.5

    ac_lines = n.lines[(n.lines.bus0 + n.lines.bus1).str.contains(region)]
    nep_ac = ac_lines.query(
        "(build_year > 2025) and (@current_year - 5 < build_year <= @current_year)"
    ).index
    # Assuming the lines are already post-discretized
    ac_capacity = ac_lines.s_nom_opt
    if scope == "expanded":
        ac_capacity -= ac_lines.s_nom_min
    elif scope == "baseyear":
        # WARNING USING GLOBAL VARIABLE!!!
        # Subtracting 2020 capacity
        ac_capacity -= networks[0].lines.loc[ac_lines.index].s_nom_min

    ac_investments = ac_capacity * ac_lines.overnight_cost * 1e-9
    # International ac_projects are only accounted with half the costs
    ac_investments[
        ~(ac_lines.bus0.str.contains(region) & ac_lines.bus1.str.contains(region))
    ] *= 0.5

    # Blindleistungskompensation
    # https://www.netzentwicklungsplan.de/sites/default/files/2023-07/NEP_2037_2045_V2023_2_Entwurf_Teil1_1.pdf
    # Tabelle 30, Abbildung 70, Kostenannahmen NEP + eigene Berechnungen, gerundet
    year = n.generators.build_year.max()
    reactive_power_compensation = (
        pd.Series(
            {
                2020: 0,
                2025: 4.4,
                2030: 8,
                2035: 15,
                2040: 10,
                2045: 1.5,
            }
        )
        / EUR20TOEUR23
    )
    var[var_name + "AC|Übernahme|Reactive Power Compensation"] = (
        reactive_power_compensation.get(year, 0) / 5
    )
    var[var_name + "AC|Übernahme|Startnetz Delta"] = 0

    var[var_name + "AC|Onshore"] = ac_investments.sum() / 5
    var[var_name + "AC|NEP|Onshore"] = ac_investments[nep_ac].sum() / 5
    var[var_name + "DC|Onshore"] = dc_investments.sum() / 5
    var[var_name + "DC|NEP|Onshore"] = dc_investments[nep_dc].sum() / 5

    for key in ["Onshore", "NEP|Onshore", "Offshore", "NEP|Offshore"]:
        var[var_name + f"{key}"] = (
            var[var_name + f"AC|{key}"] + var[var_name + f"DC|{key}"]
        )

    var[var_name + "AC"] = (
        var[var_name + "AC|Onshore"]
        + var[var_name + "AC|Offshore"]
        + var[var_name + "AC|Übernahme|Reactive Power Compensation"]
    )
    var[var_name + "AC|NEP"] = (
        var[var_name + "AC|NEP|Onshore"]
        + var[var_name + "AC|NEP|Offshore"]
        + var[var_name + "AC|Übernahme|Reactive Power Compensation"]
    )
    var[var_name + "DC"] = var[var_name + "DC|Onshore"] + var[var_name + "DC|Offshore"]
    var[var_name + "DC|NEP"] = (
        var[var_name + "DC|NEP|Onshore"] + var[var_name + "DC|NEP|Offshore"]
    )
    var["Investment|Energy Supply|Electricity|Transmission"] = (
        var[var_name + "AC"] + var[var_name + "DC"]
    )
    var[var_name + "NEP"] = var[var_name + "AC|NEP"] + var[var_name + "DC|NEP"]

    distribution_grid = n.links[
        (n.links.carrier == "electricity distribution grid")
        & n.links.bus0.str.contains(region)
        & ~n.links.reversed
    ]

    year = distribution_grid.build_year.max()
    year_pre = (year - 5) if year > 2020 else 2020

    dg_capacity = distribution_grid.p_nom_opt.sum()
    if scope == "expanded":
        dg_capacity -= distribution_grid[
            distribution_grid.build_year <= year_pre
        ].p_nom_opt.sum()
    elif scope == "baseyear":
        dg_capacity -= distribution_grid[
            distribution_grid.build_year <= 2020
        ].p_nom_opt.sum()

    dg_investment = (
        dg_capacity * distribution_grid.overnight_cost.unique().item() * 1e-9
    )
    var["Investment|Energy Supply|Electricity|Distribution"] = dg_investment / 5

    var["Investment|Energy Supply|Electricity|Transmission and Distribution"] = (
        var["Investment|Energy Supply|Electricity|Distribution"]
        + var["Investment|Energy Supply|Electricity|Transmission"]
    )

    h2_links = n.links[
        n.links.carrier.str.contains("H2 pipeline")
        & ~n.links.reversed
        & (n.links.bus0 + n.links.bus1).str.contains(region)
    ]
    year = h2_links.build_year.max()

    if scope == "expanded":
        new_h2_links = h2_links[
            ((year - 5) < h2_links.build_year) & (h2_links.build_year <= year)
        ]
    else:  # scope is all or baseyear
        new_h2_links = h2_links.copy()

    h2_expansion = new_h2_links.p_nom_opt
    h2_investments = h2_expansion * new_h2_links.overnight_cost * 1e-9
    # International h2_projects are only accounted with domestic_length_factor * costs
    if len(h2_links.carrier.unique()) == 1:
        dlf = domestic_length_factor(n, h2_links.carrier.unique().tolist(), region)
    else:
        dlf = np.array(
            list(
                domestic_length_factor(
                    n, h2_links.carrier.unique().tolist(), region
                ).values()
            )
        ).mean()

    h2_investments[
        ~(
            new_h2_links.bus0.str.contains(region)
            & new_h2_links.bus1.str.contains(region)
        )
    ] *= dlf

    var["Investment|Energy Supply|Hydrogen|Transmission and Distribution"] = var[
        "Investment|Energy Supply|Hydrogen|Transmission"
    ] = h2_investments.sum() / 5

    new_h2_links_kernnetz_i = new_h2_links[
        (new_h2_links.index.str.contains("kernnetz"))
    ].index

    new_h2_links_endogen_i = new_h2_links[
        ~(new_h2_links.index.str.contains("kernnetz"))
    ].index

    var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen"] = (
        h2_investments[new_h2_links_endogen_i].sum() / 5
    )
    var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz"] = (
        h2_investments[new_h2_links_kernnetz_i].sum() / 5
    )

    assert isclose(
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution"],
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen"]
        + var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz"
        ],
    )

    if "retrofitted" in new_h2_links.columns:
        new_h2_links_retrofitted_i = new_h2_links[
            (new_h2_links.retrofitted == 1.0)
            | (new_h2_links.index.str.contains("retrofitted"))
        ].index
    else:
        new_h2_links_retrofitted_i = new_h2_links[
            (new_h2_links.index.str.contains("retrofitted"))
        ].index

    new_h2_links_newbuild_i = new_h2_links.index.difference(new_h2_links_retrofitted_i)

    var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|New-build"] = (
        h2_investments[new_h2_links_newbuild_i].sum() / 5
    )
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Retrofitted"
    ] = h2_investments[new_h2_links_retrofitted_i].sum() / 5

    assert isclose(
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution"],
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|New-build"]
        + var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Retrofitted"
        ],
    )

    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen|New-build"
    ] = (
        h2_investments[
            new_h2_links_newbuild_i.intersection(new_h2_links_endogen_i)
        ].sum()
        / 5
    )
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen|Retrofitted"
    ] = (
        h2_investments[
            new_h2_links_retrofitted_i.intersection(new_h2_links_endogen_i)
        ].sum()
        / 5
    )
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|New-build"
    ] = (
        h2_investments[
            new_h2_links_newbuild_i.intersection(new_h2_links_kernnetz_i)
        ].sum()
        / 5
    )
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|Retrofitted"
    ] = (
        h2_investments[
            new_h2_links_retrofitted_i.intersection(new_h2_links_kernnetz_i)
        ].sum()
        / 5
    )

    assert isclose(
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen"],
        var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen|New-build"
        ]
        + var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Endogen|Retrofitted"
        ],
    )

    assert isclose(
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz"],
        var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|New-build"
        ]
        + var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|Retrofitted"
        ],
    )

    if "tags" in new_h2_links.columns:
        # extract infos from tags
        tags = new_h2_links.loc[new_h2_links_kernnetz_i].tags.values
        df = pd.DataFrame(tags, columns=["info"])
        df["info"] = df["info"].apply(ast.literal_eval)
        df["pci"] = df["info"].apply(lambda x: x["pci"])
        df["ipcei"] = df["info"].apply(lambda x: x["ipcei"])
        df["investment_costs (Mio. Euro)"] = df["info"].apply(
            lambda x: x["investment_costs (Mio. Euro)"]
        )
        df.index = new_h2_links_kernnetz_i

        pci_i = df[df["pci"] != "no"].index
        ipcei_i = df[df["ipcei"] != "no"].index
    else:
        pci_i = []
        ipcei_i = []

    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|PCI"
    ] = h2_investments[pci_i].sum() / 5
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|IPCEI"
    ] = h2_investments[ipcei_i].sum() / 5
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|PCI+IPCEI"
    ] = h2_investments[pci_i.union(ipcei_i)].sum() / 5
    var[
        "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|NOT-PCI+IPCEI"
    ] = (
        h2_investments[new_h2_links_kernnetz_i.difference(pci_i.union(ipcei_i))].sum()
        / 5
    )

    assert isclose(
        var["Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz"],
        var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|PCI+IPCEI"
        ]
        + var[
            "Investment|Energy Supply|Hydrogen|Transmission and Distribution|Kernnetz|NOT-PCI+IPCEI"
        ],
    )

    # var["Investment|Energy Supply|Electricity|Electricity Storage"] = \
    # var["Investment|Energy Supply|CO2 Transport and Storage"] =
    # var["Investment|Energy Supply|Other"] =
    # var["Investment|Energy Efficiency"] = \
    # var["Investment|Energy Supply|Heat"] = \
    # var["Investment|Energy Supply|Hydrogen"] = \
    # var["Investment|RnD|Energy Supply"] = \
    # var["Investment|Energy Demand|Transportation"] = \
    # var["Investment|Energy Demand|Transportation|LDV"] = \
    # var["Investment|Energy Demand|Transportation|Bus"] = \
    # var["Investment|Energy Demand|Transportation|Rail"] = \
    # var["Investment|Energy Demand|Transportation|Truck"] = \
    # var["Investment|Infrastructure|Transport"] = \
    # var["Investment|Energy Demand|Residential and Commercial"] = \
    # var["Investment|Energy Demand|Residential and Commercial|Low-Efficiency Buildings"] = \
    # var["Investment|Energy Demand|Residential and Commercial|Medium-Efficiency Buildings"] = \
    # var["Investment|Energy Demand|Residential and Commercial|High-Efficiency Buildings"] = \
    # var["Investment|Energy Demand|Residential and Commercial|Building Retrofits"] =
    # var["Investment|Energy Demand|Residential and Commercial|Space Heating"] = \
    # var["Investment|Infrastructure|Industry|Green"] = \
    # var["Investment|Infrastructure|Industry|Non-Green"] = \
    # var["Investment|Industry"] = \
    return var


def get_policy(n, investment_year):
    var = pd.Series()
    try:
        co2_limit_de = n.global_constraints.loc["co2_limit-DE", "mu"]
    except KeyError:
        co2_limit_de = 0
    var["Price|Carbon"] = -n.global_constraints.loc["CO2Limit", "mu"] - co2_limit_de

    var["Price|Carbon|EU-wide Regulation All Sectors"] = -n.global_constraints.loc[
        "CO2Limit", "mu"
    ]

    var["Price|Carbon|National Climate Target"] = -co2_limit_de

    try:
        var["Price|Policy|Electricity Import Independence"] = n.global_constraints.loc[
            "Electricity_import_limit-DE", "mu"
        ]
    except KeyError:
        var["Price|Policy|Electricity Import Independence"] = 0

    try:
        var["Price|Policy|Hydrogen Import Independence"] = n.global_constraints.loc[
            "H2_import_limit-DE", "mu"
        ]
    except KeyError:
        var["Price|Policy|Hydrogen Import Independence"] = 0

    try:
        var["Price|Policy|Methanol Import Independence"] = n.global_constraints.loc[
            "methanol_import_limit-DE", "mu"
        ]
    except KeyError:
        var["Price|Policy|Methanol Import Independence"] = 0

    try:
        var["Price|Policy|Hydrogen Export Ban"] = n.global_constraints.loc[
            "H2_export_ban-DE", "mu"
        ]
    except KeyError:
        var["Price|Policy|Hydrogen Export Ban"] = 0

    try:
        var["Price|Policy|Renewable Oil Import Independence"] = (
            n.global_constraints.loc["renewable_oil_import_limit-DE", "mu"]
        )
    except KeyError:
        var["Price|Policy|Renewable Oil Import Independence"] = 0

    try:
        var["Price|Policy|Renewable Gas Import Independence"] = (
            n.global_constraints.loc["renewable_gas_import_limit-DE", "mu"]
        )
    except KeyError:
        var["Price|Policy|Renewable Gas Import Independence"] = 0

    try:
        var["Price|Policy|Hydrogen Derivate Import Independence"] = (
            n.global_constraints.loc["H2_derivate_import_limit-DE", "mu"]
        )
    except KeyError:
        var["Price|Policy|Hydrogen Derivate Import Independence"] = 0

    return var


def get_economy(n, region):
    var = pd.Series()

    def get_tsc(n, country):
        pypsa.options.params.statistics.drop_zero = False
        capex = n.statistics.capex(
            groupby=pypsa.statistics.groupers["name", "carrier"], nice_names=False
        )

        opex = n.statistics.opex(
            groupby=pypsa.statistics.groupers["name", "carrier"], nice_names=False
        )

        # filter inter country transmission lines and links
        inter_country_lines = n.lines.bus0.map(n.buses.country) != n.lines.bus1.map(
            n.buses.country
        )
        inter_country_links = n.links.bus0.map(n.buses.country) != n.links.bus1.map(
            n.buses.country
        )
        #
        transmission_carriers = get_transmission_carriers(n).get_level_values("carrier")
        transmission_lines = (
            n.lines.carrier.isin(transmission_carriers) & n.lines.active
        )
        transmission_links = (
            n.links.carrier.isin(transmission_carriers) & n.links.active
        )
        #
        country_transmission_lines = (
            (n.lines.bus0.str.contains(country)) & ~(n.lines.bus1.str.contains(country))
        ) | (
            ~(n.lines.bus0.str.contains(country)) & (n.lines.bus1.str.contains(country))
        )
        country_tranmission_links = (
            (n.links.bus0.str.contains(country)) & ~(n.links.bus1.str.contains(country))
        ) | (
            ~(n.links.bus0.str.contains(country)) & (n.links.bus1.str.contains(country))
        )
        #
        inter_country_transmission_lines = (
            inter_country_lines & transmission_lines & country_transmission_lines
        )
        inter_country_transmission_links = (
            inter_country_links & transmission_links & country_tranmission_links
        )
        inter_country_transmission_lines_i = inter_country_transmission_lines[
            inter_country_transmission_lines
        ].index
        inter_country_transmission_links_i = inter_country_transmission_links[
            inter_country_transmission_links
        ].index
        inter_country_transmission_i = inter_country_transmission_lines_i.union(
            inter_country_transmission_links_i
        )

        #
        tsc = pd.concat([capex, opex], axis=1, keys=["capex", "opex"])
        tsc = tsc.reset_index().set_index("name")
        tsc.loc[inter_country_transmission_i, ["capex", "opex"]] = (
            tsc.loc[inter_country_transmission_i, ["capex", "opex"]] / 2
        )
        tsc.rename(
            index={
                index: index + " " + country for index in inter_country_transmission_i
            },
            inplace=True,
        )
        # rename inter region links and lines
        to_rename_links = n.links[
            (n.links.bus0.str.contains(region))
            & (n.links.bus1.str.contains(region))
            & ~(n.links.index.str.contains(region))
        ].index
        to_rename_lines = n.lines[
            (n.lines.bus0.str.contains(region))
            & (n.lines.bus1.str.contains(region))
            & ~(n.lines.index.str.contains(region))
        ].index
        tsc.rename(
            index={index: index + " " + region for index in to_rename_links},
            inplace=True,
        )
        tsc.rename(
            index={index: index + " " + region for index in to_rename_lines},
            inplace=True,
        )

        tsc = (
            tsc.filter(like=country, axis=0)
            .drop("component", axis=1)
            .groupby("carrier")
            .sum()
        )

        return tsc

    def get_trade_cost(n, region, carriers):
        """
        Positive values mean cost for the domestic energy system (imports > exports)
        Negative values mean revenue for the domestic energy system (exports > imports)
        """
        export_revenue, import_cost = get_export_import(n, region, carriers, unit="€")
        return import_cost - export_revenue

    var["Total Energy System Cost|Trade|Electricity"] = (
        get_trade_cost(n, region, ["AC"]) + get_trade_cost(n, region, ["DC"])
    ) / 1e9
    var["Total Energy System Cost|Trade|Efuels"] = (
        get_trade_cost(n, region, ["renewable oil", "renewable gas", "methanol"]) / 1e9
    )
    var["Total Energy System Cost|Trade|Hydrogen"] = (
        get_trade_cost(
            n,
            region,
            ["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"],
        )
        / 1e9
    )
    var["Total Energy System Cost|Trade"] = (
        var["Total Energy System Cost|Trade|Electricity"]
        + var["Total Energy System Cost|Trade|Efuels"]
        + var["Total Energy System Cost|Trade|Hydrogen"]
    )
    # Total Energy System Cost in billion EUR2020/yr
    var["Total Energy System Cost|Non Trade"] = get_tsc(n, region).sum().sum() / 1e9

    var["Total Energy System Cost"] = (
        var["Total Energy System Cost|Non Trade"]
        + var["Total Energy System Cost|Trade"]
    )

    var["Total Energy System Cost|EU"] = (
        n.statistics.capex().sum() + n.statistics.opex().sum()
    ) / 1e9

    rev = n.statistics.revenue(groupby=["bus"]).loc["Load"]
    var["Total Energy System Cost|Load Revenue"] = (
        rev[rev.index.str.startswith(region)].sum() / 1e9
    )

    return var


def get_export_import(n, region, carriers, aggregate=True, unit="MWh"):
    # note: links can also used bidirectional if efficiency=1 (e.g. "H2 pipeline retrofitted")
    if "AC" not in carriers:
        df = n.links
        df_t = n.links_t
    if "AC" in carriers:
        if len(carriers) > 1:
            raise NotImplementedError(
                "AC lines cannot be combined with other carriers. Use carrier=['AC'] to get export_import for n.lines."
            )
        df = n.lines
        df_t = n.lines_t

    outgoing = df[
        (df.carrier.isin(carriers))
        & (df.bus0.str[:2] == region)
        & (df.bus1.str[:2] != region)
    ]

    incoming = df[
        (df.carrier.isin(carriers))
        & (df.bus0.str[:2] != region)
        & (df.bus1.str[:2] == region)
    ]
    if outgoing.empty and incoming.empty:
        logger.warning(
            f"When computing imports and export no {carriers} links found between {region} and non-{region} buses. Returning 0 or empty series."
        )
        if aggregate:
            return 0, 0
        else:
            return pd.Series(), pd.Series()
    # if p0 > 0 (=clip(lower=0)) system is withdrawing from bus0 (DE) and feeding into bus1 (non-DE) -> export
    export_outgoing = (
        df_t.p0.loc[:, outgoing.index]
        .clip(lower=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
    )
    if unit == "€":
        # bus0 is DE for outgoing links
        domestic_prices = pd.concat(
            [
                n.buses_t.marginal_price[bus].rename(link)
                for link, bus in outgoing.bus0.items()
            ],
            axis=1,
        )
        export_outgoing *= domestic_prices

    # if p1 > 0 system is withdrawing from bus1 (DE) and feeding into bus0 (non-DE) -> export
    export_incoming = (
        df_t.p1.loc[:, incoming.index]
        .clip(lower=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
    )
    if unit == "€":
        # bus1 is DE for incoming links
        domestic_prices = pd.concat(
            [
                n.buses_t.marginal_price[bus].rename(link)
                for link, bus in incoming.bus1.items()
            ],
            axis=1,
        )
        export_incoming *= domestic_prices

    exporting_p = pd.concat([export_outgoing, export_incoming], axis=1)
    if aggregate:
        exporting_p = exporting_p.values.sum()

    # if p1 < 0 (=clip(upper=0)) system is feeding into bus1 (DE) and withdrawing from bus0 (non-DE) -> import (with negative sign here)
    import_incoming = (
        df_t.p1.loc[:, incoming.index]
        .clip(upper=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
        * -1
    )
    if unit == "€":
        # bus1 is DE for incoming links
        domestic_prices = pd.concat(
            [
                n.buses_t.marginal_price[bus].rename(link)
                for link, bus in incoming.bus1.items()
            ],
            axis=1,
        )
        import_incoming *= domestic_prices

    # if p0 < 0 (=clip(upper=0)) system is feeding into bus0 (DE) and withdrawing from bus1 (non-DE) -> import (with negative sign here)
    import_outgoing = (
        df_t.p0.loc[:, outgoing.index]
        .clip(upper=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
        * -1
    )
    if unit == "€":
        # bus0 is DE for outgoing links
        domestic_prices = pd.concat(
            [
                n.buses_t.marginal_price[bus].rename(link)
                for link, bus in outgoing.bus0.items()
            ],
            axis=1,
        )
        import_outgoing *= domestic_prices

    importing_p = pd.concat([import_outgoing, import_incoming], axis=1)
    if aggregate:
        importing_p = importing_p.values.sum()

    return exporting_p, importing_p


def get_trade(n, region):
    var = pd.Series()

    # Trade|Secondary Energy|Electricity|Volume
    exports_ac, imports_ac = get_export_import(n, region, ["AC"])

    exports_dc, imports_dc = get_export_import(n, region, ["DC"])

    var["Trade|Secondary Energy|Electricity|Volume"] = (exports_ac - imports_ac) + (
        exports_dc - imports_dc
    )
    var["Trade|Secondary Energy|Electricity|Gross Import|Volume"] = (
        imports_ac + imports_dc
    )
    # var["Trade|Secondary Energy|Electricity|Volume|Exports"] = \
    #     (exporting_p_ac + exports_dc)

    # Trade|Secondary Energy|Hydrogen|Volume
    h2_carriers = ["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"]
    exports_h2, imports_h2 = get_export_import(n, region, h2_carriers)
    var["Trade|Secondary Energy|Hydrogen|Volume"] = exports_h2 - imports_h2
    var["Trade|Secondary Energy|Hydrogen|Gross Import|Volume"] = imports_h2
    # var["Trade|Secondary Energy|Hydrogen|Volume|Exports"] = \
    #     exports_h2

    # Trade|Secondary Energy|Liquids|Hydrogen|Volume

    renewable_oil_supply = (
        n.statistics.supply(bus_carrier="renewable oil", **kwargs)
        .groupby(["bus", "carrier"])
        .sum()
    )

    assert region == "DE"  # only DE is implemented at the moment

    DE_renewable_oil = renewable_oil_supply.get("DE renewable oil", pd.Series(0))
    EU_renewable_oil = renewable_oil_supply.get("EU renewable oil", pd.Series(0))

    if DE_renewable_oil.sum() == 0:
        DE_bio_fraction = 0
    else:
        DE_bio_fraction = (
            DE_renewable_oil.filter(like="bio").sum() / DE_renewable_oil.sum()
        )

    if EU_renewable_oil.sum() == 0:
        EU_bio_fraction = 0
    else:
        EU_bio_fraction = (
            EU_renewable_oil.filter(like="bio").sum() / EU_renewable_oil.sum()
        )

    exports_oil_renew, imports_oil_renew = get_export_import(
        n, region, ["renewable oil"]
    )

    var["Trade|Secondary Energy|Liquids|Biomass|Volume"] = (
        exports_oil_renew * DE_bio_fraction - imports_oil_renew * EU_bio_fraction
    )

    var["Trade|Secondary Energy|Liquids|Biomass|Gross Import|Volume"] = (
        imports_oil_renew * EU_bio_fraction
    )

    exports_meoh, imports_meoh = get_export_import(n, region, ["methanol"])

    var["Trade|Secondary Energy|Liquids|Hydrogen|Volume"] = (
        exports_oil_renew * (1 - DE_bio_fraction)
        - imports_oil_renew * (1 - EU_bio_fraction)
        + exports_meoh
        - imports_meoh
    )

    var["Trade|Secondary Energy|Liquids|Hydrogen|Gross Import|Volume"] = (
        imports_oil_renew * (1 - EU_bio_fraction) + imports_meoh
    )

    var["Trade|Secondary Energy|Methanol|Volume"] = exports_meoh - imports_meoh

    var["Trade|Secondary Energy|Methanol|Gross Import|Volume"] = imports_meoh

    # Trade|Secondary Energy|Gases|Hydrogen|Volume

    renewable_gas_supply = (
        n.statistics.supply(bus_carrier="renewable gas", **kwargs)
        .groupby(["bus", "carrier"])
        .sum()
    )

    DE_renewable_gas = renewable_gas_supply.get("DE renewable gas", pd.Series(0))
    EU_renewable_gas = renewable_gas_supply.get("EU renewable gas", pd.Series(0))

    if DE_renewable_gas.sum() == 0:
        DE_bio_fraction = 0
    else:
        DE_bio_fraction = (
            DE_renewable_gas.filter(like="bio").sum() / DE_renewable_gas.sum()
        )

    if EU_renewable_gas.sum() == 0:
        EU_bio_fraction = 0
    else:
        EU_bio_fraction = (
            EU_renewable_gas.filter(like="bio").sum() / EU_renewable_gas.sum()
        )

    assert region == "DE"  # only DE is implemented at the moment

    exports_gas_renew, imports_gas_renew = get_export_import(
        n, region, ["renewable gas"]
    )
    var["Trade|Secondary Energy|Gases|Hydrogen|Volume"] = exports_gas_renew * (
        1 - DE_bio_fraction
    ) - imports_gas_renew * (1 - EU_bio_fraction)
    var["Trade|Secondary Energy|Gases|Hydrogen|Gross Import|Volume"] = (
        imports_gas_renew * (1 - EU_bio_fraction)
    )

    var["Trade|Secondary Energy|Gases|Biomass|Volume"] = (
        exports_gas_renew * DE_bio_fraction - imports_gas_renew * EU_bio_fraction
    )
    var["Trade|Secondary Energy|Gases|Biomass|Gross Import|Volume"] = (
        imports_gas_renew * EU_bio_fraction
    )

    # Trade|Primary Energy|Coal|Volume
    # Trade|Primary Energy|Gas|Volume

    gas_fractions = _get_fuel_fractions(n, region, "gas")

    if "gas pipeline" in n.links.carrier.unique():
        exports_gas, imports_gas = get_export_import(
            n, region, ["gas pipeline", "gas pipeline new"]
        )
        var["Trade|Primary Energy|Gas|Volume"] = (
            exports_gas - imports_gas
        ) * gas_fractions["Natural Gas"]
        var["Trade|Primary Energy|Gas|Volume|Imports"] = (
            imports_gas * gas_fractions["Natural Gas"]
        )
        var["Trade|Primary Energy|Gas|Volume|Exports"] = (
            exports_gas * gas_fractions["Natural Gas"]
        )

    # Trade|Primary Energy|Oil|Volume

    # Biomass Trade

    biomass_primary_gens = n.generators.query(
        f"index.str.startswith('{region}') and index.str.endswith('solid biomass')"
    )  # use endswith to avoid the biomass transport generators
    biomass_transport_gens = n.generators.query(
        f"index.str.startswith('{region}') and index.str.endswith('solid biomass transported')"
    )

    local_biomass_potential = biomass_primary_gens.e_sum_max.sum()

    local_biomass_usage = (
        (
            n.generators_t.p[
                biomass_primary_gens.index.union(biomass_transport_gens.index)
            ]
        )
        .sum(axis=1)
        .multiply(n.snapshot_weightings.generators)
        .sum()
    )

    biomass_imports = local_biomass_usage - local_biomass_potential

    var["Trade|Primary Energy|Biomass|Net Imports"] = biomass_imports

    logger.info(
        f"""Share of imported biomass: {
            round(biomass_imports / local_biomass_usage, 3)
        }"""
    )

    return var * MWh2TWh


def get_production(region, year):
    var = pd.Series()
    # read in the industrial production data
    years = [
        int(re.search(r"(\d{4})-modified\.csv", filename).group(1))
        for filename in snakemake.input.industrial_production_per_country_tomorrow
    ]
    index = next((idx for idx, y in enumerate(years) if y == year), None)
    production = pd.read_csv(
        snakemake.input.industrial_production_per_country_tomorrow[index], index_col=0
    ).div(1e3)  # kton/a -> Mton/a

    var["Production|Non-Metallic Minerals|Cement"] = production.loc[region, "Cement"]
    var["Production|Steel"] = production.loc[
        "DE", ["Electric arc", "Integrated steelworks", "DRI + Electric arc"]
    ].sum()
    var["Production|Steel|Primary"] = (
        var["Production|Steel"] * config_industry["St_primary_fraction"][year]
    )
    var["Production|Steel|Secondary"] = var["Production|Steel"] * (
        1 - config_industry["St_primary_fraction"][year]
    )

    # optional:
    # var[""Production|Pulp and Paper"]
    # var["Production|Chemicals|Ammonia"]
    # var["Production|Chemicals|Methanol"]
    # var["Production|Pulp and Paper"]
    # var["Production|Non-Ferrous Metals"]

    return var


def get_grid_capacity(n, region, year):
    var = pd.Series()
    ### Total Capacity
    ## Transmission Grid
    # TODO add offwind (requires changes in modify_prenetwork!)

    dc_links = n.links[
        (n.links.carrier == "DC")
        & (n.links.bus0 + n.links.bus1).str.contains(region)
        & ~n.links.reversed
    ]
    ac_lines = n.lines[(n.lines.bus0 + n.lines.bus1).str.contains(region)]

    # Count length of internationl links only half
    dc_links.loc[
        ~(dc_links.bus0.str.contains(region) & dc_links.bus1.str.contains(region)),
        "length",
    ] *= 0.5
    ac_lines.loc[
        ~(ac_lines.bus0.str.contains(region) & ac_lines.bus1.str.contains(region)),
        "length",
    ] *= 0.5

    # NEP subsets
    current_year = n.generators.build_year.max()  # noqa
    nep_dc = dc_links.query(
        "(index.str.startswith('DC') or index.str.startswith('TYNDP')) and build_year > 2025 and (@current_year - 5 < build_year <= @current_year)"
    ).index
    nep_ac = ac_lines.query(
        "(build_year > 2025) and (@current_year - 5 < build_year <= @current_year)"
    ).index
    var["Capacity|Electricity|Transmission|DC"] = (
        dc_links.eval("p_nom_opt * length").sum() * MW2GW
    )
    var["Capacity|Electricity|Transmission|DC|NEP"] = (
        dc_links.loc[nep_dc].eval("p_nom_opt * length").sum() * MW2GW
    )
    var["Capacity Additions|Electricity|Transmission|DC"] = (
        dc_links.eval("(p_nom_opt - p_nom_min) * length").sum() * MW2GW
    )
    var["Capacity Additions|Electricity|Transmission|DC|NEP"] = (
        dc_links.loc[nep_dc].eval("(p_nom_opt - p_nom_min) * length").sum() * MW2GW
    )
    var["Length Additions|Electricity|Transmission|DC"] = (
        dc_links.eval("p_nom_opt - p_nom_min")
        .floordiv(
            snakemake.params.post_discretization["link_unit_size"]["DC"]
            - 5  # To account for numerical errors subtract a small capacity
        )
        .div(2000 // (snakemake.params.post_discretization["link_unit_size"]["DC"] - 5))
        .multiply(dc_links.length)
        .sum()
    )
    var["Length Additions|Electricity|Transmission|DC|NEP"] = (
        dc_links.loc[nep_dc]
        .eval("p_nom_opt - p_nom_min")
        .floordiv(snakemake.params.post_discretization["link_unit_size"]["DC"] - 5)
        .div(2000 // (snakemake.params.post_discretization["link_unit_size"]["DC"] - 5))
        .multiply(dc_links.length)
        .sum()
    )
    var["Capacity|Electricity|Transmission|AC"] = (
        ac_lines.eval("s_nom_opt * length").sum() * MW2GW
    )
    var["Capacity|Electricity|Transmission|AC|NEP"] = (
        ac_lines.loc[nep_ac].eval("s_nom_opt * length").sum() * MW2GW
    )
    var["Capacity Additions|Electricity|Transmission|AC"] = (
        ac_lines.eval("(s_nom_opt - s_nom_min) * length").sum() * MW2GW
    )
    var["Capacity Additions|Electricity|Transmission|AC|NEP"] = (
        ac_lines.loc[nep_ac].eval("(s_nom_opt - s_nom_min) * length").sum() * MW2GW
    )
    var["Length Additions|Electricity|Transmission|AC"] = (
        ac_lines.eval("s_nom_opt - s_nom_min")
        .floordiv(
            snakemake.params.post_discretization["line_unit_size"] - 5
        )  # To account for numerical errors subtract a small capacity
        .div(5265 // (snakemake.params.post_discretization["line_unit_size"] - 5))
        .multiply(ac_lines.length)
        .sum()
    )
    var["Length Additions|Electricity|Transmission|AC|NEP"] = (
        ac_lines.loc[nep_ac]
        .eval("s_nom_opt - s_nom_min")
        .floordiv(snakemake.params.post_discretization["line_unit_size"] - 5)
        .div(5265 // (snakemake.params.post_discretization["line_unit_size"] - 5))
        .multiply(ac_lines.length)
        .sum()
    )

    var["Capacity|Electricity|Transmission"] = (
        var["Capacity|Electricity|Transmission|AC"]
        + var["Capacity|Electricity|Transmission|DC"]
    )
    var["Capacity|Electricity|Transmission|NEP"] = (
        var["Capacity|Electricity|Transmission|AC|NEP"]
        + var["Capacity|Electricity|Transmission|DC|NEP"]
    )
    var["Capacity Additions|Electricity|Transmission"] = (
        var["Capacity Additions|Electricity|Transmission|AC"]
        + var["Capacity Additions|Electricity|Transmission|DC"]
    )
    var["Capacity Additions|Electricity|Transmission|NEP"] = (
        var["Capacity Additions|Electricity|Transmission|AC|NEP"]
        + var["Capacity Additions|Electricity|Transmission|DC|NEP"]
    )
    var["Length Additions|Electricity|Transmission"] = (
        var["Length Additions|Electricity|Transmission|AC"]
        + var["Length Additions|Electricity|Transmission|DC"]
    )
    var["Length Additions|Electricity|Transmission|NEP"] = (
        var["Length Additions|Electricity|Transmission|AC|NEP"]
        + var["Length Additions|Electricity|Transmission|DC|NEP"]
    )

    ## Distribution Grid
    distr_grid = n.links[
        (n.links.carrier == "electricity distribution grid")
        & (n.links.bus0.str[:2] == region)
        & ~(n.links.index.str.contains("reversed"))
    ]

    var["Capacity|Electricity|Distribution"] = distr_grid.p_nom_opt.sum() * MW2GW
    var["Capacity Additions|Electricity|Distribution"] = (
        distr_grid.eval("(p_nom_opt - p_nom_min)").sum() * MW2GW
    )

    # Hydrogen : TW*km
    # TODO: add missing variables and make nice plot

    h2_links = n.links[
        n.links.carrier.str.contains("H2 pipeline")
        & ~n.links.reversed
        & (n.links.bus0 + n.links.bus1).str.contains(region)
    ]

    # Count length of internationl links according to domestic length factor
    if len(h2_links.carrier.unique()) == 1:
        dlf = domestic_length_factor(n, h2_links.carrier.unique().tolist(), region)
    else:
        dlf = np.array(
            list(
                domestic_length_factor(
                    n, h2_links.carrier.unique().tolist(), region
                ).values()
            )
        ).mean()
    h2_links.loc[
        ~(h2_links.bus0.str.contains(region) & h2_links.bus1.str.contains(region)),
        "length",
    ] *= dlf

    # Kernnetz
    h2_links_kern = h2_links[h2_links.index.str.contains("kernnetz")]

    # Endogeneous
    endo_carriers = ["H2 pipeline", "H2 pipeline retrofitted"]
    h2_links_endo = h2_links[h2_links.carrier.isin(endo_carriers)]

    var["Capacity|Hydrogen|Transmission"] = (
        h2_links.eval("p_nom_opt * length").sum() * MW2TW
    )
    var["Capacity|Hydrogen|Transmission|Kernnetz"] = (
        h2_links_kern.eval("p_nom_opt * length").sum() * MW2TW
    )
    # var["Capacity|Hydrogen|Transmission|Kernnetz|Newbuild"] =
    # var["Capacity|Hydrogen|Transmission|Kernnetz|Retrofitted"] =
    var["Capacity|Hydrogen|Transmission|Endogenous"] = (
        h2_links_endo.eval("p_nom_opt * length").sum() * MW2TW
    )
    # var["Capacity|Hydrogen|Transmission|Endogenous|Newbuild"] =
    # var["Capacity|Hydrogen|Transmission|Endogenous|Retrofitted"] =

    assert isclose(
        var["Capacity|Hydrogen|Transmission"],
        var["Capacity|Hydrogen|Transmission|Kernnetz"]
        + var["Capacity|Hydrogen|Transmission|Endogenous"],
    ), (
        "Hydrogen transmission capacity is not correctly split into Kernnetz and Endogenous"
    )

    year = h2_links.build_year.max()
    new_h2_links = h2_links[
        ((year - 5) < h2_links.build_year) & (h2_links.build_year <= year)
    ]
    new_h2_links_kern = new_h2_links[new_h2_links.index.str.contains("kernnetz")]
    new_h2_links_endo = new_h2_links[new_h2_links.carrier.isin(endo_carriers)]

    var["Capacity Additions|Hydrogen|Transmission"] = (
        new_h2_links.eval("(p_nom_opt - p_nom_min) * length").sum() * MW2TW
    )
    var["Capacity Additions|Hydrogen|Transmission|Kernnetz"] = (
        new_h2_links_kern.eval("(p_nom_opt - p_nom_min) * length").sum() * MW2TW
    )
    # var["Capacity Additions|Hydrogen|Transmission|Kernnetz|Newbuild"] =
    # var["Capacity Additions|Hydrogen|Transmission|Kernnetz|Retrofitted"] =
    var["Capacity Additions|Hydrogen|Transmission|Endogenous"] = (
        new_h2_links_endo.eval("(p_nom_opt - p_nom_min) * length").sum() * MW2TW
    )
    # var["Capacity Additions|Hydrogen|Transmission|Endogenous|Newbuild"] =
    # var["Capacity Additions|Hydrogen|Transmission|Endogenous|Retrofitted"] =

    assert isclose(
        var["Capacity Additions|Hydrogen|Transmission"],
        var["Capacity Additions|Hydrogen|Transmission|Kernnetz"]
        + var["Capacity Additions|Hydrogen|Transmission|Endogenous"],
    ), (
        "Hydrogen transmission capacity additions are not correctly split into Kernnetz and Endogenous"
    )

    # TODO: add length additions

    return var


def hack_DC_projects(n, p_nom_start, p_nom_planned, model_year, snakemake, costs):
    logger.info(f"Hacking DC projects for year {model_year}")

    logger.info("Assuming all indices of DC projects start with 'DC' or 'TYNDP'")
    tprojs = n.links.loc[
        (n.links.index.str.startswith("DC") | n.links.index.str.startswith("TYNDP"))
        & ~n.links.reversed
    ].index

    future_projects = tprojs[n.links.loc[tprojs, "build_year"] > model_year]

    current_projects = tprojs[
        (n.links.loc[tprojs, "build_year"] > (model_year - 5))
        & (n.links.loc[tprojs, "build_year"] <= model_year)
    ]
    past_projects = tprojs[n.links.loc[tprojs, "build_year"] <= (model_year - 5)]

    for proj in tprojs:
        if not isclose(n.links.loc[proj, "p_nom"], p_nom_planned[proj]):
            logger.warning(
                f"(Post-)discretization changed p_nom of {proj} from {p_nom_planned[proj]} to {n.links.loc[proj, 'p_nom']}"
            )

    # Future projects should not have any capacity
    assert isclose(n.links.loc[future_projects, "p_nom_opt"], 0).all()

    # Setting p_nom to 0 such that n.statistics does not compute negative expanded capex or capacity additions
    # Setting p_nom_min to 0 for the grid_expansion calculation
    # This is ONLY POSSIBLE IN POST-PROCESSING
    # We pretend that the model expanded the grid endogenously
    n.links.loc[future_projects, "p_nom"] = 0
    n.links.loc[future_projects, "p_nom_min"] = 0

    # Current projects should have their p_nom_opt bigger or equal to p_nom until the year 2030 (Startnetz that we force in)
    # TODO 2030 is hard coded but should be read from snakemake config
    if model_year <= 2030:
        assert (
            n.links.loc[current_projects, "p_nom_opt"] + 0.1
            >= n.links.loc[current_projects, "p_nom"]
        ).all()

        n.links.loc[current_projects, "p_nom"] -= p_nom_start[current_projects]
        n.links.loc[current_projects, "p_nom_min"] -= p_nom_start[current_projects]
        if (snakemake.params.NEP_year == 2021) or (
            snakemake.params.NEP_transmission == "overhead"
        ):
            logger.warning("Switching DC projects to NEP23 and underground costs.")
            n.links.loc[current_projects, "overnight_cost"] = (
                n.links.loc[current_projects, "length"]
                * (
                    (1.0 - n.links.loc[current_projects, "underwater_fraction"])
                    * costs.at["HVDC underground", "investment"]
                    / 1e-9
                    + n.links.loc[current_projects, "underwater_fraction"]
                    * costs.at["HVDC submarine", "investment"]
                    / 1e-9
                )
                + costs.at["HVDC inverter pair", "investment"] / 1e-9
            )
            n.links.loc[current_projects, "capital_cost"] = (
                n.links.loc[current_projects, "length"]
                * (
                    (1.0 - n.links.loc[current_projects, "underwater_fraction"])
                    * costs.at["HVDC underground", "capital_cost"]
                    / 1e-9
                    + n.links.loc[current_projects, "underwater_fraction"]
                    * costs.at["HVDC submarine", "capital_cost"]
                    / 1e-9
                )
                + costs.at["HVDC inverter pair", "capital_cost"] / 1e-9
            )
    else:
        n.links.loc[current_projects, "p_nom"] = n.links.loc[
            current_projects, "p_nom_min"
        ]

    # Past projects should have their p_nom_opt bigger or equal to p_nom
    if model_year <= 2030 + 5:
        assert (
            n.links.loc[past_projects, "p_nom_opt"] + 0.1  # numerical error tolerance
            >= n.links.loc[past_projects, "p_nom"]
        ).all()

    return n


def hack_AC_projects(n, s_nom_start, model_year, snakemake):
    logger.info(f"Hacking AC projects for year {model_year}")

    current_projects = n.lines.query(
        "@model_year - 5 < build_year <= @model_year"
    ).index

    if snakemake.params.NEP_year == 2021:
        logger.warning("Switching AC projects to NEP23 costs post-optimization")
        n.lines.loc[current_projects, "overnight_cost"] *= 772 / 472

    # Even though the lines are available to the model from the start,
    # we pretend that the lines were expanded in the current year
    # s_nom_start is used, because the model may expand lines
    # endogenously before that or after that
    n.lines.loc[current_projects, "s_nom"] -= s_nom_start[current_projects]
    n.lines.loc[current_projects, "s_nom_min"] -= s_nom_start[current_projects]

    return n


def process_postnetworks(n, n_start, model_year, snakemake, costs):
    post_discretization = snakemake.params.post_discretization

    # logger.info("Post-Discretizing H2 pipeline")

    # assert post_discretization["link_unit_size"]["H2 pipeline"] == post_discretization["link_unit_size"]["H2 pipeline retrofitted"]
    # assert post_discretization["link_threshold"]["H2 pipeline"] == post_discretization["link_threshold"]["H2 pipeline retrofitted"]

    # _h2_lambda = lambda x: get_discretized_value(
    #     x,
    #     post_discretization["link_unit_size"]["H2 pipeline"],
    #     post_discretization["link_threshold"]["H2 pipeline"],
    # )
    # h2_links = n.links.query("carrier == 'H2 pipeline' or carrier == 'H2 pipeline retrofitted'").index
    # for attr in ["p_nom_opt", "p_nom", "p_nom_min"]:
    #     # The values  in p_nom_opt may already be discretized, here we make sure that
    #     # the same logic is applied to p_nom and p_nom_min
    #     n.links.loc[h2_links, attr] = n.links.loc[h2_links, attr].apply(_h2_lambda)
    logger.info("Adding average Kernnetz cost to carrier H2 pipeline (Kernnetz)")
    h2_links_kern = n.links.query("carrier == 'H2 pipeline (Kernnetz))'").index
    capital_costs = (
        0.7 * costs.at["H2 (g) pipeline", "capital_cost"]
        + 0.3 * costs.at["H2 (g) pipeline repurposed", "capital_cost"]
    ) * n.links.loc[h2_links_kern, "length"]
    overnight_costs = (
        0.7 * costs.at["H2 (g) pipeline", "investment"]
        + 0.3 * costs.at["H2 (g) pipeline repurposed", "investment"]
    ) * n.links.loc[h2_links_kern, "length"]
    n.links.loc[h2_links_kern, "capital_cost"] = capital_costs
    n.links.loc[h2_links_kern, "overnight_cost"] = overnight_costs

    logger.info("Post-Discretizing DC links")

    def _dc_lambda(x):
        return get_discretized_value(
            x,
            post_discretization["link_unit_size"]["DC"],
            post_discretization["link_threshold"]["DC"],
        )

    dc_links = n.links.query("carrier == 'DC'").index
    for attr in ["p_nom_opt", "p_nom", "p_nom_min"]:
        # The values  in p_nom_opt may already be discretized, here we make sure that
        # the same logic is applied to p_nom and p_nom_min
        n.links.loc[dc_links, attr] = n.links.loc[dc_links, attr].apply(_dc_lambda)

    p_nom_planned = n_start.links.loc[dc_links, "p_nom"]
    p_nom_start = n_start.links.loc[dc_links, "p_nom"].apply(_dc_lambda)

    logger.info("Post-Discretizing AC lines")

    def _ac_lambda(x):
        return get_discretized_value(
            x,
            post_discretization["line_unit_size"],
            post_discretization["line_threshold"],
        )

    for attr in ["s_nom_opt", "s_nom", "s_nom_min"]:
        # The values  in s_nom_opt may already be discretized, here we make sure that
        # the same logic is applied to s_nom and s_nom_min
        n.lines[attr] = n.lines[attr].apply(_ac_lambda)

    s_nom_start = n_start.lines["s_nom"].apply(_ac_lambda)

    n = hack_DC_projects(n, p_nom_start, p_nom_planned, model_year, snakemake, costs)
    n = hack_AC_projects(n, s_nom_start, model_year, snakemake)
    return n


def get_ariadne_var(
    n,
    industry_demand,
    energy_totals,
    sector_ratios,
    industry_production,
    costs,
    region,
    year,
    config,
    config_industry,
):
    var = pd.concat(
        [
            get_capacities(n, region),
            get_grid_capacity(n, region, year),
            # get_capacity_additions_simple(n,region),
            # get_installed_capacities(n,region),
            get_capacity_additions(n, region),
            # get_capacity_additions_nstat(n, region),
            get_production(region, year),
            get_primary_energy(n, region),
            get_secondary_energy(n, region, industry_demand),
            get_final_energy(
                n,
                region,
                industry_demand,
                energy_totals,
                sector_ratios,
                industry_production,
                config,
                config_industry,
            ),
            get_prices(n, region),
            get_emissions(n, region, energy_totals, industry_demand),
            get_policy(n, year),
            get_trade(n, region),
            get_economy(n, region),
            get_system_cost(n, region),
        ]
    )

    return var


# uses the global variables model, scenario and var2unit. For now.
def get_data(
    n,
    industry_demand,
    energy_totals,
    sector_ratios,
    industry_production,
    costs,
    region,
    year,
    config,
    config_industry,
    version="0.10",
    scenario="test",
):
    var = get_ariadne_var(
        n,
        industry_demand,
        energy_totals,
        sector_ratios,
        industry_production,
        costs,
        region,
        year,
        config,
        config_industry,
    )

    # Renaming variables

    var["Investment|Energy Supply|Electricity|Wind Onshore"] = var[
        "Investment|Energy Supply|Electricity|Wind|Onshore"
    ]

    var["Investment|Energy Supply|Electricity|Wind Offshore"] = var[
        "Investment|Energy Supply|Electricity|Wind|Offshore"
    ]

    var["Investment|Energy Supply|Electricity|Electricity Storage"] = var[
        "Investment|Energy Supply|Electricity|Storage Reservoir"
    ]

    var["Investment|Energy Supply|Heat|Heatpump"] = var[
        "Investment|Energy Supply|Heat|Heat pump"
    ]

    var["Investment|Energy Supply|Heat|Solarthermal"] = var[
        "Investment|Energy Supply|Heat|Solar thermal"
    ]

    var["Investment|Energy Supply|Hydrogen|Storage"] = var[
        "Investment|Energy Supply|Hydrogen|Reservoir"
    ]

    var["Investment|Energy Supply|Hydrogen|Electrolysis"] = var[
        "Investment|Energy Supply|Hydrogen|Electricity"
    ]

    var["Investment|Energy Supply|Hydrogen|Fossil"] = var[
        "Investment|Energy Supply|Hydrogen|Gas"
    ]
    # For internal use only and translated directly to TWh
    var["Demand|Electricity"] = var.reindex(
        [
            "Secondary Energy|Electricity|Storage Losses",
            "Secondary Energy|Electricity|Transmission Losses",
            "Secondary Energy Input|Electricity|Heat",
            "Secondary Energy Input|Electricity|Hydrogen",
            "Secondary Energy Input|Electricity|Liquids",
            "Final Energy|Industry|Electricity",
            "Final Energy|Agriculture|Electricity",
            "Final Energy|Residential and Commercial|Electricity",
            "Final Energy|Transportation|Electricity",
            "Final Energy|Carbon Dioxide Removal|Electricity",
        ]
    ).sum()

    data = []
    for v in var.index:
        try:
            unit = var2unit[v]
        except KeyError:
            unit = "NA"

        data.append(
            [
                "PyPSA-Eur " + version,
                scenario,
                region,
                v,
                unit,
                var[v],
            ]
        )

    tab = pd.DataFrame(
        data, columns=["Model", "Scenario", "Region", "Variable", "Unit", year]
    )

    return tab


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_ariadne_variables",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    config_industry = snakemake.params.config_industry
    planning_horizons = snakemake.params.planning_horizons
    post_discretization = snakemake.params.post_discretization
    ariadne_template = pd.read_excel(snakemake.input.template, sheet_name=None)
    var2unit = ariadne_template["variable_definitions"].set_index("Variable")["Unit"]
    var2unit = var2unit.str.replace("GJ", "MWh")
    var2unit = var2unit.str.replace("PJ", "TWh")

    industry_demands = [
        pd.read_csv(
            in_dem,
            index_col="TWh/a (MtCO2/a)",
        )
        .multiply(TWh2MWh)
        .rename_axis("bus")
        for in_dem in snakemake.input.industry_demands
    ]
    energy_totals = (
        pd.read_csv(
            snakemake.input.energy_totals,
            index_col=[0, 1],
        )
        .xs(
            snakemake.params.energy_totals_year,
            level="year",
        )
        .multiply(TWh2MWh)
    )

    sector_ratios = [
        pd.read_csv(
            in_sec_ratio,
            header=[0, 1],
            index_col=0,
        ).rename_axis("carrier")
        for in_sec_ratio in snakemake.input.industry_sector_ratios
    ]
    industry_production = [
        pd.read_csv(
            in_ind_prod,
            index_col="kton/a",
        ).rename_axis("country")
        for in_ind_prod in snakemake.input.industrial_production_per_country_tomorrow
    ]

    # Load data
    _networks = [pypsa.Network(fn) for fn in snakemake.input.networks]

    costs = list(
        map(
            lambda _costs: load_costs(_costs).multiply(1e-9),  # in bn €
            snakemake.input.costs,
        )
    )
    modelyears = [fn[-7:-3] for fn in snakemake.input.networks]
    # Hack the transmission projects
    networks = [
        process_postnetworks(n.copy(), _networks[0], int(my), snakemake, c)
        for n, my, c in zip(_networks, modelyears, costs)
    ]

    if "debug" == "debug":  # For debugging
        var = pd.Series()
        idx = 1
        n = networks[idx]
        c = costs[idx]
        _industry_demand = industry_demands[idx]
        industry_demand = industry_demands[idx]
        _industry_production = industry_production[idx]
        _sector_ratios = sector_ratios[idx]
        _energy_totals = energy_totals.copy()
        region = "DE"
        cap_func = n.statistics.optimal_capacity
        cap_string = "Optimal Capacity|"
        kwargs = {
            "groupby": ["bus", "carrier"],
            "at_port": True,
            "nice_names": False,
        }

    yearly_dfs = []
    for i, year in enumerate(planning_horizons):
        logger.info(f"Getting data for year {year}...")
        yearly_dfs.append(
            get_data(
                networks[i],
                industry_demands[i],
                energy_totals,
                sector_ratios[i],
                industry_production[i],
                costs[i],
                "DE",
                year=year,
                config=snakemake.config,
                config_industry=snakemake.params.config_industry,
                version=snakemake.config["version"],
                scenario=snakemake.wildcards.run,
            )
        )

    df = reduce(
        lambda left, right: pd.merge(
            left, right, on=["Model", "Scenario", "Region", "Variable", "Unit"]
        ),
        yearly_dfs,
    )

    logger.info("Gleichschaltung of AC-Startnetz with investments for AC projects")
    # In this hacky part of the code we assure that the investments for the AC projects, match those of the NEP-AC-Startnetz
    # Thus the variable 'Investment|Energy Supply|Electricity|Transmission|AC' is equal to the sum of exogeneous AC projects, endogenous AC expansion and Übernahme of NEP costs (mainly Systemdienstleistungen (Reactive Power Compensation) and lines that are below our spatial resolution)
    ac_startnetz = 14.5 / 5 / EUR20TOEUR23  # billion EUR

    ac_projects_invest = df.query(
        "Variable == 'Investment|Energy Supply|Electricity|Transmission|AC|NEP|Onshore'"
    )[planning_horizons].values.sum()

    df.loc[
        df.query(
            "Variable == 'Investment|Energy Supply|Electricity|Transmission|AC|Übernahme|Startnetz Delta'"
        ).index,
        [2025, 2030, 2035, 2040],
    ] += (ac_startnetz - ac_projects_invest) / 4

    for suffix in ["|AC|NEP", "|AC", "", " and Distribution"]:
        df.loc[
            df.query(
                f"Variable == 'Investment|Energy Supply|Electricity|Transmission{suffix}'"
            ).index,
            [2025, 2030, 2035, 2040],
        ] += (ac_startnetz - ac_projects_invest) / 4

    logger.info("Assigning mean investments of year and year + 5 to year.")
    investment_rows = df.loc[df["Variable"].str.contains("Investment")]
    average_investments = (
        investment_rows[planning_horizons]
        .add(investment_rows[planning_horizons].shift(-1, axis=1))
        .div(2)
        .fillna(0.0)
    )
    average_investments[planning_horizons[0]] += investment_rows[
        planning_horizons[0]
    ].div(2)
    average_investments[planning_horizons[-1]] += investment_rows[
        planning_horizons[-1]
    ].div(2)
    df.loc[investment_rows.index, planning_horizons] = average_investments

    df["Region"] = df["Region"].str.replace("DE", "DEU")
    df["Model"] = "PyPSA-Eur v0.10"

    with pd.ExcelWriter(snakemake.output.exported_variables_full) as writer:
        df.round(5).to_excel(writer, sheet_name="data", index=False)

    logger.info("Dropping variables which are not in the template.")
    ariadne_df = df.drop(df.loc[df["Unit"] == "NA"].index)

    meta = pd.Series(
        {
            "Model": "PyPSA-Eur v0.10",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "preliminary",
            "Internal usage within Kopernikus AG Szenarien": "yes",
            "Release for publication": "no",
        }
    )

    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        ariadne_df.round(5).to_excel(writer, sheet_name="data", index=False)
        meta.to_frame().T.to_excel(writer, sheet_name="meta", index=False)
