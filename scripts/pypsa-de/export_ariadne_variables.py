# -*- coding: utf-8 -*-
import math
import os
import re
import sys
from functools import reduce

import numpy as np
import pandas as pd
import pypsa
from numpy import isclose

paths = [
    "workflow/submodules/pypsa-eur/scripts",
    "../submodules/pypsa-eur/scripts",
    "../submodules/pypsa-eur/",
]
for path in paths:
    sys.path.insert(0, os.path.abspath(path))

from _helpers import mute_print
from prepare_sector_network import prepare_costs

# Defining global varibales

TWh2PJ = 3.6
MWh2TJ = 3.6e-3
MW2GW = 1e-3
t2Mt = 1e-6

MWh2GJ = 3.6
TWh2PJ = 3.6
MWh2PJ = 3.6e-6
toe_to_MWh = 11.630  # GWh/ktoe OR MWh/toe


def _get_oil_fossil_fraction(n, region):
    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "at_port": True,
        "nice_names": False,
    }

    if "DE" in region:
        total_oil_supply = (
            n.statistics.supply(bus_carrier="oil", **kwargs)
            .groupby("name")
            .sum()
            .reindex(
                [
                    "DE oil refining",
                    "DE renewable oil -> DE oil",
                    "EU renewable oil -> DE oil",
                ]
            )
            .dropna()
        )  # If links are not used they are dropped here

    else:
        total_oil_supply = (
            n.statistics.supply(bus_carrier="oil", **kwargs)
            .groupby("name")
            .sum()
            .reindex(
                [
                    "EU oil refining",
                    "DE renewable oil -> EU oil",
                    "EU renewable oil -> EU oil",
                ]
            )
            .dropna()
        )

    oil_fossil_fraction = (
        total_oil_supply[~total_oil_supply.index.str.contains("renewable")].sum()
        / total_oil_supply.sum()
    )

    return oil_fossil_fraction


def _get_gas_fractions(n, region):
    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "at_port": True,
        "nice_names": False,
    }

    renewable_gas_supply = (
        n.statistics.supply(bus_carrier="renewable gas", **kwargs)
        .groupby(["bus", "carrier"])
        .sum()
    )

    total_gas_supply = (
        n.statistics.supply(bus_carrier="gas", **kwargs)
        .groupby(["name", "carrier"])
        .sum()
    )

    if "DE" in region:
        domestic_gas_supply = (
            total_gas_supply.reindex(
                [
                    ("DE gas", "gas"),
                    ("DE renewable gas -> DE gas", "renewable gas"),
                    ("EU renewable gas -> DE gas", "renewable gas"),
                ]
            )
            .dropna()
            .groupby("carrier")
            .sum()
        )  # If links are not used they are dropped here
        total_imported_renewable_gas = total_gas_supply.get(
            ("EU renewable gas -> DE gas", "renewable gas"), 0
        )
        total_exported_renewable_gas = total_gas_supply.get(
            ("DE renewable gas -> EU gas", "renewable gas"), 0
        )
        domestic_renewable_gas = renewable_gas_supply.get("DE renewable gas")
        foreign_renewable_gas = renewable_gas_supply.get("EU renewable gas")
    else:
        domestic_gas_supply = (
            total_gas_supply.reindex(
                [
                    ("EU gas", "gas"),
                    ("DE renewable gas -> EU gas", "renewable gas"),
                    ("EU renewable gas -> EU gas", "renewable gas"),
                ]
            )
            .dropna()
            .groupby("carrier")
            .sum()
        )
        total_imported_renewable_gas = total_gas_supply.get(
            ("DE renewable gas -> EU gas", "renewable gas"), 0
        )
        total_exported_renewable_gas = total_gas_supply.get(
            ("EU renewable gas -> DE gas", "renewable gas"), 0
        )
        domestic_renewable_gas = renewable_gas_supply.get("EU renewable gas", 0)
        foreign_renewable_gas = renewable_gas_supply.get("DE renewable gas", 0)

    # Legacy code for dropping gas pipelines (now deactivated)
    drops = ["gas pipeline", "gas pipeline new"]
    for d in drops:
        if d in total_gas_supply.index:
            total_gas_supply.drop(d, inplace=True)

    if total_imported_renewable_gas == 0:
        imported_renewable_gas = pd.Series(
            0, index=renewable_gas_supply.index.get_level_values("carrier").unique()
        )
    else:
        imported_renewable_gas = (
            foreign_renewable_gas
            / foreign_renewable_gas.sum()
            * total_imported_renewable_gas
        )

    if total_exported_renewable_gas == 0:
        exported_renewable_gas = pd.Series(
            0, index=renewable_gas_supply.index.get_level_values("carrier").unique()
        )
    else:
        exported_renewable_gas = (
            domestic_renewable_gas
            / domestic_renewable_gas.sum()
            * total_exported_renewable_gas
        )

    renewable_gas_balance = imported_renewable_gas.add(
        domestic_renewable_gas, fill_value=0
    ).subtract(exported_renewable_gas, fill_value=0)
    # These numbers may deviate because the input data from n.statistics
    # can deviate. I don't know why exactly, but at least we can check that
    # the difference stays roughly the same after the calculation.
    assert isclose(
        domestic_gas_supply.get("renewable gas", 0) - renewable_gas_balance.sum(),
        total_gas_supply.get(
            ["DE renewable gas -> DE gas"],
            pd.Series(0)).sum() 
        + total_gas_supply.get(
            ["DE renewable gas -> EU gas"],
            pd.Series(0)
        ).sum()
        - renewable_gas_supply.get("DE renewable gas", pd.Series(0)).sum(),
        rtol=1e-3,
    )

    gas_fractions = pd.Series(
        {
            "Natural Gas": domestic_gas_supply.get("gas", 0),
            "Biomass": renewable_gas_balance.filter(like="biogas").sum(),
            "Efuel": renewable_gas_balance.get("Sabatier", 0),
        }
    ).divide(domestic_gas_supply.sum())

    assert isclose(gas_fractions.sum(), 1)

    return gas_fractions


def _get_h2_fossil_fraction(n):
    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "at_port": True,
        "nice_names": False,
    }
    total_h2_supply = (
        n.statistics.supply(bus_carrier="H2", **kwargs)
        .drop("Store", errors="ignore")
        .groupby("carrier")
        .sum()
    )

    h2_fossil_fraction = total_h2_supply.get("SMR") / total_h2_supply.sum()

    return h2_fossil_fraction


def _get_t_sum(df, df_t, carrier, region, snapshot_weightings, port):
    if type(carrier) == list:
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
    return MWh2PJ * _get_t_sum(
        n.loads,
        n.loads_t,
        carrier,
        region,
        n.snapshot_weightings.generators,
        "p",
    )


def sum_co2(n, carrier, region):
    if type(carrier) == list:
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
        print(
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


def get_investments(n, costs, region):
    def _f(**kwargs):
        return n.statistics.expanded_capex(**kwargs, cost_attribute="overnight_cost")

    return _get_capacities(
        n, region, _f, cap_string="Investment|Energy Supply|", costs=costs
    )


def get_capacity_additions_nstat(n, region):
    def _f(*args, **kwargs):
        kwargs.pop("storage", None)
        return n.statistics.expanded_capacity(*args, **kwargs)

    return _get_capacities(n, region, _f, cap_string="Capacity Additions Nstat|")


costs_dict = {
    # capacities electricity
    "BEV charger": None,
    "CCGT": "CCGT",
    "DAC": "direct air capture",
    "H2 Electrolysis": "electrolysis",
    "H2 Fuel Cell": "fuel cell",
    "OCGT": "OCGT",
    "PHS": "PHS",
    "V2G": None,
    "battery charger": "battery inverter",
    "battery discharger": "battery inverter",
    "coal": "coal",
    "gas pipeline": "CH4 (g) pipeline",
    "home battery charger": "home battery inverter",
    "home battery discharger": "home battery inverter",
    "hydro": "hydro",
    "lignite": "lignite",
    "methanolisation": "methanolisation",
    "offwind-ac": "offwind",  # TODO add grid connection cost
    "offwind-dc": "offwind",  # TODO add grid connection cost
    "offwind-float": "offwind-float",  # TODO add grid connection cost
    "oil": "oil",
    "onwind": "onwind",
    "ror": "ror",
    "rural air heat pump": "decentral air-sourced heat pump",
    "rural ground heat pump": "decentral ground-sourced heat pump",
    "rural resistive heater": "decentral resistive heater",
    "rural solar thermal": "decentral solar thermal",
    "solar": "solar-utility",
    "solar rooftop": "solar-rooftop",
    "solar-hsat": "solar-utility single-axis tracking",
    "solid biomass": "central solid biomass CHP",
    "urban central air heat pump": "central air-sourced heat pump",
    "urban central coal CHP": "central coal CHP",
    "urban central gas CHP": "central gas CHP",
    "urban central gas CHP CC": "central gas CHP",
    "urban central lignite CHP": "central coal CHP",
    "urban central oil CHP": "central gas CHP",
    "urban central resistive heater": "central resistive heater",
    "urban central solar thermal": "central solar thermal",
    "urban central solid biomass CHP": "central solid biomass CHP",
    "urban central solid biomass CHP CC": "central solid biomass CHP CC",
    "urban decentral air heat pump": "decentral air-sourced heat pump",
    "urban decentral resistive heater": "decentral resistive heater",
    "urban decentral solar thermal": "decentral solar thermal",
    "waste CHP": "waste CHP",
    "waste CHP CC": "waste CHP CC",
    # Heat capacities
    # TODO Check the units of the investments
    "DAC": "direct air capture",
    "Fischer-Tropsch": None,  #'Fischer-Tropsch' * "efficiency" ,
    "H2 Electrolysis": "electrolysis",
    "H2 Fuel Cell": "fuel cell",
    "Sabatier": "methanation",
    "methanolisation": "methanolisation",
    # 'urban central air heat pump': 'central air-sourced heat pump',
    # 'urban central coal CHP': 'central coal CHP',
    # 'urban central gas CHP': 'central gas CHP',
    # 'urban central gas CHP CC': 'central gas CHP',
    # 'urban central lignite CHP': 'central coal CHP',
    # 'urban central oil CHP': 'central gas CHP',
    # 'urban central resistive heater': 'central resistive heater',
    # 'urban central solid biomass CHP': 'central solid biomass CHP',
    # 'urban central solid biomass CHP CC': 'central solid biomass CHP CC',
    "urban central water tanks charger": "water tank charger",
    "urban central water tanks discharger": "water tank discharger",
    # 'waste CHP': 'waste CHP',
    # 'waste CHP CC': 'waste CHP CC',
    #  Decentral Heat capacities
    # TODO consider overdim_factor
    #'rural air heat pump': None,
    "rural biomass boiler": "biomass boiler",
    "rural gas boiler": "decentral gas boiler",
    #'rural ground heat pump': None,
    "rural oil boiler": "decentral oil boiler",
    "rural resistive heater": "decentral resistive heater",
    "rural water tanks charger": "water tank charger",
    "rural water tanks discharger": "water tank discharger",
    #'urban decentral air heat pump': None,
    "urban decentral biomass boiler": "biomass boiler",
    "urban decentral gas boiler": "decentral gas boiler",
    "urban decentral oil boiler": "decentral oil boiler",
    "urban decentral resistive heater": "decentral resistive heater",
    "urban decentral water tanks charger": "water tank charger",
    "urban decentral water tanks discharger": "water tank discharger",
    # Other capacities
    # 'Sabatier': 'methanation',costs.at["methanation", "fixed"]
    # * costs.at["methanation", "efficiency"]
    "biogas to gas": None,  # TODO biogas + biogas upgrading
    "biogas to gas CC": None,  # TODO costs.at["biogas CC", "fixed"]
    # + costs.at["biogas upgrading", "fixed"]
    # + costs.at["biomass CHP capture", "fixed"]
    # * costs.at["biogas CC", "CO2 stored"],
}


def _get_capacities(n, region, cap_func, cap_string="Capacity|", costs=None):

    kwargs = {
        "groupby": n.statistics.groupers.get_bus_and_carrier,
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

    # if cap_string.startswith("Investment"):
    #     technology_investments = pd.Series(
    #         [
    #             0 if costs_dict.get(key) is None else
    #             costs.at[costs_dict.get(key), "investment"]
    #             for key in capacities_electricity.index
    #         ],
    #         capacities_electricity.index
    #     )
    #     for carrier in technology_investments.index.intersection(
    #         ["onwind", "solar", "solar-hsat"]):
    #         technology_investments[carrier] += \
    #             costs.at["electricity grid connection", "investment"]

    #     for carrier in ["offwind-ac", "offwind-dc", "offwind-float"]:
    #         technology_investments[carrier] = 0.0 # TODO add grid connection cost
    #     #     apply update_wind_solar_costs(n, costs)

    #     capacities_electricity = \
    #         capacities_electricity.div(5).multiply(technology_investments)
    # else:
    #     capacities_electricity = \
    #         capacities_electricity.multiply(MW2GW)

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

    var[cap_string + "Electricity|Gas|CC"] = capacities_electricity.get("CCGT")

    var[cap_string + "Electricity|Gas|OC"] = capacities_electricity.get("OCGT")

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
    var[cap_string + "Electricity|Hydrogen"] = var[
        cap_string + "Electricity|Hydrogen|FC"
    ] = capacities_electricity.reindex(
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

    var[cap_string + "Electricity|Storage Converter|Stationary Batteries"] = (
        capacities_electricity.get("battery discharger", 0)
        + capacities_electricity.get("home battery discharger", 0)
    )

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

    if cap_string.startswith("Investment"):
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

    var[cap_string + "Electricity|Storage Reservoir|Stationary Batteries"] = pd.Series(
        {c: storage_capacities.get(c) for c in ["battery", "home battery"]}
    ).sum()

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

    # Test if we forgot something
    #
    # Unconvenient at the moment, requires further changes to n.statistics
    #
    # _drop_idx = [
    #     col for col in [
    #         "PHS",
    #         "battery discharger",
    #         "home battery discharger",
    #         "V2G",
    #     ] if col in capacities_electricity.index
    # ]
    # assert isclose(
    #     var[cap_string + "Electricity"],
    #     capacities_electricity.drop(_drop_idx).sum(),
    # )

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

    var[cap_string + "Heat|Storage Converter"] = capacities_central_heat.filter(
        like="water tanks discharger"
    ).sum()

    var[cap_string + "Heat|Storage Reservoir"] = storage_capacities.filter(
        like="water tanks"
    ).sum()

    var[cap_string + "Heat"] = (
        var[cap_string + "Heat|Solar thermal"]
        + var[cap_string + "Heat|Resistive heater"]
        + var[cap_string + "Heat|Biomass"]
        + var[cap_string + "Heat|Oil"]
        + var[cap_string + "Heat|Gas"]
        + var[cap_string + "Heat|Processes"]
        +
        # var[cap_string + "Heat|Hydrogen"] +
        var[cap_string + "Heat|Heat pump"]
        + var[cap_string + "Heat|Non-Renewable Waste"]
    )

    # This check requires further changes to n.statistics
    # assert isclose(
    #     var[cap_string + "Heat"],
    #     capacities_central_heat[
    #         # exclude storage converters (i.e., dischargers)
    #         ~capacities_central_heat.index.str.contains("discharger|DAC")
    #     ].sum()
    # )

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

    var[cap_string + "Decentral Heat|Solar thermal"] = capacities_decentral_heat.filter(
        like="solar thermal"
    ).sum()

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

    var[cap_string + "Hydrogen|Electricity"] = capacities_h2.get("H2 Electrolysis", 0)

    var[cap_string + "Hydrogen"] = (
        var[cap_string + "Hydrogen|Electricity"] + var[cap_string + "Hydrogen|Gas"]
    )

    # This check requires further changes to n.statistics
    #
    # assert isclose(
    #     var[cap_string + "Hydrogen"],
    #     capacities_h2.reindex([
    #         "H2 Electrolysis",
    #         "SMR",
    #         "SMR CC",
    #     ]).sum(), # if technology not build, reindex returns NaN
    # )

    var[cap_string + "Hydrogen|Reservoir"] = storage_capacities.get("H2", 0)

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
    var[cap_string + "Liquids|Hydrogen"] = var[cap_string + "Liquids"] = (
        capacities_liquids.get("Fischer-Tropsch", 0)
    )

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
        print(
            "Warning: carrier `methanol` not found in network.links.carrier! Assuming 0 capacities."
        )
        var[cap_string + "Methanol"] = 0

    if cap_string.startswith("Investment"):
        var = var.div(MW2GW).mul(1e-9).div(5).round(3)  # in bn € / year
    return var


def get_CHP_E_and_H_usage(n, bus_carrier, region, fossil_fraction=1):
    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }

    usage = (
        n.statistics.withdrawal(
            bus_carrier=bus_carrier,
            **kwargs,
        )
        .filter(like=region)
        .filter(like="CHP")
        .multiply(MWh2PJ)
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
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }

    var = pd.Series()

    oil_fossil_fraction = _get_oil_fossil_fraction(n, region)
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
        .multiply(MWh2PJ)
    )

    oil_CHP_E_usage, oil_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "oil", region, fossil_fraction=oil_fossil_fraction
    )

    ## Primary Energy

    var["Primary Energy|Oil|Heat"] = (
        oil_usage.filter(like="urban central oil boiler").sum() + oil_CHP_H_usage
    ) * primary_oil_factor

    var["Primary Energy|Oil|Electricity"] = (
        oil_usage.get("oil", 0) + oil_CHP_E_usage
    ) * primary_oil_factor

    var["Primary Energy|Oil"] = oil_usage.sum() * primary_oil_factor

    # At the moment, everyting that is not electricity or heat is counted as liquid fuel
    var["Primary Energy|Oil|Liquids"] = (
        var["Primary Energy|Oil"]
        - var["Primary Energy|Oil|Electricity"]
        - var["Primary Energy|Oil|Heat"]
    )

    gas_fractions = _get_gas_fractions(n, region)

    # Eventhough biogas gets routed through the EU gas bus,
    # it should be counted separately as Primary Energy|Biomass
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
        .multiply(MWh2PJ)
    )

    gas_CHP_E_usage, gas_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "gas", region, fossil_fraction=gas_fractions["Natural Gas"]
    )

    var["Primary Energy|Gas|Heat"] = (
        gas_usage.filter(like="urban central gas boiler").sum() + gas_CHP_H_usage
    )

    var["Primary Energy|Gas|Electricity"] = (
        gas_usage.reindex(
            [
                "CCGT",
                "OCGT",
            ],
        ).sum()
        + gas_CHP_E_usage
    )

    var["Primary Energy|Gas|Hydrogen"] = gas_usage.filter(like="SMR").sum()

    var["Primary Energy|Gas"] = gas_usage.sum()

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
        .multiply(MWh2PJ)
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
        .groupby("carrier")
        .sum()
        .multiply(MWh2PJ)
    )

    biomass_CHP_E_usage, biomass_CHP_H_usage = get_CHP_E_and_H_usage(
        n, "solid biomass", region
    )

    var["Primary Energy|Biomass|Gases"] = biomass_usage.filter(
        like="biogas to gas"
    ).sum()

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

    var["Primary Energy|Biomass"] = biomass_usage.sum()

    var["Primary Energy|Nuclear"] = (
        n.statistics.withdrawal(
            bus_carrier=["uranium"],
            **kwargs,
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(MWh2PJ)
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
        .multiply(MWh2PJ)
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
        .multiply(MWh2PJ)
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

    return var


def get_secondary_energy(n, region, _industry_demand):
    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }
    var = pd.Series()

    electricity_supply = (
        n.statistics.supply(bus_carrier=["low voltage", "AC"], **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .multiply(MWh2PJ)
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

    var["Secondary Energy|Electricity|Fossil"] = (
        var["Secondary Energy|Electricity|Gas"]
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
        electricity_supply.get("biogas")
    )
    var["Secondary Energy|Electricity|Biomass"] = (
        var["Secondary Energy|Electricity|Biomass|w/o CCS"]
        + var["Secondary Energy|Electricity|Biomass|w/ CCS"]
    )

    var["Secondary Energy|Electricity|Hydro"] = electricity_supply.get(
        "hydro"
    ) + electricity_supply.get("ror")
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
        .multiply(MWh2PJ)
        .values.sum()
    )

    var["Secondary Energy|Electricity|Storage Losses"] = (
        n.statistics.withdrawal(bus_carrier=["AC", "low voltage"], **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .reindex(
            [
                "BEV charger",
                "battery charger",
                "home battery charger",
                "PHS",
            ]
        )
        .subtract(
            n.statistics.supply(bus_carrier=["AC", "low voltage"], **kwargs)
            .filter(like=region)
            .groupby(["carrier"])
            .sum()
            .reindex(
                [
                    "V2G",
                    "battery discharger",
                    "home battery discharger",
                    "PHS",
                ]
            )
        )
        .multiply(MWh2PJ)
        .sum()
    )

    # TODO Compute transmission losses via links_t
    # var["Secondary Energy|Electricity|Transmission Losses"] = \
    #     n.statistics.withdrawal(
    #         bus_carrier=["AC", "low voltage"], **kwargs
    #     ).filter(like=region).groupby(["carrier"]).sum().get(
    #         ["AC", "DC", "electricity distribution grid"]
    #     ).subtract(
    #         n.statistics.supply(
    #             bus_carrier=["AC", "low voltage"], **kwargs
    #         ).filter(like=region).groupby(["carrier"]).sum().get(
    #             ["AC", "DC", "electricity distribution grid"]
    #         )
    #     ).multiply(MWh2PJ).sum()

    # supply - withdrawal
    # var["Secondary Energy|Electricity|Storage"] = \
    var["Secondary Energy|Electricity"] = (
        var["Secondary Energy|Electricity|Fossil"]
        + var["Secondary Energy|Electricity|Biomass"]
        + var["Secondary Energy|Electricity|Non-Biomass Renewables"]
        + var["Secondary Energy|Electricity|Nuclear"]
        # + var["Secondary Energy|Electricity|Transmission Losses"]
        # + var["Secondary Energy|Electricity|Storage Losses"]
        + var["Secondary Energy|Electricity|Hydrogen"]
        + var["Secondary Energy|Electricity|Waste"]
    )

    assert isclose(
        electricity_supply[
            ~electricity_supply.index.str.contains(
                "PHS" "|battery discharger" "|home battery discharger" "|V2G"
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
        .multiply(MWh2PJ)
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
        .multiply(MWh2PJ)
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
    )

    oil_fossil_fraction = _get_oil_fossil_fraction(n, region)

    oil_fuel_usage = (
        n.statistics.withdrawal(bus_carrier="oil", **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(MWh2PJ)
        .reindex(
            [
                "agriculture machinery oil",
                "kerosene for aviation",
                "land transport oil",
                "naphtha for industry",
                "shipping oil",
                # TODO is decentral heating oil considered a fuel?
                "rural oil boiler",
                "urban decentral oil boiler",
            ]
        )
    )

    total_oil_fuel_usage = oil_fuel_usage.sum()
    var["Secondary Energy|Liquids|Fossil"] = var["Secondary Energy|Liquids|Oil"] = (
        total_oil_fuel_usage * oil_fossil_fraction
    )
    var["Secondary Energy|Liquids|Hydrogen"] = total_oil_fuel_usage * (
        1 - oil_fossil_fraction
    )

    var["Secondary Energy|Liquids"] = (
        var["Secondary Energy|Liquids|Oil"] + var["Secondary Energy|Liquids|Hydrogen"]
    )

    try:
        methanol_production = (
            n.statistics.supply(bus_carrier="methanol", **kwargs)
            .filter(like=region)
            .groupby(["carrier"])
            .sum()
            .multiply(MWh2PJ)
        )

        assert (
            methanol_production.size <= 2
        )  # only methanolisation and methanol imports

        # var["Production|Chemicals|Methanol"] = \ # here units are Mt/year
        var["Secondary Energy|Methanol"] = methanol_production.get("methanolisation", 0)
        # Remeber to specify that Other Carrier == Methanol in Comments Tab
    except KeyError:
        var["Secondary Energy|Methanol"] = 0

    gas_fuel_usage = (
        n.statistics.withdrawal(bus_carrier="gas", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .multiply(MWh2PJ)
        .reindex(
            [
                "gas for industry",
                "gas for industry CC",
                "rural gas boiler",
                "urban decentral gas boiler",
            ]
        )
    )  # Building, Transport and Industry sectors

    total_gas_fuel_usage = gas_fuel_usage.sum()

    gas_fractions = _get_gas_fractions(n, region)

    # Fraction supplied by Hydrogen conversion
    var["Secondary Energy|Gases|Hydrogen"] = total_gas_fuel_usage * gas_fractions.get(
        "Efuel"
    )

    var["Secondary Energy|Gases|Biomass"] = total_gas_fuel_usage * gas_fractions.get(
        "Biomass"
    )

    var["Secondary Energy|Gases|Natural Gas"] = (
        total_gas_fuel_usage * gas_fractions.get("Natural Gas")
    )

    var["Secondary Energy|Gases"] = (
        var["Secondary Energy|Gases|Hydrogen"]
        + var["Secondary Energy|Gases|Biomass"]
        + var["Secondary Energy|Gases|Natural Gas"]
    )

    assert isclose(var["Secondary Energy|Gases"], total_gas_fuel_usage)

    industry_demand = _industry_demand.filter(
        like=region,
        axis=0,
    ).sum()
    mwh_coal_per_mwh_coke = 1.366
    # Coke is added as a coal demand, so we need to convert back to units of coke for secondary energy
    var["Secondary Energy|Solids|Coal"] = var["Secondary Energy|Solids"] = (
        industry_demand.get("coke", 0) / mwh_coal_per_mwh_coke
    )

    electricity_withdrawal = (
        n.statistics.withdrawal(bus_carrier=["low voltage", "AC"], **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .multiply(MWh2PJ)
    )

    var["Secondary Energy Input|Electricity|Hydrogen"] = electricity_withdrawal.get(
        "H2 Electrolysis", 0
    )

    var["Secondary Energy Input|Electricity|Heat"] = electricity_withdrawal.filter(
        like="urban central"
    ).sum()

    hydrogen_withdrawal = (
        n.statistics.withdrawal(bus_carrier="H2", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .multiply(MWh2PJ)
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

    var["Secondary Energy Input|Hydrogen|Liquids"] = hydrogen_withdrawal.get(
        "Fischer-Tropsch", 0
    )

    var["Secondary Energy"] = (
        var["Secondary Energy|Electricity"]
        + var["Secondary Energy|Heat"]
        + var["Secondary Energy|Hydrogen"]
        + var["Secondary Energy|Gases"]
        + var["Secondary Energy|Liquids"]
        + var["Secondary Energy|Solids"]
        + var["Secondary Energy|Methanol"]
    )

    return var


def get_final_energy(
    n, region, _industry_demand, _energy_totals, _sector_ratios, _industry_production
):

    var = pd.Series()

    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }
    h2_fossil_fraction = _get_h2_fossil_fraction(n)
    oil_fossil_fraction = _get_oil_fossil_fraction(n, region)

    if config["industry"]["ammonia"]:
        # MWh/a
        Haber_Bosch_NH3 = (
            n.statistics.supply(bus_carrier="NH3", **kwargs)
            .groupby("carrier")
            .sum()["Haber-Bosch"]
        )

        CH4_for_NH3 = (
            Haber_Bosch_NH3
            * h2_fossil_fraction
            * config["industry"]["MWh_CH4_per_tNH3_SMR"]
            / config["industry"]["MWh_NH3_per_tNH3"]
            * MWh2PJ
        )
        H2_for_NH3 = (
            Haber_Bosch_NH3
            * (1 - h2_fossil_fraction)
            / config["industry"]["MWh_H2_per_tNH3_electrolysis"]
            * MWh2PJ
        )
        subcategories = ["HVC", "Methanol", "Chlorine"]

    else:
        CH4_for_NH3 = 0
        H2_for_NH3 = 0
        subcategories = ["HVC", "Methanol", "Chlorine", "Ammonia"]

    carrier = ["hydrogen", "methane", "naphtha"]

    ip = _industry_production.loc[region, subcategories]  # kt/a
    sr = _sector_ratios["DE"].loc[carrier, subcategories]  # MWh/tMaterial
    non_energy = sr.multiply(ip).sum(axis=1) * 1e3 * MWh2PJ

    # write var
    var["Final Energy|Non-Energy Use|Gases"] = non_energy.methane + CH4_for_NH3

    gas_fractions = _get_gas_fractions(n, region)
    for gas_type in gas_fractions.index:
        var[f"Final Energy|Non-Energy Use|Gases|{gas_type}"] = (
            var["Final Energy|Non-Energy Use|Gases"] * gas_fractions[gas_type]
        )

    var["Final Energy|Non-Energy Use|Liquids"] = non_energy.naphtha

    var["Final Energy|Non-Energy Use|Liquids|Petroleum"] = (
        non_energy.naphtha * oil_fossil_fraction
    )
    var["Final Energy|Non-Energy Use|Liquids|Efuel"] = non_energy.naphtha * (
        1 - oil_fossil_fraction
    )
    var["Final Energy|Non-Energy Use|Liquids|Biomass"] = 0

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

    var["Final Energy|Industry|Electricity"] = industry_demand.get("electricity")
    # or use: sum_load(n, "industry electricity", region)
    # electricity is not used for non-energy purposes
    var["Final Energy|Industry excl Non-Energy Use|Electricity"] = var[
        "Final Energy|Industry|Electricity"
    ]

    var["Final Energy|Industry|Heat"] = industry_demand.get("low-temperature heat")
    # heat is not used for non-energy purposes
    var["Final Energy|Industry excl Non-Energy Use|Heat"] = var[
        "Final Energy|Industry|Heat"
    ]

    # var["Final Energy|Industry|Solar"] = \
    # !: Included in |Heat

    # var["Final Energy|Industry|Geothermal"] = \
    # Not implemented

    var["Final Energy|Industry|Gases"] = industry_demand.get("methane")

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
            var[f"Final Energy|Industry excl Non-Energy Use|Gases"]
            * gas_fractions[gas_type]
        )

    # var["Final Energy|Industry|Power2Heat"] = \
    # Q: misleading description

    var["Final Energy|Industry|Hydrogen"] = industry_demand.get("hydrogen")
    # subtract non-energy used hydrogen from total hydrogen demand
    var["Final Energy|Industry excl Non-Energy Use|Hydrogen"] = (
        var["Final Energy|Industry|Hydrogen"]
        - var["Final Energy|Non-Energy Use|Hydrogen"]
    )

    var["Final Energy|Industry|Liquids|Petroleum"] = (
        sum_load(n, "naphtha for industry", region) * oil_fossil_fraction
    )

    # subtract non-energy used petroleum from total petroleum demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids|Petroleum"] = (
        var["Final Energy|Industry|Liquids|Petroleum"]
        - var["Final Energy|Non-Energy Use|Liquids|Petroleum"]
    )

    var["Final Energy|Industry|Liquids|Efuel"] = sum_load(
        n, "naphtha for industry", region
    ) * (1 - oil_fossil_fraction)
    # subtract non-energy used efuels from total efuels demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids|Efuel"] = (
        var["Final Energy|Industry|Liquids|Efuel"]
        - var["Final Energy|Non-Energy Use|Liquids|Efuel"]
    )

    var["Final Energy|Industry|Liquids"] = sum_load(n, "naphtha for industry", region)
    # subtract non-energy used liquids from total liquid demand
    var["Final Energy|Industry excl Non-Energy Use|Liquids"] = (
        var["Final Energy|Industry|Liquids"]
        - var["Final Energy|Non-Energy Use|Liquids"]
    )

    # TODO This is plastics not liquids for industry! Look in industry demand!

    # var["Final Energy|Industry|Other"] = \

    var["Final Energy|Industry|Solids|Biomass"] = industry_demand.get("solid biomass")
    var["Final Energy|Industry excl Non-Energy Use|Solids|Biomass"] = var[
        "Final Energy|Industry|Solids|Biomass"
    ]

    mwh_coal_per_mwh_coke = 1.366
    # Coke is added as a coal demand, so we need to convert back to units of coke for final energy
    var["Final Energy|Industry|Solids|Coal"] = (
        industry_demand.get("coal")
        + industry_demand.get("coke") / mwh_coal_per_mwh_coke
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
        .multiply(MWh2PJ)
    )

    var["Final Energy|Residential and Commercial|Electricity"] = (
        low_voltage_electricity[
            # carrier does not contain one of the following substrings
            ~low_voltage_electricity.index.str.contains(
                "urban central|industry|agriculture|charger|distribution"
                # Excluding chargers (battery and EV)
            )
        ].sum()
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
                "rural water tanks charger",
            ],
            errors="ignore",
        )
        .multiply(MWh2PJ)
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
        .multiply(MWh2PJ)
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
        .multiply(MWh2PJ)
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
        .multiply(MWh2PJ)
    )

    var["Final Energy|Residential and Commercial|Liquids"] = oil_usage.get(
        "urban decentral oil boiler", 0
    ) + oil_usage.get("rural oil boiler", 0)

    var["Final Energy|Residential and Commercial|Liquids|Petroleum"] = (
        var["Final Energy|Residential and Commercial|Liquids"] * oil_fossil_fraction
    )

    # var["Final Energy|Residential and Commercial|Other"] = \
    # var["Final Energy|Residential and Commercial|Solids|Coal"] = \
    # ! Not implemented

    biomass_usage = (
        n.statistics.withdrawal(bus_carrier="solid biomass", **kwargs)
        .filter(like=region)
        .groupby(["carrier"])
        .sum()
        .multiply(MWh2PJ)
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

    var["Final Energy|Transportation|Electricity"] = sum_load(
        n, "land transport EV", region
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

    oil_fossil_fraction = _get_oil_fossil_fraction(n, region)

    var["Final Energy|Transportation|Liquids"] = (
        sum_load(n, "land transport oil", region)
        + (
            sum_load(n, "kerosene for aviation", region)
            * (1 - international_aviation_fraction)
        )
        + (
            sum_load(n, "shipping oil", region)
            * (1 - international_navigation_fraction)
        )
    )

    var["Final Energy|Transportation|Methanol"] = sum_load(
        n, "shipping methanol", region
    ) * (1 - international_navigation_fraction)

    # var["Final Energy|Transportation|Liquids|Biomass"] = \
    # var["Final Energy|Transportation|Liquids|Synthetic Fossil"] = \
    var["Final Energy|Transportation|Liquids|Petroleum"] = (
        var["Final Energy|Transportation|Liquids"] * oil_fossil_fraction
    )

    var["Final Energy|Transportation|Liquids|Efuel"] = var[
        "Final Energy|Transportation|Liquids"
    ] * (1 - oil_fossil_fraction)

    var["Final Energy|Bunkers|Aviation"] = var[
        "Final Energy|Bunkers|Aviation|Liquids"
    ] = (sum_load(n, "kerosene for aviation", region) * international_aviation_fraction)

    # TODO Navigation hydrogen

    var["Final Energy|Bunkers|Navigation|Liquids"] = (
        sum_load(n, "shipping oil", region) * international_navigation_fraction
    )
    var["Final Energy|Bunkers|Navigation|Methanol"] = (
        sum_load(n, "shipping methanol", region) * international_navigation_fraction
    )

    var["Final Energy|Bunkers|Navigation"] = (
        var["Final Energy|Bunkers|Navigation|Liquids"]
        + var["Final Energy|Bunkers|Navigation|Methanol"]
    )

    # var["Final Energy|Bunkers|Navigation|Gases"] = \
    # ! Not implemented
    # var["Final Energy|Bunkers|Navigation|Hydrogen"] = \
    # ! Not used

    var["Final Energy|Bunkers"] = (
        var["Final Energy|Bunkers|Navigation"] + var["Final Energy|Bunkers|Aviation"]
    )

    var["Final Energy|Transportation"] = (
        var["Final Energy|Transportation|Electricity"]
        + var["Final Energy|Transportation|Liquids"]
        + var["Final Energy|Transportation|Hydrogen"]
        + var["Final Energy|Transportation|Methanol"]
    )

    var["Final Energy|Agriculture|Electricity"] = sum_load(
        n, "agriculture electricity", region
    )
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
    # !!! TODO this assert is temporarily disbaled because of https://github.com/PyPSA/pypsa-eur/issues/985

    var["Final Energy|Electricity"] = (
        var["Final Energy|Agriculture|Electricity"]
        + var["Final Energy|Residential and Commercial|Electricity"]
        + var["Final Energy|Transportation|Electricity"]
        + var["Final Energy|Industry excl Non-Energy Use|Electricity"]
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
        var["Final Energy|Liquids"] * oil_fossil_fraction
    )

    var["Final Energy|Liquids|Efuel"] = var["Final Energy|Liquids"] * (
        1 - oil_fossil_fraction
    )

    var["Final Energy|Heat"] = (
        var["Final Energy|Agriculture|Heat"]
        + var["Final Energy|Residential and Commercial|Heat"]
        + var["Final Energy|Industry excl Non-Energy Use|Heat"]
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
        .multiply(MWh2PJ)
    )

    var["Final Energy|Waste"] = waste_withdrawal.filter(like="waste CHP").sum()

    var["Final Energy incl Non-Energy Use incl Bunkers"] = (
        var["Final Energy|Industry"]
        + var["Final Energy|Residential and Commercial"]
        + var["Final Energy|Agriculture"]
        + var["Final Energy|Transportation"]
        + var["Final Energy|Bunkers"]
    )

    var["Final Energy"] = (
        var["Final Energy|Industry excl Non-Energy Use"]
        + var["Final Energy|Residential and Commercial"]
        + var["Final Energy|Agriculture"]
        + var["Final Energy|Transportation"]
    )

    # The general problem with final energy is that for most of these categories
    # feedstocks shouls be excluded (i.e., non-energy use)
    # However this is hard to do in PyPSA.
    # TODO nevertheless it would be nice to do exactly that

    return var


def get_emissions(n, region, _energy_totals):

    energy_totals = _energy_totals.loc[region[0:2]]

    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
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

    var["Emissions|CO2"] = co2_emissions.sum() - co2_atmosphere_withdrawal.sum()

    co2_storage = (
        n.statistics.supply(bus_carrier="co2 stored", **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .multiply(t2Mt)
    )

    # Assert neglible numerical errors / leakage in stored CO2
    assert co2_storage.get("co2 stored", 0) < 1.0
    co2_storage.drop("co2 stored", inplace=True, errors="ignore")

    CHP_emissions = (
        n.statistics.supply(bus_carrier="co2", **kwargs)
        .filter(like=region)
        .filter(like="CHP")
        .multiply(t2Mt)
    )

    # exclude waste CHPs because they are accounted separately
    CHP_emissions = CHP_emissions[
        ~CHP_emissions.index.get_level_values("carrier").str.contains("waste")
    ]

    ## Account for carbon neutral fuels (e-fuels, biogas, ...)

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

    oil_fossil_fraction = _get_oil_fossil_fraction(n, region)
    # Assuming that efuel emissions are generated at the production site
    var["Emissions|CO2|Energy|Production|From Liquids"] = co2_emissions.loc[
        co2_emissions.index.isin(oil_techs)
    ].sum() * (1 - oil_fossil_fraction)

    # Fossil fuel emissions are generated where they get burned

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

    gas_fractions = _get_gas_fractions(n, region)

    var["Emissions|CO2|Energy|Production|From Gases"] = co2_emissions.loc[
        co2_emissions.index.isin(gas_techs)
    ].sum() * (1 - gas_fractions["Natural Gas"])
    co2_emissions.loc[co2_emissions.index.isin(gas_techs)] *= gas_fractions[
        "Natural Gas"
    ]
    CHP_emissions.loc[
        CHP_emissions.index.get_level_values("carrier").isin(gas_techs)
    ] *= gas_fractions["Natural Gas"]

    # TODO Methanol

    CHP_E_to_H = (
        n.links.loc[CHP_emissions.index.get_level_values("name")].efficiency
        / n.links.loc[CHP_emissions.index.get_level_values("name")].efficiency2
    )

    CHP_E_fraction = CHP_E_to_H * (1 / (CHP_E_to_H + 1))

    ccs_carriers = co2_storage.index.intersection(co2_atmosphere_withdrawal.index)
    fossil_cc_carriers = co2_storage.index.difference(co2_atmosphere_withdrawal.index)

    negative_cc = co2_storage.reindex(ccs_carriers)
    fossil_cc = co2_storage.reindex(fossil_cc_carriers)

    assert isclose(fossil_cc.sum() + negative_cc.sum(), co2_storage.sum())

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

    negative_ccs = total_ccs - fossil_cc.sum()

    co2_negative_emissions = negative_cc.multiply(negative_ccs / negative_cc.sum())

    if negative_ccs < 0:
        co2_negative_emissions *= 0
        # If not enough CO2 is captured, than additional emissions occur
        fossil_cc_emissions = -fossil_cc.multiply(negative_ccs / fossil_cc.sum())
        # All captured fossil should be sequestered for e-fuels to be carbon neutral
        # If this warning appears repeatedly we will need to add a hard constraint
        print("WARNING! Not all CO2 capture from fossil sources is captured!!!")
        print("total_ccs - fossil_cc: ", total_ccs - fossil_cc.sum())
        # TODO what to do with fossil_cc_emissions???
        if (n.links.build_year.max() == 2045) and (
            total_ccs - fossil_cc.sum() > 1
        ):  # > 1 for numerical errors
            raise Exception("Not enough CCS in 2045!")

    if not co2_atmosphere_withdrawal.get("urban central solid biomass CHP CC"):
        biomass_CHP_correction_factor = 0
    else:
        biomass_CHP_correction_factor = min(
            1,  # Can not be > 1, taking minimum to avoid numerical errors
            co2_negative_emissions.get("urban central solid biomass CHP CC")
            / co2_atmosphere_withdrawal.get("urban central solid biomass CHP CC"),
        )

    negative_CHP_emissions = (
        n.statistics.withdrawal(bus_carrier="co2", **kwargs)
        .filter(like=region)
        .filter(like="solid biomass CHP CC")
        .multiply(t2Mt)
        .multiply(  # Correcting for actual negative emissions
            biomass_CHP_correction_factor
        )
    )

    negative_CHP_E_to_H = (
        n.links.loc[negative_CHP_emissions.index.get_level_values("name")].efficiency
        / n.links.loc[negative_CHP_emissions.index.get_level_values("name")].efficiency2
    )

    negative_CHP_E_fraction = negative_CHP_E_to_H * (1 / (negative_CHP_E_to_H + 1))

    var["Carbon Sequestration"] = total_ccs

    var["Carbon Sequestration|DACCS"] = var["Carbon Sequestration"] * (
        co2_storage.filter(like="DAC").sum() / co2_storage.sum()
    )

    var["Carbon Sequestration|BECCS"] = var["Carbon Sequestration"] * (
        co2_storage.filter(like="bio").sum() / co2_storage.sum()
    )

    var["Carbon Sequestration|Other"] = (
        var["Carbon Sequestration"]
        - var["Carbon Sequestration|DACCS"]
        - var["Carbon Sequestration|BECCS"]
    )

    # ! LULUCF should also be subtracted (or added??), we get from REMIND,
    # TODO how to consider it here?

    # Make sure these values are about right
    var["Emissions|CO2|Industrial Processes"] = co2_emissions.reindex(
        [
            "process emissions",
            "process emissions CC",
        ]
    ).sum() + co2_emissions.get("industry methanol", 0)
    # process emissions is mainly cement, methanol is used for chemicals

    var["Emissions|CO2|Energy|Demand|Industry"] = co2_emissions.reindex(
        ["gas for industry", "gas for industry CC", "coal for industry"]
    ).sum() - co2_negative_emissions.get(
        "solid biomass for industry CC",
        0,
    )

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

    var["Emissions|CO2|Energy|Demand|Transportation"] = (
        co2_emissions.get("land transport oil", 0)
        + (
            co2_emissions.get("kerosene for aviation")
            * (1 - international_aviation_fraction)
        )
        + (
            co2_emissions.filter(like="shipping").sum()
            * (1 - international_navigation_fraction)
        )
    )

    var["Emissions|CO2|Energy|Demand|Transportation|Domestic Aviation"] = (
        co2_emissions.get("kerosene for aviation")
        * (1 - international_aviation_fraction)
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

    var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"] = (
        co2_emissions.reindex(
            [
                "OCGT",
                "CCGT",
                "coal",
                "lignite",
                "oil",
            ],
        ).sum()
        + CHP_emissions.multiply(CHP_E_fraction).values.sum()
    )

    var["Emissions|CO2|Energy|Supply|Electricity"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"]
        - negative_CHP_emissions.multiply(negative_CHP_E_fraction).values.sum()
    )

    var["Emissions|Gross Fossil CO2|Energy|Supply|Heat"] = (
        co2_emissions.filter(like="urban central")
        .filter(like="boiler")  # in 2020 there might be central oil boilers?!
        .sum()
        + CHP_emissions.multiply(1 - CHP_E_fraction).values.sum()
    )

    var["Emissions|CO2|Energy|Supply|Heat"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Heat"]
        - negative_CHP_emissions.multiply(1 - negative_CHP_E_fraction).values.sum()
    )

    var["Emissions|CO2|Energy|Supply|Electricity and Heat"] = (
        var["Emissions|CO2|Energy|Supply|Heat"]
        + var["Emissions|CO2|Energy|Supply|Electricity"]
    )

    var["Emissions|CO2|Energy|Supply|Hydrogen"] = var[
        "Emissions|Gross Fossil CO2|Energy|Supply|Hydrogen"
    ] = co2_emissions.filter(like="SMR").sum()

    var["Emissions|CO2|Energy|Supply|Gases"] = (-1) * co2_negative_emissions.filter(
        like="biogas to gas"
    ).sum()

    var["Emissions|CO2|Supply|Non-Renewable Waste"] = co2_emissions.reindex(
        [
            "HVC to air",
            "waste CHP",
            "waste CHP CC",
        ]
    ).sum()

    var["Emissions|CO2|Energy|Supply|Liquids and Gases"] = var[
        "Emissions|CO2|Energy|Supply|Liquids"
    ] = co2_emissions.get("oil refining", 0)

    # var["Emissions|CO2|Energy|Supply|Gases"] + \

    var["Emissions|CO2|Energy|Supply"] = (
        var["Emissions|CO2|Energy|Supply|Gases"]
        + var["Emissions|CO2|Energy|Supply|Hydrogen"]
        + var["Emissions|CO2|Energy|Supply|Electricity and Heat"]
        + var["Emissions|CO2|Energy|Supply|Liquids"]
    )

    # var["Emissions|CO2|Energy|Supply|Other Sector"] = \
    # var["Emissions|CO2|Energy|Supply|Solids"] = \

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

    var["Emissions|Gross Fossil CO2|Energy|Supply"] = (
        var["Emissions|Gross Fossil CO2|Energy|Supply|Electricity"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Heat"]
        + var["Emissions|Gross Fossil CO2|Energy|Supply|Hydrogen"]
    )

    emission_difference = var["Emissions|CO2"] - (
        var["Emissions|CO2|Energy and Industrial Processes"]
        + var["Emissions|CO2|Energy|Demand|Bunkers"]
        + var["Emissions|CO2|Supply|Non-Renewable Waste"]
        - co2_negative_emissions.get("DAC", 0)
        + var["Emissions|CO2|Energy|Production|From Liquids"]
        + var["Emissions|CO2|Energy|Production|From Gases"]
        - co2_atmosphere_withdrawal.subtract(co2_negative_emissions).sum()
    )
    print(
        "Differences in accounting for CO2 emissions:",
        emission_difference,
    )

    assert emission_difference < 1e-2  # Improve numerical stability

    return var


# functions for prices
def get_nodal_flows(n, bus_carrier, region, query="index == index or index != index"):
    """
    Get the nodal flows for a given bus carrier and region.

    Parameters:
        n (pypsa.Network): The PyPSA network object.
        bus_carrier (str): The bus carrier for which to retrieve the nodal flows.
        region (str): The region for which to retrieve the nodal flows.
        query (str, optional): A query string to filter the nodal flows. Defaults to 'index == index or index != index'.

    Returns:
        pandas.DataFrame: The nodal flows for the specified bus carrier and region.
    """

    groupby = n.statistics.groupers.get_name_bus_and_carrier

    result = (
        n.statistics.withdrawal(
            bus_carrier=bus_carrier,
            groupby=groupby,
            aggregate_time=False,
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

    Parameters:
        n (pypsa.Network): The PyPSA network object.
        bus_carrier (str): The bus carrier for which to retrieve the nodal flows.
        region (str): The region for which to retrieve the nodal flows.
        query (str, optional): A query string to filter the nodal flows. Defaults to 'index == index or index != index'.

    Returns:
        pandas.DataFrame: The nodal flows for the specified bus carrier and region.
    """

    groupby = n.statistics.groupers.get_name_bus_and_carrier

    result = (
        n.statistics.supply(
            bus_carrier=bus_carrier,
            groupby=groupby,
            aggregate_time=False,
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

    Parameters:
    - n (pandas.DataFrame): The network model.
    - load_carrier (str): The load carrier to calculate the price for.
    - region (str): The region to calculate the price in.

    Returns:
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

    Parameters:
    - n (pandas.DataFrame): The network model.
    - region (str): The region to consider.
    - carrier (str): The carrier of the generators.

    Returns:
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

    Parameters:
        n (pypsa.Network): The PyPSA network object.
        region (str): The region to consider for the links.
        carrier (str): The carrier of the links.
        gen_bus (str, optional): The bus where the main generation of the link takes place. Defaults to "p1".

    Returns:
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

    Parameters:
    - n (PyPSa network): The Ariadne model scenario output.
    - region (str): The region for which the prices are calculated.

    Returns:
    - var (pandas.Series): A series containing the calculated prices.

    This function calculates the prices of different energy sources in the Ariadne model
    based on the nodal flows and marginal prices of the model. The calculated prices are
    stored in a pandas Series object and returned.
    """

    var = pd.Series()

    kwargs = {
        "groupby": n.statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }
    try: 
        co2_limit_de = n.global_constraints.loc["co2_limit-DE", "mu"]
    except KeyError:
        co2_limit_de = 0


    # co2 additions
    co2_price = (
        -n.global_constraints.loc["CO2Limit", "mu"]
        - co2_limit_de
    )
    # specific emissions in tons CO2/MWh according to n.links[n.links.carrier =="your_carrier].efficiency2.unique().item()
    specific_emisisons = {
        "oil": 0.2571,
        "gas": 0.198,  # OCGT
        "hard coal": 0.3361,
        "lignite": 0.4069,
    }

    nodal_flows_lv = get_nodal_flows(
        n,
        "low voltage",
        region,
        query="not carrier.str.contains('agriculture')"
        "& not carrier.str.contains('industry')"
        "& not carrier.str.contains('urban central')",
    )

    nodal_prices_lv = n.buses_t.marginal_price[nodal_flows_lv.columns]

    # electricity price at the final level in the residential sector. Prices should include the effect of carbon prices.
    var["Price|Final Energy|Residential and Commercial|Electricity"] = (
        nodal_flows_lv.mul(nodal_prices_lv).values.sum()
        / nodal_flows_lv.values.sum()
        / MWh2GJ
    )

    # vars: Tier 1, Category: energy(price)

    nodal_flows_bm = get_nodal_flows(n, "solid biomass", region)
    nodal_prices_bm = n.buses_t.marginal_price[nodal_flows_bm.columns]

    # primary energy consumption of purpose-grown bioenergy crops, crop and forestry residue bioenergy, municipal solid waste bioenergy, traditional biomass, including renewable waste
    var["Price|Primary Energy|Biomass"] = (
        nodal_flows_bm.mul(nodal_prices_bm).values.sum()
        / nodal_flows_bm.values.sum()
        / MWh2GJ
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
        coal_fraction * specific_emisisons["hard coal"] * co2_price
        + lignite_fraction * specific_emisisons["lignite"] * co2_price
    )

    var["Price|Primary Energy|Coal"] = (
        get_weighted_costs(
            [coal_price, lignite_price], [nf_coal.values.sum(), nf_lignite.values.sum()]
        )
        + co2_add_coal
    ) / MWh2GJ

    # Price|Primary Energy|Gas
    nodal_flows_gas = get_nodal_flows(n, "gas", region)
    nodal_prices_gas = n.buses_t.marginal_price[nodal_flows_gas.columns]

    # co2 part
    gas_fractions = _get_gas_fractions(n, region)
    co2_add_gas = gas_fractions["Natural Gas"] * specific_emisisons["gas"] * co2_price

    var["Price|Primary Energy|Gas"] = (
        nodal_flows_gas.mul(nodal_prices_gas).values.sum()
        / nodal_flows_gas.values.sum()
        + co2_add_gas
    ) / MWh2GJ

    # Price|Primary Energy|Oil
    # if oil bus is unravelled change "EU" into region
    nodal_flows_oil = get_nodal_flows(n, "oil", "EU")
    nodal_prices_oil = n.buses_t.marginal_price[nodal_flows_oil.columns]

    # co2 part
    oil_fossil_fraction = _get_oil_fossil_fraction(n, region)
    co2_add_oil = oil_fossil_fraction * specific_emisisons["oil"] * co2_price

    var["Price|Primary Energy|Oil"] = (
        nodal_flows_oil.mul(nodal_prices_oil).values.sum()
        / nodal_flows_oil.values.sum()
        + co2_add_oil
    ) / MWh2GJ

    # Price|Secondary Energy|Electricity
    # electricity price at the secondary level, i.e. for large scale consumers (e.g. aluminum production). Prices should include the effect of carbon prices.

    nodal_flows_ac = get_nodal_flows(
        n, "AC", region, query="not carrier.str.contains('gas')"
    )
    nodal_prices_ac = n.buses_t.marginal_price[nodal_flows_ac.columns]

    var["Price|Secondary Energy|Electricity"] = (
        nodal_flows_ac.mul(nodal_prices_ac).values.sum()
        / nodal_flows_ac.values.sum()
        / MWh2GJ
    )

    var["Price|Secondary Energy|Gases|Natural Gas"] = (
        costs_gen_generators(n, region, "gas")[0] / MWh2GJ
    )

    var["Price|Secondary Energy|Gases|Hydrogen"] = (
        costs_gen_links(n, region, "Sabatier")[0] / MWh2GJ
    )

    var["Price|Secondary Energy|Gases|Biomass"] = (
        get_weighted_costs_links(["biogas to gas", "biogas to gas CC"], n, region)
        / MWh2GJ
    )

    # Price|Secondary Energy|Gases|Efuel
    # Price for gaseous Efuels at the secondary level, i.e. for large scale consumers. Prices should include the effect of carbon prices.
    # what are gaseous Efuels?

    # Price|Secondary Energy|Hydrogen (carbon costs not yet included)
    nodal_flows_h2 = get_nodal_flows(n, "H2", region)
    nodal_prices_h2 = n.buses_t.marginal_price[nodal_flows_h2.columns]

    var["Price|Secondary Energy|Hydrogen"] = (
        nodal_flows_h2.mul(nodal_prices_h2).values.sum() / nodal_flows_h2.values.sum()
    ) / MWh2GJ

    # From PIK plots
    # "Price|Final Energy|Residential and Commercial|Hydrogen" = final energy consumption by the residential sector of hydrogen
    # do we have residential applications for hydrogen?

    nf_gas_residential = get_nodal_flows(
        n,
        "gas",
        region,
        query="carrier.str.contains('rural')"
        "or carrier.str.contains('urban decentral')",
    )
    nodal_prices_gas = n.buses_t.marginal_price[nf_gas_residential.columns]

    # !!! mv much higher: check carbon effect!
    var["Price|Final Energy|Residential and Commercial|Gases"] = (
        nf_gas_residential.mul(nodal_prices_gas).values.sum()
        / nf_gas_residential.values.sum()
        / MWh2GJ
        if nf_gas_residential.values.sum() > 0
        else np.nan
    )

    # "Price|Final Energy|Residential and Commercial|Gases|Natural Gas" ?
    # "Price|Final Energy|Residential and Commercial|Liquids|Biomass" x

    var["Price|Final Energy|Residential and Commercial|Liquids|Oil"] = (
        get_weighted_costs_links(
            ["rural oil boiler", "urban decentral oil boiler"], n, region
        )
        / MWh2GJ
    )

    var["Price|Final Energy|Residential and Commercial|Liquids"] = var[
        "Price|Final Energy|Residential and Commercial|Liquids|Oil"
    ]

    var["Price|Final Energy|Residential and Commercial|Solids|Biomass"] = (
        get_weighted_costs_links(
            ["rural biomass boiler", "urban decentral biomass boiler"], n, region
        )
        / MWh2GJ
    )

    var["Price|Final Energy|Residential and Commercial|Solids"] = var[
        "Price|Final Energy|Residential and Commercial|Solids|Biomass"
    ]

    # "Price|Final Energy|Industry|Electricity"✓

    var["Price|Final Energy|Industry|Gases"] = (
        get_weighted_costs_links(["gas for industry", "gas for industry CC"], n, region)
        / MWh2GJ
    )

    # "Price|Final Energy|Industry|Heat"✓

    var["Price|Final Energy|Industry|Liquids"] = (
        price_load(n, "naphtha for industry", region)[0] / MWh2GJ
    )

    # "Price|Final Energy|Industry|Hydrogen"✓

    var["Price|Final Energy|Industry|Solids"] = (
        get_weighted_costs_links(
            [
                "solid biomass for industry",
                "solid biomass for industry CC",
                "coal for industry",
            ],
            n,
            region,
        )
        / MWh2GJ
    )

    # Rest Tier 2
    # x
    # Price|Final Energy|Transportation|Liquids|Petroleum
    # Price|Final Energy|Transportation|Liquids|Petroleum|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Petroleum|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Petroleum|Carbon Price Component
    # Price|Final Energy|Transportation|Liquids|Petroleum|Other Taxes
    # 'land transport oil' ?

    # x
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

    # x
    # Price|Final Energy|Transportation|Liquids|Biomass
    # Price|Final Energy|Transportation|Liquids|Biomass|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Biomass|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Biomass|Other Taxes

    # Price|Final Energy|Transportation|Liquids|Efuel

    # TODO THIS SEEMS INCORRECT
    df = pd.DataFrame(
        {
            c: price_load(n, c, region)
            for c in ["kerosene for aviation", "shipping methanol", "shipping oil"]
        }
    )

    var["Price|Final Energy|Transportation|Liquids|Efuel"] = (
        (df.iloc[0] * df.iloc[1]).sum() / df.iloc[1].sum() / MWh2GJ
    )

    # Price|Final Energy|Transportation|Liquids|Efuel|Sales Margin
    # Price|Final Energy|Transportation|Liquids|Efuel|Transport and Distribution
    # Price|Final Energy|Transportation|Liquids|Efuel|Other Taxes

    # Price|Final Energy|Transportation|Gases|Efuel
    # Price|Final Energy|Transportation|Gases|Efuel|Sales Margin
    # Price|Final Energy|Transportation|Gases|Efuel|Transport and Distribution
    # Price|Final Energy|Transportation|Gases|Efuel|Other Taxes

    # Price|Final Energy|Transportation|Hydrogen
    # Price|Final Energy|Transportation|Hydrogen|Sales Margin
    # Price|Final Energy|Transportation|Hydrogen|Transport and Distribution
    # Price|Final Energy|Transportation|Hydrogen|Other Taxes

    # Price|Final Energy|Transportation|Electricity

    var["Price|Final Energy|Transportation|Electricity"] = price_load(
        n, "land transport EV", region
    )[0] / (MWh2GJ)

    # Price|Final Energy|Transportation|Electricity|Sales Margin
    # Price|Final Energy|Transportation|Electricity|Transport and Distribution
    # Price|Final Energy|Transportation|Electricity|Other Taxes

    # Price|Final Energy|Residential and Commercial|Liquids|Oil

    var["Price|Final Energy|Residential and Commercial|Liquids|Oil"] = (
        get_weighted_costs_links(
            ["rural oil boiler", "urban decentral oil boiler"], n, region
        )
        / MWh2GJ
    )

    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Sales Margin
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Carbon Price Component
    # Price|Final Energy|Residential and Commercial|Liquids|Oil|Other Taxes

    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas
    # cannot really be reasonably divided from non Natural Gas resources (at least no low hanging fruit :))
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
        nf_rc_heat.mul(np_rc_heat).values.sum() / nf_rc_heat.values.sum() / MWh2GJ
    )

    # Price|Final Energy|Residential and Commercial|Heat|Sales Margin
    # Price|Final Energy|Residential and Commercial|Heat|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Heat|Other Taxes

    # Price|Final Energy|Residential and Commercial|Liquids|Biomass
    # Price|Final Energy|Residential and Commercial|Liquids|Biomass|Sales Margin
    # Price|Final Energy|Residential and Commercial|Liquids|Biomass|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Liquids|Biomass|Other Taxes

    # Price|Final Energy|Residential and Commercial|Solids|Biomass

    var["Price|Final Energy|Residential and Commercial|Solids|Biomass"] = (
        get_weighted_costs_links(
            ["rural biomass boiler", "urban decentral biomass boiler"], n, region
        )
        / MWh2GJ
    )

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

    var["Price|Final Energy|Residential and Commercial|Electricity"] = var[
        "Price|Final Energy|Residential and Commercial|Electricity"
    ]

    # Price|Final Energy|Residential and Commercial|Electricity|Sales Margin x
    # Price|Final Energy|Residential and Commercial|Electricity|Transport and Distribution
    # Price|Final Energy|Residential and Commercial|Electricity|Other Taxes
    var["Price|Final Energy|Industry|Electricity"] = price_load(
        n, "industry electricity", region
    )[0] / (MWh2GJ)

    var["Price|Final Energy|Industry|Heat"] = price_load(
        n, "low-temperature heat for industry", region
    )[0] / (MWh2GJ)

    var["Price|Final Energy|Industry|Hydrogen"] = price_load(
        n, "H2 for industry", region
    )[0] / (MWh2GJ)

    var["Price|Final Energy|Industry|Solids|Coal"] = price_load(
        n, "coal for industry", region
    )[0] / (MWh2GJ)

    # Price|Final Energy|Industry|Solids|Coal|Sales Margin x
    # Price|Final Energy|Industry|Solids|Coal|Transport and Distribution
    # Price|Final Energy|Industry|Solids|Coal|Carbon Price Component
    # Price|Final Energy|Industry|Solids|Coal|Other Taxes

    # var["Price|Final Energy|Industry|Gases|Natural Gas"] ?

    # Price|Final Energy|Industry|Gases|Natural Gas|Sales Margin x
    # Price|Final Energy|Industry|Gases|Natural Gas|Transport and Distribution
    # Price|Final Energy|Industry|Gases|Natural Gas|Carbon Price Component
    # Price|Final Energy|Industry|Gases|Natural Gas|Other Taxes

    # Price|Final Energy|Industry|Heat|Sales Margin x
    # Price|Final Energy|Industry|Heat|Transport and Distribution
    # Price|Final Energy|Industry|Heat|Other Taxes

    # Price|Final Energy|Industry|Liquids|Biomass x
    # Price|Final Energy|Industry|Liquids|Biomass|Sales Margin
    # Price|Final Energy|Industry|Liquids|Biomass|Transport and Distribution
    # Price|Final Energy|Industry|Liquids|Biomass|Other Taxes

    var["Price|Final Energy|Industry|Solids|Biomass"] = price_load(
        n, "solid biomass for industry", region
    )[0] / (MWh2GJ)

    # Price|Final Energy|Industry|Solids|Biomass|Sales Margin x
    # Price|Final Energy|Industry|Solids|Biomass|Transport and Distribution
    # Price|Final Energy|Industry|Solids|Biomass|Other Taxes

    # Price|Final Energy|Industry|Gases|Biomass ?
    # Price|Final Energy|Industry|Gases|Biomass|Sales Margin
    # Price|Final Energy|Industry|Gases|Biomass|Transport and Distribution
    # Price|Final Energy|Industry|Gases|Biomass|Other Taxes

    var["Price|Final Energy|Industry|Liquids|Efuel"] = var[
        "Price|Final Energy|Industry|Liquids"
    ]

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

    # Rest Tier3
    nodal_flows_gas = get_nodal_flows(
        n,
        "gas",
        region,
        query="not carrier.str.contains('pipeline')"
        "& not carrier == 'gas'"
        "& not carrier.str.contains('rural')"
        "& not carrier.str.contains('urban decentral')",
    )
    nodal_prices_gas = n.buses_t.marginal_price[nodal_flows_gas.columns]
    # TODO renewable gas prices
    var["Price|Secondary Energy|Gases"] = (
        nodal_flows_gas.mul(nodal_prices_gas).values.sum()
        / nodal_flows_gas.values.sum()
        / MWh2GJ
    )

    nodal_flows_oil = get_nodal_flows(
        n,
        "oil",
        "EU",
        query="not carrier.str.contains('rural')"
        "& not carrier.str.contains('urban decentral')",
    )
    nodal_prices_oil = n.buses_t.marginal_price[nodal_flows_oil.columns]

    var["Price|Secondary Energy|Liquids"] = (
        nodal_flows_oil.mul(nodal_prices_oil).values.sum()
        / nodal_flows_oil.values.sum()
        / MWh2GJ
    )

    # Price|Final Energy|Transportation|Freight|Electricity x
    # Price|Final Energy|Transportation|Freight|Gases
    # Price|Final Energy|Transportation|Freight|Hydrogen

    var["Price|Final Energy|Transportation|Freight|Liquids"] = var[
        "Price|Final Energy|Transportation|Liquids|Efuel"
    ]

    # Price|Final Energy|Transportation|Freight|Solids x

    var["Price|Final Energy|Transportation|Passenger|Electricity"] = var[
        "Price|Final Energy|Transportation|Electricity"
    ]

    # Price|Final Energy|Transportation|Passenger|Gases
    # Price|Final Energy|Transportation|Passenger|Hydrogen
    var["Price|Final Energy|Transportation|Passenger|Liquids"] = var[
        "Price|Final Energy|Transportation|Liquids|Efuel"
    ]

    # Price|Final Energy|Transportation|Passenger|Solids x

    # Price|Final Energy|Residential and Commercial|Hydrogen x
    # Price|Final Energy|Residential and Commercial|Gases|Natural Gas ?
    # Price|Final Energy|Residential and Commercial|Solids|Coal x

    # Price|Final Energy|Transportation|Electricity|Carbon Price Component ?
    # Price|Final Energy|Transportation|Gases|Carbon Price Component
    # Price|Final Energy|Transportation|Hydrogen|Carbon Price Component
    # Price|Final Energy|Transportation|Liquids|Carbon Price Component

    return var


def get_discretized_value(value, disc_int, build_threshold=0.3):

    if value == 0.0:
        return value

    add = value - value % disc_int
    value = value % disc_int
    discrete = disc_int if value > build_threshold * disc_int else 0.0

    return add + discrete


def get_grid_investments(n, costs, region, dg_cost_factor=1.0, length_factor=1.0):
    # TODO gap between years should be read from config
    # TODO Discretization units should be read from config
    var = pd.Series()

    offwind = n.generators.filter(like="offwind", axis=0).filter(like="DE", axis=0)
    offwind_connection_overnight_cost = (
        (offwind.p_nom_opt - offwind.p_nom) * offwind.connection_overnight_cost
    ) * 1e-9
    offwind_connection_ac = offwind_connection_overnight_cost.filter(like="ac")
    offwind_connection_dc = offwind_connection_overnight_cost.filter(regex="dc|float")
    var["Investment|Energy Supply|Electricity|Transmission|Offwind-AC"] = (
        offwind_connection_ac.sum() / 5
    )
    var["Investment|Energy Supply|Electricity|Transmission|Offwind-DC"] = (
        offwind_connection_dc.sum() / 5
    )
    dc_links = n.links[
        (n.links.carrier == "DC")
        & (n.links.bus0 + n.links.bus1).str.contains(region)
        & ~n.links.index.str.contains("reversed")
    ]
    dc_expansion = dc_links.p_nom_opt.apply(
        lambda x: get_discretized_value(x, 2000)
    ) - n.links.loc[dc_links.index].p_nom_min.apply(
        lambda x: get_discretized_value(x, 2000)
    )

    dc_new = (dc_expansion > 0) & (n.links.loc[dc_links.index].p_nom_min > 10)

    dc_investments = (
        dc_links.length
        * length_factor
        * (
            (1 - dc_links.underwater_fraction)
            * dc_expansion
            * costs.at["HVDC overhead", "investment"]
            + dc_links.underwater_fraction
            * dc_expansion
            * costs.at["HVDC submarine", "investment"]
        )
        + dc_new * costs.at["HVDC inverter pair", "investment"]
    )

    ac_lines = n.lines[(n.lines.bus0 + n.lines.bus1).str.contains(region)]
    ac_expansion = ac_lines.s_nom_opt.apply(
        lambda x: get_discretized_value(x, 1700)
    ) - n.lines.loc[ac_lines.index].s_nom_min.apply(
        lambda x: get_discretized_value(x, 1700)
    )
    ac_investments = (
        ac_lines.length
        * length_factor
        * ac_expansion
        * costs.at["HVAC overhead", "investment"]
    )
    var["Investment|Energy Supply|Electricity|Transmission|AC"] = (
        ac_investments.sum() + offwind_connection_ac.sum()
    ) / 5
    var["Investment|Energy Supply|Electricity|Transmission|DC"] = (
        dc_investments.sum() + offwind_connection_dc.sum()
    ) / 5
    var["Investment|Energy Supply|Electricity|Transmission"] = (
        var["Investment|Energy Supply|Electricity|Transmission|AC"]
        + var["Investment|Energy Supply|Electricity|Transmission|DC"]
    )

    distribution_grid = n.links[n.links.carrier.str.contains("distribution")].filter(
        like="DE", axis=0
    )

    year = distribution_grid.build_year.max()
    year_pre = (year - 5) if year > 2020 else 2020

    dg_expansion = (
        distribution_grid.p_nom_opt.sum()
        - distribution_grid[distribution_grid.build_year <= year_pre].p_nom_opt.sum()
    )
    dg_investment = (
        dg_expansion
        * costs.at["electricity distribution grid", "investment"]
        * dg_cost_factor
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
    year = n.links.build_year.max()
    new_h2_links = h2_links[
        ((year - 5) < h2_links.build_year) & (h2_links.build_year <= year)
    ]
    h2_costs = (
        new_h2_links.length
        * new_h2_links.p_nom_opt.apply(lambda x: get_discretized_value(x, 1500))
        * costs.at["H2 pipeline", "investment"]
    )
    # TODO add retrofitted costs!!

    var["Investment|Energy Supply|Hydrogen|Transmission"] = h2_costs.sum() / 5

    if "gas pipeline" in n.links.carrier.unique():
        gas_links = n.links[
            (
                ((n.links.carrier == "gas pipeline") & (n.links.build_year > 2020))
                | (n.links.carrier == "gas pipeline new")
            )
            & ~n.links.reversed
            & (n.links.bus0 + n.links.bus1).str.contains(region)
        ]
        year = n.links.build_year.max()
        new_gas_links = gas_links[
            ((year - 5) < gas_links.build_year) & (gas_links.build_year <= year)
        ]
        gas_costs = (
            new_gas_links.length
            * new_gas_links.p_nom_opt.apply(lambda x: get_discretized_value(x, 1200))
            * costs.at["CH4 (g) pipeline", "investment"]
        )

        var["Investment|Energy Supply|Gas|Transmission"] = gas_costs.sum() / 5

    # var["Investment|Energy Supply|Electricity|Electricity Storage"] = \
    # var["Investment|Energy Supply|Hydrogen|Fossil"] = \
    # var["Investment|Energy Supply|Hydrogen|Biomass"] = \
    # var["Investment|Energy Supply|Hydrogen|Electrolysis"] =
    # var["Investment|Energy Supply|Hydrogen|Other"] = \
    # var["Investment|Energy Supply|Liquids"] = \
    # var["Investment|Energy Supply|Liquids|Oil"] = \
    # var["Investment|Energy Supply|Liquids|Coal and Gas"] = \
    # var["Investment|Energy Supply|Liquids|Biomass"] = \
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

    # add carbon component to fossil fuels if specified
    if investment_year in snakemake.params.co2_price_add_on_fossils.keys():
        co2_price_add_on = snakemake.params.co2_price_add_on_fossils[investment_year]
    else:
        co2_price_add_on = 0.0
    try: 
        co2_limit_de = n.global_constraints.loc["co2_limit-DE", "mu"]
    except KeyError:
        co2_limit_de = 0
    var["Price|Carbon"] = (
        -n.global_constraints.loc["CO2Limit", "mu"]
        - co2_limit_de
        + co2_price_add_on
    )

    var["Price|Carbon|EU-wide Regulation All Sectors"] = (
        -n.global_constraints.loc["CO2Limit", "mu"] + co2_price_add_on
    )

    # Price|Carbon|EU-wide Regulation Non-ETS

    var["Price|Carbon|National Climate Target"] = -co2_limit_de

    # Price|Carbon|National Climate Target Non-ETS

    return var


def get_trade(n, region):
    var = pd.Series()

    def get_export_import_links(n, region, carriers):
        exporting = n.links.index[
            (n.links.carrier.isin(carriers))
            & (n.links.bus0.str[:2] == region)
            & (n.links.bus1.str[:2] != region)
        ]

        importing = n.links.index[
            (n.links.carrier.isin(carriers))
            & (n.links.bus0.str[:2] != region)
            & (n.links.bus1.str[:2] == region)
        ]

        exporting_p = (
            n.links_t.p0.loc[:, exporting]
            .clip(lower=0)
            .multiply(n.snapshot_weightings.generators, axis=0)
            .values.sum()
            - n.links_t.p0.loc[:, importing]
            .clip(upper=0)
            .multiply(n.snapshot_weightings.generators, axis=0)
            .values.sum()
        )
        importing_p = (
            n.links_t.p0.loc[:, importing]
            .clip(lower=0)
            .multiply(n.snapshot_weightings.generators, axis=0)
            .values.sum()
            - n.links_t.p0.loc[:, exporting]
            .clip(upper=0)
            .multiply(n.snapshot_weightings.generators, axis=0)
            .values.sum()
        )

        return exporting_p, importing_p

    # Trade|Primary Energy|Biomass|Volume
    # Trade|Secondary Energy|Electricity|Volume
    exporting_ac = n.lines.index[
        (n.lines.carrier == "AC")
        & (n.lines.bus0.str[:2] == region)
        & (n.lines.bus1.str[:2] != region)
    ]

    importing_ac = n.lines.index[
        (n.lines.carrier == "AC")
        & (n.lines.bus0.str[:2] != region)
        & (n.lines.bus1.str[:2] == region)
    ]
    exporting_p_ac = (
        n.lines_t.p0.loc[:, exporting_ac]
        .clip(lower=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
        .values.sum()
        - n.lines_t.p0.loc[:, importing_ac]
        .clip(upper=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
        .values.sum()
    )
    importing_p_ac = (
        n.lines_t.p0.loc[:, importing_ac]
        .clip(lower=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
        .values.sum()
        - n.lines_t.p0.loc[:, exporting_ac]
        .clip(upper=0)
        .multiply(n.snapshot_weightings.generators, axis=0)
        .values.sum()
    )

    exports_dc, imports_dc = get_export_import_links(n, region, ["DC"])

    var["Trade|Secondary Energy|Electricity|Volume"] = (
        (exporting_p_ac - importing_p_ac) + (exports_dc - imports_dc)
    ) * MWh2PJ
    var["Trade|Secondary Energy|Electricity|Gross Import|Volume"] = (
        importing_p_ac + imports_dc
    ) * MWh2PJ
    # var["Trade|Secondary Energy|Electricity|Volume|Exports"] = \
    #     (exporting_p_ac + exports_dc) * MWh2PJ

    # Trade|Secondary Energy|Hydrogen|Volume
    h2_carriers = ["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"]
    exports_h2, imports_h2 = get_export_import_links(n, region, h2_carriers)
    var["Trade|Secondary Energy|Hydrogen|Volume"] = (exports_h2 - imports_h2) * MWh2PJ
    var["Trade|Secondary Energy|Hydrogen|Gross Import|Volume"] = imports_h2 * MWh2PJ
    # var["Trade|Secondary Energy|Hydrogen|Volume|Exports"] = \
    #     exports_h2 * MWh2PJ

    # Trade|Secondary Energy|Liquids|Hydrogen|Volume
    exports_oil_renew, imports_oil_renew = get_export_import_links(
        n, region, ["renewable oil"]
    )
    var["Trade|Secondary Energy|Liquids|Hydrogen|Volume"] = (
        exports_oil_renew - imports_oil_renew
    ) * MWh2PJ
    var["Trade|Secondary Energy|Liquids|Hydrogen|Gross Import|Volume"] = (
        imports_oil_renew * MWh2PJ
    )
    # var["Trade|Secondary Energy|Liquids|Hydrogen|Volume|Exports"] = \
    #     exports_oil_renew * MWh2PJ

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
        DE_bio_fraction = DE_renewable_gas.filter(like="biogas to gas").sum() / DE_renewable_gas.sum()

    if EU_renewable_gas.sum() == 0:
        EU_bio_fraction = 0
    else:
        EU_bio_fraction = EU_renewable_gas.filter(like="biogas to gas").sum() / EU_renewable_gas.sum()

    assert region == "DE" # only DE is implemented at the moment

    exports_gas_renew, imports_gas_renew = get_export_import_links(
        n, region, ["renewable gas"]
    )
    var["Trade|Secondary Energy|Gases|Hydrogen|Volume"] = (
        exports_gas_renew * (1 - DE_bio_fraction) 
        - imports_gas_renew * (1 - EU_bio_fraction)
    ) * MWh2PJ
    var["Trade|Secondary Energy|Gases|Hydrogen|Gross Import|Volume"] = (
        imports_gas_renew * (1 - EU_bio_fraction) * MWh2PJ
    )

    var["Trade|Secondary Energy|Gases|Biomass|Volume"] = (
        exports_gas_renew * DE_bio_fraction
        - imports_gas_renew * EU_bio_fraction
    ) * MWh2PJ
    var["Trade|Secondary Energy|Gases|Biomass|Gross Import|Volume"] = (
        imports_gas_renew * EU_bio_fraction * MWh2PJ
    )


    # TODO add methanol trade, renewable gas trade

    # Trade|Primary Energy|Coal|Volume
    # Trade|Primary Energy|Gas|Volume

    gas_fractions = _get_gas_fractions(n, region)

    if "gas pipeline" in n.links.carrier.unique():
        exports_gas, imports_gas = get_export_import_links(
            n, region, ["gas pipeline", "gas pipeline new"]
        )
        var["Trade|Primary Energy|Gas|Volume"] = (
            (exports_gas - imports_gas) * MWh2PJ
        ) * gas_fractions["Natural Gas"]
        var["Trade|Primary Energy|Gas|Volume|Imports"] = (
            imports_gas * MWh2PJ * gas_fractions["Natural Gas"]
        )
        var["Trade|Primary Energy|Gas|Volume|Exports"] = (
            exports_gas * MWh2PJ * gas_fractions["Natural Gas"]
        )

    # Trade|Primary Energy|Oil|Volume

    return var


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
    ).div(
        1e3
    )  # kton/a -> Mton/a

    var["Production|Non-Metallic Minerals|Cement"] = production.loc[region, "Cement"]
    var["Production|Steel"] = production.loc[
        "DE", ["Electric arc", "Integrated steelworks", "DRI + Electric arc"]
    ].sum()
    var["Production|Steel|Primary"] = (
        var["Production|Steel"] * config["industry"]["St_primary_fraction"][year]
    )
    var["Production|Steel|Secondary"] = var["Production|Steel"] * (
        1 - config["industry"]["St_primary_fraction"][year]
    )

    # optional:
    # var[""Production|Pulp and Paper"]
    # var["Production|Chemicals|Ammonia"]
    # var["Production|Chemicals|Methanol"]
    # var["Production|Pulp and Paper"]
    # var["Production|Non-Ferrous Metals"]

    return var


def get_operational_and_capital_costs(year):
    """
    ' This function reads in the cost data from the costs.csv file and brings
    it into the database format.
    """
    var = pd.Series()
    ind = planning_horizons.index(year)
    costs = prepare_costs(
        snakemake.input.costs[ind],
        snakemake.params.costs,
        nyears=1,
    )

    sector_dict = {
        "BEV charger": "Electricity",
        "CCGT": "Electricity",
        "DAC": "Gases",
        "H2 Electrolysis": "Hydrogen",
        "H2 Fuel Cell": "Electricity",
        "OCGT": "Electricity",
        "PHS": "Electricity",
        "V2G": "Electricity",
        "battery charger": "Electricity",
        "battery discharger": "Electricity",
        "coal": "Electricity",
        "Fischer-Tropsch": "Liquids",
        "gas pipeline": "Gases",
        "home battery charger": "Electricity",
        "home battery discharger": "Electricity",
        "hydro": "Electricity",
        "lignite": "Electricity",
        "methanolisation": "Liquids",
        "offwind-ac": "Electricity",
        "offwind-dc": "Electricity",
        "offwind-float": "Electricity",
        "oil": "Electricity",
        "onwind": "Electricity",
        "ror": "Electricity",
        "rural air heat pump": "Heat",
        "rural biomass boiler": "Heat",
        "rural gas boiler": "Heat",
        "rural oil boiler": "Heat",
        "rural ground heat pump": "Heat",
        "rural resistive heater": "Heat",
        "rural solar thermal": "Heat",
        "rural water tanks charger": "Heat",
        "rural water tanks discharger": "Heat",
        "Sabatier": "Gases",
        "solar": "Electricity",
        "solar rooftop": "Electricity",
        "solar-hsat": "Electricity",
        "solid biomass": "Heat",
        "urban central air heat pump": "Heat",
        "urban central coal CHP": "Heat",
        "urban central gas CHP": "Heat",
        "urban central gas CHP CC": "Heat",
        "urban central lignite CHP": "Heat",
        "urban central oil CHP": "Heat",
        "urban central water tanks charger": "Heat",
        "urban central water tanks discharger": "Heat",
        "urban central resistive heater": "Heat",
        "urban central solar thermal": "Heat",
        "urban central solid biomass CHP": "Heat",
        "urban central solid biomass CHP CC": "Heat",
        "urban decentral air heat pump": "Heat",
        "urban decentral resistive heater": "Heat",
        "urban decentral solar thermal": "Heat",
        "urban decentral biomass boiler": "Heat",
        "urban decentral gas boiler": "Heat",
        "urban decentral oil boiler": "Heat",
        "urban decentral water tanks charger": "Heat",
        "urban decentral water tanks discharger": "Heat",
        "waste CHP": "Heat",
        "waste CHP CC": "Heat",
    }

    grid_connection = [
        "offwind-ac",
        "offwind-dc",
        "offwind-float",
        "solar",
        "solar-hsat",
    ]

    for key, tech in costs_dict.items():
        if tech is None:
            continue
        sector = sector_dict[key]

        FOM = "OM Cost|Fixed" + "|" + sector + "|" + tech
        VOM = "OM Cost|Variable" + "|" + sector + "|" + tech
        capital = "Capital Cost" + "|" + sector + "|" + tech

        var[FOM] = costs.at[tech, "fixed"] / 1e3  # EUR/MW -> EUR/kW
        var[VOM] = costs.at[tech, "VOM"] / MWh2GJ  # EUR/MWh -> EUR/GJ
        var[capital] = costs.at[tech, "investment"] / 1e3  # EUR/MW -> EUR/kW

        if key in grid_connection:
            var[FOM] += costs.at["electricity grid connection", "fixed"] / 1e3
            var[capital] += costs.at["electricity grid connection", "investment"] / 1e3

    return var


def get_ariadne_var(
    n,
    industry_demand,
    energy_totals,
    sector_ratios,
    industry_production,
    costs,
    region,
    year,
):

    var = pd.concat(
        [
            get_capacities(n, region),
            # get_capacity_additions_simple(n,region),
            # get_installed_capacities(n,region),
            get_capacity_additions(n, region),
            get_investments(n, costs, region),
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
            ),
            get_prices(n, region),
            get_emissions(n, region, energy_totals),
            get_grid_investments(
                n,
                costs,
                region,
                dg_cost_factor=snakemake.params.dg_cost_factor,
                length_factor=snakemake.params.length_factor,
            ),
            get_policy(n, year),
            get_trade(n, region),
            # get_operational_and_capital_costs(year),
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
    )

    data = []
    for v in var.index:
        try:
            unit = var2unit[v]
        except KeyError:
            print("Warning: Variable '", v, "' not in Ariadne Database", sep="")
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
        import os
        import sys

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "export_ariadne_variables",
            simpl="",
            clusters=22,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Bal_v4",
        )

    storage_costs_dict = {
        "H2": "hydrogen storage underground",
        "EV battery": None,  # 0 i think
        "PHS": None,  #'PHS', accounted already as generator??
        "battery": "battery storage",
        "biogas": None,  # not a typical store, 0 i think
        "co2 sequestered": snakemake.params.co2_sequestration_cost,  # TODO how to consider the co2_sequestration_lifetime here
        "co2 stored": "CO2 storage tank",
        "gas": "gas storage",
        "home battery": "home battery storage",
        "hydro": None,  # `hydro`, , accounted already as generator??
        "oil": 0.02,
        "rural water tanks": "decentral water tank storage",
        "solid biomass": None,  # not a store, but a potential, 0 i think
        "urban central water tanks": "central water tank storage",
        "urban decentral water tanks": "decentral water tank storage",
    }

    config = snakemake.config
    planning_horizons = snakemake.params.planning_horizons
    ariadne_template = pd.read_excel(snakemake.input.template, sheet_name=None)
    var2unit = ariadne_template["variable_definitions"].set_index("Variable")["Unit"]
    industry_demands = [
        pd.read_csv(
            in_dem,
            index_col="TWh/a (MtCO2/a)",
        )
        .multiply(TWh2PJ)
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
        .multiply(TWh2PJ)
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

    nhours = int(snakemake.params.hours[:-1])
    nyears = nhours / 8760

    costs = list(
        map(
            lambda _costs: prepare_costs(
                _costs,
                snakemake.params.costs,
                nyears,
            ).multiply(
                1e-9
            ),  # in bn €
            snakemake.input.costs,
        )
    )

    networks = [pypsa.Network(n) for n in snakemake.input.networks]

    if "debug" == "debug":  # For debugging
        var = pd.Series()
        idx = 2
        n = networks[idx]
        c = costs[idx]
        _industry_demand = industry_demands[idx]
        _energy_totals = energy_totals.copy()
        region = "DE"
        cap_func = n.statistics.optimal_capacity
        cap_string = "Optimal Capacity|"
        dg_cost_factor = snakemake.params.dg_cost_factor
        kwargs = {
            "groupby": n.statistics.groupers.get_bus_and_carrier,
            "at_port": True,
            "nice_names": False,
        }

    yearly_dfs = []
    for i, year in enumerate(planning_horizons):
        print("Getting data for year {year}...".format(year=year))
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
                version=config["version"],
                scenario=snakemake.wildcards.run,
            )
        )

    df = reduce(
        lambda left, right: pd.merge(
            left, right, on=["Model", "Scenario", "Region", "Variable", "Unit"]
        ),
        yearly_dfs,
    )

    print("Assigning mean investments of year and year + 5 to year.")
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

    print(
        "Dropping variables which are not in the template:",
        *df.loc[df["Unit"] == "NA"]["Variable"],
        sep="\n",
    )
    df.drop(df.loc[df["Unit"] == "NA"].index, inplace=True)

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
        df.round(5).to_excel(writer, sheet_name="data", index=False)
        meta.to_frame().T.to_excel(writer, sheet_name="meta", index=False)
