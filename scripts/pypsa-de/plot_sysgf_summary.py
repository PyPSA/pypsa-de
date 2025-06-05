#!/usr/bin/env python3
# SPDX-FileCopyrightText: : 2024 PyPSA-DE authors
#
# SPDX-License-Identifier: MIT

"""
Generate summary plots and statistics for system analysis.

This script processes PyPSA networks and generates plots and metrics
for energy system analysis and comparison across scenarios.
"""

import logging
import os
from pathlib import Path
import sys
import os

sys.path.append(os.getcwd())

import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import numpy as np
import pandas as pd
import pypsa
import seaborn as sns
import yaml

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def calc_ptes_cycles(n, mean=True):
    """Calculate the number of cycles for PTES (water pits) systems."""
    pits_de = n.stores.filter(regex="DE0.*water pits", axis=0).query("e_nom_opt > 0")
    if pits_de.empty:
        return 0

    discharge = (
        n.links_t.p0.filter(regex="DE.*pits discharger")
        .mul(n.snapshot_weightings.generators, axis=0)
        .sum()
    )
    discharge.index = discharge.index.str.replace(" discharger", "")
    no_cycles = discharge.div(pits_de.e_nom_opt)
    if mean:
        return no_cycles.mean()
    else:
        return no_cycles


def calc_average_dh_price(n):
    """Calculate average district heating price in EUR/MWh."""
    dh_loads = n.loads_t.p.filter(
        regex="DE.*(urban central heat|low-temperature heat for industry)"
    )
    if dh_loads.empty:
        return 0

    # Rename columns to standardize names
    dh_loads.columns = dh_loads.columns.str.replace(
        "low-temperature heat for industry", "urban central heat"
    )
    # Aggregate columns with same name
    dh_loads = dh_loads.groupby(axis=1, level=0).sum()

    prices = n.buses_t.marginal_price.filter(regex="DE0.*urban central heat")
    if prices.empty or dh_loads.sum().sum() == 0:
        return 0

    average_dh_price = dh_loads.mul(prices).sum().sum() / dh_loads.sum().sum()
    return average_dh_price


def calc_average_elec_price(n):
    """Calculate average electricity price in EUR/MWh."""
    elec_mps = n.buses_t.marginal_price.filter(regex="DE0 \\d+$")
    elec_demand = n.loads_t.p.filter(regex="DE\\d.*(\\d+|electricity|EV)$")

    if elec_mps.empty or elec_demand.empty:
        return 0

    elec_demand.columns = elec_demand.columns.str.split(" ").str[:2].str.join(" ")
    elec_demand = elec_demand.groupby(elec_demand.columns, axis=1).sum()
    elec_costs = (elec_mps * elec_demand).sum().sum()

    if elec_demand.sum().sum() == 0:
        return 0

    return elec_costs / elec_demand.sum().sum()


def calc_average_dh_price_t_ordered(n):
    """Calculate time-ordered average district heating price."""
    loads = n.loads_t.p.filter(
        regex="DE\\d.*(urban central|low-temperature) heat"
    ).clip(lower=0)
    prices = n.buses_t.marginal_price.filter(regex="DE\\d.*urban central heat")

    if loads.empty or prices.empty or loads.sum(axis=1).isnull().any():
        return pd.Series()

    weighted_average_price_t = loads.mul(prices).sum(axis=1).div(loads.sum(axis=1))
    return weighted_average_price_t


def calc_average_electricity_price_t_ordered(n):
    """Calculate time-ordered average electricity price."""
    loads = n.buses_t.p.filter(regex="DE\\d \\d$").clip(upper=0).mul(-1)
    prices = n.buses_t.marginal_price.filter(regex="DE\\d \\d$")

    if loads.empty or prices.empty or loads.sum(axis=1).isnull().any():
        return pd.Series()

    weighted_average_price_t = loads.mul(prices).sum(axis=1).div(loads.sum(axis=1))
    return weighted_average_price_t


def calc_curtailment_de(n):
    """Calculate curtailment in TWh."""
    try:
        curtailment = (
            n.statistics.curtailment(
                groupby=n.statistics.groupers.get_bus_and_carrier, nice_names=False
            )
            .xs("Generator", level=0)
            .filter(regex="DE.*wind|solar")
            .div(1e6)
            .sum()
        )
        return curtailment
    except:
        return 0


def calc_heat_venting_de(n):
    """Calculate heat venting in TWh."""
    try:
        heat_venting = (
            n.snapshot_weightings.generators
            @ n.generators_t.p.filter(regex="DE0.*heat vent")
        ).sum() / 1e6
        return heat_venting
    except:
        return 0


def get_component_mask(lines_or_links, country, other_countries, bus=0):
    """Create mask for components connecting a country to others."""
    if bus == 0:
        # For links, check both buses
        return lines_or_links.bus0.str.contains(
            country
        ) & lines_or_links.bus1.str.contains("|".join(other_countries))
    elif bus == 1:
        # For links, check the second bus
        return lines_or_links.bus1.str.contains(
            country
        ) & lines_or_links.bus0.str.contains("|".join(other_countries))
    else:
        return (
            lines_or_links.bus0.str.contains(country)
            & lines_or_links.bus1.str.contains("|".join(other_countries))
        ) | (
            lines_or_links.bus0.str.contains("|".join(other_countries))
            & lines_or_links.bus0.str.contains(country)
        )


def calc_ic_capex_correction(n, country):
    """Calculate correction value for interconnection CAPEX assuming equal shares of connected countries."""
    other_countries = n.buses.country.unique()
    other_countries = other_countries[
        (other_countries != "") & (other_countries != country)
    ]

    ic_links_bus0_mask = get_component_mask(n.links, country, other_countries, bus=0)
    ic_links_bus0 = n.links.loc[ic_links_bus0_mask]
    capex_ic_links_bus0 = ic_links_bus0.p_nom_opt.sub(ic_links_bus0.p_nom).mul(
        ic_links_bus0.capital_cost
    )

    ic_links_bus1_mask = get_component_mask(n.links, country, other_countries, bus=1)
    ic_links_bus1 = n.links.loc[ic_links_bus1_mask]
    capex_ic_links_bus1 = ic_links_bus1.p_nom_opt.sub(ic_links_bus1.p_nom).mul(
        ic_links_bus1.capital_cost
    )

    ic_links = pd.concat([capex_ic_links_bus0, capex_ic_links_bus1], axis=0)
    ic_links["capex"] = 0.5 * (
        capex_ic_links_bus1.reindex(ic_links.index, fill_value=0)
        - capex_ic_links_bus0.reindex(ic_links.index, fill_value=0)
    )

    ic_lines_bus0_mask = get_component_mask(n.lines, country, other_countries, bus=0)
    ic_lines_bus0 = n.lines.loc[ic_lines_bus0_mask]
    capex_ic_links_bus0 = ic_lines_bus0.s_nom_opt.sub(ic_lines_bus0.s_nom).mul(
        ic_lines_bus0.capital_cost
    )

    ic_lines_bus1_mask = get_component_mask(n.lines, country, other_countries, bus=1)
    ic_lines_bus1 = n.lines.loc[ic_lines_bus1_mask]
    capex_ic_links_bus1 = ic_lines_bus1.s_nom_opt.sub(ic_lines_bus1.s_nom).mul(
        ic_lines_bus1.capital_cost
    )

    ic_lines = pd.concat([capex_ic_links_bus0, capex_ic_links_bus1], axis=0)
    ic_lines["capex"] = 0.5 * (
        capex_ic_links_bus1.reindex(ic_lines.index, fill_value=0)
        - capex_ic_links_bus0.reindex(ic_lines.index, fill_value=0)
    )

    return pd.Series(
        {
            "AC": ic_lines.capex.sum(),
            "DC": ic_links.capex.sum(),
        },
        name="capex",
    )


def calc_system_costs_country(n, country):
    """Calculate system costs for a country."""
    s = n.statistics
    capex_country = s.expanded_capex(
        groupby=n.statistics.groupers.get_bus_and_carrier
    ).filter(regex=country, axis=0)
    opex_country = s.opex(groupby=n.statistics.groupers.get_bus_and_carrier).filter(
        regex=country, axis=0
    )

    system_costs_country = capex_country.sum().sum() + opex_country.sum().sum()
    ic_costs_correction = calc_ic_capex_correction(n, country).sum()

    return system_costs_country + ic_costs_correction


def calc_h2_store_capacity(n):
    """Calculate hydrogen storage capacity in TWh."""
    h2_stores = n.stores.filter(regex="DE0.*H2 Store", axis=0)
    return h2_stores.e_nom_opt.sum() / 1e6  # TWh


def calc_costs_per_tech_country(n, country):
    """Calculate costs per technology for a specific country."""
    capex = n.statistics.expanded_capex(
        groupby=n.statistics.groupers.get_bus_and_carrier, nice_names=False
    ).filter(regex=country, axis=0)

    opex = n.statistics.opex(
        groupby=n.statistics.groupers.get_bus_and_carrier, nice_names=False
    ).filter(regex=country, axis=0)

    union_index = capex.index.union(opex.index)
    capex = capex.reindex(union_index).fillna(0)
    opex = opex.reindex(union_index).fillna(0)

    # Add Index level for the costs
    total = pd.concat([capex, opex], axis=1, keys=["capex", "opex"])
    total = total.droplevel(["component", "bus"])

    # Sum up the costs for each technology
    total = total.groupby(total.index).sum()

    ic_capex_correction = calc_ic_capex_correction(n, country)
    total.loc[ic_capex_correction.index, "capex"] + ic_capex_correction

    return total


def extract_part(scenario):
    """Extract parameter part from scenario name."""
    for prefix in ["Low", "High", "0.5", "2"]:
        if scenario.startswith(prefix):
            return scenario[len(prefix) :]
    return scenario


def get_all_colors(override_colors=None):
    """Get color mapping for plots."""
    # Default colors for different technologies
    default_colors = {
        # Base colors
        "wind power": "skyblue",
        "onshore wind": "lightblue",
        "offshore wind": "deepskyblue",
        "PV": "yellow",
        "solar": "yellow",
        "battery": "grey",
        "power grid": "navy",
        "AC": "navy",
        "DC": "navy",
        # Heat technologies
        "urban central water pits": "blue",
        "urban central water tanks": "cornflowerblue",
        "urban central water pits charge": "blue",
        "urban central water pits discharge": "blue",
        "urban central water tanks charge": "cornflowerblue",
        "urban central water tanks discharge": "cornflowerblue",
        "urban central air heat pump": "gold",
        "urban central geothermal heat": "khaki",
        "urban central geothermal heat pump": "khaki",
        "urban central geothermal heat direct utilisation": "olive",
        "urban central heat vent": "red",
        "urban central resistive heater": "mediumaquamarine",
        # Decentral heat
        "decentral heat pump": "#556B2F",
        "decentral ground heat pump": "#556B2F",
        "decentral air heat pump": "#556B2F",
        "decentral biomass boiler": "#6B8E23",
        "decentral resistive heater": "mediumaquamarine",
        # CHP and conventional
        "CHP": "brown",
        "urban central solid biomass CHP": "forestgreen",
        "urban central solid biomass CHP CC": "palegreen",
        "urban central lignite CHP": "grey",
        "urban central coal CHP": "dimgrey",
        "urban central oil CHP": "darkolivegreen",
        "urban central H2 CHP": "purple",
        # Hydrogen
        "H2 OCGT": "#8A2BE2",
        "H2 Store": "purple",
        "H2 Electrolysis": "darkviolet",
        # Others
        "low-temperature heat for industry": "darkgrey",
        "waste CHP": "plum",
        "oil": "black",
        "Fischer-Tropsch": "slateblue",
        "net electricity trade": "turquoise",
        "net hydrogen trade": "pink",
        "other technologies": "dimgrey",
        "neighbour countries": "#D3D3D3",
    }

    # Override with any user-defined colors
    if override_colors:
        default_colors.update(override_colors)

    return default_colors


def clean_and_aggregate_columns(costs_agg):
    """Clean and aggregate technology columns for better visualization."""
    # Aggregate grid technologies
    costs_agg["power grid"] = costs_agg[
        ["AC", "DC", "electricity distribution grid"]
    ].sum(axis=1)
    costs_agg.drop(["AC", "DC", "electricity distribution grid"], axis=1, inplace=True)

    # Replace residential, services, rural, and urban decentral from the column names
    costs_agg.columns = costs_agg.columns.str.replace("residential ", "")
    costs_agg.columns = costs_agg.columns.str.replace("services ", "")
    costs_agg.columns = costs_agg.columns.str.replace("rural ", "decentral ")
    costs_agg.columns = costs_agg.columns.str.replace("urban decentral ", "decentral ")
    costs_agg = costs_agg.groupby(costs_agg.columns, axis=1).sum()

    # Aggregate heat pumps
    sum_before = costs_agg.sum().sum()
    costs_agg["decentral heat pump"] = costs_agg.filter(
        regex="decentral (ground|air) heat pump"
    ).sum(axis=1)
    costs_agg.drop(
        costs_agg.filter(regex="decentral (ground|air) heat pump").columns,
        axis=1,
        inplace=True,
    )
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate H2 pipeline technologies
    sum_before = costs_agg.sum().sum()
    costs_agg["H2 pipeline"] = costs_agg.filter(regex="H2 pipeline").sum(axis=1)
    costs_agg.drop(
        costs_agg.filter(regex="H2 pipeline").columns.drop(
            "H2 pipeline", errors="ignore"
        ),
        axis=1,
        inplace=True,
    )
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate wind power
    sum_before = costs_agg.sum().sum()
    costs_agg["wind power"] = costs_agg.filter(regex="(on|off)wind").sum(axis=1)
    costs_agg.drop(costs_agg.filter(regex="(on|off)wind").columns, axis=1, inplace=True)
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate oil technologies (but not CHP)
    sum_before = costs_agg.sum().sum()
    costs_agg["oil"] = costs_agg.filter(regex="(^| )oil(?!.*CHP)").sum(axis=1)
    costs_agg.drop(
        costs_agg.filter(regex="(^| )oil(?!.*CHP)").columns.drop(
            "oil", errors="ignore"
        ),
        axis=1,
        inplace=True,
    )
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate solar technologies
    sum_before = costs_agg.sum().sum()
    costs_agg["PV"] = costs_agg.filter(regex="solar").sum(axis=1)
    costs_agg.drop(costs_agg.filter(regex="solar").columns, axis=1, inplace=True)
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate biomass technologies
    sum_before = costs_agg.sum().sum()
    costs_agg["solid biomass"] = costs_agg.filter(regex="solid biomass(?!.*CHP)").sum(
        axis=1
    )
    cols_to_drop = costs_agg.filter(regex="solid biomass(?!.*CHP)").columns
    cols_to_drop = cols_to_drop[cols_to_drop != "solid biomass"]
    costs_agg.drop(cols_to_drop, axis=1, inplace=True)
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate biogas
    sum_before = costs_agg.sum().sum()
    costs_agg["biogas"] = costs_agg.filter(regex="biogas").sum(axis=1)
    costs_agg.drop(
        costs_agg.filter(regex="biogas").columns.drop("biogas", errors="ignore"),
        axis=1,
        inplace=True,
    )
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate battery technologies
    sum_before = costs_agg.sum().sum()
    costs_agg["battery"] = costs_agg.filter(regex="battery").sum(axis=1)
    costs_agg.drop(
        costs_agg.filter(regex="battery").columns.drop("battery", errors="ignore"),
        axis=1,
        inplace=True,
    )
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    # Aggregate CHP technologies
    sum_before = costs_agg.sum().sum()
    costs_agg["CHP"] = costs_agg.filter(like="CHP").sum(axis=1)
    costs_agg.drop(
        costs_agg.filter(like="CHP").columns.drop("CHP", errors="ignore"),
        axis=1,
        inplace=True,
    )
    assert abs(sum_before - costs_agg.sum().sum()) < 1

    return costs_agg


def create_summary_df(networks):
    """Create summary dataframe with key metrics from networks."""
    df = pd.DataFrame(
        columns=[
            "scenario",
            "year",
            "total_system_costs_bnEUR",
            "total_system_costs_DE_bnEUR",
            "PTES_capacity_TWh",
            "PTES_capacity_GW",
            "PTES_no_cycles",
            "TTES_capacity_TWh",
            "TTES_capacity_GW",
            "H2_store_TWh",
            "dh_price_EUR_per_MWh",
            "electricity_price_EUR_per_MWh",
            "peak_electricity_price_EUR_per_MWh",
            "peak_dh_price_EUR_per_MWh",
            "curtailment_TWh",
            "heat_venting_TWh",
        ]
    )

    for scenario, years in networks.items():
        networks_scenario = networks[scenario]
        for year, n in networks_scenario.items():
            # Extract metrics from network
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "scenario": scenario,
                            "year": year,
                            "total_system_costs_bnEUR": n.objective / 1e9,
                            "total_system_costs_DE_bnEUR": calc_system_costs_country(
                                n, "DE"
                            )
                            / 1e9,
                            "PTES_capacity_TWh": n.stores.filter(
                                regex="DE.*urban central water pits", axis=0
                            )
                            .e_nom_opt.div(1e6)
                            .sum(),
                            "PTES_capacity_GW": n.links.filter(
                                regex="DE.*urban central water pits discharger", axis=0
                            )
                            .p_nom_opt.div(1e3)
                            .sum(),
                            "PTES_no_cycles": calc_ptes_cycles(n),
                            "TTES_capacity_TWh": n.stores.filter(
                                regex="DE.*urban central water tanks", axis=0
                            )
                            .e_nom_opt.div(1e6)
                            .sum(),
                            "TTES_capacity_GW": n.links.filter(
                                regex="DE.*urban central water tanks discharger", axis=0
                            )
                            .p_nom_opt.div(1e3)
                            .sum(),
                            "H2_store_TWh": calc_h2_store_capacity(n),
                            "dh_price_EUR_per_MWh": calc_average_dh_price(n),
                            "electricity_price_EUR_per_MWh": calc_average_elec_price(n),
                            "peak_electricity_price_EUR_per_MWh": n.buses_t.marginal_price.filter(
                                regex="DE0 \\d+$"
                            )
                            .max()
                            .max(),
                            "peak_dh_price_EUR_per_MWh": n.buses_t.marginal_price.filter(
                                regex="DE0.*urban central heat"
                            )
                            .max()
                            .max(),
                            "curtailment_TWh": calc_curtailment_de(n),
                            "heat_venting_TWh": calc_heat_venting_de(n),
                        },
                        index=[0],
                    ),
                ],
                ignore_index=True,
            )

    # Add parameter extraction for sensitivity analysis
    df["parameter"] = df.scenario.apply(extract_part)
    df["case"] = df.apply(lambda x: x.scenario.replace(x.parameter, ""), axis=1)

    return df


def create_cost_aggregation(networks):
    """Create cost aggregation dataframe by technology."""
    costs_agg = pd.DataFrame()

    for scenario, years in networks.items():
        networks_scenario = networks[scenario]
        for year, n in networks_scenario.items():
            costs_de = calc_costs_per_tech_country(n, "DE").sum(axis=1)
            costs_de["neighbour countries"] = n.objective - costs_de.sum()
            costs_de["scenario"] = scenario
            costs_de = costs_de.to_frame().T.set_index("scenario")
            costs_agg = pd.concat([costs_agg, costs_de])

    # Add parameter extraction for sensitivity analysis
    costs_agg["parameter"] = costs_agg.index.to_series().apply(extract_part)
    costs_agg["case"] = costs_agg.apply(
        lambda x: x.name.replace(x.parameter, ""), axis=1
    )
    costs_agg.set_index(["parameter", "case"], inplace=True)

    # Clean and aggregate technology columns
    costs_agg = clean_and_aggregate_columns(costs_agg)

    costs_agg.index = costs_agg.index.get_level_values(
        1
    ) + costs_agg.index.get_level_values(0)

    return costs_agg


def create_price_duration_curves(networks):
    """Create price duration curves for electricity and district heating."""
    price_curves = {}

    for scenario, n in networks.items():
        networks_scenario = networks[scenario]
        for year, n in networks_scenario.items():
            elec_prices = calc_average_electricity_price_t_ordered(n)
            if not elec_prices.empty:
                elec_price_duration = (
                    elec_prices.loc[
                        np.repeat(elec_prices.index, n.snapshot_weightings.generators)
                    ]
                    .sort_values(ascending=False)
                    .reset_index(drop=True)
                )
                elec_price_duration.index = elec_price_duration.index / (
                    len(elec_price_duration) - 1
                )
                price_curves[f"{scenario}_{year}_electricity"] = elec_price_duration

            # District heating price duration curve
            dh_prices = calc_average_dh_price_t_ordered(n)
            if not dh_prices.empty:
                dh_price_duration = (
                    dh_prices.loc[
                        np.repeat(dh_prices.index, n.snapshot_weightings.generators)
                    ]
                    .sort_values(ascending=False)
                    .reset_index(drop=True)
                )
                dh_price_duration.index = dh_price_duration.index / (
                    len(dh_price_duration) - 1
                )
                price_curves[f"{scenario}_{year}_dh"] = dh_price_duration

    return price_curves


def process_networks(run_name, scenarios, planning_horizons):
    """Process networks and extract data for analysis."""
    networks = {}
    networks_path = os.path.join("results", run_name)
    networks_path = "/home/cpschau/Code/dev/pypsa-de/results/20250406_dhsubnodes_eem"

    logger.info(f"Processing networks for run: {run_name}")
    logger.info(f"Scenarios: {scenarios}")
    logger.info(f"Planning horizons: {planning_horizons}")

    for scenario in scenarios:
        scenario_path = os.path.join(networks_path, scenario, "networks")
        if not os.path.exists(scenario_path):
            logger.warning(f"Scenario path {scenario_path} does not exist, skipping")
            continue

        for year in planning_horizons:
            # Find network file for this year
            network_file = None
            for file in os.listdir(scenario_path):
                if file.endswith(f"{year}.nc"):
                    network_file = os.path.join(scenario_path, file)
                    break

            if not network_file:
                logger.warning(
                    f"No network file found for scenario {scenario} and year {year}"
                )
                continue

            # Load network and calculate metrics
            try:
                logger.info(f"Loading network for scenario {scenario} and year {year}")
                n = pypsa.Network(network_file)
                networks[scenario] = {}
                networks[scenario][year] = n
                logger.info(
                    f"Successfully loaded network for scenario {scenario} and year {year}"
                )
            except Exception as e:
                logger.error(
                    f"Error processing network for scenario {scenario} and year {year}: {e}"
                )

    if not networks:
        logger.error("No networks could be loaded")
        return None, None, None

    # Create dataframes
    summary_df = create_summary_df(networks)
    costs_agg = create_cost_aggregation(networks)
    price_curves = create_price_duration_curves(networks)

    return networks, summary_df, costs_agg, price_curves


def plot_dual_comparison(costs_agg, scenario_A, scenario_B, colors, output_path):
    """Plot comparison between two scenarios (typically with/without PTES)."""
    logger.info(f"Plotting dual comparison between {scenario_A} and {scenario_B}")

    if scenario_A not in costs_agg.index.get_level_values(
        0
    ) or scenario_B not in costs_agg.index.get_level_values(0):
        logger.error(f"Scenarios {scenario_A} or {scenario_B} not found in data")
        return

    fig, ax = plt.subplots(ncols=2, figsize=(12, 6))

    # Extract Baseline scenario data
    costs_baseline = (
        costs_agg.loc[(scenario_A, slice(None))].squeeze().sort_values(ascending=False)
    )
    costs_baseline = costs_baseline[(costs_baseline != 0)]  # Drop zero-cost entries

    # Extract differences between scenarios
    df_diff = (
        costs_agg.loc[(scenario_B, slice(None))]
        .squeeze()
        .sub(costs_baseline.reindex(costs_agg.columns).fillna(0))
    )
    df_diff = df_diff.loc[(df_diff != 0)]  # Drop zero-difference entries

    # Group small contributions into "other technologies"
    other_indices = df_diff[df_diff.abs() < 0.02 * df_diff.abs().sum()].index
    df_diff["other technologies"] = df_diff[other_indices].sum()
    df_diff.drop(other_indices, inplace=True)

    # Handle other entries in baseline scenario
    other_indices_baseline = costs_baseline.index.intersection(other_indices).union(
        costs_baseline.index.difference(df_diff.index)
    )
    costs_baseline["other technologies"] = costs_baseline[other_indices_baseline].sum()
    costs_baseline.drop(other_indices_baseline, inplace=True)

    # Order by magnitude
    costs_baseline = costs_baseline.loc[
        costs_baseline.drop("neighbour countries")
        .abs()
        .sort_values(ascending=False)
        .index.append(pd.Index(["neighbour countries"]))
    ]

    df_diff = df_diff.loc[
        df_diff.drop("neighbour countries")
        .abs()
        .sort_values(ascending=False)
        .index.append(pd.Index(["neighbour countries"]))
    ]

    # Plot Baseline scenario
    costs_baseline.to_frame().T.div(1e9).plot.bar(
        stacked=True,
        ax=ax[0],
        color=costs_baseline.index.map(colors).fillna("black"),
        legend=False,
    )
    ax[0].set_ylabel("Costs per Technology [bn€]", fontsize=12)
    ax[0].set_xlabel("")
    ax[0].set_title(f"Total System Costs in {scenario_A} Scenario")
    ax[0].set_ylim(0, 140)

    # Plot differences
    df_diff.to_frame().T.div(1e9).plot.bar(
        stacked=True,
        ax=ax[1],
        color=df_diff.index.map(colors).fillna("black"),
        legend=False,
    )

    # Add horizontal line at 0
    ax[1].hlines(
        y=0,
        xmin=ax[1].get_xlim()[0],
        xmax=ax[1].get_xlim()[1],
        color="black",
        linewidth=0.5,
        zorder=3,
    )

    # Add markers for the total and German system cost changes
    markers_total = df_diff.sum() / 1e9
    ax[1].hlines(
        y=markers_total,
        xmin=-0.42,
        xmax=0.42,
        color="black",
        linewidth=3,
        zorder=2,
        label="net difference",
        path_effects=[patheffects.withStroke(linewidth=3)],
    )

    markers_de = df_diff.drop("neighbour countries").sum() / 1e9
    ax[1].hlines(
        y=markers_de,
        xmin=-0.42,
        xmax=0.42,
        color="black",
        linewidth=1.5,
        linestyle="--",
        zorder=2,
        label="German system difference",
        path_effects=[patheffects.withStroke(linewidth=3)],
    )

    # Annotate differences
    ref_sum = costs_baseline.sum() / 1e9
    ax[1].annotate(
        (
            f"+{markers_total:.1f}\n(+{markers_total/ref_sum*100:.2f}%)"
            if markers_total > 0
            else f"{markers_total:.1f}\n({markers_total/ref_sum*100:.2f}%)"
        ),
        (0, -0.8),
        textcoords="offset points",
        xytext=(0, -70),
        ha="center",
        va="bottom",
        path_effects=[patheffects.withStroke(linewidth=3, foreground="white")],
        bbox=dict(boxstyle="round,pad=0", edgecolor="none", facecolor="white", alpha=0),
        fontsize=14,
    )

    ref_sum_de = costs_baseline.drop("neighbour countries").sum() / 1e9
    ax[1].annotate(
        (
            f"+{markers_de:.1f}\n(+{markers_de/ref_sum_de*100:.2f}%)"
            if markers_de > 0
            else f"{markers_de:.1f}\n({markers_de/ref_sum_de*100:.2f}%)"
        ),
        (0, -0.8),
        textcoords="offset points",
        xytext=(0, 100),
        ha="center",
        va="top",
        path_effects=[patheffects.withStroke(linewidth=3, foreground="white")],
        bbox=dict(boxstyle="round,pad=0", edgecolor="none", facecolor="white", alpha=0),
        fontsize=14,
    )

    # Add titles and labels
    ax[1].set_ylabel("Cost Difference [bn€]", fontsize=12)
    ax[1].set_title(f"Difference in Investments ({scenario_A} - {scenario_B})")
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_xlabel("")

    # Increase tick size
    ax[0].tick_params(labelsize=12)
    ax[1].tick_params(labelsize=12)

    # Add legend
    handles, labels = ax[1].get_legend_handles_labels()
    labels = [
        label.replace("urban central ", "district heating ")
        .replace("water pits", "PTES")
        .replace("water tanks", "TTES")
        .replace("OCGT", "open-cycle gas turbine")
        .replace("Electolysis", "electrolysis")
        .replace("Store", " storage")
        for label in labels
    ]

    from matplotlib.lines import Line2D

    second_legend_lines = [
        Line2D([0], [0], color="black", linewidth=3),
        Line2D([0], [0], color="black", linestyle="dashed", linewidth=3),
    ]
    second_legend_labels = ["Total System Costs", "German System Costs"]

    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.05, 0.6),
        loc="center left",
        ncol=1,
        frameon=False,
        title="Technology",
    )

    fig.legend(
        second_legend_lines,
        second_legend_labels,
        bbox_to_anchor=(1.05, 0.15),
        loc="center left",
        ncol=1,
        frameon=False,
        title="Net difference in",
        fontsize=12,
    )

    # Save figure
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_path, f"dual_comparison_{scenario_A}_{scenario_B}.pdf"),
        bbox_inches="tight",
        pad_inches=0.1,
    )

    logger.info(f"Dual comparison plot saved to {output_path}")


def plot_sensitivity_analysis(
    costs_agg, sensitivities, reference_scenario, colors, output_path
):
    """Plot sensitivity analysis showing changes compared to reference scenario."""
    if reference_scenario not in costs_agg.index.get_level_values(0):
        logger.error(f"Reference scenario {reference_scenario} not found in data")
        return

    logger.info(f"Plotting sensitivity analysis for {sensitivities}")

    # Initialize figure with 2 columns and 2 rows (1 subplot per parameter)
    num_rows = (len(sensitivities) + 1) // 2
    fig, axs = plt.subplots(nrows=num_rows, ncols=2, figsize=(12, 5 * num_rows))

    # Flatten the axs array for easier iteration
    axs = axs.flatten() if len(sensitivities) > 1 else [axs]

    # Iterate over sensitivities
    for i, sensitivity in enumerate(sensitivities):
        if i >= len(axs):
            break

        ax = axs[i]

        # Prepare data for both scenarios (Low and High) of the current sensitivity
        scenarios = [f"Low{sensitivity}", f"High{sensitivity}"]
        df_diff_scenarios = []

        for scenario in scenarios:
            if scenario not in costs_agg.index.get_level_values(0):
                logger.warning(
                    f"Scenario {scenario} not found for sensitivity {sensitivity}"
                )
                continue

            # Calculate difference from reference
            df_diff_scenario = (
                costs_agg.loc[(scenario, slice(None))]
                .squeeze()
                .sub(costs_agg.loc[(reference_scenario, slice(None))].squeeze())
            )

            # Group small contributors into "other technologies"
            other_indices = df_diff_scenario[
                abs(df_diff_scenario) < 0.02 * abs(df_diff_scenario).sum()
            ].index
            df_diff_scenario["other technologies"] = df_diff_scenario[
                other_indices
            ].sum()
            df_diff_scenario.drop(other_indices, inplace=True)

            # Append to the list for concatenation
            df_diff_scenarios.append(
                df_diff_scenario.to_frame(name=scenario.replace(sensitivity, "")).T
            )

        if not df_diff_scenarios:
            continue

        # Concatenate Low and High data for the current sensitivity
        df_diff_combined = pd.concat(df_diff_scenarios).div(1e9)

        # Sort by magnitude
        df_diff_combined = df_diff_combined[
            df_diff_combined.drop("neighbour countries", axis=1)
            .abs()
            .sum()
            .sort_values(ascending=False)
            .index.append(pd.Index(["neighbour countries"]))
        ]

        # Plot stacked bar chart
        df_diff_combined.plot.bar(
            stacked=True,
            ax=ax,
            color=df_diff_combined.columns.map(colors).fillna("black"),
            legend=False,
            width=0.8,
        )

        # Add horizontal line at 0
        ax.hlines(0, -0.5, len(df_diff_combined) - 0.5, color="black", linewidth=0.5)

        # Add markers and annotations for total system cost changes
        markers_total = df_diff_combined.sum(axis=1)

        for x, (idx, y) in enumerate(markers_total.items()):
            ax.hlines(
                y=y,
                xmin=x - 0.4,
                xmax=x + 0.4,
                color="black",
                linewidth=3,
                zorder=3,
                linestyles="solid",
            )

            # Add absolute value annotation
            ax.text(
                x,
                y + 0.1 if y > 0 else y - 0.1,
                f"+{y:.2f}" if y > 0 else f"{y:.2f}",
                color="black",
                ha="center",
                va="bottom" if y > 0 else "top",
                fontsize=12,
                path_effects=[patheffects.withStroke(linewidth=3, foreground="white")],
            )

            # Add relative change annotation
            reference_costs = (
                costs_agg.loc[(reference_scenario, slice(None))].squeeze().sum() / 1e9
            )
            relative_change = (y / reference_costs) * 100
            ax.text(
                x,
                y + 0.3 if y > 0 else y - 0.3,
                f"({relative_change:+.2f}%)",
                color="black",
                ha="center",
                va="bottom" if y > 0 else "top",
                fontsize=12,
                path_effects=[patheffects.withStroke(linewidth=3, foreground="white")],
            )

        # Add titles and labels
        ax.set_title(f"{sensitivity} Parameter Sensitivity")
        ax.set_ylabel("Cost Difference [bn EUR]")
        ax.set_xticklabels(df_diff_combined.index, rotation=0)

        # Set y-axis limits consistently
        y_max = max(2.5, df_diff_combined.abs().max().max() * 1.2)
        ax.set_ylim(-y_max, y_max)

    # Add legend
    handles, labels = zip(*[ax.get_legend_handles_labels() for ax in axs])
    handles = [item for sublist in handles for item in sublist]
    labels = [item for sublist in labels for item in sublist]
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=min(len(labels), 4),
        frameon=False,
    )

    # Remove unused subplots
    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])

    # Adjust layout
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(
        os.path.join(output_path, f"sensitivity_analysis.pdf"),
        bbox_inches="tight",
        pad_inches=0.1,
    )

    logger.info(f"Sensitivity analysis plot saved to {output_path}")


def plot_price_duration_curves(price_curves, output_path):
    """Plot price duration curves for electricity and district heating."""
    # Plot electricity price duration curves
    fig_elec, ax_elec = plt.subplots(figsize=(10, 6))
    for scenario, data in price_curves.items():
        if "electricity" in scenario:
            data.plot(ax=ax_elec, label=scenario.replace("_electricity", ""))

    ax_elec.set_xlabel("Duration Proportion")
    ax_elec.set_ylabel("Price [EUR/MWh]")
    ax_elec.set_title("Electricity Price Duration Curve")
    ax_elec.legend(loc="upper right")
    ax_elec.grid(True, alpha=0.3)

    fig_elec.savefig(
        os.path.join(output_path, "electricity_price_duration.pdf"),
        bbox_inches="tight",
    )

    # Plot district heating price duration curves
    fig_dh, ax_dh = plt.subplots(figsize=(10, 6))
    for scenario, data in price_curves.items():
        if "dh" in scenario:
            data.plot(ax=ax_dh, label=scenario.replace("_dh", ""))

    ax_dh.set_xlabel("Duration Proportion")
    ax_dh.set_ylabel("Price [EUR/MWh]")
    ax_dh.set_title("District Heating Price Duration Curve")
    ax_dh.legend(loc="upper right")
    ax_dh.grid(True, alpha=0.3)

    fig_dh.savefig(
        os.path.join(output_path, "dh_price_duration.pdf"),
        bbox_inches="tight",
    )


def plot_summary_metrics(summary_df, output_path):
    """Plot key summary metrics across scenarios."""
    # Select key metrics to plot
    key_metrics = [
        "total_system_costs_bnEUR",
        "PTES_capacity_TWh",
        "TTES_capacity_TWh",
        "electricity_price_EUR_per_MWh",
        "dh_price_EUR_per_MWh",
        "curtailment_TWh",
    ]

    # Create a pivot table for plotting
    plot_data = summary_df.groupby(["scenario", "year"])[key_metrics].sum()

    # Plot the data
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, metric in enumerate(key_metrics):
        ax = axes[i]
        plot_data[metric].sort_values().plot(kind="bar", ax=ax)
        ax.set_title(metric.replace("_", " ").title())
        ax.set_ylabel(metric.split("_")[-1])
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(
        os.path.join(output_path, "summary_metrics.pdf"),
        bbox_inches="tight",
    )


def generate_summary_pdf(
    summary_df, costs_agg, reference_scenario, output_path, run_name
):
    """Generate a summary PDF with key findings."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 15))

    # Title and general info
    fig.suptitle(f"System Analysis Summary - Run: {run_name}", fontsize=16)

    # Plot 1: System costs comparison
    ax = axes[0]
    costs_data = summary_df.set_index("scenario")["total_system_costs_bnEUR"]
    costs_data.sort_values().plot(kind="bar", ax=ax)
    ax.set_title("Total System Costs")
    ax.set_ylabel("Billion EUR")
    ax.grid(True, axis="y", alpha=0.3)

    # Plot 2: PTES Capacity comparison
    ax = axes[1]
    ptes_data = summary_df.set_index("scenario")["PTES_capacity_TWh"]
    ptes_data.sort_values(ascending=False).plot(kind="bar", ax=ax)
    ax.set_title("PTES Capacity")
    ax.set_ylabel("TWh")
    ax.grid(True, axis="y", alpha=0.3)

    # Plot 3: Price comparison
    ax = axes[2]
    price_data = pd.DataFrame(
        {
            "Electricity": summary_df.set_index("scenario")[
                "electricity_price_EUR_per_MWh"
            ],
            "District Heating": summary_df.set_index("scenario")[
                "dh_price_EUR_per_MWh"
            ],
        }
    )
    price_data.plot(kind="bar", ax=ax)
    ax.set_title("Average Energy Prices")
    ax.set_ylabel("EUR/MWh")
    ax.grid(True, axis="y", alpha=0.3)

    # Save figure
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()


def main(snakemake):
    """Main function to generate plots from network data."""
    # Configure logging
    configure_logging(snakemake)

    # Get parameters from snakemake
    run_name = snakemake.params.run
    scenarios = snakemake.params.scenarios
    planning_horizons = snakemake.params.planning_horizons
    reference_scenario = snakemake.params.reference_scenario
    sensitivity_runs = snakemake.params.sensitivity_runs

    # Get color overrides from config if available
    try:
        override_colors = snakemake.params.plotting["override_tech_colors"]
    except:
        override_colors = {}

    # Create output directory
    output_path = os.path.dirname(snakemake.output.sysgf_summary)
    os.makedirs(output_path, exist_ok=True)

    logger.info(f"Generating system analysis for run: {run_name}")

    # Process networks and collect data
    networks, summary_df, costs_agg, price_curves = process_networks(
        run_name, scenarios, planning_horizons
    )

    if not networks:
        logger.error("No networks could be loaded")
        return

    # Get color mapping
    colors = get_all_colors(override_colors)

    # Generate plots

    # 1. Plot price duration curves
    plot_price_duration_curves(price_curves, output_path)

    # 2. Plot summary metrics
    plot_summary_metrics(summary_df, output_path)

    # 3. Plot dual comparison if configured
    if (
        "dual_comparison" in snakemake.params.plotting
        and snakemake.params.plotting["dual_comparison"]["enable"]
    ):
        scenario_A = snakemake.params.plotting["dual_comparison"]["scenario_A"]
        scenario_B = snakemake.params.plotting["dual_comparison"]["scenario_B"]
        if scenario_A in costs_agg.index.get_level_values(
            0
        ) and scenario_B in costs_agg.index.get_level_values(0):
            plot_dual_comparison(costs_agg, scenario_A, scenario_B, colors, output_path)

    # 4. Plot sensitivity analysis if configured
    if sensitivity_runs:
        plot_sensitivity_analysis(
            costs_agg, sensitivity_runs, reference_scenario, colors, output_path
        )

    # 5. Generate summary PDF
    generate_summary_pdf(
        summary_df,
        costs_agg,
        reference_scenario,
        snakemake.output.sysgf_summary,
        run_name,
    )

    # Save data for further analysis
    summary_df.to_csv(os.path.join(output_path, "summary_metrics.csv"))

    logger.info("System analysis completed successfully")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_sysgf_summary",
        )
    main(snakemake)
