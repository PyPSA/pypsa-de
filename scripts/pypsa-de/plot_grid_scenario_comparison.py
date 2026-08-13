# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Compare grid-topology (LT) or portfolio-evaluation (ST) network variants of
the stochastic-grid analysis (see rules/pypsa-de/stochastic_grid.smk).

Ported from a sibling repo's system_plots_scenario_comparison.py. That
script iterates `networks[scenario][year]` across cross-run scenarios (same
network, different `run`) for a full myopic pathway; here all variants
share the same `run`/Scenario name and are standalone one-year networks, so
the same plot functions are driven by `networks[network_id]` /
`variables[network_id]` instead - no `scenario`/`year` nesting needed.
`utils.py` (this directory) carries over tech_colors/sector_colors/df_to_png/
aggregate_by_keywords etc. unchanged.
"""

import logging
import os
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from pypsa.statistics import get_transmission_carriers

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utils import (
    aggregate_by_keywords,
    capa_groups,
    df_to_png,
    scenario_color,
    scenario_label,
    sector_colors,
    storage_energy_capa_groups,
    tech_colors,
)

# Regional capacity-diff maps reuse system_plots.py's map-drawing machinery
# (same directory, plain top-level import - see module docstring there).
from system_plots import (
    group_capacity_by_tech,
    plot_group_capacity_diff_maps,
    power_capacity_by_bus_carrier,
    storage_energy_capacity_by_bus_carrier,
)

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

plt.switch_backend("Agg")


def _load_get_export_import():
    """Import get_export_import from export_ariadne_variables.py.

    That module lives in the same (hyphenated) `pypsa-de` directory, so it
    can't be reached via a normal dotted import.
    """
    import importlib.util

    path = Path(__file__).parent / "export_ariadne_variables.py"
    spec = importlib.util.spec_from_file_location("export_ariadne_variables", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.get_export_import


get_export_import = _load_get_export_import()


BUS_CARRIER_TO_SECTOR = {
    "Electricity": [
        "AC",
        "low voltage",
        "battery",
        "iron-air battery",
        "EV battery",
        "home battery",
        "DSM",
    ],
    "Heat": [
        "urban central heat",
        "urban central water tanks",
        "urban central water pits",
        "rural heat",
        "rural water tanks",
        "urban decentral heat",
        "urban decentral water tanks",
    ],
    "H2": ["H2"],
    "Gas": ["gas", "gas primary", "biogas", "gas for industry", "renewable gas"],
    "Coal": ["lignite", "coal", "coal for industry"],
    "Fuels": [
        "oil",
        "oil primary",
        "land transport oil",
        "shipping oil",
        "kerosene for aviation",
        "agriculture machinery oil",
        "renewable oil",
        "naphtha for industry",
        "methanol",
        "industry methanol",
        "shipping methanol",
    ],
    "Biomass": ["solid biomass", "solid biomass for industry"],
    "CO2": ["co2", "co2 stored", "co2 sequestered", "process emissions"],
}

CARRIER_GROUPS = {
    "electricity": {
        "bus_carrier": ["AC", "low voltage"],
        "drop_stores": True,
        "drop_carriers": ["AC", "DC"],
        "threshold": 1.0,  # TWh
    },
    "heat": {
        "bus_carrier": ["urban decentral heat", "rural heat", "urban central heat"],
        "drop_stores": False,
        "drop_carriers": [],
        "threshold": 1.0,
    },
    "H2": {
        "bus_carrier": ["H2"],
        "drop_stores": False,
        "drop_carriers": [],
        "threshold": 1.0,
    },
    "oil": {
        "bus_carrier": ["oil"],
        "drop_stores": False,
        "drop_carriers": [],
        "threshold": 1.0,
    },
    "gas": {
        "bus_carrier": ["gas"],
        "drop_stores": False,
        "drop_carriers": [],
        "threshold": 1.0,
    },
    "co2": {
        "bus_carrier": ["co2"],
        "drop_stores": False,
        "drop_carriers": [],
        "threshold": 1.0,  # Mt
    },
    "solid_biomass": {
        "bus_carrier": ["solid biomass"],
        "drop_stores": False,
        "drop_carriers": [],
        "threshold": 1.0,
    },
}


# ---------------------------------------------------------------------------
# Ported ~1:1 from system_plots_scenario_comparison.py, `networks[scenario][year]`
# replaced by `networks[network_id]` (see module docstring).
# ---------------------------------------------------------------------------


def get_generation_consumption(networks, network_ids, carrier_name, config, region="DE", kwargs=None):
    if kwargs is None:
        kwargs = {
            "groupby": pypsa.statistics.groupers["name", "bus", "carrier"],
            "at_port": True,
            "nice_names": False,
        }

    gen_data = {}
    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        result = (
            n.statistics.supply(bus_carrier=config["bus_carrier"], **kwargs)
            .filter(like=region)
        )
        if config["drop_stores"]:
            result = result.drop(["Store"], errors="ignore")
        gen_data[network_id] = (
            result.groupby(["carrier"]).sum().drop(config["drop_carriers"], errors="ignore").div(1e6)
        )

    con_data = {}
    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        result = (
            n.statistics.withdrawal(bus_carrier=config["bus_carrier"], **kwargs)
            .filter(like=region)
        )
        if config["drop_stores"]:
            result = result.drop(["Store"], errors="ignore")
        con_data[network_id] = (
            result.groupby(["carrier"]).sum().drop(config["drop_carriers"], errors="ignore").div(1e6)
        )

    gen_df = pd.concat([gen_data[nid] for nid in gen_data], axis=1)
    gen_df.columns = list(gen_data.keys())
    con_df = pd.concat([con_data[nid] for nid in con_data], axis=1)
    con_df.columns = list(con_data.keys())
    return gen_df, con_df


def plot_storage_timeseries(networks, network_ids, carrier, region="DE", output_path=None):
    """Plot storage timeseries for a specific carrier."""
    fig, ax = plt.subplots(figsize=(10, 6))

    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        stores_i = n.stores[
            (n.stores.carrier == carrier) & (n.stores.bus.str.contains(region))
        ].index

        if len(stores_i) > 0:
            energy_capacity = n.stores.loc[stores_i].e_nom_opt.sum() / 1e6
            (n.stores_t.e[stores_i].sum(axis=1) / 1e6).plot(
                label=f"{scenario_label(network_id)} (Capacity: {energy_capacity:.2f} TWh)",
                ax=ax,
            )

    ax.set_xlabel("Time [h]")
    ax.set_ylabel("Energy [TWh]")
    ax.set_title(f"{carrier} Storage")
    ax.legend(fontsize=8)

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

    return fig, ax


def plot_marginal_prices(networks, network_ids, carriers, plot_dir, tag, ylim=(-50, 400)):
    """
    Plot marginal prices time series and price duration curves for specified carriers.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    fig2, ax2 = plt.subplots(1, 1, figsize=(12, 5))

    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        buses = n.buses.index[
            (n.buses.index.str[:2] == "DE") & (n.buses.carrier.isin(carriers))
        ]
        prices = n.buses_t.marginal_price[buses]

        ax.plot(prices.mean(axis=1).to_numpy(), label=scenario_label(network_id), ls=":")

        sorted_prices = np.sort(prices.mean(axis=1).values)
        pdc = np.arange(1, len(sorted_prices) + 1) / len(sorted_prices)
        ax2.plot(pdc, sorted_prices, label=scenario_label(network_id), ls=":")

    ax.set_title(f"Marginal Prices ({', '.join(carriers)})")
    ax.set_xlabel("Time")
    ax.set_ylabel("Price (EUR/MWh)")
    ax.legend(loc="upper right", frameon=True, fancybox=True, shadow=True, fontsize=8)
    ax.grid()
    ax.set_ylim(ylim)
    fig.savefig(
        f"{plot_dir}/marginal_prices_{tag}_{'_'.join(carriers)}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
    )

    ax2.set_title(f"Price Duration Curve ({', '.join(carriers)})")
    ax2.set_xlabel("PDC")
    ax2.set_ylabel("Price (EUR/MWh)")
    ax2.legend(loc="upper right", frameon=True, fancybox=True, shadow=True, fontsize=8)
    ax2.grid()
    ax2.set_ylim(ylim)
    fig2.savefig(
        f"{plot_dir}/marginal_prices_pdc_{tag}_{'_'.join(carriers)}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
    )

    plt.close(fig)
    plt.close(fig2)


def get_tsc(n, country, region="DE"):
    pypsa.options.set_option("params.statistics.drop_zero", False)
    capex = n.statistics.capex(
        groupby=pypsa.statistics.groupers["name", "carrier"], nice_names=False
    )
    opex = n.statistics.opex(
        groupby=pypsa.statistics.groupers["name", "carrier"], nice_names=False
    )

    inter_country_lines = n.lines.bus0.map(n.buses.country) != n.lines.bus1.map(n.buses.country)
    inter_country_links = n.links.bus0.map(n.buses.country) != n.links.bus1.map(n.buses.country)

    transmission_carriers = get_transmission_carriers(n).get_level_values("carrier")
    transmission_lines = n.lines.carrier.isin(transmission_carriers) & n.lines.active
    transmission_links = n.links.carrier.isin(transmission_carriers) & n.links.active

    country_transmission_lines = (
        (n.lines.bus0.str.contains(country)) & ~(n.lines.bus1.str.contains(country))
    ) | (~(n.lines.bus0.str.contains(country)) & (n.lines.bus1.str.contains(country)))
    country_transmission_links = (
        (n.links.bus0.str.contains(country)) & ~(n.links.bus1.str.contains(country))
    ) | (~(n.links.bus0.str.contains(country)) & (n.links.bus1.str.contains(country)))

    inter_country_transmission_lines = (
        inter_country_lines & transmission_lines & country_transmission_lines
    )
    inter_country_transmission_links = (
        inter_country_links & transmission_links & country_transmission_links
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

    tsc = pd.concat([capex, opex], axis=1, keys=["capex", "opex"])
    tsc = tsc.reset_index().set_index("name")
    tsc.loc[inter_country_transmission_i, ["capex", "opex"]] = (
        tsc.loc[inter_country_transmission_i, ["capex", "opex"]] / 2
    )
    tsc.rename(
        index={index: index + " " + country for index in inter_country_transmission_i},
        inplace=True,
    )

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
    tsc.rename(index={index: index + " " + region for index in to_rename_links}, inplace=True)
    tsc.rename(index={index: index + " " + region for index in to_rename_lines}, inplace=True)

    tsc = tsc.filter(like=country, axis=0)

    tsc_sum = (
        tsc.filter(like=country, axis=0)
        .drop("component", axis=1)
        .groupby("carrier")
        .sum()
    )

    return tsc_sum, tsc


def get_trade_cost(n, region, carriers):
    """
    Positive values mean cost for the domestic energy system (imports > exports)
    Negative values mean revenue for the domestic energy system (exports > imports)
    """
    export_revenue, import_cost = get_export_import(n, region, carriers, unit="EUR")
    return import_cost - export_revenue


def get_component_bus(row, n):
    """Get the bus for a component based on its type."""
    name = row.name
    component = row["component"]

    if component == "Generator":
        return n.generators.loc[name, "bus"] if name in n.generators.index else None
    elif component == "Store":
        return n.stores.loc[name, "bus"] if name in n.stores.index else None
    elif component == "StorageUnit":
        return n.storage_units.loc[name, "bus"] if name in n.storage_units.index else None
    elif component == "Link":
        return (
            n.links.loc[name, "bus1"]
            if name in n.links.index
            else n.links[n.links.carrier == row.carrier].bus1.head(1).values[0]
        )
    elif component == "Load":
        return n.loads.loc[name, "bus"] if name in n.loads.index else None
    elif component == "Line":
        return "DE0 0"
    else:
        return None


def map_bus_to_sector(bus_name, n, bus_carrier_to_sector):
    """Map a bus to its sector based on carrier."""
    if pd.isna(bus_name):
        return "Other"

    if bus_name in n.buses.index:
        bus_carrier = n.buses.loc[bus_name, "carrier"]
    else:
        return "Other"

    for sector, carriers in bus_carrier_to_sector.items():
        if bus_carrier in carriers:
            return sector

    return "Other"


def plot_price_duration_curves(
    networks,
    network_ids,
    carriers=("AC", "low voltage"),
    regions=("DE",),
    output_dir=None,
    ylim=(-50, 500),
    table_smaller=1,
    table_bigger=300,
    close_plot=True,
):
    """Plot electricity price duration curves comparison across network_ids."""
    fig, ax = plt.subplots(figsize=(9, 5))
    network_data = {}

    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        buses = n.buses[
            n.buses.carrier.isin(carriers) & n.buses.index.str.startswith(tuple(regions))
        ].index
        lmps = n.buses_t.marginal_price[buses].values.flatten()
        lmps_sorted = np.sort(lmps)[::-1]
        pct = np.arange(len(lmps_sorted)) / len(lmps_sorted) * 100
        ax.plot(
            pct,
            lmps_sorted,
            label=scenario_label(network_id),
            color=scenario_color(network_id),
            linewidth=2,
        )
        network_data[network_id] = {
            "avg": lmps_sorted.mean(),
            "std": lmps_sorted.std(),
            f"<{table_smaller}": (lmps_sorted < table_smaller).mean() * 100,
            f">{table_bigger}": (lmps_sorted > table_bigger).mean() * 100,
        }

    if not network_data:
        plt.close(fig)
        return

    y_range = ylim[1] - ylim[0]
    label_y_position = ylim[0] + 0.02 * y_range

    quantiles = [0, 10, 25, 50, 75, 100]
    for q in quantiles:
        ax.axvline(q, ls="--", lw=0.8, color="grey", alpha=0.5)
        ax.text(q + 0.5, label_y_position, f"{q}%", rotation=90, fontsize=8, alpha=0.7)

    ax.set(
        ylim=ylim,
        xlim=(-5, 105),
        xlabel="Percentage of time",
        ylabel="EUR/MWh",
    )
    ax.grid(alpha=0.3, linestyle=":")

    table_data = []
    for nid, m in network_data.items():
        table_data.append(
            [
                scenario_label(nid),
                f"{m['avg']:.2f}",
                f"{m['std']:.2f}",
                f"{m[f'<{table_smaller}']:.0f}%",
                f"{m[f'>{table_bigger}']:.0f}%",
            ]
        )

    col_labels = ["Scenario", "Avg", "Std", f"<{table_smaller}", f">{table_bigger}"]
    # First column needs to fit the (abbreviated) scenario labels.
    name_width = min(0.55, max(0.25, 0.014 * max(len(scenario_label(nid)) for nid in network_data)))
    other_width = (1 - name_width) / (len(col_labels) - 1)
    col_widths = [name_width] + [other_width] * (len(col_labels) - 1)
    table = plt.table(
        cellText=table_data,
        colLabels=col_labels,
        colWidths=col_widths,
        loc="upper right",
        cellLoc="center",
        colColours=["#34495e"] * len(col_labels),
        fontsize=7,
    )

    table.auto_set_font_size(False)
    table.set_fontsize(6.5)
    table.scale(1.0, 1.2)

    for i, nid in enumerate(network_data, start=1):
        table[i, 0].set_facecolor(scenario_color(nid))
        table[i, 0].set_text_props(weight="bold", color="white")

    for j in range(len(col_labels)):
        table[0, j].set_text_props(weight="bold", color="white", fontsize=6.5)

    plt.tight_layout()
    if output_dir:
        plt.savefig(output_dir, bbox_inches="tight", dpi=300)
    if close_plot:
        plt.close()
    else:
        plt.show()


def bar_plot_variables(
    variables,
    network_ids,
    tech_colors,
    plot_vars=None,
    sign_flip_vars=("Electricity", "Hydrogen", "eFuels"),
    title=None,
    ylabel="TWh/a",
    output_dir=None,
    filename_prefix="trade_volume",
):
    """Plot multiple trade variables across network_ids in one plot."""

    if plot_vars is None:
        plot_vars = {
            "Electricity": "Trade|Secondary Energy|Electricity|Volume",
            "Gas": "Primary Energy|Gas",
            "Oil": "Primary Energy|Oil",
            "Hydrogen": "Trade|Secondary Energy|Hydrogen|Volume",
            "eFuels": "Trade|Secondary Energy|Efuels|Volume",
            "Biomass": "Trade|Primary Energy|Biomass|Net Imports",
        }

    fig, ax = plt.subplots(figsize=(max(10, len(network_ids) * 1.6), 5))

    data = {var_name: [] for var_name in plot_vars.keys()}

    for network_id in network_ids:
        df = variables[network_id]
        for var_name, var_path in plot_vars.items():
            try:
                value = df.loc[var_path]
                if hasattr(value, "values"):
                    value = value.values[0]
                value = float(value)
                if var_name in sign_flip_vars:
                    value = value * -1
            except (KeyError, IndexError):
                value = 0.0

            data[var_name].append(value)

    x = np.arange(len(network_ids))
    n_vars = len(plot_vars)
    width = 0.8 / max(n_vars, 1)

    total_width = width * n_vars
    start_offset = -total_width / 2 + width / 2

    for i, var_name in enumerate(plot_vars.keys()):
        offset = start_offset + width * i
        color = tech_colors.get(var_name if var_name != "Hydrogen" else "H2", "gray")
        bars = ax.bar(
            x + offset,
            data[var_name],
            width,
            label=var_name,
            color=color,
            edgecolor="black",
            linewidth=1.2,
            alpha=0.85,
        )

        for bar in bars:
            height = bar.get_height()
            if abs(height) > 0.5:
                va = "bottom" if height >= 0 else "top"
                y_pos = height
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    y_pos,
                    f"{height:.1f}",
                    ha="center",
                    va=va,
                    fontsize=8,
                    fontweight="bold",
                )

    ax.set_xticks(x)
    ax.set_xticklabels([scenario_label(nid) for nid in network_ids], fontsize=9, fontweight="bold", rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(loc="best", fontsize=9, framealpha=0.9)
    ax.grid(axis="y", alpha=0.3, linestyle=":", linewidth=1)
    ax.axhline(y=0, color="black", linewidth=1.5)

    if title:
        ax.set_title(title, fontsize=13, fontweight="bold")

    plt.tight_layout()
    if output_dir:
        plt.savefig(output_dir / f"{filename_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_curtailment(networks, network_ids, tech_colors, output_dir=None, region="DE"):
    """Plot curtailment for wind and solar technologies across network_ids."""

    wind_techs = ["onwind", "offwind-ac", "offwind-dc"]
    solar_techs = ["solar", "solar rooftop", "solar-hsat"]
    kwargs = {
        "groupby": pypsa.statistics.groupers["bus", "carrier"],
        "nice_names": False,
    }

    network_data = {nid: {"wind": {}, "solar": {}} for nid in network_ids}
    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        curtailment = (
            n.statistics.curtailment(bus_carrier=["AC", "low voltage"], **kwargs)
            .filter(like=region)
            .groupby("carrier")
            .sum()
        )
        for tech in wind_techs:
            network_data[network_id]["wind"][tech] = (
                curtailment[tech] / 1e6 if tech in curtailment.index else 0.0
            )
        for tech in solar_techs:
            network_data[network_id]["solar"][tech] = (
                curtailment[tech] / 1e6 if tech in curtailment.index else 0.0
            )

    wind_total_all = sum(v for sc in network_ids for v in network_data[sc]["wind"].values())
    solar_total_all = sum(v for sc in network_ids for v in network_data[sc]["solar"].values())
    if wind_total_all < 1 and solar_total_all < 1:
        return

    fig, ax = plt.subplots(figsize=(max(8, len(network_ids) * 1.3), 5))

    x = np.arange(len(network_ids))
    width = 0.35

    bottom_wind = np.zeros(len(network_ids))
    for tech in wind_techs:
        values = [network_data[sc]["wind"][tech] for sc in network_ids]
        ax.bar(
            x - width / 2, values, width, bottom=bottom_wind,
            label=tech, color=tech_colors[tech], edgecolor="black", linewidth=0.5,
        )
        bottom_wind += values

    bottom_solar = np.zeros(len(network_ids))
    for tech in solar_techs:
        values = [network_data[sc]["solar"][tech] for sc in network_ids]
        ax.bar(
            x + width / 2, values, width, bottom=bottom_solar,
            label=tech, color=tech_colors[tech], edgecolor="black", linewidth=0.5,
        )
        bottom_solar += values

    max_value = max(max(bottom_wind), max(bottom_solar))

    if max_value < 1.0:
        plt.close(fig)
        return

    for i, network_id in enumerate(network_ids):
        wind_total = sum(network_data[network_id]["wind"].values())
        solar_total = sum(network_data[network_id]["solar"].values())

        if wind_total > 0:
            ax.text(x[i] - width / 2, wind_total + 1, f"{wind_total:.1f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
        if solar_total > 0:
            ax.text(x[i] + width / 2, solar_total + 1, f"{solar_total:.1f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([scenario_label(nid) for nid in network_ids], fontsize=9, fontweight="bold", rotation=30, ha="right")
    ax.set_ylabel("TWh/a", fontsize=12)
    ax.set_ylim(0, max_value * 1.15)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper right", fontsize=9, framealpha=0.9, ncol=2)

    ax.grid(axis="y", alpha=0.3, linestyle=":", linewidth=1)

    plt.tight_layout()

    if output_dir:
        plt.savefig(output_dir / "curtailment.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def get_capacities(networks, network_ids, years, year_of=None, max_usage=None):
    """Extract and group capacities for all network_ids and years.

    Parameters
    ----------
    year_of : dict, optional
        Maps network_id -> the single year its network actually represents.
        A (network_id, year) cell only gets real data when
        year_of[network_id] == year; every other cell is filled with
        all-zero groups. Used for the LT pathway-year rows, where most
        network_ids only have data for the target year, not the earlier
        pathway years (e.g. 2025) mixed into the same comparison - without
        this, every network_id's single network would be (wrongly) reused
        for every row. Defaults to every network_id belonging to every
        year (single shared year, the pre-pathway-year behaviour).
    max_usage : list, optional
        List of carrier names for which to use maximum usage instead of capacity
    """
    if max_usage is None:
        max_usage = []

    kwargs = {"groupby": ["bus", "carrier"], "at_port": True, "nice_names": False}
    all_data = {}
    zero_groups = {
        group_name: dict.fromkeys(techs, 0.0) for group_name, techs in capa_groups.items()
    }

    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        for year in years:
            if year_of is not None and year_of.get(network_id) != year:
                all_data[(network_id, year)] = zero_groups
                continue
            caps = (
                n.statistics.optimal_capacity(
                    bus_carrier=["AC", "low voltage"], **kwargs
                )
                .filter(like="DE")
                .groupby("carrier")
                .sum()
                .drop(["AC", "DC", "electricity distribution grid"], errors="ignore")
            )

            if max_usage:
                ct = "DE"
                buses = n.buses.index[(n.buses.index.str[:2] == ct)].drop("DE", errors="ignore")
                balance = (
                    n.statistics.energy_balance(
                        aggregate_time=False,
                        nice_names=False,
                        groupby=["bus", "carrier", "bus_carrier"],
                    )
                    .loc[:, buses, :, :]
                    .droplevel(["component", "bus"])
                )
                carriers = ["AC", "low voltage"]
                mask = balance.index.get_level_values("bus_carrier").isin(carriers)
                nb = balance[mask].groupby("carrier").sum().T

                for carrier in max_usage:
                    if carrier in nb.columns:
                        max_val = nb[carrier].abs().max()  # in MW, matches caps units
                        if carrier in caps.index:
                            caps.loc[carrier] = max_val

            grouped = {}
            for group_name, techs in capa_groups.items():
                grouped[group_name] = {}
                for tech_name, carrier_list in techs.items():
                    val = caps[caps.index.isin(carrier_list)].sum()
                    grouped[group_name][tech_name] = abs(val) / 1000  # GW

            all_data[(network_id, year)] = grouped

    return all_data


def plot_capacity_comparison(data, network_ids, years, tech_colors, max_usage=None):
    """Create stacked bar chart comparing network_ids across years.

    Parameters
    ----------
    max_usage : list, optional
        List of technology names where max usage is plotted instead of capacity
    """
    if max_usage is None:
        max_usage = []

    fig, axes = plt.subplots(
        len(years),
        4,
        figsize=(16, 4 * len(years)),
        gridspec_kw={"wspace": 0.3, "hspace": 0.4},
    )
    if len(years) == 1:
        axes = axes.reshape(1, -1)

    group_names = list(capa_groups.keys())
    x = np.arange(len(network_ids))
    width = 0.6

    for year_idx, year in enumerate(years):
        for group_idx, group_name in enumerate(group_names):
            ax = axes[year_idx, group_idx]
            tech_names = list(capa_groups[group_name].keys())
            bottoms = np.zeros(len(network_ids))

            for tech_name in tech_names:
                values = [
                    data[(nid, year)][group_name].get(tech_name, 0) for nid in network_ids
                ]
                color = tech_colors.get(tech_name, "#CCCCCC")
                hatch = None
                if "CHP" in tech_name:
                    hatch = "//"
                elif "decentral" in tech_name:
                    hatch = "\\\\"

                label = tech_name + "*" if tech_name in max_usage else tech_name

                bars = ax.bar(
                    x,
                    values,
                    width,
                    bottom=bottoms,
                    color=color,
                    label=label,
                    hatch=hatch,
                    edgecolor="white" if hatch else None,
                    linewidth=0.5,
                )

                for i, (bar, val) in enumerate(zip(bars, values)):
                    if val > 1:
                        ax.text(
                            bar.get_x() + bar.get_width() / 2,
                            bottoms[i] + val / 2,
                            f"{round(val)}",
                            ha="center",
                            va="center",
                            fontsize=9,
                            color="white",
                            weight="bold",
                        )
                bottoms += values

            for i, total in enumerate(bottoms):
                if total > 0:
                    ax.text(
                        i,
                        total + max(bottoms) * 0.02,
                        f"{round(total)}",
                        ha="center",
                        va="bottom",
                        fontsize=10,
                        weight="bold",
                    )

            ax.set_xticks(x)
            ax.set_xticklabels(
                [scenario_label(nid) for nid in network_ids],
                rotation=45 if len(network_ids) > 5 else 0,
                ha="right" if len(network_ids) > 5 else "center",
            )
            ax.set_ylabel("Installed capacity (GW)", fontsize=10)
            ax.set_ylim(0, max(bottoms) * 1.15 if max(bottoms) > 0 else 1)
            ax.grid(axis="y", alpha=0.3)

            if year_idx == 0:
                ax.set_title(group_name, fontsize=12, weight="bold")
            if group_idx == len(group_names) - 1:
                ax.text(
                    1.05,
                    0.5,
                    str(year),
                    transform=ax.transAxes,
                    rotation=270,
                    va="center",
                    fontsize=14,
                    weight="bold",
                )
            if len(network_ids) <= 5:
                ax.set_xlabel("Network", fontsize=10)

            if year_idx == len(years) - 1:
                handles, labels = ax.get_legend_handles_labels()

                ax.legend(
                    handles,
                    labels,
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.15),
                    ncol=1,
                    fontsize=9,
                    frameon=False,
                )

                if group_idx == len(group_names) - 1 and max_usage:
                    ax.text(
                        0.1,
                        -0.65,
                        "* Maximum capacity usage instead of installed capacity",
                        transform=ax.transAxes,
                        ha="center",
                        fontsize=10,
                        style="italic",
                    )

    plt.tight_layout()
    return fig


def plot_energy_system_cost_comparison(
    df,
    network_ids,
    sectors=("Electricity", "Heat", "H2", "Fuels", "Gas", "Biomass", "Other"),
    colors=None,
    save_path=None,
    figsize=(9, 5),
):
    """
    Plot total energy system cost comparison across network_ids.

    Left bars show costs stacked by type (CAPEX, OPEX, Trade).
    Right bars show costs stacked by sector within each cost type.
    """
    if colors is None:
        colors = sector_colors

    capex = np.array([df.filter(like="CAPEX", axis=0)[nid].sum() for nid in network_ids])
    opex = np.array([df.filter(like="OPEX", axis=0)[nid].sum() for nid in network_ids])
    trade = np.array([df.filter(like="Trade", axis=0)[nid].sum() for nid in network_ids])
    total = capex + opex + trade

    sector_capex = np.zeros((len(network_ids), len(sectors)))
    sector_opex = np.zeros((len(network_ids), len(sectors)))
    sector_trade = np.zeros((len(network_ids), len(sectors)))

    for i, nid in enumerate(network_ids):
        for j, sector in enumerate(sectors):
            sector_capex[i, j] = df[df.index.str.contains(f"CAPEX.*{sector}")][nid].sum()
            sector_opex[i, j] = df[df.index.str.contains(f"OPEX.*{sector}")][nid].sum()
            sector_trade[i, j] = df[df.index.str.contains(f"Trade.*{sector}")][nid].sum()

    base_idx = network_ids.index("Base") if "Base" in network_ids else None
    if base_idx is not None:
        base_total = total[base_idx]
        pct_diff = (total - base_total) / base_total * 100
    else:
        pct_diff = np.zeros(len(network_ids))

    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(network_ids))
    width = 0.35

    ax.bar(x - width / 2, capex, width, color="thistle")
    ax.bar(x - width / 2, opex, width, bottom=capex, color="pink")
    ax.bar(x - width / 2, trade, width, bottom=capex + opex, color="navajowhite")

    for i in range(len(network_ids)):
        ax.text(i - width / 2, capex[i] / 2, "CAPEX", ha="center", va="center", fontsize=9, fontweight="bold", color="white")
        ax.text(i - width / 2, capex[i] + opex[i] / 2, "OPEX", ha="center", va="center", fontsize=9, fontweight="bold", color="white")
        ax.text(i - width / 2, capex[i] + opex[i] + trade[i] / 2, "Trade", ha="center", va="center", fontsize=9, fontweight="bold", color="white")

    for i in range(len(network_ids)):
        ax.text(i - width / 2, total[i], f"{total[i]:.1f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    bottom = np.zeros(len(network_ids))

    for j, sector in enumerate(sectors):
        ax.bar(x + width / 2, sector_capex[:, j], width, bottom=bottom, color=colors[sector], label=sector, edgecolor="white", linewidth=0.5)
        bottom += sector_capex[:, j]

    capex_top = bottom.copy()

    for j, sector in enumerate(sectors):
        ax.bar(x + width / 2, sector_opex[:, j], width, bottom=bottom, color=colors[sector], edgecolor="white", linewidth=0.5)
        bottom += sector_opex[:, j]

    opex_top = bottom.copy()

    for j, sector in enumerate(sectors):
        ax.bar(x + width / 2, sector_trade[:, j], width, bottom=bottom, color=colors[sector], edgecolor="white", linewidth=0.5)
        bottom += sector_trade[:, j]

    for i in range(len(network_ids)):
        ax.hlines(capex_top[i], i + width / 2 - width / 2, i + width / 2 + width / 2, colors="black", linewidth=2.5, zorder=10)
        ax.hlines(opex_top[i], i + width / 2 - width / 2, i + width / 2 + width / 2, colors="black", linewidth=2.5, zorder=10)

    if base_idx is not None:
        for i, network_id in enumerate(network_ids):
            if network_id != "Base":
                sign = "+" if pct_diff[i] > 0 else ""
                ax.text(i + width / 2, total[i], f"{sign}{pct_diff[i]:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_ylabel("Cost (billion EUR)", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels([scenario_label(nid) for nid in network_ids], rotation=30, ha="right")
    ax.legend(fontsize=9, ncol=4, loc="upper right", title="Sectors")
    ax.grid(axis="y", alpha=0.3)
    if max(total) > 0:
        ax.set_ylim(0, max(total) * 1.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return fig, ax


def make_capacity_tables(networks, network_ids, output_dir, kwargs):
    """Save capacity_table_{bus_carrier}.png for a set of network_ids.

    Split out of the main comparison body so it can also be run, with a
    single network_id, for each earlier myopic-pathway year mixed into the
    LT comparison (see the "optimal_YYYY" pseudo network_ids and the
    pathway-year tables section in `__main__`).
    """
    carriers_sets = [
        ["AC", "low voltage"],
        ["urban central heat", "rural heat", "urban decentral heat"],
        ["H2"],
        ["co2 stored"],
        ["oil"],
        ["renewable oil"],
        ["gas"],
        ["renewables gas"],
        ["solid biomass"],
    ]

    for bc in carriers_sets:
        capa = {}
        for network_id in network_ids:
            n = networks.get(network_id)
            if n is None:
                continue
            stats = n.statistics.optimal_capacity(bus_carrier=bc, **kwargs).filter(like="DE")
            if stats.empty:
                capa[network_id] = pd.Series(dtype=float)
            else:
                capa[network_id] = stats.groupby(["carrier"]).sum().div(1e3)  # GW

        if not capa:
            continue
        capa_df = pd.concat(capa, axis=1)
        capa_df.columns = list(capa.keys())

        if not capa_df.empty:
            df = round(capa_df[capa_df.gt(1).any(axis=1)], 2)
            if not df.empty:
                df_to_png(
                    df.rename(columns=scenario_label),
                    f"{output_dir}/capacity_table_{'-'.join(bc)}.png",
                )


def make_generation_consumption_tables(networks, network_ids, output_dir, kwargs_gc, region="DE"):
    """Save generation/consumption/heat/re-fuel tables for a set of
    network_ids (see make_capacity_tables docstring for why this is split
    out of the main comparison body)."""
    for carrier_name, config in CARRIER_GROUPS.items():
        gen_df, con_df = get_generation_consumption(
            networks, network_ids, carrier_name, config, region, kwargs_gc
        )

        if not gen_df.empty and gen_df.gt(config["threshold"]).any().any():
            df = round(gen_df[gen_df.gt(config["threshold"]).any(axis=1)], 2)
            df_to_png(
                df.rename(columns=scenario_label),
                os.path.join(output_dir, f"{carrier_name}_generation.png"),
            )

        if not con_df.empty and con_df.gt(config["threshold"]).any().any():
            df = round(con_df[con_df.gt(config["threshold"]).any(axis=1)], 2)
            df_to_png(
                df.rename(columns=scenario_label),
                os.path.join(output_dir, f"{carrier_name}_consumption.png"),
            )

    heat_levels = CARRIER_GROUPS["heat"]["bus_carrier"]
    for bc in heat_levels:
        heat_gen = {}
        for network_id in network_ids:
            n = networks.get(network_id)
            if n is None:
                continue
            heat_gen[network_id] = (
                n.statistics.supply(bus_carrier=[bc], **kwargs_gc)
                .filter(like=region)
                .groupby(["carrier"])
                .sum()
                .div(1e6)
            )
        if not heat_gen:
            continue
        heat_gen_df = pd.concat(heat_gen, axis=1)
        heat_gen_df.columns = list(heat_gen.keys())

        if not heat_gen_df.empty and heat_gen_df.gt(0.1).any().any():
            df = round(heat_gen_df[heat_gen_df.gt(0.1).any(axis=1)], 2)
            filename = bc.replace(" ", "_")
            df_to_png(
                df.rename(columns=scenario_label),
                os.path.join(output_dir, f"heat_gen_{filename}.png"),
            )

    re_oil_config = {"bus_carrier": ["renewable oil"], "drop_stores": False, "drop_carriers": [], "threshold": 1.0}
    re_oil_gen_df, _ = get_generation_consumption(networks, network_ids, "renewable_oil", re_oil_config, region, kwargs_gc)
    if not re_oil_gen_df.empty and re_oil_gen_df.gt(1.0).any().any():
        df = round(re_oil_gen_df[re_oil_gen_df.gt(1.0).any(axis=1)], 2)
        df_to_png(df.rename(columns=scenario_label), os.path.join(output_dir, "re_oil_generation.png"))

    re_gas_config = {"bus_carrier": ["renewable gas"], "drop_stores": False, "drop_carriers": [], "threshold": 1.0}
    re_gas_gen_df, _ = get_generation_consumption(networks, network_ids, "renewable_gas", re_gas_config, region, kwargs_gc)
    if not re_gas_gen_df.empty and re_gas_gen_df.gt(1.0).any().any():
        df = round(re_gas_gen_df[re_gas_gen_df.gt(1.0).any(axis=1)], 2)
        df_to_png(df.rename(columns=scenario_label), os.path.join(output_dir, "re_gas_generation.png"))


def make_system_cost_tables(networks, network_ids, output_dir, region="DE"):
    """Compute the system-cost breakdown for network_ids and save
    system_cost.png / system_cost_plotted.png / system_cost_{sector}.png
    (see make_capacity_tables docstring for why this is split out of the
    main comparison body). Returns (var_all, var_plot, tsc_all) for callers
    that also want the total_energy_system_cost_by_scenario.png bar chart.
    """
    var_all = pd.DataFrame(columns=network_ids)
    tsc_all = {}

    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        var = pd.Series(dtype=float)

        var["Total Energy System Cost|Trade|Electricity"] = (
            get_trade_cost(n, region, ["AC"]) + get_trade_cost(n, region, ["DC"])
        ) / 1e9
        var["Total Energy System Cost|Trade|eFuels"] = (
            get_trade_cost(n, region, ["renewable oil", "renewable gas", "methanol"]) / 1e9
        )
        var["Total Energy System Cost|Trade|H2"] = (
            get_trade_cost(
                n, region, ["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"]
            )
            / 1e9
        )
        tsc_sum, tsc = get_tsc(n, region)

        var["Total Energy System Cost|EU"] = (
            n.statistics.capex().sum() + n.statistics.opex().sum()
        ) / 1e9

        biomass_import_ind = tsc[tsc.index.str.contains("biomass transported")].index
        var["Total Energy System Cost|Trade|Biomass"] = tsc.loc[biomass_import_ind, "opex"].sum().sum() / 1e9
        tsc.loc[biomass_import_ind, "opex"] = 0

        gas_import_ind = tsc[tsc.carrier == "gas primary"].index
        var["Total Energy System Cost|Trade|Gas"] = tsc.loc[gas_import_ind, "opex"].sum().sum() / 1e9
        tsc.loc[gas_import_ind, "opex"] = 0

        oil_import_ind = tsc[tsc.carrier == "oil primary"].index
        var["Total Energy System Cost|Trade|Oil"] = tsc.loc[oil_import_ind, "opex"].sum().sum() / 1e9
        tsc.loc[oil_import_ind, "opex"] = 0

        tsc_all[network_id] = tsc

        lignite_ind = n.links[(n.links.carrier == "lignite") & (n.links.bus1.str.contains("DE"))].index
        lignite_cost = (
            (n.links_t.p0[lignite_ind].sum(axis=1) * n.buses_t.marginal_price["EU lignite"]).sum() / 1e9
            if len(lignite_ind)
            else 0.0
        )

        coal_ind = n.links[(n.links.carrier == "coal") & (n.links.bus1.str.contains("DE"))].index
        coal_cost = (
            (n.links_t.p0[coal_ind].sum(axis=1) * n.buses_t.marginal_price["EU coal"]).sum() / 1e9
            if len(coal_ind)
            else 0.0
        )

        var["Total Energy System Cost|Trade|Coal"] = lignite_cost + coal_cost

        var["Total Energy System Cost"] = (
            tsc_sum.sum().sum() / 1e9
            + var["Total Energy System Cost|Trade|Electricity"]
            + var["Total Energy System Cost|Trade|eFuels"]
            + var["Total Energy System Cost|Trade|H2"]
            + var["Total Energy System Cost|Trade|Coal"]
        )

        var["Total Energy System Cost|Trade"] = (
            var["Total Energy System Cost|Trade|Electricity"]
            + var["Total Energy System Cost|Trade|eFuels"]
            + var["Total Energy System Cost|Trade|H2"]
            + var["Total Energy System Cost|Trade|Biomass"]
            + var["Total Energy System Cost|Trade|Gas"]
            + var["Total Energy System Cost|Trade|Oil"]
            + var["Total Energy System Cost|Trade|Coal"]
        )

        var["Total Energy System Cost|Trade|Fuels"] = (
            var["Total Energy System Cost|Trade|eFuels"] + var["Total Energy System Cost|Trade|Oil"]
        )

        var["Total Energy System Cost|Non Trade"] = tsc[["capex", "opex"]].sum().sum() / 1e9

        assert (
            abs(
                var["Total Energy System Cost"]
                - (var["Total Energy System Cost|Non Trade"] + var["Total Energy System Cost|Trade"])
            )
            < 1e-6
        )

        df = tsc
        df["bus"] = df.apply(lambda row: get_component_bus(row, n), axis=1)
        df["sector"] = df["bus"].apply(lambda x: map_bus_to_sector(x, n, BUS_CARRIER_TO_SECTOR))

        sector_costs = df.groupby("sector")[["capex", "opex"]].sum()

        for sector in ["Electricity", "Heat", "H2", "Fuels", "Gas", "Biomass"]:
            var["Total Energy System Cost|CAPEX|" + sector] = sector_costs["capex"].get(sector, 0.0) / 1e9
            var["Total Energy System Cost|OPEX|" + sector] = sector_costs["opex"].get(sector, 0.0) / 1e9

        other = [s for s in sector_costs.index if s not in ["Electricity", "Heat", "H2", "Fuels", "Gas", "Biomass"]]
        var["Total Energy System Cost|CAPEX|Other"] = sector_costs.loc[other, "capex"].sum() / 1e9
        var["Total Energy System Cost|OPEX|Other"] = sector_costs.loc[other, "opex"].sum() / 1e9

        var["Total Energy System Cost|Trade|Other"] = var["Total Energy System Cost|Trade|Coal"]

        var_all[network_id] = var

    df_to_png(round(var_all, 2).rename(columns=scenario_label), output_dir / "system_cost.png")

    index_keep = [
        "Total Energy System Cost|CAPEX|Electricity",
        "Total Energy System Cost|CAPEX|Heat",
        "Total Energy System Cost|CAPEX|H2",
        "Total Energy System Cost|CAPEX|Fuels",
        "Total Energy System Cost|CAPEX|Gas",
        "Total Energy System Cost|CAPEX|Biomass",
        "Total Energy System Cost|CAPEX|Other",
        "Total Energy System Cost|OPEX|Electricity",
        "Total Energy System Cost|OPEX|Heat",
        "Total Energy System Cost|OPEX|H2",
        "Total Energy System Cost|OPEX|Fuels",
        "Total Energy System Cost|OPEX|Gas",
        "Total Energy System Cost|OPEX|Biomass",
        "Total Energy System Cost|OPEX|Other",
        "Total Energy System Cost|Trade|Electricity",
        "Total Energy System Cost|Trade|H2",
        "Total Energy System Cost|Trade|Fuels",
        "Total Energy System Cost|Trade|Gas",
        "Total Energy System Cost|Trade|Biomass",
        "Total Energy System Cost|Trade|Other",
    ]
    var_plot = pd.DataFrame(index=index_keep, columns=network_ids)

    for s in network_ids:
        var_plot[s] = var_all[s]

    assert np.isclose(var_plot.sum(), var_all.loc["Total Energy System Cost"]).all()

    df_to_png(round(var_plot, 2).rename(columns=scenario_label), output_dir / "system_cost_plotted.png")

    for sector in ["Electricity", "Heat", "H2", "Fuels", "Gas", "Biomass", "Other"]:
        res = pd.DataFrame()
        for network_id in network_ids:
            if network_id not in tsc_all:
                continue
            df = tsc_all[network_id]
            res[network_id] = df[df.sector == sector].groupby(by="carrier").sum()[["capex"]]

        df = round(res / 1e9, 2)
        if not df.empty:
            df_filtered = df[(df > 1).any(axis=1)]
            if not df_filtered.empty:
                df_to_png(
                    df_filtered.rename(columns=scenario_label), output_dir / f"system_cost_{sector}.png"
                )

    return var_all, var_plot, tsc_all


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_grid_scenario_comparison",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
            variant="LT",
        )

    configure_logging(snakemake)

    all_network_ids = list(snakemake.params.network_ids)
    year = int(snakemake.wildcards.planning_horizons)
    variant = snakemake.wildcards.variant
    output_dir = Path(snakemake.output[0])
    output_dir.mkdir(parents=True, exist_ok=True)

    # "optimal_YYYY" pseudo network_ids (see GRID_EXPORT_NETWORK_IDS_LT in
    # stochastic_grid.smk) are the plain, suffix-less pathway-reference
    # networks for years earlier than this comparison's own target `year` -
    # mixed into the LT comparison so the buildout leading up to `year` is
    # visible, but kept out of the main (single-year) tables/plots below and
    # given their own pathway-year tables instead (see the end of this
    # script). ST comparisons never contain these, so pathway_ids is empty
    # and nid_to_year trivially maps every id to `year`.
    nid_to_year = {
        nid: (int(nid.removeprefix("optimal_")) if nid.startswith("optimal_") else year)
        for nid in all_network_ids
    }
    pathway_ids = [nid for nid in all_network_ids if nid.startswith("optimal_")]
    network_ids = [nid for nid in all_network_ids if nid not in pathway_ids]

    logger.info(f"Comparing {variant} network variants: {network_ids}")
    if pathway_ids:
        logger.info(f"Pathway-year reference networks: {pathway_ids}")

    variables = {}
    for path, network_id in zip(snakemake.input.exported_variables, all_network_ids):
        _df = pd.read_excel(path, index_col=list(range(5)), sheet_name="data").droplevel(
            ["Model", "Scenario", "Region"]
        )
        variables[network_id] = _df.droplevel("Unit")[nid_to_year[network_id]]

    logger.info("Loading networks...")
    networks = {}
    for path, network_id in zip(snakemake.input.networks, all_network_ids):
        try:
            networks[network_id] = pypsa.Network(path)
        except Exception as e:
            logger.warning(f"Could not load network for {network_id} ({path}): {e}")

    logger.info(f"Loaded {len(networks)} network variants")

    kwargs = {
        "groupby": pypsa.statistics.groupers["bus", "carrier"],
        "nice_names": False,
    }

    ####### PLOTTING #########

    ## Capacities (bar plots) ##
    carrier_sets = [
        ["AC", "low voltage"],
        ["urban central heat", "rural heat", "urban decentral heat"],
        ["H2"],
        ["oil"],
        ["gas"],
        ["co2 stored"],
        ["co2 sequestered"],
    ]
    region = "DE"

    for bc in carrier_sets:
        try:
            fig, ax = plt.subplots(1, 1, figsize=(10, 15))
            capas = {}

            for network_id in network_ids:
                n = networks.get(network_id)
                if n is None:
                    continue
                capas[network_id] = (
                    n.statistics.optimal_capacity(bus_carrier=bc, **kwargs)
                    .filter(like=region)
                    .groupby("carrier")
                    .sum()
                    .div(1e3)
                )

            if not capas:
                plt.close(fig)
                continue

            df = pd.concat(capas, axis=1)
            df.columns = list(capas.keys())

            if df.empty:
                logger.info(f"Skipping {bc}: No data")
                plt.close(fig)
                continue

            df.plot(kind="barh", title=f"Optimal capacity ({bc})", ylabel="GW", ax=ax)
            plt.grid()
            plt.tight_layout()
            plt.savefig(f"{output_dir}/capacity_{'_'.join(bc)}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)

        except (IndexError, KeyError) as e:
            logger.info(f"Skipping {bc}: {e}")
            plt.close(fig)
            continue

    # CAPACITIES TABLES
    make_capacity_tables(networks, network_ids, output_dir, kwargs)

    # CAPACITY COMPARISON (stacked, grouped by tech category) - one row per
    # pathway year (just `year` for ST, `year` plus every earlier pathway
    # year for LT), each row only populated for the network_id(s) that
    # actually represent that year (see get_capacities' year_of parameter).
    logger.info("Plotting capacity comparison for DE")
    capacity_years = sorted(set(nid_to_year.values()))
    capacity_data = get_capacities(
        networks, all_network_ids, capacity_years, year_of=nid_to_year, max_usage=["BEV charger", "V2G"]
    )
    fig = plot_capacity_comparison(
        capacity_data,
        all_network_ids,
        capacity_years,
        tech_colors,
        max_usage=["Vehicle-to-grid", "BEV charging"],
    )
    fig.savefig(output_dir / "capacity_comparison.png", dpi=300, bbox_inches="tight")

    # REGIONAL CAPACITY-DIFF MAPS (LT only) - each of the topology-*
    # variants against the unconstrained "optimal" network, plus each
    # stochastic (topology-stochastic__{gs}) variant against its own
    # deterministic counterpart (topology-{gs}), for every capa_groups /
    # storage_energy_capa_groups category (see utils.py).
    if variant == "LT":
        logger.info("Plotting regional capacity-diff maps...")
        onshore_regions = gpd.read_file(snakemake.input.regions_onshore[0]).set_index("name")

        diff_pairs = []
        if "optimal" in network_ids:
            diff_pairs += [(nid, "optimal") for nid in network_ids if nid != "optimal"]
        for nid in network_ids:
            if nid.startswith("topology-stochastic__"):
                det_id = "topology-" + nid.removeprefix("topology-stochastic__")
                if det_id in network_ids:
                    diff_pairs.append((nid, det_id))

        for nid_a, nid_b in diff_pairs:
            n_a, n_b = networks.get(nid_a), networks.get(nid_b)
            if n_a is None or n_b is None:
                continue
            label_a, label_b = scenario_label(nid_a), scenario_label(nid_b)

            power_a = power_capacity_by_bus_carrier(n_a, region)
            power_b = power_capacity_by_bus_carrier(n_b, region)
            energy_a = storage_energy_capacity_by_bus_carrier(n_a, region)
            energy_b = storage_energy_capacity_by_bus_carrier(n_b, region)

            for group_name, tech_dict in capa_groups.items():
                group_slug = group_name.lower().replace(" + ", "_").replace(" ", "_")
                plot_group_capacity_diff_maps(
                    group_capacity_by_tech(n_a, tech_dict, power_a),
                    group_capacity_by_tech(n_b, tech_dict, power_b),
                    group_name,
                    onshore_regions,
                    output_dir / f"capacity_diff_map_{group_slug}_{label_a}_vs_{label_b}.png",
                    label_a,
                    label_b,
                    region=region,
                    unit="GW",
                )

            for group_name, tech_dict in storage_energy_capa_groups.items():
                group_slug = group_name.lower().replace(" + ", "_").replace(" ", "_")
                plot_group_capacity_diff_maps(
                    group_capacity_by_tech(n_a, tech_dict, energy_a),
                    group_capacity_by_tech(n_b, tech_dict, energy_b),
                    group_name,
                    onshore_regions,
                    output_dir / f"capacity_diff_map_{group_slug}_{label_a}_vs_{label_b}.png",
                    label_a,
                    label_b,
                    region=region,
                    unit="GWh",
                )

    ### GENERATION & CONSUMPTION ###

    kwargs_gc = {
        "groupby": pypsa.statistics.groupers["name", "bus", "carrier"],
        "at_port": True,
        "nice_names": False,
    }

    make_generation_consumption_tables(networks, network_ids, output_dir, kwargs_gc, region)

    # STORAGE PLOTS
    heat_stores_carriers = [
        "urban central water tanks",
        "urban central water pits",
        "urban decentral water tanks",
        "urban decentral water pits",
        "rural water tanks",
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    for network_id in network_ids:
        n = networks.get(network_id)
        if n is None:
            continue
        stores_i = n.stores[
            (n.stores.carrier.isin(heat_stores_carriers)) & (n.stores.bus.str.contains(region))
        ].index

        if len(stores_i) > 0:
            energy_capa = n.stores.loc[stores_i].e_nom_opt.sum() / 1e6
            (n.stores_t.e[stores_i].sum(axis=1) / 1e6).plot(
                label=f"{scenario_label(network_id)} (Capacity: {energy_capa:.2f} TWh)", ax=ax
            )

    ax.set_xlabel("Time [h]")
    ax.set_ylabel("Energy [TWh]")
    ax.set_title("Heat Storage")
    ax.legend(fontsize=8)
    plt.savefig(os.path.join(output_dir, "heat_storage.png"), dpi=300, bbox_inches="tight")
    plt.close()

    plot_storage_timeseries(networks, network_ids, "H2 Store", region, os.path.join(output_dir, "h2_storage.png"))
    plot_storage_timeseries(networks, network_ids, "oil", region, os.path.join(output_dir, "oil_storage.png"))
    plot_storage_timeseries(networks, network_ids, "gas", region, os.path.join(output_dir, "gas_storage.png"))

    logger.info(f"All generation and consumption plots saved to: {output_dir}")

    ### PRICES ###
    carriers_sets = [
        ["AC", "low voltage"],
        ["urban central heat", "rural heat", "urban decentral heat"],
        ["H2"],
        ["co2 stored"],
        ["oil"],
        ["renewable oil"],
        ["gas"],
        ["renewables gas"],
        ["solid biomass"],
    ]

    for carriers in carriers_sets:
        try:
            plot_marginal_prices(networks, network_ids, carriers, output_dir, variant)
        except Exception as e:
            logger.warning(f"Failed to plot marginal prices for carriers: {carriers} - {e}")

    ### PRICES - PAPER ###

    price_carriers_to_plot = [
        ["AC", "low voltage"],
        ["urban central heat"],
        ["urban decentral heat"],
        ["rural heat"],
        ["urban central heat", "urban decentral heat", "rural heat"],
        ["H2"],
        ["gas"],
        ["renewable gas"],
        ["methanol"],
        ["renewable oil"],
        ["oil"],
    ]

    for c in price_carriers_to_plot:
        plot_price_duration_curves(
            networks,
            network_ids,
            carriers=c,
            regions=["DE"],
            output_dir=output_dir / f"pdc_{'_'.join(c)}.png",
            ylim=(-50, 500),
            table_smaller=1,
            table_bigger=200,
        )

    ### SYSTEM COSTS ###
    logger.info("Plotting system costs...")

    var_all, var_plot, tsc_all = make_system_cost_tables(networks, network_ids, output_dir, region)

    var_plot_renamed = var_plot.rename(
        index=lambda i: i.replace("Total Energy System Cost|", "").replace("|", " ")
    )
    fig, ax = plot_energy_system_cost_comparison(
        var_plot_renamed,
        network_ids=network_ids,
        save_path=output_dir / "total_energy_system_cost_by_scenario.png",
        colors=sector_colors,
    )

    ### TRADE ###
    logger.info("Plotting trade...")
    plot_vars = {
        "Electricity": "Trade|Secondary Energy|Electricity|Volume",
        "Gas": "Primary Energy|Gas",
        "Oil": "Primary Energy|Oil",
        "Hydrogen": "Trade|Secondary Energy|Hydrogen|Volume",
        "eFuels": "Trade|Secondary Energy|Efuels|Volume",
        "Biomass": "Trade|Primary Energy|Biomass|Net Imports",
    }

    bar_plot_variables(
        variables,
        network_ids=network_ids,
        tech_colors=dict(sector_colors, **tech_colors),
        plot_vars=plot_vars,
        sign_flip_vars=["Oil", "Gas", "Biomass"],
        output_dir=output_dir,
        filename_prefix="trade_volume",
    )

    plot_vars = {
        "eFuels": "Trade|Secondary Energy|Efuels|Volume",
        "Renewable Gas": "Trade|Secondary Energy|Efuels|Renewable Gas|Volume",
        "Renewable Oil": "Trade|Secondary Energy|Efuels|Renewable Oil|Volume",
        "Methanol": "Trade|Secondary Energy|Efuels|Methanol|Volume",
    }

    bar_plot_variables(
        variables,
        network_ids=network_ids,
        tech_colors=dict(sector_colors, **tech_colors),
        plot_vars=plot_vars,
        sign_flip_vars=["Oil"],
        output_dir=output_dir,
        filename_prefix="trade_volume_efuels",
    )

    plot_vars = {
        "Electricity": "Total Energy System Cost|Trade|Electricity",
        "Gas": "Total Energy System Cost|Trade|Gas",
        "Oil": "Total Energy System Cost|Trade|Oil",
        "Hydrogen": "Total Energy System Cost|Trade|Hydrogen",
        "eFuels": "Total Energy System Cost|Trade|Efuels",
        "Biomass": "Total Energy System Cost|Trade|Biomass",
    }

    bar_plot_variables(
        variables,
        network_ids=network_ids,
        tech_colors=dict(sector_colors, **tech_colors),
        plot_vars=plot_vars,
        sign_flip_vars=[],
        output_dir=output_dir,
        filename_prefix="trade_cost",
    )

    plot_vars = {
        "eFuels": "Total Energy System Cost|Trade|Efuels",
        "Renewable Gas": "Total Energy System Cost|Trade|Efuels|Renewable Gas",
        "Renewable Oil": "Total Energy System Cost|Trade|Efuels|Renewable Oil",
        "Methanol": "Total Energy System Cost|Trade|Efuels|Methanol",
    }

    bar_plot_variables(
        variables,
        network_ids=network_ids,
        tech_colors=dict(sector_colors, **tech_colors),
        plot_vars=plot_vars,
        sign_flip_vars=[],
        output_dir=output_dir,
        filename_prefix="trade_cost_efuels",
    )

    ### CURTAILMENT ###
    logger.info("Plotting curtailment...")
    plot_curtailment(networks, network_ids, tech_colors, output_dir)

    ### PATHWAY-YEAR TABLES ###
    # Earlier myopic-pathway years (e.g. 2025 relative to a 2030 target)
    # aren't mixed into the tables above (their absolute capacity/cost
    # levels aren't directly comparable to the target year's), but do get
    # their own, single-network table set in a same-named subfolder - see
    # "optimal_YYYY" in GRID_EXPORT_NETWORK_IDS_LT (stochastic_grid.smk).
    # capacity_comparison.png is the one exception: it already shows every
    # pathway year as its own row (see the capacity_years block above).
    for pid in pathway_ids:
        py = nid_to_year[pid]
        py_dir = output_dir / str(py)
        py_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Building {py} pathway-year tables ({pid}) -> {py_dir}")
        make_capacity_tables(networks, [pid], py_dir, kwargs)
        make_generation_consumption_tables(networks, [pid], py_dir, kwargs_gc, region)
        make_system_cost_tables(networks, [pid], py_dir, region)

    logger.info(f"All plots saved to: {output_dir}")
