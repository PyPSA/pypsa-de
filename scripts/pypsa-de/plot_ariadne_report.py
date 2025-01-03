# -*- coding: utf-8 -*-
import logging
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from datetime import datetime
from itertools import compress
from multiprocessing import Pool

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from matplotlib.lines import Line2D
from pypsa.plot import add_legend_lines

from scripts._helpers import configure_logging, mock_snakemake, set_scenario_config
from scripts.export_ariadne_variables import get_discretized_value, process_postnetworks
from scripts.plot_power_network import load_projection
from scripts.plot_summary import preferred_order, rename_techs
from scripts.prepare_sector_network import prepare_costs
from scripts.pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

logger = logging.getLogger(__name__)

####### definitions #######
THRESHOLD = 5  # GW

CARRIER_GROUPS = {
    "electricity": ["AC", "low voltage"],
    "heat": [
        "urban central heat",
        "urban decentral heat",
        "rural heat",
        "residential urban decentral heat",
        "residential rural heat",
        "services urban decentral heat",
        "services rural heat",
    ],
    # "hydrogen": "H2",
    # "oil": "oil",
    # "methanol": "methanol",
    # "ammonia": "NH3",
    # "biomass": ["solid biomass", "biogas"],
    # "CO2 atmosphere": "co2",
    # "CO2 stored": "co2 stored",
    # "methane": "gas",
}

year_colors = [
    "dimgrey",
    "darkorange",
    "seagreen",
    "cadetblue",
    "hotpink",
    "darkviolet",
]
markers = [
    "v",
    "^",
    "<",
    ">",
    "1",
    "2",
    "3",
    "4",
    "*",
    "+",
    "d",
    "o",
    "|",
    "s",
    "P",
    "p",
    "h",
]

date_format = "%Y-%m-%d %H:%M:%S"
reduced_date_format = "%Y-%m-%d"


resistive_heater = [
    "urban central resistive heater",
    "rural resistive heater",
    "urban decentral resistive heater",
]
gas_boiler = [
    "urban central gas boiler",
    "rural gas boiler",
    "urban decentral gas boiler",
]
air_heat_pump = [
    "urban central air heat pump",
    "rural air heat pump",
    "urban decentral air heat pump",
]
water_tanks_charger = [
    "urban central water tanks charger",
    "rural water tanks charger",
    "urban decentral water tanks charger",
]
water_tanks_discharger = [
    "urban central water tanks discharger",
    "rural water tanks discharger",
    "urban decentral water tanks discharger",
]
solar_thermal = [
    "urban decentral solar thermal",
    "urban central solar thermal",
    "rural solar thermal",
]
carrier_renaming = {
    "urban central solid biomass CHP CC": "biomass CHP CC",
    "urban central solid biomass CHP": "biomass CHP",
    "urban central gas CHP": "gas CHP",
    "urban central gas CHP CC": "gas CHP CC",
    "urban central air heat pump": "air heat pump",
    "urban central resistive heater": "resistive heater",
}
carrier_renaming_reverse = {
    "biomass CHP CC": "urban central solid biomass CHP CC",
    "biomass CHP": "urban central solid biomass CHP",
    "gas CHP": "urban central gas CHP",
    "gas CHP CC": "urban central gas CHP CC",
    "air heat pump": "urban central air heat pump",
    "resistive heater": "urban central resistive heater",
}
c1_groups = [
    resistive_heater,
    gas_boiler,
    air_heat_pump,
    water_tanks_charger,
    water_tanks_discharger,
    solar_thermal,
]
c1_groups_name = [
    "resistive heater",
    "gas boiler",
    "air heat pump",
    "water tanks charger",
    "water tanks discharger",
    "solar thermal",
]
solar = [
    "solar rooftop",
    "solar-hsat",
    "solar",
]
electricity_load = [
    "electricity",
    "industry electricity",
    "agriculture electricity",
]
electricity_imports = [
    "AC",
    "DC",
]


####### functions #######
def get_condense_sum(df, groups, groups_name, return_original=False):
    """
    return condensed df, that has been groupeb by condense groups
    Arguments:
        df: df you want to condense (carriers have to be in the columns)
        groups: group labels you want to condense on
        groups_name: name of the new grouped column
        return_original: boolean to specify if the original df should also be returned
    Returns:
        condensed df
    """
    result = df

    for group, name in zip(groups, groups_name):
        # check if carrier are in columns
        bool = [c in df.columns for c in group]
        # updated to carriers within group that are in columns
        group = list(compress(group, bool))

        result[name] = df[group].sum(axis=1)
        result.drop(group, axis=1, inplace=True)

    if return_original:
        return result, df

    return result


def plot_nodal_balance(
    network,
    nodal_balance,
    tech_colors,
    savepath,
    carriers=["AC", "low voltage"],
    start_date="2019-01-01 00:00:00",
    end_date="2019-12-31 00:00:00",
    regions=["DE"],
    model_run="Model run",
    c1_groups=c1_groups,
    c1_groups_name=c1_groups_name,
    loads=["electricity", "industry electricity", "agriculture electricity"],
    plot_lmps=True,
    plot_loads=True,
    resample=None,
    nice_names=False,
    threshold=1e-3,  # in GWh
    condense_groups=None,
    condense_names=None,
    ylabel="total electricity balance [GW]",
    title="Electricity balance",
):
    carriers = carriers
    loads = loads
    start_date = start_date
    end_date = end_date
    regions = regions
    period = network.generators_t.p.index[
        (network.generators_t.p.index >= start_date)
        & (network.generators_t.p.index <= end_date)
    ]
    # ToDo find out why this is overwriting itself
    rename = {}

    mask = nodal_balance.index.get_level_values("bus_carrier").isin(carriers)
    nb = balance[mask].groupby("carrier").sum().div(1e3).T.loc[period]
    if plot_loads:
        df_loads = abs(nb[loads].sum(axis=1))
    # condense groups (summarise carriers to groups)
    nb = get_condense_sum(nb, c1_groups, c1_groups_name)
    # rename unhandy column names
    nb.rename(columns=carrier_renaming, inplace=True)
    # summarise some carriers if specified
    if condense_groups is not None:
        nb = get_condense_sum(nb, condense_groups, condense_names)

    ## summaris low contributing carriers according to their sum over the period (threshold in GWh)
    techs_below_threshold = nb.columns[nb.abs().sum() < threshold].tolist()
    if techs_below_threshold:
        other = {tech: "other" for tech in techs_below_threshold}
        rename.update(other)
        tech_colors["other"] = "grey"

    if rename:
        nb = nb.T.groupby(nb.columns.map(lambda a: rename.get(a, a))).sum().T

    if resample is not None:
        nb = nb.resample(resample).mean()

    df = nb
    # split into df with positive and negative values
    df_neg, df_pos = df.clip(upper=0), df.clip(lower=0)
    df_pos = df_pos[df_pos.sum().sort_values(ascending=False).index]
    df_neg = df_neg[df_neg.sum().sort_values().index]
    # get colors
    c_neg, c_pos = (
        [tech_colors[col] for col in df_neg.columns],
        [tech_colors[col] for col in df_pos.columns],
    )

    fig, ax = plt.subplots(figsize=(14, 8))

    # plot positive values
    ax = df_pos.plot.area(ax=ax, stacked=True, color=c_pos, linewidth=0.0)

    # rename negative values that are also present on positive side, so that they are not shown and plot negative values
    f = lambda c: "out_" + c
    cols = [f(c) if (c in df_pos.columns) else c for c in df_neg.columns]
    cols_map = dict(zip(df_neg.columns, cols))
    ax = df_neg.rename(columns=cols_map).plot.area(
        ax=ax, stacked=True, color=c_neg, linewidth=0.0
    )

    # plot lmps
    if plot_lmps:
        lmps = network.buses_t.marginal_price[
            network.buses[network.buses.carrier.isin(carriers)].index
        ].mean(axis=1)[period]
        ax2 = lmps.plot(
            style="--", color="black", label="lmp (mean over buses)", secondary_y=True
        )
        ax2.grid(False)
        # set limits of secondary y-axis
        ax2.set_ylim(
            [
                -1.5
                * lmps.max()
                * abs(df_neg.sum(axis=1).min())
                / df_pos.sum(axis=1).max(),
                1.5 * lmps.max(),
            ]
        )
        ax2.legend(title="Legend for right y-axis", loc="upper right")
        ax2.set_ylabel("lmp [€/MWh]")

    # plot loads
    if plot_loads:
        df_loads.plot(style=":", color="black", label="loads")

    # explicitly filter out duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    filtered_handles_labels = [
        (h, l) for h, l in zip(handles, labels) if not l.startswith("out_")
    ]
    handles, labels = zip(*filtered_handles_labels)

    if nice_names:
        nice_names_dict = network.carriers.nice_name.to_dict()
        labels = [nice_names_dict.get(l, l) for l in labels]

    # rescale the y-axis
    ax.set_ylim([1.05 * df_neg.sum(axis=1).min(), 1.05 * df_pos.sum(axis=1).max()])
    ax.legend(
        handles,
        labels,
        ncol=1,
        loc="upper center",
        bbox_to_anchor=(1.22 if plot_lmps else 1.13, 1.01),
        title="Legend for left y-axis" if plot_lmps else "Legend",
    )

    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.set_title(
        f"{title} ({model_run}; {datetime.strptime(start_date, date_format).strftime(reduced_date_format)} - {datetime.strptime(end_date, date_format).strftime(reduced_date_format)})",
        fontsize=16,
        pad=15,
    )
    ax.grid(True)
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()

    return fig


def plot_stacked_area_steplike(ax, df, colors={}):
    if isinstance(colors, pd.Series):
        colors = colors.to_dict()

    df_cum = df.cumsum(axis=1)

    previous_series = np.zeros_like(df_cum.iloc[:, 0].values)

    for col in df_cum.columns:
        ax.fill_between(
            df_cum.index,
            previous_series,
            df_cum[col],
            step="pre",
            linewidth=0,
            color=colors.get(col, "grey"),
            label=col,
        )
        previous_series = df_cum[col].values


def plot_energy_balance_timeseries(
    df,
    time=None,
    ylim=None,
    resample=None,
    rename={},
    preferred_order=[],
    ylabel="",
    colors={},
    threshold=0,
    dir="",
):
    if time is not None:
        df = df.loc[time]

    timespan = df.index[-1] - df.index[0]
    long_time_frame = timespan > pd.Timedelta(weeks=5)

    techs_below_threshold = df.columns[df.abs().max() < threshold].tolist()

    if techs_below_threshold:
        other = {tech: "other" for tech in techs_below_threshold}
        rename.update(other)
        colors["other"] = "grey"

    if rename:
        df = df.T.groupby(df.columns.map(lambda a: rename.get(a, a))).sum().T

    if resample is not None:
        # upsampling to hourly resolution required to handle overlapping block
        df = df.resample("1h").ffill().resample(resample).mean()

    order = (df / df.max()).var().sort_values().index
    if preferred_order:
        order = preferred_order.intersection(order).append(
            order.difference(preferred_order)
        )
    df = df.loc[:, order]

    # fillna since plot_stacked_area_steplike cannot deal with NaNs
    pos = df.where(df > 0).fillna(0.0)
    neg = df.where(df < 0).fillna(0.0)

    fig, ax = plt.subplots(figsize=(10, 4), layout="constrained")

    plot_stacked_area_steplike(ax, pos, colors)
    plot_stacked_area_steplike(ax, neg, colors)

    plt.xlim((df.index[0], df.index[-1]))

    if not long_time_frame:
        # Set major ticks every Monday
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MONDAY))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%e\n%b"))
        # Set minor ticks every day
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%e"))
    else:
        # Set major ticks every first day of the month
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%e\n%b"))
        # Set minor ticks every 15th of the month
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%e"))

    ax.tick_params(axis="x", which="minor", labelcolor="grey")
    ax.grid(axis="y")

    # half the labels because pos and neg create duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    half = int(len(handles) / 2)
    fig.legend(handles=handles[:half], labels=labels[:half], loc="outside right upper")

    ax.axhline(0, color="grey", linewidth=0.5)

    if ylim is None:
        # ensure y-axis extent is symmetric around origin in steps of 100 units
        ylim = np.ceil(max(-neg.sum(axis=1).min(), pos.sum(axis=1).max()) / 100) * 100
    plt.ylim([-ylim, ylim])

    is_kt = any(s in ylabel.lower() for s in ["co2", "steel", "hvc"])
    unit = "kt/h" if is_kt else "GW"
    plt.ylabel(f"{ylabel} balance [{unit}]")

    if not long_time_frame:
        # plot frequency of snapshots on top x-axis as ticks
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(df.index)
        ax2.grid(False)
        ax2.tick_params(axis="x", length=2)
        ax2.xaxis.set_tick_params(labelbottom=False)
        ax2.set_xticklabels([])

    if resample is None:
        resample = f"native-{time}"
    fn = f"ts-balance-{ylabel.replace(' ', '_')}-{resample}"
    # plt.savefig(dir + "/" + fn + ".pdf")
    plt.savefig(dir + "/" + fn + ".png")
    plt.close()


def plot_storage(
    network,
    tech_colors,
    savepath,
    model_run="Model run",
    start_date="2019-01-01 00:00:00",
    end_date="2019-12-31 00",
    regions=["DE"],
):
    # State of charge [per unit of max] (all stores and storage units)
    # Ratio of total generation of max state of charge

    n = network

    # storage carriers
    st_carriers = [
        "battery",
        "EV battery",
        "PHS",
        "hydro",
        "H2 Store",
    ]  # "battery", "Li ion",
    # generation carriers
    carriers = [
        "battery discharger",
        "V2G",
        "PHS",
        "hydro",
        "H2",
    ]  # "battery discharger", "V2G",
    period = n.generators_t.p.index[
        (n.generators_t.p.index >= start_date) & (n.generators_t.p.index <= end_date)
    ]

    stor_res = pd.DataFrame(
        index=st_carriers,
        columns=[
            "max_charge",
            "gen_charge_ratio",
            "max_stor_cap",
            "gen_charge_cap_ratio",
            "gen_sum",
        ],
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 10), height_ratios=[1, 1.5])

    for i, c in enumerate(st_carriers):
        if c in n.storage_units.carrier.unique().tolist():
            charge = n.storage_units_t.state_of_charge.loc[
                :, n.storage_units.carrier == c
            ].sum(axis=1)[period]
            gen_sum = (
                n.storage_units_t.p_dispatch.loc[
                    period, n.storage_units.carrier == carriers[i]
                ]
                .sum()
                .sum()
            )
            gen = n.storage_units_t.p_dispatch.loc[
                period, n.storage_units.carrier == carriers[i]
            ].sum(axis=1)
            index = n.storage_units[n.storage_units.carrier == c].index
            max_stor_cap = (n.storage_units.max_hours * n.storage_units.p_nom_opt)[
                index
            ].sum()
            stor_res.loc[c, "max_charge"] = charge.max()
            stor_res.loc[c, "gen_charge_ratio"] = gen_sum / charge.max()
            stor_res.loc[c, "max_stor_cap"] = max_stor_cap
            stor_res.loc[c, "gen_charge_cap_ratio"] = gen_sum / max_stor_cap
            stor_res.loc[c, "gen_sum"] = gen_sum

        elif c in n.stores.carrier.unique().tolist():
            # state of charge (sum over different stores at same location)
            charge = n.stores_t.e.loc[:, n.stores.carrier == c].sum(axis=1)[period]
            gen_sum = (
                -n.links_t.p1.loc[period, n.links.carrier == carriers[i]].sum().sum()
            )
            max_stor_cap = n.stores.e_nom_opt[
                n.stores[n.stores.carrier == c].index
            ].sum()
            stor_res.loc[c, "max_charge"] = charge.max()
            stor_res.loc[c, "gen_charge_ratio"] = gen_sum / charge.max()
            stor_res.loc[c, "max_stor_cap"] = max_stor_cap
            stor_res.loc[c, "gen_charge_cap_ratio"] = gen_sum / max_stor_cap
            stor_res.loc[c, "gen_sum"] = gen_sum

        if c in ["battery", "EV battery", "Li ion"]:
            ax1.plot(
                charge / max_stor_cap,
                label=c,
                color=tech_colors[c],
                marker=markers[i],
                markevery=[0],
                mfc="white",
                mec="black",
            )

        else:
            ax2.plot(
                charge / max_stor_cap,
                label=c,
                color=tech_colors[c],
                marker=markers[i],
                markevery=[0],
                mfc="white",
                mec="black",
            )

    ax1.set_title(f"State of charge of short-term storage technologies ({model_run})")
    ax2.set_title(
        f"State of charge of mid- and long-term storage technologies({model_run})"
    )
    ax1.set_ylabel("State of charge [per unit of max storage capacity]")
    ax2.set_ylabel("State of charge [per unit of max storage capacity]")
    ax1.legend(loc="lower right")
    ax2.legend(loc="lower right")

    fig.tight_layout(pad=3)
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()

    return fig


def plot_price_duration_curve(
    networks,
    year_colors,
    savepath,
    years,
    carriers=["AC", "low voltage"],
    aggregate=True,
    model_run="Model run",
    regions=["DE"],
    y_lim_values=[-50, 300],
    language="english",
):
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8, 6))

    for i, n in enumerate(networks.values()):
        buses = n.buses[
            (n.buses.carrier.isin(carriers))
            & n.buses.index.str.startswith(tuple(regions))
        ].index
        lmps = pd.DataFrame(n.buses_t.marginal_price[buses])
        if aggregate:
            lmps = pd.DataFrame(lmps.mean(axis=1))
        else:
            lmps = pd.DataFrame(lmps.values.flatten())

        lmps.columns = ["lmp"]
        lmps.sort_values(by="lmp", ascending=False, inplace=True)
        lmps["percentage"] = np.arange(len(lmps)) / len(lmps) * 100
        ax.plot(lmps["percentage"], lmps["lmp"], label=years[i], color=year_colors[i])

        ax.set_ylim(y_lim_values)
        # # add corridor which contains 75 % of the generation around the median
        # ax.hlines(df["lmp"].loc[df["lmp"][df["gen_cumsum_norm"] > 0.125].index[0]], 0, 1, color=year_colors[i], ls="--", lw=1)
        # ax.hlines(df["lmp"].loc[df["lmp"][df["gen_cumsum_norm"] > 0.875].index[0]], 0, 1,  color=year_colors[i], ls="--", lw =1)

        ax.set_ylabel(
            "Strompreis [$€/MWh_{el}$]"
            if language == "german"
            else "Electricity Price [$€/MWh_{el}$]"
        )
        ax.set_xlabel(
            "Zeitanteil [%]" if language == "german" else "Fraction of time [%]"
        )
        ax.set_title(
            (
                f"Strompreisdauerlinien (Modell: {model_run})"
                if language == "german"
                else f"Electricity price duration curves {model_run}"
            ),
            fontsize=16,
        )
        ax.legend()
        ax.grid(True)

    fig.tight_layout()
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()

    return fig


def plot_price_duration_hist(
    networks,
    year_colors,
    savepath,
    years,
    carriers=["AC", "low voltage"],
    aggregate=True,
    model_run="Model run",
    regions=["DE"],
    x_lim_values=[-50, 300],
):
    fig, axes = plt.subplots(ncols=1, nrows=len(years), figsize=(8, 3 * len(years)))
    axes = axes.flatten()

    for i, n in enumerate(networks.values()):
        buses = n.buses[
            (n.buses.carrier.isin(carriers))
            & n.buses.index.str.startswith(tuple(regions))
        ].index
        lmps = pd.DataFrame(n.buses_t.marginal_price[buses])
        if aggregate:
            lmps = pd.DataFrame(lmps.mean(axis=1))
        else:
            lmps = pd.DataFrame(lmps.values.flatten())

        lmps.columns = ["lmp"]
        axes[i].hist(
            lmps,
            bins=100,
            color=year_colors[i],
            alpha=0.5,
            label=years[i],
            range=x_lim_values,
        )
        axes[i].legend()

    axes[i].set_xlabel("Electricity Price [$€/MWh_{el}$")
    plt.suptitle(f"Electricity prices ({model_run})", fontsize=16, y=0.99)
    fig.tight_layout()
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()

    return fig


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def group_pipes(df, drop_direction=False):
    """
    Group pipes which connect same buses and return overall capacity.
    """
    df = df.copy()
    if drop_direction:
        positive_order = df.bus0 < df.bus1
        df_p = df[positive_order]
        swap_buses = {"bus0": "bus1", "bus1": "bus0"}
        df_n = df[~positive_order].rename(columns=swap_buses)
        df = pd.concat([df_p, df_n])

    # there are pipes for each investment period rename to AC buses name for plotting
    df["index_orig"] = df.index
    df.index = df.apply(
        lambda x: f"H2 pipeline {x.bus0.replace(' H2', '')} -> {x.bus1.replace(' H2', '')}",
        axis=1,
    )
    return df.groupby(level=0).agg(
        {"p_nom_opt": "sum", "bus0": "first", "bus1": "first", "index_orig": "first"}
    )


def plot_h2_map(n, regions, savepath, only_de=False):
    logger.info("Plotting H2 map")
    logger.info("Assigning location")
    assign_location(n)

    h2_storage = n.stores[n.stores.carrier.isin(["H2", "H2 Store"])]
    regions["H2"] = (
        h2_storage.rename(index=h2_storage.bus.map(n.buses.location))
        .e_nom_opt.groupby(level=0)
        .sum()
        .div(1e6)
    )  # TWh
    regions["H2"] = regions["H2"].where(regions["H2"] > 0.1)

    bus_size_factor = 1e5
    linewidth_factor = 4e3
    # MW below which not drawn
    line_lower_threshold = 0

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    carriers = ["H2 Electrolysis", "H2 Fuel Cell"]

    elec = n.links[n.links.carrier.isin(carriers)].index

    bus_sizes = (
        n.links.loc[elec, "p_nom_opt"].groupby([n.links["bus0"], n.links.carrier]).sum()
        / bus_size_factor
    )

    # make a fake MultiIndex so that area is correct for legend
    bus_sizes.rename(index=lambda x: x.replace(" H2", ""), level=0, inplace=True)
    # drop all links which are not H2 pipelines
    n.links.drop(
        n.links.index[~n.links.carrier.str.contains("H2 pipeline")], inplace=True
    )

    h2_new = n.links[n.links.carrier == "H2 pipeline"]
    h2_retro = n.links[n.links.carrier == "H2 pipeline retrofitted"]
    h2_kern = n.links[n.links.carrier == "H2 pipeline (Kernnetz)"]

    # sum capacitiy for pipelines from different investment periods
    logger.info("Grouping pipes")
    h2_new = group_pipes(h2_new)

    if not h2_retro.empty:
        h2_retro = (
            group_pipes(h2_retro, drop_direction=True).reindex(h2_new.index).fillna(0)
        )

    if not h2_kern.empty:
        h2_kern = (
            group_pipes(h2_kern, drop_direction=True).reindex(h2_new.index).fillna(0)
        )

    h2_total = n.links.p_nom_opt.groupby(level=0).sum()
    link_widths_total = h2_total / linewidth_factor

    # drop all reversed pipe
    n.links.drop(n.links.index[n.links.index.str.contains("reversed")], inplace=True)
    n.links.rename(index=lambda x: x.split("-2")[0], inplace=True)
    n.links = n.links.groupby(level=0).agg(
        {
            **{
                col: "first" for col in n.links.columns if col != "p_nom_opt"
            },  # Take first value for all columns except 'p_nom_opt'
            "p_nom_opt": "sum",  # Sum values for 'p_nom_opt'
        }
    )
    link_widths_total = link_widths_total.reindex(n.links.index).fillna(0.0)
    link_widths_total[n.links.p_nom_opt < line_lower_threshold] = 0.0

    carriers_pipe = ["H2 pipeline", "H2 pipeline retrofitted", "H2 pipeline (Kernnetz)"]
    total = n.links.p_nom_opt.where(n.links.carrier.isin(carriers_pipe), other=0.0)

    retro = n.links.p_nom_opt.where(
        n.links.carrier == "H2 pipeline retrofitted", other=0.0
    )

    kern = n.links.p_nom_opt.where(
        n.links.carrier == "H2 pipeline (Kernnetz)", other=0.0
    )

    link_widths_total = total / linewidth_factor
    link_widths_total[n.links.p_nom_opt < line_lower_threshold] = 0.0

    link_widths_retro = retro / linewidth_factor
    link_widths_retro[n.links.p_nom_opt < line_lower_threshold] = 0.0

    link_widths_kern = kern / linewidth_factor
    link_widths_kern[n.links.p_nom_opt < line_lower_threshold] = 0.0

    n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    regions = regions.to_crs(proj.proj4_init)
    logger.info("Plotting map")
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": proj})

    color_h2_pipe = "#b3f3f4"
    color_retrofit = "#499a9c"
    color_kern = "#6b3161"

    bus_colors = {"H2 Electrolysis": "#ff29d9", "H2 Fuel Cell": "#805394"}

    n.plot(
        geomap=True,
        bus_sizes=bus_sizes,
        bus_colors=bus_colors,
        link_colors=color_h2_pipe,
        link_widths=link_widths_total,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n.plot(
        geomap=True,
        bus_sizes=0,
        link_colors=color_retrofit,
        link_widths=link_widths_retro,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n.plot(
        geomap=True,
        bus_sizes=0,
        link_colors=color_kern,
        link_widths=link_widths_kern,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    regions.plot(
        ax=ax,
        column="H2",
        cmap="Blues",
        linewidths=0,
        legend=True,
        vmax=6,
        vmin=0,
        legend_kwds={
            "label": "Hydrogen Storage [TWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    sizes = [50, 10]
    labels = [f"{s} GW" for s in sizes]
    sizes = [s / bus_size_factor * 1e3 for s in sizes]
    logger.info("Adding legend")
    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1),
        labelspacing=0.8,
        handletextpad=0,
        frameon=False,
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    sizes = [30, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.23, 1),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw,
    )

    colors = [bus_colors[c] for c in carriers] + [
        color_h2_pipe,
        color_retrofit,
        color_kern,
    ]
    labels = carriers + [
        "H2 pipeline (new)",
        "H2 pipeline (repurposed)",
        "H2 pipeline (Kernnetz)",
    ]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.13),
        ncol=2,
        frameon=False,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)

    ax.set_facecolor("white")

    fig.savefig(savepath, bbox_inches="tight")
    plt.close()


def plot_h2_map_de(n, regions, tech_colors, savepath, specify_buses=None):
    assign_location(n)

    h2_storage = n.stores[n.stores.carrier.isin(["H2", "H2 Store"])]
    regions["H2"] = (
        h2_storage.rename(index=h2_storage.bus.map(n.buses.location))
        .e_nom_opt.groupby(level=0)
        .sum()
        .div(1e6)
    )  # TWh
    regions["H2"] = regions["H2"].where(regions["H2"] > 0.1)

    linewidth_factor = 4e3
    # MW below which not drawn
    line_lower_threshold = 0

    # buses and size
    if specify_buses is None:
        bus_size_factor = 1e5
        carriers = ["H2 Electrolysis", "H2 Fuel Cell"]
        elec = n.links[
            (n.links.carrier.isin(carriers)) & (n.links.index.str.contains("DE"))
        ].index
        bus_sizes = (
            n.links.loc[elec, "p_nom_opt"]
            .groupby([n.links["bus0"], n.links.carrier])
            .sum()
            / bus_size_factor
        )
    if specify_buses == "production":
        bus_size_factor = 2e8
        h2_producers = n.links.index[
            n.links.index.str.startswith("DE")
            & (n.links.bus1.map(n.buses.carrier) == "H2")
        ]
        carriers = h2_producers.map(n.links.carrier).unique().tolist()
        production = -n.links_t.p1[h2_producers].multiply(
            n.snapshot_weightings.generators, axis=0
        )
        bus_sizes = (
            production.sum()
            .groupby(
                [
                    production.sum().index.map(n.links.bus1),
                    production.sum().index.map(n.links.carrier),
                ]
            )
            .sum()
            / bus_size_factor
        )

    if specify_buses == "consumption":
        bus_size_factor = 2e8
        # links
        h2_consumers_links = n.links.index[
            n.links.index.str.startswith("DE")
            & (n.links.bus0.map(n.buses.carrier) == "H2")
        ]
        consumption_links = n.links_t.p0[h2_consumers_links].multiply(
            n.snapshot_weightings.generators, axis=0
        )
        bus_sizes_links = (
            consumption_links.sum()
            .groupby(
                [
                    consumption_links.sum().index.map(n.links.bus0),
                    consumption_links.sum().index.map(n.links.carrier),
                ]
            )
            .sum()
            / bus_size_factor
        )
        # loads
        h2_consumers_loads = n.loads.index[
            n.loads.bus.str.startswith("DE")
            & (n.loads.bus.map(n.buses.carrier) == "H2")
        ]
        consumption_loads = n.loads_t.p[h2_consumers_loads].multiply(
            n.snapshot_weightings.generators, axis=0
        )
        bus_sizes_loads = (
            consumption_loads.sum()
            .groupby(
                [
                    consumption_loads.sum().index.map(n.loads.bus),
                    consumption_loads.sum().index.map(n.loads.carrier),
                ]
            )
            .sum()
            / bus_size_factor
        )

        bus_sizes = pd.concat([bus_sizes_links, bus_sizes_loads])

        def rename_carriers(carrier):
            if "H2" in carrier and "OCGT" in carrier:
                return "H2 OCGT"
            elif "H2" in carrier and "CHP" in carrier:
                return "H2 CHP"
            else:
                return carrier

        bus_sizes = bus_sizes.rename(index=lambda x: rename_carriers(x), level=1)
        bus_sizes = bus_sizes.groupby(level=[0, 1]).sum()
        tech_colors["H2 CHP"] = "darkorange"
        # only select 4 most contributing carriers and summarise rest as other
        others = (
            (bus_sizes.groupby(level=[1]).sum() / bus_sizes.sum())
            .sort_values(ascending=False)[5:]
            .index.tolist()
        )
        replacement_dict = {value: "other" for value in others}
        bus_sizes = bus_sizes.rename(index=replacement_dict, level=1)
        bus_sizes = bus_sizes.groupby(level=[0, 1]).sum()
        carriers = bus_sizes.index.get_level_values(1).unique().tolist()

    # make a fake MultiIndex so that area is correct for legend
    bus_sizes.rename(index=lambda x: x.replace(" H2", ""), level=0, inplace=True)

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)
    # drop all links which are not H2 pipelines or not in germany
    n.links.drop(
        n.links.index[
            ~(
                (n.links["carrier"].str.contains("H2 pipeline"))
                & (n.links.index.str.contains("DE"))
            )
        ],
        inplace=True,
    )

    h2_new = n.links[n.links.carrier == "H2 pipeline"]
    h2_retro = n.links[n.links.carrier == "H2 pipeline retrofitted"]
    h2_kern = n.links[n.links.carrier == "H2 pipeline (Kernnetz)"]

    # sum capacitiy for pipelines from different investment periods
    h2_new = group_pipes(h2_new)

    if not h2_retro.empty:
        h2_retro = (
            group_pipes(h2_retro, drop_direction=True).reindex(h2_new.index).fillna(0)
        )

    if not h2_kern.empty:
        h2_kern = (
            group_pipes(h2_kern, drop_direction=True).reindex(h2_new.index).fillna(0)
        )

    h2_total = n.links.p_nom_opt.groupby(level=0).sum()
    link_widths_total = h2_total / linewidth_factor

    # drop all reversed pipe
    n.links.drop(n.links.index[n.links.index.str.contains("reversed")], inplace=True)
    n.links.rename(index=lambda x: x.split("-2")[0], inplace=True)
    n.links = n.links.groupby(level=0).agg(
        {
            **{
                col: "first" for col in n.links.columns if col != "p_nom_opt"
            },  # Take first value for all columns except 'p_nom_opt'
            "p_nom_opt": "sum",  # Sum values for 'p_nom_opt'
        }
    )
    link_widths_total = link_widths_total.reindex(n.links.index).fillna(0.0)
    link_widths_total[n.links.p_nom_opt < line_lower_threshold] = 0.0

    carriers_pipe = ["H2 pipeline", "H2 pipeline retrofitted", "H2 pipeline (Kernnetz)"]
    total = n.links.p_nom_opt.where(n.links.carrier.isin(carriers_pipe), other=0.0)

    retro = n.links.p_nom_opt.where(
        n.links.carrier == "H2 pipeline retrofitted", other=0.0
    )

    kern = n.links.p_nom_opt.where(
        n.links.carrier == "H2 pipeline (Kernnetz)", other=0.0
    )

    link_widths_total = total / linewidth_factor
    link_widths_total[n.links.p_nom_opt < line_lower_threshold] = 0.0

    link_widths_retro = retro / linewidth_factor
    link_widths_retro[n.links.p_nom_opt < line_lower_threshold] = 0.0

    link_widths_kern = kern / linewidth_factor
    link_widths_kern[n.links.p_nom_opt < line_lower_threshold] = 0.0

    n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    regions = regions.to_crs(proj.proj4_init)

    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": proj})

    color_h2_pipe = "#b3f3f4"
    color_retrofit = "#499a9c"
    color_kern = "#6b3161"

    n.plot(
        geomap=True,
        bus_sizes=bus_sizes,
        bus_colors=tech_colors,
        link_colors=color_h2_pipe,
        link_widths=link_widths_total,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n.plot(
        geomap=True,
        bus_sizes=0,
        link_colors=color_retrofit,
        link_widths=link_widths_retro,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n.plot(
        geomap=True,
        bus_sizes=0,
        link_colors=color_kern,
        link_widths=link_widths_kern,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    regions.plot(
        ax=ax,
        column="H2",
        cmap="Blues",
        linewidths=0,
        legend=True,
        vmax=6,
        vmin=0,
        legend_kwds={
            "label": "Hydrogen Storage [TWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    # Set geographic extent for Germany
    ax.set_extent([5.5, 15.5, 48.0, 56], crs=ccrs.PlateCarree())  # Germany bounds

    if specify_buses is None:
        sizes = [5, 1]
        labels = [f"{s} GW" for s in sizes]
        sizes = [s / bus_size_factor * 1e3 for s in sizes]
        n_cols = 2
    elif specify_buses == "production":
        sizes = [10, 1]
        labels = [f"{s} TWh" for s in sizes]
        sizes = [s / bus_size_factor * 1e6 for s in sizes]
        n_cols = 2
    elif specify_buses == "consumption":
        sizes = [10, 1]
        labels = [f"{s} TWh" for s in sizes]
        sizes = [s / bus_size_factor * 1e6 for s in sizes]
        n_cols = 3

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1),
        labelspacing=0.8,
        handletextpad=0,
        frameon=False,
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    sizes = [30, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.23, 1),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw,
    )

    colors = [tech_colors[c] for c in carriers] + [
        color_h2_pipe,
        color_retrofit,
        color_kern,
    ]
    labels = carriers + [
        "H2 pipeline (new)",
        "H2 pipeline (repurposed)",
        "H2 pipeline (Kernnetz)",
    ]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.13),
        ncol=n_cols,
        frameon=False,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)

    ax.set_facecolor("white")
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()


### electricity transmission


def plot_elec_map_de(
    network,
    base_network,
    tech_colors,
    regions_de,
    savepath,
    expansion_case="total-expansion",
):
    m = network.copy()
    m.mremove("Bus", m.buses[m.buses.x == 0].index)
    m.buses.drop(m.buses.index[m.buses.carrier != "AC"], inplace=True)

    m_base = base_network.copy()

    # storage as cmap on map
    battery_storage = m.stores[m.stores.carrier.isin(["battery"])]
    regions_de["battery"] = (
        battery_storage.rename(
            index=battery_storage.bus.str.replace(" battery", "").map(m.buses.location)
        )
        .e_nom_opt.groupby(level=0)
        .sum()
        .div(1e3)
    )  # GWh
    regions_de["battery"] = regions_de["battery"].where(regions_de["battery"] > 0.1)

    # buses
    bus_size_factor = 0.5e6
    carriers = ["onwind", "offwind-ac", "offwind-dc", "solar", "solar-hsat"]
    elec = m.generators[
        (m.generators.carrier.isin(carriers)) & (m.generators.bus.str.contains("DE"))
    ].index
    bus_sizes = (
        m.generators.loc[elec, "p_nom_opt"]
        .groupby([m.generators.bus, m.generators.carrier])
        .sum()
        / bus_size_factor
    )
    replacement_dict = {
        "onwind": "Onshore Wind",
        "offwind-ac": "Offshore Wind",
        "offwind-dc": "Offshore Wind",
        "solar": "Solar",
        "solar-hsat": "Solar",
    }
    bus_sizes = bus_sizes.rename(index=replacement_dict, level=1)
    bus_sizes = bus_sizes.groupby(level=[0, 1]).sum()
    carriers = bus_sizes.index.get_level_values(1).unique().tolist()

    # lines
    linew_factor = 1e3
    linkw_factor = 0.5e3

    # line widths
    startnetz_i = m.lines[m.lines.build_year != 0].index
    total_exp_linew = m.lines.s_nom_opt - m_base.lines.s_nom_min
    total_exp_linew[startnetz_i] = m.lines.s_nom_opt[startnetz_i]
    total_exp_noStart_linew = total_exp_linew.copy()
    total_exp_noStart_linew.loc[startnetz_i] = 0
    startnetz_linew = m.lines.s_nom_opt.loc[startnetz_i]

    # link widths
    tprojs = m.links.loc[
        (m.links.index.str.startswith("DC") | m.links.index.str.startswith("TYNDP"))
        & ~m.links.reversed
    ].index
    links_i = m.links.index[m.links.carrier == "DC"]
    total_exp_linkw = (m.links.p_nom_opt - m_base.links.p_nom).loc[links_i]
    total_exp_linkw[tprojs] = m.links.p_nom_opt[tprojs]
    total_exp_noStart_linkw = total_exp_linkw.copy()
    total_exp_noStart_linkw.loc[tprojs] = 0
    startnetz_linkw = m.links.p_nom_opt[tprojs]

    if expansion_case == "total-expansion":
        line_widths = total_exp_linew / linew_factor
        link_widths = total_exp_linkw / linkw_factor
    elif expansion_case == "startnetz":
        line_widths = startnetz_linew / linew_factor
        link_widths = startnetz_linkw / linkw_factor
    elif expansion_case == "pypsa":
        line_widths = total_exp_noStart_linew / linew_factor
        link_widths = total_exp_noStart_linkw / linkw_factor
    else:
        line_widths = None
        link_widths = None

    regions_de = regions_de.to_crs(proj.proj4_init)
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": proj})

    m.plot(
        ax=ax,
        margin=0.06,
        bus_sizes=bus_sizes,
        bus_colors=tech_colors,
        line_widths=line_widths,
        line_colors=tech_colors["AC"],
        link_widths=link_widths,
        link_colors=tech_colors["DC"],
    )

    regions_de.plot(
        ax=ax,
        column="battery",
        cmap="Oranges",
        linewidths=0,
        legend=True,
        legend_kwds={
            "label": "Battery Storage [GWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    # Set geographic extent for Germany
    ax.set_extent([5.5, 15.5, 47, 56], crs=ccrs.PlateCarree())

    sizes = [10, 5]
    labels = [f"{s} GW" for s in sizes]
    sizes = [s / bus_size_factor * 1e3 for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1),
        labelspacing=0.8,
        handletextpad=0,
        frameon=True,
        facecolor="white",
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=m.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    # AC
    sizes_ac = [10, 5]
    labels_ac = [f"HVAC ({s} GW)" for s in sizes_ac]
    scale = 1e3 / linew_factor
    sizes_ac = [s * scale for s in sizes_ac]

    # DC
    sizes_dc = [5, 2]
    labels_dc = [f"HVDC ({s} GW)" for s in sizes_dc]
    scale = 1e3 / linkw_factor
    sizes_dc = [s * scale for s in sizes_dc]

    sizes = sizes_ac + sizes_dc
    labels = labels_ac + labels_dc
    colors = [tech_colors["AC"]] * len(sizes_ac) + [tech_colors["DC"]] * len(sizes_dc)

    legend_kw = dict(
        loc=[0.2, 0.9],
        frameon=True,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=13,
        ncol=2,
        facecolor="white",
    )

    add_legend_lines(ax, sizes, labels, colors=colors, legend_kw=legend_kw)

    colors = [tech_colors[c] for c in carriers]
    labels = carriers
    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 0.9),
        ncol=2,
        frameon=True,
        facecolor="white",
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_ariadne_report",
            simpl="",
            clusters=49,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Bal_v4",
        )

    configure_logging(snakemake)

    ### Modify postnetworks (this might be moved to a separate script)

    # Load costs (needed for modification)
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

    # Load data
    _networks = [pypsa.Network(fn) for fn in snakemake.input.networks]
    modelyears = [fn[-7:-3] for fn in snakemake.input.networks]

    # Hack the transmission projects
    networks = [
        process_postnetworks(n.copy(), _networks[0], int(my), snakemake, costs)
        for n, my in zip(_networks, modelyears)
    ]
    del _networks
    ###

    # ensure output directory exist
    for dir in snakemake.output[2:]:
        if not os.path.exists(dir):
            os.makedirs(dir)

    # configs
    config = snakemake.config
    planning_horizons = snakemake.params.planning_horizons
    tech_colors = snakemake.params.plotting["tech_colors"]

    # update tech_colors
    colors_update = (
        networks[0].carriers.color.rename(networks[0].carriers.nice_name).to_dict()
    )
    colors_update = {k: v for k, v in colors_update.items() if v != ""}
    tech_colors.update(colors_update)

    # define possible renaming and grouping of carriers
    c_g = [solar, electricity_load, electricity_imports]
    c_n = ["Solar", "Electricity load", "Electricity trade"]

    # add colors for renaming and condensed groups
    for old_name, new_name in carrier_renaming.items():
        if old_name in tech_colors:
            tech_colors[new_name] = tech_colors[old_name]
    for name in c1_groups_name:
        tech_colors[name] = tech_colors[f"urban central {name}"]
    for name in c1_groups_name:
        tech_colors[name] = tech_colors[f"urban central {name}"]

    # carrier names manual
    tech_colors["urban central oil CHP"] = tech_colors["oil"]
    tech_colors["Solar"] = tech_colors["solar"]
    tech_colors["Electricity load"] = tech_colors["electricity"]
    tech_colors["Electricity trade"] = tech_colors["AC"]
    tech_colors["Offshore Wind"] = tech_colors["offwind-ac"]
    tech_colors["urban decentral heat"] = tech_colors["urban central heat"]
    tech_colors["urban decentral biomass boiler"] = tech_colors["biomass boiler"]
    tech_colors["rural biomass boiler"] = tech_colors["biomass boiler"]
    tech_colors["urban decentral oil boiler"] = tech_colors["oil boiler"]
    tech_colors["rural oil boiler"] = tech_colors["oil boiler"]
    tech_colors["rural ground heat pump"] = tech_colors["ground heat pump"]

    ### plotting
    for year in planning_horizons:
        network = networks[planning_horizons.index(year)].copy()
        ct = "DE"
        buses = network.buses.index[(network.buses.index.str[:2] == ct)].drop("DE")
        balance = (
            network.statistics.energy_balance(
                aggregate_time=False,
                nice_names=False,
                groupby=network.statistics.groupers.get_bus_and_carrier_and_bus_carrier,
            )
            .loc[:, buses, :, :]
            .droplevel("bus")
        )

        # electricity supply and demand
        logger.info("Plotting electricity supply and demand for year %s", year)
        plot_nodal_balance(
            network=network,
            nodal_balance=balance,
            tech_colors=tech_colors,
            start_date="2019-01-01 00:00:00",
            end_date="2019-12-31 00:00:00",
            savepath=f"{snakemake.output.elec_balances}/elec-all-year-DE-{year}.png",
            model_run=snakemake.wildcards.run,
            resample="D",
            plot_lmps=False,
            plot_loads=False,
            nice_names=True,
            threshold=1e2,  # in GWh as sum over period
            condense_groups=c_g,
            condense_names=c_n,
            title="Strombilanz",
        )

        plot_nodal_balance(
            network=network,
            nodal_balance=balance,
            tech_colors=tech_colors,
            start_date="2019-01-01 00:00:00",
            end_date="2019-01-31 00:00:00",
            savepath=f"{snakemake.output.elec_balances}/elec-Jan-DE-{year}.png",
            model_run=snakemake.wildcards.run,
            nice_names=True,
            threshold=1e2,
            condense_groups=[electricity_load, electricity_imports],
            condense_names=["Electricity load", "Electricity trade"],
            title="Strombilanz",
        )

        plot_nodal_balance(
            network=network,
            nodal_balance=balance,
            tech_colors=tech_colors,
            start_date="2019-05-01 00:00:00",
            end_date="2019-05-31 00:00:00",
            savepath=f"{snakemake.output.elec_balances}/elec-May-DE-{year}.png",
            model_run=snakemake.wildcards.run,
            nice_names=True,
            threshold=1e2,
            condense_groups=[electricity_load, electricity_imports],
            condense_names=["Electricity load", "Electricity trade"],
            title="Strombilanz",
        )

        # heat supply and demand
        logger.info("Plotting heat supply and demand")
        for carriers in ["urban central heat", "urban decentral heat", "rural heat"]:
            plot_nodal_balance(
                network=network,
                nodal_balance=balance,
                tech_colors=tech_colors,
                start_date="2019-01-01 00:00:00",
                end_date="2019-12-31 00:00:00",
                savepath=f"{snakemake.output.heat_balances}/heat-all-year-DE-{carriers}-{year}.png",
                model_run=snakemake.wildcards.run,
                resample="D",
                plot_lmps=False,
                plot_loads=False,
                nice_names=True,
                threshold=1e1,  # in GWh as sum over period
                condense_groups=c_g,
                condense_names=c_n,
                carriers=[carriers],
                ylabel="Heat [GW]",
                title=f"{carriers} balance",
            )

            plot_nodal_balance(
                network=network,
                nodal_balance=balance,
                tech_colors=tech_colors,
                start_date="2019-01-01 00:00:00",
                end_date="2019-01-31 00:00:00",
                savepath=f"{snakemake.output.heat_balances}/heat-Jan-DE-{carriers}-{year}.png",
                model_run=snakemake.wildcards.run,
                plot_lmps=False,
                plot_loads=False,
                nice_names=True,
                threshold=1e1,
                carriers=[carriers],
                ylabel="Heat [GW]",
                title=f"{carriers} balance",
            )

            plot_nodal_balance(
                network=network,
                nodal_balance=balance,
                tech_colors=tech_colors,
                start_date="2019-05-01 00:00:00",
                end_date="2019-05-31 00:00:00",
                savepath=f"{snakemake.output.heat_balances}/heat-May-DE-{carriers}-{year}.png",
                model_run=snakemake.wildcards.run,
                plot_lmps=False,
                plot_loads=False,
                nice_names=True,
                threshold=1e1,
                carriers=[carriers],
                ylabel="Heat [GW]",
                title=f"{carriers} balance",
            )

        # storage
        logger.info("Plotting storage")
        plot_storage(
            network=network,
            tech_colors=tech_colors,
            start_date="2019-01-01 00:00:00",
            end_date="2019-12-31 00:00:00",
            savepath=f"{snakemake.output.results}/storage-DE-{year}.png",
            model_run=snakemake.wildcards.run,
        )

    ## price duration
    logger.info("Plotting price duration curve")
    networks_dict = {int(my): n for n, my in zip(networks, modelyears)}
    plot_price_duration_curve(
        networks=networks_dict,
        year_colors=year_colors,
        savepath=snakemake.output.elec_price_duration_curve,
        model_run=snakemake.wildcards.run,
        years=planning_horizons,
        language="german",
    )

    plot_price_duration_hist(
        networks=networks_dict,
        year_colors=year_colors,
        savepath=snakemake.output.elec_price_duration_hist,
        model_run=snakemake.wildcards.run,
        years=planning_horizons,
    )

    ## hydrogen transmission
    logger.info("Plotting hydrogen transmission")
    map_opts = snakemake.params.plotting["map"]
    snakemake.params.plotting["projection"] = {"name": "EqualEarth"}
    proj = load_projection(snakemake.params.plotting)
    regions = gpd.read_file(snakemake.input.regions_onshore_clustered).set_index("name")

    for year in planning_horizons:
        network = networks[planning_horizons.index(year)].copy()
        logger.info(f"Plotting hydrogen transmission for {year}")
        plot_h2_map(
            network,
            regions,
            savepath=f"{snakemake.output.h2_transmission}/h2_transmission_all-regions_{year}.png",
        )

        regions_de = regions[regions.index.str.startswith("DE")]
        logger.info(f"Plotting hydrogen transmission for {year} in DE")
        del network
        for sb in ["production", "consumption"]:
            network = networks[planning_horizons.index(year)].copy()
            plot_h2_map_de(
                network,
                regions_de,
                tech_colors=tech_colors,
                specify_buses=sb,
                savepath=f"{snakemake.output.h2_transmission}/h2_transmission_DE_{sb}_{year}.png",
            )
            del network

    ## electricity transmission
    logger.info("Plotting electricity transmission")
    for year in planning_horizons:
        scenarios = ["total-expansion", "startnetz", "pypsa"]
        for s in scenarios:
            plot_elec_map_de(
                networks[planning_horizons.index(year)],
                networks[planning_horizons.index(2020)],
                tech_colors,
                regions_de,
                savepath=f"{snakemake.output.elec_transmission}/elec-transmission-DE-{s}-{year}.png",
                expansion_case=s,
            )

    ## nodal balances general (might not be very robust)
    logger.info("Plotting nodal balances")
    plt.style.use(["bmh", snakemake.input.rc])

    year = 2045
    network = networks[planning_horizons.index(year)].copy()
    n = network

    months = pd.date_range(freq="ME", **snakemake.config["snapshots"]).map(
        lambda x: x.strftime("%Y-%m")
    )

    balance = n.statistics.energy_balance(aggregate_time=False)

    # only DE
    ct = "DE"
    buses = n.buses.index[(n.buses.index.str[:2] == ct)].drop("DE")
    balance = (
        n.statistics.energy_balance(
            aggregate_time=False,
            groupby=n.statistics.groupers.get_bus_and_carrier_and_bus_carrier,
        )
        .loc[:, buses, :, :]
        .droplevel("bus")
    )

    n.carriers.color.update(snakemake.config["plotting"]["tech_colors"])
    n.carriers.color.update(tech_colors)
    colors = n.carriers.color.rename(n.carriers.nice_name)
    # replace empty values TODO add empty values with colors to plotting config
    colors[colors.values == ""] = "lightgrey"

    # wrap in function for multiprocessing
    def process_group(group, carriers, balance, months, colors):
        if not isinstance(carriers, list):
            carriers = [carriers]

        mask = balance.index.get_level_values("bus_carrier").isin(carriers)
        df = balance[mask].groupby("carrier").sum().div(1e3).T

        # daily resolution for each carrier
        plot_energy_balance_timeseries(
            df,
            resample="D",
            ylabel=group,
            colors=colors,
            threshold=THRESHOLD,
            dir=dir,
        )

        # monthly resolution for each carrier
        plot_energy_balance_timeseries(
            df,
            resample="M",
            ylabel=group,
            colors=colors,
            threshold=THRESHOLD,
            dir=dir,
        )

        # native resolution for each month and carrier
        for month in months:
            plot_energy_balance_timeseries(
                df,
                time=month,
                ylabel=group,
                colors=colors,
                threshold=THRESHOLD,
                dir=dir,
            )

    args = [
        (group, carriers, balance, months, colors)
        for group, carriers in CARRIER_GROUPS.items()
    ]
    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_group, args)
