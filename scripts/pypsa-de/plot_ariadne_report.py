# -*- coding: utf-8 -*-
import logging
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import locale
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
from export_ariadne_variables import get_discretized_value, process_postnetworks
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

from scripts._helpers import configure_logging, mock_snakemake, set_scenario_config
from scripts.plot_power_network import load_projection
from scripts.plot_summary import preferred_order, rename_techs
from scripts.prepare_sector_network import prepare_costs

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

backup_techs = {
    "Gas": [
        "OCGT",
        "CCGT",
        "biogas",
        "urban central gas CHP",
        "urban central gas CHP CC",
    ],
    # "Öl": ["oil", "urban central oil CHP"],
    "Kohle": ["lignite", "coal", "urban central coal CHP", "urban central lignite CHP"],
    "Wasserkraft": ["PHS", "hydro"],
    "Batterie": ["battery discharger", "home battery discharger"],
    "Wasserstoff": [
        "H2 OCGT",
        "H2 retrofit OCGT",
        "urban central H2 CHP",
        "urban central H2 retrofit CHP",
    ],
    # "Müllverbrennung": ["waste CHP", "waste CHP CC"],
    "Biomasse": [
        "solid biomass",
        "urban central solid biomass CHP",
        "urban central solid biomass CHP CC",
    ],
}

vre_gens = [
    "onwind",
    "offwind-ac",
    "offwind-dc",
    "solar",
    "solar-hsat",
    "solar rooftop",
    "ror",
]

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

carriers_in_german = {
    "DC": "Gleichstrom",
    "AC": "Wechselstrom",
    "hydro": "Wasserkraft (Reservoir & Damm)",
    "offwind-ac": "Offshore-Wind (AC)",
    "offwind-dc": "Offshore-Wind (DC)",
    "solar": "Solar",
    "solar-hsat": "Solar (HSAT)",
    "onwind": "Onshore-Wind",
    "PHS": "Pumpspeicherkraftwerk",
    "ror": "Laufwasserkraft",
    "": "",
    "lignite": "Braunkohle",
    "coal": "Steinkohle",
    "oil": "Öl",
    "uranium": "Uran",
    "none": "keine",
    "co2": "CO2",
    "co2 stored": "CO2 gespeichert",
    "co2 sequestered": "CO2 sequestriert",
    "gas": "Gas",
    "H2": "Wasserstoffspeicher",
    "battery": "Batteriespeicher",
    "EV battery": "Elektrofahrzeug-Batterie",
    "urban central heat": "Zentrale städtische Heizung",
    "urban central water tanks": "Zentrale städtische Wassertanks",
    "urban central solar thermal": "Zentrale städtische Solarthermie",
    "biogas": "Biogas",
    "solid biomass": "Biomasse",
    "methanol": "Methanol",
    "home battery": "Hausbatterie",
    "rural heat": "Ländliche Wärme",
    "rural solar thermal": "Ländliche Solarthermie",
    "rural water tanks": "Ländliche Wassertanks",
    "urban decentral heat": "Dezentrale städtische Heizung",
    "urban decentral solar thermal": "Dezentrale städtische Solarthermie",
    "urban decentral water tanks": "Dezentrale städtische Wassertanks",
    "oil primary": "Primäröl",
    "solar rooftop": "Solar-Dach",
    "urban central heat vent": "Zentrale städtische Wärmeentlüftung",
    "gas primary": "Primärgas",
    "solid biomass for industry": "Biomasse (Industrie)",
    "shipping oil": "Schiffsöl",
    "agriculture machinery oil": "Landwirtschaftsmaschinenöl",
    "naphtha for industry": "Naphtha für die Industrie",
    "land transport oil": " Öl (Transport)",
    "kerosene for aviation": "Kerosin für die Luftfahrt",
    "non-sequestered HVC": "Nicht-sequestriertes HVC",
    "shipping methanol": "Methanol für Schifffahrt",
    "gas for industry": "Gas (Industrie)",
    "low voltage": "Niederspannung",
    "industry methanol": "Methanol (Industrie)",
    "coal for industry": "Kohle (Industrie)",
    "process emissions": "Prozessemissionen",
    "industry electricity": "Industrieelektrizität",
    "low-temperature heat for industry": "Niedertemperaturwärme (Industrie)",
    "agriculture electricity": "Landwirtschaftliche Elektrizität",
    "electricity": "Elektrizität",
    "land transport EV": "Elektrofahrzeuge (Transport)",
    "agriculture heat": "Landwirtschaftliche Wärme",
    "urban central air heat pump": "Zentrale städtische Luftwärmepumpe",
    "electricity distribution grid": "Stromverteilungsnetz",
    "battery charger": "Batterie Laden",
    "waste CHP": "Müll-KWK",
    "urban central water tanks discharger": "Entladung zentraler städtischer Wassertanks",
    "rural water tanks discharger": "Entladung ländlicher Wassertanks",
    "urban decentral biomass boiler": "Dezentrale städtische Biomassekessel",
    "BEV charger": "E-Fahrzeug Laden",
    "rural ground heat pump": "Erdwärmepumpe",
    "Fischer-Tropsch": "Fischer-Tropsch",
    "urban decentral water tanks discharger": "Entladung dezentraler städtischer Wassertanks",
    "urban central gas CHP": "Gas-KWK",
    "Sabatier": "Sabatier-Prozess",
    "gas compressing": "Gasverdichtung",
    "home battery charger": "Hausbatterie Laden",
    "battery discharger": "Batterie Entladung",
    "H2 pipeline retrofitted": "Nachgerüstete Wasserstoffpipeline",
    "rural resistive heater": "Ländlicher Widerstandsheizer",
    "urban central water tanks charger": "Ladung zentraler städtischer Wassertanks",
    "urban central solid biomass CHP CC": "Biomasse-KWK mit CO2-Abscheidung",
    "rural air heat pump": "Ländliche Luftwärmepumpe",
    "DAC": "Direkte CO2-Abscheidung",
    "urban decentral air heat pump": "Dezentrale städtische Luftwärmepumpe",
    "waste CHP CC": "Müll-KWK mit CO2-Abscheidung",
    "biomass to liquid CC": "Biomasse zu Flüssigkeit mit CO2-Abscheidung",
    "HVC to air": "HVC in die Luft",
    "methanolisation": "Methanolisation",
    "gas for industry CC": "Gas mit CO2-Abscheidung (Industrie)",
    "H2 pipeline": "Wasserstoffpipeline",
    "solid biomass for industry CC": "Biomasse mit CO2-Abscheidung (Industrie)",
    "urban decentral gas boiler": "Dezentrale städtische Gaskessel",
    "urban central gas CHP CC": "Gas-KWK mit CO2-Abscheidung",
    "rural gas boiler": "Ländlicher Gaskessel",
    "process emissions CC": "Prozessemissionen mit CO2-Abscheidung",
    "urban decentral water tanks charger": "Ladung dezentraler städtischer Wassertanks",
    "biogas to gas CC": "Biogas zu Gas mit CO2-Abscheidung",
    "urban decentral resistive heater": "Dezentrale städtische Widerstandsheizer",
    "biogas to gas": "Biogas zu Gas",
    "OCGT": "Gas (OCGT)",
    "urban central solid biomass CHP": "Biomasse-KWK",
    "urban central gas boiler": "Zentraler städtischer Gaskessel",
    "urban central resistive heater": "Zentraler städtischer Widerstandsheizer",
    "oil refining": "Ölraffinierung",
    "rural biomass boiler": "Ländlicher Biomassekessel",
    "SMR": "Dampfreformierung",
    "biomass to liquid": "Biomasse zu Flüssigkeit",
    "rural water tanks charger": "Ladung ländlicher Wassertanks",
    "home battery discharger": "Hausbatterie Entladung",
    "H2 Store": "Wasserstoffspeicher",
    "renewable oil": "Erneuerbares Öl",
    "renewable gas": "Erneuerbares Gas",
    "CCGT": "Gas (CCGT)",
    "nuclear": "Kernenergie",
    "gas CHP": "Gas-KWK",
    "gas CHP CC": "Gas KWK mit CO2-Abscheidung",
    "urban central coal CHP": "Steinkohle-KWK",
    "Electricity trade": "Stromhandel",
    "urban central gas CHP": "Gas-KWK",
    "urban central biomass CHP": "Biomasse-KWK",
    "biomass CHP": "Biomasse-KWK",
    "air heat pump": "Luftwärmepumpe",
    "Electricity load": "Stromlast",
    "methanolisation": "Methanolisierung",
    "resistive heater": "Widerstandsheizung",
    "gas boiler": "Gaskessel",
    "H2 Electrolysis": "Elektrolyse",
    "H2 Fuel Cell": "Brennstoffzelle (Strom)",
    "H2 for industry": "H2 für Industrie",
    "H2 OCGT": "Wasserstoff (OCGT)",
    "H2 CHP": "H2 KWK",
    "land transport fuel cell": "Brennstoffzelle (Verkehr)",
    "other": "Sonstige",
    "SMR CC": "Dampfreformierung mit CCS",
    "H2 pipeline (new)": "H2 Pipeline (Neubau)",
    "H2 pipeline (repurposed)": "H2 Pipeline (Umstellung)",
    "H2 pipeline (Kernnetz)": "H2 Pipeline (Kernnetz)",
    "heat pump": "Wärmepumpe",
    "urban central H2 CHP": "Wasserstoff-KWK",
    "urban central H2 retrofit CHP": "Wasserstoff-KWK (Umrüstung)",
    "urban central oil CHP": "Öl-KWK",
    "urban central lignite CHP": "Braunkohle-KWK",
    "H2 retrofit OCGT": "Wasserstoff (OCGT;Umrüstung)",
}


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


def plot_nodal_elec_balance(
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

    if resample == "D" and network.snapshots.size < 365:
        # code is not working at low resolution!
        logger.error(
            "Temporal resolution does not allow for daily resampling! Please use hihger resolution results or change the 'resample' flag."
        )
        return

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

    ## summaris low contributing carriers acccording to their sum over the period (threshold in GWh)
    techs_below_threshold = nb.columns[nb.abs().sum() < threshold].tolist()
    if techs_below_threshold:
        other = {tech: "other" for tech in techs_below_threshold}
        rename.update(other)
        tech_colors["other"] = "grey"

    if rename:
        nb = nb.T.groupby(nb.columns.map(lambda a: rename.get(a, a))).sum().T

    if resample is not None:
        nb = nb.resample(resample).mean()

    # plot positive values
    preferred_order_pos = [
        "solar",
        "onwind",
        "offwind-ac",
        "offwind-dc",
        "ror",
        "H2 OCGT",
        "H2 CHP",
        "urban central H2 retrofit CHP",
        "CCGT",
        "gas CHP",
        "waste CHP CC",
        "urban central biomass CHP",
        "hydro",
        "PHS",
        "battery discharger",
        "Stromimport",
        "other",
    ]
    preferred_order_neg = [
        "Electricity load",
        "electricity distribution grid",
        "BEV charger",
        "air heat pump",
        "rural ground heat pump",
        "resistive heater",
        "battery charger",
        "Stromexport",
        "H2 Electrolysis",
        "methanolisation",
        "Sonstige",
    ]
    pos_c = {
        "solar": "#f9d002",
        "onwind": "#235ebc",
        "offwind-ac": "#6895dd",
        "offwind-dc": "#74c6f2",
        "ror": "#81a3de",
        "H2 OCGT": "#ff0000",
        "H2 CHP": "#bf3737",
        "urban central H2 retrofit CHP": "#ff8282",
        "CCGT": "#edc566",
        "gas CHP": "#e67c12",
        "waste CHP CC": "#a85522",
        "urban central biomass CHP": "#9d9042",
        "hydro": "#5379ad",
        "PHS": "#6999db",
        "battery discharger": "#76e388",
        "Stromimport": "#97ad8c",
        "other": "#8f9c9a",
    }
    neg_c = {
        "Electricity load": "#110d63",
        "electricity distribution grid": "#baf238",
        "air heat pump": "#e3d3ff",
        "rural ground heat pump": "#9f6df7",
        "resistive heater": "#493173",
        "BEV charger": "#81a3de",
        "battery charger": "#76e388",
        "Stromexport": "#97ad8c",
        "H2 Electrolysis": "#ff8282",
        "methanolisation": "#872f2f",
        "Sonstige": "#8f9c9a",
    }

    df = nb
    # split into df with positive and negative values
    df_neg, df_pos = df.clip(upper=0), df.clip(lower=0)

    df_pos["solar"] = (
        df_pos.get("solar", 0)
        + df_pos.get("Solar", 0)
        + df_pos.get("solar-hsat", 0)
        + df_pos.get("solar rooftop", 0)
        + df_pos.get("solar thermal", 0)
    )

    columns_to_drop = ["solar-hsat", "Solar", "solar rooftop", "solar thermal"]
    df_pos = df_pos.drop(
        columns=[col for col in columns_to_drop if col in df_pos.columns]
    )

    df_neg = df_neg[df_neg.sum().sort_values().index]

    fig, ax = plt.subplots(figsize=(14, 12))
    # Reorder the DataFrame columns based on the preferred order
    df_pos["Stromimport"] = (
        df["Electricity trade"].where(df["Electricity trade"] > 0).fillna(0)
    )
    df_pos = df_pos.drop(columns=["Electricity trade"], errors="ignore")
    df_pos = df_pos.rename(columns={"urban central H2 CHP": "H2 CHP"})
    df_pos["other"] = df_pos.drop(columns=preferred_order_pos, errors="ignore").sum(
        axis=1
    )

    df_pos = df_pos.reindex(columns=preferred_order_pos)
    ax = df_pos.plot.area(ax=ax, stacked=True, color=pos_c, linewidth=0.0)

    # rename negative values that are also present on positive side, so that they are not shown and plot negative values
    f = lambda c: "out_" + c
    cols = [f(c) if (c in df_pos.columns) else c for c in df_neg.columns]
    cols_map = dict(zip(df_neg.columns, cols))
    df_neg["Stromexport"] = (
        df["Electricity trade"].where(df["Electricity trade"] < 0).fillna(0)
    )
    df_neg = df_neg.drop(columns=["Electricity trade"], errors="ignore")
    df_neg["Sonstige"] = df_neg.drop(columns=preferred_order_neg, errors="ignore").sum(
        axis=1
    )
    df_neg = df_neg.reindex(columns=preferred_order_neg)
    ax = df_neg.plot.area(ax=ax, stacked=True, color=neg_c, linewidth=0.0)

    # plot lmps
    if plot_lmps:
        lmps = network.buses_t.marginal_price[
            network.buses[network.buses.carrier.isin(carriers)].index
        ].mean(axis=1)[period]
        ax2 = lmps.plot(
            style="--",
            color="black",
            label="Knotenpreise (gemittelt)",
            secondary_y=True,
        )
        ax2.grid(False)

        # Manually remove "(right)" from the legend
        handles, labels = ax2.get_legend_handles_labels()
        # Remove '(right)' from the labels, if present
        labels = [label.replace(" (right)", "") for label in labels]

        # Update the legend with the cleaned labels
        ax2.legend(handles, labels, loc="upper left")

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
        ax2.legend(loc="upper right")
        ax2.set_ylabel("Knotenpreise [€/MWh]")

    # explicitly filter out duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    filtered_handles_labels = [
        (h, l) for h, l in zip(handles, labels) if not l.startswith("out_")
    ]
    handles, labels = zip(*filtered_handles_labels)

    if nice_names & (not german_carriers):
        nice_names_dict = network.carriers.nice_name.to_dict()
        labels = [nice_names_dict.get(l, l) for l in labels]

    # rescale the y-axis
    ax.set_ylim([1.05 * df_neg.sum(axis=1).min(), 1.05 * df_pos.sum(axis=1).max()])

    # Get all handles and labels
    handles, labels = ax.get_legend_handles_labels()
    if german_carriers:
        labels = [carriers_in_german.get(l, l) for l in labels]

    # Split into "Erzeugung" and "Verbrauch"
    erzeugung_handles = handles[:17]
    erzeugung_labels = labels[:17]
    verbrauch_handles = handles[17:]
    verbrauch_labels = labels[17:]
    subtitle_erzeugung = Patch(color="none", label="Erzeugung")
    subtitle_verbrauch = Patch(color="none", label="Verbrauch")

    # Combine all handles and labels
    combined_handles = (
        [subtitle_erzeugung]
        + erzeugung_handles
        + [subtitle_verbrauch]
        + verbrauch_handles
    )
    combined_labels = (
        ["Erzeugung"] + erzeugung_labels + ["Verbrauch"] + verbrauch_labels
    )

    legend = ax.legend(
        combined_handles,
        combined_labels,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.2),
    )
    for text in legend.get_texts():
        text.set_fontsize(10)  # Adjust font size (e.g., 14)
    for patch in legend.legend_handles:
        patch.set_width(20)  # Adjust rectangle width
        patch.set_height(10)

    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%B"))

    german_months = [
        "Jan",
        "Feb",
        "März",
        "Apr",
        "Mai",
        "Juni",
        "Juli",
        "Aug",
        "Sept",
        "Okt",
        "Nov",
        "Dez",
    ]

    # Custom formatter function
    def format_month(x, pos):
        month = mdates.num2date(x).month  # Extract the month as an integer (1-12)
        return german_months[month - 1]  # Get the German month name

    # Apply custom formatter
    ax.xaxis.set_major_formatter(FuncFormatter(format_month))

    # Manually set ticks to avoid duplicating December
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[:-1])
    ax.set_xlabel("")
    ax.set_title(
        f"{title} {model_run}",
        fontsize=16,
        pad=15,
    )
    ax.grid(True)
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()

    return fig


def plot_nodal_heat_balance(
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
    german_carriers=False,
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

        if c in ["battery", "EV battery", "Li ion", "PHS", "hydro"]:
            ax1.plot(
                charge / 1e6,  # max_stor_cap,
                label=c,
                color=tech_colors[c],
                marker=markers[i],
                markevery=[0],
                mfc="white",
                mec="black",
            )

        else:
            ax2.plot(
                charge / 1e6,  # max_stor_cap,
                label=ger[c],
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
    # ax2.set_ylabel("Speicherstad [-]")
    ax2.set_ylabel("TWh")
    ax1.legend(loc="lower right")
    ax2.legend(loc="lower right")
    ax2.set_xticklabels(["Jan", "Mrz", "Mai", "Jul", "Sept", "Nov", "Jan"])

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
                f"Strompreisdauerlinien"
                if languange == "german"
                else f"Electricity price duration curves"
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

    # only plot 2030 onwards
    years = years[2:]
    networks = dict(islice(networks.items(), 2, None))
    year_colors = year_colors[2:]
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

    axes[i].set_xlabel("Strompreis [$€/MWh_{el}$]")
    plt.suptitle(f"Strompreise", fontsize=16, y=0.99)
    fig.tight_layout()
    fig.savefig(savepath, bbox_inches="tight")
    plt.close()

    return fig


def plot_backup_capacity(
    networks, tech_colors, savepath, backup_techs, vre_gens, region="DE"
):

    kwargs = {
        "groupby": networks[2020].statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }

    df_all = pd.DataFrame()

    for year in np.arange(2020, 2050, 5):

        n = networks[year]

        electricity_cap = (
            n.statistics.optimal_capacity(bus_carrier=["low voltage", "AC"], **kwargs)
            .filter(like=region)
            .groupby(["carrier"])
            .sum()
            .drop(
                ["AC", "DC", "electricity distribution grid"],
                errors="ignore",
            )
        )

        df = round(electricity_cap.sort_values(ascending=False) / 1e3, 2)
        non_vre_gens = df[~df.index.isin(vre_gens)].index
        df = df.loc[non_vre_gens]
        df = df[df > 1e-3]

        df_all = pd.concat([df_all, df], axis=1)

    df_all.columns = np.arange(2020, 2050, 5)

    tech_colors["coal"] = "black"

    # Create figure
    plt.figure(figsize=(18, 5))

    # Track x-axis positions
    x_positions = []
    x_labels = []
    x_offset = 0

    df = df_all

    # Plot for each group
    for group, techs_in_group in backup_techs.items():
        # Find matching techs in the dataframe
        matching_techs = [tech for tech in techs_in_group if tech in df.index]

        # Prepare data for this group
        group_data = df.loc[matching_techs]

        # Prepare bottom for stacking
        bottom = np.zeros(len(group_data.columns))

        # Plot each technology in the group
        for j, tech in enumerate(matching_techs):
            values = group_data.loc[tech].fillna(0)
            plt.bar(
                np.arange(len(group_data.columns)) + x_offset,
                values,
                1,
                bottom=bottom,
                label=tech,
                color=tech_colors[tech],
            )
            bottom += values

        # Store x-axis positions for this group
        x_positions.append(x_offset + len(group_data.columns) / 2)
        x_labels.append(group)

        # Update x offset with extra spacing
        x_offset += len(group_data.columns) + 1  # Add extra space between groups

    plt.ylim(0, 100)

    # Customize the plot
    plt.title("Kapazität Backup-Kraftwerke (Strom)", fontsize=22)
    plt.ylabel("GW", fontsize=16)

    # Create custom x-tick labels with group names below years
    x_ticks = np.concatenate(
        [
            np.arange(len(df.columns)) + i
            for i in range(0, x_offset, len(df.columns) + 1)
        ]
    )
    x_tick_labels = np.tile(df.columns, len(backup_techs))

    plt.xticks(x_ticks, x_tick_labels, rotation=45)

    # Add group labels below x-axis ticks
    for pos, label in zip(x_positions, x_labels):
        plt.text(
            pos,
            plt.gca().get_ylim()[0] - plt.gca().get_ylim()[1] * 0.15,
            label,
            horizontalalignment="center",
            verticalalignment="top",
            fontweight="bold",
            fontsize=16,
        )

    plt.grid(axis="y")

    # Replace legend labels with German names from carriers_in_german
    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels = [
        carriers_in_german.get(label, label) for label in labels
    ]  # Replace labels if in dict
    plt.legend(handles, new_labels, loc="upper center", ncol=7)

    plt.tight_layout()
    plt.savefig(savepath, bbox_inches="tight")


def plot_backup_generation(
    networks, tech_colors, savepath, backup_techs, vre_gens, region="DE"
):

    tech_colors["coal"] = "black"

    kwargs = {
        "groupby": networks[2020].statistics.groupers.get_name_bus_and_carrier,
        "nice_names": False,
    }

    df_all = pd.DataFrame()

    for year in np.arange(2020, 2050, 5):
        n = networks[year]

        electricity_supply_de = (
            n.statistics.supply(bus_carrier=["low voltage", "AC"], **kwargs)
            .filter(like=region)
            .groupby(["carrier"])
            .sum()
            .drop(
                ["AC", "DC", "electricity distribution grid"],
                errors="ignore",
            )
        )

        df = round(electricity_supply_de.sort_values(ascending=False) / 1e6, 2)
        non_vre_gens = df[~df.index.isin(vre_gens)].index
        df = df.loc[non_vre_gens]
        df = df[df > 0.01]
        df_all = pd.concat([df_all, df], axis=1)

    df_all.columns = np.arange(2020, 2050, 5)

    # Create figure
    plt.figure(figsize=(18, 5))

    # Track x-axis positions
    x_positions = []
    x_labels = []
    x_offset = 0

    df = df_all

    # Plot for each group
    for group, techs_in_group in backup_techs.items():
        # Find matching techs in the dataframe
        matching_techs = [tech for tech in techs_in_group if tech in df.index]

        # Prepare data for this group
        group_data = df.loc[matching_techs]

        # Prepare bottom for stacking
        bottom = np.zeros(len(group_data.columns))

        # Plot each technology in the group
        for j, tech in enumerate(matching_techs):
            values = group_data.loc[tech].fillna(0)
            plt.bar(
                np.arange(len(group_data.columns)) + x_offset,
                values,
                1,
                bottom=bottom,
                label=tech,
                color=tech_colors[tech],
            )
            bottom += values

        # Store x-axis positions for this group
        x_positions.append(x_offset + len(group_data.columns) / 2)
        x_labels.append(group)

        # Update x offset with extra spacing
        x_offset += len(group_data.columns) + 1  # Add extra space between groups

    plt.ylim(0, 200)

    # Customize the plot
    plt.title("Versorgung Backup-Kraftwerke (Strom)", fontsize=22)
    plt.ylabel("TWh", fontsize=16)

    # Create custom x-tick labels with group names below years
    x_ticks = np.concatenate(
        [
            np.arange(len(df.columns)) + i
            for i in range(0, x_offset, len(df.columns) + 1)
        ]
    )
    x_tick_labels = np.tile(df.columns, len(backup_techs))

    plt.xticks(x_ticks, x_tick_labels, rotation=45)

    # Add group labels below x-axis ticks
    for pos, label in zip(x_positions, x_labels):
        plt.text(
            pos,
            plt.gca().get_ylim()[0] - plt.gca().get_ylim()[1] * 0.15,
            label,
            horizontalalignment="center",
            verticalalignment="top",
            fontweight="bold",
            fontsize=16,
        )

    plt.grid(axis="y")

    # Replace legend labels with German names from carriers_in_german
    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels = [
        carriers_in_german.get(label, label) for label in labels
    ]  # Replace labels if in dict
    plt.legend(handles, new_labels, loc="upper center", ncol=7)

    plt.tight_layout()
    plt.savefig(savepath, bbox_inches="tight")


def plot_elec_prices_spatial(
    network, tech_colors, savepath, onshore_regions, year="2045", region="DE"
):

    # onshore_regions = gpd.read_file("/home/julian-geis/repos/pypsa-ariadne-1/resources/20241203-force-onwind-south-49cl-disc/KN2045_Bal_v4/regions_onshore_base_s_49.geojson")
    # onshore_regions = onshore_regions.set_index('name')

    n = network
    buses = n.buses[n.buses.carrier == "AC"].index

    df = onshore_regions
    df["elec_price"] = n.buses_t.marginal_price[buses].mean()

    # Netzentgelte, Annuität NEP 2045 - Annuität PyPSA 2045 / Stromverbrauch Pypsa 2045
    pypsa_netzentgelt = (15.82 - 6.53) / 1.237
    elec_price_de = df["elec_price"][df.index.str.contains("DE")]
    max_above_mean = elec_price_de.max() - elec_price_de.mean()

    mean_with_netzentgelt = elec_price_de.mean() + pypsa_netzentgelt
    # Calculate the difference from the mean_with_netzentgelt
    df["elec_price_diff"] = mean_with_netzentgelt - df["elec_price"]

    crs = ccrs.PlateCarree()

    # Calculate aspect ratio based on geographic extent
    extent = [5, 16, 47, 56]  # Germany bounds
    aspect_ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])

    # Set figure size dynamically based on aspect ratio
    fig_width = 16  # You can adjust this value
    fig_height = fig_width / aspect_ratio

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(
        1, 2, subplot_kw={"projection": crs}, figsize=(fig_width, 0.5 * fig_height)
    )

    # First subplot: elec_price
    mvs = df["elec_price"][df.index.str.contains("DE")]
    vmin = np.nanmin(mvs)
    vmax = np.nanmax(mvs)

    ax1.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
    ax1.coastlines(edgecolor="black", linewidth=0.5)
    ax1.set_facecolor("white")
    ax1.add_feature(cartopy.feature.OCEAN, color="azure")

    df[df.index.str.contains("DE")].to_crs(crs.proj4_init).plot(
        column="elec_price",
        ax=ax1,
        cmap=plt.get_cmap("cividis_r"),
        linewidth=0.05,
        edgecolor="grey",
        legend=True,
        vmin=vmin,
        vmax=vmax,
        legend_kwds={
            "label": "Strompreise ($€/MWh$)",
            "orientation": "vertical",
            "shrink": 0.9,
        },
    )

    # Set geographic extent for Germany
    ax1.set_extent(extent, crs)  # Germany bounds
    ax1.set_aspect(aspect_ratio)

    # Second subplot: elec_price_diff
    mvs_diff = df["elec_price_diff"][df.index.str.contains("DE")]

    ax2.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
    ax2.coastlines(edgecolor="black", linewidth=0.5)
    ax2.set_facecolor("white")
    ax2.add_feature(cartopy.feature.OCEAN, color="azure")

    df[df.index.str.contains("DE")].to_crs(crs.proj4_init).plot(
        column="elec_price_diff",
        ax=ax2,
        cmap=plt.get_cmap("Greens"),
        linewidth=0.05,
        edgecolor="grey",
        vmax=np.nanmax(mvs_diff) + np.nanmin(mvs_diff),
        vmin=0,
        legend=True,
        legend_kwds={
            "label": "Durchschnittliche Preisreduktion für Endkunden($€/MWh$)",
            "orientation": "vertical",
            "shrink": 0.9,
        },
    )

    # Set geographic extent for Germany
    ax2.set_extent(extent, crs)  # Germany bounds
    ax2.set_aspect(aspect_ratio)

    # Adjust layout to place legends on the right

    plt.suptitle(
        "Nodale Strompreise und durchschnittliche Preisreduktion",
        ha="center",
        fontsize=16,
        y=0.98,
        x=0.5,
    )
    plt.subplots_adjust(right=0.85)
    fig.tight_layout()
    plt.show()

    fig.savefig(savepath, bbox_inches="tight")


def plot_elec_prices_spatial_new(
    network, tech_colors, savepath, onshore_regions, year="2045", region="DE"
):

    # onshore_regions = gpd.read_file("/home/julian-geis/repos/pypsa-ariadne-1/resources/20241203-force-onwind-south-49cl-disc/KN2045_Bal_v4/regions_onshore_base_s_49.geojson")
    # onshore_regions = onshore_regions.set_index('name')

    n = network
    buses = n.buses[n.buses.carrier == "AC"].index

    df = onshore_regions
    df["elec_price"] = n.buses_t.marginal_price[buses].mean()

    # Netzentgelte, Annuität NEP 2045 - Annuität PyPSA 2045 / Stromverbrauch Pypsa 2045
    pypsa_netzentgelt = (6.53 + 27.51) / 1.237
    nep_netzentgelt = (15.82 + 27.51) / 1.237
    elec_price_de = df["elec_price"][df.index.str.contains("DE")]
    max_above_mean = elec_price_de.max() - elec_price_de.mean()

    # Calculate the difference from the mean_with_netzentgelt

    df["elec_price_nep"] = elec_price_de.mean() + nep_netzentgelt
    df["elec_price_pypsa"] = df["elec_price"] + pypsa_netzentgelt
    df["elec_price_diff"] = df["elec_price_nep"] - df["elec_price_pypsa"]

    crs = ccrs.PlateCarree()

    # Calculate aspect ratio based on geographic extent
    extent = [5, 16, 47, 56]  # Germany bounds
    aspect_ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])

    # Set figure size dynamically based on aspect ratio
    fig_width = 14  # You can adjust this value
    fig_height = fig_width / aspect_ratio

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(
        1, 2, subplot_kw={"projection": crs}, figsize=(fig_width, 0.55 * fig_height)
    )

    # First subplot: elec_price
    mvs = df["elec_price_pypsa"][df.index.str.contains("DE")]
    vmin = np.nanmin(mvs)
    vmax = max(np.nanmax(mvs), df["elec_price_nep"].unique().item())

    ax1.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
    ax1.coastlines(edgecolor="black", linewidth=0.5)
    ax1.set_facecolor("white")
    ax1.add_feature(cartopy.feature.OCEAN, color="azure")
    ax1.set_title("Durchschnittspreis, NEP Ausbau", pad=15)
    img1 = (
        df[df.index.str.contains("DE")]
        .to_crs(crs.proj4_init)
        .plot(
            column="elec_price_nep",
            ax=ax1,
            linewidth=0.05,
            edgecolor="grey",
            legend=False,
            vmin=vmin,
            vmax=vmax,
            cmap="viridis_r",
        )
    )

    # Set geographic extent for Germany
    ax1.set_extent(extent, crs)  # Germany bounds
    ax1.set_aspect(aspect_ratio)

    # Second subplot: elec_price_diff
    ax2.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
    ax2.coastlines(edgecolor="black", linewidth=0.5)
    ax2.set_facecolor("white")
    ax2.add_feature(cartopy.feature.OCEAN, color="azure")
    ax2.set_title("Regionale Preiszonen, $PyPSA$-$DE$ Ausbau", pad=15)

    img2 = (
        df[df.index.str.contains("DE")]
        .to_crs(crs.proj4_init)
        .plot(
            column="elec_price_diff",
            ax=ax2,
            cmap="viridis",
            linewidth=0.05,
            edgecolor="grey",
            vmax=vmax - vmin,
            vmin=0,
            legend=False,
        )
    )

    # Set geographic extent for Germany
    ax2.set_extent(extent, crs)  # Germany bounds
    ax2.set_aspect(aspect_ratio)

    # Create a new axis for the colorbar at the bottom
    cax2 = fig.add_axes([0.15, 0.14, 0.6, 0.03])  # [left, bottom, width, height]
    cax1 = fig.add_axes([0.15, 0.00, 0.6, 0.03])  # [left, bottom, width, height]

    # Add colorbar to the new axis
    cbar1 = fig.colorbar(
        img1.get_figure().get_axes()[0].collections[0],
        cax=cax1,
        orientation="horizontal",
    )
    cbar1.set_label("Börsenstrompreis zzgl. durchschnittlichem Netzentgelt [$€/MWh$]")
    cbar1.set_ticklabels(np.linspace(vmax, vmin, 6).round(1))
    cbar1.ax.invert_xaxis()

    cbar2 = fig.colorbar(
        img2.get_figure().get_axes()[0].collections[0],
        cax=cax2,
        orientation="horizontal",
    )
    cbar2.set_label("Durchschnittliche Preisreduktion für Endkunden [$€/MWh$]")
    cbar2.set_ticklabels(np.linspace(0, vmax - vmin, 6).round(1))

    plt.subplots_adjust(right=0.75, bottom=0.22)
    # plt.show()

    fig.savefig(savepath, bbox_inches="tight")


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
        bus_size_factor = 4e8
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
        bus_size_factor = 4e8
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
    ax.set_extent([5.5, 15.5, 47, 56], crs=ccrs.PlateCarree())

    if specify_buses is None:
        sizes = [5, 1]
        labels = [f"{s} GW" for s in sizes]
        sizes = [s / bus_size_factor * 1e3 for s in sizes]
        n_cols = 2
        title = ""
    elif specify_buses == "production":
        sizes = [50, 25, 5]
        labels = [f"{s} TWh" for s in sizes]
        sizes = [s / bus_size_factor * 1e6 for s in sizes]
        n_cols = 2
        title = "Wasserstoffinfrastruktur (Produktion)"
        loc_patches = (0.8, -0.09)
    elif specify_buses == "consumption":
        sizes = [50, 25, 5]
        labels = [f"{s} TWh" for s in sizes]
        sizes = [s / bus_size_factor * 1e6 for s in sizes]
        n_cols = 2
        title = "Wasserstoffinfrastruktur (Verbrauch)"
        loc_patches = (0.75, -0.16)

    legend_kw_circles = dict(
        loc="lower center",
        bbox_to_anchor=(0.1, -0.15),
        labelspacing=1.5,
        handletextpad=0.5,
        frameon=True,
        facecolor="white",
        fontsize=10,
        ncol=1,
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw_circles,
    )
    legend = ax.get_legend()
    legend.get_frame().set_boxstyle("square, pad=0.7")

    sizes = [30, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw_lines = dict(
        loc="lower center",
        bbox_to_anchor=(0.3, -0.07),
        frameon=True,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=10,
        ncol=1,
        facecolor="white",
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw_lines,
    )
    legend = ax.get_legend()
    legend.get_frame().set_boxstyle("square, pad=0.7")

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

    labels = [carriers_in_german.get(c, c) for c in labels]

    legend_kw_patches = dict(
        loc="lower center",
        bbox_to_anchor=loc_patches,
        ncol=n_cols,
        frameon=True,
        facecolor="white",
        fontsize=10,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw_patches)

    ax.set_facecolor("white")
    ax.set_title(title, fontsize=16, pad=20)
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
        title = "Stromnetzausbau (gesamt)"
    elif expansion_case == "startnetz":
        line_widths = startnetz_linew / linew_factor
        link_widths = startnetz_linkw / linkw_factor
        title = "Stromnetzausbau (Startnetz)"
    elif expansion_case == "pypsa":
        line_widths = total_exp_noStart_linew / linew_factor
        link_widths = total_exp_noStart_linkw / linkw_factor
        title = "Stromnetzausbau (nur modellseitig)"
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

    sizes = [40, 20, 10]
    labels = [f"{s} GW" for s in sizes]
    sizes = [s / bus_size_factor * 1e3 for s in sizes]

    legend_kw_circles = dict(
        loc="lower center",
        bbox_to_anchor=(0.15, -0.2),
        labelspacing=0.8,
        handletextpad=0.5,
        frameon=True,
        facecolor="white",
        fontsize=14,
        ncol=1,
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=m.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw_circles,
    )
    # ensure cirlce is not outside the box
    legend = ax.get_legend()
    legend.get_frame().set_boxstyle("square, pad=0.7")

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

    legend_kw_lines = dict(
        loc="lower center",
        bbox_to_anchor=(0.65, -0.12),  # 0.7
        frameon=True,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=14,
        ncol=2,
        facecolor="white",
    )

    add_legend_lines(ax, sizes, labels, colors=colors, legend_kw=legend_kw_lines)

    colors = [tech_colors[c] for c in carriers]
    labels = carriers

    legend_kw_patches = dict(
        loc="lower center",
        bbox_to_anchor=(0.66, -0.24),
        ncol=2,
        frameon=True,
        facecolor="white",
        fontsize=14,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw_patches)

    ax.set_title(title, fontsize=16, pad=20)

    fig.savefig(savepath, bbox_inches="tight")
    plt.close()


# electricity capacity map
def plot_cap_map_de(
    network,
    tech_colors,
    regions_de,
    savepath,
):

    m = network.copy()
    m.mremove("Bus", m.buses[m.buses.x == 0].index)
    m.buses.drop(m.buses.index[m.buses.carrier != "AC"], inplace=True)

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
    carriers_links = [
        "H2 Fuel Cell",
        "urban central H2 CHP",
        "H2 OCGT",
        "urban central solid biomass CHP",
        "urban central solid biomass CHP CC",
        "waste CHP",
        "waste CHP CC",
        "CCGT",
        "urban central gas CHP",
        "urban central gas CHP CC",
        "OCGT",
    ]
    elec = m.generators[
        (m.generators.carrier.isin(carriers)) & (m.generators.bus.str.contains("DE"))
    ].index
    elec_links = m.links[
        (m.links.carrier.isin(carriers_links)) & (m.links.index.str.contains("DE"))
    ].index
    bus_sizes = (
        m.generators.loc[elec, "p_nom_opt"]
        .groupby([m.generators.bus, m.generators.carrier])
        .sum()
        / bus_size_factor
    )
    bus_sizes_links = (
        m.links.loc[elec_links, "p_nom_opt"]
        .groupby([m.links.bus1, m.links.carrier])
        .sum()
        / bus_size_factor
    )
    bus_sizes = pd.concat([bus_sizes, bus_sizes_links])
    replacement_dict = {
        "onwind": "Onshore Wind",
        "offwind-ac": "Offshore Wind",
        "offwind-dc": "Offshore Wind",
        "solar": "Solar",
        "solar-hsat": "Solar",
        "H2 Fuel Cell": "H2",
        "urban central H2 CHP": "H2",
        "H2 OCGT": "H2",
        "urban central solid biomass CHP": "Biomasse",
        "urban central solid biomass CHP CC": "Biomasse",
        "waste CHP": "Abfall",
        "waste CHP CC": "Abfall",
        "CCGT": "Gas",
        "urban central gas CHP": "Gas",
        "urban central gas CHP CC": "Gas",
        "OCGT": "Gas",
    }
    tech_colors["Biomasse"] = tech_colors["biomass"]
    tech_colors["Abfall"] = tech_colors["waste"]
    tech_colors["Gas"] = tech_colors["gas"]
    bus_sizes = bus_sizes.rename(index=replacement_dict, level=1)
    bus_sizes = bus_sizes.groupby(level=[0, 1]).sum()
    carriers = bus_sizes.index.get_level_values(1).unique().tolist()

    regions_de = regions_de.to_crs(proj.proj4_init)
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": proj})

    m.plot(
        ax=ax,
        margin=0.06,
        bus_sizes=bus_sizes,
        bus_colors=tech_colors,
        line_alpha=0,
        link_alpha=0,
    )

    regions_de.plot(
        ax=ax,
        column="battery",
        cmap="Oranges",
        linewidths=0,
        legend=True,
        legend_kwds={
            "label": "Batteriespeicher [GWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    # Set geographic extent for Germany
    ax.set_extent([5.5, 15.5, 47, 56], crs=ccrs.PlateCarree())

    sizes = [50, 15, 5]
    labels = [f"  {s} GW" for s in sizes]
    sizes = [s / bus_size_factor * 1e3 for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1),
        borderaxespad=0.5,
        labelspacing=1.2,
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
    legend = ax.get_legend()
    legend.get_frame().set_boxstyle("square, pad=0.7")

    legend_kw = dict(
        loc=[0.2, 0.9],
        frameon=True,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=13,
        ncol=2,
        facecolor="white",
    )

    colors = [tech_colors[c] for c in carriers]
    labels = carriers
    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.5, 0),
        ncol=3,
        frameon=True,
        facecolor="white",
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)
    ax.set_title(f"Installierte Leistung Stromsektor {year}")
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
    tech_colors["H2 OCGT"] = "#3b4cc0"
    tech_colors["H2 retrofit OCGT"] = "#9abbff"
    tech_colors["urban central H2 CHP"] = "#c9d7f0"
    tech_colors["urban central H2 retrofit CHP"] = "#edd1c2"

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
        plot_nodal_elec_balance(
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

        plot_nodal_elec_balance(
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

        plot_nodal_elec_balance(
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
            plot_nodal_heat_balance(
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

            plot_nodal_heat_balance(
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

            plot_nodal_heat_balance(
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

    plot_backup_capacity(
        networks=networks_dict,
        tech_colors=tech_colors,
        savepath=snakemake.output.backup_capacity,
        backup_techs=backup_techs,
        vre_gens=vre_gens,
        region="DE",
    )

    plot_backup_generation(
        networks=networks_dict,
        tech_colors=tech_colors,
        savepath=snakemake.output.backup_generation,
        backup_techs=backup_techs,
        vre_gens=vre_gens,
        region="DE",
    )

    # load regions
    regions = gpd.read_file(snakemake.input.regions_onshore_clustered).set_index("name")

    year = 2045
    plot_elec_prices_spatial_new(
        network=networks[planning_horizons.index(year)].copy(),
        tech_colors=tech_colors,
        onshore_regions=regions,
        savepath=snakemake.output.elec_prices_spatial_de,
        region="DE",
        year=year,
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
        plot_cap_map_de(
            networks[planning_horizons.index(year)],
            tech_colors,
            regions_de,
            savepath=f"{snakemake.output.elec_transmission}/elec-cap-DE-{year}.png",
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
