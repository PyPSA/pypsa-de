# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts._helpers import mock_snakemake

TWh2PJ = 3.6


def plot_trade(
    df,
    savepath,
):
    # WARNING CODE DUPLICATION
    # WARNING i just asked COPILOT to reformat the two plots below
    # load data and convert to TWh
    h2_balance = (
        df.loc["Trade|Secondary Energy|Hydrogen|Volume"] / TWh2PJ * -1
    )  # exports-imports
    h2_import = df.loc["Trade|Secondary Energy|Hydrogen|Gross Import|Volume"] / TWh2PJ
    h2_export = h2_balance - h2_import

    elec_balance = (
        df.loc["Trade|Secondary Energy|Electricity|Volume"] / TWh2PJ * -1
    )  # exports-imports
    elec_import = (
        df.loc["Trade|Secondary Energy|Electricity|Gross Import|Volume"] / TWh2PJ
    )
    elec_export = elec_balance - elec_import

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Plot Electricity Trade
    ax1.bar(
        elec_import.columns, elec_import.loc["PJ/yr"], color="#70af1d", label="Import"
    )
    ax1.bar(
        elec_export.columns, elec_export.loc["PJ/yr"], color="#3f630f", label="Export"
    )
    ax1.scatter(
        elec_balance.columns,
        elec_balance.loc["PJ/yr"],
        color="black",
        marker="x",
        label="Netto",
    )
    ax1.set_xticks(elec_balance.columns)
    ax1.set_xticklabels(elec_balance.columns)
    ax1.axhline(0, color="black", linewidth=0.5)
    ax1.set_ylabel("TWh")
    ax1.set_title("Stromaustausch")
    ax1.legend(loc="upper left")

    # Plot Hydrogen Trade
    ax2.bar(h2_import.columns, h2_import.loc["PJ/yr"], color="#f081dc", label="Import")
    ax2.bar(h2_export.columns, h2_export.loc["PJ/yr"], color="#6b3161", label="Export")
    ax2.scatter(
        h2_balance.columns,
        h2_balance.loc["PJ/yr"],
        color="black",
        marker="x",
        label="Netto",
    )
    ax2.set_xticks(h2_balance.columns)
    ax2.set_xticklabels(h2_balance.columns)
    ax2.axhline(0, color="black", linewidth=0.5)
    ax2.set_ylabel("TWh")
    ax2.set_title("Wasserstoffaustausch")
    ax2.legend()

    plt.tight_layout()
    fig.savefig(savepath, bbox_inches="tight")


def plot_h2_trade(
    df,
    savepath,
):
    # load data and convert to TWh
    h2_balance = (
        df.loc["Trade|Secondary Energy|Hydrogen|Volume"] / TWh2PJ * -1
    )  # exports-imports
    h2_import = df.loc["Trade|Secondary Energy|Hydrogen|Gross Import|Volume"] / TWh2PJ
    h2_export = h2_balance - h2_import

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(h2_import.columns, h2_import.loc["PJ/yr"], color="#f081dc", label="Import")
    ax.bar(h2_export.columns, h2_export.loc["PJ/yr"], color="#6b3161", label="Export")
    ax.scatter(
        h2_balance.columns,
        h2_balance.loc["PJ/yr"],
        color="black",
        marker="x",
        label="Netto",
    )

    ax.set_xticks(h2_balance.columns)
    ax.set_xticklabels(h2_balance.columns)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylabel("TWh")
    ax.set_title("Wasserstoff Import/Export Deutschland")
    ax.legend()

    plt.tight_layout()
    fig.savefig(savepath, bbox_inches="tight")


def plot_elec_trade(
    df,
    savepath,
):
    # load data and convert to TWh
    elec_balance = (
        df.loc["Trade|Secondary Energy|Electricity|Volume"] / TWh2PJ * -1
    )  # exports-imports
    elec_import = (
        df.loc["Trade|Secondary Energy|Electricity|Gross Import|Volume"] / TWh2PJ
    )
    elec_export = elec_balance - elec_import

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(
        elec_import.columns, elec_import.loc["PJ/yr"], color="#70af1d", label="Import"
    )
    ax.bar(
        elec_export.columns, elec_export.loc["PJ/yr"], color="#3f630f", label="Export"
    )
    ax.scatter(
        elec_balance.columns,
        elec_balance.loc["PJ/yr"],
        color="black",
        marker="x",
        label="Netto",
    )

    ax.set_xticks(elec_balance.columns)
    ax.set_xticklabels(elec_balance.columns)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylabel("Twh")
    ax.set_title("Strom Import/Export Deutschland")
    ax.legend(loc="upper left")

    plt.tight_layout()
    fig.savefig(savepath, bbox_inches="tight")


def plot_Kernnetz(df, savepath=None, currency_year=2020):
    key = "Investment|Energy Supply|Hydrogen|Transmission and Distribution|"

    data = {
        "Kategorie": ["FNB", "PyPSA"],
        "Kernnetz-Zubau": [
            12.3 + 0.6,
            df.loc[key + "Kernnetz|New-build"].values.sum() * 5,
        ],
        "Kernnetz-Umstellung": [
            3.2 + 0.2,
            df.loc[key + "Kernnetz|Retrofitted"].values.sum() * 5,
        ],
        "Endogen-Zubau": [
            None,
            df.loc[key + "Endogen|New-build"].values.sum() * 5,
        ],
        "Endogen-Umstellung": [
            None,
            df.loc[key + "Endogen|Retrofitted"].values.sum() * 5,
        ],
        "PCI+IPCEI": [
            7.8,
            df.loc[key + "Kernnetz|PCI+IPCEI"].values.sum() * 5,
        ],
        "Not PCI+IPCEI": [
            8.4,
            df.loc[key + "Kernnetz|NOT-PCI+IPCEI"].values.sum() * 5,
        ],
    }

    plotframe = pd.DataFrame(data)
    plotframe.set_index("Kategorie", inplace=True)

    if currency_year == 2023:
        plotframe.loc["PyPSA"] *= 1.1076
    elif currency_year == 2020:
        plotframe.loc["FNB"] /= 1.1076
    else:
        raise ValueError("Currency year not supported")

    # Set up the plot
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))

    # Create bars
    bar_width = 0.35
    x = np.arange(2)  # Three groups

    group1 = [
        "Kernnetz-Zubau",
        "Kernnetz-Umstellung",
        "Endogen-Zubau",
        "Endogen-Umstellung",
    ]
    group2 = ["PCI+IPCEI", "Not PCI+IPCEI"]

    colors_dict = {
        # Group 1 - Using turquoise to pink spectrum
        "Kernnetz-Zubau": "#00A4B4",  # Turquoise
        "Kernnetz-Umstellung": "#20CFB4",  # Bright turquoise-mint
        "Endogen-Zubau": "#FF69B4",  # Hot pink
        "Endogen-Umstellung": "#BA55D3",  # Medium orchid (purple-pink)
        # Group 2 - Keeping warm colors
        "PCI+IPCEI": "#D95F02",  # Orange
        "Not PCI+IPCEI": "#E7B031",  # Golden yellow
    }
    # Create grouped bars
    y_limit = plotframe[group1].sum(axis=1).max() + 3
    plotframe[group1].plot(
        kind="bar",
        stacked=True,
        ax=ax[0],
        ylim=(0, y_limit),
        ylabel="Investment (Mrd. €)",
        color=[colors_dict.get(x, "#333333") for x in group1],
    )
    plotframe[group2].plot(
        kind="bar",
        stacked=True,
        ax=ax[1],
        ylim=(0, y_limit),
        color=[colors_dict.get(x, "#333333") for x in group2],
    )

    plt.suptitle("Investitionen ins Wasserstoffnetz")

    # Adjust layout
    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, bbox_inches="tight")
    else:
        plt.show()


def plot_NEP_Trassen(df, savepath=None, gleichschaltung=True):
    NEP_Trassen = {
        "NEP-DC": {
            # Zu-/Umbeseilung + Ersatz-/Parallelneubau + Neubau
            "Startnetz": 321 + 560 + 2466,
            "Zubaunetz": 0 + 179 + 4396,
        },
        "NEP-AC": {"Startnetz": 919 + 2081 + 599, "Zubaunetz": 2279 + 3846 + 1714},
    }

    key = "Length Additions|Electricity|Transmission|"

    NEP_Trassen["PyPSA-DC"] = {
        "exogen": df.loc[key + "DC|NEP"].values.sum(),
        "endogen": (df.loc[key + "DC"].values - df.loc[key + "DC|NEP"].values).sum(),
    }
    NEP_Trassen["PyPSA-AC"] = {
        "exogen": df.loc[key + "AC|NEP"].values.sum(),
        "endogen": (df.loc[key + "AC"].values - df.loc[key + "AC|NEP"].values).sum(),
    }

    data = {
        "Kategorie": ["DC", "AC"],
        "Startnetz": [
            NEP_Trassen["NEP-DC"]["Startnetz"],
            NEP_Trassen["NEP-AC"]["Startnetz"],
        ],
        "Zubaunetz": [
            NEP_Trassen["NEP-DC"]["Zubaunetz"],
            NEP_Trassen["NEP-AC"]["Zubaunetz"],
        ],
        "exogen": [
            NEP_Trassen["PyPSA-DC"]["exogen"],
            NEP_Trassen["PyPSA-AC"]["exogen"],
        ],
        "Übernahme": [
            0,
            NEP_Trassen["NEP-AC"]["Startnetz"] - NEP_Trassen["PyPSA-AC"]["exogen"],
        ],
        "endogen": [
            NEP_Trassen["PyPSA-DC"]["endogen"],
            NEP_Trassen["PyPSA-AC"]["endogen"],
        ],
    }

    plotframe = pd.DataFrame(data)

    # Define the width of the bars
    bar_width = 0.35
    indices = np.arange(len(plotframe))  # Bar positions

    plt.clf()

    plt.bar(indices, plotframe["Startnetz"], bar_width, label="Startnetz")
    plt.bar(
        indices,
        plotframe["Zubaunetz"],
        bar_width,
        bottom=plotframe["Startnetz"],
        label="Zubaunetz",
    )
    plt.bar(indices + bar_width, plotframe["exogen"], bar_width, label="exogen")
    bottom = plotframe["exogen"].copy()
    if gleichschaltung:
        plt.bar(
            indices + bar_width,
            plotframe["Übernahme"],
            bar_width,
            bottom=plotframe["exogen"],
            label="Übernahme",
            color="darkgreen",
        )
        bottom += plotframe["Übernahme"]
    plt.bar(
        indices + bar_width,
        plotframe["endogen"],
        bar_width,
        bottom=bottom,
        label="endogen",
    )

    plt.xlabel("Kategorie")
    plt.ylabel("km")
    plt.title("Trassenlänge Übertragungsnetz Onshore")

    # Adjust the x-ticks to be between the two bars
    plt.xticks(indices + bar_width / 2, plotframe["Kategorie"])
    plt.legend()
    if savepath:
        plt.savefig(savepath, bbox_inches="tight")
    else:
        plt.show()

    plotframe["NEP-Total"] = plotframe["Startnetz"] + plotframe["Zubaunetz"]
    plotframe["PyPSA-Total"] = (
        plotframe["exogen"] + plotframe["endogen"] + plotframe["Übernahme"]
    )
    plotframe.to_csv(snakemake.output.trassenlaenge_csv)


def plot_NEP(df, savepath=None, gleichschaltung=True, currency_year=2020):
    key = "Investment|Energy Supply|Electricity|Transmission|"

    data = {
        "Kategorie": ["DC", "AC", "System-\ndienstleistungen", "Onshore", "Offshore"],
        "Startnetz": [
            26,
            14.5,
            9.4,
            None,
            12.4,
        ],
        "Zubaunetz": [
            46.2,
            30.5,
            29.5,
            None,
            145.1,
        ],
        "exogen": [
            df.loc[key + "DC|NEP|Onshore"].values.sum() * 5,
            df.loc[key + "AC|NEP|Onshore"].values.sum() * 5,
            0,  # see "Übernahme"
            None,
            df.loc[key + "NEP|Offshore"].values.sum() * 5,
        ],
        "endogen": [
            (
                df.loc[key + "DC|Onshore"].values
                - df.loc[key + "DC|NEP|Onshore"].values
            ).sum()
            * 5,
            (
                df.loc[key + "AC|Onshore"].values
                - df.loc[key + "AC|NEP|Onshore"].values
            ).sum()
            * 5,
            0,
            None,
            (
                df.loc[key + "Offshore"].values - df.loc[key + "NEP|Offshore"].values
            ).sum()
            * 5,
        ],
        "Übernahme": [
            0,
            df.loc[key + "AC|Übernahme|Startnetz Delta"].values.sum() * 5,
            df.loc[key + "AC|Übernahme|Reactive Power Compensation"].values.sum() * 5,
            None,
            0,
        ],
    }

    plotframe = pd.DataFrame(data)

    if currency_year == 2023:
        plotframe["exogen"] *= 1.1076
        plotframe["endogen"] *= 1.1076
        plotframe["Übernahme"] *= 1.1076

    elif currency_year == 2020:
        plotframe["Startnetz"] /= 1.1076
        plotframe["Zubaunetz"] /= 1.1076
    else:
        raise ValueError("Currency year not supported")

    plotframe.set_index("Kategorie", inplace=True)
    plotframe.loc["Onshore"] = plotframe.loc[
        ["AC", "DC", "System-\ndienstleistungen"]
    ].sum()
    # Define the width of the bars
    bar_width = 0.35
    indices = np.arange(len(plotframe))  # Bar positions

    plt.clf()

    plt.bar(indices, plotframe["Startnetz"], bar_width, label="Startnetz")
    plt.bar(
        indices,
        plotframe["Zubaunetz"],
        bar_width,
        bottom=plotframe["Startnetz"],
        label="Zubaunetz",
    )
    plt.bar(indices + bar_width, plotframe["exogen"], bar_width, label="exogen")
    bottom = plotframe["exogen"].copy()
    if gleichschaltung:
        plt.bar(
            indices + bar_width,
            plotframe["Übernahme"],
            bar_width,
            bottom=plotframe["exogen"],
            label="Übernahme",
            color="darkgreen",
        )
        bottom += plotframe["Übernahme"]
    plt.bar(
        indices + bar_width,
        plotframe["endogen"],
        bar_width,
        bottom=bottom,
        label="endogen",
    )

    plt.xlabel("Kategorie")
    plt.ylabel(f"Milliarden EUR{currency_year}")
    plt.title(f"Investitionen ins Übertragungsnetz in EUR{currency_year}")

    # Adjust the x-ticks to be between the two bars
    plt.xticks(indices + bar_width / 2, plotframe.index)
    plt.legend()

    # Rename the category to remove the format string again
    plotframe.rename(
        index={"System-\ndienstleistungen": "Systemdienstleistungen"}, inplace=True
    )
    plotframe["NEP-Total"] = plotframe["Startnetz"] + plotframe["Zubaunetz"]
    plotframe["PyPSA-Total"] = (
        plotframe["exogen"] + plotframe["endogen"] + plotframe["Übernahme"]
    )

    # Add total costs annotations

    plt.text(
        1.05,
        0.95,
        f"NEP: {round(plotframe.loc["Onshore","NEP-Total"] + plotframe.loc["Offshore","NEP-Total"],1)}",
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.5),
    )
    plt.text(
        1.05,
        0.85,
        f"PyPSA: {round(plotframe.loc["Onshore","PyPSA-Total"] + plotframe.loc["Offshore","PyPSA-Total"],1)}",
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.5),
    )

    if savepath:
        plt.savefig(savepath, bbox_inches="tight")
    else:
        plt.show()

    plotframe.to_csv(snakemake.output.transmission_investment_csv)


def secondary_energy_plot(ddf, name="Secondary Energy"):
    # Get Secondary Energy data
    df = ddf[ddf.index.get_level_values("Variable").str.startswith(name)]
    # Get Unit and delete Multiindex
    unit = df.index.get_level_values("Unit").unique().dropna().item()
    df.index = df.index.droplevel("Unit")

    # Simplify variable names
    df.index = pd.Index(
        map(
            lambda x: x[(x.find("|") + 1) :],
            df.index,
        ),
        name=df.index.names[0],
    )

    # Get detailed data
    exclude = ("Fossil", "Renewables", "Losses")  # include Losses once fixed
    detailed = df[(df.index.str.count("[|]") == 1) & ~df.index.str.endswith(exclude)]

    ax = detailed.T.plot.area(ylabel=unit, title="Detailed" + name)
    ax.legend(bbox_to_anchor=(1, 1))

    coarse = df[(df.index.str.count("[|]") == 0)]

    ax = coarse.T.plot.area(ylabel=unit, title=name)
    ax.legend(bbox_to_anchor=(1, 1))


def ariadne_subplot(
    df,
    ax,
    title,
    select_regex=r"",
    drop_regex=r"",
    stacked=True,
    unit=None,
):
    df = df.T.copy()

    if select_regex:
        df = df.filter(
            regex=select_regex,
        )
    if drop_regex:
        df = df.filter(
            regex=drop_regex,
        )

    # Check that all values have the same Unit
    if not unit:
        unit = df.columns.get_level_values("Unit").unique().dropna().item()

    df.columns = df.columns.droplevel("Unit")

    # Simplify variable names
    df.columns = pd.Index(
        map(
            lambda x: x[(x.find("|") + 1) :],
            df.columns,
        ),
        name=df.columns.names[0],
    )

    if df.empty:
        # Create an empty plot if DataFrame is empty
        ax.plot([], [])
        ax.set_title("Ooops! Empty DataFrame")
        return ax

    return df.plot.area(ax=ax, title=title, legend=False, stacked=stacked, ylabel=unit)


def side_by_side_plot(df, dfhybrid, title, savepath, rshift=1.25, **kwargs):
    idx = df.index.union(dfhybrid.index, sort=False)

    df = df.reindex(idx)
    dfhybrid = dfhybrid.reindex(idx)

    fig, axes = plt.subplots(ncols=2, sharey=True)
    ax = ariadne_subplot(df, axes[0], "PyPSA-Eur", **kwargs)
    ax2 = ariadne_subplot(dfhybrid, axes[1], "REMIND-EU v1.1", **kwargs)

    handles, labels = ax.get_legend_handles_labels()
    labels2 = ax2.get_legend_handles_labels()[1]
    assert labels == labels2

    fig.legend(
        list(reversed(handles)), list(reversed(labels)), bbox_to_anchor=(rshift, 0.9)
    )
    fig.suptitle(title)
    title = title.replace(" ", "_")
    fig.savefig(savepath, bbox_inches="tight")
    return fig


def within_plot(
    df,
    df2,
    title,
    savepath,
    select_regex=r"",
    drop_regex=r"",
    write_sum=False,
    unit="EUR_2020/GJ",
    **kwargs,
):
    df = df.T.copy()
    if select_regex:
        df = df.filter(
            regex=select_regex,
        )
    if drop_regex:
        df = df.filter(
            regex=drop_regex,
        )

    n = df.shape[1]
    if n == 0:
        print(
            f"Warning! Apparently the variables required for this plot (({title}) are missing."
        )
        fig = plt.figure()
        plt.title(
            f"Warning! Apparently the variables required for this plot ({title}) are missing."
        )
        fig.savefig(savepath, bbox_inches="tight")
        return fig
    rows = n // 2 + n % 2

    fig, axes = plt.subplots(rows, 2, figsize=(10, 5 * rows))
    axes = axes.flatten()

    for i, var in enumerate(df.columns.get_level_values("Variable")):
        axes[i].plot(df.xs(var, axis=1, level=0), label="PyPSA-Eur")
        if var in df2.index.get_level_values("Variable"):
            axes[i].plot(df2.T.xs(var, axis=1, level=0), label="REMIND-EU")
        axes[i].set_title(var)
        axes[i].legend()

        if write_sum:
            sum_df1 = round(df.xs(var, axis=1, level=0).sum().values.item(), 2)
            if var in df2.index.get_level_values("Variable"):
                sum_df2 = round(df2.T.xs(var, axis=1, level=0).sum().values.item(), 2)
            else:
                sum_df2 = np.nan
            # Annotate plot with the sum of variables
            sum_text = f"Sum: \nPyPSA-Eur = {sum_df1},\nREMIND-EU = {sum_df2}"
            axes[i].annotate(
                sum_text,
                xy=(0, 1),
                xycoords="axes fraction",
                fontsize=12,
                xytext=(5, -5),
                textcoords="offset points",
                ha="left",
                va="top",
            )

    # Remove the last subplot if there's an odd number of plots
    if n % 2 != 0:
        fig.delaxes(axes[-1])

    plt.suptitle(f"{title} in ({unit})", fontsize="xx-large", y=1.0)
    plt.tight_layout()
    plt.close()
    fig.savefig(savepath, bbox_inches="tight")

    return fig


def elec_val_plot(df, savepath):
    # electricity validation for 2020
    elec_capacities = pd.DataFrame(
        index=[
            "ror",
            "hydro",
            "battery",
            "biomass",
            "nuclear",
            "lignite",
            "coal",
            "oil",
            "gas",
            "wind_onshore",
            "wind_offshore",
            "solar",
        ]
    )
    elec_generation = pd.DataFrame(
        index=[
            "net exports",
            "ror",
            "hydro",
            "battery",
            "biomass",
            "nuclear",
            "lignite",
            "coal",
            "oil",
            "gas",
            "wind",
            "solar",
        ]
    )

    elec_capacities["real"] = [
        4.94,
        9.69,
        2.4,
        8.72,
        8.11,
        20.86,
        23.74,
        4.86,
        32.54,
        54.25,
        7.86,
        54.36,
    ]  # https://energy-charts.info/charts/installed_power/chart.htm?l=en&c=DE&year=2020

    elec_capacities["pypsa"] = [
        0,
        df.loc[("Capacity|Electricity|Hydro", "GW"), "2020"],
        0,
        df.loc[("Capacity|Electricity|Biomass", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Nuclear", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Coal|Lignite", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Coal|Hard Coal", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Oil", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Gas", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Wind|Onshore", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Wind|Offshore", "GW"), "2020"],
        df.loc[("Capacity|Electricity|Solar", "GW"), "2020"],
    ]

    elec_generation["real (gross)"] = [
        -18.9,  # net exports",
        np.nan,  # ror
        18.7,  # hydro
        np.nan,  # battery
        45,  # biomass
        64,  # nuclear
        91,  # lignite
        43,  # coal
        4.7,  # oil
        95,  # gas
        132,  # wind
        50,  # solar
    ]
    # https://www.destatis.de/DE/Themen/Branchen-Unternehmen/Energie/Erzeugung/Tabellen/bruttostromerzeugung.html
    # https://www.bdew.de/media/documents/Bruttostromerz_D_Entw_10J_online_o_dw2x_jaehrlich_FS_05042024_nlA6lUa.pdf

    # obtained so that it fits how pypsa models the energy sector
    elec_generation["real (net, pypsa representation)"] = [
        -18.9,  # 34 TWH in 2019             "net exports",
        0,  # "ror" due to reporting
        18.54,  # "hydro"
        np.nan,  # "battery"
        44.85,  # "biomass"
        60.91,  # "nuclear"
        84.5,  # "lignite"  (82.13 + 2.37 (industrial own production)
        38.7,  # "coal"
        3.71,  #  "oil"
        91.7,  # "gas" 34.6 (industry self consumption) + 57.1 (fossil gas grid feed in)
        129.64,  # "wind"
        48.5,  # "solar"
    ]
    # https://energy-charts.info/charts/energy_pie/chart.htm?l=de&c=DE&interval=year&year=2020
    # https://www.bundesnetzagentur.de/SharedDocs/Pressemitteilungen/DE/2021/20210102_smard.html
    # https://energy-charts.info/charts/energy/chart.htm?l=en&c=DE&interval=year&year=2020&source=total

    elec_generation["pypsa (net)"] = [
        -df.loc[("Trade|Secondary Energy|Electricity|Volume", "PJ/yr"), "2020"] / 3.6,
        0,
        df.loc[("Secondary Energy|Electricity|Hydro", "PJ/yr"), "2020"] / 3.6,
        0,
        df.loc[("Secondary Energy|Electricity|Biomass", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Nuclear", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Coal|Lignite", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Coal|Hard Coal", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Oil", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Gas", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Wind", "PJ/yr"), "2020"] / 3.6,
        df.loc[("Secondary Energy|Electricity|Solar", "PJ/yr"), "2020"] / 3.6,
    ]

    # elec_generation.loc["sum/10"] = elec_generation.sum().div(10)
    elec_generation.loc["sum_real_gross-sum_pypsa", "sum_real_gross-sum_pypsa"] = (
        elec_generation.sum()["real (gross)"] - elec_generation.sum()["pypsa (net)"]
    )
    elec_generation.loc["sum_real_net-sum_pypsa", "sum_real_net-sum_pypsa"] = (
        elec_generation.sum()["real (net, pypsa representation)"]
        - elec_generation.sum()["pypsa (net)"]
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), width_ratios=[1, 1.5])
    elec_capacities.plot(kind="bar", ax=axes[0])
    axes[0].set_ylabel("GW")
    axes[0].set_title("Installed Capacities Germany 2020")

    elec_generation.plot(kind="bar", ax=axes[1])
    axes[1].set_ylabel("TWh")
    axes[1].set_title("Electricity Generation Germany 2020")

    plt.tight_layout()
    plt.close()
    fig.savefig(savepath, bbox_inches="tight")

    return fig


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_ariadne_variables",
            simpl="",
            clusters=49,
            opts="",
            ll="v1.2",
            sector_opts="None",
            planning_horizons="2045",
            run="KN2045_Bal_v4",
            # configfiles="config/config.public.yaml"
        )

    df = (
        pd.read_excel(
            snakemake.input.exported_variables_full,
            index_col=list(range(5)),
            # index_col=["Model", "Scenario", "Region", "Variable", "Unit"],
            sheet_name="data",
        )
        .groupby(["Variable", "Unit"], dropna=False)
        .sum()
    ).round(5)

    df.columns = df.columns.astype(str)
    leitmodell = "REMIND-EU v1.1"

    dfremind = pd.read_csv(
        snakemake.input.ariadne_database,
        index_col=["model", "scenario", "region", "variable", "unit"],
    ).loc[leitmodell, snakemake.params.fallback_reference_scenario, "Deutschland"][
        df.columns
    ]
    dfremind.index.names = df.index.names

    side_by_side_plot(
        df,
        dfremind,
        "Primary Energy in PJ_yr",
        savepath=snakemake.output.primary_energy,
        select_regex=r"Primary Energy\|[^|]*$",
        drop_regex=r"^(?!.*(Fossil|Price)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Primary Energy in PJ_yr",
        savepath=snakemake.output.primary_energy_detailed,
        select_regex=r"Primary Energy\|[^|]*\|[^|]*$",
        drop_regex=r"^(?!.*(CCS|Price|Volume)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Secondary Energy in PJ_yr",
        savepath=snakemake.output.secondary_energy,
        select_regex=r"Secondary Energy\|[^|]*$",
        drop_regex=r"^(?!.*(Price)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Secondary Energy in PJ_yr",
        savepath=snakemake.output.secondary_energy_detailed,
        # Secondary Energy|Something|Something (exactly two pipes)
        select_regex=r"Secondary Energy\|[^|]*\|[^|]*$",
        # Not ending in Fossil or Renewables (i.e., categories)
        drop_regex=r"^(?!.*(Fossil|Renewables|Losses|Price|Volume)).+",
    )

    if df.loc["Final Energy|Industry excl Non-Energy Use|Hydrogen", "2020"].item() < 0:
        val = df.loc["Final Energy|Industry excl Non-Energy Use|Hydrogen", "2020"]
        df.loc["Final Energy|Industry excl Non-Energy Use|Hydrogen", "2020"] = 0
        df.loc["Final Energy|Hydrogen", "2020"] = 0
        print("WARNING! NEGATIVE HYDROGEN DEMAND IN INDUSTRY IN 2020! ", val)
    if df.loc["Final Energy|Industry excl Non-Energy Use|Hydrogen", "2025"].item() < 0:
        val = df.loc["Final Energy|Industry excl Non-Energy Use|Hydrogen", "2025"]
        df.loc["Final Energy|Industry excl Non-Energy Use|Hydrogen", "2025"] = 0
        df.loc["Final Energy|Hydrogen", "2025"] = 0
        print("WARNING! NEGATIVE HYDROGEN DEMAND IN INDUSTRY IN 2025! ", val)
    side_by_side_plot(
        df,
        dfremind,
        "Final Energy in PJ_yr",
        savepath=snakemake.output.final_energy,
        select_regex=r"Final Energy\|[^|]*$",
        rshift=1.45,
        drop_regex=r"^(?!.*(Price|Non-Energy Use)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Final Energy in PJ_yr",
        savepath=snakemake.output.final_energy_detailed,
        select_regex=r"Final Energy\|[^|]*\|[^|]*$",
        rshift=1.7,
        drop_regex=r"^(?!.*(Price|\|Solids\||Non-Energy Use\|)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Capacity in GW",
        savepath=snakemake.output.capacity,
        select_regex=r"Capacity\|[^|]*$",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Capacity in GW",
        savepath=snakemake.output.capacity_detailed,
        select_regex=r"Capacity\|[^|]*\|[^|]*$",
        drop_regex=r"^(?!.*(Reservoir|Converter)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Demand Emissions in Mt",
        savepath=snakemake.output.energy_demand_emissions,
        select_regex=r"Emissions\|CO2\|Energy\|Demand\|[^|]*$",
        stacked=False,
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Supply Emissions in Mt",
        savepath=snakemake.output.energy_supply_emissions,
        select_regex=r"Emissions\|CO2\|Energy\|Supply\|[^|]*$",
        stacked=False,
        drop_regex=r"^(?!.*(and)).+",
    )

    side_by_side_plot(
        df,
        dfremind,
        "Detailed Emissions in Mt",
        savepath=snakemake.output.co2_emissions,
        select_regex=r"Emissions\|CO2\|[^|]*$",
        stacked=False,
        # drop_regex=r"^(?!.*(and)).+",
        unit="Mt CO2-equiv/yr",
    )

    within_plot(
        df,
        dfremind,
        title="Price|Primary Energy",
        savepath=snakemake.output.primary_energy_price,
        select_regex=r"Price\|Primary Energy\|[^|]*$",
    )

    within_plot(
        df[
            df.index.get_level_values("Variable").str.startswith(
                "Price|Secondary Energy"
            )
        ],
        dfremind,
        title="Price|Secondary Energy",
        savepath=snakemake.output.secondary_energy_price,
    )

    # within_plot(
    #     df[df.index.get_level_values("Variable").str.startswith("Price|Final Energy|Residential")],
    #     dfremind,
    #     title = "Price|Final Energy|Residential",
    #     savepath=snakemake.output.final_energy_residential_price,
    #     #select_regex=r"Price\|Final Energy\|Residential\|[^|]*$"
    # )

    within_plot(
        df[
            df.index.get_level_values("Variable").str.startswith(
                "Price|Final Energy|Industry"
            )
        ],
        dfremind,
        title="Price|Final Energy|Industry",
        savepath=snakemake.output.final_energy_industry_price,
        # select_regex=r"Price\|Final Energy\|Industry\|[^|]*$"
    )

    within_plot(
        df[
            df.index.get_level_values("Variable").str.startswith(
                "Price|Final Energy|Transportation"
            )
        ],
        dfremind,
        title="Price|Final Energy|Transportation",
        savepath=snakemake.output.final_energy_transportation_price,
        # select_regex=r"Price\|Final Energy\|Industry\|[^|]*$"
    )

    within_plot(
        df[
            df.index.get_level_values("Variable").str.startswith(
                "Price|Final Energy|Residential and Commercial"
            )
        ],
        dfremind,
        title="Price|Final Energy|Residential and Commercial",
        savepath=snakemake.output.final_energy_residential_commercial_price,
        # select_regex=r"Price\|Final Energy\|Industry\|[^|]*$"
    )

    within_plot(
        df[df.index.get_level_values("Variable").str.startswith("Price")],
        dfremind,
        title="All prices",
        savepath=snakemake.output.all_prices,
    )

    within_plot(
        df[df.index.get_level_values("Variable").str.startswith("Price|Carbon")],
        dfremind,
        title="Price of carbon",
        savepath=snakemake.output.policy_carbon,
        unit="EUR/tCO2",
    )

    within_plot(
        df[
            df.index.get_level_values("Variable").str.startswith(
                "Investment|Energy Supply"
            )
        ],
        dfremind,
        title="Investment in Energy Supply",
        savepath=snakemake.output.investment_energy_supply,
        unit="Billionen EUR",
        write_sum=True,
    )

    elec_val_plot(df, savepath=snakemake.output.elec_val_2020)

    within_plot(
        df[df.index.get_level_values("Variable").str.startswith("Trade")],
        dfremind,
        title="Trade",
        savepath=snakemake.output.trade,
        unit="PJ/yr",
    )

    plot_NEP(df, savepath=snakemake.output.NEP_plot)
    plot_NEP_Trassen(df, savepath=snakemake.output.NEP_Trassen_plot)

    plot_Kernnetz(
        df, savepath=snakemake.output.Kernnetz_Investment_plot, currency_year=2020
    )

    plot_elec_trade(
        df,
        savepath=snakemake.output.elec_trade,
    )

    plot_h2_trade(
        df,
        savepath=snakemake.output.h2_trade,
    )

    plot_trade(
        df,
        savepath=snakemake.output.trade_balance,
    )
