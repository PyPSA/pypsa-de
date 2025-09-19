import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import collections
import itertools
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, mock_snakemake

groups = {
    "gas": ["gas CHP", "OCGT", "CCGT", "gas"],
    "heat vent": ["heat vent"],
    "water tanks": ["water tank", "water pit"],
    "heat pump": ["heat pump"],
    "resistive heater": ["resistive heater"],
    "biomass": ["biomass"],
    "lignite": ["lignite"],
    "coal": ["coal"],
    "oil": ["oil"],
    "waste": ["waste"],
    "solar": ["solar"],
    "offwind": ["offwind"],
}


def aggregate_by_keywords(opex_comp_agg, groups):
    """
    Aggregate rows in opex_comp_agg according to keyword groups.

    Parameters
    ----------
    opex_comp_agg : pd.DataFrame
        DataFrame with row index as technology names.
    groups : dict
        Keys = new aggregated name,
        Values = list of substrings to match in the index.

    Returns
    -------
    pd.DataFrame
    """
    df_out = opex_comp_agg.copy()
    for new_name, keywords in groups.items():
        mask = df_out.index.to_series().str.contains("|".join(keywords))
        if mask.any():
            summed = df_out.loc[mask].sum()
            df_out = df_out.drop(df_out.index[mask])
            df_out.loc[new_name] = summed
    return df_out


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "regret_plots",
        )

    configure_logging(snakemake)
    config = snakemake.config
    planning_horizons = snakemake.params.planning_horizons
    scenarios = ["HighDemand", "LowDemand"]
    tech_colors = snakemake.params.plotting["tech_colors"]

    # Nested dict: networks[year][scenario][decision] = Network
    networks = collections.defaultdict(lambda: collections.defaultdict(dict))

    for fn in snakemake.input.regret_networks:
        parts = fn.split(os.sep)

        # scenario is the folder name 2 levels up
        scenario = parts[-3]
        if scenario not in scenarios:
            raise ValueError(
                f"Unexpected scenario '{scenario}' in {fn}. Allowed: {scenarios}"
            )

        # extract year (4 digits before .nc)
        m = re.search(r"_(\d{4})\.nc$", fn)
        if not m:
            raise ValueError(f"Could not parse year from {fn}")
        year = int(m.group(1))

        # extract decision_* (string until the 2nd underscore in filename)
        filename = parts[-1]
        decision = "_".join(filename.split("_")[:2])
        if not decision.startswith("decision"):
            raise ValueError(f"Unexpected decision string in {filename}")

        # load and store
        # print(f"Loading {fn} ...")
        # print(f"  scenario: {scenario}, year: {year}, decision: {decision}")
        networks[year][scenario][decision] = pypsa.Network(fn)

    # ensure output directory exist
    if not os.path.exists(snakemake.params.output_dir):
        os.makedirs(snakemake.params.output_dir)

    # Plot electricity price duration curves

    fig, ax = plt.subplots(
        figsize=(10, 5 * len(planning_horizons)), nrows=len(planning_horizons), ncols=1
    )
    ax = ax.flatten()

    decisions = ["decision_HighDemand", "decision_LowDemand"]

    for i, year in enumerate(planning_horizons):
        for scenario, decision in itertools.product(scenarios, decisions):
            n = networks[year][scenario][decision]
            lmps = n.buses_t.marginal_price.loc[
                :, (n.buses.carrier == "AC") & (n.buses.index.str.startswith("DE"))
            ]
            lmps_sorted = pd.DataFrame(
                lmps.values.flatten(), columns=["lmp"]
            ).sort_values(by="lmp", ascending=False)
            lmps_sorted["percentage"] = (
                np.arange(len(lmps_sorted)) / len(lmps_sorted) * 100
            )

            ax[i].plot(
                lmps_sorted["percentage"],
                lmps_sorted["lmp"],
                label=f"{scenario}_{decision} (avg: {lmps_sorted['lmp'].mean():.2f})",
            )

        ax[i].set_ylim(-50, 300)
        ax[i].legend()
        ax[i].set_xlabel("Percentage of time")
        ax[i].set_ylabel("€/MWh")
        ax[i].set_title(f"Price duration curves {year}")

    plt.tight_layout()
    plt.savefig(snakemake.output.elec_price_comp_de, bbox_inches="tight")
    plt.close()

    # Print CO2 prices

    # for i, year in enumerate(years):
    #     for scenario, decision in itertools.product(scenarios, decisions):

    #         n = networks[year][scenario][decision]

    #         print(f"CO2 price for {year}, {scenario}, {decision}: {n.global_constraints.loc["CO2Limit", "mu"] + n.global_constraints.loc["co2_limit-DE", "mu"]}")

    # Plot OPEX

    kwargs = {
        "groupby": ["bus", "carrier"],
        "at_port": True,
        "nice_names": False,
    }

    fig, axes = plt.subplots(
        nrows=len(planning_horizons), ncols=1, figsize=(12, 6 * len(planning_horizons))
    )
    axes = axes.flatten()

    for i, year in enumerate(planning_horizons):
        opex_comp = pd.DataFrame(
            columns=["_".join(tup) for tup in itertools.product(scenarios, decisions)]
        )

        # Collect OPEX for all scenario-decision combinations
        for scenario, decision in itertools.product(scenarios, decisions):
            n = networks[year][scenario][decision]

            opex = (
                n.statistics.opex(**kwargs)
                .filter(like="DE")
                .groupby("carrier")
                .sum()
                .multiply(1e-9)  # to billion €
            )
            opex_comp[f"{scenario}_{decision}"] = opex

        # Aggregate cost components with less than 0.1 (100 Mio €) as "Other"
        opex_comp_agg = aggregate_by_keywords(opex_comp, groups)
        small_rows = opex_comp_agg.abs().max(axis=1) < 0.1
        other_row = opex_comp_agg[small_rows].sum(axis=0)
        opex_comp_agg = opex_comp_agg.loc[~small_rows]
        opex_comp_agg.loc["Other"] = other_row

        # Prepare labels with line breaks
        labels = [col.replace("_", "\n") for col in opex_comp_agg.columns]

        # Plot stacked bar
        ax = axes[i]
        bottom = np.zeros(len(opex_comp_agg.columns))

        for tech in opex_comp_agg.index:
            values = opex_comp_agg.loc[tech].values
            ax.bar(
                labels,
                values,
                bottom=bottom,
                color=tech_colors.get(tech, "#333333"),
                label=tech,
            )

            # Add numbers in the middle, except for 'Other'
            if tech != "Other":
                for j, val in enumerate(values):
                    if val > 0:  # only if positive
                        ax.text(
                            j,
                            bottom[j] + val / 2,  # middle of the segment
                            f"{val:.2f}",
                            ha="center",
                            va="center",
                            fontsize=8,
                            color="white",
                        )

            bottom += values

        # Add total sum labels on top of bars
        totals = opex_comp_agg.sum(axis=0)
        for j, total in enumerate(totals):
            ax.text(
                j,
                total + total * 0.02,
                f"{total:.2f}",
                ha="center",
                va="bottom",
                fontsize=10,
            )

        # Adjust y-limit
        ax.set_ylim(0, max(totals) * 1.08)
        ax.set_ylabel("OPEX [billion €]")
        ax.set_title(f"Stacked OPEX composition by technology, {year}")

    # Legend outside
    axes[-1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.savefig(snakemake.params.output_dir + "/opex_comp_de.png", bbox_inches="tight")
    plt.close()
