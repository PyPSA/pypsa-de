import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))  # Adds 'scripts/' to path
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)  # Adds repo root

import re
from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, mock_snakemake
from matplotlib.patches import Rectangle

groups = {
    "gas (+CHP)": ["gas CHP", "OCGT", "CCGT"],
    "heat pump": ["heat pump"],
    "resistive heater": ["resistive heater"],
    "biomass (+ CHP)": ["biomass"],
    "coal (+ CHP)": ["coal"],
    "oil (+ CHP)": ["oil"],
    "waste CHP": ["waste"],
    "solar": ["solar"],
    "offwind": ["offwind"],
}


def aggregate_by_keywords(df, groups):
    """
    Aggregate rows in df according to keyword groups.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with row index as technology names.
    groups : dict
        Keys = new aggregated name,
        Values = list of substrings to match in the index.

    Returns
    -------
    pd.DataFrame
    """
    df_out = df.copy()
    for new_name, keywords in groups.items():
        mask = df_out.index.to_series().str.contains("|".join(keywords))
        if mask.any():
            summed = df_out.loc[mask].sum()
            df_out = df_out.drop(df_out.index[mask])
            df_out.loc[new_name] = summed
    return df_out


def plot_capacity_comparison(
    df,
    scenarios=("AriadneDemand", "LowDemand"),
    tech_colors=None,
    plot_diff=False,
    title="Electricity capacities",
    ylabel="GW",
    save_path=None,
    figsize=(12, 6),
    hatch_second="//",  # hatch for second scenario in non-diff plot
    show_dummy_legend=True,  # legend with empty + hatched boxes
):
    """
    Plot electricity capacities by carrier for multiple scenarios, or their difference.

    Parameters
    ----------
    df : pd.DataFrame            # index: carriers, columns: scenario names
    scenarios : tuple(str, str)  # (base, compare)
    tech_colors : dict           # carrier -> color
    plot_diff : bool             # if True, plot (compare - base) as signed bars
    """
    if tech_colors is None:
        tech_colors = {}

    base_name, cmp_name = scenarios
    df = df.copy()

    # --- Establish supply/demand membership AND ORDER from base scenario only ---
    base = df[base_name]
    supply_order = base[base >= 0].sort_values(ascending=False).index
    demand_order = (
        base[base < 0].abs().sort_values(ascending=False).index
    )  # sort by |capacity|

    # Build plotting table
    if plot_diff:
        # signed difference (base - compare)
        diff = (df[base_name].abs() - df[cmp_name].abs()).rename("Diff").to_frame()
        df_plot_supply = diff.loc[supply_order.intersection(diff.index)]
        df_plot_demand = diff.loc[demand_order.intersection(diff.index)]
        plot_columns = ["Diff"]
    else:
        # normal comparison: keep both scenarios
        df_plot_supply = df[[base_name, cmp_name]].loc[
            supply_order.intersection(df.index)
        ]
        df_plot_demand = df[[base_name, cmp_name]].loc[
            demand_order.intersection(df.index)
        ]
        # For the demand side, show magnitudes (upward bars). Clip to avoid negative heights.
        df_plot_demand = df_plot_demand.apply(lambda s: (-s).clip(lower=0))
        plot_columns = [base_name, cmp_name]

    # Positions
    n_supply = len(df_plot_supply)
    n_demand = len(df_plot_demand)
    pos_supply = np.arange(n_supply)
    pos_demand = np.arange(n_demand) + n_supply + 1  # +1 gap
    bar_width = 0.35

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Supply bars
    for j, col in enumerate(plot_columns):
        if n_supply:
            vals = df_plot_supply[col].values
            ax.bar(
                pos_supply + (j - 0.5) * bar_width
                if len(plot_columns) > 1
                else pos_supply,
                vals,
                width=bar_width if len(plot_columns) > 1 else 0.8,
                color=[tech_colors.get(t, "#1f77b4") for t in df_plot_supply.index],
                hatch=(
                    ""
                    if (not plot_diff and col == base_name)
                    else (hatch_second if not plot_diff else "")
                ),
                alpha=0.9,
                linewidth=0.0,
            )

    # Demand bars
    for j, col in enumerate(plot_columns):
        if n_demand:
            vals = (
                df_plot_demand[col].values
                if not plot_diff
                else (df_plot_demand[col].values)
            )
            # In diff mode we plot signed values; in normal mode they are already positive magnitudes
            ax.bar(
                pos_demand + (j - 0.5) * bar_width
                if len(plot_columns) > 1
                else pos_demand,
                vals,
                width=bar_width if len(plot_columns) > 1 else 0.8,
                color=[tech_colors.get(t, "#1f77b4") for t in df_plot_demand.index],
                hatch=(
                    ""
                    if (not plot_diff and col == base_name)
                    else (hatch_second if not plot_diff else "")
                ),
                alpha=0.9,
                linewidth=0.0,
            )

    # Divider line between supply and demand
    if n_supply:
        ax.axvline(x=n_supply - 0.5, color="black", linestyle="--", linewidth=1.0)

    # Baseline for signed diffs
    if plot_diff:
        ax.axhline(0, color="black", linewidth=1.0)

    # X labels
    tech_labels = list(df_plot_supply.index) + list(df_plot_demand.index)
    all_positions = np.concatenate([pos_supply, pos_demand])
    ax.set_xticks(all_positions)
    ax.set_xticklabels(tech_labels, rotation=45, ha="right")

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)

    # Legend
    if plot_diff:
        handles = [
            mpatches.Patch(
                facecolor="white",
                edgecolor="black",
                label=f"{base_name} - {cmp_name}",
            )
        ]
        ax.legend(handles=handles, title="Difference")
    elif show_dummy_legend:
        handles = [
            mpatches.Patch(
                facecolor="white", edgecolor="black", hatch="", label=base_name
            ),
            mpatches.Patch(
                facecolor="white", edgecolor="black", hatch=hatch_second, label=cmp_name
            ),
        ]
        ax.legend(handles=handles, title="Scenario")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "regret_plots_lt",
        )

    configure_logging(snakemake)
    config = snakemake.config
    planning_horizons = snakemake.params.planning_horizons
    scenarios = ["AriadneDemand", "LowDemand"]
    tech_colors = snakemake.params.plotting["tech_colors"]

    # Load networks
    networks = defaultdict(dict)

    for fn in snakemake.input.networks:
        scenario = fn.split(os.sep)[-3]
        year = int(re.search(r"_(\d{4})\.nc$", fn).group(1))
        networks[scenario][year] = pypsa.Network(fn)

    # Load variables
    vars_dict = {}

    for fn in snakemake.input.regret_variables:
        df = (
            pd.read_excel(
                fn,
                index_col=list(range(5)),
                # index_col=["Model", "Scenario", "Region", "Variable", "Unit"],
                sheet_name="data",
            )
            .groupby(["Variable", "Unit"], dropna=False)
            .sum()
        ).round(5)

        if "AriadneDemand" in fn:
            vars_dict["AriadneDemand"] = df
        elif "LowDemand" in fn:
            vars_dict["LowDemand"] = df

    # ensure output directory exist
    if not os.path.exists(snakemake.output[-1]):
        os.makedirs(snakemake.output[-1])

    # Capacity plot DE

    tech_colors["gas (+CHP)"] = tech_colors["OCGT"]
    tech_colors["biomass (+ CHP)"] = tech_colors["biomass"]
    tech_colors["coal (+ CHP)"] = tech_colors["coal"]
    tech_colors["oil (+ CHP)"] = tech_colors["oil"]
    tech_colors["heat pump"] = tech_colors["heat pump"]
    tech_colors["waste CHP"] = tech_colors["waste"]

    kwargs = {
        "groupby": ["bus", "carrier"],
        "at_port": True,
        "nice_names": False,
    }

    capa_comp = pd.DataFrame(columns=scenarios)

    for year in planning_horizons:
        for scenario in scenarios:
            n = networks[scenario][year]

            capacities = (
                n.statistics.optimal_capacity(
                    bus_carrier=["AC", "low voltage"],
                    **kwargs,
                )
                .filter(like="DE")
                .groupby("carrier")
                .sum()
                .drop(
                    ["AC", "DC", "electricity distribution grid"],
                    errors="ignore",
                )
                .multiply(1e-3)  # MW → GW
            )
            capa_comp[scenario] = capacities

        capa_comp_agg = aggregate_by_keywords(capa_comp, groups)
        # drop capa with less than 100 MW in both scenarios
        capa_comp_agg = capa_comp_agg[(capa_comp_agg.abs() >= 0.1).any(axis=1)]

        plot_capacity_comparison(
            df=capa_comp_agg,
            scenarios=["AriadneDemand", "LowDemand"],
            tech_colors=tech_colors,
            plot_diff=False,
            title=f"Electricity capacities in DE: {year}",
            save_path=snakemake.output.elec_capa_comp_de_2025
            if year == 2025
            else snakemake.output[-1] + f"/elec_capa_comp_de_{year}.png",
        )

        plot_capacity_comparison(
            df=capa_comp_agg,
            scenarios=["AriadneDemand", "LowDemand"],
            tech_colors=tech_colors,
            plot_diff=True,
            title=f"Difference of electricity capacities in DE: {year}",
            save_path=snakemake.output[-1] + f"/elec_capa_diff_de_{year}.png",
        )

    # Capacity plot outside DE

    for year in planning_horizons:
        for scenario in scenarios:
            n = networks[scenario][year]

            capacities = (
                n.statistics.optimal_capacity(
                    bus_carrier=["AC", "low voltage"],
                    **kwargs,
                )
                .loc[
                    lambda df: ~df.index.get_level_values("bus").str.contains(
                        "DE", regex=False
                    )
                ]
                .groupby("carrier")
                .sum()
                .drop(["AC", "DC", "electricity distribution grid"], errors="ignore")
                .multiply(1e-3)  # MW → GW
            )
            capa_comp[scenario] = capacities

        capa_comp_agg = aggregate_by_keywords(capa_comp, groups)
        # drop capa with less than 100 MW in both scenarios
        capa_comp_agg = capa_comp_agg[(capa_comp_agg.abs() >= 0.1).any(axis=1)]

        plot_capacity_comparison(
            df=capa_comp_agg,
            scenarios=["AriadneDemand", "LowDemand"],
            tech_colors=tech_colors,
            plot_diff=False,
            title=f"Electricity capacities in EU (outside DE): {year}",
            save_path=snakemake.output[-1] + f"/elec_capa_comp_eu_{year}.png",
        )

        plot_capacity_comparison(
            df=capa_comp_agg,
            scenarios=["AriadneDemand", "LowDemand"],
            tech_colors=tech_colors,
            plot_diff=True,
            title=f"Difference of electricity capacities in EU (outside DE): {year}",
            save_path=snakemake.output[-1] + f"/elec_capa_diff_eu_{year}.png",
        )

    # Electricity demand as bar plot

    demand_comp = pd.DataFrame(columns=scenarios)

    for year in planning_horizons:
        for scenario in scenarios:
            n = networks[scenario][year]

            electricity_withdrawal = (
                n.statistics.withdrawal(bus_carrier=["low voltage", "AC"], **kwargs)
                .filter(like="DE")
                .groupby(["carrier"])
                .sum()
                .multiply(1e-6)  # MWh → TWh
            )

            demand_comp[scenario] = electricity_withdrawal

        demand_comp_agg = aggregate_by_keywords(demand_comp, groups)
        # drop capa with less than 100 MW in both scenarios
        demand_comp_agg = demand_comp_agg[(demand_comp_agg.abs() >= 0.1).any(axis=1)]

        plot_capacity_comparison(
            df=demand_comp_agg,
            scenarios=["AriadneDemand", "LowDemand"],
            tech_colors=tech_colors,
            plot_diff=False,
            title=f"Electricity demand in DE: {year}",
            ylabel="TWh",
            save_path=snakemake.output[-1] + f"/elec_demand_comp_de_{year}.png",
        )

        plot_capacity_comparison(
            df=demand_comp_agg,
            scenarios=["AriadneDemand", "LowDemand"],
            tech_colors=tech_colors,
            plot_diff=True,
            title=f"Difference of electricity demand in DE: {year}",
            ylabel="TWh",
            save_path=snakemake.output[-1] + f"/elec_demand_diff_de_{year}.png",
        )

        # ToDo

        # Electricity demand temporal (+ import)
        # Split electricity demand of distribution grid

        # System cost CAPEX
        capex_data = {}

        for scenario in scenarios:
            # Extract CAPEX data
            df = vars_dict[scenario]
            capex = df[
                df.index.get_level_values("Variable").str.startswith(
                    "System Cost|CAPEX"
                )
            ]
            capex_top = capex[
                capex.index.get_level_values("Variable").str.count(r"\|") == 2
            ]

            # Reset index and prepare data
            capex_reset = capex_top.reset_index()
            capex_reset["Category"] = capex_reset["Variable"].str.split("|").str[2]

            # Get year columns on first iteration
            if scenario == scenarios[0]:
                year_cols = [
                    col
                    for col in capex_reset.columns
                    if str(col) in ["2025", "2030", "2035"]
                ]

            # Store processed data
            capex_data[scenario] = capex_reset.set_index("Category")[year_cols].fillna(
                0
            )

        # Extract for easier access
        capex_low_plot = capex_data["LowDemand"]
        capex_ariadne_plot = capex_data["AriadneDemand"]

        # Set up the plot
        fig, ax = plt.subplots(figsize=(12, 8))

        # Define years and categories
        years = [str(col) for col in year_cols]  # Convert to strings for labels
        categories = capex_low_plot.index.tolist()

        # Define colors for each category using tech_colors
        category_color_map = {
            "Electricity": "gold",
            "Gases": tech_colors["gas"],
            "Heat": tech_colors["heat"],
            "Hydrogen": tech_colors["H2"],
            "Liquids": tech_colors["oil"],
            "Methanol": tech_colors["methanol"],
        }

        # Create colors list based on categories
        colors = [category_color_map.get(cat, "#808080") for cat in categories]

        # Set up bar positions
        x = np.arange(len(years))
        width = 0.35

        # Create stacked bars
        bottom_low = np.zeros(len(years))
        bottom_ariadne = np.zeros(len(years))

        for i, category in enumerate(categories):
            # Get values for each year
            low_values = [capex_low_plot.loc[category, col] for col in year_cols]
            ariadne_values = [
                capex_ariadne_plot.loc[category, col] for col in year_cols
            ]

            # Plot bars - LowDemand with hatching, AriadneDemand without
            ax.bar(
                x - width / 2,
                low_values,
                width,
                bottom=bottom_low,
                label=category,
                color=colors[i],
                alpha=0.8,
                hatch="///",
            )
            ax.bar(
                x + width / 2,
                ariadne_values,
                width,
                bottom=bottom_ariadne,
                color=colors[i],
                alpha=0.8,
            )

            # Add category values inside bars (only if value > threshold for readability)
            for j in range(len(years)):
                # LowDemand category values
                if low_values[j] > 0.5:  # Only show if value is significant
                    ax.text(
                        j - width / 2,
                        bottom_low[j] + low_values[j] / 2,
                        f"{low_values[j]:.1f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        fontweight="bold",
                        color="white",
                    )

                # AriadneDemand category values
                if ariadne_values[j] > 0.5:  # Only show if value is significant
                    ax.text(
                        j + width / 2,
                        bottom_ariadne[j] + ariadne_values[j] / 2,
                        f"{ariadne_values[j]:.1f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        fontweight="bold",
                        color="white",
                    )

            # Update bottom for stacking
            bottom_low += low_values
            bottom_ariadne += ariadne_values

        # Customize the plot
        ax.set_xlabel("Year", fontsize=12)
        ax.set_ylabel("billion €", fontsize=12)
        ax.set_title(
            "System Cost CAPEX Comparison: LowDemand vs AriadneDemand",
            fontsize=14,
            fontweight="bold",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(years)

        # Category legend elements
        category_elements = [
            Rectangle((0, 0), 1, 1, facecolor=colors[i], alpha=0.8, label=cat)
            for i, cat in enumerate(categories)
        ]

        # Scenario legend elements
        scenario_elements = [
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor="gray",
                alpha=0.8,
                hatch="///",
                label="LowDemand",
            ),
            Rectangle((0, 0), 1, 1, facecolor="gray", alpha=0.8, label="AriadneDemand"),
        ]

        # Create separate legends
        category_legend = ax.legend(
            handles=category_elements,
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            title="Categories",
        )
        scenario_legend = ax.legend(
            handles=scenario_elements,
            loc="upper left",
            bbox_to_anchor=(1.05, 0.6),
            title="Scenarios",
        )

        # Add both legends to the plot
        ax.add_artist(category_legend)

        # Add total CAPEX values on top of bars
        for i, year in enumerate(years):
            ax.text(
                i - width / 2,
                bottom_low[i] + max(max(bottom_low), max(bottom_ariadne)) * 0.02,
                f"{bottom_low[i]:.1f}",
                ha="center",
                fontweight="bold",
                fontsize=9,
            )
            ax.text(
                i + width / 2,
                bottom_ariadne[i] + max(max(bottom_low), max(bottom_ariadne)) * 0.02,
                f"{bottom_ariadne[i]:.1f}",
                ha="center",
                fontweight="bold",
                fontsize=9,
            )

        # Add grid
        ax.grid(True, alpha=0.3, axis="y")
        ax.set_axisbelow(True)

        # Adjust y-axis limits to accommodate top labels
        y_max = max(max(bottom_low), max(bottom_ariadne))
        ax.set_ylim(0, y_max * 1.1)
        plt.savefig(snakemake.output[-1] + "/capex_de.png", bbox_inches="tight")
        plt.close()
