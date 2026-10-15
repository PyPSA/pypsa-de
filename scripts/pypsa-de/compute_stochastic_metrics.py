# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Rolls up the stochastic-grid-scenario analysis into the WS/SP/EEV/EVPI/ECIU
comparison from notebooks/stochastic_grid_uncertainty.ipynb, using only the
already-solved topology and dispatch-only evaluation networks (no re-solving
here):

- WS  (wait-and-see): expected cost if the grid outcome were known before
  investing = sum_s probability_s * total_cost(solved topology-s)
- SP  (stochastic solution): the stochastic portfolio re-dispatched on each real
  grid (cost matrix "stochastic" row = capex + correct-physics expected opex),
  NOT total_cost(topology-stochastic): the native scenario-dimension objective
  prices dispatch with phantom-cycle physics (per-scenario line reactance is not
  honoured), so it is not comparable to the deterministic WS/EEV below
- EEV (expected-value planner): capex(topology-eev) + expected opex of
  dispatching the eev portfolio under each real topology
- EVPI = SP - WS   (value of perfect information about the grid)
- ECIU = EEV - SP  (value of stochastic planning over naive EV planning)

`total_cost(n) = n.objective + n.objective_constant`, not `n.objective` alone
(see `compute_metrics`'s docstring for why this matters).

The same five metrics are also computed for Germany only (`*_de` outputs). A
whole-network objective is not decomposable per country, so the DE cost reuses
the per-country accounting from `export_ariadne_variables`: `get_tsc` (capex +
opex of DE-located assets, cross-border transmission split 50/50) plus the net
trade cost (imports minus exports valued at domestic marginal prices). To stay
comparable with the network-wide metrics - which treat the grid as exogenous and
omit its capex - DE transmission capex is likewise excluded (`de_total_cost`).
"""

import logging
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from pypsa.descriptors import nominal_attrs
from pypsa.statistics import get_transmission_carriers

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from export_ariadne_variables import get_trade_cost, get_tsc

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

REGION = "DE"

# Carrier groups traded across the DE border, as in export_ariadne_variables's
# get_economy (AC lines and DC links must each be queried on their own).
TRADE_CARRIER_GROUPS = [
    ["AC"],
    ["DC"],
    ["renewable oil", "renewable gas", "methanol"],
    ["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"],
]


def get_one_scenario(n: pypsa.Network) -> pypsa.Network:
    """One scenario's slice of `n`, as a plain (non-scenario) network.

    `n.get_scenario()` assumes every component's static index carries a
    `(scenario, name)` MultiIndex, which doesn't hold for `sub_networks`
    (topology metadata repopulated by `n.determine_network_topology()`, not
    itself scenario-varying) and raises `TypeError: Index must be a
    MultiIndex` (or, once emptied, `KeyError` from the same assumption in
    get_scenario()'s empty-component branch). Not needed for capacity/capex
    aggregation below, so it's replaced with a properly (empty)
    MultiIndexed frame on a copy first, to route around the bug entirely.
    """
    n = n.copy()
    empty_mi = pd.MultiIndex.from_arrays([[], []], names=["scenario", "name"])
    n.c.sub_networks.static = n.c.sub_networks.static.reindex(empty_mi)
    return n.get_scenario(n.scenarios[0])


def capacities_by_carrier(n: pypsa.Network) -> pd.Series:
    """Installed capacity by carrier, taking one scenario slice if `n` has scenarios."""
    if n.has_scenarios:
        n = get_one_scenario(n)
    caps = n.statistics.optimal_capacity()
    return caps.groupby("carrier").sum()


def de_capacities_by_carrier(n: pypsa.Network, region: str = REGION) -> pd.Series:
    """Installed capacity by carrier for assets located in `region`."""
    if n.has_scenarios:
        n = get_one_scenario(n)
    caps = n.statistics.optimal_capacity(groupby=["bus", "carrier"])
    caps = caps[caps.index.get_level_values("bus").str.startswith(region)]
    return caps.groupby("carrier").sum()


def portfolio_capex(n: pypsa.Network) -> float:
    """
    Capex of a portfolio's capacities, counted once. Restricted to
    originally-extendable components (generators/storage/stores/non-grid
    links) - the fixed grid Lines/DC-Links (non-extendable, given as the
    scenario's topology) are excluded, matching the notebook's convention
    that grid capex is exogenous and not part of the compared portfolio.
    """
    if n.has_scenarios:
        n = get_one_scenario(n)
    total = 0.0
    for component, attr in nominal_attrs.items():
        static = n.components[component].static
        opt_attr = f"{attr}_opt"
        if not {"capital_cost", opt_attr, f"{attr}_extendable"}.issubset(static.columns):
            continue
        ext = static.query(f"{attr}_extendable")
        total += (ext["capital_cost"] * ext[opt_attr]).sum()
    return total


def de_trade_cost(n: pypsa.Network, region: str = REGION) -> float:
    """Net DE trade cost (imports minus exports at domestic prices); positive = net import cost."""
    return float(sum(get_trade_cost(n, region, c) for c in TRADE_CARRIER_GROUPS))


def de_total_cost(n: pypsa.Network, region: str = REGION) -> float:
    """DE total system cost: capex + opex of DE assets (grid capex excluded) plus net trade cost.

    Mirrors the network-wide convention of treating the (exogenous) grid as a
    given: the capex of transmission carriers is dropped, so WS/SP/EEV compare
    portfolios and dispatch, not the scenario's fixed grid buildout.
    """
    tsc = get_tsc(n, region)
    transmission = get_transmission_carriers(n).get_level_values("carrier")
    grid_capex = tsc.loc[tsc.index.isin(transmission), "capex"].sum()
    return tsc.values.sum() - grid_capex + de_trade_cost(n, region)


def build_capacities_table(
    topologies: dict[str, pypsa.Network], cap_func=capacities_by_carrier
) -> pd.DataFrame:
    return pd.DataFrame(
        {name: cap_func(n) for name, n in topologies.items()}
    ).fillna(0)


def build_cost_matrix(
    topologies: dict[str, pypsa.Network],
    evaluations: dict[tuple[str, str], pypsa.Network],
    grid_scenario_names: list[str],
    probabilities: dict[str, float],
) -> pd.DataFrame:
    """Rows = portfolios, columns = realized topologies (+ expected), cell = capex + opex."""
    rows = {}
    for portfolio, n_portfolio in topologies.items():
        capex = portfolio_capex(n_portfolio)
        row = {}
        for grid_scenario in grid_scenario_names:
            n_eval = evaluations[(portfolio, grid_scenario)]
            opex = n_eval.objective
            row[grid_scenario] = capex + opex
        rows[portfolio] = row
    cost_matrix = pd.DataFrame(rows).T
    cost_matrix["expected"] = sum(
        probabilities[s] * cost_matrix[s] for s in grid_scenario_names
    )
    return cost_matrix


def build_cost_matrix_de(
    evaluations: dict[tuple[str, str], pypsa.Network],
    grid_scenario_names: list[str],
    portfolios: list[str],
    probabilities: dict[str, float],
    region: str = REGION,
) -> pd.DataFrame:
    """DE cost matrix: each cell is `de_total_cost` of the portfolio re-dispatched on a realized grid."""
    rows = {}
    for portfolio in portfolios:
        rows[portfolio] = {
            s: de_total_cost(evaluations[(portfolio, s)], region)
            for s in grid_scenario_names
        }
    cost_matrix = pd.DataFrame(rows).T
    cost_matrix["expected"] = sum(
        probabilities[s] * cost_matrix[s] for s in grid_scenario_names
    )
    return cost_matrix


def total_cost(n: pypsa.Network) -> float:
    """True total system cost of a solved network: `n.objective + n.objective_constant`,
    not `n.objective` alone.

    Per PyPSA's own convention (`pypsa.Network.objective`'s docstring), `n.objective` is
    only the LP's *variable* part - `objective_constant` corrects for capital cost charged
    on the *full* `s_nom_opt`/`p_nom_opt` of extendable transmission-like assets whose
    optimization is bounded below by an already-installed baseline (`s_nom_min`/
    `p_nom_min`), so as not to double-charge that sunk baseline portion inside the
    objective PyPSA actually optimizes. `objective` alone therefore *omits* that baseline's
    annuitized cost - it isn't "excluded exogenous grid cost" the way `portfolio_capex()`
    deliberately excludes fixed (non-extendable) grid Lines/DC-Links; PyPSA's LP builder
    already never reads `capital_cost` for non-extendable components at all (true sunk
    cost, zero contribution either way), so that exclusion needs no special handling here.

    Confirmed empirically against a solved run: `objective_constant` was ~9.36bn EUR,
    identical across all four portfolios (so it isn't tracking each portfolio's own,
    genuinely different, fixed-grid capex) - and omitting it from SP alone very nearly
    reproduced the entire previously-reported ECIU, i.e. the "value of stochastic planning"
    metric was almost entirely this omitted baseline-capex term, not a real cost difference.
    """
    return n.objective + n.objective_constant


def compute_metrics(
    cost_matrix: pd.DataFrame,
    topologies: dict[str, pypsa.Network],
    grid_scenario_names: list[str],
    probabilities: dict[str, float],
    total_cost_fn=total_cost,
) -> pd.Series:
    ws = sum(
        probabilities[s] * total_cost_fn(topologies[s]) for s in grid_scenario_names
    )
    # SP from the cost matrix (stochastic portfolio re-dispatched on each real grid,
    # correct physics), on the same basis as EEV - not total_cost(topologies[
    # "stochastic"]), whose native scenario-dimension objective prices dispatch with
    # phantom-cycle physics and can sit above EEV. The stochastic *plan* is still
    # sized by that degraded solve, so a residual ECIU < 0 is a genuine signal (the
    # plan underperforms the naive one), not a measurement artefact.
    sp = cost_matrix.loc["stochastic", "expected"]
    eev = cost_matrix.loc["eev", "expected"]
    evpi = sp - ws
    eciu = eev - sp

    tol = 1e-3
    negative = {
        name: val
        for name, val in {"EVPI (SP-WS)": evpi, "ECIU (EEV-SP)": eciu}.items()
        if val < -tol
    }
    if negative:
        detail = ", ".join(f"{k} = {v / 1e9:+.3f} bn EUR" for k, v in negative.items())
        logger.warning(
            f"Negative stochastic metric(s), expected >= 0: {detail} "
            f"(WS={ws / 1e9:.3f}, SP={sp / 1e9:.3f}, EEV={eev / 1e9:.3f} bn EUR). "
            "A negative ECIU means the stochastic plan underperforms the naive EEV "
            "plan - expected when the phantom-cycle physics of the native scenario "
            "solve dominate (see stochastic_grid_uncertainty.ipynb). For the DE-only "
            "metrics a subset cost need not satisfy the network-wide relations."
        )

    return pd.Series(
        {"WS": ws, "SP": sp, "EEV": eev, "EVPI": evpi, "ECIU": eciu},
        name="value",
    )


def plot_metrics(metrics: pd.Series, path: str, title: str) -> None:
    """Expected system cost (WS/SP/EEV) and information/stochastic value (EVPI/ECIU)."""
    m = metrics / 1e9
    fig, axd = plt.subplot_mosaic(
        [["levels", "deltas"]],
        figsize=(9, 4),
        layout="constrained",
        width_ratios=[3, 2],
    )
    axd["levels"].bar(
        ["WS", "SP", "EEV"],
        [m["WS"], m["SP"], m["EEV"]],
        color=["#4477aa", "#228833", "#ee6677"],
    )
    axd["levels"].set_ylabel("expected system cost [bn EUR]")
    axd["levels"].set_title("Expected cost")
    axd["deltas"].bar(
        ["EVPI\n(SP-WS)", "ECIU\n(EEV-SP)"],
        [m["EVPI"], m["ECIU"]],
        color=["#66ccee", "#aa3377"],
    )
    axd["deltas"].axhline(0, color="k", lw=0.8)
    axd["deltas"].set_ylabel("value [bn EUR]")
    axd["deltas"].set_title("Value of information")
    fig.suptitle(title)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_cost_matrix(cost_matrix: pd.DataFrame, path: str, title: str) -> None:
    """Heatmap of total cost [bn EUR]: rows = portfolios, columns = realized grid."""
    data = cost_matrix / 1e9
    fig, ax = plt.subplots(
        figsize=(1.8 + 1.1 * data.shape[1], 1.5 + 0.55 * data.shape[0]),
        layout="constrained",
    )
    im = ax.imshow(data.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(data.shape[1]), labels=data.columns, rotation=45, ha="right")
    ax.set_yticks(range(data.shape[0]), labels=data.index)
    ax.set_xlabel("realized grid")
    ax.set_ylabel("portfolio")
    vmid = (data.values.max() + data.values.min()) / 2
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data.values[i, j]
            ax.text(
                j, i, f"{val:.1f}",
                ha="center", va="center", fontsize=8,
                color="white" if val < vmid else "black",
            )
    fig.colorbar(im, ax=ax, label="total cost [bn EUR]")
    ax.set_title(title)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_capacities(capacities: pd.DataFrame, path: str, title: str) -> None:
    """Table of installed capacity [GW] by carrier (rows) and portfolio (columns).

    All-zero carriers are dropped from the figure for legibility (the CSV keeps them).
    """
    data = (capacities / 1e3).round(1)
    data = data.loc[(data != 0).any(axis=1)]
    fig, ax = plt.subplots(
        figsize=(1.5 + 1.3 * data.shape[1], 0.5 + 0.3 * data.shape[0]),
        layout="constrained",
    )
    ax.axis("off")
    tbl = ax.table(
        cellText=data.values,
        rowLabels=data.index,
        colLabels=data.columns,
        loc="center",
        cellLoc="right",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    ax.set_title(f"{title} - installed capacity [GW]")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "compute_stochastic_metrics",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
        )

    configure_logging(snakemake)

    grid_scenario_names = snakemake.params.grid_scenario_names
    grid_scenario_values = snakemake.params.grid_scenario_values
    probabilities = {
        name: cfg["probability"] for name, cfg in snakemake.params.scenarios.items()
    }

    # unpack()-flattened named inputs: "topology_<grid_scenario>" and
    # "eval_<portfolio>__<grid_scenario>" (see rules/pypsa-de/stochastic_grid.smk).
    topologies = {}
    evaluations = {}
    for key, path in snakemake.input.items():
        if key.startswith("topology_"):
            topologies[key[len("topology_") :]] = pypsa.Network(path)
        elif key.startswith("eval_"):
            portfolio, grid_scenario = key[len("eval_") :].split("__", 1)
            evaluations[(portfolio, grid_scenario)] = pypsa.Network(path)

    title = f"Stochastic grid {snakemake.wildcards.planning_horizons}"
    portfolios = list(topologies)

    # Network-wide
    capacities = build_capacities_table(topologies)
    capacities.to_csv(snakemake.output.capacities)
    plot_capacities(capacities, snakemake.output.capacities_png, title)

    cost_matrix = build_cost_matrix(
        topologies, evaluations, grid_scenario_names, probabilities
    )
    cost_matrix.to_csv(snakemake.output.cost_matrix)
    plot_cost_matrix(cost_matrix, snakemake.output.cost_matrix_png, title)

    metrics = compute_metrics(
        cost_matrix, topologies, grid_scenario_names, probabilities
    )
    metrics.to_csv(snakemake.output.metrics)
    plot_metrics(metrics, snakemake.output.metrics_png, title)
    logger.info(f"Stochastic-grid metrics (network):\n{metrics}")

    # Germany only
    capacities_de = build_capacities_table(topologies, de_capacities_by_carrier)
    capacities_de.to_csv(snakemake.output.capacities_de)
    plot_capacities(capacities_de, snakemake.output.capacities_de_png, f"{title} (DE)")

    cost_matrix_de = build_cost_matrix_de(
        evaluations, grid_scenario_names, portfolios, probabilities
    )
    cost_matrix_de.to_csv(snakemake.output.cost_matrix_de)
    plot_cost_matrix(cost_matrix_de, snakemake.output.cost_matrix_de_png, f"{title} (DE)")

    metrics_de = compute_metrics(
        cost_matrix_de, topologies, grid_scenario_names, probabilities, de_total_cost
    )
    metrics_de.to_csv(snakemake.output.metrics_de)
    plot_metrics(metrics_de, snakemake.output.metrics_de_png, f"{title} (DE)")
    logger.info(f"Stochastic-grid metrics (DE):\n{metrics_de}")
