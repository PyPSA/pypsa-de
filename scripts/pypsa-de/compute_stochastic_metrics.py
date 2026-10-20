# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Rolls up the stochastic-grid-scenario analysis into the WS/SP/EEV/EVPI/ECIU
comparison from notebooks/stochastic_grid_uncertainty.ipynb, using only the
already-solved topology and dispatch-only evaluation networks (no re-solving
here):

- WS  (wait-and-see): expected cost if the grid outcome were known before
  investing = sum_s probability_s * cost_matrix[s, s], the matrix diagonal -
  portfolio s (the plan optimised for grid s) evaluated on grid s
- SP  (stochastic solution): the stochastic portfolio re-dispatched on each real
  grid (cost matrix "stochastic" row = capex + correct-physics expected opex),
  NOT the native scenario-dimension objective of topology-stochastic, which
  prices dispatch with phantom-cycle physics (per-scenario line reactance is not
  honoured) and is not comparable to the deterministic WS/EEV
- EEV (expected-value planner): capex(topology-eev) + expected opex of
  dispatching the eev portfolio under each real topology
- EVPI = SP - WS   (value of perfect information about the grid)
- ECIU = EEV - SP  (value of stochastic planning over naive EV planning)

WS, SP and EEV all read from the one objective-basis cost matrix - portfolio
capex (new investment, counted once) plus dispatch opex (`n.objective` of the
re-dispatch). They therefore share a basis and EVPI/ECIU are differences on equal
footing. The sunk baseline capex (`objective_constant`) is excluded from every
term alike, so it cannot leak into EVPI/ECIU.

The same five metrics are also computed for Germany only (`*_de` outputs). A
whole-network objective is not decomposable per country, so the DE cost reuses
the per-country accounting from `export_ariadne_variables`: `get_tsc` (capex +
opex of DE-located assets, cross-border transmission split 50/50) plus the net
trade cost (imports minus exports valued at domestic marginal prices). To stay
comparable with the network-wide metrics - which treat the grid as exogenous and
omit its capex - DE transmission capex is likewise excluded (`de_total_cost`).

Two caveats when comparing these numbers to other outputs:
- Basis: the network-wide WS/SP/EEV use the objective-basis cost matrix
  (portfolio capex + dispatch opex, sunk `objective_constant` excluded). The `*_de` metrics and the
  `total_energy_system_cost_by_scenario.png` plot use the statistics basis
  (`capex + opex` of all assets, incl. already-built capacity). The two can
  differ substantially, so a network-wide metric is NOT a DE metric scaled up.
- Scope: the `*_de` metrics are Germany only; the network-wide metrics cover the
  whole modelled system. The DE share and its complement reconstruct the
  whole-network statistics total (DE + outside-DE); the country-by-country
  build-up is shown in notebooks/system_cost_decomposition.ipynb.
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


def foreign_portfolio_capex(n: pypsa.Network, region: str = REGION) -> float:
    """Capex of the originally-extendable, non-grid capacity located outside `region`.

    This is the set `freeze_outside_de_capacities` fixes to a source network
    (generators/storage/stores and non-DC links with a foreign bus). For the
    freeze source it is its own optimised foreign buildout; added to every other
    (frozen) portfolio's cost-matrix row it accounts for that same buildout -
    which the freeze turned non-extendable, so `portfolio_capex` silently drops
    it - making the objective-basis matrix comparable across rows. Sunk/brownfield
    foreign capacity (non-extendable in the source too) and the exogenous grid
    (Lines/Transformers, DC Links) stay excluded, exactly as in `portfolio_capex`.
    """
    if n.has_scenarios:
        n = get_one_scenario(n)
    country = n.buses.country
    foreign = set(country[(country != "") & (country != region)].unique())
    total = 0.0
    for component, attr in nominal_attrs.items():
        if component in ("Line", "Transformer"):
            continue
        static = n.components[component].static
        opt_attr = f"{attr}_opt"
        if static.empty or not {"capital_cost", opt_attr, f"{attr}_extendable"}.issubset(static.columns):
            continue
        if {"bus0", "bus1"}.issubset(static.columns):
            c0, c1 = static.bus0.map(country), static.bus1.map(country)
            mask = ~c0.eq(region) & ~c1.eq(region) & (c0.isin(foreign) | c1.isin(foreign))
            if component == "Link":
                mask &= static.carrier != "DC"
        elif "bus" in static.columns:
            mask = static.bus.map(country).isin(foreign)
        else:
            continue
        sel = static[mask & static[f"{attr}_extendable"]]
        total += (sel["capital_cost"] * sel[opt_attr]).sum()
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


def de_total_cost_with_grid(n: pypsa.Network, region: str = REGION) -> float:
    """Like `de_total_cost` but KEEPING the DE transmission-grid capex, so each
    realized topology carries the cost of its own grid buildout."""
    return get_tsc(n, region).values.sum() + de_trade_cost(n, region)


def total_cost_stats(n: pypsa.Network) -> float:
    """Whole-network system cost on the statistics basis: capex + opex of ALL
    assets, incl. already-built capacity and the grid. Differs from the objective
    basis used by the cost matrix, which omits the sunk-capacity baseline
    (`objective_constant`) and the exogenous grid capex."""
    return n.statistics.capex().sum() + n.statistics.opex().sum()


def outside_de_non_trade_cost(n: pypsa.Network, region: str = REGION) -> float:
    """Capex+opex of everything outside `region`: the whole-network statistics
    total minus the DE share (`get_tsc`). This is the exact complement, so
    ``get_tsc(region) + outside_de_non_trade_cost == total_cost_stats``. It covers
    every other country plus the shared EU-level primary-energy buses that
    `get_tsc` attributes to no country. Non-trade (trade nets out internally).

    The country-by-country build-up of this same quantity, and that it sums to
    the whole-network total, is demonstrated in
    notebooks/system_cost_decomposition.ipynb.
    """
    return total_cost_stats(n) - get_tsc(n, region).values.sum()


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
    extra_capex: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Rows = portfolios, columns = realized topologies (+ expected), cell = capex + opex.

    `extra_capex` adds a per-portfolio capex offset (the frozen foreign buildout
    under `freeze_outside_de_capacities`, see `foreign_portfolio_capex`); empty
    when no freeze applies.
    """
    extra_capex = extra_capex or {}
    rows = {}
    for portfolio, n_portfolio in topologies.items():
        capex = portfolio_capex(n_portfolio) + extra_capex.get(portfolio, 0.0)
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


# Wait-and-see optimality tolerance (relative). The property is exact in theory,
# so any deviation is solver-gap noise (observed <~0.1%); the margin stays well
# below gross accounting errors such as a freeze source whose foreign capex is
# counted while every other portfolio's is dropped (see portfolio_capex).
WAIT_AND_SEE_RTOL = 1e-2


def assert_wait_and_see_optimality(
    cost_matrix: pd.DataFrame,
    grid_scenario_names: list[str],
    rtol: float = WAIT_AND_SEE_RTOL,
) -> None:
    """Assert the objective-basis cost matrix obeys wait-and-see optimality.

    Each grid scenario's own portfolio must be the cheapest cell in its realized
    grid column (no other portfolio beats the one optimised for that grid), and
    the stochastic portfolio must be the cheapest cell in the `expected` column.
    Only `build_cost_matrix`'s objective-basis matrix satisfies this exactly; the
    statistics-basis (`cost_matrix_total`/`_outside_de`) and DE matrices use a
    different cost definition and are deliberately not checked.

    Raises
    ------
    AssertionError
        If any named portfolio exceeds its column minimum by more than `rtol`.
    """
    checks = [(s, s) for s in grid_scenario_names]
    if "stochastic" in cost_matrix.index:
        checks.append(("stochastic", "expected"))
    violations = []
    for portfolio, column in checks:
        col = cost_matrix[column]
        best = col.min()
        cost = col[portfolio]
        if cost > best * (1 + rtol):
            cheapest = col.idxmin()
            violations.append(
                f"column '{column}': expected '{portfolio}' cheapest, but "
                f"'{cheapest}' is cheaper ({best / 1e9:.3f} < {cost / 1e9:.3f} bn "
                f"EUR, {portfolio} is {cost / best - 1:.2%} above the minimum)"
            )
    assert not violations, (
        "Wait-and-see optimality violated in the objective-basis cost matrix - a "
        "portfolio is beaten on the grid it was optimised for (or the stochastic "
        "plan is not expected-cheapest). This signals an inconsistent cost basis "
        "across rows, e.g. freeze_outside_de_capacities dropping frozen foreign "
        "capex from portfolio_capex for every portfolio except the freeze source. "
        "Details:\n  " + "\n  ".join(violations)
    )


def build_eval_cost_matrix(
    evaluations: dict[tuple[str, str], pypsa.Network],
    grid_scenario_names: list[str],
    portfolios: list[str],
    probabilities: dict[str, float],
    cost_fn,
) -> pd.DataFrame:
    """Cost matrix (rows = portfolios, columns = realized grids + `expected`)
    where each cell is `cost_fn(n)` of the portfolio re-dispatched on a realized
    grid (the `_st` evaluation network, i.e. the ST realization of that grid)."""
    rows = {
        portfolio: {
            s: cost_fn(evaluations[(portfolio, s)]) for s in grid_scenario_names
        }
        for portfolio in portfolios
    }
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
    return build_eval_cost_matrix(
        evaluations,
        grid_scenario_names,
        portfolios,
        probabilities,
        lambda n: de_total_cost(n, region),
    )


def compute_metrics(
    cost_matrix: pd.DataFrame,
    grid_scenario_names: list[str],
    probabilities: dict[str, float],
) -> pd.Series:
    # All three levels read from the one objective-basis cost matrix, so EVPI/ECIU
    # are equal-basis differences. WS (wait-and-see, grid known before investing):
    # each scenario gets its own optimal plan, i.e. the diagonal cost_matrix[s, s]
    # (portfolio s on grid s). SP: the stochastic portfolio re-dispatched with
    # correct physics (the "stochastic" row's expected column) - not the native
    # scenario-dimension objective, whose phantom-cycle physics can sit above EEV;
    # the stochastic *plan* is still sized by that degraded solve, so a residual
    # ECIU < 0 is a genuine signal (plan underperforms the naive one), not an
    # artefact. EEV: the expected-value planner's "eev" row.
    ws = sum(probabilities[s] * cost_matrix.loc[s, s] for s in grid_scenario_names)
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

    # Under freeze_outside_de_capacities=<grid_scenario>, that source portfolio is
    # the only one left with its foreign non-grid capacity extendable; every other
    # portfolio has it frozen non-extendable, so portfolio_capex drops it from their
    # rows (but keeps it for the source). Add the source's foreign buildout back to
    # every frozen row so the objective-basis matrix is comparable across rows. The
    # statistics-basis (cost_matrix_total/_outside_de) and DE matrices count frozen
    # capacity via n.statistics already, so they need no offset.
    freeze_source = snakemake.params.freeze_outside_de_capacities
    extra_capex = {}
    if freeze_source and freeze_source in topologies:
        foreign_capex = foreign_portfolio_capex(topologies[freeze_source])
        extra_capex = {p: foreign_capex for p in topologies if p != freeze_source}
        logger.info(
            f"freeze_outside_de_capacities={freeze_source}: adding its foreign "
            f"non-grid capex ({foreign_capex / 1e9:.3f} bn EUR) to every frozen "
            "portfolio's objective-basis cost-matrix row."
        )

    # Network-wide
    capacities = build_capacities_table(topologies)
    capacities.to_csv(snakemake.output.capacities)
    plot_capacities(capacities, snakemake.output.capacities_png, title)

    cost_matrix = build_cost_matrix(
        topologies, evaluations, grid_scenario_names, probabilities, extra_capex
    )
    cost_matrix.to_csv(snakemake.output.cost_matrix)
    plot_cost_matrix(cost_matrix, snakemake.output.cost_matrix_png, title)
    assert_wait_and_see_optimality(cost_matrix, grid_scenario_names)

    metrics = compute_metrics(cost_matrix, grid_scenario_names, probabilities)
    metrics.to_csv(snakemake.output.metrics)
    plot_metrics(metrics, snakemake.output.metrics_png, title)
    logger.info(f"Stochastic-grid metrics (network):\n{metrics}")

    # Whole-network total cost on the statistics basis (incl. existing-capacity
    # capex and the grid), complementing the objective-basis matrix above.
    cost_matrix_total = build_eval_cost_matrix(
        evaluations, grid_scenario_names, portfolios, probabilities, total_cost_stats
    )
    cost_matrix_total.to_csv(snakemake.output.cost_matrix_total)
    plot_cost_matrix(
        cost_matrix_total, snakemake.output.cost_matrix_total_png, f"{title} (total incl. existing capex)"
    )

    # Outside-DE, on the same per-country get_tsc accounting as the DE cost.
    cost_matrix_outside_de = build_eval_cost_matrix(
        evaluations, grid_scenario_names, portfolios, probabilities,
        outside_de_non_trade_cost,
    )
    cost_matrix_outside_de.to_csv(snakemake.output.cost_matrix_outside_de)
    plot_cost_matrix(
        cost_matrix_outside_de, snakemake.output.cost_matrix_outside_de_png, f"{title} (outside DE)"
    )

    # Germany only
    capacities_de = build_capacities_table(topologies, de_capacities_by_carrier)
    capacities_de.to_csv(snakemake.output.capacities_de)
    plot_capacities(capacities_de, snakemake.output.capacities_de_png, f"{title} (DE)")

    cost_matrix_de = build_cost_matrix_de(
        evaluations, grid_scenario_names, portfolios, probabilities
    )
    cost_matrix_de.to_csv(snakemake.output.cost_matrix_de)
    plot_cost_matrix(cost_matrix_de, snakemake.output.cost_matrix_de_png, f"{title} (DE)")

    metrics_de = compute_metrics(cost_matrix_de, grid_scenario_names, probabilities)
    metrics_de.to_csv(snakemake.output.metrics_de)
    plot_metrics(metrics_de, snakemake.output.metrics_de_png, f"{title} (DE)")
    logger.info(f"Stochastic-grid metrics (DE):\n{metrics_de}")

    # DE cost INCLUDING the realized topology's own grid capex (de_total_cost
    # excludes it). Uses the ST evaluation networks, i.e. each column's grid is
    # that scenario's realized single-grid topology.
    cost_matrix_de_grid = build_eval_cost_matrix(
        evaluations, grid_scenario_names, portfolios, probabilities,
        de_total_cost_with_grid,
    )
    cost_matrix_de_grid.to_csv(snakemake.output.cost_matrix_de_grid)
    plot_cost_matrix(
        cost_matrix_de_grid, snakemake.output.cost_matrix_de_grid_png, f"{title} (DE incl. grid)"
    )
