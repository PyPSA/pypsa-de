# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Rolls up the stochastic-grid-scenario analysis into the WS/SP/EEV/EVPI/ECIU
comparison from notebooks/stochastic_grid_uncertainty.ipynb, using only the
already-solved topology and dispatch-only evaluation networks (no re-solving
here):

- WS  (wait-and-see): expected cost if the grid outcome were known before
  investing = sum_s probability_s * objective(solved topology-s)
- SP  (stochastic solution): objective(solved topology-stochastic), PyPSA's
  own expected-cost objective for a scenario-enabled network
- EEV (expected-value planner): capex(topology-eev) + expected opex of
  dispatching the eev portfolio under each real topology
- EVPI = SP - WS   (value of perfect information about the grid)
- ECIU = EEV - SP  (value of stochastic planning over naive EV planning)
"""

import logging

import pandas as pd
import pypsa
from pypsa.descriptors import nominal_attrs

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


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


def build_capacities_table(topologies: dict[str, pypsa.Network]) -> pd.DataFrame:
    return pd.DataFrame(
        {name: capacities_by_carrier(n) for name, n in topologies.items()}
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


def compute_metrics(
    cost_matrix: pd.DataFrame,
    topologies: dict[str, pypsa.Network],
    grid_scenario_names: list[str],
    probabilities: dict[str, float],
) -> pd.Series:
    ws = sum(
        probabilities[s] * topologies[s].objective for s in grid_scenario_names
    )
    sp = topologies["stochastic"].objective
    eev = cost_matrix.loc["eev", "expected"]
    evpi = sp - ws
    eciu = eev - sp

    if not (ws <= sp + 1e-3 <= eev + 1e-3):
        logger.warning(
            f"Expected WS <= SP <= EEV, got WS={ws:.3e}, SP={sp:.3e}, EEV={eev:.3e}. "
            "This may indicate a bug in the scenario-awareness of a custom constraint "
            "(see scripts/solve_network.py / additional_functionality.py)."
        )

    return pd.Series(
        {"WS": ws, "SP": sp, "EEV": eev, "EVPI": evpi, "ECIU": eciu},
        name="value",
    )


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

    capacities = build_capacities_table(topologies)
    capacities.to_csv(snakemake.output.capacities)

    cost_matrix = build_cost_matrix(
        topologies, evaluations, grid_scenario_names, probabilities
    )
    cost_matrix.to_csv(snakemake.output.cost_matrix)

    metrics = compute_metrics(
        cost_matrix, topologies, grid_scenario_names, probabilities
    )
    metrics.to_csv(snakemake.output.metrics)

    logger.info(f"Stochastic-grid metrics:\n{metrics}")
