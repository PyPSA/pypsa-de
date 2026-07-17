# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Dispatch-only re-solve: freezes one grid-topology portfolio's optimized
capacities (generators/storage/stores/non-grid links) onto a *different*
topology's fixed grid, then re-solves dispatch only. This is what makes the
WS/SP/EEV/EVPI/ECIU comparison in
notebooks/stochastic_grid_uncertainty.ipynb possible: a portfolio's own
objective only reflects the topology it was optimized for, so evaluating it
under a different realized topology requires this frozen-capacity re-solve.

`n.optimize.fix_optimal_capacities()` only re-fixes a network's *own*
already-solved capacities (it reads `*_opt` from the same network), so the
one extra step here is transplanting the portfolio's `*_opt` values onto the
canvas network before calling it.
"""

import logging

import pypsa
from pypsa.descriptors import nominal_attrs

from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    scenario_slice,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts import solve_network
from scripts.solve_network import collect_kwargs, create_optimization_model, prepare_network

logger = logging.getLogger(__name__)


def transplant_optimal_capacities(canvas: pypsa.Network, portfolio: pypsa.Network) -> None:
    """Copy `portfolio`'s optimized capacities onto `canvas`'s own extendable components.

    Non-anticipativity: first-stage (extendable-capacity) decisions are
    identical across scenarios, so an arbitrary single scenario's slice is
    representative. Using `scenario_slice()` per component here rather than
    `portfolio.get_scenario(...)`, since the latter iterates *every*
    component including `sub_networks` - topology metadata repopulated by
    `n.determine_network_topology()` at solve time, which doesn't carry a
    scenario dimension even though the rest of the network does, and trips
    `get_scenario()`'s blanket assumption that every component's static
    index is a `(scenario, name)` MultiIndex.
    """
    for component, attr in nominal_attrs.items():
        static = canvas.components[component].static
        ext_i = static.query(f"{attr}_extendable").index
        if ext_i.empty:
            continue
        portfolio_static = scenario_slice(portfolio.components[component].static)
        static.loc[ext_i, f"{attr}_opt"] = portfolio_static.loc[ext_i, f"{attr}_opt"]


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "evaluate_grid_portfolio",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
            portfolio="reduced_40",
            grid_scenario="reduced_10",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # extra_functionality() in solve_network.py reads this as a module global,
    # same as when solve_network.py itself runs as __main__.
    solve_network.snakemake = snakemake

    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    canvas = pypsa.Network(snakemake.input.canvas)
    portfolio = pypsa.Network(snakemake.input.portfolio)

    transplant_optimal_capacities(canvas, portfolio)
    canvas.optimize.fix_optimal_capacities()

    prepare_network(
        canvas,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
    )

    model_kwargs, solve_kwargs = collect_kwargs(
        snakemake.config,
        snakemake.params.solving,
        planning_horizons,
        log_fn=snakemake.log.solver,
        mode="single",
    )

    create_optimization_model(
        canvas,
        config=snakemake.config,
        params=snakemake.params,
        model_kwargs=model_kwargs,
        solve_kwargs=solve_kwargs,
        planning_horizons=planning_horizons,
    )
    status, condition = canvas.optimize.solve_model(**solve_kwargs)
    if status != "ok":
        raise RuntimeError(f"Dispatch-only solve failed: {status}, {condition}")

    canvas.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    canvas.export_to_netcdf(snakemake.output.network)
