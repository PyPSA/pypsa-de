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

import copy
import logging

import pandas as pd
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

# Load-shedding is deliberately excluded here: it's a CO2-accounting sink
# (loads are non-positive - "emitted", not consumed), not a real energy
# service, so shedding it would be physically meaningless and could mask
# emissions rather than unmet demand.
LOAD_SHEDDING_EXCLUDED_CARRIERS = {"process emissions"}


def with_load_shedding_carriers(solve_opts: dict, n: pypsa.Network) -> dict:
    """`solve_opts` with `load_shedding.carriers` filled in from `n`'s own loads.

    Carriers found on `n` fall back to `default_cost`, but any per-carrier
    cost already set in config's `load_shedding.carriers` takes precedence -
    a static carrier list would otherwise go stale if the sector
    configuration changes what's connected to this network, so the carrier
    *set* is derived fresh from `n` while letting config override individual
    costs. `all_carriers` is forced off: left on (the config default),
    `add_load_balance_components` (solve_network.py) adds shedding
    generators to every bus *not* in `carriers` too, which would silently
    reintroduce `LOAD_SHEDDING_EXCLUDED_CARRIERS`.
    """
    solve_opts = copy.deepcopy(solve_opts)
    load_shedding = solve_opts.get("load_shedding", {})
    if not load_shedding.get("enable"):
        return solve_opts

    load_bus_carriers = set(n.buses.loc[n.loads.bus.unique(), "carrier"].unique())
    load_bus_carriers -= LOAD_SHEDDING_EXCLUDED_CARRIERS
    default_cost = load_shedding.get("default_cost", 1e5)
    configured_costs = load_shedding.get("carriers", {})
    load_shedding["carriers"] = {
        carrier: configured_costs.get(carrier, default_cost) for carrier in load_bus_carriers
    }
    load_shedding["all_carriers"] = False
    solve_opts["load_shedding"] = load_shedding
    return solve_opts


def check_load_shedding(n: pypsa.Network, threshold: float = 0.05) -> float:
    """Log warnings if shed load exceeds `threshold` of served demand, both
    network-wide and broken down by bus carrier and by individual bus.

    A little shedding at a single tight snapshot is expected numerical
    slack; a meaningful share means this portfolio/topology combination
    couldn't actually meet demand there. The network-wide aggregate alone
    can hide this: a small carrier or a single grid-constrained node can be
    almost entirely unserved while staying under `threshold` of total
    system-wide demand - exactly the failure mode a transmission-topology
    comparison is meant to surface. Returns the network-wide shed fraction
    for the caller to record in `n.meta`; the per-carrier/per-bus breakdowns
    are logged only, not returned.
    """
    shed_i = n.generators.query("carrier == 'load'").index
    if shed_i.empty:
        return 0.0

    weights = n.snapshot_weightings.generators
    shed_by_gen = n.generators_t.p[shed_i].clip(lower=0).mul(weights, axis=0).sum()
    served_by_load = (
        n.get_switchable_as_dense("Load", "p_set").clip(lower=0).mul(weights, axis=0).sum()
    )

    shed_by_bus = shed_by_gen.groupby(n.generators.loc[shed_i, "bus"]).sum()
    served_by_bus = served_by_load.groupby(n.loads.bus).sum()
    bus_i = served_by_bus.index.union(shed_by_bus.index)
    shed_by_bus = shed_by_bus.reindex(bus_i, fill_value=0.0)
    served_by_bus = served_by_bus.reindex(bus_i, fill_value=0.0)

    total_shed = shed_by_bus.sum()
    total_served = served_by_bus.sum()
    fraction = total_shed / total_served if total_served else 0.0
    if fraction > threshold:
        logger.warning(
            f"Load shedding used {fraction:.2%} of total demand - this "
            "portfolio/topology combination likely couldn't fully meet "
            "demand; treat its cost as a capacity-adequacy shortfall "
            "signal, not a normal opex comparison."
        )

    bus_carrier = n.buses.loc[bus_i, "carrier"]
    carrier_shed = shed_by_bus.groupby(bus_carrier).sum()
    carrier_served = served_by_bus.groupby(bus_carrier).sum()
    carrier_fraction = (carrier_shed / carrier_served.replace(0, pd.NA)).dropna()
    for carrier, carrier_frac in carrier_fraction[carrier_fraction > threshold].items():
        logger.warning(
            f"Load shedding for carrier '{carrier}' used {carrier_frac:.2%} of "
            "that carrier's demand, even though this may be masked in the "
            "network-wide aggregate."
        )

    bus_fraction = (shed_by_bus / served_by_bus.replace(0, pd.NA)).dropna()
    breached = bus_fraction[bus_fraction > threshold].sort_values(ascending=False)
    if not breached.empty:
        worst = ", ".join(f"{bus} ({frac:.0%})" for bus, frac in breached.head(10).items())
        more = f", and {len(breached) - 10} more" if len(breached) > 10 else ""
        logger.warning(
            f"Load shedding exceeded {threshold:.0%} of local demand at "
            f"{len(breached)} individual bus(es) - worst: {worst}{more}"
        )

    return fraction


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

    solve_opts = with_load_shedding_carriers(snakemake.params.solving["options"], canvas)
    prepare_network(
        canvas,
        solve_opts=solve_opts,
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

    shed_fraction = check_load_shedding(canvas)

    canvas.meta = dict(
        snakemake.config,
        wildcards=dict(snakemake.wildcards),
        load_shed_fraction=shed_fraction,
    )
    canvas.export_to_netcdf(snakemake.output.network)
