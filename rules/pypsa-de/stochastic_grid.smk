# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC BY 4.0
#
# Stochastic optimization over uncertain grid-expansion (AC line / DC link)
# capacities, myopic sector-coupled pathway only. See
# notebooks/stochastic_grid_uncertainty.ipynb for the underlying methodology.
#
# This is an opt-in analysis (enabled via config["stochastic_grid_scenarios"]),
# not part of the default `rule all` - request `stochastic_grid_analysis` (or
# any of its intermediate outputs) explicitly.

_sgs = config.get("stochastic_grid_scenarios", {})
GRID_SCENARIO_NAMES = list(_sgs.get("scenarios", {}).keys())
GRID_SCENARIO_VALUES = GRID_SCENARIO_NAMES + ["eev", "stochastic"]


wildcard_constraints:
    grid_scenario="|".join(GRID_SCENARIO_VALUES) if GRID_SCENARIO_VALUES else "(?!)",
    portfolio="|".join(GRID_SCENARIO_VALUES) if GRID_SCENARIO_VALUES else "(?!)",


def _scenario_csv_paths(wildcards, names):
    cfg = config["stochastic_grid_scenarios"]["scenarios"]
    paths = []
    for name in names:
        for key in ("links", "lines"):
            path = cfg.get(name, {}).get(key)
            if path:
                paths.append(path.format(clusters=wildcards.clusters))
    return paths


def get_grid_scenario_csvs(wildcards):
    names = (
        [wildcards.grid_scenario]
        if wildcards.grid_scenario in GRID_SCENARIO_NAMES
        else GRID_SCENARIO_NAMES
    )
    return _scenario_csv_paths(wildcards, names)


def stochastic_grid_solving(wildcards):
    """`solving` config adjusted for PyPSA scenario-dimension gaps.

    Two PyPSA limitations with the `n.set_scenarios()` scenario dimension
    used by the "stochastic" grid_scenario/portfolio, both worked around
    here for all grid-topology variants (not just the scenario-bearing one)
    for consistency across the comparison:

    - `assign_all_duals`: dual-assignment for custom (non-native)
      constraints doesn't support the scenario dimension (raises inside
      `pypsa.optimization.optimize.assign_duals`). Nothing in this analysis
      (WS/SP/EEV/EVPI/ECIU metrics) reads shadow prices, so this is free to
      switch off.
    - `transmission_losses`: post-processing of piecewise-linear line losses
      doesn't support the scenario dimension either (raises inside
      `pypsa.optimization.optimize.post_processing`; reported upstream at
      https://github.com/PyPSA/PyPSA/issues, not yet fixed as of PyPSA
      master). Disabled here rather than blocking on an upstream fix - at
      this pipeline's 27-cluster resolution and given this analysis compares
      *differences* in cost across grid-topology scenarios rather than
      absolute system cost, the lost fidelity is a reasonable tradeoff.
    """
    import copy

    solving = copy.deepcopy(config["solving"])
    solving["options"]["assign_all_duals"] = False
    solving["options"]["transmission_losses"] = 0
    return solving


def evaluate_grid_portfolio_solving(wildcards):
    """`stochastic_grid_solving` plus load-shedding, for evaluate_grid_portfolio only.

    A portfolio sized for one grid topology can genuinely have no feasible
    dispatch when frozen and re-solved on a worse one (e.g. reduced_10's own
    portfolio evaluated against reduced_40's grid: less transmission, same
    fixed capacities) - a real capacity-adequacy shortfall, not a bug.
    Load-shedding turns that into a large-but-finite cost instead of a hard
    solver failure, so one infeasible cell doesn't take down the whole
    compute_stochastic_metrics comparison.

    Deliberately not added to `stochastic_grid_solving` (used by
    solve_grid_topology_network too): letting the capacity-expansion solves
    lean on cheap load-shedding instead of building adequate capacity would
    quietly distort the very capacities this analysis is comparing.
    """
    import copy

    solving = stochastic_grid_solving(wildcards)
    solving = copy.deepcopy(solving)
    solving["options"]["load_shedding"]["enable"] = True
    return solving


rule build_grid_topology:
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_final.nc"
        ),
        csvs=get_grid_scenario_csvs,
    output:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}.nc"
        ),
    log:
        logs(
            "build_grid_topology_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{grid_scenario}.log"
        ),
    threads: 1
    resources:
        mem_mb=4000,
    params:
        stochastic_grid_scenarios=config_provider("stochastic_grid_scenarios"),
        grid_scenario=lambda w: w.grid_scenario,
    message:
        "Building grid-topology variant '{wildcards.grid_scenario}' for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizon"
    script:
        scripts("pypsa-de/build_grid_topology.py")


rule solve_grid_topology_network:
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}.nc"
        ),
        co2_totals_name=resources("co2_totals.csv"),
        energy_totals=resources("energy_totals.csv"),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}.yaml",
        model=(
            RESULTS
            + "models/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}.nc"
            if config["solving"]["options"]["store_model"]
            else []
        ),
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_grid_topology_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}"
        )
    shadow:
        shadow_config
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        solving=stochastic_grid_solving,
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        energy_year=config_provider("energy", "energy_totals_year"),
    message:
        "Solving grid-topology variant '{wildcards.grid_scenario}' for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizon"
    script:
        scripts("solve_network.py")


rule evaluate_grid_portfolio:
    input:
        portfolio=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{portfolio}.nc",
        canvas=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-{grid_scenario}.nc"
        ),
        co2_totals_name=resources("co2_totals.csv"),
        energy_totals=resources("energy_totals.csv"),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_portfolio-{portfolio}_on-{grid_scenario}_st.nc",
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_portfolio-{portfolio}_on-{grid_scenario}_st_solver.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_portfolio-{portfolio}_on-{grid_scenario}_st_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/evaluate_grid_portfolio/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_portfolio-{portfolio}_on-{grid_scenario}"
        )
    shadow:
        shadow_config
    threads: 4
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        options=config_provider("solving", "options"),
        solving=evaluate_grid_portfolio_solving,
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        energy_year=config_provider("energy", "energy_totals_year"),
    message:
        "Evaluating portfolio '{wildcards.portfolio}' capacities dispatched on grid topology '{wildcards.grid_scenario}'"
    script:
        scripts("pypsa-de/evaluate_grid_portfolio.py")


def get_topology_networks(wildcards):
    return {
        f"topology_{grid_scenario}": RESULTS
        + f"networks/base_s_{wildcards.clusters}_{wildcards.opts}_{wildcards.sector_opts}_{wildcards.planning_horizons}_topology-{grid_scenario}.nc"
        for grid_scenario in GRID_SCENARIO_VALUES
    }


def get_evaluation_networks(wildcards):
    return {
        f"eval_{portfolio}__{grid_scenario}": RESULTS
        + f"networks/base_s_{wildcards.clusters}_{wildcards.opts}_{wildcards.sector_opts}_{wildcards.planning_horizons}_portfolio-{portfolio}_on-{grid_scenario}_st.nc"
        for portfolio in GRID_SCENARIO_VALUES
        for grid_scenario in GRID_SCENARIO_NAMES
    }


rule compute_stochastic_metrics:
    input:
        # unpack() flattens each function's returned dict into named inputs
        # directly (snakemake.input.topology_reduced_40,
        # snakemake.input.eval_reduced_40__reduced_10, ...) - matches the
        # existing `unpack(input_network_year)` pattern in solve_perfect.smk.
        unpack(get_topology_networks),
        unpack(get_evaluation_networks),
    output:
        capacities=RESULTS
        + "stochastic_grid_capacities_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        cost_matrix=RESULTS
        + "stochastic_grid_cost_matrix_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        metrics=RESULTS
        + "stochastic_grid_metrics_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    log:
        logs(
            "compute_stochastic_metrics_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    threads: 1
    resources:
        mem_mb=4000,
    params:
        grid_scenario_names=GRID_SCENARIO_NAMES,
        grid_scenario_values=GRID_SCENARIO_VALUES,
        scenarios=config_provider("stochastic_grid_scenarios", "scenarios"),
    message:
        "Computing WS/SP/EEV/EVPI/ECIU stochastic-grid metrics for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizon"
    script:
        scripts("pypsa-de/compute_stochastic_metrics.py")


rule stochastic_grid_analysis:
    input:
        expand(
            RESULTS
            + "stochastic_grid_metrics_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
            run=config["run"]["name"],
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=_sgs.get("planning_horizons", []),
        ),
        # Unconstrained reference network (no grid_scenario override, so grid
        # Lines/DC-Links stay freely extendable rather than fixed to a
        # reduced-capacity CSV) - "no delay, grid can expand to its own
        # optimum". Not part of the WS/SP/EEV/EVPI/ECIU comparison (it isn't
        # capex-comparable - grid buildout cost is endogenous here but
        # excluded-as-exogenous everywhere else), just produced alongside it
        # for manual inspection.
        expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            run=config["run"]["name"],
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=_sgs.get("planning_horizons", []),
        ),
    message:
        "Collected all stochastic-grid-scenario analysis outputs."
