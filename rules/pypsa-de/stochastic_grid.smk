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

# Run scope - trade completeness for speed:
# - "full": LT topology variants (incl. the native stochastic solve), the ST
#   portfolio evaluations, the WS/SP/EEV/EVPI/ECIU metrics, and the LT + ST
#   scenario-comparison plots.
# - "lt_only": only the LT topology variants (deterministic, eev AND the
#   stochastic solve) with their system / LT-comparison plots. No ST portfolio
#   evaluations and no stochastic metrics (both need the ST re-dispatch runs).
# - "lt_deterministic": like lt_only but the native stochastic solve is skipped
#   too, so no PyPSA scenario dimension is built at all - only the deterministic
#   + eev single-scenario topologies. The fastest path to topology results (and
#   it needs no PyPSA #1951 / grid-master env, since nothing is scenario-indexed).
_GRID_RUN_MODE = _sgs.get("run_mode", "full")
assert _GRID_RUN_MODE in {"full", "lt_only", "lt_deterministic"}, (
    "stochastic_grid_scenarios.run_mode must be one of "
    f"'full'/'lt_only'/'lt_deterministic', got {_GRID_RUN_MODE!r}"
)
_GRID_INCLUDE_STOCHASTIC = _GRID_RUN_MODE != "lt_deterministic"
_GRID_INCLUDE_ST = _GRID_RUN_MODE == "full"

GRID_SCENARIO_VALUES = (
    GRID_SCENARIO_NAMES + ["eev"] + (["stochastic"] if _GRID_INCLUDE_STOCHASTIC else [])
)

# Investigation year(s) the grid analysis is built and solved for. Single source
# of truth = scenario.planning_horizons (for the standalone single-year analysis,
# set that to one year); an explicit stochastic_grid_scenarios.planning_horizons
# still overrides it with a subset if given.
_GRID_HORIZONS = _sgs.get("planning_horizons") or config.get("scenario", {}).get(
    "planning_horizons", []
)

# First myopic-pathway horizon - the pristine, never-solved prenetwork used as
# the "grid stays put" baseline when auto-generating `*_exogen` scenario CSVs.
_GRID_FIRST_HORIZON = (
    str(min(config["scenario"]["planning_horizons"]))
    if config.get("scenario", {}).get("planning_horizons")
    else None
)


# "topology-stochastic" itself is excluded: it's the actual two-stage
# stochastic solve, kept on PyPSA's native scenario dimension (MultiIndex
# over the grid_scenario branches) rather than resolved into a single
# network. The Ariadne-variable export (like `assign_all_duals`/
# `transmission_losses` elsewhere in this file, see `stochastic_grid_solving`
# above) doesn't support that dimension. Its branches are unstacked into
# plain single-scenario networks instead (see
# unstack_stochastic_grid_topology below) and exported as
# "topology-stochastic__{grid_scenario}". All ST portfolio evaluations -
# including the "stochastic" portfolio's capacities dispatched onto one
# canvas - and the other LT topology variants are already plain
# single-scenario networks and are exported directly.
GRID_EXPORT_NETWORK_IDS_LT_TOPOLOGY = [
    f"topology-{gs}" for gs in GRID_SCENARIO_NAMES + ["eev"]
] + (
    [f"topology-stochastic__{gs}" for gs in GRID_SCENARIO_NAMES]
    if _GRID_INCLUDE_STOCHASTIC
    else []
)


def _grid_pathway_years(target_year):
    """Configured myopic-pathway years up to and including `target_year`,
    oldest first - e.g. [2025, 2030] for target_year=2030. Used to mix the
    plain (suffix-less) `base_s_..._{year}.nc` pathway networks into the LT
    comparison alongside the target year's grid-topology variants, so the
    buildout leading up to the stochastic-grid analysis year is visible.
    """
    return sorted(y for y in config["scenario"]["planning_horizons"] if y <= target_year)


# "optimal" aliases whichever year a given rule instance's own
# {planning_horizons} wildcard is bound to; "optimal_YYYY" is a fixed
# earlier pathway year mixed into that (otherwise single-year) instance.
# The two never overlap for one target year, but the wildcard-constraint
# regex below is a static superset across every stochastic_grid_scenarios
# target year (usually just one), so both forms are listed.
_GRID_PATHWAY_YEARS_ALL = sorted(
    {
        y
        for target in _GRID_HORIZONS
        for y in _grid_pathway_years(target)
        if y != target
    }
)

# The endogenous "optimal"/"optimal_YYYY" reference networks are only built,
# solved and compared when opted into (see solve_endogenous_optimal); otherwise
# the comparison is purely over the grid-topology variants.
_GRID_SOLVE_OPTIMAL = _sgs.get("solve_endogenous_optimal", False)
_GRID_OPTIMAL_IDS = (
    ["optimal"] + [f"optimal_{y}" for y in _GRID_PATHWAY_YEARS_ALL]
    if _GRID_SOLVE_OPTIMAL
    else []
)
GRID_EXPORT_NETWORK_IDS_LT = GRID_EXPORT_NETWORK_IDS_LT_TOPOLOGY + _GRID_OPTIMAL_IDS
GRID_EXPORT_NETWORK_IDS_ST = (
    [
        f"portfolio-{p}_on-{gs}_st"
        for p in GRID_SCENARIO_VALUES
        for gs in GRID_SCENARIO_NAMES
    ]
    if _GRID_INCLUDE_ST
    else []
)
GRID_EXPORT_NETWORK_IDS = GRID_EXPORT_NETWORK_IDS_LT + GRID_EXPORT_NETWORK_IDS_ST


wildcard_constraints:
    grid_scenario="|".join(GRID_SCENARIO_VALUES) if GRID_SCENARIO_VALUES else "(?!)",
    portfolio="|".join(GRID_SCENARIO_VALUES) if GRID_SCENARIO_VALUES else "(?!)",
    network_id="|".join(GRID_EXPORT_NETWORK_IDS) if GRID_EXPORT_NETWORK_IDS else "(?!)",


def _year_for_network_id(network_id, fallback_year):
    """Actual data year for one network_id: "optimal_YYYY" pins its own
    year, everything else (including plain "optimal") uses the enclosing
    rule instance's own {planning_horizons} wildcard value."""
    if network_id.startswith("optimal_"):
        return network_id.removeprefix("optimal_")
    return fallback_year


def _grid_network_path(clusters, opts, sector_opts, network_id, fallback_year):
    """.nc path for one network_id - suffix-less for "optimal"/"optimal_YYYY"
    (the plain pathway-reference network), `_{network_id}` suffixed for the
    regular topology-*/portfolio-*_st grid-scenario variants."""
    year = _year_for_network_id(network_id, fallback_year)
    suffix = "" if network_id == "optimal" or network_id.startswith("optimal_") else f"_{network_id}"
    return RESULTS + f"networks/base_s_{clusters}_{opts}_{sector_opts}_{year}{suffix}.nc"


def _grid_prenetwork_path(clusters, opts, sector_opts, network_id, fallback_year):
    """Pre-solve counterpart of _grid_network_path's network_id, for the
    grid-topology comparison maps (system_plots.py's
    plot_grid_topology_maps) - the network as handed to the solver, with
    whatever grid CSV override (if any) is already applied but before
    capacity/dispatch optimization:
    - "optimal"/"optimal_YYYY": the plain prepared network (no override).
    - "topology-{gs}"/"topology-stochastic__{gs}": {gs}'s own
      build_grid_topology output (the joint stochastic solve's
      "topology-stochastic__{gs}" branch starts from the exact same grid
      CSV override as the deterministic "topology-{gs}" one).
    - "portfolio-{p}_on-{gs}_st": {gs}'s build_grid_topology output too -
      that's the frozen "canvas" grid topology being evaluated against.
    """
    year = _year_for_network_id(network_id, fallback_year)
    if network_id.startswith("topology-stochastic__"):
        gs = network_id.removeprefix("topology-stochastic__")
    elif network_id.startswith("topology-"):
        gs = network_id.removeprefix("topology-")
    elif network_id.startswith("portfolio-") and network_id.endswith("_st"):
        gs = network_id.rsplit("_on-", 1)[1].removesuffix("_st")
    else:
        gs = None
    suffix = f"_topology-{gs}" if gs is not None else "_final"
    template = resources(
        "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}" + suffix + ".nc"
    )
    return (
        template.replace("{clusters}", clusters)
        .replace("{opts}", opts)
        .replace("{sector_opts}", sector_opts)
        .replace("{planning_horizons}", year)
    )


def _grid_scenario_generated_csv(wildcards, name, key):
    """Per-run resources path of an auto-generated scenario CSV (see
    build_grid_scenario_csvs). Used whenever a scenario has no explicit CSV
    path in config, so overrides are re-derived from the current topology
    instead of a hand-maintained file that goes stale on reclustering."""
    template = resources(
        "grid_scenarios/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_"
        + f"{name}_{key}.csv"
    )
    path = (
        template.replace("{clusters}", wildcards.clusters)
        .replace("{opts}", wildcards.opts)
        .replace("{sector_opts}", wildcards.sector_opts)
        .replace("{planning_horizons}", wildcards.planning_horizons)
    )
    # resources() leaves the {run} scenario wildcard for snakemake to expand,
    # which it does for input strings but not inside the nested params dict this
    # feeds (params.scenario_csvs in build_grid_topology), so resolve it here.
    return path.replace("{run}", wildcards.run) if "{run}" in path else path


def _grid_scenario_csv_map(wildcards):
    """{scenario: {'links': path, 'lines': path}} for the scenarios this
    build_grid_topology instance needs. Explicit config paths win (manual
    override); otherwise the auto-generated per-run resources path is used."""
    names = (
        [wildcards.grid_scenario]
        if wildcards.grid_scenario in GRID_SCENARIO_NAMES
        else GRID_SCENARIO_NAMES
    )
    cfg = config["stochastic_grid_scenarios"]["scenarios"]
    result = {}
    for name in names:
        entry = {}
        for key in ("links", "lines"):
            explicit = cfg.get(name, {}).get(key)
            entry[key] = (
                explicit.format(clusters=wildcards.clusters)
                if explicit
                else _grid_scenario_generated_csv(wildcards, name, key)
            )
        result[name] = entry
    return result


def get_grid_scenario_csvs(wildcards):
    return [p for entry in _grid_scenario_csv_map(wildcards).values() for p in entry.values()]


def _grid_scenario_csv_source_networks(wildcards):
    """Networks build_grid_scenario_csvs reads: the canvas (target-year
    prenetwork the CSV is applied to) plus either the pristine first-horizon
    prenetwork (kind 'exogen') or the solved postnetwork of the scenario's own
    year (kind 'optimal'). Name is '<year>_<kind>[_<carrier>][_<direction>]';
    only year and kind pick the source network, so trailing filters are ignored
    here (see build_grid_scenario_csvs.parse_grid_scenario)."""
    year, kind = wildcards.grid_scenario.split("_")[:2]
    ins = {
        "canvas": resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_final.nc"
        )
    }
    if kind == "optimal":
        ins["solved"] = (
            RESULTS
            + f"networks/base_s_{wildcards.clusters}_{wildcards.opts}_{wildcards.sector_opts}_{year}.nc"
        )
    else:
        ins["pristine"] = (
            resources(
                "networks/base_s_{clusters}_{opts}_{sector_opts}_"
                + f"{_GRID_FIRST_HORIZON}_final.nc"
            )
            .replace("{clusters}", wildcards.clusters)
            .replace("{opts}", wildcards.opts)
            .replace("{sector_opts}", wildcards.sector_opts)
        )
    return ins


def _deep_update(base, override):
    """Recursively merge ``override`` into ``base`` (in place) and return it."""
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_update(base[k], v)
        else:
            base[k] = v
    return base


def stochastic_grid_solving(wildcards):
    """`solving` config adjusted for PyPSA scenario-dimension gaps.

    Two PyPSA limitations with the `n.set_scenarios()` scenario dimension used
    by the "stochastic" grid_scenario:

    - `assign_all_duals`: dual-assignment for custom (non-native) constraints
      doesn't support the scenario dimension (raises inside
      `pypsa.optimization.optimize.assign_duals`). Nothing in this analysis
      (WS/SP/EEV/EVPI/ECIU metrics) reads those shadow prices, so it is forced
      off for the stochastic solve ONLY; the single-scenario variants
      (deterministic, eev, portfolio evaluations) keep the configured
      `solving.options.assign_all_duals` - it works fine without a scenario
      dimension and it does not change the optimum, only whether duals are
      stored.
    - `transmission_losses`: post-processing of piecewise-linear line losses
      doesn't support the scenario dimension either (raises inside
      `pypsa.optimization.optimize.post_processing`; reported upstream at
      https://github.com/PyPSA/PyPSA/issues, not yet fixed as of PyPSA
      master). Disabled here for *all* variants (unlike assign_all_duals): it
      changes the physics/optimum, so keeping it uniform across the comparison
      matters, and this analysis compares *differences* in cost across
      grid-topology scenarios rather than absolute system cost.
    """
    import copy

    solving = copy.deepcopy(config["solving"])
    solving["options"]["transmission_losses"] = 0
    if wildcards.grid_scenario == "stochastic":
        # `assign_all_duals` is forced off for the scenario-dimension solve
        # only (it raises there); single-scenario variants keep the configured
        # solving.options.assign_all_duals. Then apply any stochastic-only
        # `solving` overrides (e.g. solver_options, options.*), deep-merged for
        # this solve alone - free-form, so new solver settings need no rule
        # change (add them under `stochastic_grid_scenarios: stochastic_solving`).
        solving["options"]["assign_all_duals"] = False
        _deep_update(solving, _sgs.get("stochastic_solving") or {})
    return solving


def evaluate_grid_portfolio_solving(wildcards):
    """`stochastic_grid_solving` plus load-shedding, for evaluate_grid_portfolio only.

    A portfolio sized for one grid topology can genuinely have no feasible
    dispatch when frozen and re-solved on a worse one (e.g. 2035_exogen's own
    portfolio evaluated against 2025_exogen's grid: less transmission, same
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


rule build_grid_scenario_csvs:
    input:
        unpack(_grid_scenario_csv_source_networks),
    output:
        lines=resources(
            "grid_scenarios/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{grid_scenario}_lines.csv"
        ),
        links=resources(
            "grid_scenarios/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{grid_scenario}_links.csv"
        ),
    log:
        logs(
            "build_grid_scenario_csvs_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{grid_scenario}.log"
        ),
    threads: 1
    resources:
        mem_mb=4000,
    params:
        scope=config_provider(
            "stochastic_grid_scenarios", "scope", default="de_and_interconnectors"
        ),
    message:
        "Generating grid-scenario CSVs for '{wildcards.grid_scenario}' at {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizon"
    script:
        scripts("pypsa-de/build_grid_scenario_csvs.py")


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
        mem_mb=16000,
    params:
        # Only the `scenarios` sub-dict, NOT the whole stochastic_grid_scenarios
        # block: sibling keys (run_mode, solve_endogenous_optimal,
        # stochastic_solve_resources, ...) don't affect the built topology
        # network, so passing them here would needlessly rerun this rule (and
        # the expensive downstream solve) via snakemake's `params` rerun-trigger
        # whenever one is toggled - e.g. switching run_mode lt_only -> full.
        grid_scenarios=config_provider("stochastic_grid_scenarios", "scenarios"),
        outside_de_grid=config_provider(
            "stochastic_grid_scenarios", "outside_de_grid", default="endogenous"
        ),
        grid_scenario=lambda w: w.grid_scenario,
        scenario_csvs=_grid_scenario_csv_map,
    message:
        "Building grid-topology variant '{wildcards.grid_scenario}' for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizon"
    script:
        scripts("pypsa-de/build_grid_topology.py")


def grid_topology_solve_mem_mb(wildcards):
    """Solver memory for solve_grid_topology_network.

    The "stochastic" variant carries PyPSA's native scenario dimension, so its LP is
    much larger - and scales super-linearly in the number of scenarios - than the
    single-scenario variants. Its memory is set explicitly via
    `stochastic_grid_scenarios: stochastic_solve_resources: mem_mb`; unset (or any
    other variant) falls back to the general `solving: mem_mb`.
    """
    override = (_sgs.get("stochastic_solve_resources") or {}).get("mem_mb")
    if wildcards.grid_scenario == "stochastic" and override is not None:
        return override
    return config_provider("solving", "mem_mb")(wildcards)


def grid_topology_solve_runtime(wildcards):
    """Runtime for solve_grid_topology_network, analogous to
    grid_topology_solve_mem_mb: the "stochastic" variant may set
    `stochastic_grid_scenarios: stochastic_solve_resources: runtime`, otherwise the
    general `solving: runtime` applies (as it does for every other variant).
    """
    override = (_sgs.get("stochastic_solve_resources") or {}).get("runtime")
    if wildcards.grid_scenario == "stochastic" and override is not None:
        return override
    return config_provider("solving", "runtime", default="6h")(wildcards)


def grid_topology_solve_threads(wildcards):
    """Thread count for solve_grid_topology_network, analogous to
    grid_topology_solve_mem_mb: the "stochastic" variant may set
    `stochastic_grid_scenarios: stochastic_solve_resources: threads`, otherwise
    the general `solver_threads` applies (as for every other variant)."""
    override = (_sgs.get("stochastic_solve_resources") or {}).get("threads")
    if wildcards.grid_scenario == "stochastic" and override is not None:
        return override
    return solver_threads(wildcards)


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
    threads: grid_topology_solve_threads
    resources:
        mem_mb=grid_topology_solve_mem_mb,
        runtime=grid_topology_solve_runtime,
    params:
        solving=stochastic_grid_solving,
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        energy_year=config_provider("energy", "energy_totals_year"),
    message:
        "Solving grid-topology variant '{wildcards.grid_scenario}' for {wildcards.clusters} clusters, "
        "{wildcards.planning_horizons} planning horizon "
        "(threads={threads}, mem_mb={resources.mem_mb}, runtime={resources.runtime} min)"
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


rule unstack_stochastic_grid_topology:
    """Extract one branch of the joint two-stage stochastic solve.

    `topology-stochastic` keeps PyPSA's native scenario dimension across
    both grid_scenario branches at once; this pulls one branch out as a
    plain single-scenario network so it can go through the same
    Ariadne-variable export as everything else (see network_id
    "topology-stochastic__{grid_scenario}" above).
    """
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-stochastic.nc",
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_topology-stochastic__{grid_scenario}.nc",
    log:
        logs(
            "unstack_stochastic_grid_topology_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{grid_scenario}.log"
        ),
    threads: 1
    resources:
        mem_mb=16000,
    message:
        "Unstacking scenario '{wildcards.grid_scenario}' from the joint stochastic grid-topology network"
    script:
        scripts("pypsa-de/unstack_stochastic_grid_topology.py")


rule export_ariadne_variables_grid_scenario:
    """Ariadne-variable export for one grid-topology/portfolio network variant.

    Reuses export_ariadne_variables.py with a single-network input (instead
    of the full myopic-pathway network list), since each stochastic-grid
    variant is a standalone one-year network. All variants share the same
    `run` wildcard/scenario name, so - unlike the main export rule - the
    output is distinguished by `network_id` (e.g. "topology-stochastic" or
    "portfolio-eev_on-2025_exogen_st") rather than by scenario.
    """
    input:
        template="data/template_ariadne_database.xlsx",
        industry_demands=lambda w: [
            resources(
                "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
            )
            .replace("{clusters}", w.clusters)
            .replace("{planning_horizons}", _year_for_network_id(w.network_id, w.planning_horizons))
        ],
        networks=lambda w: [
            _grid_network_path(w.clusters, w.opts, w.sector_opts, w.network_id, w.planning_horizons)
        ],
        costs=lambda w: [
            resources("costs_{planning_horizons}_processed.csv").replace(
                "{planning_horizons}", _year_for_network_id(w.network_id, w.planning_horizons)
            )
        ],
        industrial_production_per_country_tomorrow=lambda w: [
            resources(
                "industrial_production_per_country_tomorrow_{planning_horizons}-modified.csv"
            ).replace(
                "{planning_horizons}", _year_for_network_id(w.network_id, w.planning_horizons)
            )
        ],
        industry_sector_ratios=lambda w: [
            resources("industry_sector_ratios_{planning_horizons}.csv").replace(
                "{planning_horizons}", _year_for_network_id(w.network_id, w.planning_horizons)
            )
        ],
        industrial_production=resources("industrial_production_per_country.csv"),
        energy_totals=resources("energy_totals.csv"),
    output:
        exported_variables=RESULTS
        + "ariadne/exported_variables_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{network_id}.xlsx",
        exported_variables_full=RESULTS
        + "ariadne/exported_variables_full_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{network_id}.xlsx",
    log:
        logs(
            "export_ariadne_variables_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{network_id}.log"
        ),
    resources:
        mem_mb=16000,
    params:
        planning_horizons=lambda w: [int(_year_for_network_id(w.network_id, w.planning_horizons))],
        config_industry=config_provider("industry"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        post_discretization=config_provider("solving", "options", "post_discretization"),
        NEP_year=lambda w: config_provider("costs", "custom_cost_fn")(w)[-8:-4],
        NEP_transmission=config_provider("costs", "transmission"),
    message:
        "Exporting Ariadne variables for grid-topology network variant '{wildcards.network_id}'"
    script:
        scripts("pypsa-de/export_ariadne_variables.py")


rule system_plots_grid_scenario:
    """system_plots.py for one grid-topology/portfolio network variant.

    Single-network/single-year counterpart of `system_plots` (which covers
    the full myopic pathway for a `run`), analogous to how
    export_ariadne_variables_grid_scenario relates to export_ariadne_variables.
    One output folder per `network_id`, since - like the export rule - all
    variants share the same `run`/scenario name.
    """
    params:
        planning_horizons=lambda w: [int(_year_for_network_id(w.network_id, w.planning_horizons))],
        plotting=config_provider("plotting"),
        output_dir=RESULTS + "system/plots/{network_id}",
    input:
        networks=lambda w: [
            _grid_network_path(w.clusters, w.opts, w.sector_opts, w.network_id, w.planning_horizons)
        ],
        prenetworks=lambda w: [
            _grid_prenetwork_path(w.clusters, w.opts, w.sector_opts, w.network_id, w.planning_horizons)
        ],
        regions_onshore=expand(
            resources("regions_onshore_base_s_{clusters}.geojson"),
            allow_missing=True,
        ),
    output:
        flag=touch(
            RESULTS
            + "system/plots/{network_id}/.system_plots_complete_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.flag"
        ),
    resources:
        mem_mb=40000,
    log:
        logs(
            "system_plots_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{network_id}.log"
        ),
    message:
        "Plotting system plots for grid-topology network variant '{wildcards.network_id}'"
    script:
        scripts("pypsa-de/system_plots.py")


def _lt_network_ids(wildcards):
    """LT comparison set for one rule instance: the grid-topology variants
    (all at the instance's own {planning_horizons} year), plus - only when
    solve_endogenous_optimal is on - one suffix-less endogenous pathway network
    per configured year up to and including that target year ("optimal" for the
    target year itself, "optimal_YYYY" for each earlier one)."""
    if not _GRID_SOLVE_OPTIMAL:
        return GRID_EXPORT_NETWORK_IDS_LT_TOPOLOGY
    target_year = int(wildcards.planning_horizons)
    pathway_ids = [
        "optimal" if y == target_year else f"optimal_{y}"
        for y in _grid_pathway_years(target_year)
    ]
    return GRID_EXPORT_NETWORK_IDS_LT_TOPOLOGY + pathway_ids


def _comparison_network_ids(wildcards):
    return _lt_network_ids(wildcards) if wildcards.variant == "LT" else GRID_EXPORT_NETWORK_IDS_ST


rule plot_grid_scenario_comparison:
    """Compare LT topology or ST portfolio-evaluation network variants.

    All variants share the same `run`/Scenario name (see
    export_ariadne_variables_grid_scenario above), so the comparison is
    keyed by `network_id` instead - one directory of plots per `variant`
    (LT topology-* networks, or ST portfolio-*_on-*_st networks).
    """
    input:
        networks=lambda w: [
            _grid_network_path(w.clusters, w.opts, w.sector_opts, nid, w.planning_horizons)
            for nid in _comparison_network_ids(w)
        ],
        exported_variables=lambda w: [
            RESULTS
            + f"ariadne/exported_variables_full_base_s_{w.clusters}_{w.opts}_{w.sector_opts}_{w.planning_horizons}_{nid}.xlsx"
            for nid in _comparison_network_ids(w)
        ],
        regions_onshore=expand(
            resources("regions_onshore_base_s_{clusters}.geojson"),
            allow_missing=True,
        ),
    output:
        directory(
            RESULTS
            + "scenario_comparison/{variant}/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    log:
        logs(
            "plot_grid_scenario_comparison_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{variant}.log"
        ),
    threads: 1
    resources:
        mem_mb=50000,
        runtime=480,  
    params:
        network_ids=lambda w: _comparison_network_ids(w),
        plotting=config_provider("plotting"),
    wildcard_constraints:
        variant="LT|ST",
    message:
        "Plotting {wildcards.variant} grid-scenario comparison"
    script:
        scripts("pypsa-de/plot_grid_scenario_comparison.py")


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
        # directly (snakemake.input.topology_2025_exogen,
        # snakemake.input.eval_2035_exogen__2025_exogen, ...) - matches the
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
        mem_mb=80000,
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
        # WS/SP/EEV/EVPI/ECIU metrics - only in the "full" run_mode (they need
        # the ST portfolio-evaluation re-dispatch runs).
        (
            expand(
                RESULTS
                + "stochastic_grid_metrics_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
                run=config["run"]["name"],
                clusters=config["scenario"]["clusters"],
                opts=config["scenario"]["opts"],
                sector_opts=config["scenario"]["sector_opts"],
                planning_horizons=_GRID_HORIZONS,
            )
            if _GRID_INCLUDE_ST
            else []
        ),
        # Endogenous "optimal" reference network (no grid_scenario override, so
        # grid Lines/DC-Links stay freely extendable rather than fixed to a
        # reduced-capacity CSV) - "no delay, grid can expand to its own optimum".
        # Opt-in via `stochastic_grid_scenarios: solve_endogenous_optimal`; not
        # part of the WS/SP/EEV/EVPI/ECIU comparison (grid buildout cost is
        # endogenous here but excluded-as-exogenous everywhere else), just
        # produced alongside it for manual inspection.
        (
            expand(
                RESULTS
                + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
                run=config["run"]["name"],
                clusters=config["scenario"]["clusters"],
                opts=config["scenario"]["opts"],
                sector_opts=config["scenario"]["sector_opts"],
                planning_horizons=_GRID_HORIZONS,
            )
            if _sgs.get("solve_endogenous_optimal", False)
            else []
        ),
        # Ariadne-variable export (LT topology + ST portfolio-evaluation
        # networks), one xlsx per network variant since they all share the
        # same `run`/scenario name.
        expand(
            RESULTS
            + "ariadne/exported_variables_full_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{network_id}.xlsx",
            run=config["run"]["name"],
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=_GRID_HORIZONS,
            network_id=GRID_EXPORT_NETWORK_IDS,
        ),
        # Scenario-comparison plots - always LT; ST only in the "full" run_mode.
        expand(
            RESULTS
            + "scenario_comparison/{variant}/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=_GRID_HORIZONS,
            variant=(["LT", "ST"] if _GRID_INCLUDE_ST else ["LT"]),
        ),
        # Per-network-variant system plots (storage/capacity maps, energy
        # balances, capacity comparison) - one folder per network variant.
        expand(
            RESULTS
            + "system/plots/{network_id}/.system_plots_complete_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.flag",
            run=config["run"]["name"],
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=_GRID_HORIZONS,
            network_id=GRID_EXPORT_NETWORK_IDS,
        ),
    message:
        "Collected all stochastic-grid-scenario analysis outputs."
