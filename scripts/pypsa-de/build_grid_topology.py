# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Builds one grid-topology variant of a network for stochastic optimization
over uncertain grid-expansion (AC line / DC link) capacities, see
notebooks/stochastic_grid_uncertainty.ipynb for the underlying methodology.

Three kinds of `grid_scenario` wildcard value, all built from the same input
network:

- a named scenario (e.g. "2025_exogen"): apply that scenario's CSV overrides
  onto a deterministic (non-scenario) network - "perfect information" for
  that one topology.
- "eev": apply the probability-weighted average of all scenarios' CSV
  overrides onto a single deterministic network - the naive planner who
  ignores uncertainty and builds for the "expected" grid.
- "stochastic": call n.set_scenarios(...) and apply every scenario's own
  CSV overrides onto its own (scenario, name) slice - the joint two-stage
  stochastic program PyPSA solves natively.

Every variant scales a reduced AC line's reactance/resistance with its capacity
(x, r ~ 1/num_parallel ~ 1/s_nom, so a smaller corridor is more reactive and
more resistive). The deterministic and eev (single-scenario) variants also drop
lines delayed to zero capacity from the network - and thus from the Kirchhoff
cycle - rather than leave a phantom cycle (see _set_line_capacity). The
stochastic variant scales num_parallel per (scenario, name) slice, so each
scenario carries its own impedance; the zeroed lines are instead removed from
each scenario's cycle at solve time by additional_functionality.add_scenario_kvl,
which replaces PyPSA's native scenario KVL (that one reads impedance from
scenarios[0] only).
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

NOM_ATTRS = {"Link": "p_nom", "Line": "s_nom"}
CSV_KEYS = {"Link": "links", "Line": "lines"}


def _read_csv(path, clusters):
    df = pd.read_csv(path.format(clusters=clusters), index_col=0)
    df.index = df.index.astype(str)  # component names are always strings in PyPSA
    return df


def _check_names_exist(df, network_index, component, csv_path):
    missing = df.index.difference(network_index)
    if len(missing):
        raise ValueError(
            f"{component} name(s) {list(missing)} from {csv_path} not found in the "
            f"network. Check that this CSV matches the `clusters` resolution in use."
        )


def _apply_overrides(static, df, component, csv_path):
    nom_attr = NOM_ATTRS[component]
    _check_names_exist(df, static.index, component, csv_path)
    static.loc[df.index, df.columns] = df.values
    if f"{nom_attr}_extendable" not in df.columns:
        static.loc[df.index, f"{nom_attr}_extendable"] = False


def _set_line_capacity(n, s_nom_new):
    """
    Override AC-line ``s_nom`` and co-scale the physical line parameters.

    A line's series reactance and resistance scale inversely with the number of
    parallel circuits (``x = x_per_length * length / num_parallel``, likewise
    ``r``), hence inversely with ``s_nom``. Every line carries a standard
    ``type``, so ``n.calculate_dependent_values()`` recomputes ``x``/``r`` from
    ``num_parallel`` at solve time: scaling ``num_parallel`` by the
    scenario/canvas capacity ratio is therefore the single lever that makes a
    delayed corridor correctly more reactive *and* more resistive.

    A corridor delayed to ``s_nom=0`` is absent, so it is removed from the
    network - and thus from the Kirchhoff cycle - rather than left as a
    phantom cycle that would distort flows on the rest of its loop.

    ``n`` must still hold the canvas (fully-built) values when called; only
    deterministic/eev single-scenario networks are passed here. The stochastic
    network scales impedance per scenario slice via
    _scale_scenario_line_capacity instead (it cannot drop lines per scenario).

    Parameters
    ----------
    n : pypsa.Network
        Single-scenario network whose lines carry canvas values.
    s_nom_new : pandas.Series
        Target ``s_nom`` per line name (canvas value for undelayed lines, a
        reduced value for delayed ones, 0 for lines not yet built).
    """
    lines = n.components["Line"].static
    frac = s_nom_new.where(s_nom_new > 0) / lines.loc[s_nom_new.index, "s_nom"]
    live = frac.dropna().index
    dropped = s_nom_new.index.difference(live)

    lines.loc[live, "num_parallel"] *= frac[live]
    lines.loc[live, "s_nom"] = s_nom_new[live]
    lines.loc[live, "s_nom_extendable"] = False
    if len(dropped):
        n.remove("Line", list(dropped))
    logger.info(
        "Lines: scaled num_parallel (-> x/r) on %d reduced corridor(s), dropped "
        "%d zeroed corridor(s).",
        len(live),
        len(dropped),
    )


def build_deterministic_topology(n, scenario_cfg, clusters):
    """Apply a single scenario's overrides onto a plain deterministic network."""
    for component, key in CSV_KEYS.items():
        path = scenario_cfg.get(key)
        if not path:
            continue
        df = _read_csv(path, clusters)
        if component == "Line":
            _check_names_exist(df, n.components["Line"].static.index, "Line", path)
            _set_line_capacity(n, df["s_nom"])
        else:
            _apply_overrides(n.components[component].static, df, component, path)
    return n


def build_eev_topology(n, scenarios_cfg, clusters):
    """Probability-weighted blend of all scenarios' overrides onto one deterministic network."""
    for component, key in CSV_KEYS.items():
        nom_attr = NOM_ATTRS[component]
        static = n.components[component].static

        dfs = {}
        for name, cfg in scenarios_cfg.items():
            path = cfg.get(key)
            if path:
                df = _read_csv(path, clusters)
                _check_names_exist(df, static.index, component, path)
                dfs[name] = df

        if not dfs:
            continue

        touched = pd.Index(sorted(set().union(*(df.index for df in dfs.values()))))
        baseline = static.loc[touched, nom_attr]
        blended = pd.Series(0.0, index=touched)
        for name, cfg in scenarios_cfg.items():
            df = dfs.get(name)
            values = df[nom_attr] if df is not None else pd.Series(dtype=float)
            values = values.reindex(touched).fillna(baseline)
            blended += cfg["probability"] * values

        if component == "Line":
            _set_line_capacity(n, blended)
        else:
            static.loc[touched, nom_attr] = blended
            static.loc[touched, f"{nom_attr}_extendable"] = False
    return n


def _scale_scenario_line_capacity(static, idx, s_nom_new):
    """Scale one scenario's line ``num_parallel`` (-> ``x``/``r``) by capacity.

    Per-scenario counterpart of :func:`_set_line_capacity` for one
    ``(scenario, name)`` slice of the joint stochastic network: scale
    ``num_parallel`` (and so ``x``/``r``) by ``s_nom_new / canvas s_nom`` for
    corridors reduced to a positive capacity and set the new ``s_nom``
    (non-extendable). Corridors delayed to ``s_nom=0`` keep their canvas
    impedance (the ratio is undefined) but are removed from this scenario's
    Kirchhoff cycle at solve time by
    :func:`additional_functionality.add_scenario_kvl`, so they cannot angle-lock
    the grid. Unlike the single-scenario case, lines cannot be dropped per
    scenario, so every line stays in the shared frame.
    """
    s_nom_new = s_nom_new.copy()
    s_nom_new.index = idx
    frac = s_nom_new.where(s_nom_new > 0) / static.loc[idx, "s_nom"]
    live = frac.dropna().index
    static.loc[live, "num_parallel"] *= frac[live]
    static.loc[idx, "s_nom"] = s_nom_new
    static.loc[idx, "s_nom_extendable"] = False


def build_stochastic_topology(n, scenarios_cfg, clusters):
    """Set up the joint two-stage stochastic network (PyPSA-native scenario dimension).

    Each scenario gets its own line impedance (``_scale_scenario_line_capacity``),
    so it must be solved with ``additional_functionality.add_scenario_kvl``, which
    builds one Kirchhoff cycle per scenario from that scenario's reactances and
    drops its ``s_nom=0`` corridors. PyPSA's native scenario KVL would read
    impedance from the first scenario only.
    """
    probabilities = {name: cfg["probability"] for name, cfg in scenarios_cfg.items()}
    n.set_scenarios(probabilities)

    for name, cfg in scenarios_cfg.items():
        for component, key in CSV_KEYS.items():
            path = cfg.get(key)
            if not path:
                continue
            df = _read_csv(path, clusters)
            static = n.components[component].static
            own_names = static.xs(name, level="scenario").index
            _check_names_exist(df, own_names, component, path)
            idx = pd.MultiIndex.from_product([[name], df.index])
            if component == "Line":
                _scale_scenario_line_capacity(static, idx, df["s_nom"])
            else:
                _apply_overrides(static, df.set_axis(idx), component, path)
    return n


def _foreign_countries(n):
    """Real (2-letter) non-DE country codes present in the network - excludes
    DE itself and the empty-string country of EU-level fossil/CO2 buses."""
    return {c for c in n.buses.country.dropna().unique() if c and c != "DE"}


def apply_outside_de_grid(n, spec):
    """
    Fix the transmission grid *outside* Germany to a committed exogen year, or
    leave it endogenously extendable.

    ``spec`` is ``"endogenous"`` (foreign grid stays as prepared, i.e. freely
    optimisable) or ``"<year>_exogen"`` (foreign AC lines and DC links are
    floored to the capacity committed by ``<year>``, i.e. ``build_year <= year``,
    and made non-extendable; not-yet-committed foreign branches are dropped). Only
    foreign-foreign branches are touched - DE and interconnector branches are the
    uncertain grid set by the scenario override and are left to
    build_grid_topology. Must run on the plain network before set_scenarios.
    """
    if spec in (None, "endogenous"):
        return n
    tokens = spec.split("_")
    if len(tokens) != 2 or tokens[1] != "exogen":
        raise ValueError(
            f"outside_de_grid {spec!r} must be 'endogenous' or '<year>_exogen'."
        )
    year = int(tokens[0])

    country = n.buses.country
    foreign = _foreign_countries(n)

    def foreign_mask(static):
        c0, c1 = static.bus0.map(country), static.bus1.map(country)
        return ~c0.eq("DE") & ~c1.eq("DE") & (c0.isin(foreign) | c1.isin(foreign))

    lines = n.components["Line"].static
    fl = lines.index[foreign_mask(lines)]
    s_nom_new = lines.loc[fl, "s_nom"].where(lines.loc[fl, "build_year"] <= year, 0.0)
    _set_line_capacity(n, s_nom_new)

    links = n.components["Link"].static
    fk = links.index[foreign_mask(links) & (links.carrier == "DC")]
    p_nom_new = links.loc[fk, "p_nom"].where(links.loc[fk, "build_year"] <= year, 0.0)
    links.loc[fk, "p_nom"] = p_nom_new
    links.loc[fk, "p_nom_extendable"] = False

    logger.info(
        "outside_de_grid=%s: fixed %d foreign AC line(s) and %d foreign DC link(s) "
        "to the committed grid, non-extendable.",
        spec,
        len(fl),
        len(fk),
    )
    return n


def build_grid_topology(n, grid_scenario, scenarios_cfg, clusters):
    if grid_scenario == "stochastic":
        return build_stochastic_topology(n, scenarios_cfg, clusters)
    elif grid_scenario == "eev":
        return build_eev_topology(n, scenarios_cfg, clusters)
    else:
        return build_deterministic_topology(n, scenarios_cfg[grid_scenario], clusters)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_grid_topology",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2035",
            grid_scenario="2025_exogen",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    grid_scenario = snakemake.params.grid_scenario
    clusters = snakemake.wildcards.clusters

    # Merge each scenario's probability (config) with its resolved CSV paths
    # (params.scenario_csvs, auto-generated by build_grid_scenario_csvs unless
    # overridden by an explicit config path); see rules/pypsa-de/stochastic_grid.smk.
    scenario_csvs = snakemake.params.scenario_csvs
    scenarios_cfg = {
        name: {**cfg, **scenario_csvs.get(name, {})}
        for name, cfg in snakemake.params.stochastic_grid_scenarios["scenarios"].items()
    }

    apply_outside_de_grid(n, snakemake.params.get("outside_de_grid", "endogenous"))

    build_grid_topology(n, grid_scenario, scenarios_cfg, clusters)

    n.export_to_netcdf(snakemake.output.network)
