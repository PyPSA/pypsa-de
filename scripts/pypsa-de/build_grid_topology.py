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


def build_deterministic_topology(n, scenario_cfg, clusters):
    """Apply a single scenario's overrides onto a plain deterministic network."""
    for component, key in CSV_KEYS.items():
        path = scenario_cfg.get(key)
        if not path:
            continue
        df = _read_csv(path, clusters)
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

        static.loc[touched, nom_attr] = blended
        static.loc[touched, f"{nom_attr}_extendable"] = False
    return n


def build_stochastic_topology(n, scenarios_cfg, clusters):
    """Set up the joint two-stage stochastic network (PyPSA-native scenario dimension)."""
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
            _apply_overrides(static, df.set_axis(idx), component, path)
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

    scenarios_cfg = snakemake.params.stochastic_grid_scenarios["scenarios"]
    grid_scenario = snakemake.params.grid_scenario
    clusters = snakemake.wildcards.clusters

    build_grid_topology(n, grid_scenario, scenarios_cfg, clusters)

    n.export_to_netcdf(snakemake.output.network)
