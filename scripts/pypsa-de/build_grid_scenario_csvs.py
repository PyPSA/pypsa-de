# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Generate a grid-topology scenario's capacity-override CSVs (one for AC lines,
one for DC links) directly from the networks the workflow already builds, so
they never go stale when a clustering or base-network change renumbers
components. This reproduces the `*_exogen`/`*_optimal` export logic of
notebooks/grid_topology.ipynb as a workflow rule (see build_grid_topology.py,
which consumes the output).

The `grid_scenario` name encodes a year and a kind, `<year>_<kind>`:

- ``exogen``  - committed (NEP/TYNDP/manual) transmission projects only, zero
  endogenous expansion. Taken from the pristine first-horizon prenetwork by
  flooring every not-yet-due branch (``build_year > year``) to 0.
- ``optimal`` - capacities read off the solved postnetwork of ``year``.

Both are scoped (default ``de_and_interconnectors``) and validated against the
canvas network build_grid_topology applies the CSV to.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

NOM_ATTR = {"Line": "s_nom", "Link": "p_nom"}


def _branches(n, component):
    return n.lines if component == "Line" else n.links.query("carrier == 'DC'")


def branch_table(n, component, horizon, gate_lines=False):
    """
    One row per branch: effective capacity, build_year, endpoint countries.

    Capacity is the solved value if the network is solved (any ``*_nom_opt >
    0``), else the nominal value. For an unsolved network, not-yet-due branches
    (``build_year > horizon``) are floored to 0 - always for DC Links, and for
    AC Lines only when ``gate_lines`` is set (matching grid_topology.ipynb).
    """
    nom = NOM_ATTR[component]
    df = _branches(n, component)
    opt = df[f"{nom}_opt"]
    if (opt > 0).any():
        capacity = opt
    elif component == "Link" or gate_lines:
        not_yet_due = df.build_year > horizon
        capacity = df[nom].where(~not_yet_due, 0.0)
    else:
        capacity = df[nom]
    return pd.DataFrame(
        {
            "capacity_mw": capacity,
            "country0": n.buses.loc[df.bus0, "country"].to_numpy(),
            "country1": n.buses.loc[df.bus1, "country"].to_numpy(),
        },
        index=df.index,
    )


def scope_mask(table, scope):
    de0, de1 = table.country0 == "DE", table.country1 == "DE"
    if scope == "de_internal":
        return de0 & de1
    if scope == "de_and_interconnectors":
        return de0 | de1
    if scope == "whole_network":
        return pd.Series(True, index=table.index)
    raise ValueError(f"unknown scope {scope!r}")


def scenario_table(source, canvas, component, year, kind, scope):
    """Scoped, canvas-validated (name -> nom) override table, rounded to 1 dp."""
    nom = NOM_ATTR[component]
    table = branch_table(source, component, year, gate_lines=(component == "Line"))
    table = table.loc[scope_mask(table, scope)]
    missing = table.index.difference(_branches(canvas, component).index)
    if len(missing):
        raise ValueError(
            f"{component} name(s) {list(missing)} not found in the canvas network. "
            "Pristine/solved and canvas networks are out of sync."
        )
    out = table[["capacity_mw"]].rename(columns={"capacity_mw": nom})
    out.index.name = "name"
    return out.round(1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_grid_scenario_csvs",
            clusters="49",
            opts="",
            sector_opts="none",
            planning_horizons="2035",
            grid_scenario="2025_exogen",
        )

    configure_logging(snakemake)

    scope = snakemake.params.scope
    year_str, kind = snakemake.wildcards.grid_scenario.rsplit("_", 1)
    year = int(year_str)

    canvas = pypsa.Network(snakemake.input.canvas)
    if kind == "optimal":
        source = pypsa.Network(snakemake.input.solved)
    elif kind == "exogen":
        source = pypsa.Network(snakemake.input.pristine)
    else:
        raise ValueError(
            f"Unknown grid-scenario kind {kind!r} in "
            f"{snakemake.wildcards.grid_scenario!r}; expected 'exogen' or 'optimal'."
        )

    for component, out_path in [
        ("Line", snakemake.output.lines),
        ("Link", snakemake.output.links),
    ]:
        table = scenario_table(source, canvas, component, year, kind, scope)
        table.to_csv(out_path)
        logger.info(f"wrote {out_path} ({len(table)} rows)")
