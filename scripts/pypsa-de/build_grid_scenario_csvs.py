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

The `grid_scenario` name encodes a year, a kind and optional filters,
``<year>_<kind>[_<carrier>][_<direction>]``:

- ``year``      - build-year cutoff (any integer, not just a planning horizon).
- ``exogen``    - committed (NEP/TYNDP/manual) transmission projects only, zero
  endogenous expansion. Taken from the pristine first-horizon prenetwork by
  flooring every not-yet-due branch (``build_year > year``) to 0.
- ``optimal``   - capacities read off the solved postnetwork of ``year``.
- ``carrier``   - optional ``AC``/``DC``: delay only that carrier; the other is
  reset to its fully built-out (canvas) capacity.
- ``direction`` - optional ``NS``/``WE``: delay only branches running mostly
  north-south / west-east (by endpoint orientation); the rest stay built out.

Both qualifiers are optional and order-free. Everything is scoped (default
``de_and_interconnectors``) and validated against the canvas network
build_grid_topology applies the CSV to.
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

NOM_ATTR = {"Line": "s_nom", "Link": "p_nom"}
CARRIERS = {"AC", "DC"}
DIRECTIONS = {"NS", "WE"}


def parse_grid_scenario(name):
    """
    Parse a ``grid_scenario`` name into ``(year, kind, carrier, direction)``.

    Format ``<year>_<kind>[_<carrier>][_<direction>]``; the two qualifiers are
    optional and order-free, and ``None`` means "no restriction on that axis".
    """
    tokens = name.split("_")
    if len(tokens) < 2:
        raise ValueError(
            f"grid_scenario {name!r} must be '<year>_<kind>[_<carrier>][_<direction>]'."
        )
    year = int(tokens[0])
    kind = tokens[1]
    if kind not in ("exogen", "optimal"):
        raise ValueError(
            f"unknown kind {kind!r} in {name!r}; expected 'exogen' or 'optimal'."
        )
    carrier = direction = None
    for tok in tokens[2:]:
        if tok in CARRIERS and carrier is None:
            carrier = tok
        elif tok in DIRECTIONS and direction is None:
            direction = tok
        else:
            raise ValueError(
                f"unknown or duplicate qualifier {tok!r} in {name!r}; expected one "
                f"of {sorted(CARRIERS)} and/or {sorted(DIRECTIONS)}."
            )
    return year, kind, carrier, direction


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


def _is_north_south(source, component):
    """
    True where a branch runs more north-south than west-east.

    Longitude spans are scaled by cos(latitude) so the comparison is in physical
    distance rather than raw degrees.
    """
    df = _branches(source, component)
    y0 = source.buses.y.loc[df.bus0].to_numpy()
    y1 = source.buses.y.loc[df.bus1].to_numpy()
    x0 = source.buses.x.loc[df.bus0].to_numpy()
    x1 = source.buses.x.loc[df.bus1].to_numpy()
    dlat = y1 - y0
    dlon = (x1 - x0) * np.cos(np.deg2rad((y0 + y1) / 2.0))
    return pd.Series(np.abs(dlat) >= np.abs(dlon), index=df.index)


def delay_mask(source, component, carrier, direction):
    """
    Boolean over a component's branches: True where the scenario delays them.

    ``carrier`` (AC/DC) and ``direction`` (NS/WE) narrow the selection; ``None``
    on an axis means no restriction there.
    """
    df = _branches(source, component)
    comp_carrier = "AC" if component == "Line" else "DC"
    if carrier is not None and carrier != comp_carrier:
        return pd.Series(False, index=df.index)
    mask = pd.Series(True, index=df.index)
    if direction is not None:
        ns = _is_north_south(source, component)
        mask &= ns if direction == "NS" else ~ns
    return mask


def scenario_table(source, canvas, component, year, kind, scope, carrier=None, direction=None):
    """
    Scoped, canvas-validated (name -> nom) override table, rounded to 1 dp.

    Only branches selected by the optional ``carrier``/``direction`` filters are
    delayed; every other in-scope branch is reset to the canvas (fully built-out
    target-year) capacity, so the filter narrows *what the scenario delays*.
    """
    nom = NOM_ATTR[component]
    table = branch_table(source, component, year, gate_lines=(component == "Line"))
    table = table.loc[scope_mask(table, scope)]

    missing = table.index.difference(_branches(canvas, component).index)
    if len(missing):
        raise ValueError(
            f"{component} name(s) {list(missing)} not found in the canvas network. "
            "Pristine/solved and canvas networks are out of sync."
        )

    delay = delay_mask(source, component, carrier, direction).reindex(
        table.index, fill_value=False
    )
    undelayed = table.index[~delay]
    if len(undelayed):
        canvas_nom = _branches(canvas, component)[nom].reindex(undelayed)
        table.loc[undelayed, "capacity_mw"] = canvas_nom.to_numpy()

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
    year, kind, carrier, direction = parse_grid_scenario(
        snakemake.wildcards.grid_scenario
    )

    canvas = pypsa.Network(snakemake.input.canvas)
    source = pypsa.Network(
        snakemake.input.solved if kind == "optimal" else snakemake.input.pristine
    )

    for component, out_path in [
        ("Line", snakemake.output.lines),
        ("Link", snakemake.output.links),
    ]:
        table = scenario_table(
            source, canvas, component, year, kind, scope, carrier, direction
        )
        table.to_csv(out_path)
        logger.info(f"wrote {out_path} ({len(table)} rows)")
