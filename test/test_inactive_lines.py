# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests that inactive AC lines stay out of the Kirchhoff cycles in
`scripts/solve_network.py`.

PyPSA builds the cycle basis from all lines and only drops the flows of
inactive ones, so a cycle through an inactive line becomes a phantom
constraint. `detach_inactive_lines` removes them for the solve and
`reattach_lines` restores them for the exported network.
"""

import sys

import pandas as pd
import pypsa
import pytest

sys.path.insert(0, ".")

from scripts.solve_network import detach_inactive_lines, reattach_lines  # noqa: E402


def _triangle(inactive=True):
    """
    Triangle A-B-C with a cheap generator at A, an expensive one at B and load
    at C. The line A-B is inactive (or absent), so A and B each have a radial
    path to C and the cheap generator serves the whole load.
    """
    n = pypsa.Network()
    n.set_snapshots(range(2))
    n.add("Bus", ["A", "B", "C"], v_nom=380.0)
    n.add("Generator", "cheap", bus="A", p_nom=200, marginal_cost=10)
    n.add("Generator", "dear", bus="B", p_nom=200, marginal_cost=100)
    n.add("Load", "load", bus="C", p_set=100)
    lines = {"AC": ("A", "C"), "BC": ("B", "C")}
    if inactive:
        lines["AB"] = ("A", "B")
    for name, (b0, b1) in lines.items():
        n.add("Line", name, bus0=b0, bus1=b1, x=0.1, r=0.01, s_nom=500)
    if inactive:
        n.lines.loc["AB", "active"] = False
        n.lines_t.s_max_pu = pd.DataFrame({"AB": [0.7, 0.8]}, index=n.snapshots)
    return n


def _solve(n):
    n.optimize(solver_name="highs")
    return n.objective


def test_inactive_line_creates_phantom_cycle_without_detach():
    assert _solve(_triangle()) > _solve(_triangle(inactive=False)) + 1


def test_detached_solve_matches_network_without_the_line():
    n = _triangle()
    removed = detach_inactive_lines(n)
    assert _solve(n) == pytest.approx(_solve(_triangle(inactive=False)))
    reattach_lines(n, *removed)


def test_reattach_restores_lines_with_zero_outputs():
    n = _triangle()
    before = n.lines.copy()
    removed = detach_inactive_lines(n)
    _solve(n)
    reattach_lines(n, *removed)

    assert set(n.lines.index) == set(before.index)
    cols = before.columns.drop(["s_nom_opt", "sub_network"])
    pd.testing.assert_series_equal(n.lines.loc["AB", cols], before.loc["AB", cols])
    assert n.lines.at["AB", "s_nom_opt"] == n.lines.at["AB", "s_nom"]
    assert (n.lines_t.p0["AB"] == 0).all()
    assert n.lines_t.s_max_pu["AB"].tolist() == [0.7, 0.8]


def test_no_inactive_lines_is_a_no_op():
    n = _triangle(inactive=False)
    before = n.lines.copy()
    reattach_lines(n, *detach_inactive_lines(n))
    pd.testing.assert_frame_equal(n.lines, before)
