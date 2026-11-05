# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for scripts/pypsa-de/build_grid_scenario_csvs.py."""

import importlib.util
import sys

import pypsa
import pytest

sys.path.insert(0, ".")

_spec = importlib.util.spec_from_file_location(
    "build_grid_scenario_csvs", "scripts/pypsa-de/build_grid_scenario_csvs.py"
)
csvs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(csvs)


def _network():
    """DE corridor whose line `old` is replaced by the upgrade `new` in 2030."""
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], country="DE", x=[0.0, 1.0], y=[0.0, 0.0])
    n.add("Line", "old", bus0="a", bus1="b", x=1.0, s_nom=100.0, build_year=1990, lifetime=40)
    n.add("Line", "new", bus0="a", bus1="b", x=1.0, s_nom=200.0, build_year=2030)
    return n


@pytest.mark.parametrize("year, old, new", [(2025, 100.0, 0.0), (2029, 100.0, 0.0), (2030, 0.0, 200.0)])
def test_branch_table_retires_replaced_lines(year, old, new):
    table = csvs.branch_table(_network(), "Line", year, gate_lines=True)
    assert table.capacity_mw.to_dict() == {"old": old, "new": new}


def test_undelayed_line_keeps_canvas_retirement():
    canvas = _network()
    canvas.lines.loc["old", "active"] = False  # retired in the canvas year
    table = csvs.scenario_table(_network(), canvas, "Line", 2028, "exogen", "de_internal", carrier="DC")
    assert table.s_nom.to_dict() == {"old": 0.0, "new": 200.0}


def test_branch_table_solved_network_ignores_inactive_lines():
    n = _network()
    n.lines["s_nom_opt"] = n.lines.s_nom  # solved: non-extendable s_nom_opt = s_nom
    n.lines.loc["old", "active"] = False
    table = csvs.branch_table(n, "Line", 2035)
    assert table.capacity_mw.to_dict() == {"old": 0.0, "new": 200.0}
