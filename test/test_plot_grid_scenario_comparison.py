# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for scripts/pypsa-de/plot_grid_scenario_comparison.py."""

import importlib.util
import sys

import pypsa
import pytest

sys.path.insert(0, ".")

_spec = importlib.util.spec_from_file_location(
    "plot_grid_scenario_comparison", "scripts/pypsa-de/plot_grid_scenario_comparison.py"
)
pgsc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pgsc)


def test_grid_capacity_volume_ignores_inactive_branches():
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], country="DE")
    n.add("Line", ["old", "new"], bus0="a", bus1="b", x=1.0, length=10.0,
          s_nom_opt=[100.0, 200.0], active=[False, True])
    n.add("Link", ["dc_old", "dc_new"], bus0="a", bus1="b", carrier="DC", length=10.0,
          p_nom_opt=[300.0, 400.0], active=[False, True])

    out = pgsc.grid_capacity_volume(n)

    assert out.loc[("AC line", "DE internal"), "capacity_GW"] == pytest.approx(0.2)
    assert out.loc[("DC link", "DE internal"), "capacity_GW"] == pytest.approx(0.4)


def test_grid_capacity_volume_counts_dc_reversed_twin_once():
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], country="DE")
    n.add("Link", ["dc", "dc-reversed"], bus0=["a", "b"], bus1=["b", "a"],
          carrier="DC", length=10.0, p_nom_opt=400.0)

    out = pgsc.grid_capacity_volume(n)

    assert out.loc[("DC link", "DE internal"), "capacity_GW"] == pytest.approx(0.4)
    assert out.loc[("DC link", "DE internal"), "n"] == 1
