# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests scripts/pypsa-de/compute_stochastic_metrics.py, in particular the
wait-and-see optimality guard on the objective-basis cost matrix: each grid
scenario's own portfolio must be cheapest on its realized grid, and the
stochastic portfolio must be cheapest in expectation. A gross violation (e.g.
freeze_outside_de_capacities dropping frozen foreign capex from portfolio_capex
for every portfolio except the freeze source) must raise.
"""

import importlib.util
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytest

sys.path.insert(0, ".")
sys.path.append("./scripts")

_spec = importlib.util.spec_from_file_location(
    "compute_stochastic_metrics", "scripts/pypsa-de/compute_stochastic_metrics.py"
)
csm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(csm)

GRID_SCENARIO_NAMES = ["2025_exogen", "2030_exogen", "2035_exogen"]
COLUMNS = GRID_SCENARIO_NAMES + ["expected"]
PORTFOLIOS = GRID_SCENARIO_NAMES + ["eev", "stochastic"]


def _matrix(rows: dict[str, list[float]]) -> pd.DataFrame:
    return pd.DataFrame(rows, index=COLUMNS).T.reindex(PORTFOLIOS)


def _valid_matrix() -> pd.DataFrame:
    """
    Each diagonal cell is the column min; stochastic is min of `expected`.
    The 2030 diagonal is beaten by stochastic by 0.005% - within tolerance.
    """
    return _matrix(
        {
            "2025_exogen": [100.0, 130.0, 140.0, 115.0],
            "2030_exogen": [110.0, 100.005, 125.0, 118.0],  # stochastic beats by 0.005%
            "2035_exogen": [120.0, 115.0, 100.0, 116.0],
            "eev": [135.0, 120.0, 118.0, 121.0],
            "stochastic": [101.0, 100.0, 101.0, 99.0],
        }
    )


def test_valid_matrix_passes():
    csm.assert_wait_and_see_optimality(_valid_matrix(), GRID_SCENARIO_NAMES)


def test_matrix_without_stochastic_passes():
    m = _valid_matrix().drop(index="stochastic")
    csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)


def test_inflated_diagonal_raises():
    """
    Freeze-style bug: the 2035_exogen row is inflated on every column, so its
    diagonal is the column max, not the min.
    """
    m = _valid_matrix()
    m.loc["2035_exogen"] = [349.0, 366.0, 349.0, 355.0]
    with pytest.raises(AssertionError, match="2035_exogen"):
        csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)


def test_stochastic_not_expected_cheapest_raises():
    m = _valid_matrix()
    m.loc["stochastic", "expected"] = 130.0
    with pytest.raises(AssertionError, match="expected"):
        csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)


@pytest.mark.parametrize("factor, ok", [(0.5, True), (2.0, False)])
def test_violation_tolerance(factor, ok):
    """eev beats the 2025 diagonal by `factor` x WAIT_AND_SEE_RTOL."""
    m = _valid_matrix()
    diag = m.loc["2025_exogen", "2025_exogen"]
    m.loc["eev", "2025_exogen"] = diag / (1 + factor * csm.WAIT_AND_SEE_RTOL)
    if ok:
        csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)
    else:
        with pytest.raises(AssertionError, match="2025_exogen"):
            csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)


def test_compute_metrics_ws_uses_diagonal():
    """WS is the probability-weighted cost-matrix diagonal (each scenario's own
    portfolio on its own grid), SP/EEV the stochastic/eev `expected` cells."""
    probs = {"2025_exogen": 0.1, "2030_exogen": 0.5, "2035_exogen": 0.4}
    m = _matrix(
        {
            "2025_exogen": [100.0, 130.0, 140.0, 0.0],
            "2030_exogen": [110.0, 101.0, 120.0, 0.0],
            "2035_exogen": [115.0, 112.0, 102.0, 0.0],
            "eev": [118.0, 115.0, 110.0, 150.0],
            "stochastic": [116.0, 111.0, 108.0, 140.0],
        }
    )
    met = csm.compute_metrics(m, GRID_SCENARIO_NAMES, probs)
    ws = 0.1 * 100.0 + 0.5 * 101.0 + 0.4 * 102.0
    assert met["WS"] == pytest.approx(ws)
    assert met["SP"] == pytest.approx(140.0)
    assert met["EEV"] == pytest.approx(150.0)
    assert met["EVPI"] == pytest.approx(140.0 - ws)
    assert met["ECIU"] == pytest.approx(150.0 - 140.0)


def _freeze_network() -> pypsa.Network:
    """One DE bus, two FR buses, with a mix of domestic/foreign and extendable/
    frozen non-grid assets plus the exogenous grid. Only the foreign, non-grid,
    extendable capex (FR generator 20*100 + FR H2 pipeline 3*100 = 2300) counts."""
    n = pypsa.Network()
    n.add("Bus", "DE1", country="DE")
    n.add("Bus", "FR1", country="FR")
    n.add("Bus", "FR2", country="FR")
    n.add("Generator", "g_de", bus="DE1", p_nom_extendable=True, capital_cost=100.0)
    n.add("Generator", "g_fr", bus="FR1", p_nom_extendable=True, capital_cost=100.0)
    n.add("Generator", "g_fr_sunk", bus="FR1", p_nom_extendable=False, capital_cost=100.0)
    n.add("Link", "dc", bus0="DE1", bus1="FR1", carrier="DC", p_nom_extendable=True, capital_cost=100.0)
    n.add("Link", "h2", bus0="FR1", bus1="FR2", carrier="H2 pipeline", p_nom_extendable=True, capital_cost=100.0)
    n.generators.loc[["g_de", "g_fr", "g_fr_sunk"], "p_nom_opt"] = [10.0, 20.0, 5.0]
    n.links.loc[["dc", "h2"], "p_nom_opt"] = [50.0, 3.0]
    return n


def test_foreign_portfolio_capex_counts_only_foreign_extendable_non_grid():
    assert csm.foreign_portfolio_capex(_freeze_network()) == pytest.approx(2300.0)


def test_extra_capex_restores_optimality():
    """A freeze-inflated source row (all foreign capex on the source only) is made
    comparable once every frozen row gets the same foreign-capex offset."""
    m = _valid_matrix()
    m.loc["2035_exogen"] += 138.0  # source carries the foreign buildout, others don't
    with pytest.raises(AssertionError):
        csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)
    for p in ["2025_exogen", "2030_exogen", "eev", "stochastic"]:
        m.loc[p] += 138.0
    csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)


def _lt_network(opex: float) -> pypsa.Network:
    """LT network with portfolio capex 20 (2 MW x 10) and objective capex + opex."""
    n = pypsa.Network()
    n.add("Bus", "b")
    n.add("Generator", "g", bus="b", p_nom_extendable=True, capital_cost=10.0)
    n.generators.loc["g", "p_nom_opt"] = 2.0
    n._objective = 20.0 + opex  # as stored by a solve; the property has no setter
    return n


def _lt_st_inputs(st_opex: dict[tuple[str, str], float]):
    names = ["2025_exogen", "2030_exogen"]
    probs = {"2025_exogen": 0.4, "2030_exogen": 0.6}
    topologies = {s: _lt_network(100.0) for s in names}
    topologies["stochastic"] = _lt_network(100.0)
    evaluations = {
        (p, s): SimpleNamespace(objective=st_opex.get((p, s), 100.0))
        for p in topologies
        for s in names
    }
    return topologies, evaluations, names, probs


@pytest.mark.parametrize(
    "st_opex, ok",
    [
        ({}, True),
        ({("2025_exogen", "2030_exogen"): 150.0}, True),  # off-diagonal is free
        ({("2030_exogen", "2030_exogen"): 100.005}, True),  # within rtol
        ({("2030_exogen", "2030_exogen"): 101.0}, False),  # loose ST solve
        ({("stochastic", "2025_exogen"): 102.5}, False),  # E[ST] 101 vs SP 100
    ],
)
def test_lt_st_consistency(st_opex, ok):
    args = _lt_st_inputs(st_opex)
    if ok:
        csm.assert_lt_st_consistency(*args)
    else:
        with pytest.raises(AssertionError, match="does not reproduce"):
            csm.assert_lt_st_consistency(*args)


@pytest.mark.parametrize(
    "lt_meta, lt_expected",
    [
        ({"relative_gap": 2e-4}, "LT solve <= 2.0e-04 (<= 40.0 M EUR)"),
        ({}, "LT solve not recorded"),
    ],
)
def test_solver_gap_note_reports_largest_gap(lt_meta, lt_expected):
    topologies = {"a": SimpleNamespace(objective=2e11, meta=lt_meta)}
    evaluations = {
        ("a", "a"): SimpleNamespace(objective=2e11, meta={"relative_gap": 1e-5}),
        ("a", "b"): SimpleNamespace(objective=1e11, meta={"relative_gap": 3e-5}),
    }
    note = csm.solver_gap_note(topologies, evaluations)
    assert lt_expected in note
    assert "ST cell <= 3.0e-05 (<= 3.0 M EUR)" in note


def test_solver_gap_note_counts_missing_bounds():
    evaluations = {
        ("a", "a"): SimpleNamespace(objective=2e11, meta={"relative_gap": 1e-5}),
        ("a", "b"): SimpleNamespace(objective=1e11, meta={"relative_gap": np.nan}),
    }
    note = csm.solver_gap_note({}, evaluations)
    assert "ST cell <= 1.0e-05 (<= 2.0 M EUR), no bound for 1" in note
