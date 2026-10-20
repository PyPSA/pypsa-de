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
    The 2030 diagonal is beaten by stochastic by 0.03% - within tolerance.
    """
    return _matrix(
        {
            "2025_exogen": [100.0, 130.0, 140.0, 115.0],
            "2030_exogen": [110.0, 100.03, 125.0, 118.0],  # stochastic beats by 0.03%
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


def test_violation_within_tolerance_passes():
    m = _valid_matrix()
    m.loc["eev", "2025_exogen"] = m.loc["2025_exogen", "2025_exogen"] * (1 + 5e-3)
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
