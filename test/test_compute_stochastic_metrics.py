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


def test_small_beating_is_reported_not_raised(caplog):
    """eev beats the 2025 diagonal by 0.5%: within the guard, but reported."""
    m = _valid_matrix()
    m.loc["eev", "2025_exogen"] = 100.0 / 1.005
    csm.assert_wait_and_see_optimality(m, GRID_SCENARIO_NAMES)
    row = csm.wait_and_see_table(m, GRID_SCENARIO_NAMES).loc["2025_exogen"]
    assert (row.planned, row.cheapest) == ("2025_exogen", "eev")
    assert row.excess_pct == pytest.approx(0.5)
    assert "not cheapest" in caplog.text


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


def _network(base_extendable: bool) -> pypsa.Network:
    """Load 1 MW then 3 MW; `base` (mc 1, capex 50) and a 10 MW `peak` (mc 30).

    As LT plan `base` settles at 1 MW, so `peak` covers 2 MW at the 3 MW
    hour and sets the price to 30 there.
    """
    n = pypsa.Network()
    n.set_snapshots(range(2))
    n.add("Bus", "b")
    n.add("Load", "d", bus="b", p_set=[1.0, 3.0])
    n.add("Generator", "base", bus="b", marginal_cost=1.0, capital_cost=50.0,
          p_nom=1.0, p_nom_extendable=base_extendable)
    n.add("Generator", "peak", bus="b", marginal_cost=30.0, p_nom=10.0)
    return n


def _lt(probabilities: dict[str, float] | None = None) -> pypsa.Network:
    n = _network(base_extendable=True)
    if probabilities:
        n.set_scenarios(probabilities)
    n.optimize(solver_name="highs", include_objective_constant=False)
    n.model.solver_model = None  # as loaded from disk; get_one_scenario copies `n`
    return n


def _st(shed_cost: float | None) -> pypsa.Network:
    """ST re-dispatch of the LT plan (`base` fixed at 1 MW), optionally with shedding."""
    n = _network(base_extendable=False)
    if shed_cost is not None:
        n.add("Generator", "shed", bus="b", carrier="load", p_nom=10.0,
              marginal_cost=shed_cost)
    n.optimize(solver_name="highs", include_objective_constant=False)
    n.meta = {"co2_slack_used": {}}
    return n


@pytest.mark.parametrize(
    "shed_cost, delta, bound",
    [
        (None, 0.0, 0.0),  # no relaxation: LT and ST must match
        (40.0, 0.0, 0.0),  # shedding dearer than the peak plant stays unused
        (20.0, -20.0, 20.0),  # sheds 2 MW at 20 instead of 30: bound is tight
    ],
)
def test_lt_st_table_bounds_shedding_gain(shed_cost, delta, bound):
    names = ["a", "b"]
    probs = {"a": 0.4, "b": 0.6}
    st = _st(shed_cost)
    topologies = {"a": _lt(), "b": _lt(), "stochastic": _lt(probs)}
    evaluations = {(p, s): st for p in topologies for s in names}

    table = csm.lt_st_table(topologies, evaluations, names, probs)

    assert table.delta.tolist() == pytest.approx([delta] * 3, abs=1e-6)
    assert table.gain_bound.tolist() == pytest.approx([bound] * 3, abs=1e-6)
    csm.assert_lt_st_consistency(table)


def test_shed_snapshots_ignore_solver_residue():
    st = _st(20.0)  # sheds 2 MW in the second hour
    st.generators_t.p.loc[0, "shed"] = 1e-6
    assert csm.shedding_stats(st).shed_snapshots == 1


@pytest.mark.parametrize("slack, unknown", [(1e-4, False), (5.0, True)])
def test_co2_slack_without_lt_price(slack, unknown):
    """Solver residue counts as unused; real slack without LT price gives NaN."""
    st = _st(None)
    st.meta = {
        "co2_slack_used": {"co2_slack-CO2Limit": slack},
        "solving": {"options": {"co2_cap_slack": {"penalty": 1e4}}},
    }
    assert np.isnan(csm.co2_slack_value(_lt(), st)) == unknown


@pytest.mark.parametrize(
    "delta, bound, gap_tol, ok",
    [
        (0.005, 0.0, 0.0, True),  # within rtol
        (0.5, 0.0, 0.0, False),  # ST dearer than LT: loose solve
        (0.5, 0.0, 1.0, True),  # ... but within the reported solver gaps
        (-20.005, 20.0, 0.0, True),  # gain at the bound, within rtol
        (-21.0, 20.0, 0.0, False),  # gain beyond what the relaxations can save
        (-21.0, 20.0, 2.0, True),  # ... but within the reported solver gaps
        (-21.0, np.nan, 0.0, True),  # unknown bound: only ST dearer than LT is checked
        (0.5, np.nan, 0.0, False),
    ],
)
def test_assert_lt_st_consistency(delta, bound, gap_tol, ok):
    table = pd.DataFrame(
        {
            "lt_opex": [100.0],
            "delta": [delta],
            "delta_pct": [delta],
            "gain_bound": [bound],
            "gap_tol": [gap_tol],
        },
        index=["a"],
    )
    if ok:
        csm.assert_lt_st_consistency(table)
    else:
        with pytest.raises(AssertionError, match="dearer than the LT dispatch"):
            csm.assert_lt_st_consistency(table)


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


def test_lt_st_table_gap_tol_sums_reported_gaps():
    names = ["a", "b"]
    probs = {"a": 0.4, "b": 0.6}
    st = _st(None)
    st.meta["relative_gap"] = 1e-3
    topologies = {"a": _lt(), "b": _lt(), "stochastic": _lt(probs)}
    topologies["a"].meta = {"relative_gap": 1e-2}
    evaluations = {(p, s): st for p in topologies for s in names}

    table = csm.lt_st_table(topologies, evaluations, names, probs)

    lt_gap = 1e-2 * abs(topologies["a"].objective)
    assert table.gap_tol["a"] == pytest.approx(lt_gap + 1e-3 * abs(st.objective))
    assert table.gap_tol["stochastic"] == pytest.approx(1e-3 * abs(st.objective))


GAP_CHECK = {
    "lt": {"max_rel_gap": 5e-5, "fail_rel_gap": None},
    "st": {"max_rel_gap": 5e-5, "fail_rel_gap": 1e-3},
}


@pytest.mark.parametrize(
    "lt_gap, st_gap, rough",
    [
        (1e-5, 1e-5, []),
        (4e-4, 1e-5, [("LT", "a", "")]),
        (np.nan, 2e-4, [("ST", "a", "b")]),
    ],
)
def test_solver_gap_table_flags_rough_solves(lt_gap, st_gap, rough):
    topologies = {"a": SimpleNamespace(meta={"relative_gap": lt_gap})}
    evaluations = {
        ("a", "a"): SimpleNamespace(meta={"relative_gap": 1e-5}),
        ("a", "b"): SimpleNamespace(meta={"relative_gap": st_gap}),
    }
    table = csm.solver_gap_table(topologies, evaluations, GAP_CHECK)
    assert table.index[table.verdict == "rough"].tolist() == rough


LT_CELL = ("a", "a")


def test_cost_matrix_lt_cell_matches_st_cell():
    """An LT network standing in for a diagonal ST cell counts opex only, not its capex twice."""
    topologies = {"a": _lt()}
    matrices = [
        csm.build_cost_matrix(topologies, {LT_CELL: n}, ["a"], {"a": 1.0})
        for n in (_st(None), _lt())
    ]
    pd.testing.assert_frame_equal(*matrices)


def test_lt_cells_skip_lt_st_check_and_co2_slack():
    names = ["a", "b"]
    topologies = {"a": _lt(), "b": _lt()}
    evaluations = {(p, s): _st(20.0) for p in topologies for s in names}
    evaluations[LT_CELL] = _lt()

    table = csm.lt_st_table(topologies, evaluations, names, {}, {LT_CELL})
    shedding = csm.shedding_table(evaluations, {LT_CELL})

    assert table.index.tolist() == ["b"]
    assert np.isnan(shedding.loc[LT_CELL, "co2_slack_t"])
    assert shedding.loc[("a", "b"), "co2_slack_t"] == 0


def test_solver_gap_table_judges_lt_cell_by_lt_standard():
    evaluations = {LT_CELL: SimpleNamespace(meta={"relative_gap": 2e-3})}
    table = csm.solver_gap_table({}, evaluations, GAP_CHECK, {LT_CELL})
    assert table.loc[("ST", "a", "a")].tolist() == [2e-3, "rough", "LT"]
