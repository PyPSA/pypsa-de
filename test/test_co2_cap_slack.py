# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the elastic CO2-cap slack used to stabilise the ST grid-dispatch solves:
`add_co2_cap_constraint` (scripts/_helpers.py) adds a hard cap unless the slack
is enabled for that constraint, and the penalized slack stays at zero whenever
the cap is reachable, so `check_co2_slack`
(scripts/pypsa-de/evaluate_grid_portfolio.py) passes. A cap the frozen
portfolio cannot meet forces the slack open and must trip the guard.
"""

import importlib.util
import sys
from types import SimpleNamespace

import linopy
import pytest

sys.path.insert(0, ".")
sys.path.append("./scripts")

from scripts._helpers import add_co2_cap_constraint  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "evaluate_grid_portfolio", "scripts/pypsa-de/evaluate_grid_portfolio.py"
)
egp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(egp)


def _network(enable, floor):
    """Single-variable linopy model wrapped as a minimal `n` for the helpers.

    `x >= floor` is the forced emission; the cap is `x <= 8` (plus slack). A
    floor above the cap can only be met by opening the slack.
    """
    m = linopy.Model()
    x = m.add_variables(lower=0, name="x")
    m.add_constraints(x >= floor, name="floor")
    m.objective = 1.0 * x  # minimise emissions down onto the cap
    cfg = {
        "enable": enable,
        "penalty": 10000,
        "max_overshoot": 0.02,
        "assert_threshold": 0.005,
        "constraints": ["CO2Limit"],
    }
    n = SimpleNamespace(
        model=m,
        has_scenarios=False,
        config={"solving": {"options": {"co2_cap_slack": cfg}}},
    )
    return n, x


def _solve_cap(n, x, limit=8.0):
    add_co2_cap_constraint(n, "CO2Limit", 1.0 * x, limit)
    n.model.solve(solver_name="highs")


def test_hard_constraint_when_disabled():
    n, x = _network(enable=False, floor=5.0)
    _solve_cap(n, x)
    assert "co2_slack-CO2Limit" not in n.model.variables
    assert egp.check_co2_slack(n) == {}


def test_slack_stays_zero_when_cap_reachable():
    n, x = _network(enable=True, floor=5.0)
    _solve_cap(n, x)
    assert "co2_slack-CO2Limit" in n.model.variables
    used = egp.check_co2_slack(n)
    assert used == pytest.approx({"co2_slack-CO2Limit": 0.0})


def test_guard_raises_when_cap_unreachable():
    # floor 8.1 vs cap 8 forces slack 0.1: above assert_threshold (0.04) yet
    # within the max_overshoot bound (0.16), so the solve stays feasible.
    n, x = _network(enable=True, floor=8.1)
    _solve_cap(n, x)
    assert n.model["co2_slack-CO2Limit"].solution.sum() == pytest.approx(0.1)
    with pytest.raises(AssertionError, match="CO2 cap overshoot"):
        egp.check_co2_slack(n)


class _FakeNet:
    """Records the solver_options of each solve and replays scripted results."""

    def __init__(self, results):
        self._results = iter(results)
        self.calls = []
        self.optimize = self

    def solve_model(self, **kwargs):
        self.calls.append(kwargs["solver_options"])
        return next(self._results)


_KW = {"solver_name": "gurobi", "solver_options": {"threads": 16}}


def test_fallback_stops_on_first_optimal():
    n = _FakeNet([("ok", "optimal")])
    status, condition = egp.solve_with_fallback(n, _KW, [{"BarHomogeneous": 1}])
    assert (status, condition) == ("ok", "optimal")
    assert n.calls == [{"threads": 16}]  # tier 1 only, no robust knobs paid


def test_fallback_escalates_then_succeeds():
    n = _FakeNet([("ok", "suboptimal"), ("ok", "optimal")])
    status, condition = egp.solve_with_fallback(n, _KW, [{"BarHomogeneous": 1}])
    assert (status, condition) == ("ok", "optimal")
    assert n.calls == [{"threads": 16}, {"threads": 16, "BarHomogeneous": 1}]


def test_fallback_exhausted_returns_last_nonoptimal():
    ladder = [{"BarHomogeneous": 1}, {"BarHomogeneous": 1, "crossover": 1}]
    n = _FakeNet([("ok", "suboptimal")] * 3)
    status, condition = egp.solve_with_fallback(n, _KW, ladder)
    assert condition == "suboptimal"
    assert len(n.calls) == 3  # tier 1 + two fallback tiers
