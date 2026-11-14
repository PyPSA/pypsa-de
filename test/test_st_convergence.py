# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the convergence check of the ST grid-dispatch solves
(`solve_with_fallback`, scripts/pypsa-de/evaluate_grid_portfolio.py).

Barrier without crossover can end 'optimal' while gurobi reports no finite
`ObjBound`. Such a solve is accepted with a NaN gap.
"""

import importlib.util
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, ".")
sys.path.append("./scripts")

_spec = importlib.util.spec_from_file_location(
    "evaluate_grid_portfolio", "scripts/pypsa-de/evaluate_grid_portfolio.py"
)
egp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(egp)

MAX_REL_GAP = 5e-5
OBJ = 2e11


def _network(bound: float, dual_vio: float, condition: str) -> SimpleNamespace:
    sm = SimpleNamespace(
        ObjVal=OBJ,
        ObjBound=bound,
        DualVio=dual_vio,
    )
    return SimpleNamespace(
        model=SimpleNamespace(solver_model=sm),
        optimize=SimpleNamespace(solve_model=lambda **kw: ("ok", condition)),
    )


@pytest.mark.parametrize(
    "bound, dual_vio, condition, rejected",
    [
        (OBJ * (1 - 1e-6), 1e-8, "optimal", False),
        (OBJ * (1 - 1e-3), 1e-8, "optimal", True),
        (-np.inf, 4e-6, "optimal", False),
        (OBJ, 1e-8, "suboptimal", True),
    ],
)
def test_solve_with_fallback_rejects_unconverged(bound, dual_vio, condition, rejected):
    n = _network(bound, dual_vio, condition)
    _, _, gap = egp.solve_with_fallback(n, {}, [], MAX_REL_GAP)
    assert (gap > MAX_REL_GAP) == rejected


def test_missing_bound_gives_nan_gap():
    n = _network(-np.inf, 4e-6, "optimal")
    _, _, gap = egp.solve_with_fallback(n, {}, [], MAX_REL_GAP)
    assert np.isnan(gap)
