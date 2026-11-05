# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for scripts/pypsa-de/modify_prenetwork.py."""

import importlib.util
import sys

import numpy as np
import pypsa
import pytest

sys.path.insert(0, ".")

_spec = importlib.util.spec_from_file_location(
    "modify_prenetwork", "scripts/pypsa-de/modify_prenetwork.py"
)
mp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mp)


@pytest.mark.parametrize(
    "built_until, retired",
    [
        (2024, []),
        (2025, ["old_2025"]),
        (2030, ["old_2025", "old_2030"]),
        (np.inf, ["old_2025", "old_2030", "old_2040"]),
    ],
)
def test_retire_replaced_lines(built_until, retired):
    n = pypsa.Network()
    n.add("Bus", ["a", "b"])
    n.add("Line", "base", bus0="a", bus1="b", x=1.0)
    for end in (2025, 2030, 2040):
        n.add("Line", f"old_{end}", bus0="a", bus1="b", x=1.0, build_year=1990, lifetime=end - 1990)

    mp.retire_replaced_lines(n, built_until)

    assert sorted(n.lines.index[~n.lines.active]) == retired
