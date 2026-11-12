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


def _cutoff_network() -> pypsa.Network:
    """DE and FR buses with projects due 2029 and 2031 on every branch type."""
    n = pypsa.Network()
    n.add("Bus", ["DE1", "DE2"], country="DE")
    n.add("Bus", ["FR1", "FR2"], country="FR")
    for y in (2029, 2031):
        n.add("Line", f"de_{y}", bus0="DE1", bus1="DE2", x=1.0, build_year=y)
        n.add("Line", f"defr_{y}", bus0="DE1", bus1="FR1", x=1.0, build_year=y)
        n.add("Line", f"fr_{y}", bus0="FR1", bus1="FR2", x=1.0, build_year=y)
        n.add("Link", f"dc_{y}", bus0="DE1", bus1="DE2", carrier="DC", build_year=y)
        n.add("Link", f"h2_{y}", bus0="DE1", bus1="DE2", carrier="H2 pipeline", build_year=y)
    n.add("Line", "old_2029", bus0="DE1", bus1="DE2", x=1.0, build_year=1990, lifetime=39)
    return n


@pytest.mark.parametrize(
    "de_cutoff, inactive",
    [
        (2030, ["dc_2031", "de_2031", "defr_2031", "fr_2031", "h2_2031"]),
        (
            2028,
            ["dc_2029", "dc_2031", "de_2029", "de_2031", "defr_2029", "defr_2031", "fr_2031", "h2_2031"],
        ),
    ],
)
def test_deactivate_early_transmission_expansion_de_cutoff(de_cutoff, inactive):
    """The DE cutoff delays only AC lines and DC links touching DE."""
    n = _cutoff_network()
    mp.deactivate_early_transmission_expansion(n, 2030, de_cutoff)
    off = n.lines.index[~n.lines.active].union(n.links.index[~n.links.active])
    assert sorted(off) == inactive


@pytest.mark.parametrize("de_cutoff, retired", [(2030, True), (2028, False)])
def test_replaced_line_follows_de_cutoff(de_cutoff, retired):
    """A DE line replaced by a 2029 upgrade stays in service while the upgrade is delayed."""
    n = _cutoff_network()
    mp.retire_replaced_lines(n, mp.transmission_cutoff(n, n.lines, 2030, de_cutoff))
    assert (not n.lines.at["old_2029", "active"]) == retired
