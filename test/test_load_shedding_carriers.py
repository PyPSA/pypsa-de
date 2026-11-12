# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests which bus carriers get load shedding in the ST grid-dispatch solves
(`with_load_shedding_carriers`, scripts/pypsa-de/evaluate_grid_portfolio.py).
"""

import importlib.util
import sys

import pypsa
import pytest

sys.path.insert(0, ".")
sys.path.append("./scripts")

_spec = importlib.util.spec_from_file_location(
    "evaluate_grid_portfolio", "scripts/pypsa-de/evaluate_grid_portfolio.py"
)
egp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(egp)


@pytest.fixture
def network():
    n = pypsa.Network()
    for carrier in [
        "low voltage",
        "urban central heat",
        "coal for industry",
        "process emissions",
        "gas",
    ]:
        n.add("Bus", carrier, carrier=carrier)
    for bus in [
        "low voltage",
        "urban central heat",
        "coal for industry",
        "process emissions",
    ]:
        n.add("Load", bus, bus=bus, p_set=1.0)
    return n


@pytest.mark.parametrize(
    ("all_carriers", "expected"),
    [
        (
            True,
            {"low voltage": 4000, "urban central heat": 1e5, "coal for industry": 1e5},
        ),
        (False, {"low voltage": 4000}),
    ],
)
def test_load_shedding_carriers(network, all_carriers, expected):
    opts = {
        "load_shedding": {
            "enable": True,
            "default_cost": 1e5,
            "all_carriers": all_carriers,
            "carriers": {"low voltage": 4000, "AC": 4000},
        }
    }
    shedding = egp.with_load_shedding_carriers(opts, network)["load_shedding"]
    assert shedding["carriers"] == expected
    assert shedding["all_carriers"] is False
