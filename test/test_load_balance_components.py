# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the load-shedding generators of `add_load_balance_components`
(scripts/solve_network.py), also on stochastic networks (`n.set_scenarios`).
"""

import sys

import pypsa
import pytest

sys.path.insert(0, ".")
sys.path.append("./scripts")

from scripts.solve_network import add_load_balance_components  # noqa: E402


@pytest.mark.parametrize("stochastic", [False, True])
@pytest.mark.parametrize(
    ("all_carriers", "expected"),
    [
        (False, {"AC bus load shedding": 4000}),
        (True, {"AC bus load shedding": 4000, "gas bus load shedding": 1e5}),
    ],
)
def test_load_shedding_generators(stochastic, all_carriers, expected):
    n = pypsa.Network()
    n.add("Bus", "AC bus", carrier="AC")
    n.add("Bus", "gas bus", carrier="gas")
    if stochastic:
        n.set_scenarios({"s1": 0.4, "s2": 0.6})
    cfg = {"default_cost": 1e5, "all_carriers": all_carriers, "carriers": {"AC": 4000}}

    add_load_balance_components(n, cfg)

    gens = n.generators.query("carrier == 'load'")
    if stochastic:
        assert set(gens.index.unique("scenario")) == {"s1", "s2"}
        gens = gens.xs("s1", level="scenario")
    assert gens.marginal_cost.to_dict() == expected
