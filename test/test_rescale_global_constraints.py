# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the row scaling of large global constraints for solver numerics:
`rescale_global_constraints` (scripts/_helpers.py) divides a cap row by a
factor before the solve, `unscale_global_constraint_duals` divides the shadow
price back afterwards. The optimum and the shadow prices must be unchanged,
both for a PyPSA-built cap (CO2Limit) and for a custom cap added through
`add_co2_cap_constraint` (co2_limit-DE).
"""

import sys

import pypsa
import pytest

sys.path.insert(0, ".")

from scripts._helpers import (  # noqa: E402
    add_co2_cap_constraint,
    rescale_global_constraints,
    unscale_global_constraint_duals,
)

LIMIT = 3e3  # tCO2; coal-only dispatch emits 4500 t, so the cap binds


def _network(cname, factor):
    """Coal vs gas on one bus; the CO2 cap forces gas in at 20/0.6 EUR/t."""
    n = pypsa.Network()
    n.set_snapshots(range(3))
    n.add("Carrier", "coal", co2_emissions=1.0)
    n.add("Carrier", "gas", co2_emissions=0.4)
    n.add("Bus", "b")
    n.add("Load", "l", bus="b", p_set=[1e3, 2e3, 1.5e3])
    n.add("Generator", "coal", bus="b", carrier="coal", p_nom=1e4, marginal_cost=10)
    n.add("Generator", "gas", bus="b", carrier="gas", p_nom=1e4, marginal_cost=30)
    if cname == "CO2Limit":
        n.add(
            "GlobalConstraint",
            cname,
            type="primary_energy",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=LIMIT,
        )
    cfg = {"enable": factor != 1, "factor": factor, "constraints": [cname]}
    n.config = {"solving": {"options": {"rescale_global_constraints": cfg}}}
    return n


def _solve(cname, factor):
    """Return (objective, raw shadow price, shadow price after unscaling)."""
    n = _network(cname, factor)

    def extra(n, snapshots):
        if cname != "CO2Limit":
            p = n.model["Generator-p"]
            lhs = (p.loc[:, "coal"] + 0.4 * p.loc[:, "gas"]).sum()
            add_co2_cap_constraint(n, cname, lhs, LIMIT)
            n.add("GlobalConstraint", cname, constant=LIMIT, sense="<=", type="")
        rescale_global_constraints(n)

    n.optimize(solver_name="highs", extra_functionality=extra)
    raw = n.global_constraints.at[cname, "mu"]
    unscale_global_constraint_duals(n)
    return n.objective, raw, n.global_constraints.at[cname, "mu"]


@pytest.mark.parametrize("cname", ["CO2Limit", "co2_limit-DE"])
@pytest.mark.parametrize("factor", [1e3, 1e6])
def test_rescaling_keeps_optimum_and_shadow_price(cname, factor):
    obj_ref, _, mu_ref = _solve(cname, 1)
    obj, raw, mu = _solve(cname, factor)
    assert abs(mu_ref) == pytest.approx(100 / 3)
    assert raw == pytest.approx(factor * mu_ref)
    assert obj == pytest.approx(obj_ref, rel=1e-9)
    assert mu == pytest.approx(mu_ref, rel=1e-6)
