# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Extracts one branch of the joint two-stage stochastic grid-topology network
(`topology-stochastic`, see build_grid_topology.py) as a plain, single-
scenario network - e.g. for Ariadne-variable export, which doesn't support
PyPSA's native scenario dimension (a MultiIndex over the grid_scenario
branches).

Network.get_scenario() does this natively, but crashes on SubNetwork: that
component isn't scenario-indexed by n.set_scenarios() to begin with (it's
bookkeeping derived from determine_network_topology(), not solved data), so
it's dropped here and recomputed for the extracted network instead.

Also rescales buses_t.marginal_price by 1/scenario_weight: PyPSA's
stochastic objective weights each scenario's cost by its probability, but
the per-scenario nodal balance constraint is added unweighted, so its dual
(marginal_price) comes out scaled by that same probability - confirmed
intended upstream ("expected shadow price", matching how CAPEX/OPEX are
expected costs too), not a bug. Left as-is, a branch unstacked from e.g. a
50/50 two-scenario network would show ~half the price of an equivalent
deterministic solve. Rescaling here (rather than only where it's plotted)
also keeps it correct for anything else downstream that reads
marginal_price off these networks (e.g. the trade-cost calculations in
plot_grid_scenario_comparison.py).
"""

import logging

import pypsa

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

# Nominal-capacity attribute per component type, used to detect the components
# a scenario branch zeroes out (see _fill_zero_capacity_dispatch).
NOM_ATTRS = {
    "Generator": "p_nom",
    "Line": "s_nom",
    "Link": "p_nom",
    "StorageUnit": "p_nom",
    "Store": "e_nom",
    "Transformer": "s_nom",
}


def _fill_zero_capacity_dispatch(n):
    """Replace all-NaN dispatch time series of zero-capacity components with 0.

    The joint two-stage stochastic solve does not emit per-scenario operational
    values for components a given branch zeroes out (e.g. not-yet-built DC
    projects with ``p_nom == 0`` in the "2025_exogen" grid topology). After the
    ``xs(scenario)`` slice in :func:`unstack_scenario` those columns come out
    all-NaN, which then propagates into anything summing dispatch downstream
    (trade-cost / energy-balance calculations). A component with no capacity
    carries no flow, so filling the missing values with 0 is exact.
    """
    for c in n.components.values():
        nom_attr = NOM_ATTRS.get(c.name)
        if nom_attr is None or c.static.empty:
            continue
        opt_attr = f"{nom_attr}_opt"
        capacity = c.static[opt_attr] if opt_attr in c.static else c.static[nom_attr]
        zero_i = capacity.index[capacity.fillna(0.0) == 0.0]
        if zero_i.empty:
            continue
        for df in c.dynamic.values():
            cols = df.columns.intersection(zero_i)
            if len(cols):
                df.loc[:, cols] = df.loc[:, cols].fillna(0.0)


def unstack_scenario(n, scenario):
    weight = n.scenario_weightings.loc[scenario, "weight"]

    n2 = n.copy()
    n2.name = f"{n2.name} - Scenario '{scenario}'"
    n2._scenarios_data = n2._scenarios_data.iloc[:0]

    for c in n2.components.values():
        if c.name == "SubNetwork":
            continue
        if not c.static.empty:
            c.static = c.static.xs(scenario, level="scenario", axis=0)
        else:
            c.static.index = c.static.index.droplevel("scenario")
        for k, v in c.dynamic.items():
            if not v.empty:
                c.dynamic[k] = v.xs(scenario, level="scenario", axis=1)
            else:
                c.dynamic[k].columns = c.dynamic[k].columns.droplevel("scenario")

    _fill_zero_capacity_dispatch(n2)

    if len(n2.sub_networks):
        n2.remove("SubNetwork", n2.sub_networks.index.tolist())
    n2.determine_network_topology()

    if not n2.buses_t.marginal_price.empty:
        n2.buses_t.marginal_price = n2.buses_t.marginal_price / weight

    return n2


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "unstack_stochastic_grid_topology",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2035",
            grid_scenario="2025_exogen",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    n_scenario = unstack_scenario(n, snakemake.wildcards.grid_scenario)
    n_scenario.export_to_netcdf(snakemake.output.network)
