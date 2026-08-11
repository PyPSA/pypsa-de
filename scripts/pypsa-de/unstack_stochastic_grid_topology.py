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
"""

import logging

import pypsa

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def unstack_scenario(n, scenario):
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

    if len(n2.sub_networks):
        n2.remove("SubNetwork", n2.sub_networks.index.tolist())
    n2.determine_network_topology()

    return n2


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "unstack_stochastic_grid_topology",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
            grid_scenario="reduced_40",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    n_scenario = unstack_scenario(n, snakemake.wildcards.grid_scenario)
    n_scenario.export_to_netcdf(snakemake.output.network)
