# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Solves linear optimal dispatch in hourly resolution using the capacities of
previous capacity expansion in rule :mod:`solve_network`.
"""

import logging

import numpy as np
import pypsa
import logging
import os
import sys
from _benchmark import memory_logger

sys.path.append(os.path.abspath(os.path.dirname(__file__)))  # Adds 'scripts/' to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))  # Adds repo root

from solve_network import prepare_network, solve_network

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_operations_network_myopic",
            simpl="",
            clusters=1,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2025",
            run="KN2045_Bal_v4_24H",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    logger.info(f"Solving again with fixed capacities")

    solve_opts = snakemake.params.options

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)

    # co2 constraint can become infeasible bc of numerical issues (2 ways to handle it)
    
    # (1) round co2 constraint to make feasible (rounding with <= softens the constraint)
    # n.global_constraints.loc["CO2Limit" , "constant"] = round(n.global_constraints.loc["CO2Limit" , "constant"])

    # (2) multiply co2 store e_nom_opt by 2
    n.stores.loc[n.stores.carrier == "co2", "e_nom_opt"] *= 2

    # gas for industry load cannot be satisfied bc of numerical issues
    if  round(n.loads[n.loads.carrier == "gas for industry"].p_set.iloc[0], 5) > n.links[n.links.carrier.isin(["gas for industry" , "gas for industry CC"])].p_nom_opt.sum():
        n.links.loc[n.links.carrier == "gas for industry", "p_nom_opt"] += n.loads[n.loads.carrier == "gas for industry"].p_set.iloc[0]- n.links[n.links.carrier.isin(["gas for industry" , "gas for industry CC"])].p_nom_opt.sum()

    n.optimize.fix_optimal_capacities()

    n = prepare_network(
        n,
        solve_opts=solve_opts,
        config=snakemake.config,
        foresight=snakemake.params.foresight,
        planning_horizons=snakemake.params.planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        snakemake=snakemake,
    )

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        n = solve_network(
            n,
            config=snakemake.config,
            params=snakemake.params,
            solving=snakemake.params.solving,
            log_fn=snakemake.log.solver,
            snakemake=snakemake,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)
    