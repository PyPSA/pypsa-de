import importlib.util
import logging
import pathlib
import re

import numpy as np
import pypsa

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import solve_network

_spec_path = pathlib.Path(__file__).resolve().parent / "modify_prenetwork.py"
_spec = importlib.util.spec_from_file_location(
    "scripts.pypsa_de.modify_prenetwork", _spec_path
)
_modify_prenetwork = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_modify_prenetwork)
remove_flexibility_options = _modify_prenetwork.remove_flexibility_options


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "solve_st_low_res_network",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
            decision="LowDemand",
            run="HighDemand",
            sensitivity="base",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.st_low_res_prenetwork)

    planning_horizons = snakemake.wildcards.get("planning_horizons", None)
    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )
    np.random.seed(snakemake.params.solving["options"].get("seed", 123))

    if "no_flex" in snakemake.params.st_sensitivity:
        logger.info(
            "Running sensitivity of the short term model with less flexibility options."
        )
        remove_flexibility_options(n)

    gas_price = re.findall(r"gas_price_(\d{2,3})", snakemake.params.st_sensitivity)
    if gas_price:
        gas_price = int(gas_price[0])
        logger.info(
            f"Running sensitivity of the short term model with gas price set to {gas_price} €/MWh."
        )
        n.generators.loc[n.generators.carrier == "gas primary", "marginal_cost"] = (
            gas_price
        )

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        solve_network(
            n,
            config=snakemake.config,
            params=snakemake.params,
            solving=snakemake.params.solving,
            planning_horizons=planning_horizons,
            rule_name=snakemake.rule,
            log_fn=snakemake.log.solver,
            snakemake=snakemake,
        )
    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    n.export_to_netcdf(snakemake.output.st_low_res_network)
