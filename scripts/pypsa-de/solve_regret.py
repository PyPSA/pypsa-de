import logging

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

logger = logging.getLogger(__name__)


def fix_capacities(realization, decision):
    realization.optimize.fix_optimal_capacities()  # also sets p_nom = p_nom_opt
    decision.optimize.fix_optimal_capacities()

    nominal_attrs = {
        "generators": "p_nom",
        "lines": "s_nom",
        "links": "p_nom",
        "stores": "e_nom",
    }

    realization.links.index.intersection(decision.links.index)

    for name, attr in nominal_attrs.items():
        real = getattr(realization, name)
        deci = getattr(decision, name)

        common = real.index.intersection(deci.index)
        if not real.index.equals(deci.index):
            logger.warning(
                f"Indices of {name} in realization and decision networks do not match. "
                "This may lead to unexpected results."
            )

        real.loc[common, attr] = deci.loc[common, attr]

        if name == "links":
            virtual_links = [
                "oil refining",
                "gas compressing",
                "BEV charger",
                "land transport oil",
                "land transport fuel cell",
                "unsustainable bioliquids",
                "solid biomass for industry",
                "gas for industry",
                "industry methanol",
                "naphtha for industry",
                "process emissions",
                "coal for industry",
                "H2 for industry",
                "shipping methanol",
                "shipping oil",
                "kerosene for aviation",
                "agriculture machinery oil",
                "renewable oil",
                "methanol",
                "renewable gas",
            ]
            real.loc[real.carrier.isin(virtual_links), "p_nom_extendable"] = True
            real.loc[real.carrier.isin(virtual_links), "p_nom_min"] = real.loc[
                real.carrier.isin(virtual_links), "p_nom"
            ]

            real.loc[real.carrier == "SMR", "p_nom_extendable"] = True
            real.loc[real.carrier == "SMR", "p_nom_min"] = real.loc[
                real.carrier == "SMR", "p_nom"
            ]

            real.loc[real.carrier == "waste CHP", "p_nom_extendable"] = True
            real.loc[real.carrier == "waste CHP", "p_nom_min"] = real.loc[
                real.carrier == "waste CHP", "p_nom"
            ]

            real.loc[
                real.carrier == "electricity distribution grid", "p_nom_extendable"
            ] = True  # either this or load shedding?

            real.loc[
                real.carrier == "electricity distribution grid", "p_nom_extendable"
            ] = True  # either this or load shedding?
            real.loc[real.carrier == "electricity distribution grid", "p_nom_min"] = (
                real.loc[real.carrier == "electricity distribution grid", "p_nom"]
            )

        if name == "generators":
            fuels_and_vents = [
                "lignite",
                "coal",
                "oil primary",
                "uranium",
                "gas primary",
                "urban central heat vent",
                "rural heat vent",
                "urban decentral heat vent",
            ]
            real.loc[real.carrier.isin(fuels_and_vents), "p_nom_extendable"] = True

    return realization


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "solve_regret",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
            decision="KN2045_Mix",
            run="LowDemand",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    # if snakemake.input.realization == snakemake.input.decision:
    #     import os
    #     import sys

    #     src = os.path.abspath(snakemake.input.realization)
    #     dst = os.path.abspath(snakemake.output.regret_network)
    #     os.symlink(src, dst)
    #     sys.exit(0)

    logger.info("Loading realization and decision networks")

    import pathlib

    # Touch output file to ensure it exists
    pathlib.Path(snakemake.output.regret_network).touch()

    realization = pypsa.Network(snakemake.input.realization)
    decision = pypsa.Network(snakemake.input.decision)

    planning_horizons = snakemake.wildcards.get("planning_horizons", None)
    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )

    n = fix_capacities(realization, decision)

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
    n.export_to_netcdf(snakemake.output.regret_network)
