import logging
import pathlib

import numpy as np
import pypsa

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import prepare_network, solve_network

logger = logging.getLogger(__name__)


def check_matching_components(new, deci, name, attr):
    if not new.index.symmetric_difference(deci.index).empty:
        logger.error(
            f"Indices of {name} in realization and decision networks do not match. "
            "This may lead to unexpected results."
            f"Offending indices are: {new.index.symmetric_difference(deci.index).tolist()}"
        )
        assert (
            new.query(f"{attr}_extendable")
            .index.symmetric_difference(deci.query(f"{attr}_extendable").index)
            .empty
        ), (
            f"Indices of {name} with {attr}_extendable in realization and decision networks do not match."
        )


def _unfix_bottlenecks(new, deci, name, extendable_i):
    if name == "links":
        virtual_links = [
            "land transport oil",
            "land transport fuel cell",
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
            # "co2 sequestered",
        ]

        _idx = new.loc[new.carrier.isin(virtual_links)].index.intersection(extendable_i)
        # Virtual links are free
        new.loc[_idx, "p_nom_extendable"] = True

        bottleneck_links = [
            "SMR",
            "rural gas boiler",
            "urban decentral gas boiler",
            # For 2035
            # "rural biomass boiler",
            # "urban decentral biomass boiler",
        ]
        # Bottleneck links can be extended, but not reduced
        _idx = new.loc[new.carrier.isin(bottleneck_links)].index.intersection(
            extendable_i
        )
        new.loc[_idx, "p_nom_extendable"] = True
        new.loc[_idx, "p_nom_min"] = deci.loc[_idx, "p_nom_opt"]

        # No limits for waste burning outside DE
        _idx = new.query(
            "carrier == 'HVC to air' and not index.str.startswith('DE')"
        ).index.intersection(extendable_i)
        new.loc[_idx, "p_nom_extendable"] = True
        new.loc[_idx, "p_nom_min"] = deci.loc[_idx, "p_nom_opt"]
        # Tight limits inside DE
        _idx = new.query(
            "carrier == 'waste CHP' and index.str.startswith('DE')"
        ).index.intersection(extendable_i)
        new.loc[_idx, "p_nom_extendable"] = True
        new.loc[_idx, "p_nom_min"] = deci.loc[_idx, "p_nom_opt"]
        # new.loc[_idx, "p_nom_max"] = deci.loc[_idx, "p_nom_opt"] + 0.1

        _idx = new.query(
            "carrier == 'electricity distribution grid'"
        ).index.intersection(extendable_i)
        new.loc[_idx, "p_nom_extendable"] = True
        new.loc[_idx, "p_nom_min"] = deci.loc[_idx, "p_nom_opt"]
        # new.loc[_idx, "p_nom_max"] = deci.loc[_idx, "p_nom_opt"] + 0.1

    # if name == "lines":
    #     # For lines allow only a minimal extension, to avoid powerflow constraint issues
    #     new.loc[extendable_i, "s_nom_extendable"] = True

    #     new.loc[extendable_i, "s_nom_min"] = deci.loc[extendable_i, "s_nom_opt"]
    #     new.loc[extendable_i, "s_nom_max"] = deci.loc[extendable_i, "s_nom_opt"] + 10

    if name == "generators":
        fuels = [
            "lignite",
            "coal",
            "oil primary",
            "uranium",
            "gas primary",
        ]
        vents = [
            "urban central heat vent",
            "rural heat vent",
            "urban decentral heat vent",
        ]
        _idx = new.loc[new.carrier.isin(fuels + vents)].index.intersection(extendable_i)
        new.loc[_idx, "p_nom_extendable"] = True

    # if name == "stores":
    # For co2 sequestered no limits are needed since the global constraints is active
    # _idx = new.query("carrier == 'co2 sequestered'").index.intersection(
    #     extendable_i
    # )
    # new.loc[_idx, "e_nom_extendable"] = True
    # new.loc[_idx, "e_nom_max"] = deci.loc[_idx, "e_nom_opt"] + 1

    # _idx = new.query("carrier == 'co2 stored'").index.intersection(
    #     extendable_i
    # )
    # new.loc[_idx, "e_nom_extendable"] = True

    return


def fix_capacities(realization, decision, scope="", strict=False):
    n = realization.copy()

    nominal_attrs = {
        "generators": "p_nom",
        "lines": "s_nom",
        "links": "p_nom",
        "stores": "e_nom",
    }

    for name, attr in nominal_attrs.items():
        new = getattr(n, name)
        deci = getattr(decision, name)

        greater = deci.query(f"{attr}_opt > {attr}_max")
        if not greater.empty:
            logger.error(
                f"Decision network have {name} with {attr}_opt > {attr}_max. "
                f"These assets are: {greater.index.tolist()}"
            )
        smaller = deci.query(f"{attr}_min > {attr}_opt")
        if not smaller.empty:
            logger.error(
                f"Decision network have {name} with {attr}_min > {attr}_opt. "
                "This may lead to unexpected results."
                f"These assets are: {smaller.index.tolist()}"
            )

        check_matching_components(new, deci, name, attr)

        common_i = new.query("carrier != 'BEV charger'").index.intersection(
            deci.query("carrier != 'BEV charger'").index
        )  # Exclude BEV chargers from modification because the p_nom are taken from UBA
        extendable_i = new.query(
            f"{attr}_extendable and index.str.startswith('{scope}')"
        ).index

        # Copy all assets (maybe this should only be done for current planning horizon)
        new.loc[common_i, attr] = deci.loc[common_i, attr]
        new.loc[common_i, attr + "_opt"] = deci.loc[common_i, attr + "_opt"]
        new.loc[common_i, attr + "_min"] = deci.loc[common_i, attr + "_min"]
        new.loc[common_i, attr + "_max"] = deci.loc[common_i, attr + "_max"]

        # Fix everything...
        new.loc[extendable_i, attr + "_extendable"] = False
        new.loc[extendable_i, attr] = new.loc[extendable_i, attr + "_opt"]

        if not strict:
            _unfix_bottlenecks(new, deci, name, extendable_i)

        if name == "stores":
            # there is only one co2 atmosphere store which should always be extendable, hence no intersection with extendable_i needed
            _idx = new.query("carrier == 'co2'").index
            new.loc[_idx, "e_nom_extendable"] = True

        # Above several assets are switched to extendable again, for these the p_nom value is restored, i.e. set to the value from the decision network
        _idx = new.query(f"{attr}_extendable").index.intersection(extendable_i)
        new.loc[_idx, attr] = deci.loc[_idx, attr]

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "solve_regret",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2035",
            decision="LowDemand",
            run="AriadneDemand",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    logger.info("Loading realization and decision networks")

    # Touch output file to ensure it exists
    pathlib.Path(snakemake.output.regret_network).touch()

    realization = pypsa.Network(snakemake.input.realization)
    decision = pypsa.Network(snakemake.input.decision)

    planning_horizons = snakemake.wildcards.get("planning_horizons", None)
    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )

    n = fix_capacities(realization, decision, scope="DE", strict=False)

    n_pre = n.copy()

    snakemake.params.solving["options"]["noisy_costs"] = False

    # TODO remove this again
    # snakemake.params.solving["options"]["load_shedding"] = 100

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
        regret_run=True,  # snakemake.params.get("regret_run", False),
    )

    # n.add(
    #     "Generator",
    #     "co2 atmosphere",
    #     bus="co2 atmosphere",
    #     p_min_pu=-1,
    #     p_max_pu=0,
    #     p_nom_extendable=True,
    #     marginal_cost=realization.global_constraints.loc["CO2Limit", "mu"],
    # )

    # n.global_constraints.drop("CO2Limit", inplace=True)

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

    # logger.info((n.lines.s_nom_opt - decision.lines.s_nom_opt).sort_values())

    logger.info(
        (
            n.links.query("carrier == 'electricity distribution grid'").p_nom_opt
            - decision.links.query(
                "carrier == 'electricity distribution grid'"
            ).p_nom_opt
        ).sort_values()
    )

    logger.info(
        (decision.global_constraints.mu - n.global_constraints.mu).round().sort_values()
    )
