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


def _unfix_bottlenecks(new, deci, name, extendable_i):
    if name == "links":
        # Links that have 0-cost and are extendable
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
            "co2 sequestered",
        ]

        _idx = new.loc[new.carrier.isin(virtual_links)].index.intersection(extendable_i)
        new.loc[_idx, "p_nom_extendable"] = True

        # Bottleneck links can be extended, but not reduced to fix infeasibilities due to numerical inconsistencies
        bottleneck_links = [
            "electricity distribution grid",
            "waste CHP",
            "SMR",
            # Boilers create bottlenecks AND should be extendable for fixed_profile_scaling constraints to be applied correctly
            "rural gas boiler",
            "urban decentral gas boiler",
            # Biomass for 2035 when gas is banned
            "rural biomass boiler",
            "urban decentral biomass boiler",
        ]
        _idx = new.loc[new.carrier.isin(bottleneck_links)].index.intersection(
            extendable_i
        )
        new.loc[_idx, "p_nom_extendable"] = True
        new.loc[_idx, "p_nom_min"] = deci.loc[_idx, "p_nom_opt"]

        # Waste outside DE can also be burned directly
        _idx = new.query(
            "carrier == 'HVC to air' and not index.str.startswith('DE')"
        ).index.intersection(extendable_i)
        new.loc[_idx, "p_nom_extendable"] = True
        new.loc[_idx, "p_nom_min"] = deci.loc[_idx, "p_nom_opt"]

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

    return


def fix_capacities(realization, decision, scope="DE", strict=False):
    logger.info(f"Fixing all capacities for scope: {scope}")
    if scope == "EU":
        scope = ""
    if not strict:
        logger.info("Freeing virtual links, bottlenecks and fossil generators.")

    # Copy all existing assets from the decision network
    n = decision.copy()

    # The constraints and loads are taken from the realization network
    n.global_constraints = realization.global_constraints.copy()
    n.loads = realization.loads.copy()

    nominal_attrs = {
        "generators": "p_nom",
        "lines": "s_nom",
        "links": "p_nom",
        "stores": "e_nom",
    }

    for name, attr in nominal_attrs.items():
        new = getattr(n, name)
        deci = getattr(decision, name)
        real = getattr(realization, name)

        # Scenario specific assets are taken from the realization network
        _idx = real.query("carrier in ['BEV charger', 'V2G', 'EV battery']").index
        new.loc[_idx, attr] = real.loc[_idx, attr]

        # Start with fixing everything...
        extendable_i = new.query(
            f"{attr}_extendable and index.str.startswith('{scope}')"
        ).index
        new.loc[extendable_i, attr + "_extendable"] = False
        new.loc[extendable_i, attr] = new.loc[extendable_i, attr + "_opt"]

        # Some links should be extendable to avoid infeasibilities or allow burning more fossil fuels
        if not strict:
            _unfix_bottlenecks(new, deci, name, extendable_i)

        # The CO2 constraints on atmosphere and sequestration need extendable stores to work correctly
        if name == "stores":
            logger.info("Freeing co2 atmosphere and sequestered stores.")
            # there is only one co2 atmosphere store which should always be extendable, hence no intersection with extendable_i needed
            _idx = new.query("carrier == 'co2'").index
            new.loc[_idx, "e_nom_extendable"] = True
            # co2 sequestered stores from previous planning horizons should not be extendable
            _idx = new.query("carrier == 'co2 sequestered'").index.intersection(
                extendable_i
            )
            new.loc[_idx, "e_nom_extendable"] = True

        # Above several assets are switched to extendable again, for these the p_nom value is restored to the value from the decision network

        _idx = new.query(f"{attr}_extendable")

        if not _idx.difference(extendable_i).empty:
            raise ValueError(
                "Assets that are not extendable in the decision network have been set to extendable. This should not happen. Aborting."
            )

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
            decision="AriadneDemand",
            run="AriadneDemand",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Touch output file to ensure it exists
    pathlib.Path(snakemake.output.regret_network).touch()

    realization = pypsa.Network(snakemake.input.realization)
    decision = pypsa.Network(snakemake.input.decision)

    planning_horizons = snakemake.wildcards.get("planning_horizons", None)
    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )
    solve_opts = snakemake.params.solving["options"]
    assert solve_opts["noisy_costs"] == False, (
        "Noisy costs should not be used in regret runs."
    )
    np.random.seed(solve_opts.get("seed", 123))

    n = fix_capacities(
        realization, decision, scope=snakemake.params.scope_to_fix, strict=False
    )

    if solve_opts["post_discretization"].get("enable") and not solve_opts.get(
        "skip_iterations"
    ):
        # Undo the last lines of optimize transmission expansion iteratively
        n.lines.s_nom_extendable = False
        n.lines.s_nom = n.lines.s_nom_opt

        discretized_links = n.links.query(
            f"carrier in {list(solve_opts['post_discretization'].get('link_unit_size').keys())}"
        ).index
        n.links.loc[discretized_links, "p_nom_extendable"] = False
        n.links.loc[discretized_links, "p_nom"] = n.links.loc[
            discretized_links, "p_nom_opt"
        ]

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
        regret_run=True,
    )

    if snakemake.params.scope_to_fix == "EU":
        logger.info(
            f"Fixing Scope 'EU' chosen. Setting the CO2 price to the price from the realization network to avoid infeasibilities: {realization.global_constraints.loc['CO2Limit', 'mu']} €/t_CO2"
        )
        logger.warning(
            "Please make sure that the long-term run with unchanged demand is consistent with the short-term run."
        )
        n.add(
            "Generator",
            "co2 atmosphere",
            bus="co2 atmosphere",
            p_min_pu=-1,
            p_max_pu=0,
            p_nom_extendable=True,
            carrier="co2",
            marginal_cost=realization.global_constraints.loc["CO2Limit", "mu"],
        )
        n.global_constraints.drop("CO2Limit", inplace=True)

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

    constraint_diff = (
        (decision.global_constraints.mu - n.global_constraints.mu)
        .round(2)
        .sort_values()
    )

    logger.info(
        "Difference in global constraints (decision - regret_network): %s",
        constraint_diff,
    )

    if snakemake.input.realization == snakemake.input.decision:
        if abs(constraint_diff["CO2Limit"]) > 1:
            logger.error(
                "Difference in CO2 price between long-term and short-term model is too high: %s",
                constraint_diff["CO2Limit"],
            )

    n.export_to_netcdf(snakemake.output.regret_network)
