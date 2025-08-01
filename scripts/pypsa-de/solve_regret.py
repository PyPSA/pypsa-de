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


def fix_capacities(realization, decision):
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
            logger.warning(
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
        extendable_i = new.query(f"{attr}_extendable").index

        # Copy all assets (maybe this should only be done for current planning horizon)
        new.loc[common_i, attr] = deci.loc[common_i, attr]
        new.loc[common_i, attr + "_opt"] = deci.loc[common_i, attr + "_opt"]
        new.loc[common_i, attr + "_min"] = deci.loc[common_i, attr + "_min"]
        new.loc[common_i, attr + "_max"] = deci.loc[common_i, attr + "_max"]

        # Ideally nothing should be extendable
        new.loc[extendable_i, attr + "_extendable"] = False
        # In this case p_nom = p_nom_opt
        new.loc[extendable_i, attr] = new.loc[extendable_i, attr + "_opt"]

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
            ]
            # TODO double check if these are all needed
            essential_links = [
                "SMR",
                "waste CHP",
                # "biogas to gas",
            ]
            # TODO double check if these are all needed
            bottleneck_links = [
                "electricity distribution grid",
                "methanolisation",
                "rural gas boiler",
                "urban decentral gas boiler",
                "urban central gas CHP",
                # For 2035
                "rural biomass boiler",
                "urban decentral biomass boiler",
                "DAC",
            ]
            links_to_free = virtual_links + essential_links + bottleneck_links
            _idx = new.loc[new.carrier.isin(links_to_free)].index.intersection(
                extendable_i
            )
            new.loc[_idx, "p_nom_extendable"] = True

            # For essential links and bottleneck links allow more, but not less
            links_to_limit = essential_links + bottleneck_links
            _idx = new.loc[new.carrier.isin(links_to_limit)].index.intersection(
                extendable_i
            )
            new.loc[_idx, "p_nom_min"] = new.loc[_idx, "p_nom_opt"]

            # For DC fix everything to p_nom_opt
            _idx = new.loc[new.carrier == "DC"].index.intersection(extendable_i)
            new.loc[_idx, "p_nom_min"] = new.loc[_idx, "p_nom_opt"]
            new.loc[_idx, "p_nom_max"] = new.loc[_idx, "p_nom_opt"]
            new.loc[_idx, "p_nom_extendable"] = True

        if name == "lines":
            # For lines fix everything to s_nom_opt
            _idx = new.index.intersection(extendable_i)

            new.loc[_idx, "s_nom_min"] = new.loc[_idx, "s_nom_opt"]
            new.loc[_idx, "s_nom_max"] = new.loc[_idx, "s_nom_opt"]
            new.loc[_idx, "s_nom_extendable"] = True

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
            _idx = new.loc[new.carrier.isin(fuels + vents)].index.intersection(
                extendable_i
            )
            new.loc[_idx, "p_nom_extendable"] = True

        if name == "stores":
            # there is only one co2 atmosphere store which should always be extendable, hence no intersection with extendable_i needed
            _idx = new.query("carrier == 'co2'").index
            new.loc[_idx, "e_nom_extendable"] = True

            # TODO double check if these are all needed
            co2_stores = ["co2 stored", "co2 sequestered"]
            _idx = new.loc[new.carrier.isin(co2_stores)].index.intersection(
                extendable_i
            )
            new.loc[_idx, "e_nom_extendable"] = True

            # Allow less, but not more
            if (
                min(
                    decision.global_constraints.mu.get("co2_sequestration_limit"),
                    realization.global_constraints.mu.get("co2_sequestration_limit"),
                )
                < 1
            ):
                new.loc[_idx, "e_nom_max"] = new.loc[_idx, "e_nom_opt"]

        # Above several assets were switch to extendable again, for these the p_nom value is restored, i.e. set to the value from the decision network
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
            planning_horizons="2030",
            decision="LowDemand",
            run="LowDemand",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

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

    n_pre = n.copy()

    logger.info("Adding CO2 removal service outside DE.")
    n.add("Carrier", "CO2 removal service")

    # pypsa calculates with CO2-tonnes-equivalent not single units of CO2 -> marginal cost in €/tCO2

    n.add(
        "Generator",
        "CO2 removal service",
        bus="co2 atmosphere",
        carrier="CO2 removal service",
        p_nom=1e6,
        marginal_cost=350,
        p_nom_extendable=False,
        p_min_pu=0,
        p_max_pu=1.0,
        sign=-1,
    )

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
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
    n.export_to_netcdf(snakemake.output.regret_network)
