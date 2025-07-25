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

        if not new.index.symmetric_difference(deci.index).empty:
            logger.error(
                f"Indices of {name} in realization and decision networks do not match. "
                "This may lead to unexpected results."
            )
            assert (
                new.query(f"{attr}_extendable")
                .index.symmetric_difference(deci.query(f"{attr}_extendable").index)
                .empty
            ), (
                f"Indices of {name} with {attr}_extendable in realization and decision networks do not match."
            )
            # raise ValueError("Indices of realization and decision networks do not match.")

        common_i = new.index.intersection(deci.index).difference(
            new.query("carrier == 'BEV charger'").index
        )  # Exclude BEV chargers from modification because the p_nom are taken from UBA

        extendable_i = new.query(f"{attr}_extendable").index

        if not deci.query(f"{attr}_opt > {attr}_max").empty:
            logger.warning(
                f"Decision network have {name} with {attr}_opt > {attr}_max. "
                f"These assets are: {deci.query(f'{attr}_opt > {attr}_max').index.tolist()}"
            )
            _idx = deci.query(f"{attr}_opt > {attr}_max").index
            deci.loc[_idx, attr + "_opt"] = deci.loc[_idx, attr + "_max"]
            logger.warning(
                f"Setting {name} with {attr}_opt > {attr}_max to {attr}_max."
            )

        if not deci.query(f"{attr}_min > {attr}_opt").empty:
            ValueError(
                f"Decision network have {name} with {attr}_min > {attr}_opt. "
                "This may lead to unexpected results."
                f"These assets are: {deci.query(f'{attr}_min > {attr}_opt').index.tolist()}"
            )
            _idx = deci.query(f"{attr}_min > {attr}_opt").index
            deci.loc[_idx, attr + "_opt"] = deci.loc[_idx, attr + "_min"]

        new.loc[common_i, attr] = deci.loc[common_i, attr]
        new.loc[common_i, attr + "_opt"] = deci.loc[common_i, attr + "_opt"]
        new.loc[common_i, attr + "_min"] = deci.loc[common_i, attr + "_min"]
        new.loc[common_i, attr + "_max"] = deci.loc[common_i, attr + "_max"]

        # Ideally nothing should be extendable
        new.loc[extendable_i, attr + "_extendable"] = False

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
                "biogas to gas",
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
            # For fuels allow more, but not less
            _idx = new.loc[new.carrier.isin(fuels)].index.intersection(extendable_i)
            new.loc[_idx, "p_nom_min"] = new.loc[_idx, "p_nom_opt"]

        if name == "stores":
            # there is only one co2 atmosphere store which is always extendable, hence no intersection with extendable_i needed
            _idx = new.query("carrier == 'co2'").index
            new.loc[_idx, "e_nom_extendable"] = True

            # Just making sure that the defaults are set correctly
            assert (new.loc[_idx, "e_nom_min"] == 0).all()
            assert (new.loc[_idx, "e_nom_max"] == np.inf).all()
            assert (new.loc[_idx, "e_nom"] == 0).all()

            # TODO double check if these are all needed
            co2_stores = ["co2 stored", "co2 sequestered"]
            _idx = new.loc[new.carrier.isin(co2_stores)].index.intersection(
                extendable_i
            )
            new.loc[_idx, "e_nom_extendable"] = True
            # TODO this should probably be active -  Allow less, but not more
            # new.loc[_idx, "e_nom_max"] = new.loc[_idx, "e_nom_opt"]

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
            run="AriadneDemand",
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
    # TODO remove this attempt at a hotfix
    n.links.loc[n.links.carrier == "methanolisation", "p_min_pu"] = 0.0
    # n.links.loc[n.links.carrier.str.contains('vent'), 'p_max_pu'] = 1.0
    n_pre = n.copy()
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
