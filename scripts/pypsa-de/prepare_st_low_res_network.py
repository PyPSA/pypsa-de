import importlib.util
import logging
import pathlib

import numpy as np
import pypsa

from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import prepare_network


def _load_attr_from_file(filename: str, attr_name: str) -> object:
    """
    Load attribute attr_name from a local python file given by filename (including '.py').
    """
    if not filename.endswith(".py"):
        raise ValueError("filename must include the '.py' extension")
    module_stem = pathlib.Path(filename).stem
    _spec_path = pathlib.Path(__file__).resolve().parent / filename
    _spec = importlib.util.spec_from_file_location(
        f"scripts.pypsa_de.{module_stem}", _spec_path
    )
    _mod = importlib.util.module_from_spec(_spec)
    assert _spec is not None and _spec.loader is not None
    _spec.loader.exec_module(_mod)
    return getattr(_mod, attr_name)


_unfix_bottlenecks = _load_attr_from_file(
    "prepare_regret_network.py", "_unfix_bottlenecks"
)
remove_flexibility_options = _load_attr_from_file(
    "modify_prenetwork.py", "remove_flexibility_options"
)


logger = logging.getLogger(__name__)

eeg_targets = {
    2030: {
        "solar": 215000,
        "onwind": 115000,
        "offwind": 30000,
    }
}

co2_prices = {
    2030: 220,
}


def scale_new_res_to_target(n, targets, year, ratio=1.0):
    for tech, target in targets[year].items():
        logger.info(
            f"Scaling installed capacity of {tech} in DE to target of {target * ratio} MW."
        )
        gens = n.generators.query(
            f"carrier.str.contains('{tech}') and not carrier.str.contains('solar thermal') and index.str.startswith('DE')"
        )
        existing_gens = gens[gens.build_year < year]
        new_gens = gens[gens.build_year == year]
        existing_cap = existing_gens.p_nom_opt.sum()
        new_cap = new_gens.p_nom_opt.sum()
        assert np.isclose(existing_cap + new_cap, gens.p_nom_opt.sum())

        scale_factor = (target * ratio - existing_cap) / new_cap

        n.generators.loc[new_gens.index, "p_nom_opt"] *= scale_factor
        n.generators.loc[new_gens.index, "p_nom"] *= scale_factor


def fix_capacities(n_lt, no_flex=False):
    n = n_lt.copy()

    nominal_attrs = {
        "generators": "p_nom",
        "lines": "s_nom",
        "links": "p_nom",
        "stores": "e_nom",
    }

    for name, attr in nominal_attrs.items():
        new = getattr(n, name)
        lt = getattr(n_lt, name)

        extendable_i = new.query(f"{attr}_extendable").index

        new.loc[extendable_i, attr + "_extendable"] = False
        new.loc[extendable_i, attr] = new.loc[extendable_i, attr + "_opt"]

        _unfix_bottlenecks(new, lt, name, extendable_i)

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

        _idx = new.query(f"{attr}_extendable").index

        new.loc[_idx, attr] = lt.loc[_idx, attr]

    if no_flex:
        logger.info("Realization network is from a run without flexibility.")
        remove_flexibility_options(n)
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "prepare_st_low_res_network",
            clusters=27,
            opts="",
            sector_opts="none",
            st_years="2030",
            run="HighDemand",
        )

    configure_logging(snakemake)

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n_lt = pypsa.Network(snakemake.input.network)
    st_years = snakemake.wildcards.get("st_years", None)
    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )
    solve_opts = snakemake.params.solving["options"]
    assert solve_opts["noisy_costs"] == False, (
        "Noisy costs should not be used in regret runs."
    )
    np.random.seed(solve_opts.get("seed", 123))

    strict = snakemake.params["strict"]
    scope_to_fix = snakemake.params["scope_to_fix"]
    h2_vent = snakemake.params["h2_vent"]

    n = fix_capacities(n_lt, snakemake.params.get("no_flex_lt_run", False))

    scale_new_res_to_target(n, eeg_targets, int(st_years), ratio=1.0)

    if h2_vent:
        logger.info("H2 venting activated for short-term run.")
        n.add("Carrier", "H2 vent", color="#dd2e23", nice_name="H2 vent")

        n.add(
            "Generator",
            n.buses.query("carrier=='H2'").index,
            " vent",
            bus=n.buses.query("carrier=='H2'").index,
            carrier="H2 vent",
            sign=-1,
            marginal_cost=1,
            p_nom=1e6,
        )

    # Manipulating the global constraints
    to_keep = [
        "biomass limit",
        "unsustainable biomass limit",
        "co2_sequestration_limit",
        "CO2Limit",
        "co2_limit-DE",
    ]

    to_drop = n.global_constraints.index.difference(to_keep)
    logger.info("Short-term run detected. Dropping the following constraints:")
    logger.info(to_drop)

    n.global_constraints.drop(to_drop, inplace=True)

    # If running with post-discretization the last lines of optimize_transmission_expansion_iteratively have to be undone for the operational run
    if solve_opts["post_discretization"].get("enable") and not solve_opts.get(
        "skip_iterations"
    ):
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
        planning_horizons=st_years,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
        regret_run=True,
    )

    n.add(
        "Generator",
        "co2 atmosphere",
        bus="co2 atmosphere",
        p_min_pu=-1,
        p_max_pu=0,
        p_nom_extendable=True,
        carrier="co2",
        marginal_cost=co2_prices[int(st_years)],
    )
    logger.info("Adding negative CO2 generator and dropping co2 limits.")
    n.global_constraints.drop("CO2Limit", inplace=True)
    n.global_constraints.drop("co2_limit-DE", inplace=True)

    n.export_to_netcdf(snakemake.output.st_low_res_prenetwork)
