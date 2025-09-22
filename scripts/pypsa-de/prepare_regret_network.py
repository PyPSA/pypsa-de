# Import the function dynamically since the folder name contains a hyphen which is invalid in a module name.
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

_spec_path = pathlib.Path(__file__).resolve().parent / "modify_prenetwork.py"
_spec = importlib.util.spec_from_file_location(
    "scripts.pypsa_de.modify_prenetwork", _spec_path
)
_modify_prenetwork = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_modify_prenetwork)
remove_flexibility_options = _modify_prenetwork.remove_flexibility_options

logger = logging.getLogger(__name__)

uc_params_custom = {
    "OCGT": {
        "p_min_pu": 0.2,
        "start_up_cost": 20,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 1,
    },
    "CCGT": {
        "p_min_pu": 0.45,
        "start_up_cost": 80,
        "min_up_time": 3,
        "min_down_time": 2,
        "ramp_limit_up": 1,
    },
    "coal": {
        "p_min_pu": 0.325,
        "start_up_cost": 60,
        "min_up_time": 6,
        "min_down_time": 6,
        "ramp_limit_up": 1,
    },
    "lignite": {
        "p_min_pu": 0.325,
        "start_up_cost": 80,
        "min_up_time": 10,
        "min_down_time": 10,
        "ramp_limit_up": 1,
    },
    "nuclear": {
        "p_min_pu": 0.5,
        "start_up_cost": 100,
        "min_up_time": 8,
        "min_down_time": 10,
    },
    "oil": {
        "p_min_pu": 0.2,
        "start_up_cost": 30,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 1,
    },
    "urban central solid biomass CHP": {
        "p_min_pu": 0.38,
        "start_up_cost": 50,
        "min_up_time": 2,
        "min_down_time": 2,
    },
}

uc_params_optimistic = {
    "OCGT": {
        "p_min_pu": 0.15,
        "start_up_cost": 20,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 1,
        "ramp_limit_start_up": 0.3,
        "ramp_limit_shut_down": 0.3,
    },
    "CCGT": {
        "p_min_pu": 0.4,
        "start_up_cost": 100,
        "min_up_time": 2,
        "min_down_time": 2,
        "ramp_limit_up": 1,
        "ramp_limit_start_up": 0.5,
        "ramp_limit_shut_down": 0.5,
    },
    "coal": {
        "p_min_pu": 0.3,
        "start_up_cost": 80,
        "min_up_time": 4,
        "min_down_time": 4,
        "ramp_limit_up": 0.8,
        "ramp_limit_start_up": 0.3,
        "ramp_limit_shut_down": 0.3,
    },
    "lignite": {
        "p_min_pu": 0.3,
        "start_up_cost": 120,
        "min_up_time": 5,
        "min_down_time": 5,
        "ramp_limit_up": 0.7,
        "ramp_limit_start_up": 0.3,
        "ramp_limit_shut_down": 0.3,
    },
    "nuclear": {
        "p_min_pu": 0.45,
        "start_up_cost": 200,
        "min_up_time": 6,
        "min_down_time": 8,
        "ramp_limit_up": 0.2,
        "ramp_limit_start_up": 0.2,
        "ramp_limit_shut_down": 0.2,
    },
    "oil": {
        "p_min_pu": 0.2,
        "start_up_cost": 30,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 1,
        "ramp_limit_start_up": 0.3,
        "ramp_limit_shut_down": 0.3,
    },
    "urban central solid biomass CHP": {
        "p_min_pu": 0.35,
        "start_up_cost": 50,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 0.9,
        "ramp_limit_start_up": 0.4,
        "ramp_limit_shut_down": 0.4,
    },
}

uc_params_average = {
    "OCGT": {
        "p_min_pu": 0.2,
        "start_up_cost": 40,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 1,
        "ramp_limit_start_up": 0.2,
        "ramp_limit_shut_down": 0.2,
    },
    "CCGT": {
        "p_min_pu": 0.45,
        "start_up_cost": 150,
        "min_up_time": 3,
        "min_down_time": 2,
        "ramp_limit_up": 1,
        "ramp_limit_start_up": 0.45,
        "ramp_limit_shut_down": 0.45,
    },
    "coal": {
        "p_min_pu": 0.325,
        "start_up_cost": 120,
        "min_up_time": 5,
        "min_down_time": 6,
        "ramp_limit_up": 0.7,
        "ramp_limit_start_up": 0.38,
        "ramp_limit_shut_down": 0.38,
    },
    "lignite": {
        "p_min_pu": 0.325,
        "start_up_cost": 150,
        "min_up_time": 7,
        "min_down_time": 6,
        "ramp_limit_up": 0.6,
        "ramp_limit_start_up": 0.4,
        "ramp_limit_shut_down": 0.4,
    },
    "nuclear": {
        "p_min_pu": 0.5,
        "start_up_cost": 250,
        "min_up_time": 6,
        "min_down_time": 10,
        "ramp_limit_up": 0.3,
        "ramp_limit_start_up": 0.3,
        "ramp_limit_shut_down": 0.3,
    },
    "oil": {
        "p_min_pu": 0.2,
        "start_up_cost": 50,
        "min_up_time": 1,
        "min_down_time": 1,
        "ramp_limit_up": 0.8,
        "ramp_limit_start_up": 0.2,
        "ramp_limit_shut_down": 0.2,
    },
    "urban central solid biomass CHP": {
        "p_min_pu": 0.38,
        "start_up_cost": 80,
        "min_up_time": 2,
        "min_down_time": 2,
        "ramp_limit_up": 0.7,
        "ramp_limit_start_up": 0.38,
        "ramp_limit_shut_down": 0.38,
    },
}

uc_params_conservative = {
    "OCGT": {
        "p_min_pu": 0.25,
        "start_up_cost": 60,
        "min_up_time": 2,
        "min_down_time": 2,
        "ramp_limit_up": 0.9,
        "ramp_limit_start_up": 0.15,
        "ramp_limit_shut_down": 0.15,
    },
    "CCGT": {
        "p_min_pu": 0.5,
        "start_up_cost": 200,
        "min_up_time": 4,
        "min_down_time": 3,
        "ramp_limit_up": 0.8,
        "ramp_limit_start_up": 0.35,
        "ramp_limit_shut_down": 0.35,
    },
    "coal": {
        "p_min_pu": 0.35,
        "start_up_cost": 160,
        "min_up_time": 6,
        "min_down_time": 8,
        "ramp_limit_up": 0.5,
        "ramp_limit_start_up": 0.25,
        "ramp_limit_shut_down": 0.25,
    },
    "lignite": {
        "p_min_pu": 0.35,
        "start_up_cost": 200,
        "min_up_time": 8,
        "min_down_time": 10,
        "ramp_limit_up": 0.4,
        "ramp_limit_start_up": 0.2,
        "ramp_limit_shut_down": 0.2,
    },
    "nuclear": {
        "p_min_pu": 0.55,
        "start_up_cost": 400,
        "min_up_time": 10,
        "min_down_time": 12,
        "ramp_limit_up": 0.15,
        "ramp_limit_start_up": 0.15,
        "ramp_limit_shut_down": 0.15,
    },
    "oil": {
        "p_min_pu": 0.25,
        "start_up_cost": 80,
        "min_up_time": 2,
        "min_down_time": 2,
        "ramp_limit_up": 0.6,
        "ramp_limit_start_up": 0.15,
        "ramp_limit_shut_down": 0.15,
    },
    "urban central solid biomass CHP": {
        "p_min_pu": 0.4,
        "start_up_cost": 120,
        "min_up_time": 3,
        "min_down_time": 3,
        "ramp_limit_up": 0.5,
        "ramp_limit_start_up": 0.3,
        "ramp_limit_shut_down": 0.3,
    },
}


def add_unit_commitment(
    n,
    uc_params=uc_params_average,
    carriers=["OCGT", "coal", "lignite", "urban central solid biomass CHP"],
    regions=["DE"],
):
    """
    Add unit commitment parameters to links in the network based on a UC parameter dictionary.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network.

    uc_params : dict
        Nested dict with carrier names as keys and dict of UC parameters as values.
        Example:
        {
            "OCGT": {
                "p_min_pu": 0.2,
                "start_up_cost": 40,
                "min_up_time": 1,
                "min_down_time": 1,
                "ramp_limit_up": 1,
                "ramp_limit_start_up": 0.2,
                "ramp_limit_shut_down": 0.2,
            },
            ...
        }

    carriers : list, optional
        List of carriers to process (default = all carriers in uc_params).

    regions : list
        List of region codes to filter buses (default = ["DE"]).
    """

    def get_filtered_links(carrier_list):
        carrier_mask = n.links.carrier.isin(carrier_list)
        region_mask = n.links.bus0.str.contains(
            "|".join(regions), na=False
        ) | n.links.bus1.str.contains("|".join(regions), na=False)
        return n.links[carrier_mask & region_mask].index

    # If no carriers specified, use all available in uc_params
    carriers_to_process = carriers if carriers is not None else list(uc_params.keys())

    available_carriers = set(carriers_to_process) & set(n.links.carrier.unique())

    for carrier in available_carriers:
        links_i = get_filtered_links([carrier])
        if len(links_i) == 0:
            continue

        # apply UC parameters from dict
        for param, value in uc_params[carrier].items():
            if param in n.links.columns:
                n.links.loc[links_i, param] = value

        # ensure committable flag
        n.links.loc[links_i, "committable"] = True


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


def fix_capacities(realization, decision, scope="DE", strict=False, no_flex=False):
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
    # Copy the whole realization network, because copying loads_t directly is not type stable
    n.loads_t = realization.copy().loads_t

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

        _idx = new.query(f"{attr}_extendable").index

        new.loc[_idx, attr] = deci.loc[_idx, attr]

    if no_flex:
        logger.info("Realization network is from a run without flexibility.")
        remove_flexibility_options(n)
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "prepare_regret_network",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
            decision="LowDemand",
            run="HighDemand",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Touch output file to ensure it exists
    pathlib.Path(snakemake.output.regret_prenetwork).touch()

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

    strict = snakemake.params["strict"]
    scope_to_fix = snakemake.params["scope_to_fix"]
    h2_vent = snakemake.params["h2_vent"]

    n = fix_capacities(
        realization,
        decision,
        scope=scope_to_fix,
        strict=strict,
        no_flex=decision.meta.get("iiasa_database").get("no_flex_lt_run", False),
    )

    unit_commitment = snakemake.params.get("unit_commitment")

    if unit_commitment["enable"]:
        logger.info(
            f"Add unit commitment in {unit_commitment['regions']} for carriers {unit_commitment['carriers']} with parameter set '{unit_commitment['params']}' to the network."
        )

        uc_params_str = unit_commitment["params"]
        if uc_params_str == "custom":
            uc_params = uc_params_custom
        elif uc_params_str == "optimistic":
            uc_params = uc_params_optimistic
        elif uc_params_str == "conservative":
            uc_params = uc_params_conservative
        elif uc_params_str == "average":
            uc_params = uc_params_average

        add_unit_commitment(
            n,
            uc_params,
            carriers=unit_commitment["carriers"],
            regions=unit_commitment["regions"],
        )

    if strict:
        logger.info(
            "Strict regret run chosen. No capacities are extendable. Activating load shedding to prevent infeasibilites."
        )
        snakemake.params.solving["options"]["load_shedding"] = True

    if h2_vent:
        logger.info("H2 venting activated for regret run.")
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
        # n.generators_t.p[n.generators.query("carrier == 'H2 vent'").index].T.mul(n.snapshot_weightings.generators).T.sum()

    # Manipulating the global constraints
    to_keep = [
        "biomass limit",
        "unsustainable biomass limit",
        "co2_sequestration_limit",
        "CO2Limit",
        "co2_limit-DE",
    ]

    to_drop = n.global_constraints.index.difference(to_keep)

    logger.info("Regret run detected. Dropping the following constraints:")
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
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
        regret_run=True,
    )

    # These constraints have to be changed AFTER prepare_network

    if scope_to_fix == "EU":
        logger.info(
            f"Fixing Scope 'EU' chosen. Setting the EU CO2 price to the sum of the EU and DE CO2 prices from the realization network: {realization.global_constraints.loc['CO2Limit', 'mu'] + realization.global_constraints.loc['co2_limit-DE', 'mu']} €/t_CO2"
        )
        n.add(
            "Generator",
            "co2 atmosphere",
            bus="co2 atmosphere",
            p_min_pu=-1,
            p_max_pu=0,
            p_nom_extendable=True,
            carrier="co2",
            marginal_cost=(
                realization.global_constraints.loc["CO2Limit", "mu"]
                + realization.global_constraints.loc["co2_limit-DE", "mu"]
            ),
        )
        logger.info("Adding negative CO2 generator and dropping co2 limits.")
        n.global_constraints.drop("CO2Limit", inplace=True)
        n.global_constraints.drop("co2_limit-DE", inplace=True)

    n.export_to_netcdf(snakemake.output.regret_prenetwork)
