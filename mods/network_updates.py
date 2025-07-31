"""Functions to update the industry demand in the PyPSA-AT model."""

# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
from logging import getLogger

import pandas as pd
import pypsa
from snakemake.script import Snakemake

logger = getLogger(__name__)


def modify_austrian_transmission_capacities(
    n: pypsa.Network, austrian_transmission_capacities: str
):
    """
    Update transmission capacities for Austria.

    The function is expected to run on clustered pre-networks. It
    Will read capacities provided in a data file and update the
    respective values.

    Parameters
    ----------
    n
        The pre-network to update during rule `modify_prenetwork`.

    austrian_transmission_capacities
        The path to the data file used to update the capacities.

    Returns
    -------
    :
    """
    logger.info("Modifying grid capacities for Austria.")

    # transmission_carrier = get_transmission_carriers(n)
    # to_concat = []
    # for component, carrier in transmission_carrier:
    #     capacity_column = f"{'p' if component == 'Link' else 's'}_nom"
    #     to_concat.append(
    #         n.static(component).query(f"carrier == @carrier "
    #                                   f"& (bus0.str.startswith('AT') "
    #                                   f"| bus1.str.startswith('AT'))")[["bus0", "bus1", capacity_column]]
    #     )
    # template = pd.concat(to_concat).sort_index()
    # template.to_csv(austrian_grid_capacities)

    capacities = pd.read_csv(austrian_transmission_capacities, index_col=0).sort_index()

    for c in n.branch_components:
        p = f"{'p' if c == 'Link' else 's'}_nom"
        overwrite = capacities[["bus0", "bus1", p]].dropna(subset=[p])
        n.static(c).update(overwrite)

    # todo: test if 2020 capacities are in result network
    # todo: support all years. currently only 2020 is possible


def modify_austrian_industry_demand(existing_industry, year):
    """Update the industry demand in the PyPSA-AT model for Austria."""

    logger.info("Updating industry demand for Austria.")

    return existing_industry


def modify_austrian_gas_storage_capacities():
    """Update gas and H2 storage capacities for Austria."""


def modify_biomass_potentials():
    """Update biomass potentials."""


def modify_heat_demand():
    """Update heat demands."""


def electricity_base_load_split(n: pypsa.Network, snakemake: Snakemake):
    """Split electricity base load to sectoral loads."""


def unravel_gas_import_and_production(
    n: pypsa.Network, snakemake: Snakemake, costs: pd.DataFrame
):
    """
    Differentiate LNG, pipeline and production gas generators.

    Production is cheaper than pipeline gas and LNG is
    more expensive than pipeline gas.

    Parameters
    ----------
    n
        The network before optimisation.
    snakemake
        The snakemake workflow object.
    costs
        The costs data for the current planning horizon.

    Returns
    -------
    :
    """
    config = snakemake.config
    gas_generators = n.static("Generator").query("carrier == 'gas'")
    if gas_generators.empty and config.get("gas_compression_losses", 0):
        logger.debug(
            "Skipping unravel gas generators because "
            "industry.gas_compression_losses is set."
        )
        return

    if not config.get("mods", {}).get("unravel_natural_gas_imports", {}).get("enable"):
        logger.debug(
            "Skipping unravel natural gas imports because "
            "the modification was not requested."
        )
        return

    logger.info("Unravel gas import types.")
    gas_input_nodes = pd.read_csv(
        snakemake.input.gas_input_nodes_simplified, index_col=0
    )

    # remove combined gas generators
    n.remove("Generator", gas_generators.index)
    ariadne_gas_fuel_price = costs.at["gas", "fuel"]
    cost_factors = config["mods"]["unravel_natural_gas_imports"]

    for import_type in ("lng", "pipeline", "production"):
        cost_factor = cost_factors[import_type]
        p_nom = gas_input_nodes[import_type].dropna()
        p_nom.rename(lambda x: x + " gas", inplace=True)
        nodes = p_nom.index
        suffix = (
            " production" if import_type == "production" else f" {import_type} import"
        )
        carrier = f"{import_type} gas"
        marginal_cost = ariadne_gas_fuel_price * cost_factor
        n.add(
            "Generator",
            nodes,
            suffix=suffix,
            bus=nodes,
            carrier=carrier,
            p_nom_extendable=False,
            marginal_cost=marginal_cost,
            p_nom=p_nom,
        )

    # make sure that the total gas generator capacity was not changed by this modification
    old_p_nom = gas_generators["p_nom"].sum()
    new_p_nom = (
        n.static("Generator").query("carrier.str.endswith(' gas')")["p_nom"].sum()
    )
    assert old_p_nom.round(8) == new_p_nom.round(8), (
        f"Unraveling imports changed total capacities: old={old_p_nom}, new={new_p_nom}."
    )


def unravel_electricity_base_load(n: pypsa.Network, snakemake: Snakemake) -> None:
    """
    Split electricity baseload into sectoral loads.

    Parameters
    ----------
    n
    snakemake

    Returns
    -------
    :
    """
    # config = snakemake.config
    # print(config)

    # electricity base load is from: https://nbviewer.org/github/Open-Power-System-Data/datapackage_timeseries/blob/2020-10-06/main.ipynb
    # total load=total generation−auxilary/self−consumption in power plants+imports−exports−consumption by storages
    # base_load = n.static("Load").query("carrier == 'electricity'")
    # print(base_load)

    # energy_totals.csv:
    # contains load data for sectors:
    #  - residential
    #  - services
    #  - transport (road, international & national navigation & aviation)
    # by energy carrier: electricity, heat, fuel
    #

    # todo: households and services
    # todo: electricity transport rail
    # todo: electricity industry
