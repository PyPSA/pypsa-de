"""Functions to update the industry demand in the PyPSA-AT model."""

from logging import getLogger

import pandas as pd
import pypsa

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
