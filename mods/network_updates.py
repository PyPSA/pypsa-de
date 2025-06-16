"""Functions to update the industry demand in the PyPSA-AT model."""

from logging import getLogger

logger = getLogger(__name__)


def update_austrian_grid_capacities(n, snakemake):
    """
    Update transmission capacities for Austria.

    The function is expected to run on clustered pre-networks. It
    Will read capacities provided in a data file and update the
    respective values.

    Parameters
    ----------
    n
        The pre-network to update during rule `modify_prenetwork`.
    snakemake
        The Snakemake object.

    Returns
    -------
    :
    """
    logger.info("Updating grid capacities for Austria.")
    # for every year:
    # transmission carrier:
    # 'AC',
    # 'gas pipeline new'
    # 'gas pipeline',
    # 'DC',
    # 'H2 pipeline',
    # 'H2 pipeline (Kernnetz)',
    # 'H2 pipeline retrofitted'


def update_austrian_industry_demand(existing_industry, year):
    """Update the industry demand in the PyPSA-AT model for Austria."""

    logger.info("Updating industry demand for Austria.")

    return existing_industry
