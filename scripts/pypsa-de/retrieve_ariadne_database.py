import logging

import pyam

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("retrieve_ariadne_database")

    configure_logging(snakemake)
    logger.info("Retrieving from IIASA database 'ariadne2'.")

    db = pyam.read_iiasa("ariadne2")

    logger.info("Successfully retrieved database.")
    db.timeseries().to_csv(snakemake.output.data)
