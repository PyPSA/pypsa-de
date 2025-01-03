# -*- coding: utf-8 -*-
import logging

import pyam

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("retrieve_ariadne_database")

    configure_logging(snakemake)
    logger.info(
        f"Retrieving from IIASA database {snakemake.params.db_name}\nmodels {list(snakemake.params.leitmodelle.values())}\nscenarios {snakemake.params.scenarios}"
    )

    db = pyam.read_iiasa(
        snakemake.params.db_name,
        model=snakemake.params.leitmodelle.values(),
        scenario=snakemake.params.scenarios,
        # Download only the most recent iterations of scenarios
    )

    logger.info(f"Successfully retrieved database.")
    db.timeseries().to_csv(snakemake.output.data)
