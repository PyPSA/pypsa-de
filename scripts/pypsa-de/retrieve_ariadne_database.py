import logging

import pandas as pd
import pyam

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("retrieve_ariadne_database")

    configure_logging(snakemake)
    if snakemake.params.get("source") == "primary":
        logger.info("Retrieving from IIASA database 'ariadne2'.")

        db = pyam.read_iiasa("ariadne2")

        logger.info("Successfully retrieved database.")
        db.timeseries().to_csv(snakemake.output.data)

    elif snakemake.params.get("source") == "archive":
        # Read all sheets first; then select the one called "data".
        sheets = pd.read_excel(snakemake.input.raw_xlsx, sheet_name=None)
        df = sheets["data"]
        df.loc[df["region"] == "DEU", "region"] = "Deutschland"
        df.to_csv(snakemake.output.data, index=False)
        template = sheets["variables"]
        template = template.rename(
            columns={
                "unit": "Unit",
                "variable": "Variable",
            }
        )
        template.to_csv(snakemake.output.template, index=False)
