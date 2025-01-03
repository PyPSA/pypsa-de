# -*- coding: utf-8 -*-
import logging

import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "modify_existing_heating",
            run="KN2045_Bal_v4",
        )

    configure_logging(snakemake)

    leitmodell = snakemake.params.leitmodelle["buildings"]
    logger.info(f"Using {leitmodell} for heating demand modification.")

    existing_heating = pd.read_csv(snakemake.input.existing_heating, index_col=0)

    ariadne = pd.read_csv(
        snakemake.input.ariadne,
        index_col=["model", "scenario", "region", "variable", "unit"],
    ).loc[
        leitmodell,
        snakemake.params.fallback_reference_scenario,
        "Deutschland",
        :,
        "million",
    ]

    logger.info("Heating demand before modification:{existing_heating.loc['Germany']}")

    mapping = {
        "gas boiler": "Gas Boiler",
        "oil boiler": "Oil Boiler",
        "air heat pump": "Heat Pump|Electrical|Air",
        "ground heat pump": "Heat Pump|Electrical|Ground",
        "biomass boiler": "Biomass Boiler",
    }

    year = "2020"
    for tech in mapping:
        stock = ariadne.at[
            f"Stock|Space Heating|{mapping[tech]}",
            year,
        ]

        peak = (
            stock
            * existing_heating.loc["Germany"].sum()
            / ariadne.at[f"Stock|Space Heating", year]
        )
        existing_heating.at["Germany", tech] = peak

    logger.info(
        f"Heating demand after modification with {leitmodell}: {existing_heating.loc['Germany']}"
    )

    existing_heating.to_csv(snakemake.output.existing_heating)
