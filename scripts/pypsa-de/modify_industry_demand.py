# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script modifies the industrial production values to match the FORECAST
model.

This includes
- Production|Non-Metallic Minerals|Cement
- Production|Steel
- Production|Chemicals|Ammonia
- Production|Chemicals|Methanol
- Production|Pulp and Paper
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)
if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "modify_industry_demand",
            simpl="",
            clusters=22,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Mix",
            planning_horizons=2020,
        )

    configure_logging(snakemake)

    year = snakemake.input.industrial_production_per_country_tomorrow.split("_")[
        -1
    ].split(".")[0]
    existing_industry = pd.read_csv(
        snakemake.input.industrial_production_per_country_tomorrow, index_col=0
    )

    # read in ariadne database
    if year == "2020":
        logger.info(
            "For 2020, using hardcoded values from FORECAST from the ariadne2-internal database",
        )
        ariadne = pd.DataFrame(
            {
                "Production|Chemicals|Ammonia": 2.891851,
                "Production|Chemicals|Methanol": 1.359,
                "Production|Non-Metallic Minerals": 68.635925,
                "Production|Non-Metallic Minerals|Cement": 34.966,
                "Production|Pulp and Paper": 40.746,
                "Production|Steel": 40.621,
            },
            index=[year],
        ).T.multiply(1000)
    else:
        # leitmodell for industry demand
        leitmodell = "FORECAST v1.0"
        ariadne = (
            pd.read_csv(
                snakemake.input.ariadne,
                index_col=["model", "scenario", "region", "variable", "unit"],
            )
            .loc[
                leitmodell,
                snakemake.params.reference_scenario,
                "Deutschland",
                :,
                "Mt/yr",
            ]
            .multiply(1000)
        )

    logger.info(
        "German industry demand before modification",
    )
    logger.info(
        existing_industry.loc[
            "DE",
            [
                "Cement",
                "Electric arc",
                "Integrated steelworks",
                "DRI + Electric arc",
                "Ammonia",
                "Methanol",
                "Pulp production",
                "Paper production",
                "Ceramics & other NMM",
            ],
        ],
    )

    # write Cement, Ammonia and Methanol directly to dataframe
    existing_industry.loc["DE", "Cement"] = ariadne.loc[
        "Production|Non-Metallic Minerals|Cement", year
    ]
    existing_industry.loc["DE", "Ammonia"] = ariadne.loc[
        "Production|Chemicals|Ammonia", year
    ]
    existing_industry.loc["DE", "Methanol"] = ariadne.loc[
        "Production|Chemicals|Methanol", year
    ]

    # get ratio of pulp and paper production
    pulp_ratio = existing_industry.loc["DE", "Pulp production"] / (
        existing_industry.loc["DE", "Pulp production"]
        + existing_industry.loc["DE", "Paper production"]
    )

    existing_industry.loc["DE", "Pulp production"] = (
        ariadne.loc["Production|Pulp and Paper", year] * pulp_ratio
    )
    existing_industry.loc["DE", "Paper production"] = ariadne.loc[
        "Production|Pulp and Paper", year
    ] * (1 - pulp_ratio)

    # non-metallic minerals
    existing_industry.loc["DE", "Ceramics & other NMM"] = (
        ariadne.loc["Production|Non-Metallic Minerals", year]
        - ariadne.loc["Production|Non-Metallic Minerals|Cement", year]
    )

    # get steel ratios from existing_industry
    steel = existing_industry.loc[
        "DE", ["Electric arc", "Integrated steelworks", "DRI + Electric arc"]
    ]
    ratio = steel / steel.sum()

    # multiply with steel production including primary and secondary steel since distinguishing is taken care of later
    existing_industry.loc[
        "DE", ["Electric arc", "Integrated steelworks", "DRI + Electric arc"]
    ] = ratio * ariadne.loc["Production|Steel", year]

    logger.info("German demand after modification")
    logger.info(
        existing_industry.loc[
            "DE",
            [
                "Cement",
                "Electric arc",
                "Integrated steelworks",
                "DRI + Electric arc",
                "Ammonia",
                "Methanol",
                "Pulp production",
                "Paper production",
                "Ceramics & other NMM",
            ],
        ],
    )

    existing_industry.to_csv(
        snakemake.output.industrial_production_per_country_tomorrow
    )
