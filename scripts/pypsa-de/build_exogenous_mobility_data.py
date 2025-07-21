import logging

import pandas as pd

from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def get_mobility_data(
    db,
    year,
    non_land_liquids,
    ageb_for_mobility=True,
    uba_for_mobility=False,
):
    """
    Retrieve the German mobility demand from the transport model.

    Sum over the subsectors Bus, LDV, Rail, and Truck for the fuels
    electricity, hydrogen, and synthetic fuels.
    """
    subsectors = ["Bus", "LDV", "Rail", "Truck"]
    fuels = ["Electricity", "Hydrogen", "Liquids"]

    mobility_data = pd.Series(0.0, index=fuels)

    if year == "2020":
        logger.info(
            "For 2020, using hard-coded transport data from the Ariadne2-internal database."
        )

        mobility_data = pd.Series(
            {
                "Electricity": 0.0 + 17.0 + 35.82 + 0.0,
                "Hydrogen": 0.0 + 0.0 + 0.0 + 0.0,
                "Liquids": 41.81 + 1369.34 + 11.18 + 637.23,
            }
        )

        mobility_data = mobility_data.div(3.6e-6)  # convert PJ to MWh
        mobility_data["million_EVs"] = 0.658407 + 0.120261  # BEV + PHEV

        if ageb_for_mobility or uba_for_mobility:
            if uba_for_mobility:
                logger.warning(
                    "For 2020, using historical AGEB and KBA data instead of UBA projections."
                )
            # AGEB 2020, https://ag-energiebilanzen.de/daten-und-fakten/bilanzen-1990-bis-2030/?_jahresbereich-bilanz=2011-2020
            mobility_data = pd.Series(
                {
                    "Electricity": 39129.0 + 2394.0,  # Schiene + Straße
                    "Hydrogen": 0.0,
                    "Liquids": 140718.0
                    + 1261942.0
                    + 10782.0
                    + 638820.0,  # Bio Strasse + Diesel Strasse + Diesel Schiene + Otto Strasse
                }
            )
            mobility_data = mobility_data.div(3.6e-3)  # convert PJ to MWH
            # https://www.kba.de/DE/Statistik/Produktkatalog/produkte/Fahrzeuge/fz27_b_uebersicht.html
            # FZ27_202101, table FZ 27.2, 1. January 2021:
            mobility_data["million_EVs"] = 0.358498 + 0.280149

    elif year == "2025" and uba_for_mobility:
        # https://www.umweltbundesamt.de/sites/default/files/medien/11850/publikationen/projektionsbericht_2025.pdf, Abbildung 64 & 59,
        mobility_data = pd.Series(
            {
                "Electricity": 21.0,
                "Hydrogen": 0.0,
                "Liquids": 524.0 + 51.0,
            }
        )
        mobility_data["Liquids"] -= non_land_liquids[
            int(year)
        ]  # remove domestic navigation and aviation from UBA data to avoid double counting
        mobility_data = mobility_data.mul(1e6)  # convert TWh to MWh
        mobility_data["million_EVs"] = 2.7 + 1.2  # BEV + PHEV

    elif year == "2030" and uba_for_mobility:
        mobility_data = pd.Series(
            {
                "Electricity": 57.0,
                "Hydrogen": 14.0,
                "Liquids": 418.0 + 34.0 + 1.0,
            }
        )
        mobility_data["Liquids"] -= non_land_liquids[int(year)]
        mobility_data = mobility_data.mul(1e6)
        mobility_data["million_EVs"] = 8.7 + 1.8

    elif year == "2035" and uba_for_mobility:
        mobility_data = pd.Series(
            {
                "Electricity": 117.0,
                "Hydrogen": 36.0,
                "Liquids": 237.0 + 26.0 + 1.0,
            }
        )
        mobility_data["Liquids"] -= non_land_liquids[int(year)]
        mobility_data = mobility_data.mul(1e6)
        mobility_data["million_EVs"] = 18.9 + 1.8

    else:
        if uba_for_mobility:
            logger.error(
                f"Year {year} is not supported for UBA mobility projections. Please use only 2020, 2025, 2030, 2035."
            )
            raise NotImplementedError(
                f"Year {year} is not supported for UBA mobility projections. Please use only 2020, 2025, 2030, 2035."
            )

        df = db[year].loc[snakemake.params.leitmodelle["transport"]]

        for fuel in fuels:
            for subsector in subsectors:
                key = f"Final Energy|Transportation|{subsector}|{fuel}"
                mobility_data.loc[fuel] += df.get((key, "TWh/yr"), 0.0)

        mobility_data = mobility_data.mul(1e6)  # convert TWh to MWh
        mobility_data["million_EVs"] = (
            df.loc["Stock|Transportation|LDV|BEV", "million"]
            + df.loc["Stock|Transportation|LDV|PHEV", "million"]
        )

    return mobility_data


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_exogenous_mobility_data",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2020",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    db = pd.read_csv(
        snakemake.input.ariadne,
        index_col=["model", "scenario", "region", "variable", "unit"],
    ).loc[
        :,
        snakemake.params.reference_scenario,
        "Deutschland",
        :,
        :,
    ]

    energy_totals = (
        pd.read_csv(
            snakemake.input.energy_totals,
            index_col=[0, 1],
        )
        .xs(
            snakemake.params.energy_totals_year,
            level="year",
        )
        .loc["DE"]
    )

    domestic_aviation = energy_totals.loc["total domestic aviation"] * pd.Series(
        snakemake.params.aviation_demand_factor
    )

    domestic_navigation = energy_totals.loc["total domestic navigation"] * pd.Series(
        snakemake.params.shipping_oil_share
    )

    non_land_liquids = domestic_aviation + domestic_navigation

    logger.info(
        f"Retrieving German mobility demand from {snakemake.params.leitmodelle['transport']} transport model."
    )
    # get mobility_data data
    mobility_data = get_mobility_data(
        db,
        snakemake.wildcards.planning_horizons,
        non_land_liquids,
        ageb_for_mobility=snakemake.params.ageb_for_mobility,
        uba_for_mobility=snakemake.params.uba_for_mobility,
    )

    mobility_data.to_csv(snakemake.output.mobility_data, header=False)
