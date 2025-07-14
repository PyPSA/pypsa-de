import logging

import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def get_transport_data(
    db,
    year,
    non_land_liquids,
    ageb_for_mobility=True,
    uba_for_mobility=False,
):
    """
    Retrieve the German mobility demand from the transport_data model.

    Sum over the subsectors Bus, LDV, Rail, and Truck for the fuels
    electricity, hydrogen, and synthetic fuels.
    """
    subsectors = ["Bus", "LDV", "Rail", "Truck"]
    fuels = ["Electricity", "Hydrogen", "Liquids"]

    transport_demand = pd.Series(0.0, index=fuels)

    if year == "2020":
        logger.info(
            "For 2020, using hard-coded transport data from the Ariadne2-internal database."
        )

        transport_demand = pd.Series(
            {
                "Electricity": 0.0 + 17.0 + 35.82 + 0.0,
                "Hydrogen": 0.0 + 0.0 + 0.0 + 0.0,
                "Liquids": 41.81 + 1369.34 + 11.18 + 637.23,
            }
        )

        transport_demand = transport_demand.div(3.6e-6)  # convert PJ to MWh
        transport_demand["number_of_cars"] = 0.658407 + 0.120261  # BEV + PHEV

        if ageb_for_mobility or uba_for_mobility:
            if uba_for_mobility:
                logger.warning(
                    "For 2020, using historical AGEB and KBA data instead of UBA projections."
                )
            # AGEB 2020, https://ag-energiebilanzen.de/daten-und-fakten/bilanzen-1990-bis-2030/?_jahresbereich-bilanz=2011-2020
            transport_demand = pd.Series(
                {
                    "Electricity": 39129 + 2394,  # Schiene + Straße
                    "Hydrogen": 0,
                    "Liquids": 140718
                    + 1261942
                    + 10782
                    + 638820,  # Bio Strasse + Diesel Strasse + Diesel Schiene + Otto Strasse
                }
            )
            transport_demand = transport_demand.div(3.6e-3)  # convert PJ to MWH
            # https://www.kba.de/DE/Statistik/Produktkatalog/produkte/Fahrzeuge/fz27_b_uebersicht.html
            # FZ27_202101, table FZ 27.2, 1. January 2021:
            transport_demand["number_of_cars"] = 0.358498 + 0.280149

    elif year == "2025" and uba_for_mobility:
        # https://www.umweltbundesamt.de/sites/default/files/medien/11850/publikationen/projektionsbericht_2025.pdf, Abbildung 64 & 59,
        transport_demand = pd.Series(
            {
                "Electricity": 21,
                "Hydrogen": 0.0,
                "Liquids": 524 + 51,
            }
        )
        transport_demand["Liquids"] -= non_land_liquids[
            int(year)
        ]  # remove domestic navigation and aviation from UBA data to avoid double counting
        transport_demand = transport_demand.mul(1e6)  # convert TWh to MWh
        transport_demand["number_of_cars"] = 2.7 + 1.2  # BEV + PHEV

    elif year == "2030" and uba_for_mobility:
        transport_demand = pd.Series(
            {
                "Electricity": 57,
                "Hydrogen": 14,
                "Liquids": 418 + 34 + 1,
            }
        )
        transport_demand["Liquids"] -= non_land_liquids[int(year)]
        transport_demand = transport_demand.mul(1e6)
        transport_demand["number_of_cars"] = 8.7 + 1.8

    elif year == "2035" and uba_for_mobility:
        transport_demand = pd.Series(
            {
                "Electricity": 117,
                "Hydrogen": 36,
                "Liquids": 237 + 26 + 1,
            }
        )
        transport_demand["Liquids"] -= non_land_liquids[int(year)]
        transport_demand = transport_demand.mul(1e6)
        transport_demand["number_of_cars"] = 18.9 + 1.8

    else:
        if uba_for_mobility:
            logger.error(
                f"Year {year} is not supported for UBA mobility projections. Please use only 2020, 2025, 2030, 2035."
            )

        df = db[year].loc[snakemake.params.leitmodelle["transport"]]

        for fuel in fuels:
            for subsector in subsectors:
                key = f"Final Energy|Transportation|{subsector}|{fuel}"
                transport_demand.loc[fuel] += df.get((key, "TWh/yr"), 0.0)

        transport_demand = transport_demand.mul(1e6)  # convert TWh to MWh
        transport_demand["number_of_cars"] = (
            df.loc["Stock|Transportation|LDV|BEV", "million"]
            + df.loc["Stock|Transportation|LDV|PHEV", "million"]
        )

    return transport_demand


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_exogenous_mobility_demand",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2030",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

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
    # get transport_data data
    transport_data = get_transport_data(
        db,
        snakemake.wildcards.planning_horizons,
        non_land_liquids,
        ageb_for_mobility=snakemake.params.ageb_for_mobility,
        uba_for_mobility=snakemake.params.uba_for_mobility,
    )

    # get German mobility weighting
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    # only get German data
    pop_layout = pop_layout[pop_layout.ct == "DE"].fraction

    mobility_demand = pd.DataFrame(
        pop_layout.values[:, None] * transport_data.values,
        index=pop_layout.index,
        columns=transport_data.index,
    )

    mobility_demand.to_csv(snakemake.output.mobility_demand)
