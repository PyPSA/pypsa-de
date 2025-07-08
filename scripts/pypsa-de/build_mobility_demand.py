import logging

import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def get_transport_data(db, year, ageb_for_transport=False):
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

        transport_demand = pd.Series()
        # if 2020
        transport_demand["Electricity"] = 0.0 + 17.0 + 35.82 + 0.0
        transport_demand["Hydrogen"] = 0.0 + 0.0 + 0.0 + 0.0
        transport_demand["Liquids"] = 41.81 + 1369.34 + 11.18 + 637.23
        transport_demand = transport_demand.div(3.6e-6)  # convert PJ to MWh
        transport_demand["number_of_cars"] = 0.658407 + 0.120261  # BEV + PHEV

        if ageb_for_transport:
            # AGEB 2020, https://ag-energiebilanzen.de/daten-und-fakten/bilanzen-1990-bis-2030/?_jahresbereich-bilanz=2011-2020
            transport_demand["Electricity"] = 39129 + 2394  # Schiene + Straße
            transport_demand["Hydrogen"] = 0
            transport_demand["Liquids"] = (
                140718 + 1261942 + 10782 + 638820
            )  # Bio Strasse + Diesel Strasse + Diesel Schiene + Otto Strasse
            transport_demand = transport_demand.div(3.6e-3)  # convert TJ to MWH
            # https://www.kba.de/DE/Statistik/Produktkatalog/produkte/Fahrzeuge/fz27_b_uebersicht.html
            # FZ27_202101, table FZ 27.2, 1. January 2021:
            transport_demand["number_of_cars"] = 0.358498 + 0.280149

    elif year == "2025" and ageb_for_transport:
        # AGEB2024 for train demand 25, linear extrapolation with AGEB2024 + AGEB2023 for EVs
        transport_demand["Electricity"] = 39761 + 2 * 21270 - 16180
        transport_demand["Hydrogen"] = 0
        # AGEB2024 for Liquids demand 25
        transport_demand["Liquids"] = 116323 + 9650 + 1158250 + 702618
        transport_demand = transport_demand.div(3.6e-3)
        # FZ27_202504, 202404, table FZ 27.8,
        # linear extrapolation to 1. January 2026: "1. January 2025" + ("1. January 2025" - "1. January 2024")
        # 2 * (1,810,815 + 968,734) - (1,555,265 + 922,876) = 3080957
        # rounded upwards
        transport_demand["number_of_cars"] = 3.1  # million, BEV + PHEV

    else:
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
            "build_mobility_demand",
            simpl="",
            clusters=22,
            opts="",
            ll="vopt",
            sector_opts="none",
            planning_horizons="2020",
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

    logger.info(
        f"Retrieving German mobility demand from {snakemake.params.leitmodelle['transport']} transport model."
    )
    # get transport_data data
    transport_data = get_transport_data(
        db, snakemake.wildcards.planning_horizons, snakemake.params.ageb_for_transport
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
