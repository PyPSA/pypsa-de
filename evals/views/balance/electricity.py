"""Module for electricity evaluations."""

from pathlib import Path

from evals.constants import DataModel, Group, TradeTypes
from evals.fileio import Exporter
from evals.plots import ESMBarChart, ESMGroupedBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import filter_by, rename_aggregate


def view_balance_electricity(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the electricity production & demand by country and year.

    Electricity production is sum of power supply from
      - "Generator" components,
      - "Link" supply with AC bus_carrier, and
      - Pump-Hydro-Storage and Run-Of-River inflows
    """
    bus_carrier = ["AC", "low voltage", "EV battery"]
    transmission_carrier = ["AC", "DC"]
    storage_carrier = ["PHS", "BEV charger", "EV battery", "V2G"]

    # todo: read csvs and compare with calculated results to test the code
    # supply = read_pypsa_csv(result_path, "nodal_supply", index_cols=4) #.drop(transmission_carrier, level=DataModel.CARRIER)
    # demand = read_pypsa_csv(result_path, "nodal_withdrawal", index_cols=4)# .drop(transmission_carrier, level=DataModel.CARRIER)
    # print(filter_by(supply, year="2050", bus_carrier="AC").sum() / 1e6)
    # print(filter_by(demand, year="2050", bus_carrier="AC").sum() / 1e6)
    #
    # balance = read_pypsa_csv(result_path, "nodal_energy_balance", index_cols=4)
    # print(filter_by(balance, year="2050", bus_carrier="AC").drop(transmission_carrier, level=DataModel.CARRIER).sum() / 1e6)

    supply = (
        collect_myopic_statistics(
            networks,
            statistic="supply",
            bus_carrier=bus_carrier,
        )
        .drop(transmission_carrier, level=DataModel.CARRIER)
        # .pipe(rename_aggregate, dict.fromkeys(transmission_carrier, "Import"))
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_out))
    )
    print(filter_by(supply, year="2050", bus_carrier="AC").sum() / 1e6)

    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=bus_carrier,
        )
        .drop(transmission_carrier, level=DataModel.CARRIER)
        # .pipe(rename_aggregate, dict.fromkeys(transmission_carrier, "Export"))
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_in))
        .mul(-1)
    )

    trade_statistics = []
    for scope, direction, alias in [
        (TradeTypes.FOREIGN, "import", Group.import_foreign),
        (TradeTypes.FOREIGN, "export", Group.export_foreign),
        (TradeTypes.DOMESTIC, "import", Group.import_domestic),
        (TradeTypes.DOMESTIC, "export", Group.export_domestic),
    ]:
        trade = (
            collect_myopic_statistics(
                networks,
                statistic="trade_energy",
                scope=scope,
                direction=direction,
                bus_carrier=bus_carrier,
            )
            # the trade statistic wrongly finds transmission between EU -> country buses.
            # Those are dropped by the filter_by statement.
            .pipe(filter_by, carrier=transmission_carrier)
            .pipe(rename_aggregate, alias)
        )
        trade.attrs["unit"] = "MWh_el"
        trade_statistics.append(trade)

    exporter = Exporter(
        statistics=[supply, demand] + trade_statistics,
        view_config=config["view"],
    )

    exporter.defaults.plotly.chart = ESMGroupedBarChart
    exporter.defaults.plotly.xaxis_title = ""
    # exporter.defaults.plotly.chart = ESMBarChart
    exporter.defaults.plotly.pattern = dict.fromkeys(
        [
            Group.import_foreign,
            Group.export_foreign,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )

    if exporter.defaults.plotly.chart == ESMBarChart:
        # combine bus carrier to export netted technologies, although
        # they have difference bus_carrier in index , e.g.
        # electricity distribution grid, (AC, low voltage)
        exporter.statistics[0] = rename_aggregate(demand, "AC", level="bus_carrier")
        exporter.statistics[1] = rename_aggregate(supply, "AC", level="bus_carrier")

    # todo: electricity load split

    exporter.export(result_path, subdir)
