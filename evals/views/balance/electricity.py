"""Module for electricity evaluations."""

from pathlib import Path

from evals.constants import DataModel, Group, TradeTypes
from evals.fileio import Exporter
from evals.plots.facetbars import ESMGroupedBarChart
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
    # (
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="trade_energy",
    #         scope=TradeTypes.DOMESTIC,
    #         direction="import",
    #         bus_carrier="AC",
    #     )
    #     # the trade statistic wrongly finds transmission between EU -> country buses.
    #     # Those are dropped by the filter_by statement.
    #     .pipe(filter_by, carrier=["AC", "DC"])
    #     .pipe(rename_aggregate, "domestic import")
    # )

    transmission_carrier = ["AC", "DC"]

    supply = collect_myopic_statistics(
        networks,
        statistic="supply",
        bus_carrier=["AC", "low voltage"],
    ).drop(transmission_carrier, level=DataModel.CARRIER)

    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=["AC", "low voltage"],
        )
        .drop(transmission_carrier, level=DataModel.CARRIER)
        .mul(-1)
    )

    trade_statistics = []
    for scope, direction, alias in [
        (TradeTypes.FOREIGN, "import", Group.import_foreign),
        (TradeTypes.FOREIGN, "export", Group.export_foreign),
        (TradeTypes.DOMESTIC, "import", Group.import_domestic),
        (TradeTypes.DOMESTIC, "export", Group.export_domestic),
    ]:
        trade_statistics.append(
            collect_myopic_statistics(
                networks,
                statistic="trade_energy",
                scope=scope,
                direction=direction,
                bus_carrier="AC",
            )
            # the trade statistic wrongly finds transmission between EU -> country buses.
            # Those are dropped by the filter_by statement.
            .pipe(filter_by, carrier=["AC", "DC"]).pipe(rename_aggregate, alias)
        )

    # import_foreign = (
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="trade_energy",
    #         scope=TradeTypes.FOREIGN,
    #         direction="import",
    #         bus_carrier="AC",
    #     )
    #     # the trade statistic wrongly finds transmission between EU -> country buses.
    #     # Those are dropped by the filter_by statement.
    #     .pipe(filter_by, carrier=["AC", "DC"])
    #     .pipe(rename_aggregate, Group.import_foreign)
    # )
    # export_foreign = (
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="trade_energy",
    #         scope=TradeTypes.FOREIGN,
    #         direction="export",
    #         bus_carrier="AC",
    #     )
    #     # the trade statistic wrongly finds transmission between EU -> country buses.
    #     # Those are dropped by the filter_by statement.
    #     .pipe(filter_by, carrier=["AC", "DC"])
    #     .pipe(rename_aggregate, Group.export_foreign)
    # )
    # import_domestic = (
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="trade_energy",
    #         scope=TradeTypes.DOMESTIC,
    #         direction="import",
    #         bus_carrier="AC",
    #     )
    #     # the trade statistic wrongly finds transmission between EU -> country buses.
    #     # Those are dropped by the filter_by statement.
    #     .pipe(filter_by, carrier=["AC", "DC"])
    #     .pipe(rename_aggregate, Group.import_foreign)
    # )

    exporter = Exporter(
        statistics=[supply, demand] + trade_statistics,
        statistics_unit="MWh",
        view_config=config["view"],
    )

    exporter.defaults.plotly.chart = ESMGroupedBarChart
    exporter.defaults.plotly.xaxis_title = ""
    exporter.defaults.plotly.pattern = dict.fromkeys(
        [
            Group.import_foreign,
            Group.export_foreign,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )

    # todo: split PHS input and output
    # todo: split Transport input and output

    exporter.export(result_path, subdir)
