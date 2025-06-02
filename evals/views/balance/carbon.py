"""Module for carbon dioxide nodal balances."""

from pathlib import Path

from evals.constants import (
    DataModel,
    Group,
    TradeTypes,
)
from evals.fileio import Exporter
from evals.plots import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import rename_aggregate


def view_balance_carbon(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the carbon balance.

    Returns
    -------
    :
    """
    unit = "t_co2"
    # storage and sequestration
    bus_carrier = ["co2", "co2 sequestered", "co2 stored"]
    # bus_carrier = "co2"

    supply = collect_myopic_statistics(
        networks,
        statistic="supply",
        bus_carrier=bus_carrier,
    )
    pipelines = supply.filter(like="pipeline", axis=0).index.unique(DataModel.CARRIER)
    supply = supply.drop(pipelines, level=DataModel.CARRIER)
    supply.attrs["unit"] = unit

    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=bus_carrier,
        )
        .mul(-1)
        .drop(pipelines, level=DataModel.CARRIER)
    )
    demand.attrs["unit"] = unit

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
            .filter(like="pipeline", axis=0)
            .pipe(rename_aggregate, alias)
        )
        trade.attrs["unit"] = unit
        trade_statistics.append(trade)

    exporter = Exporter(
        statistics=[supply, demand] + trade_statistics,
        view_config=config["view"],
    )
    # todo: split storage in and storage out

    # exporter.defaults.plotly.chart = ESMGroupedBarChart
    exporter.defaults.plotly.chart = ESMBarChart
    exporter.defaults.plotly.xaxis_title = ""
    exporter.defaults.plotly.pattern = dict.fromkeys(
        [
            Group.export_foreign,
            Group.import_foreign,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )

    exporter.export(result_path, subdir)
