"""Module for hydrogen nodal balances."""

from pathlib import Path

from evals.constants import (
    BusCarrier,
    DataModel,
    Group,
    TradeTypes,
)
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import rename_aggregate


def view_balance_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the Hydrogen balance.

    Returns
    -------
    :

    Notes
    -----
    See eval module docstring for parameter description.
    """
    h2_supply = collect_myopic_statistics(
        networks,
        statistic="supply",
        bus_carrier=BusCarrier.H2,
    )
    pipelines = h2_supply.filter(like="pipeline", axis=0).index.unique(
        DataModel.CARRIER
    )
    h2_supply = h2_supply.drop(pipelines, level=DataModel.CARRIER)

    h2_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=BusCarrier.H2,
        )
        .mul(-1)
        .drop(pipelines, level=DataModel.CARRIER)
    )

    h2_import_foreign = (
        collect_myopic_statistics(
            networks,
            statistic="trade_energy",
            scope=TradeTypes.FOREIGN,
            direction="import",
            bus_carrier=BusCarrier.H2,
        )
        .filter(like="pipeline", axis=0)
        .pipe(rename_aggregate, Group.import_foreign)
    )

    h2_export_foreign = (
        collect_myopic_statistics(
            networks,
            statistic="trade_energy",
            scope=TradeTypes.FOREIGN,
            direction="export",
            bus_carrier=BusCarrier.H2,
        )
        .filter(like="pipeline", axis=0)
        .pipe(rename_aggregate, Group.export_foreign)
    )

    h2_import_domestic = (
        collect_myopic_statistics(
            networks,
            statistic="trade_energy",
            scope=TradeTypes.DOMESTIC,
            direction="import",
            bus_carrier=BusCarrier.H2,
        )
        .filter(like="pipeline", axis=0)
        .pipe(rename_aggregate, Group.import_domestic)
    )

    h2_export_domestic = (
        collect_myopic_statistics(
            networks,
            statistic="trade_energy",
            scope=TradeTypes.DOMESTIC,
            direction="export",
            bus_carrier=BusCarrier.H2,
        )
        .filter(like="pipeline", axis=0)
        .pipe(rename_aggregate, Group.export_domestic)
    )

    metric = Exporter(
        statistics=[
            h2_supply,
            h2_demand,
            h2_import_foreign,
            h2_export_foreign,
            h2_import_domestic,
            h2_export_domestic,
        ],
        view_config=config["view"],
        statistics_unit="MWh",
    )

    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.pattern = dict.fromkeys(
        [
            Group.export_foreign,
            Group.import_foreign,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )

    metric.export(result_path, subdir)
