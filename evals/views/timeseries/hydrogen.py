# -*- coding: utf-8 -*-
"""Export time series views for hydrogen."""

from pathlib import Path

from evals.constants import (
    BusCarrier,
    DataModel,
    Group,
    TradeTypes,
)
from evals.fileio import Exporter
from evals.plots.timeseries import ESMTimeSeriesChart
from evals.statistic import collect_myopic_statistics
from evals.utils import rename_aggregate


def view_timeseries_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """Evaluate the Hydrogen balance time series."""
    h2_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    )
    pipelines = h2_production.filter(like="pipeline", axis=0).index.unique(
        DataModel.CARRIER
    )
    h2_production = h2_production.drop(pipelines, level=DataModel.CARRIER)
    h2_production = rename_aggregate(
        h2_production, {"H2 Store": Group.storage_out}, level=DataModel.CARRIER
    )

    h2_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=BusCarrier.H2,
            aggregate_time=False,
        )
        .drop(pipelines, level=DataModel.CARRIER)
        .mul(-1.0)
    )
    h2_demand = rename_aggregate(
        h2_demand, {"H2 Store": Group.storage_in}, level=DataModel.CARRIER
    )

    trade_saldo = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=(TradeTypes.FOREIGN, TradeTypes.DOMESTIC),
        direction="saldo",
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    )
    # todo: trade statistic finds transformation Links connected to EU buses, e.g. "EU NH3" or "EU renewable gas"
    trade_saldo = trade_saldo.filter(like="pipeline", axis=0)

    trade_saldo = rename_aggregate(trade_saldo, trade_saldo.attrs["name"])

    metric = Exporter(
        statistics=[h2_production, h2_demand, trade_saldo],
        statistics_unit="MWh",
        view_config=config["view"],
    )

    # view specific settings
    metric.defaults.excel.chart = None
    metric.defaults.plotly.chart = ESMTimeSeriesChart
    metric.defaults.plotly.plotby = [DataModel.YEAR, DataModel.LOCATION]
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]
    metric.defaults.plotly.xaxis_title = ""
    # metric.defaults.plotly.legend_header = "Production/Demand"
    metric.defaults.plotly.cutoff = 0.001  # MWh
    metric.defaults.plotly.pattern = dict.fromkeys(
        [Group.export_net, Group.import_net, Group.import_global], "/"
    )

    metric.export(result_path, subdir)
