# -*- coding: utf-8 -*-
"""Module for hydrogen nodal balances."""

from pathlib import Path

from evals.constants import (
    BusCarrier,
    DataModel,
    Group,
    TradeTypes,
)
from evals.fileio import Metric
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import make_evaluation_result_directories, rename_aggregate


def view_balance_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:  # numpydoc ignore=PR01
    """
    Evaluate the Hydrogen balance.

    Returns
    -------
    :

    Notes
    -----
    See eval module docstring for parameter description.
    """
    h2_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        bus_carrier=BusCarrier.H2,
    )
    pipelines = h2_production.filter(like="pipeline", axis=0).index.unique(
        DataModel.CARRIER
    )
    h2_production = h2_production.drop(pipelines, level=DataModel.CARRIER)

    h2_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        # comps=["Link", "Load"],
        bus_carrier=BusCarrier.H2,
    ).mul(-1)
    h2_demand = h2_demand.drop(pipelines, level=DataModel.CARRIER)

    h2_import_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="import",
        bus_carrier=BusCarrier.H2,
    ).filter(like="pipeline", axis=0)
    # todo: find missing supply energy
    h2_import_foreign = rename_aggregate(h2_import_foreign, Group.import_european)

    h2_export_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="export",
        bus_carrier=BusCarrier.H2,
    ).filter(like="pipeline", axis=0)
    h2_export_foreign = rename_aggregate(h2_export_foreign, Group.export_european)

    h2_import_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        bus_carrier=BusCarrier.H2,
    ).filter(like="pipeline", axis=0)
    h2_import_domestic = rename_aggregate(h2_import_domestic, Group.import_domestic)

    h2_export_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="export",
        bus_carrier=BusCarrier.H2,
    ).filter(like="pipeline", axis=0)
    h2_export_domestic = rename_aggregate(h2_export_domestic, Group.export_domestic)

    metric = Metric(
        statistics=[
            h2_production,
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
            Group.export_european,
            Group.import_european,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export(output_path, config["view"]["export"])
    metric.consistency_checks(config["view"])
