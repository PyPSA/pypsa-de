"""Module for gas storage evaluations."""

from pathlib import Path

from evals.constants import BusCarrier
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics


def view_gas_storage_capacities(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate optimal storage capacities for CH4 and H2.

    Returns
    -------
    :
    """
    gas_stores = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        storage=True,
        bus_carrier=[BusCarrier.H2],
    )

    metric = Exporter(
        statistics=[gas_stores],
        statistics_unit="MWh",
        view_config=config["view"],
    )

    metric.defaults.plotly.chart = ESMBarChart
    # prevent dropping empty years
    metric.defaults.plotly.cutoff_drop = False

    metric.export(result_path, subdir)
