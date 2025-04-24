"""Export electricity storage Volumes."""

from pathlib import Path

from evals.constants import BusCarrier
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics


def view_capacity_ac_storage(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """Evaluate electricity Stores."""
    ac_storage = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        storage=True,
        bus_carrier=BusCarrier.ac_stores(),
    )

    metric = Exporter(
        statistics=[ac_storage], statistics_unit="MWh", view_config=config["view"]
    )

    metric.defaults.plotly.chart = ESMBarChart
    # prevent dropping empty years
    metric.defaults.plotly.cutoff_drop = False

    # metric.defaults.plotly.category_orders = (
    #     Group.reservoir,
    #     Group.battery_storage,
    #     Group.phs,
    # )

    # metric.defaults.plotly.footnotes = (
    #     " Further potential power storage volumes resulting from vehicle-to-grid "
    #     "options are not included here. <br> Storage volumes smaller than 1 TWh "
    #     "are not displayed. <br> Storages are located according to their market "
    #     "connection.",
    #     "",
    # )

    metric.export(result_path, subdir)
