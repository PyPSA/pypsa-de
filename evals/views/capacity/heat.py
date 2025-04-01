# -*- coding: utf-8 -*-
"""
Create a view for optimal heat production capacities.

The view shows one stacked bar per year for different groups of
technologies.
"""

from pathlib import Path

from evals.constants import BusCarrier
from evals.fileio import Metric
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics


def view_capacity_heat(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the optimal heat capacities to produce heat.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See eval module docstring for parameter description.
    """
    heat_capacity = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        bus_carrier=BusCarrier.HEAT_URBAN_CENTRAL,
    ).clip(lower=0)

    metric = Metric(
        statistics=[heat_capacity],
        statistics_unit=heat_capacity.attrs["unit"],
        view_config=config["view"],
    )

    # constant view specific settings
    metric.defaults.plotly.chart = ESMBarChart

    metric.export(result_path, subdir)
