# -*- coding: utf-8 -*-
"""Create a view for optimal heat production capacities.

The view shows one stacked bar per year for different groups of
technologies.
"""

from pathlib import Path

from evals.constants import BusCarrier, Group
from evals.metric import Metric
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import make_evaluation_result_directories


def view_heat_capacity(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the optimal heat capacities to produce heat.

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

    view_config = config["view"]

    metric = Metric(
        statistics=[heat_capacity],
        statistics_unit=heat_capacity.attrs["unit"],
        view_config=view_config,
    )

    # constant view specific settings:
    metric.defaults.plotly.chart = ESMBarChart

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export(output_path, view_config["export"])
    metric.consistency_checks(view_config["checks"])
