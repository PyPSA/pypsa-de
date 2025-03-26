# -*- coding: utf-8 -*-
"""Create a view for optimal heat production capacities.

The view shows one stacked bar per year for different groups of
technologies.
"""

from pathlib import Path

from evals.constants import TITLE_SUFFIX, BusCarrier, Group
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
        metric_name=view_config["name"],
        is_unit=heat_capacity.attrs["unit"],
        to_unit=view_config["unit"],
        statistics=[heat_capacity],
    )

    title = view_config["name"] + TITLE_SUFFIX
    metric.cfg.mapping = view_config["categories"]
    metric.cfg.excel.title = title
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = view_config["file_name"]
    metric.cfg.plotly.cutoff = view_config["cutoff"]
    metric.cfg.plotly.category_orders = (
        Group.storage_out,
        Group.solar_thermal,
        Group.ft,
        Group.chp_biomass,
        Group.fuel_cell_heat,
        Group.resistive_heater,
        Group.ch4_boiler,
        Group.chp_ch4,
        Group.chp_coal,
        Group.heat_pump,
        # demand in reversed order:
        Group.storage_in,
        Group.grid_losses,
        Group.hh_and_services_heat,
        Group.industry,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export(output_path, view_config["export"])
    metric.consistency_checks(view_config["checks"])
