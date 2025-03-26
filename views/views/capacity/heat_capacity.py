# -*- coding: utf-8 -*-
"""Create a view for optimal heat production capacities.

The view shows one stacked bar per year for different groups of
technologies.
"""

from importlib import resources
from pathlib import Path

from esmtools.constants import TITLE_SUFFIX, BusCarrier, Carrier, Group
from esmtools.metric import Metric
from esmtools.plots.barchart import ESMBarChart
from esmtools.statistic import collect_myopic_statistics
from esmtools.utils import make_evaluation_result_directories

MAPPING = {
    # Carrier.dac: Group.dac,
    # Carrier.ft_1: Group.ft,
    # Carrier.ft_2: Group.ft,
    "DAC": "Direct Air Capture",
    "Fischer-Tropsch 1": "Fischer-Tropsch",
    "Fischer-Tropsch 2": "Fischer-Tropsch",
    "H2 Fuel Cell": "Fuel Cell (Heat)",
    "urban central air heat pump": "Heat Pump",
    "urban central coal CHP CC heat": "Coal CHP",
    "urban central coal CHP heat": "Coal CHP",
    "urban central gas CHP CC heat": "Gas CHP",
    "urban central gas CHP heat": "Gas CHP",
    "urban central gas boiler": "Gas Boiler",
    "urban central lignite CHP CC heat": "Coal CHP",
    "urban central lignite CHP heat": "Coal CHP",
    "urban central resistive heater": "Resistive Heater",
    "urban central solar thermal": "Solar Thermal",
    "urban central solid biomass CHP": "Biomass CHP",
    "urban central solid biomass CHP CC": "Biomass CHP",
    "urban central water tanks charger": "Storage In",
    "urban central water tanks discharger": "Storage Out",
}


def view_heat_capacity(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the optimal heat capacities to produce heat.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See esmtools.eval module docstring for parameter description.
    """
    heat_capacity = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        bus_carrier=BusCarrier.HEAT_URBAN_CENTRAL,
    ).clip(lower=0)

    metric = Metric(
        metric_name="District Heat Production Capacities",
        is_unit=heat_capacity.attrs["unit"],
        to_unit="GW",
        statistics=[heat_capacity],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.mapping = MAPPING
    metric.cfg.excel.title = title
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "heatKapas_{location}"
    metric.cfg.plotly.cutoff = 0.0001
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

    # from esmtools.utils import get_mapping
    #
    # mapping_all = get_mapping(metric.cfg.mapping)
    # mapping_used = {
    #     k: v for k, v in mapping_all.items() if k in metric.df.index.unique("carrier")
    # }

    output_path = make_evaluation_result_directories(result_path, subdir)
    # metric.export(output_path, config)
    metric.export_plotly(output_path)
    metric.assertion_checks(view_config["checks"])
