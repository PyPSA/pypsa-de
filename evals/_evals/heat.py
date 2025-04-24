"""Module for heat evaluations."""

from pathlib import Path

import pandas as pd
from constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel, Group
from metric import Metric
from plots.barchart import ESMBarChart
from plots.timeseries import ESMTimeSeriesChart
from statistic import collect_myopic_statistics
from utils import make_evaluation_result_directories, operations_override


def eval_district_heat_balance(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """
    Evaluate the district heat energy balance per country.

    Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See eval docstring for parameter description.
    """
    with operations_override(networks, "Load", "p_set"):
        # fixme: why p_set instead of p?!
        heat = collect_myopic_statistics(
            networks,
            statistic="energy_balance",
            bus_carrier=BusCarrier.HEAT_URBAN_CENTRAL,
        )

    heat_loss = get_heat_loss_factor(networks)
    heat = split_heat_hh_service_losses(heat, heat_loss)

    excluded_technologies = [
        Carrier.grid_losses,
        Carrier.water_tanks_charger_urban_central,
        Carrier.water_tanks_discharger_urban_central,
    ]
    heat = heat.drop(excluded_technologies, level=DataModel.CARRIER, errors="ignore")
    heat.columns = ["Energy Balance (MWh)"]

    metric = Metric(
        "Heat Balance", is_unit="MWh", to_unit="TWh", statistics=[heat.squeeze()]
    )

    metric.defaults.mapping = "district_heat"
    metric.defaults.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "heat_balance_{location}"
    metric.defaults.plotly.category_orders = (
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
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_district_heat_balance_ts(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the district heat energy balance time series.

    Writes 2 Excel files and 1 TimeSeriesChart per country and
    country.

    Notes
    -----
    See eval docstring for parameter description.
    """
    heat = collect_myopic_statistics(
        networks,
        statistic="energy_balance",
        bus_carrier=BusCarrier.HEAT_URBAN_CENTRAL,
        aggregate_time=False,
    )

    heat_loss = get_heat_loss_factor(networks)
    heat = split_heat_hh_service_losses(heat, heat_loss)

    metric = Metric(
        "District Heat Production and Demand",
        is_unit="MWh",
        to_unit="MWh",
        statistics=[heat],
    )

    metric.defaults.mapping = "district_heat"
    metric.defaults.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.plotly.chart = ESMTimeSeriesChart
    metric.defaults.plotly.plotby = [DataModel.LOCATION, DataModel.YEAR]
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]
    metric.defaults.plotly.file_name_template = "heat_prod_dem_time_{year}_{location}"
    metric.defaults.plotly.legend_header = "Production/Demand"

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_plotly(output_path)
