# -*- coding: utf-8 -*-
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

    metric.cfg.mapping = "district_heat"
    metric.cfg.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "heat_balance_{location}"
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

    metric.cfg.mapping = "district_heat"
    metric.cfg.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.plotly.chart = ESMTimeSeriesChart
    metric.cfg.plotly.plotby = [DataModel.LOCATION, DataModel.YEAR]
    metric.cfg.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]
    metric.cfg.plotly.file_name_template = "heat_prod_dem_time_{year}_{location}"
    metric.cfg.plotly.legend_header = "Production/Demand"

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_plotly(output_path)


def split_heat_hh_service_losses(
    df: pd.DataFrame | pd.Series, heat_loss: int
) -> pd.DataFrame:
    """Split urban heat amounts by a heat loss factor.

    Amounts for urban central and decentral heat contain
    distribution losses. However, the evaluation shows final demands
    in the results. Therefore, heat network distribution losses need
    to be separated from the total amounts because grid distribution
    losses do not arrive at the metering endpoint.

    Parameters
    ----------
    df
        The input data frame with values for urban central heat
        technologies.
    heat_loss
        The heat loss factor from the configuration file.

    Returns
    -------
    :
        The data frame with split heat amounts for end user demand
        (hh and service), distribution grid losses (grid losses) and
        anything else from the input data frame
        (not urban central heat).
    """
    if isinstance(df, pd.Series):
        df = df.to_frame(f"{df.attrs['name']} ({df.attrs['unit']})")

    # distinguish between grid losses and final consumption
    loss_factor = heat_loss / (1 + heat_loss)
    urban_heat_bus_carrier = [
        BusCarrier.HEAT_URBAN_CENTRAL,
        # fixme: Are urban_services and urban_residential correct?
        #  I thought, that decentral technologies do not have a
        #  distribution grid, hence no grid losses.
        BusCarrier.HEAT_URBAN_SERVICES,
        BusCarrier.HEAT_URBAN_RESIDENTIAL,
    ]
    central_heat = df.query(f"{DataModel.CARRIER}.isin(@urban_heat_bus_carrier)")

    grid_losses = central_heat.mul(loss_factor)
    grid_losses_mapper = dict.fromkeys(urban_heat_bus_carrier, Carrier.grid_losses)
    grid_losses = grid_losses.rename(grid_losses_mapper, level=DataModel.CARRIER)

    hh_services = central_heat.mul(1 - loss_factor)
    hh_services_mapper = dict.fromkeys(urban_heat_bus_carrier, Carrier.hh_and_services)
    hh_services = hh_services.rename(hh_services_mapper, level=DataModel.CARRIER)

    heat_no_urban = df.drop(
        urban_heat_bus_carrier, level=DataModel.CARRIER, errors="ignore"
    )

    return pd.concat([heat_no_urban, hh_services, grid_losses]).sort_index()


def get_heat_loss_factor(networks: dict) -> int:
    """Return the heat loss factor for district heating from the config.

    Parameters
    ----------
    networks
        The loaded networks.

    Returns
    -------
    The heat loss factor for district heating networks.
    """
    heat_loss_factors = {
        n.meta["sector"]["district_heating_loss"] for n in networks.values()
    }
    assert len(heat_loss_factors) == 1, "Varying loss factors are not supported."
    return heat_loss_factors.pop()
