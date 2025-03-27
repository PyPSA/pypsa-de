# -*- coding: utf-8 -*-
"""Industry Evaluations."""

from pathlib import Path

import pandas as pd
from constants import TITLE_SUFFIX, DataModel, Group
from fileio import prepare_industry_demand
from metric import Metric
from plots.barchart import ESMBarChart
from plots.facetbars import ESMGroupedBarChart
from utils import make_evaluation_result_directories, rename_aggregate

INDUSTRY_BUS_CARRIER = {
    "electricity": "Electricity",
    "methane": "Methane",
    "low-temperature heat": "District Heat",
    "hard coal": "Coal",
    "Fischer-Tropsch": "Oil",
    "solid biomass": "Biomass",
    "hydrogen": "Hydrogen",
    "naphtha": "Oil",
    "gas feedstock": "Methane",
    "other feedstock": "Oil",
}

INDUSTRY_SECTORS_SHORT = {
    "Iron and Steel": "Iron/Steel",
    "Non-metallic Minerals": "N.M. Minerals",
    "Chemical": "Chemicals",
    "Paper, Pulp and Print": "Paper",
    "Non-ferrous Metals": "N.F. Metals",
    "Food and Tabacco": "Food",
    "Transport Equipment": "Transport",
    "Machinery": "Machinery",
    "Textile and Leather": "Textile",
    "Wood and Wood Products": "Wood",
    "Non-specified Industry": "Other",
}


def eval_industry_total(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the energy demand by carrier.

    Results are grouped by bus_carrier (and not by carrier
    technology as usual) to show the source energy mix.

    Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See pacakge docstring for parameter description.
    """
    industry_demand = fetch_industry_demand_statistics(networks, result_path)

    metric = Metric(
        "Industrial Demand by Carrier Total",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[industry_demand],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.mapping = INDUSTRY_BUS_CARRIER
    metric.cfg.excel.title = title
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "ind_sec_tot_{location}"
    metric.cfg.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]
    metric.cfg.plotly.category_orders = (
        Group.oil,
        Group.ch4,
        Group.h2,
        Group.electrictiy,
        Group.heat_district,
        Group.coal,
        Group.biomass,
    )
    metric.cfg.plotly.footnotes = (
        "Methane demand includes additional demand for Carbon Capture if used "
        "by the model. <br>Industry demand includes energetic and non-energetic "
        "demands.",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_industry_sectoral(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the energy demand by carrier and sector.

    Results are grouped by bus_carrier (and not by carrier
    technology as usual) to show the source energy mix.

    Writes 2 Excel files and 1 GroupedBarChart per country.

    Notes
    -----
    See eval docstring for parameter description.
    """
    industry_demand = fetch_industry_demand_statistics(networks, result_path)
    # use abbreviated sector names
    industry_demand = rename_aggregate(
        industry_demand, INDUSTRY_SECTORS_SHORT, level=DataModel.BUS_CARRIER
    )

    metric = Metric(
        "Industrial Demand by Carrier Sectoral",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[industry_demand],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.mapping = INDUSTRY_BUS_CARRIER

    metric.cfg.excel.chart_title = title
    metric.cfg.excel.pivot_index = [DataModel.LOCATION, DataModel.CARRIER]
    metric.cfg.excel.pivot_columns = [DataModel.BUS_CARRIER, DataModel.YEAR]
    metric.cfg.excel.axis_labels = [DataModel.YEAR.title(), metric.df.attrs["unit"]]

    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMGroupedBarChart
    metric.cfg.plotly.file_name_template = "ind_sec_dem_{location}"
    metric.cfg.plotly.cutoff = 0.0  # None in original Toolbox
    metric.cfg.plotly.xaxis_title = ""
    metric.cfg.plotly.category_orders = (
        Group.oil,
        Group.ch4,
        Group.h2,
        Group.electrictiy,
        Group.biomass,
        Group.heat_district,
        Group.coal,
    )
    metric.cfg.plotly.footnotes = (
        "N.F. Metals = Non-ferrous Metals; N.M. Minerals = Non-metallic Minerals; "
        "Paper = Pulp, Paper and Printing; Textile = Textile and Leather; Transport"
        " = Transport Equipment; Wood = Wood and Wood Products",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def fetch_industry_demand_statistics(
    networks: dict, result_path: str | Path
) -> pd.Series:
    """Calculate the industry demand statistics.

    fixme: AT Hydrogen Demand in 2020 and sector Industry is -0.06 TWh
           in old and in new evaluation. Could be a problem in the
           interpolated data file.

    Parameters
    ----------
    networks
        The loaded networks.
    result_path
        The path to the results directory, needed to locate resource
        files.

    Returns
    -------
    :
        A collection of industry demand statistics.
    """
    non_demand = [
        "current electricity",
        "process emission",
        "process emission from feedstock",
    ]
    industry_demand = prepare_industry_demand(
        result_path, networks, sub_directory="interpolated_data"
    ).drop(non_demand, level=DataModel.CARRIER)

    return industry_demand
