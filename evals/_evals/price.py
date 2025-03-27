# -*- coding: utf-8 -*-
"""Evaluate energy prices."""

from pathlib import Path

from constants import TITLE_SUFFIX, DataModel
from metric import Metric
from plots.barchart import ESMBarChart
from statistic import collect_myopic_statistics
from utils import make_evaluation_result_directories


def eval_capex(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the capital expenditure."""
    capex = collect_myopic_statistics(networks, statistic="capex").drop(
        "", level=DataModel.CARRIER, errors="ignore"
    )

    metric = Metric(
        metric_name="Capex",
        is_unit=capex.attrs["unit"],
        to_unit=capex.attrs["unit"],
        statistics=[capex],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "capex_{location}"
    metric.defaults.plotly.cutoff = 0.1

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_opex(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the operational expenditure."""
    opex = collect_myopic_statistics(networks, statistic="opex").drop(
        "", level=DataModel.CARRIER, errors="ignore"
    )

    metric = Metric(
        metric_name="Opex",
        is_unit=opex.attrs["unit"],
        to_unit=opex.attrs["unit"],
        statistics=[opex],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "opex_{location}"
    metric.defaults.plotly.cutoff = 0.1

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_market_value(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate market values."""
    market_value = collect_myopic_statistics(networks, statistic="market_value").drop(
        "", level=DataModel.CARRIER, errors="ignore"
    )

    metric = Metric(
        metric_name="MarketValue",
        is_unit=market_value.attrs["unit"],
        to_unit=market_value.attrs["unit"],
        statistics=[market_value],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "market_value_{location}"
    metric.defaults.plotly.cutoff = 0.1

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_revenue(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate revenues."""
    revenue = collect_myopic_statistics(networks, statistic="revenue").drop(
        "", level=DataModel.CARRIER, errors="ignore"
    )

    metric = Metric(
        metric_name="Revenue",
        is_unit=revenue.attrs["unit"],
        to_unit=revenue.attrs["unit"],
        statistics=[revenue],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "revenue_{location}"
    metric.defaults.plotly.cutoff = 0.1

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
