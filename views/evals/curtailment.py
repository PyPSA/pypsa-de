"""Module for curtailment evaluations."""

from pathlib import Path

from esmtools.constants import TITLE_SUFFIX, DataModel
from esmtools.metric import Metric
from esmtools.plots.barchart import ESMBarChart
from esmtools.statistic import collect_myopic_statistics
from esmtools.utils import make_evaluation_result_directories


def eval_capacity_factor(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate capacity factors by country and year."""
    capacity_factor = collect_myopic_statistics(networks, statistic="capacity_factor")

    metric = Metric(
        metric_name="CapacityFactor",
        is_unit="MWh",
        to_unit="GWh",
        statistics=[capacity_factor],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.excel.chart_title = title
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "capacity_factor_{location}"
    metric.cfg.plotly.cutoff = 0.1

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_curtailment(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the curtailment by country and year.

    Curtailment ("Abregelung") is the amount of energy that was not
    dispatched by generators, although it could have been dispatched.
    I.e. it is defined as the difference between "p_max_pu * p_nom_opt"
    and the "p" time series.

    Returns
    -------
    :
        Writes X Excel files and X XXXX Chart per country.

    Notes
    -----
    See esmtools.eval docstring for parameter description.
    """
    curtailment = collect_myopic_statistics(networks, statistic="curtailment")

    # Exclude carriers with very large values because they
    # distort plots otherwise
    curtailment = curtailment.drop(
        "urban central solar thermal", level=DataModel.CARRIER, errors="ignore"
    )

    metric = Metric(
        metric_name="Curtailment",
        is_unit="MWh",
        to_unit="GWh",
        statistics=[curtailment],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.excel.chart_title = title
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "curtailment_{location}"
    metric.cfg.plotly.cutoff = 0.1

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
