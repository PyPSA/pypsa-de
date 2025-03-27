# -*- coding: utf-8 -*-
from pathlib import Path

from evals.constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel
from evals.metric import Metric
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import make_evaluation_result_directories


def view_electricity_capacities(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the optimal capacity for AC technologies.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.
    """
    ac_generation_and_storage = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        comps=("Generator", "Store", "StorageUnit"),
        bus_carrier=BusCarrier.AC,
    ).drop("value of lost load", level=DataModel.CARRIER, errors="ignore")

    transmission_or_storage_links = ["", "DC", Carrier.v2g]
    ac_production = (
        collect_myopic_statistics(
            networks,
            statistic="optimal_capacity",
            comps="Link",
            bus_carrier=BusCarrier.AC,
        )
        .clip(lower=0)
        .drop(transmission_or_storage_links, level=DataModel.CARRIER, errors="ignore")
    )
    e_caps = [ac_generation_and_storage, ac_production]
    view_config = config["view"]
    metric = Metric(
        metric_name=view_config["name"],
        is_unit=e_caps[0].attrs["unit"],
        to_unit=view_config["unit"],
        statistics=e_caps,
    )
    title = view_config["name"] + TITLE_SUFFIX
    metric.cfg.mapping = view_config["categories"]
    metric.cfg.excel.title = title
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = view_config["file_name"]
    metric.cfg.plotly.cutoff = view_config["cutoff"]
    metric.cfg.plotly.category_orders = view_config["legend_order"]

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export(output_path, view_config["export"])
    metric.consistency_checks(view_config["checks"])
