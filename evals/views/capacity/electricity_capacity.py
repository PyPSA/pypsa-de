# -*- coding: utf-8 -*-
from pathlib import Path

import pandas as pd

from evals.constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel, Group
from evals.metric import Metric
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import make_evaluation_result_directories


def view_electricity_capacities(
    result_path: str | Path,
    networks: dict,
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

    metric = Metric(
        metric_name="Optimal Capacity Electricity",
        is_unit="MW",
        to_unit="GW",
        statistics=[ac_generation_and_storage, ac_production],
    )

    metric.cfg.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "eKapas_{location}"
    # metric.cfg.plotly.cutoff = 0.0
    metric.cfg.plotly.category_orders = (
        Group.battery_storage,
        Group.nuclear_power,
        Group.pp_thermal,
        Group.phs,
        Group.reservoir,
        Group.ror,
        Group.wind,
        Group.pv,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
