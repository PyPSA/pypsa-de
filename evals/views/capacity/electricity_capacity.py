# -*- coding: utf-8 -*-
from pathlib import Path

import pandas as pd

from evals.constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel, Group
from evals.fileio import Metric
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
    )

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
    e_capas = [ac_generation_and_storage, ac_production]
    view_config = config["view"]
    metric = Metric(
        statistics_unit=e_capas[0].attrs["unit"],
        view_config=view_config,
        statistics=e_capas,
    )

    metric.defaults.plotly.chart = ESMBarChart

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export(output_path, view_config["export"])
    metric.consistency_checks(config["view"])
