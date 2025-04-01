# -*- coding: utf-8 -*-
from pathlib import Path

from evals.constants import BusCarrier, Carrier, DataModel
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics


def view_electricity_capacities(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """
    Evaluate the optimal capacity for AC technologies.

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
    metric = Exporter(
        statistics=e_capas,
        statistics_unit=ac_production.attrs["unit"],
        view_config=view_config,
    )

    metric.defaults.plotly.chart = ESMBarChart
    metric.export(result_path, subdir)
