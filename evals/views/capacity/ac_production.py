# -*- coding: utf-8 -*-
"""Export Electricity production capacities."""

from pathlib import Path

from evals.constants import BusCarrier, DataModel
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics


def view_electricity_capacities(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the optimal capacity for AC technologies.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.
    """
    ac_capacity = (
        collect_myopic_statistics(
            networks,
            statistic="optimal_capacity",
            comps=("Generator", "Link", "StorageUnit"),
            bus_carrier=BusCarrier.AC,
            # aggregate_components=None,
        ).clip(lower=0)
        # drop Links connected to StorageUnits and DC Links
        .drop(
            [
                "battery discharger",
                "DC",
                "battery charger",
                "electricity distribution grid",
            ],
            level=DataModel.CARRIER,
            errors="ignore",
        )
    )
    ac_capacity = ac_capacity[ac_capacity > 0]

    # solar-hsat: correct, map to PV
    # lignite + coal + biogas + biomass: correct, produces AC from fuel
    # carr = "OCGT"
    # networks["2040"].static("Link").query("carrier == @carr").filter(like="bus").T
    # networks["2040"].static("Link").query("carrier == @carr").filter(like="eff").T

    # ac_capacity.to_frame().query("component == 'Link'")  # & carrier == 'AC'
    # ac_capacity.to_frame().query("carrier == 'DC'")  # & carrier == 'AC'

    #
    # transmission_or_storage_links = ["", "DC", Carrier.v2g]
    # ac_production = (
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="optimal_capacity",
    #         comps="Link",
    #         bus_carrier=BusCarrier.AC,
    #     )
    #     .clip(lower=0)
    #     .drop(transmission_or_storage_links, level=DataModel.CARRIER, errors="ignore")
    # )

    metric = Exporter(
        # statistics=[ac_generation_and_storage, ac_production],
        statistics=[ac_capacity],
        statistics_unit="MWh",
        view_config=config["view"],
    )

    metric.defaults.plotly.chart = ESMBarChart
    metric.export(result_path, subdir)
