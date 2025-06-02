"""
Create a view for optimal heat production capacities.

The view shows one stacked bar per year for different groups of
technologies.
"""

from pathlib import Path

from evals.constants import BusCarrier
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics


def view_capacity_heat_production(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the optimal heat capacities to produce heat.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See eval module docstring for parameter description.
    """
    heat_capacity = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        bus_carrier=BusCarrier.HEAT_URBAN_CENTRAL,
    ).clip(lower=0)

    # correct unit to Power
    heat_capacity.attrs["unit"] = heat_capacity.attrs["unit"].replace("MWh", "MW")

    exporter = Exporter(statistics=[heat_capacity], view_config=config["view"])

    # constant view specific settings
    exporter.defaults.plotly.chart = ESMBarChart

    exporter.export(result_path, subdir)
