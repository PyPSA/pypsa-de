"""Module for primary energy evaluations."""

from functools import partial
from pathlib import Path

import pandas as pd

from evals.constants import BusCarrier, DataModel
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics, get_location
from evals.utils import drop_from_multtindex_by_regex, filter_by


def view_heat_primary_energy(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the energy required for heat production and generation.

    Results are grouped by bus_carrier and not by carrier
    as usual to show the input energy carrier mix.

    Returns
    -------
    :

    Notes
    -----
    See eval docstring for parameter description.
    """
    link_energy = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
    )

    # drop non energy CO2 rows because they have mass unit and not energy
    # todo: is this really justified?
    link_energy = link_energy.drop(["co2", "co2 stored"], level=DataModel.BUS_CARRIER)

    # only keep Links that have at least one heat bus_carrier connected at one of their branches
    carrier_with_heat_buses = []
    heat_buses = BusCarrier.heat_buses()
    for carrier, data in link_energy.groupby(DataModel.CARRIER):
        if len(data.index.unique(DataModel.BUS_CARRIER).intersection(heat_buses)):
            carrier_with_heat_buses.append(carrier)
    heat_links = filter_by(link_energy, carrier=carrier_with_heat_buses)

    # drop storage technologies and DAC (direct air capture withdraws heat)
    heat_production = drop_from_multtindex_by_regex(heat_links, "DAC|water tanks")

    def _fuel_split(df):
        withdrawal = df[df.lt(0)]
        supply = df[df.ge(0)]
        heat_supply = filter_by(supply, bus_carrier=heat_buses).sum()
        return withdrawal * heat_supply / supply.sum()

    fuel_withdrawal = heat_production.groupby(
        [DataModel.YEAR, DataModel.LOCATION, DataModel.CARRIER], group_keys=False
    ).apply(_fuel_split)

    assert "EU" not in fuel_withdrawal.index.unique(DataModel.LOCATION)

    generator_supply = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.heat_buses(),
    )

    exporter = Exporter(
        statistics=[fuel_withdrawal.mul(-1), generator_supply],
        statistics_unit="MWh",
        view_config=config["view"],
    )

    # view specific static settings:
    exporter.defaults.plotly.chart = ESMBarChart
    exporter.defaults.excel.pivot_index = [DataModel.LOCATION, DataModel.BUS_CARRIER]
    exporter.defaults.plotly.plot_category = DataModel.BUS_CARRIER
    exporter.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]

    exporter.export(result_path, subdir=subdir)
