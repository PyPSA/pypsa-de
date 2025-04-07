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
    to_concat = []
    eu_carrier = set()
    for location_port in ("0", "1"):
        # some Links draw from buses with EU location at port 0, others
        # from port 1. We calculate both, drop EU from both results, merge
        # the results and drop duplicates.
        # todo: shouldn't this be done most of the time or always even?
        #   idea -> grouper that returns the location but with EU fixed to localized locations
        _link_energy_at_location_port = collect_myopic_statistics(
            networks,
            comps="Link",
            statistic="energy_balance",
            groupby=[
                partial(get_location, location_port=location_port),
                "carrier",
                "bus_carrier",
            ],
        )
        to_concat.append(_link_energy_at_location_port)
        _eu_carrier_at_location_port = filter_by(
            _link_energy_at_location_port, location="EU"
        ).index.unique("carrier")
        eu_carrier.union(_eu_carrier_at_location_port)

    link_energy = (
        pd.concat(to_concat).drop_duplicates().drop("EU", level=DataModel.LOCATION)
    )

    # test if some carrier were lost by dropping EU previously
    # eu_carrier = eu_carrier_port_0.union(eu_carrier_port_1)
    assert not any(eu_carrier.difference(link_energy.index.unique("carrier")))

    # drop non energy CO2 rows because they have mass unit and not energy
    link_energy = link_energy.drop(["co2", "co2 stored"], level="bus_carrier")

    # only keep Links that have at least one heat bus_carrier connected at one of their branches
    carrier_with_heat_buses = []
    heat_buses = BusCarrier.heat_buses()
    for carrier, data in link_energy.groupby("carrier"):
        if len(data.index.unique("bus_carrier").intersection(heat_buses)):
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
