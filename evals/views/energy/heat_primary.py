"""Module for primary energy evaluations."""

from functools import partial
from pathlib import Path

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
    # we might be able to drop the input_energy statistics as a whole:
    # collect_myopic_statistics(networks, comps="Link", statistic="energy_balance").filter(like="methanolisation")
    # Out[11]:
    # year  location  carrier          bus_carrier
    # 2045  XK        methanolisation
    # AC                   -160316.433175
    # H2                   -673358.210901
    # co2 stored           -146742.357560
    # methanol              591679.859919
    # urban central heat     14792.579326
    #
    # do this and filter for carrier, that have one of the heat buses in the bus_carrier
    # index level. é voila, we should have everything we need?
    link_energy = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
        groupby=[partial(get_location, location_port=2), "carrier", "bus_carrier"],
    )

    carrier_with_heat_buses = []
    heat_buses = BusCarrier.heat_buses()
    for carrier, data in link_energy.groupby("carrier"):
        if len(data.index.unique("bus_carrier").intersection(heat_buses)):
            carrier_with_heat_buses.append(carrier)

    heat_links = filter_by(link_energy, carrier=carrier_with_heat_buses)

    heat_production = drop_from_multtindex_by_regex(heat_links, "DAC|water tanks")

    filter_by(heat_production, year="2045", location="EU").index.unique("carrier")

    link_input = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_input",
        bus_carrier=BusCarrier.heat_buses(),
    )
    link_input = drop_from_multtindex_by_regex(link_input, "DAC|water tanks")

    generator_input = collect_myopic_statistics(
        networks,
        statistic="energy_input",
        comps="Generator",
        bus_carrier=BusCarrier.heat_buses(),
    )
    generator_input = drop_from_multtindex_by_regex(generator_input, "heat vent")

    exporter = Exporter(
        statistics=[link_input, generator_input],
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
