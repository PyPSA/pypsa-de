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

    # for location_port in ("0", "1"):
    link_energy_port0 = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
        groupby=[
            partial(get_location, location_port="0"),
            "carrier",
            "bus_carrier",
            "unit",
        ],
    )
    eu_carrier_port_0 = filter_by(link_energy_port0, location="EU").index.unique(
        "carrier"
    )

    link_energy_port1 = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
        groupby=[
            partial(get_location, location_port="1"),
            "carrier",
            "bus_carrier",
            "unit",
        ],
    )
    eu_carrier_port_1 = filter_by(link_energy_port1, location="EU").index.unique(
        "carrier"
    )

    link_energy = (
        pd.concat([link_energy_port0, link_energy_port1])
        .drop_duplicates()
        .drop("EU", level=DataModel.LOCATION)
    )

    # test if some carrier were lost by dropping EU previously
    eu_carrier = eu_carrier_port_0.union(eu_carrier_port_1)
    assert not any(eu_carrier.difference(link_energy.index.unique("carrier")))

    # drop non energy rows, such as CO2
    energy_units = [u for u in link_energy.index.unique("unit") if u.startswith("MWh")]
    link_energy = filter_by(link_energy, unit=energy_units).droplevel("unit")

    # only keep Links that have at least one heat bus_carrier connected at one of their branches
    carrier_with_heat_buses = []
    heat_buses = BusCarrier.heat_buses()
    for carrier, data in link_energy.groupby("carrier"):
        if len(data.index.unique("bus_carrier").intersection(heat_buses)):
            carrier_with_heat_buses.append(carrier)
    heat_links = filter_by(link_energy, carrier=carrier_with_heat_buses)

    # drop storage technology Links
    heat_production = drop_from_multtindex_by_regex(heat_links, "water tanks")

    # # need to drop bus_carrier with non energy units, e.g. CO2, to prevent mixing them with MWh
    # heat_production = drop_from_multtindex_by_regex(
    #     heat_production, "co2", level=DataModel.BUS_CARRIER
    # )

    # # collect excluded carrier to assert no energies are lost in the final result.
    # eu_carrier.update(
    #     set(filter_by(heat_production, location="EU").index.unique("carrier"))
    # )

    # assert "EU" not in heat_production.index.unique("location"), (
    #     f"Must localize energy withdrawal for bus_carriers: {filter_by(heat_production, location='EU').index.unique('carrier')}"
    # )

    # assert all(c in link_energy.index.unique("carrier") for c in eu_carrier), (
    #     "Some carrier were included with EU location, but are missing in the result: ..."
    # )

    # for every nodal Link that has heat supply, calculate the share of withdrawal carrier energies
    # todo: refactor to apply(func) for performance
    to_concat = []
    for idx, df in heat_production.groupby(["year", "location", "carrier"]):
        # cases:
        # one input -> multiple outputs:
        #       * calculate efficiency share for outputs (for heat)
        #       * multiply by input. done
        # multiple inputs -> multiple outputs:
        #       * calculate efficiency share for outputs (for heat)
        #       * multiply by all inputs. done
        withdrawal = df[df.lt(0)]
        supply = df[df.ge(0)]

        if len(supply.index.unique("bus_carrier").intersection(heat_buses)):
            heat_supply = filter_by(supply, bus_carrier=heat_buses)
            heat_share = heat_supply.item() / supply.sum()
            withdrawal_for_heat = withdrawal * heat_share
            to_concat.append(withdrawal_for_heat)

    primary_energy = pd.concat(to_concat)

    assert not ("EU" in primary_energy.index.unique("location"))

    # link_input = collect_myopic_statistics(
    #     networks,
    #     comps="Link",
    #     statistic="energy_input",
    #     bus_carrier=BusCarrier.heat_buses(),
    # )
    # link_input = drop_from_multtindex_by_regex(link_input, "DAC|water tanks")

    generator_input = collect_myopic_statistics(
        networks,
        statistic="energy_input",
        comps="Generator",
        bus_carrier=BusCarrier.heat_buses(),
    )
    generator_input = drop_from_multtindex_by_regex(generator_input, "heat vent")

    exporter = Exporter(
        statistics=[heat_production, generator_input],
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
