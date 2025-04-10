"""Evaluate nodal heat balances."""

from pathlib import Path

from evals.constants import BusCarrier, DataModel
from evals.fileio import Exporter
from evals.plots.facetbars import ESMGroupedBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    calculate_input_share,
    filter_for_carrier_connected_to,
    get_heat_loss_factor,
    rename_aggregate,
    split_urban_heat_losses_and_consumption,
)


def view_balance_heat(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the heat balance.

    Parameters
    ----------
    result_path
    networks
    config
    subdir

    Returns
    -------
    :
    """
    # from evals.fileio import read_networks
    # from evals.utils import filter_by
    #
    # networks = read_networks(
    #     result_path="/IdeaProjects/pypsa-at/results/20240627public_db/8Gt_Bal_v3"
    # )
    #
    # bus_carrier = "urban central heat"
    # location = "FR0"
    # year = "2020"  # no low voltage -> only solid Biomass
    #
    # meth = collect_myopic_statistics(
    #     networks,
    #     statistic="energy_balance",
    #     # bus_carrier=BusCarrier.heat_buses(),
    # ).pipe(
    #     filter_by,
    #     location=location,
    #     year=year,
    #     carrier="methanolisation",
    # )
    #
    # by_tech = collect_myopic_statistics(
    #     networks,
    #     statistic="supply",
    #     bus_carrier=BusCarrier.heat_buses(),
    # ).pipe(
    #     filter_by,
    #     location=location,
    #     year=year,
    #     bus_carrier=bus_carrier,
    #     carrier="methanolisation",
    # )
    # print(by_tech)
    #
    # by_bus_carrier = (
    #     collect_myopic_statistics(
    #         networks,
    #         comps="Link",
    #         statistic="energy_balance",
    #     ).pipe(filter_for_carrier_connected_to, bus_carrier)
    #     .pipe(calculate_input_share, bus_carrier)
    #     .pipe(filter_by, location=location, year=year, carrier="methanolisation")
    # )
    # print(by_bus_carrier)
    #
    # # I think it will never be possible to fetch input energies without unit conversions.
    #
    # too small amounts yielded:
    link_energy_balance = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
    )

    # for every heat bus, calculate the amounts of supply for heat
    heat_supply = []
    for bus_carrier in BusCarrier.heat_buses():
        supply = (
            link_energy_balance.pipe(
                filter_for_carrier_connected_to, bus_carrier, kind="supply"
            )
            .pipe(calculate_input_share, bus_carrier)
            # drop technology names in favour of input bus carrier names:
            .pipe(rename_aggregate, bus_carrier)
            .swaplevel(DataModel.BUS_CARRIER, DataModel.CARRIER)
            # .pipe(rename_aggregate, {"urban central water tanks discharger": "Storage"})
            # .mul(-1)  # need to reverse from input (=withdrawal) to bus supply
        )
        supply.index = supply.index.set_names(DataModel.YEAR_IDX_NAMES)
        heat_supply.append(supply)

    heat_supply.append(
        collect_myopic_statistics(
            networks,
            statistic="supply",
            comps="Generator",
            bus_carrier=BusCarrier.heat_buses(),
        )
    )

    # disabled: because with ambient heat a surplus is generated on the suppy side
    # heat_supply.append(
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="ambient_heat",
    #     ).pipe(rename_aggregate, "Ambient Heat")
    # )

    """
    Yields correct balances.

    Notes
    -----
    This evaluation is only correct for Links (multi-port components)
    that have one input branch. The translation from technology name
    to category in done in the categories' mapper. If a technology has
    more than one input bus_carrier, this approach is insufficient. We
    cannot rename technologies with multiple bus carrier inputs to
    "Electricity", or "Hydrogen" without hiding the underlying demands.

    For example, Electrolsyis needs AC and H2. Therefore it is not possible to
    show Electricity or Hydrogen in the supply side of this view.
    """
    # heat_supply = [
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="supply",
    #         bus_carrier=BusCarrier.heat_buses(),
    #     ).pipe(rename_aggregate, {"urban central water tanks discharger": "Storage"})
    # ]

    heat_loss_factor = get_heat_loss_factor(networks)
    heat_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=BusCarrier.heat_buses(),
        )
        .pipe(split_urban_heat_losses_and_consumption, heat_loss_factor)
        .pipe(rename_aggregate, {"urban central water tanks charger": "Storage"})
        .mul(-1)
    )

    exporter = Exporter(
        statistics=heat_supply + [heat_demand],
        view_config=config["view"],
        statistics_unit="MWh",
    )

    # todo: check if ambient heat from heat pumps is already included in the input energy shares (= in low voltage)
    # todo: check if unit conversions might cause the problem. Looking at you: MWh_LHV
    # todo: check if CHPs produce energy and if those amounts cause the problem.

    # static view settings:
    exporter.defaults.plotly.chart = ESMGroupedBarChart
    exporter.defaults.plotly.xaxis_title = ""
    exporter.defaults.plotly.pattern = {"Demand": "/"}

    exporter.export(result_path, subdir)
