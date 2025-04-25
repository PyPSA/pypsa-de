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
    link_energy_balance = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
    )

    # for every heat bus, calculate the amounts of supply for heat
    heat_supply = []
    for bus_carrier in BusCarrier.heat_buses():
        supply = (
            link_energy_balance.pipe(filter_for_carrier_connected_to, bus_carrier)
            .pipe(calculate_input_share, bus_carrier)
            # drop technology names in favour of input bus carrier names:
            .pipe(rename_aggregate, bus_carrier)
            .swaplevel(DataModel.BUS_CARRIER, DataModel.CARRIER)
            # .pipe(rename_aggregate, {"urban central water tanks discharger": "Storage"})
            # .mul(-1)  # need to reverse from input (=withdrawal) to bus supply
        )
        supply.index = supply.index.set_names(DataModel.YEAR_IDX_NAMES)
        supply.attrs["unit"] = "MWh_th"
        heat_supply.append(supply)

    generator_supply = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.heat_buses(),
    )
    generator_supply.attrs["unit"] = "MWh_th"
    heat_supply.append(generator_supply)

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

    For example, Electrolysis needs AC and H2. Therefore it is not possible to
    show Electricity or Hydrogen in the supply side of this view.
    """
    # heat_supply = [
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="supply",
    #         bus_carrier=BusCarrier.heat_buses(),
    #     )  # .pipe(rename_aggregate, {"urban central water tanks discharger": "Storage"})
    # ]

    heat_loss_factor = get_heat_loss_factor(networks)
    heat_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=BusCarrier.heat_buses(),
        )
        .pipe(split_urban_heat_losses_and_consumption, heat_loss_factor)
        # .pipe(rename_aggregate, {"urban central water tanks charger": "Storage"})
        .mul(-1)
    )

    exporter = Exporter(
        statistics=heat_supply + [heat_demand],
        view_config=config["view"],
    )

    # static view settings:
    exporter.defaults.plotly.chart = ESMGroupedBarChart
    exporter.defaults.plotly.xaxis_title = ""
    exporter.defaults.plotly.pattern = {"Demand": "/"}

    exporter.export(result_path, subdir)
