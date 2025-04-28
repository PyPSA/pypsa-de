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
            # CO2 supply are CO2 emissions that do not help heat production
            .drop(["co2", "co2 stored"], level=DataModel.BUS_CARRIER)
            .pipe(calculate_input_share, bus_carrier)
            # drop technology names in favour of input bus carrier names:
            .pipe(rename_aggregate, bus_carrier)
            .swaplevel(DataModel.BUS_CARRIER, DataModel.CARRIER)
        )
        supply.index = supply.index.set_names(DataModel.YEAR_IDX_NAMES)
        supply.attrs["unit"] = "MWh_th"
        heat_supply.append(supply)

    heat_loss_factor = get_heat_loss_factor(networks)
    heat_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=BusCarrier.heat_buses(),
        )
        .pipe(split_urban_heat_losses_and_consumption, heat_loss_factor)
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
