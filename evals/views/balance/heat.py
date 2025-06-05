"""Evaluate nodal heat balances."""

from pathlib import Path

import pandas as pd

import evals.plots as plots
from evals.constants import DataModel as DM
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
) -> None:
    """
    Evaluate the heat balance.

    Returns
    -------
    :
    """
    bus_carrier = config["view"]["bus_carrier"]
    # todo: storage links

    link_energy_balance = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_balance",
    )

    # for every heat bus, calculate the amounts of supply for heat
    to_concat = []
    for bc in bus_carrier:
        p = (
            link_energy_balance.pipe(filter_for_carrier_connected_to, bc)
            # CO2 supply are CO2 emissions that do not help heat production
            .drop(["co2", "co2 stored"], level=DM.BUS_CARRIER)
            .pipe(calculate_input_share, bc)
            # drop technology names in favour of input bus carrier names:
            .pipe(rename_aggregate, bc)
            .swaplevel(DM.BUS_CARRIER, DM.CARRIER)
        )
        p.index = p.index.set_names(DM.YEAR_IDX_NAMES)
        p.attrs["unit"] = "MWh_th"
        to_concat.append(p)

    supply = pd.concat(to_concat)

    heat_loss_factor = get_heat_loss_factor(networks)
    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=bus_carrier,
        )
        .pipe(split_urban_heat_losses_and_consumption, heat_loss_factor)
        .mul(-1)
    )

    exporter = Exporter(statistics=[supply, demand], view_config=config["view"])

    # static view settings:
    exporter.defaults.plotly.chart = ESMGroupedBarChart
    exporter.defaults.plotly.xaxis_title = ""
    exporter.defaults.plotly.pattern = {"Demand": "/"}

    exporter.export(result_path, config["global"]["subdir"])
    chart_class = getattr(plots, config["view"]["chart"])
    exporter.defaults.plotly.chart = chart_class

    if chart_class == plots.ESMGroupedBarChart:
        exporter.defaults.plotly.xaxis_title = ""
    elif chart_class == plots.ESMBarChart:
        # combine bus carrier to export netted technologies, although
        # they have difference bus_carrier in index , e.g.
        # electricity distribution grid, (AC, low voltage)
        exporter.statistics[0] = rename_aggregate(
            demand, bus_carrier[0], level=DM.BUS_CARRIER
        )
        exporter.statistics[1] = rename_aggregate(
            supply, bus_carrier[0], level=DM.BUS_CARRIER
        )

    exporter.export(result_path, config["global"]["subdir"])
