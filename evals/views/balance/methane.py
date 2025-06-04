"""Module for methane nodal balances."""

from pathlib import Path

import evals.plots as plots
from evals.constants import DataModel as DM
from evals.constants import (
    Group,
    TradeTypes,
)
from evals.fileio import Exporter
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    filter_by,
    get_storage_carriers,
    get_transmission_techs,
    rename_aggregate,
)


def view_balance_methane(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the methane balance.

    Returns
    -------
    :

    Notes
    -----
    See eval module docstring for parameter description.
    """
    bus_carrier = config["view"]["bus_carrier"]
    transmission_techs = get_transmission_techs(networks, bus_carrier)
    transmission_comps = [comp for comp, carr in transmission_techs]
    transmission_carrier = [carr for comp, carr in transmission_techs]
    storage_carrier = get_storage_carriers(networks) + config["view"].get(
        "storage_links", []
    )

    # unit = "MWh_LHV"
    # bus_carrier = ["gas", "gas primary", "biogas", "gas for industry"]

    supply = (
        collect_myopic_statistics(
            networks,
            statistic="supply",
            bus_carrier=bus_carrier,
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_out))
        .droplevel(DM.COMPONENT)
    )

    # pipelines = (
    #     supply.filter(like="pipeline", axis=0)
    #     .index.unique(DataModel.CARRIER)
    #     .pipe(rename_aggregate, {"gas Store": Group.storage_out})
    # )
    # supply = supply.drop(pipelines, level=DataModel.CARRIER, errors="ignore")
    # supply.attrs["unit"] = unit  # renewable gas lacks unit (unit is '')

    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            bus_carrier=bus_carrier,
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        # .drop(transmission_techs, level=DM.CARRIER)
        # .pipe(rename_aggregate, dict.fromkeys(transmission_carrier, "Export"))
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_in))
        .mul(-1)
        .droplevel(DM.COMPONENT)
    )

    trade_statistics = []
    if any(transmission_techs):
        for scope, direction, alias in [
            (TradeTypes.FOREIGN, "import", Group.import_foreign),
            (TradeTypes.FOREIGN, "export", Group.export_foreign),
            (TradeTypes.DOMESTIC, "import", Group.import_domestic),
            (TradeTypes.DOMESTIC, "export", Group.export_domestic),
        ]:
            trade = (
                collect_myopic_statistics(
                    networks,
                    statistic="trade_energy",
                    scope=scope,
                    direction=direction,
                    bus_carrier=bus_carrier,
                    aggregate_components=None,
                )
                # the trade statistic wrongly finds transmission between EU -> country buses.
                # Those are dropped by the filter_by statement.
                .pipe(
                    filter_by,
                    component=transmission_comps,
                    carrier=transmission_carrier,
                )
                .pipe(rename_aggregate, alias)
                .droplevel(DM.COMPONENT)
            )
            trade.attrs["unit"] = "MWh_LHV"
            trade_statistics.append(trade)

    exporter = Exporter(
        statistics=[supply, demand] + trade_statistics,
        view_config=config["view"],
    )

    exporter.defaults.plotly.chart = getattr(plots, config["view"]["chart"])
    if exporter.defaults.plotly.chart == plots.ESMGroupedBarChart:
        exporter.defaults.plotly.xaxis_title = ""
    elif exporter.defaults.plotly.chart == plots.ESMBarChart:
        # combine bus carrier to export netted technologies, although
        # they have difference bus_carrier in index , e.g.
        # electricity distribution grid, (AC, low voltage)
        exporter.statistics[0] = rename_aggregate(demand, "gas", level="bus_carrier")
        exporter.statistics[1] = rename_aggregate(supply, "gas", level="bus_carrier")

    exporter.export(result_path, config["global"]["subdir"])
