"""Common functions to generate balances."""

from evals import plots as plots
from evals.constants import DataModel as DM
from evals.constants import Group, TradeTypes
from evals.fileio import Exporter
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    filter_by,
    get_storage_carriers,
    get_transmission_techs,
    rename_aggregate,
)


def simple_bus_balance(
    networks: dict,
    config: dict,
    result_path,
) -> None:
    bus_carrier = config["view"]["bus_carrier"]
    transmission_techs = get_transmission_techs(networks, bus_carrier)
    transmission_comps = [comp for comp, carr in transmission_techs]
    transmission_carrier = [carr for comp, carr in transmission_techs]
    storage_carrier = get_storage_carriers(networks) + config["view"].get(
        "storage_links", []
    )

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
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_in))
        .mul(-1)
        .droplevel(DM.COMPONENT)
    )

    trade_statistics = []
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
        trade.attrs["unit"] = supply.attrs["unit"]
        trade_statistics.append(trade)

    exporter = Exporter(
        statistics=[supply, demand] + trade_statistics,
        view_config=config["view"],
    )

    chart_class = getattr(plots, config["view"]["chart"])
    exporter.defaults.plotly.chart = chart_class

    if chart_class == plots.ESMGroupedBarChart:
        exporter.defaults.plotly.xaxis_title = ""
    elif chart_class == plots.ESMBarChart:
        # combine bus carrier to export netted technologies, although
        # they have difference bus_carrier in index , e.g.
        # electricity distribution grid, (AC, low voltage)
        exporter.statistics[0] = rename_aggregate(
            demand, bus_carrier[0], level="bus_carrier"
        )
        exporter.statistics[1] = rename_aggregate(
            supply, bus_carrier[0], level="bus_carrier"
        )

    exporter.export(result_path, config["global"]["subdir"])
