from pathlib import Path

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
    (
        bus_carrier,
        transmission_comps,
        transmission_carrier,
        storage_carrier,
    ) = _parse_view_config_items(networks, config)

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
        exporter.statistics = [
            rename_aggregate(s, bus_carrier[0], level=DM.BUS_CARRIER)
            for s in exporter.statistics
        ]

    exporter.export(result_path, config["global"]["subdir"])


def simple_timeseries(
    networks: dict,
    config: dict,
    result_path: str | Path,
) -> None:
    """Export simple time series views."""
    (
        bus_carrier,
        transmission_comps,
        transmission_carrier,
        storage_carrier,
    ) = _parse_view_config_items(networks, config)

    supply = (
        collect_myopic_statistics(
            networks,
            statistic="supply",
            bus_carrier=bus_carrier,
            aggregate_time=False,
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
            aggregate_time=False,
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_in))
        .droplevel(DM.COMPONENT)
        .mul(-1)
    )

    trade_saldo = (
        collect_myopic_statistics(
            networks,
            statistic="trade_energy",
            scope=(TradeTypes.FOREIGN, TradeTypes.DOMESTIC),
            direction="saldo",
            bus_carrier=bus_carrier,
            aggregate_time=False,
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
        )
        .droplevel(DM.COMPONENT)
    )
    trade_saldo.attrs["unit"] = supply.attrs["unit"]
    trade_saldo = rename_aggregate(trade_saldo, trade_saldo.attrs["name"])

    exporter = Exporter(
        statistics=[supply, demand, trade_saldo],
        view_config=config["view"],
    )

    # view specific settings
    exporter.defaults.excel.chart = None  # charts bloat the xlsx file
    chart_class = getattr(plots, config["view"]["chart"])
    exporter.defaults.plotly.chart = chart_class

    exporter.defaults.plotly.plotby = [DM.YEAR, DM.LOCATION]
    exporter.defaults.plotly.pivot_index = [
        DM.YEAR,
        DM.LOCATION,
        DM.CARRIER,
    ]
    exporter.defaults.plotly.xaxis_title = ""

    exporter.export(result_path, config["global"]["subdir"])


def simple_optimal_capacity(
    networks: dict, config: dict, result_path: str | Path, kind: str = None
) -> None:
    """Export optimal capacities for production or demand or both."""
    (
        bus_carrier,
        transmission_comps,
        transmission_carrier,
        storage_carrier,
    ) = _parse_view_config_items(networks, config)

    optimal_capacity = (
        collect_myopic_statistics(
            networks,
            statistic="optimal_capacity",
            bus_carrier=bus_carrier,
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        .drop(storage_carrier, level=DM.CARRIER, errors="ignore")
        .droplevel(DM.COMPONENT)
    )

    if kind == "production":
        optimal_capacity = optimal_capacity[optimal_capacity > 0]
    elif kind == "demand":
        optimal_capacity = optimal_capacity[optimal_capacity < 0]

    # correct units for AC capacities
    optimal_capacity.attrs["unit"] = optimal_capacity.attrs["unit"].replace("MWh", "MW")

    exporter = Exporter(
        statistics=[optimal_capacity],
        view_config=config["view"],
    )

    # view specific constant settings
    chart_class = getattr(plots, config["view"]["chart"])
    exporter.defaults.plotly.chart = chart_class

    exporter.export(result_path, config["global"]["subdir"])


def _parse_view_config_items(networks: dict, config: dict) -> tuple:
    bus_carrier = config["view"]["bus_carrier"]
    transmission_techs = get_transmission_techs(networks, bus_carrier)
    transmission_comps = [comp for comp, carr in transmission_techs]
    transmission_carrier = [carr for comp, carr in transmission_techs]
    storage_carrier = get_storage_carriers(networks) + config["view"].get(
        "storage_links", []
    )
    return (
        bus_carrier,
        transmission_comps,
        transmission_carrier,
        storage_carrier,
    )
