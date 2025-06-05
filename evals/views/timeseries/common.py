"""Common functions for timeseries views."""

from pathlib import Path

import evals.plots as plots
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


def simple_timeseries(
    networks: dict,
    config: dict,
    result_path: str | Path,
    subdir: str | Path = "evaluation",
) -> None:
    """Export simple time series views."""
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
    # metric.defaults.plotly.legend_header = "Production/Demand"
    exporter.defaults.plotly.pattern = dict.fromkeys(
        [Group.export_net, Group.import_net, Group.import_global], "/"
    )

    exporter.export(result_path, subdir)
