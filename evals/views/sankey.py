# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Plot sankey diagrams."""

from pathlib import Path

from evals import plots as plots
from evals.constants import DataModel as DM
from evals.constants import Group, TradeTypes
from evals.fileio import Exporter
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    filter_by,
    rename_aggregate,
)
from evals.views.common import _parse_view_config_items


def view_sankey(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """
    Evaluate the carbon balance.

    Returns
    -------
    :
    """
    (
        _,
        transmission_comps,
        transmission_carrier,
        storage_links,
    ) = _parse_view_config_items(networks, config)

    supply = (
        collect_myopic_statistics(
            networks,
            statistic="supply",
            # bus_carrier=bus_carrier,
            aggregate_components=None,
        ).pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        # .pipe(rename_aggregate, dict.fromkeys(storage_links, Group.storage_out))
        # .droplevel(DM.COMPONENT)
    )
    supply.attrs["unit"] = config["view"]["unit"]

    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            # bus_carrier=bus_carrier,
            aggregate_components=None,
        ).pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        # .pipe(rename_aggregate, dict.fromkeys(storage_links, Group.storage_in))
        # .mul(-1)
        # .droplevel(DM.COMPONENT)
    )
    demand.attrs["unit"] = config["view"]["unit"]

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
                # bus_carrier=bus_carrier,
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
            .abs()
            # .droplevel(DM.COMPONENT)
        )
        trade.attrs["unit"] = config["view"]["unit"]
        trade_statistics.append(trade)

    exporter = Exporter(
        statistics=[supply, demand] + trade_statistics,
        view_config=config["view"],
    )

    # df = exporter.df
    # print(df)

    chart_class = getattr(plots, config["view"]["chart"])
    exporter.defaults.plotly.chart = chart_class

    exporter.defaults.plotly.plotby = [DM.YEAR, DM.LOCATION]
    exporter.defaults.plotly.pivot_index = [
        DM.COMPONENT,
        DM.YEAR,
        DM.LOCATION,
        DM.CARRIER,
        DM.BUS_CARRIER,
    ]
    # exporter.defaults.plotly.xaxis_title = ""

    exporter.export(result_path, config["global"]["subdir"])
