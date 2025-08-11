# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Plot sankey diagrams."""

from pathlib import Path

import pandas as pd

from evals import plots as plots
from evals.constants import DataModel as DM
from evals.constants import Group, TradeTypes
from evals.fileio import Exporter
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    drop_from_multtindex_by_regex,
    filter_by,
    insert_index_level,
    rename_aggregate,
)
from evals.views.common import _parse_view_config_items

IDX = ["year", "component", "location", "carrier"]
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_rows", 500)


def _process_single_input_link(
    supply: pd.Series,
    demand: pd.Series,
    bc_in: str,
):
    balance = supply.groupby(IDX).sum() + demand.groupby(IDX).sum()
    losses = balance[balance < 0]
    surplus = balance[balance > 0]
    # losses = insert_index_level(losses, bc_in, "bus_carrier", pos=4)
    losses = insert_index_level(losses, f"{bc_in} losses", "bus_carrier", pos=4).mul(-1)
    surplus = insert_index_level(surplus, "ambient heat", "bus_carrier", pos=4)
    # if not losses.empty:
    #     # need to rename the carrier to avoid mixing with supply
    #     carrier = losses.index.unique("carrier").item()
    #     losses = rename_aggregate(losses, f"{carrier} losses")
    return pd.concat([losses, surplus])


def collect_imbalances(supply, demand):
    bc_in = demand.index.unique("bus_carrier")
    # bc_out = supply.index.unique("bus_carrier")

    if len(bc_in) > 1:
        to_concat = []
        for _bc in bc_in:
            demand_bc = filter_by(demand, bus_carrier=_bc)
            demand_share = demand_bc.sum() / demand.sum()
            supply_bc = supply * demand_share
            to_concat.append(_process_single_input_link(supply_bc, demand_bc, _bc))
            mapper = {
                v: f"{v} from {_bc}" for v in supply_bc.index.unique("bus_carrier")
            }
            to_concat.append(rename_aggregate(supply_bc, mapper, level="bus_carrier"))
        return pd.concat(to_concat)
    return _process_single_input_link(supply, demand, bc_in.item())


def get_supply(networks, transmission_comps, transmission_carrier):
    supply = (
        collect_myopic_statistics(
            networks,
            statistic="supply",
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        .pipe(
            drop_from_multtindex_by_regex, "co2|process emissions", level="bus_carrier"
        )
        .pipe(
            rename_aggregate,
            {
                "hydro": "hydro supply",
                "PHS": "PHS supply",
                "H2 Store": "H2 Store supply",
                "gas": "gas supply",
            },
        )
    )
    supply.attrs["unit"] = "MWh"

    return supply


def get_demand(networks, transmission_comps, transmission_carrier, unit):
    demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            aggregate_components=None,
        )
        .pipe(
            filter_by,
            component=transmission_comps,
            carrier=transmission_carrier,
            exclude=True,
        )
        .pipe(
            drop_from_multtindex_by_regex, "co2|process emissions", level="bus_carrier"
        )
        .pipe(
            rename_aggregate,
            {
                "hydro": "hydro demand",
                "PHS": "PHS demand",
                "H2 Store": "H2 Store demand",
                "gas": "gas demand",
            },
        )
        .mul(-1)
    )
    demand.attrs["unit"] = unit

    return demand


def net_distribution_grid_losses(supply, demand):
    grid_losses = (
        filter_by(supply, carrier="electricity distribution grid")
        .groupby(IDX)
        .sum()
        .add(
            filter_by(demand, carrier="electricity distribution grid")
            .groupby(IDX)
            .sum()
        )
    )
    grid_losses = insert_index_level(grid_losses, "losses", "bus_carrier", pos=4)
    supply.drop("electricity distribution grid", level="carrier", inplace=True)
    demand.drop("electricity distribution grid", level="carrier", inplace=True)

    return grid_losses


def get_trade_statistics(networks, transmission_comps, transmission_carrier, unit):
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
                aggregate_components=None,
            )
            # the trade statistic wrongly finds transmission between EU -> country buses.
            # Those are dropped by the filter_by statement.
            .pipe(
                filter_by,
                component=transmission_comps,
                carrier=transmission_carrier,
            )
            .pipe(drop_from_multtindex_by_regex, "co2", level="bus_carrier")
            .pipe(rename_aggregate, alias)
            .abs()
        )
        trade.attrs["unit"] = unit
        trade_statistics.append(trade)

    return trade_statistics


def get_link_losses(supply, demand):
    link_losses = []
    link_supply_carrier = filter_by(supply, component="Link").index.unique("carrier")
    link_demand_carrier = filter_by(demand, component="Link").index.unique("carrier")
    link_carrier = link_supply_carrier.union(link_demand_carrier)
    for carrier in link_carrier:
        link_supply = filter_by(supply, carrier=carrier, component="Link")
        link_demand = filter_by(demand, carrier=carrier, component="Link")
        # balance = link_supply.droplevel("bus_carrier").add(link_demand.droplevel("bus_carrier")).groupby(["year", "component", "location", "carrier"]).sum()
        if link_supply.empty or link_demand.empty:
            print(f"Skipping carrier '{carrier}' due to empty supply or demand.")
            continue
        link_losses.append(collect_imbalances(link_supply, link_demand))
        # demand.drop(link_demand.index, inplace=True)

    return link_losses


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

    supply = get_supply(networks, transmission_comps, transmission_carrier)
    demand = get_demand(
        networks, transmission_comps, transmission_carrier, unit=supply.attrs["unit"]
    )

    # todo:
    #  - calculate regional oil import from regional oil demand
    #  - calculate regional NH3 Load from regional NH3 production
    #  - calculate StorageUnit losses
    #  - harmonize V2G amounts
    grid_losses = net_distribution_grid_losses(supply, demand)
    trade_statistics = get_trade_statistics(
        networks, transmission_comps, transmission_carrier, unit=supply.attrs["unit"]
    )
    link_losses = get_link_losses(supply, demand)
    # storage_systems = (
    #     "rural water tanks",
    #     "urban central water tanks",
    #     "urban decentral water tanks",
    #     "urban central water pits",
    # )
    #
    # for storage_system in storage_systems:
    #     charger = f"{storage_system} charger"
    #     charger_losses = (
    #         filter_by(supply, carrier=charger)
    #         .droplevel("bus_carrier")
    #         .add(filter_by(demand, carrier=charger).droplevel("bus_carrier"))
    #     )
    #     assert charger_losses.abs().le(1.5).all(), (
    #         f"Charger Losses detected for carrier: {storage_system}"
    #     )
    #     supply.drop(charger, level="carrier", inplace=True)
    #     demand.drop(charger, level="carrier", inplace=True)
    #     discharger = f"{storage_system} discharger"
    #     discharger_losses = (
    #         filter_by(supply, carrier=discharger)
    #         .droplevel("bus_carrier")
    #         .add(filter_by(demand, carrier=discharger).droplevel("bus_carrier"))
    #     )
    #     assert discharger_losses.abs().le(1.5).all(), (
    #         f"Storage system imbalances detected for carrier: {storage_system}"
    #     )
    #     supply.drop(discharger, level="carrier", inplace=True)
    #     demand.drop(discharger, level="carrier", inplace=True)

    for_industry_losses = []
    # for_industry_carrier = (
    #     "coal for industry",
    #     "gas for industry",
    #     "gas for industry CC",
    #     "naphtha for industry",
    #     "solid biomass for industry",
    #     "solid biomass for industry CC",
    #     "low-temperature heat for industry",
    #     "agriculture machinery oil",
    # )
    # for industry_carrier in for_industry_carrier:
    #     industry_supply = filter_by(supply, carrier=industry_carrier, component="Link")
    #     industry_demand = filter_by(demand, carrier=industry_carrier, component="Link")
    #     if industry_supply.empty and industry_demand.empty:
    #         continue
    #     demand_bus_carrier = industry_demand.index.unique("bus_carrier").item()
    #     balance = industry_supply.droplevel("bus_carrier").add(
    #         industry_demand.droplevel("bus_carrier")
    #     )
    #     if balance.le(0).all():
    #         losses = insert_index_level(
    #             balance, demand_bus_carrier, "bus_carrier", pos=4
    #         )
    #         for_industry_losses.append(losses)
    #     elif balance.abs().gt(1e-3).any():
    #         raise ValueError(
    #             f"Unexpected carrier '{industry_carrier}' supplies energy to Load bus."
    #         )
    #     supply.drop(industry_supply.index, inplace=True)
    #     demand.drop(industry_demand.index, inplace=True)

    exporter = Exporter(
        statistics=[supply, demand, grid_losses]
        + trade_statistics
        + for_industry_losses
        + link_losses,
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
