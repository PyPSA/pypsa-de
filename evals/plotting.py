import evals.plots as plots
from evals.constants import Group, TradeTypes
from evals.fileio import Exporter
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    filter_by,
    get_storage_carriers,
    get_transmission_techs,
    rename_aggregate,
)


def plot_bus_balance(
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

    # todo: read csvs and compare with calculated results to test the code
    # supply = read_pypsa_csv(result_path, "nodal_supply", index_cols=4) #.drop(transmission_carrier, level=DM.CARRIER)
    # demand = read_pypsa_csv(result_path, "nodal_withdrawal", index_cols=4)# .drop(transmission_carrier, level=DM.CARRIER)
    # print(filter_by(supply, year="2050", bus_carrier="AC").sum() / 1e6)
    # print(filter_by(demand, year="2050", bus_carrier="AC").sum() / 1e6)
    #
    # balance = read_pypsa_csv(result_path, "nodal_energy_balance", index_cols=4)
    # print(filter_by(balance, year="2050", bus_carrier="AC").drop(transmission_carrier, level=DM.CARRIER).sum() / 1e6)

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
        # .drop(transmission_carrier, level=DM.CARRIER)
        # .pipe(rename_aggregate, dict.fromkeys(transmission_carrier, "Import"))
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_out))
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
        # .drop(transmission_techs, level=DM.CARRIER)
        # .pipe(rename_aggregate, dict.fromkeys(transmission_carrier, "Export"))
        .pipe(rename_aggregate, dict.fromkeys(storage_carrier, Group.storage_in))
        .mul(-1)
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

    # exporter.defaults.plotly.pattern = dict.fromkeys(
    #     [
    #         Group.import_foreign,
    #         Group.export_foreign,
    #         Group.import_domestic,
    #         Group.export_domestic,
    #     ],
    #     "/",
    # )

    exporter.export(result_path, config["global"]["subdir"])


if __name__ == "__main__":
    from pathlib import Path

    from evals.fileio import read_networks, read_views_config
    from evals.views import (
        view_balance_electricity,
        view_balance_hydrogen,
        view_balance_methane,
    )

    _result_path = Path("results/v2025.02/KN2045_Mix")
    _networks = read_networks(_result_path)

    # evaluations = [
    #     (
    #         "Electricity",
    #         ["AC", "low voltage", "EV battery", "DC"],
    #         ["BEV charger", "V2G"],
    #     ),
    #     ("Hydrogen", ["H2"], []),
    #     (
    #         "Methane",
    #         [
    #             "gas",
    #         ],
    #     ),
    # ]

    for func in [view_balance_electricity, view_balance_methane, view_balance_hydrogen][
        :1
    ]:
        _config = read_views_config(func)
        plot_bus_balance(_networks, _config, _result_path)

    # for _name, _bus_carrier, _additional_storage_techs in evaluations:
    #     plot_bus_balance(_networks, _config, _result_path)
