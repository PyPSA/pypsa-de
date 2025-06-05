"""Module for hydrogen nodal balances."""

from pathlib import Path

from evals.plotting import simple_bus_balance


def view_balance_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the Hydrogen balance.

    Returns
    -------
    :

    Notes
    -----
    See eval module docstring for parameter description.
    """
    simple_bus_balance(networks, config, result_path)

    # supply = collect_myopic_statistics(
    #     networks,
    #     statistic="supply",
    #     bus_carrier=BusCarrier.H2,
    # )
    # pipelines = supply.filter(like="pipeline", axis=0).index.unique(DataModel.CARRIER)
    # supply = supply.drop(pipelines, level=DataModel.CARRIER).pipe(
    #     rename_aggregate, {"H2 Store": Group.storage_out}
    # )
    #
    # demand = (
    #     collect_myopic_statistics(
    #         networks,
    #         statistic="withdrawal",
    #         bus_carrier=BusCarrier.H2,
    #     )
    #     .mul(-1)
    #     .drop(pipelines, level=DataModel.CARRIER)
    #     .pipe(rename_aggregate, {"H2 Store": Group.storage_in})
    # )
    #
    # trade_statistics = []
    # for scope, direction, alias in [
    #     (TradeTypes.FOREIGN, "import", Group.import_foreign),
    #     (TradeTypes.FOREIGN, "export", Group.export_foreign),
    #     (TradeTypes.DOMESTIC, "import", Group.import_domestic),
    #     (TradeTypes.DOMESTIC, "export", Group.export_domestic),
    # ]:
    #     trade = (
    #         collect_myopic_statistics(
    #             networks,
    #             statistic="trade_energy",
    #             scope=scope,
    #             direction=direction,
    #             bus_carrier=BusCarrier.H2,
    #         )
    #         .filter(like="pipeline", axis=0)
    #         .pipe(rename_aggregate, alias)
    #     )
    #     trade.attrs["unit"] = "MWh_LHV"
    #     trade_statistics.append(trade)
    #
    # metric = Exporter(
    #     statistics=[supply, demand] + trade_statistics,
    #     view_config=config["view"],
    # )
    #
    # metric.defaults.plotly.chart = ESMBarChart
    # metric.defaults.plotly.pattern = dict.fromkeys(
    #     [
    #         Group.export_foreign,
    #         Group.import_foreign,
    #         Group.import_domestic,
    #         Group.export_domestic,
    #     ],
    #     "/",
    # )
    #
    # metric.export(result_path, subdir)
