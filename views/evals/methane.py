"""Module for methane evaluations."""

from logging import getLogger
from pathlib import Path

from esmtools.constants import (
    TITLE_SUFFIX,
    BusCarrier,
    Carrier,
    DataModel,
    Group,
    TradeTypes,
)
from esmtools.metric import Metric
from esmtools.plots.barchart import ESMBarChart
from esmtools.plots.timeseries import ESMTimeSeriesChart
from esmtools.statistic import collect_myopic_statistics
from esmtools.utils import (
    filter_by,
    make_evaluation_result_directories,
    rename_aggregate,
)


def eval_ch4_balance(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate Methane balances."""
    ch4_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        bus_carrier=BusCarrier.CH4,
    ).drop(["gas pipeline", ""], level=DataModel.CARRIER, errors="ignore")

    # fixme: why is CNG long/short in production?
    cng_production = filter_by(
        ch4_production, carrier=[Carrier.cng_long, Carrier.cng_short]
    ).sum()
    if cng_production > 10.0:  # MWh
        logger = getLogger(__name__)
        logger.warning(
            f"CNG long/short technologies supply large values to "
            f"the gas bus. The total amount of energy supplied to "
            f"the gas bus is {cng_production:.2f} MWh"
        )

    # non-EU imports are generator components connected to
    # CH4 buses in defined regions
    ch4_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.CH4,
    )
    ch4_generation = rename_aggregate(ch4_generation, Group.global_market)

    ch4_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            comps=["Link", "Load"],
            bus_carrier=BusCarrier.CH4,
        )
        .drop(["gas pipeline"], level=DataModel.CARRIER)
        .mul(-1)
    )

    # import / export
    ch4_import_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="import",
        bus_carrier=BusCarrier.CH4,
    )
    ch4_import_foreign = rename_aggregate(ch4_import_foreign, Group.import_european)

    ch4_export_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="export",
        bus_carrier=BusCarrier.CH4,
    )
    ch4_export_foreign = rename_aggregate(ch4_export_foreign, Group.export_european)

    ch4_import_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        bus_carrier=BusCarrier.CH4,
    )
    ch4_import_domestic = rename_aggregate(ch4_import_domestic, Group.import_domestic)

    ch4_export_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="export",
        bus_carrier=BusCarrier.CH4,
    )
    ch4_export_domestic = rename_aggregate(ch4_export_domestic, Group.export_domestic)

    metric = Metric(
        metric_name="Methane Balance",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[
            ch4_production,
            ch4_generation,
            ch4_import_foreign,
            ch4_demand,
            ch4_export_foreign,
            ch4_import_domestic,
            ch4_export_domestic,
        ],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX

    metric.cfg.mapping = "capacity"
    metric.cfg.excel.chart_title = title

    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "methane_balance_{location}"
    metric.cfg.plotly.cutoff = 0.1  # TWh
    metric.cfg.plotly.pattern = dict.fromkeys(
        [
            Group.export_european,
            Group.import_european,
            Group.global_market,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )
    metric.cfg.plotly.category_orders = (
        # order production from outside to zero
        Group.ch4_bio_processing,
        Group.import_european,
        Group.global_market,
        # --- zero ---
        # demand ordered from zero to outside:
        Group.industry,
        Group.chp_electricity,
        Group.heat_decentral,
        Group.heat_district,
        Group.ocgt_electricity,
        Group.transport,
        Group.smr,
        Group.export_european,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_ch4_balance_ts(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate Methane balances time series."""
    ch4_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        bus_carrier=BusCarrier.CH4,
        aggregate_time=False,
    ).drop(["gas pipeline", ""], level=DataModel.CARRIER, errors="ignore")

    cng_production = (
        filter_by(ch4_production, carrier=[Carrier.cng_long, Carrier.cng_short])
        .sum()
        .sum()
    )
    if cng_production > 10.0:  # MWh
        logger = getLogger(__name__)
        logger.warning(
            f"CNG long/short technologies supply large values to "
            f"the gas bus. The total amount of energy supplied to "
            f"the gas bus is {cng_production:.2f} MWh"
        )

    global_import = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.CH4,
        aggregate_time=False,
    )
    global_import = rename_aggregate(global_import, Group.import_global)
    # note, that me must differentiate between imports from trade via
    # pipelines and global import from non-EU countries by the Group
    # name. The EU aggregation in esmtools.utils.aggregate_eu will
    # drop import/export names to remove all zero traces.

    ch4_storage_out = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Store",
        bus_carrier=BusCarrier.CH4,
        aggregate_time=False,
    )
    ch4_storage_out = rename_aggregate(ch4_storage_out, Group.storage_out)

    ch4_demand = (
        collect_myopic_statistics(
            networks,
            statistic="withdrawal",
            comps=["Link", "Load"],
            bus_carrier=BusCarrier.CH4,
            aggregate_time=False,
        )
        .drop(["gas pipeline"], level=DataModel.CARRIER, errors="ignore")
        .mul(-1.0)
    )

    ch4_storage_in = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Store",
        bus_carrier=BusCarrier.CH4,
        aggregate_time=False,
    ).mul(-1.0)
    ch4_storage_in = rename_aggregate(ch4_storage_in, Group.storage_in)

    # import / export
    trade_saldo = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=(TradeTypes.FOREIGN, TradeTypes.DOMESTIC),
        direction="saldo",
        bus_carrier=BusCarrier.CH4,
        aggregate_time=False,
    )
    trade_saldo = rename_aggregate(trade_saldo, trade_saldo.attrs["name"])

    metric = Metric(
        metric_name="Methane Production and Demand",
        is_unit="MWh",
        to_unit="MWh",
        statistics=[
            ch4_production,
            global_import,
            ch4_storage_out,
            ch4_demand,
            ch4_storage_in,
            trade_saldo,
        ],
    )

    metric.cfg.mapping = "capacity"
    metric.cfg.excel.chart = None

    metric.cfg.plotly.title = metric.df.attrs["name"] + " {location} in {year}"
    metric.cfg.plotly.file_name_template = "methane_prod_dem_time_{year}_{location}"
    metric.cfg.plotly.chart = ESMTimeSeriesChart
    metric.cfg.plotly.plotby = [DataModel.YEAR, DataModel.LOCATION]
    metric.cfg.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]
    metric.cfg.plotly.xaxis_title = ""
    metric.cfg.plotly.legend_header = "Production/Demand"
    # The Toolbox rounds values to second digit before clipping.
    metric.cfg.plotly.cutoff = 0.001  # MWh
    metric.cfg.plotly.pattern = dict.fromkeys(
        [Group.export_net, Group.import_net, Group.import_global], "/"
    )

    metric.cfg.plotly.category_orders = (
        # list production entries from outside to zero
        Group.import_net,
        Group.import_global,
        Group.storage_out,
        Group.smr,
        Group.electrolysis,
        Group.ch4_bio_processing,
        # list demand entries from outside to zero
        Group.export_net,
        Group.storage_in,
        Group.chp_electricity,
        Group.ocgt_electricity,
        Group.heat_district,
        Group.heat_decentral,
        Group.smr,
        Group.transport,
        Group.industry,
    )
    # metric.cfg.plotly.footnotes = (
    #     "<br><b>Miscellaneous</b> includes residential and "
    #     "services rural H2-powered PEMFC.",
    #     "",
    # )

    output_path = make_evaluation_result_directories(result_path, subdir)
    # metric.export_excel(output_path)
    # metric.export_csv(output_path)
    metric.export_plotly(output_path)
