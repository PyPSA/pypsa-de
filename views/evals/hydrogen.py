"""Module for hydrogen evaluations."""

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
from esmtools.utils import make_evaluation_result_directories, rename_aggregate


def eval_h2_balance_ts(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the Hydrogen balance time series."""
    h2_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        bus_carrier=BusCarrier.H2,
        carrier=[Carrier.h2_electrolysis, Carrier.smr, Carrier.smr_cc],
        aggregate_time=False,
    )

    h2_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    )
    h2_generation = rename_aggregate(h2_generation, Group.import_global)

    h2_storage_out = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps=["Store", "StorageUnit"],
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    )
    h2_storage_out = rename_aggregate(
        h2_storage_out, Group.storage_out, level=DataModel.CARRIER
    )

    h2_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps=["Link", "Load"],
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    ).mul(-1.0)
    h2_demand = h2_demand.drop(
        [Carrier.h2_pipeline, Carrier.h2_pipeline_retro], level=DataModel.CARRIER
    )

    h2_storage_in = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps=["Store", "StorageUnit"],
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    ).mul(-1.0)
    h2_storage_in = rename_aggregate(
        h2_storage_in, Group.storage_in, level=DataModel.CARRIER
    )

    # import / export
    trade_saldo = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=(TradeTypes.FOREIGN, TradeTypes.DOMESTIC),
        direction="saldo",
        bus_carrier=BusCarrier.H2,
        aggregate_time=False,
    )
    trade_saldo = rename_aggregate(trade_saldo, trade_saldo.attrs["name"])

    metric = Metric(
        metric_name="Hydrogen Production and Demand",
        is_unit="MWh",
        to_unit="MWh",
        statistics=[
            h2_production,
            h2_generation,
            h2_storage_out,
            h2_demand,
            h2_storage_in,
            trade_saldo,
        ],
    )
    metric.cfg.mapping = "capacity"
    metric.cfg.excel.chart = None

    metric.cfg.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.plotly.file_name_template = "hydrogen_prod_dem_time_{year}_{location}"
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
        # list demand entries from outside to zero
        Group.export_net,
        Group.misc,
        Group.storage_in,
        Group.methanation,
        Group.ft,
        Group.fuel_cell,
        Group.transport,
        Group.industry,
    )
    metric.cfg.plotly.footnotes = (
        "<br><b>Miscellaneous</b> includes residential and "
        "services rural H2-powered PEMFC.",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    # metric.export_excel(output_path)
    # metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_h2_balance(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the Hydrogen balance."""
    h2_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        bus_carrier=BusCarrier.H2,
        carrier=[Carrier.h2_electrolysis, Carrier.smr, Carrier.smr_cc],
    )

    # non-EU imports are generator components connected to
    # H2 buses in defined regions
    h2_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.H2,
    )
    h2_generation = rename_aggregate(h2_generation, Group.import_global)

    h2_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps=["Link", "Load"],
        bus_carrier=BusCarrier.H2,
    ).mul(-1)
    h2_demand = h2_demand.drop(
        [Carrier.h2_pipeline, Carrier.h2_pipeline_retro],
        level=DataModel.CARRIER,
        errors="ignore",
    )

    h2_import_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="import",
        bus_carrier=BusCarrier.H2,
    )
    h2_import_foreign = rename_aggregate(h2_import_foreign, Group.import_european)

    h2_export_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="export",
        bus_carrier=BusCarrier.H2,
    )
    h2_export_foreign = rename_aggregate(h2_export_foreign, Group.export_european)

    h2_import_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        bus_carrier=BusCarrier.H2,
    )
    h2_import_domestic = rename_aggregate(h2_import_domestic, Group.import_domestic)

    h2_export_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="export",
        bus_carrier=BusCarrier.H2,
    )
    h2_export_domestic = rename_aggregate(h2_export_domestic, Group.export_domestic)

    metric = Metric(
        metric_name="Hydrogen Balance",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[
            h2_production,
            h2_generation,
            h2_demand,
            h2_import_foreign,
            h2_export_foreign,
            h2_import_domestic,
            h2_export_domestic,
        ],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX

    metric.cfg.mapping = "h2_balance"  # "capacity"
    metric.cfg.excel.chart_title = title

    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "hydrogen_balance_{location}"
    # fixme: Industry values for 2030 AT02 (Vienna) Hydrogen Balance
    #  (TWh) are exactly 0.1. The utils.apply_cutoff function truncates
    #  0.1 although the boundary is inclusive (0.1 should be kept). This
    #  is probably because of the imprecise float data type. The 0.099
    #  below is a hotfix for this problem. A better should be found!
    metric.cfg.plotly.cutoff = 0.0999999  # should be 0.1 TWh
    metric.cfg.plotly.pattern = dict.fromkeys(
        [
            Group.export_european,
            Group.import_european,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )
    metric.cfg.plotly.category_orders = (
        # production ordered from zero to outside
        Group.electrolysis,
        Group.smr,
        Group.import_domestic,
        Group.import_european,
        # demand ordered from zero to outside:
        Group.industry,
        Group.transport,
        Group.synth_fuels,
        Group.export_domestic,
        Group.export_european,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
