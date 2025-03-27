# -*- coding: utf-8 -*-
"""Module for electricity evaluations."""

from functools import partial
from pathlib import Path

from constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel, Group
from metric import Metric
from plots.barchart import ESMBarChart
from statistic import (
    collect_myopic_statistics,
    get_location_and_carrier_and_bus_carrier,
)
from utils import make_evaluation_result_directories, rename_aggregate


def eval_electricity_capacities(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the optimal capacity for AC technologies.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.
    """
    ac_generation_and_storage = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        comps=("Generator", "Store", "StorageUnit"),
        bus_carrier=BusCarrier.AC,
    ).drop("value of lost load", level=DataModel.CARRIER, errors="ignore")

    transmission_or_storage_links = ["", "DC", Carrier.v2g]
    ac_production = (
        collect_myopic_statistics(
            networks,
            statistic="optimal_capacity",
            comps="Link",
            bus_carrier=BusCarrier.AC,
        )
        .clip(lower=0)
        .drop(transmission_or_storage_links, level=DataModel.CARRIER, errors="ignore")
    )

    metric = Metric(
        metric_name="Optimal Capacity Electricity",
        is_unit="MW",
        to_unit="GW",
        statistics=[ac_generation_and_storage, ac_production],
    )

    metric.defaults.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "eKapas_{location}"
    # metric.cfg.plotly.cutoff = 0.0
    metric.defaults.plotly.category_orders = (
        Group.battery_storage,
        Group.nuclear_power,
        Group.pp_thermal,
        Group.phs,
        Group.reservoir,
        Group.ror,
        Group.wind,
        Group.pv,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_fuel_production_capacities(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the fuel production capacities per country.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See eval docstring for parameter description.
    """
    pipelines_or_biogas = [
        Carrier.gas_pipepline,
        Carrier.h2_pipeline,
        Carrier.h2_pipeline_retro,
        Carrier.biogas_to_ch4,  # biogas from energy amounts in a separate statistic
    ]
    # gas generation is excluded from this evaluation.
    gas_production = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        comps="Link",
        bus_carrier=[BusCarrier.CH4, BusCarrier.H2],
    ).clip(
        lower=0
    )  # supply only
    gas_production = gas_production.drop(
        pipelines_or_biogas, level=DataModel.CARRIER, errors="ignore"
    )

    oil_production = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        comps="Link",
        at_port=True,
        bus_carrier=[BusCarrier.FT_1, BusCarrier.FT_2],
        # use the port 0 (H2) locations for all branches
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="0"),
    ).clip(lower=0)

    # rename the FT1/2 bus_carrier to Fischer-Tropsch.
    # The problem is, that the rename_aggregate function will
    # only hit the carrier level, but the respective bus_carriers stay
    # untouched.
    # This only has an impact if the cut_off values is high
    # enough to exclude separate carrier - bus_carrier combinations,
    # but not the combined carrier sums.
    # Another way to explain, is that the cut_off is applied to
    # the carrier-bus_carrier items separately, whereas the old Toolbox
    # applies the cutoff to carrier level sums and this leads to
    # differences between the old and the new evaluation library.
    oil_production = rename_aggregate(
        oil_production,
        dict.fromkeys([BusCarrier.FT_1, BusCarrier.FT_2], Group.ft),
        level=DataModel.BUS_CARRIER,
    )

    biogas_to_gas = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        bus_carrier=BusCarrier.CH4,
        carrier=[Carrier.biogas_to_ch4],
    )
    hours_productive = 8000  # 760h / year maintenance assumed
    biogas_production = biogas_to_gas.div(hours_productive)  # MW

    metric = Metric(
        metric_name="Fuel Production Capacities",
        is_unit="MW",
        to_unit="GW",
        statistics=[
            gas_production,
            oil_production,
            biogas_production,
        ],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = "capacity"
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.file_name_template = "fKapas_{location}"
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.cutoff = 0.0001  # GW
    metric.defaults.plotly.legend_header = "Types of Production"
    metric.defaults.plotly.category_orders = (
        Group.electrolysis,
        Group.smr,
        Group.methanation,
        Group.ch4_bio_processing,
        Group.ft,
    )
    metric.defaults.plotly.footnotes = (
        f"Values in thermic units (Output)<br>Bio Methane Processing approximated "
        f"as total quantity divided by {hours_productive} hours (since only base "
        f"load, capacity is linked to quantity).",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
