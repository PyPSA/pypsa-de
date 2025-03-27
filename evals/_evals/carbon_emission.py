# -*- coding: utf-8 -*-
"""Module for CO2 Emission evaluations."""

from functools import partial
from pathlib import Path

import pandas as pd
import pypsa.statistics
from constants import TITLE_SUFFIX, Carrier, DataModel, Group
from fileio import (
    prepare_co2_emissions,
    prepare_costs,
    prepare_industry_demand,
    prepare_nodal_energy,
)
from metric import Metric
from plots.barchart import ESMBarChart
from statistic import (
    collect_myopic_statistics,
    get_location_and_carrier_and_bus_carrier,
)
from utils import (
    filter_by,
    make_evaluation_result_directories,
    rename_aggregate,
)


def eval_co2(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate Carbon Dioxide balances per country.

    Returns
    -------
    :
        Exports one BarChart per country and associated tables and JSON.

    Notes
    -----
    See pacakge docstring for parameter description.
    """
    costs = prepare_costs(result_path, n_years=25)
    co2_intensity = filter_by(costs, carrier=["oil", "coal", "gas"])
    co2_intensity = co2_intensity["CO2 intensity"].unstack("year")  # t CO2 / MWh

    nodal_energy = prepare_nodal_energy(result_path)
    oil_for_transport = filter_by(
        nodal_energy,
        carrier=[
            "Fischer-Tropsch road freight",
            "Fischer-Tropsch rail",
            "Fischer-Tropsch domestic navigation",
            "Fischer-Tropsch domestic aviation",
        ],
    )
    transport_co2 = oil_for_transport.mul(co2_intensity.loc[["oil"], :]).T.squeeze()

    industry_demand = prepare_industry_demand(result_path, networks)

    ft_for_industry = filter_by(industry_demand, carrier="Fischer-Tropsch")
    ft_co2 = (ft_for_industry * co2_intensity.loc[["oil"], :]).T.squeeze()
    ft_co2 = rename_aggregate(ft_co2, "Fischer-Tropsch industry")

    coal_for_industry = filter_by(industry_demand, carrier=["hard coal"])
    coal_co2 = coal_for_industry.mul(co2_intensity.loc[["coal"], :]).T.squeeze()
    coal_co2 = rename_aggregate(coal_co2, "hard coal industry")

    process_emissions_eu = collect_myopic_statistics(
        networks,
        "energy_balance",
        comps="Link",
        bus_carrier="co2",
        carrier=["process emissions", "process emissions CC"],
    )
    process_emissions_country = filter_by(
        industry_demand, carrier=["process emission", "process emission from feedstock"]
    )
    process_emissions = _distribute_eu_values_by_country_share(
        process_emissions_eu, process_emissions_country
    )

    # beware: the old toolbox accessed n.df("Load").p_set. We assume
    # that p_set equals the max value from the respective n.pnl("Load")
    # time series. If this assumption breaks, differences will arise.
    gas_loads_p_set = collect_myopic_statistics(
        networks,
        "withdrawal",
        comps="Load",
        bus_carrier="gas",
        carrier=["gas domestic navigation", "gas road freight"],
        aggregate_time="max",
    )
    gas_emissions = gas_loads_p_set.mul(co2_intensity.loc[["gas"], :]).T.squeeze()

    # the 'max' time aggregation skips multiplying with snapshot weights
    # fixme: question - Why max load for all hours?
    first_network = next(iter(networks.values()))
    gas_emissions *= pypsa.statistics.get_weightings(first_network, "Load").sum()

    co2_from_links = collect_myopic_statistics(
        networks,
        "energy_balance",
        comps="Link",
        bus_carrier="co2",
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="1"),
    ).drop("EU", level=DataModel.LOCATION)

    dac = collect_myopic_statistics(
        networks,
        "energy_balance",
        comps="Link",
        bus_carrier="co2",
        carrier=[Carrier.dac],
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="2"),
    )

    ft_import_eu = collect_myopic_statistics(
        networks,
        "energy_balance",
        comps="Link",
        bus_carrier="co2",
        carrier=["Fischer-Tropsch import link 1"],
    )
    ft_import = _distribute_eu_values_by_country_share(ft_import_eu, transport_co2)

    co2_vent_eu = collect_myopic_statistics(
        networks,
        "energy_balance",
        comps="Link",
        bus_carrier="co2 store",
        carrier=["co2 vent"],
    )

    # historic CO2 emissions from data file lack associated years
    co2_country_share = pd.DataFrame.from_dict(
        {
            year: prepare_co2_emissions(result_path, "THI")
            for year in co2_vent_eu.index.unique("year")
        }
    )
    co2_country_share.columns.name = "year"
    co2_vent = (
        co2_vent_eu.droplevel("location")  # location = 'EU'
        .mul(co2_country_share)
        .T.stack("location")  # location becomes country code
        .reorder_levels(DataModel.YEAR_IDX_NAMES)
    )

    metric = Metric(
        "CO2 Emissions",
        is_unit="t",
        to_unit="Mt",
        statistics=[
            transport_co2,
            ft_co2,
            coal_co2,
            ft_import,
            process_emissions,
            gas_emissions,
            dac,
            co2_from_links,
            co2_vent,
        ],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = "co2"
    metric.defaults.excel.chart_title = title

    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.pivot_index = DataModel.YEAR_IDX_NAMES[:-1]
    metric.defaults.plotly.legend_header = "Types of Emitters"
    metric.defaults.plotly.unit = f"{metric.df.attrs['unit']} CO2"
    metric.defaults.plotly.cutoff = 0.1
    metric.defaults.plotly.category_orders = (
        Group.transport,
        Group.smr,
        Group.smr_cc,
        Group.industry,
    )
    metric.defaults.plotly.file_name_template = "co2_emissions_{location}"

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_plotly(output_path)
    metric.export_excel(output_path)
    metric.export_csv(output_path)


def _distribute_eu_values_by_country_share(
    statistic_eu: pd.Series, country_base: pd.Series
) -> pd.Series:
    """Spread CO2 emission from EU level to countries.

    Parameters
    ----------
    statistic_eu
        The CO2 emissions for the EU region.
    country_base
        The base values for a country.

    Returns
    -------
    :
        CO2 emissions by country.
    """
    statistic_eu = statistic_eu.droplevel(DataModel.LOCATION)  # location must be 'EU'
    totals_country = country_base.groupby([DataModel.YEAR, DataModel.LOCATION]).sum()
    totals_eu = country_base.groupby(DataModel.YEAR).sum()
    share = totals_country / totals_eu
    statistic_country = statistic_eu * share

    return statistic_country.reorder_levels(DataModel.YEAR_IDX_NAMES)
