# -*- coding: utf-8 -*-
"""Module for electricity evaluations."""

import logging
from pathlib import Path

import pandas as pd
from constants import (
    TITLE_SUFFIX,
    BusCarrier,
    Carrier,
    DataModel,
    Group,
    TradeTypes,
)
from metric import Metric
from plots.barchart import ESMBarChart
from plots.timeseries import ESMTimeSeriesChart
from statistic import collect_myopic_statistics
from utils import (
    expand_to_time_series,
    filter_by,
    get_mapping,
    make_evaluation_result_directories,
    operations_override,
    rename_aggregate,
)


def eval_electricity_demand(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the electricity demand by country and year.

    Electricity demand is sum of power supply from
      - electricity loads split by sector,
      - "Link" component withdrawal at AC buses, and
      - electricity demand from Battery Electric Vehicles
        and Vehicle-2-Grid technologies

    Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See eval docstring for parameter description.
    """
    ac_load_split = (
        collect_myopic_statistics(networks, statistic="ac_load_split")
        .clip(upper=0)  # withdrawal only
        .mul(-1.0)  # as positive values
    )

    losses_and_hh_load = [Carrier.grid_losses, "industry + hh & services load"]
    # note, that industry and hh-services are still in the load split.
    # They are in 'industry' and 'domestic homes and trade'. Just the
    # combination "industry + hh & services load" is superfluous.
    ac_load_split = ac_load_split.drop(losses_and_hh_load, level=DataModel.CARRIER)

    bev_demand = collect_myopic_statistics(
        networks, statistic="bev_v2g", carrier=[Carrier.bev_demand]
    ).mul(-1.0)

    ac_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Link",
        bus_carrier=BusCarrier.AC,
    )

    transmission_or_storage = ["", "DC", "BEV charger", "battery charger"]
    ac_demand = ac_demand.drop(
        transmission_or_storage, level=DataModel.CARRIER, errors="ignore"
    )

    metric = Metric(
        metric_name="Power Demand",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[ac_load_split, bev_demand, ac_demand],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = "e_demand"
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "elec_demand_{location}"
    metric.defaults.plotly.category_orders = (
        Group.industry,
        Group.hh_and_services,
        Group.dac,
        Group.electrolysis,
        Group.heat_decentral,
        Group.heat_district,
        Group.transport,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_electricity_production(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the electricity production by country and year.

    Electricity production is sum of power supply from
      - "Generator" component,
      - Link supply at AC buses, and
      - Pump-Hydro-Storage and Run-Of-River inflows

    Writes 2 Excel files and 1 BarChart per country.
    """
    ac_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        bus_carrier=BusCarrier.AC,
    ).drop("", level="carrier", errors="ignore")

    ac_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.AC,
    )
    phs_hydro_dispatched_power = collect_myopic_statistics(
        networks,
        "phs_split",
        carrier=["hydro Dispatched Power", "PHS Dispatched Power from Inflow"],
    )

    metric = Metric(
        "Electricity Production",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[ac_production, ac_generation, phs_hydro_dispatched_power],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = "e_production"
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "elec_production_{location}"
    metric.defaults.plotly.category_orders = (
        Group.ror,
        Group.phs_inflow,
        Group.wind,
        Group.pv,
        Group.pp_thermal,
        Group.chp_biomass,
        Group.nuclear_power,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_residual_load(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the residual load for AC technologies.

    The residual load is defined as the part of the inflexible demand,
    that cannot be met by renewable generators, i.e. the difference
    of the inflexible demand minus PV, Wind and Hydro.

    It is a measure for loads, that must be met by backup
    generators with higher marginal costs, e.g. gas, oil, H2, etc.

    Writes 1 Excel file and 1 TimeSeries chart per country.
    """
    mapping = get_mapping("default", "external")
    renewables = [
        carrier
        for carrier, group in mapping.items()
        if group in [Group.pv, Group.wind, Group.ror]
    ]

    renewable_production = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        aggregate_time=False,
        bus_carrier=BusCarrier.AC,
        carrier=renewables,
    )

    ac_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Load",
        aggregate_time=False,
        bus_carrier=BusCarrier.AC,
    )

    # assuming that Plug-In Hybrid Electric Vehicles cannot load shift
    phev_demand = collect_myopic_statistics(
        networks,
        statistic="energy_balance",
        comps="Link",
        aggregate_time=False,
        bus_carrier=BusCarrier.AC,
        carrier=[Carrier.phev_long, Carrier.phev_short],
    )

    metric = Metric(
        metric_name="Residual Load",
        is_unit="MWh",
        to_unit="MWh",
        statistics=[renewable_production.abs(), ac_demand.abs(), phev_demand.abs()],
    )

    metric.defaults.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.plotly.chart = ESMTimeSeriesChart
    metric.defaults.plotly.plotby = [DataModel.LOCATION, DataModel.YEAR]
    metric.defaults.plotly.file_name_template = "residual_load_{year}_{location}"
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,  # contains only AC
    ]
    metric.defaults.plotly.categories = DataModel.BUS_CARRIER
    metric.defaults.plotly.footnotes = (
        "Residual load is calculated as the households & services demand,"
        " grid losses, demand from industry and from transport(excluding"
        " EV) minus the total production from photovoltaics, wind power"
        " and run-of-river, cut off at zero.",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_plotly(output_path)


def eval_electricity_balance_ts(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the electricity production and demand balance.

    Electricity production is sum of supply from
      - "StorageUnit", "Store", and "Generator" components,
      - Link supply at AC buses, and
      - electricity imports

    Electricity demand is the sum of withdrawal from
      - "Store", "Generator", and "Load" components,
      - p_store time series for "StorageUnit" components,
      - Link withdrawal at AC buses, and
      - electricity exports

    Loads for industry and rail are subtracted from "electricity"
    (Base Load) and added as separate carrier "Industry" and
    "Transport", respectively. Note, that this introduces possible
    imbalances if industry+rail is larger than the base load.

    The additional data series "Inflexible Demand" is calculated as
    the sum of loads that cannot be shut down. Contributing
    categories are "Industry", "Base Load", and "Transport"
    (excluding "BEV charger" carrier). They are plotted against the
    production for visual comparison with the same.

    Writes 1 TimeSeries chart per country per year.
    """
    output_path = make_evaluation_result_directories(result_path, subdir)

    ac_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        # no StorageUnit because it uses p_dispatch
        comps=("Store", "Generator"),
        bus_carrier=BusCarrier.AC,
        aggregate_time=False,
    )
    ac_generation = ac_generation.drop(
        "value of lost load", level="carrier", errors="ignore"
    )

    with operations_override(networks, "StorageUnit", "p_dispatch"):
        phs_supply = collect_myopic_statistics(
            networks,
            statistic="supply",
            comps="StorageUnit",
            bus_carrier=BusCarrier.AC,
            aggregate_time=False,
        )

    transmission_links = ["", "DC"]
    ac_production = collect_myopic_statistics(
        networks,
        statistic="energy_balance",
        comps="Link",
        bus_carrier=BusCarrier.AC,
        aggregate_time=False,
    ).drop(transmission_links, level=DataModel.CARRIER, errors="ignore")

    # Some technologies exist in both, supply and withdrawal as "Storage Out"
    # and "Storage In", respectively. The old mappings for
    # e_production_ts and e_demand_ts only differ for those carrier.
    # By renaming them explicitly here to Storage Out, we can combine
    # the mappings and avoid some redundancy.
    phs_supply = rename_aggregate(phs_supply, {Carrier.phs: Group.storage_out})
    ac_production = rename_aggregate(
        ac_production, {Carrier.battery_discharger: Group.storage_out}
    )

    with operations_override(networks, "StorageUnit", "p_store"):
        phs_demand = collect_myopic_statistics(
            networks,
            statistic="supply",  # p_store has positive values for demand
            comps="StorageUnit",
            bus_carrier=BusCarrier.AC,
            aggregate_time=False,
        ).mul(
            -1
        )  # demand is negative

    # replicate time series for industry and rail yearly demand
    rail_industry = collect_myopic_statistics(
        networks,
        statistic="ac_load_split",
        carrier=[Carrier.electricity_rail, Carrier.industry],
        drop_zero_rows=False,
    )
    rail_industry = expand_to_time_series(rail_industry, networks["2030"].snapshots)

    ac_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        # no Storage Unit because it uses p_store
        comps=("Store", "Generator", "Load"),
        bus_carrier=BusCarrier.AC,
        aggregate_time=False,
        drop_zero_rows=False,
    ).mul(
        -1.0
    )  # demand is negative
    ac_demand = ac_demand.drop("value of lost load", level="carrier", errors="ignore")
    ac_demand = _subtract_rail_and_industry_demands_from_electricity_base_load(
        ac_demand, rail_industry
    )

    trade_saldo = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=(TradeTypes.FOREIGN, TradeTypes.DOMESTIC),
        direction="saldo",
        bus_carrier=BusCarrier.AC,  # bus_carrier DC does not exist
        aggregate_time=False,
    )
    trade_saldo = rename_aggregate(trade_saldo, "trade saldo")

    to_concat = [
        ac_generation,
        ac_production,
        phs_supply,
        ac_demand,
        phs_demand,
        rail_industry,
        trade_saldo,
    ]
    statistics = pd.concat(to_concat)

    inflexible_demand = _calculate_inflexible_demand(statistics)

    # reverse sign to plot inflexible demand upwards along the supply
    metric = Metric(
        "Power Production and Demand",
        is_unit="MWh",
        to_unit="MWh",
        statistics=[statistics, -inflexible_demand],
    )

    metric.defaults.mapping = "electricity"
    metric.defaults.plotly.title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.plotly.chart = ESMTimeSeriesChart
    metric.defaults.plotly.file_name_template = "elec_prod_dem_time_{year}_{location}"
    metric.defaults.plotly.plotby = [DataModel.YEAR, DataModel.LOCATION]
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]
    metric.defaults.plotly.cutoff = 0.04999999  # None in the original Toolbox
    # the original toolbox discards values smaller than 0.05 from
    # production and demand time series:
    # df.apply(lambda x: np.nan if abs(x) < 0.05 else x)
    metric.defaults.plotly.xaxis_title = ""
    metric.defaults.plotly.pattern = dict.fromkeys(
        [Group.import_net, Group.export_net], "/"
    )
    metric.defaults.plotly.legend_header = "Production/Demand"
    metric.defaults.plotly.footnotes = (
        """
<b>Storage In/Out</b> includes: Pumped Hydro Storage, Reservoirs, Batteries
and V2G.
<br><b>Base load</b> includes: demand for households, demand for services,
grid losses. <b>Industry</b> includes: demand for industry, demand for
direct air capture.
<br><b>Thermal Powerplants</b> includes: Biomass CHPs, Hydrogen Fuell Cells,
Methane Powerplants and PEMFCs, Coal Powerplants, Lignite Powerplants, Oil
Powerplants.
<br><b>Inflexible Demand</b> includes: Base Load, Industry and Transport.""",
        "",
    )
    metric.defaults.plotly.category_orders = (
        "Inflexible Demand",
        Group.import_net,
        Group.storage_out,
        Group.pv,
        Group.wind,
        Group.pp_thermal,
        Group.ror,
        # --- zero ---
        Group.export_net,
        Group.storage_in,
        Group.electrolysis,
        Group.heat,
        Group.transport,
        Group.base_load,
        Group.industry,
    )

    metric.export_plotly(output_path)


def eval_electricity_balance(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the electricity production & demand by country and year.

    Electricity production is sum of power supply from
      - "Generator" components,
      - "Link" supply with AC bus_carrier, and
      - Pump-Hydro-Storage and Run-Of-River inflows

    Writes 2 Excel files and 1 BarChart per country.
    """
    ac_load_split = collect_myopic_statistics(networks, statistic="ac_load_split")
    # some technologies are included in Links, or do not count as demand
    link_carrier = [Carrier.grid_losses, "industry + hh & services load"]
    ac_load_split = ac_load_split.drop(link_carrier, level=DataModel.CARRIER)

    phs_hydro = collect_myopic_statistics(
        networks,
        "phs_split",
        carrier=[Carrier.phs_dispatched_power_inflow, Carrier.hydro_dispatched_power],
    )

    bev_demand = collect_myopic_statistics(
        networks, statistic="bev_v2g", carrier=[Carrier.bev_demand]
    )

    transmission_or_storage_links = [
        "",
        "BEV charger",
        "battery charger",
        "battery discharger",
        "DC",
        "V2G",
    ]
    ac_balance = collect_myopic_statistics(
        networks,
        statistic="energy_balance",
        comps="Link",
        bus_carrier=BusCarrier.AC,
    ).drop(transmission_or_storage_links, level=DataModel.CARRIER, errors="ignore")

    ac_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.AC,
    ).drop(Carrier.value_lost_load, level=DataModel.CARRIER, errors="ignore")

    trade_statistics = _fetch_electricity_trade_statistics(networks, net_trade=False)

    metric = Metric(
        metric_name="Power Balance",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[ac_load_split, bev_demand, ac_balance, ac_generation, phs_hydro]
        + trade_statistics,
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = "e_production"  # e_demand
    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "elec_balance_{location}"
    metric.defaults.plotly.cutoff = 0.0001
    metric.defaults.plotly.pattern = dict.fromkeys(
        [
            Group.import_foreign,
            Group.export_foreign,
            Group.import_domestic,
            Group.export_domestic,
        ],
        "/",
    )
    metric.defaults.plotly.category_orders = (
        Group.nuclear_power,
        Group.chp_biomass,
        Group.pp_thermal,
        Group.pv,
        Group.wind,
        Group.phs_inflow,
        Group.ror,
        Group.import_domestic,
        Group.import_foreign,
        # --- zero ---
        Group.industry,
        Group.hh_and_services,
        Group.dac,
        Group.electrolysis,
        Group.heat_decentral,
        Group.heat_district,
        Group.transport,
        Group.export_domestic,
        Group.export_foreign,
    )
    metric.defaults.plotly.footnotes = (  # fixme: Footnote not displayed
        "Balance does not include grid losses due to storage technologies.",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def _subtract_rail_and_industry_demands_from_electricity_base_load(
    ac_demand: pd.DataFrame,
    rail_industry: pd.DataFrame,
    keep_regions: tuple = ("AT",),
    regional_clipping: bool = True,
) -> pd.DataFrame:
    """Calculate the base load without rail and industry demands.

    Subtract rail and industry demand from electricity (Base Load),
    because they were added from dedicated metrics already. Note,
    that this is not the same as "industry + hh & services load"
    obtained from the ac_load_split, because the ac load split has
    time aggregated values whereas here we work with time series.

    There is an inconsistency between the old and the new
    evaluation. The old evaluation contains base load
    values per country. The new evaluation contains base load
    values per node. It can happen, that one node has positive
    values, but the other has not. In this case, the new
    evaluation drops positive values for the one node and keeps
    negative values (demand) for the other node.
    The old evaluation contains the sum of both nodes. If the sum
    is positive, it will be dropped.
    So, if the positive node value is bigger than the negative
    node value, the sum is positive and the country has base load
    zero in the old evaluation. In the new evaluation, the negative
    node load values prevail.
    This leads to differences between the old and the new
    evaluation, and I am not sure which is correct.

    To circumnavigate this problem, base load clipping is applied
    after location aggregation.

    New:
    nodes = ['ES0 0', 'ES2 0']
    year = "2050"
    snapshots                     2015-01-01 00:00:00  2015-01-01 03:00:00  2015-01-01 06:00:00
    2050 ES0 0    electricity AC  -2731.878352         206.106368           249.330098
         ES2 0    electricity AC  -131.025686          -66.052936           -71.928966

    Only 206.106368 will be clipped to zero. -66.052 will appear in the results.

    Old:
    demand_unfiltered["Base Load"].iloc[0::3]
    snapshot
    2015-01-01 00:00:00     2862.904037
    2015-01-01 03:00:00     -140.053436  # clipped to zero
    2015-01-01 06:00:00     -177.401128  # clipped to zero
    2015-01-01 09:00:00     4510.239152
    2015-01-01 12:00:00     7362.349636

    Parameters
    ----------
    ac_demand
        A data frame with the electricity Load.
    rail_industry
        The rail and industry demands as time series.
    keep_regions
        A collection of country codes, that should not
        be clipped at country level.
    regional_clipping
        Whether to apply the clipping at region or at
        nodal level.

    Returns
    -------
    :
        The electricity carrier demand without rail and industry
        demands.
    """
    idx_base_load = pd.IndexSlice[:, :, "electricity", :]
    idx_rail = pd.IndexSlice[:, :, "electricity rail", :]
    idx_industry = pd.IndexSlice[:, :, "industry", :]
    base_load = (
        ac_demand.loc[idx_base_load, :]
        - rail_industry.loc[idx_rail, :].to_numpy()
        - rail_industry.loc[idx_industry, :].to_numpy()
    )
    if any(base_load > 0):
        logger = logging.getLogger(__name__)
        logger.warning(
            "Dropping positive base load values. This leads to "
            "imbalanced time steps in the time series balances."
        )

    if regional_clipping:
        nodes = base_load.index.get_level_values(DataModel.LOCATION)
        regions = [
            node if node.startswith(keep_regions) else node[:2] for node in nodes
        ]

        clipped = []
        for _, region_year in base_load.groupby([regions, DataModel.YEAR]):
            country_sum = region_year.sum()
            positive_snapshots = country_sum[country_sum > 0].index
            region_year[positive_snapshots] = 0.0
            clipped.append(region_year)

        clipped_base_load = pd.concat(clipped)
        ac_demand = pd.concat(
            [
                ac_demand.drop("electricity", level=DataModel.CARRIER, errors="ignore"),
                clipped_base_load,
            ]
        )

    else:  # nodal clipping
        ac_demand.loc[idx_base_load, :] = base_load.clip(upper=0)

    return ac_demand


def _calculate_inflexible_demand(statistics: pd.DataFrame) -> pd.DataFrame:
    """Calculate the inflexible demand time series.

    Parameters
    ----------
    statistics
        The combined electricity production and demand statistics
        without location aggregation and with the original carrier
        and bus_carrier names.

    Returns
    -------
    :
        The inflexible demand statistic as a time series.
    """
    # add inflexible demand as a sum of technologies that cannot be
    # shut down (load shedding not allowed)
    mapping = get_mapping("electricity", "external")
    # fixme: if the electricity mapping changes, bugs are introduced
    #  because the inflexible carrier might change. That's very bad and
    #  must be prevented. However, this is the case in the old Toolbox
    #  too.
    inflexible_carrier = [
        carrier
        for carrier, group in mapping.items()
        if group in [Group.industry, Group.base_load, Group.transport]
    ]
    demand = statistics.mul(statistics.le(0))
    # BEV charger can be load shifted (= flexible). Need to exclude
    # them, because they are mapped to "Transport" group.
    inflexible_carrier.remove(Carrier.bev_charger)
    inflexible_demand = filter_by(demand, carrier=inflexible_carrier)
    inflexible_demand = rename_aggregate(inflexible_demand, "Inflexible Demand")

    return inflexible_demand


def _fetch_electricity_trade_statistics(
    networks: dict, aggregate_time: str | bool = "sum", net_trade: bool = True
) -> list:
    """Calculate foreign and domenstic trade statistics.

    Parameters
    ----------
    networks
        The loaded networks in a dictionary.
    aggregate_time
        The time aggregation to apply. Pass False to return time series.
    net_trade
        Whether to combine domestic and foreign trade statistics. The
        import and export sides will still be separate.

    Returns
    -------
    :
        Statistics for electricity trade.
    """
    ac_import_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="import",
        bus_carrier=BusCarrier.AC,
        aggregate_time=aggregate_time,
    )
    alias = Group.import_net if net_trade else Group.import_foreign
    ac_import_foreign = rename_aggregate(ac_import_foreign, alias)

    ac_import_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        bus_carrier=BusCarrier.AC,
        aggregate_time=aggregate_time,
    )
    alias = Group.import_net if net_trade else Group.import_domestic
    ac_import_domestic = rename_aggregate(ac_import_domestic, alias)

    ac_export_foreign = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.FOREIGN,
        direction="export",
        bus_carrier=BusCarrier.AC,
        aggregate_time=aggregate_time,
    )
    alias = Group.export_net if net_trade else Group.export_foreign
    ac_export_foreign = rename_aggregate(ac_export_foreign, alias)

    ac_export_domestic = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="export",
        bus_carrier=BusCarrier.AC,
        aggregate_time=aggregate_time,
    )
    alias = Group.export_net if net_trade else Group.export_domestic
    ac_export_domestic = rename_aggregate(ac_export_domestic, alias)

    return [
        ac_import_foreign,
        ac_import_domestic,
        ac_export_foreign,
        ac_export_domestic,
    ]
