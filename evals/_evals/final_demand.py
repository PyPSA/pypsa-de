"""Evaluate the Final Energy Demand."""

from functools import partial
from pathlib import Path

import pandas as pd
from constants import (
    TITLE_SUFFIX,
    BusCarrier,
    Carrier,
    DataModel,
    Group,
    Mapping,
)
from fileio import prepare_industry_demand, prepare_nodal_energy
from metric import Metric
from plots.barchart import ESMBarChart
from plots.facetbars import ESMGroupedBarChart
from statistic import (
    collect_myopic_statistics,
    get_location_and_carrier_and_bus_carrier,
)
from utils import (
    drop_from_multtindex_by_regex,
    filter_by,
    insert_index_level,
    make_evaluation_result_directories,
    rename_aggregate,
)

from evals.heat import get_heat_loss_factor, split_heat_hh_service_losses
from evals.industry import (
    INDUSTRY_BUS_CARRIER,
    fetch_industry_demand_statistics,
)


def eval_fed(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """
    Evaluate the final energy demand per country.

    This evaluation uses the carrier index level to track sectors.

    Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See evals docstring for parameter description.
    """
    fed_statistics = final_energy_demand_statistics(networks, result_path)

    metric = Metric(
        metric_name="Final Energy Demand Total",
        is_unit="MWh",
        to_unit="TWh",
        statistics=fed_statistics,
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = Mapping.bus_carrier

    metric.defaults.excel.title = title
    metric.defaults.excel.pivot_index = [
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]

    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.plot_category = DataModel.BUS_CARRIER
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    metric.defaults.plotly.file_name_template = "final_demand_{location}"
    metric.defaults.plotly.legend_header = "Energy Carrier"
    # metric.cfg.plotly.footnotes = (
    #     "International aviation and navigation are not included.",
    #     "",
    # )
    metric.defaults.plotly.category_orders = (
        Group.coal,
        Group.solar_thermal,
        Group.h2,
        Group.biomass,
        Group.oil,
        Group.ch4,
        Group.electrictiy,
        Group.heat_district,
    )
    metric.defaults.plotly.xaxis_title = "<b>Years</b>"
    metric.defaults.plotly.cutoff = 0.04

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_fed_sectoral(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """
    Evaluate the final energy demand per country and sector.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See eval docstring for parameter description.
    """
    fed_statistics = final_energy_demand_statistics(networks, result_path)

    metric = Metric(
        metric_name="Final Energy Demand Total by Carrier and Sector",
        is_unit="MWh",
        to_unit="TWh",
        statistics=fed_statistics,
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = Mapping.bus_carrier

    metric.defaults.excel.title = title
    metric.defaults.excel.pivot_index = [DataModel.LOCATION, DataModel.BUS_CARRIER]
    metric.defaults.excel.pivot_columns = [DataModel.BUS_CARRIER, DataModel.YEAR]
    metric.defaults.excel.axis_labels = [
        DataModel.YEAR.title(),
        metric.df.attrs["unit"],
    ]

    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMGroupedBarChart
    metric.defaults.plotly.facet_column = DataModel.CARRIER  # sector subplots
    metric.defaults.plotly.plot_category = DataModel.BUS_CARRIER  # categories
    metric.defaults.plotly.file_name_template = "final_demand_sec_{location}"
    metric.defaults.plotly.legend_header = "Energy Carrier"
    metric.defaults.plotly.xaxis_title = ""  # skip to use sector names
    metric.defaults.plotly.cutoff = 0.0  # None in original Toolbox
    metric.defaults.plotly.category_orders = (
        Group.coal,
        Group.solar_thermal,
        Group.biomass,
        Group.oil,
        Group.ch4,
        Group.h2,
        Group.electrictiy,
        Group.heat_district,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_fed_biomass(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the Final Energy Demand of Biomass per country."""
    ac_input = collect_myopic_statistics(
        networks,
        "energy_input",
        comps="Link",
        bus_carrier=[BusCarrier.AC],
        # fixme: why do we not include the full amounts of biomass
        #  needed for AC production? In the case of AC, the efficiency
        #  is applied, whereas for heat production the full amount of
        #  solid biomass is given to produce heat (including losses).
        # include_losses=True,
    )
    biomass_for_ac = filter_by(ac_input, bus_carrier=[BusCarrier.SOLID_BIOMASS])

    heat_input = collect_myopic_statistics(
        networks,
        "energy_input",
        comps="Link",
        bus_carrier=[
            BusCarrier.HEAT_RURAL_RESIDENTIAL,
            BusCarrier.HEAT_RURAL_SERVICES,
            BusCarrier.HEAT_URBAN_CENTRAL,
            # BusCarrier.HEAT_URBAN_SERVICES,
            # BusCarrier.HEAT_URBAN_RESIDENTIAL,
        ],
        include_losses=True,
    ).drop("", level=DataModel.CARRIER, errors="ignore")
    biomass_for_heat = filter_by(heat_input, bus_carrier=[BusCarrier.SOLID_BIOMASS])
    # need to rename CHP carrier to keep them separate from AC
    _heat_mapping = {
        Carrier.chp_urban_central_solid_biomass: Group.heat_district,
        Carrier.chp_urban_central_solid_biomass_cc: Group.heat_district,
    }
    biomass_for_heat = rename_aggregate(biomass_for_heat, _heat_mapping)

    _industry = prepare_industry_demand(result_path, networks)
    biomass_for_industry = filter_by(_industry, carrier=BusCarrier.SOLID_BIOMASS)

    # need to aggregate sectors to fit the expected data model
    groups = [DataModel.YEAR, DataModel.LOCATION, DataModel.CARRIER]
    biomass_for_industry = biomass_for_industry.groupby(groups).sum()
    biomass_for_industry.index.names = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    biomass_for_industry = insert_index_level(
        biomass_for_industry, Group.industry, DataModel.CARRIER, pos=-1
    )

    metric = Metric(
        metric_name="Final Demand of Biomass",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[biomass_for_ac, biomass_for_heat, biomass_for_industry],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX

    metric.defaults.mapping = "biomass"
    metric.defaults.excel.chart_title = title

    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.legend_header = "Final Biomass Demand"
    metric.defaults.plotly.file_name_template = "biomass_balance_{location}"
    metric.defaults.plotly.cutoff = 0.02  # TWh
    metric.defaults.plotly.category_orders = (
        # order production from outside to zero
        Group.industry,
        Group.electrictiy,
        Group.heat_decentral,
        Group.heat_district,
    )
    metric.defaults.plotly.footnotes = (
        "The use of biomass for heat and power in COP plants is calculated "
        "based on the respective efficiencies.",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_fed_building_heat(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the final energy demand of heat for buildings."""
    fed_district_heat = _fetch_fed_district_heat(networks)

    fed_decentral = _fetch_fed_decentral_heat(networks)

    metric = Metric(
        metric_name=(
            "Final Energy Demand for Room Heat and Warm Water Households and Services"
        ),
        is_unit="MWh",
        to_unit="TWh",
        statistics=[fed_decentral, fed_district_heat],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = Mapping.bus_carrier

    metric.defaults.excel.title = title
    metric.defaults.excel.pivot_index = [
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]

    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.plot_category = DataModel.BUS_CARRIER
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    metric.defaults.plotly.file_name_template = "building_sec_heat_{location}"
    metric.defaults.plotly.legend_header = "Energy Carrier"
    metric.defaults.plotly.cutoff = 0.04
    # metric.cfg.plotly.footnotes = (
    #     "International aviation and navigation are not included.",
    #     "",
    # )
    metric.defaults.plotly.category_orders = (
        Group.coal,
        Group.solar_thermal,
        Group.h2,
        Group.biomass,
        Group.oil,
        Group.ch4,
        Group.electrictiy,
        Group.heat_district,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_fed_hh_services(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the total final energy demand of heat."""
    fed_district_heat = _fetch_fed_district_heat(networks)
    fed_decentral = _fetch_fed_decentral_heat(networks)
    fed_hh_services = _fetch_fed_homes_and_trade(networks)

    metric = Metric(
        metric_name="Final Energy Demand Households & Services",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[fed_district_heat, fed_decentral, fed_hh_services],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.defaults.mapping = Mapping.bus_carrier

    metric.defaults.excel.title = title
    metric.defaults.excel.pivot_index = [
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]

    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.plot_category = DataModel.BUS_CARRIER
    metric.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    metric.defaults.plotly.file_name_template = "building_sec_tot_{location}"
    metric.defaults.plotly.legend_header = "Energy Carrier"
    metric.defaults.plotly.cutoff = 0.04
    # metric.cfg.plotly.footnotes = (
    #     "International aviation and navigation are not included.",
    #     "",
    # )
    metric.defaults.plotly.category_orders = (
        Group.coal,
        Group.solar_thermal,
        Group.h2,
        Group.biomass,
        Group.oil,
        Group.ch4,
        Group.electrictiy,
        Group.heat_district,
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def final_energy_demand_statistics(networks: dict, result_path: Path) -> list:
    """
    Return all final energy demand statistics.

    The same energy demand statistics are used in the total
    and in the sectoral final energy demand evaluations.

    Parameters
    ----------
    networks
        The loaded networks.
    result_path
        The path to the results directory, needed to locate resource
        files.

    Returns
    -------
    :
        A list with all final energy demand statics. The statistics have
        the sector names in the "carrier" index level.
    """
    return [
        _fetch_fed_transport(networks, result_path),
        _fetch_fed_industry(networks, result_path),
        _fetch_fed_decentral_heat(networks),
        _fetch_fed_district_heat(networks),
        _fetch_fed_homes_and_trade(networks),
    ]


def _fetch_fed_transport(networks: dict, result_path: Path) -> pd.Series | pd.DataFrame:
    """
    Calculate the final energy demand for transport technologies.

    Parameters
    ----------
    networks
        The loaded networks.
    result_path
        The path to the results directory, needed to locate resource
        files.

    Returns
    -------
    :
        The combined final energy demand for transport technologies from
        the nodal energy demand resource file, Links, road freight Loads
        and BEV passenger demand.
    """
    nodal_energy = prepare_nodal_energy(result_path)
    nodal_transport_carrier = [
        Carrier.ft_domestic_aviation,
        Carrier.ft_domestic_navigation,
        Carrier.ft_rail,
        Carrier.road_freight_ft,
        "H2 domestic aviation",
        "H2 domestic navigation",
        "H2 rail",
        Carrier.electricity_rail,
        Carrier.ch4_navigation_domestic,
    ]
    transport_nodal = filter_by(nodal_energy, carrier=nodal_transport_carrier)

    # Transport technologies are like energy "consumers", the FED is at
    # their supply side. In addition, Fischer-Tropsch buses have
    # location "EU" and not the country of consumption. We calculate
    # them using the bus location from port 1.
    fed_transport_links = collect_myopic_statistics(
        networks,
        statistic="energy_balance",  # FixMe: should be withdrawal!
        comps="Link",
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="1"),
        bus_carrier=[
            BusCarrier.FT_1,
            BusCarrier.FT_2,
            BusCarrier.AC,
            BusCarrier.CH4,
            BusCarrier.H2,
        ],
        carrier=[
            Carrier.cng_long,
            Carrier.cng_short,
            Carrier.phev_long,
            Carrier.phev_short,
            Carrier.hev_long,
            Carrier.hev_short,
            Carrier.ice_long,
            Carrier.ice_short,
            Carrier.fcev_long,
            Carrier.fcev_short,
        ],
    )

    # fixme: I am not sure why, but transport technologies have mixed
    #  signs. For example, 'CNG long' is supply and 'CNG short' is
    #  withdrawal, although both Links have the same bus ports and
    #  respective port efficiencies. We fix this using absolute values,
    #  but an explanation would be nice.
    #  See the bugfix_large-gas-supply-from-CNG.ipynb
    fed_transport_links = fed_transport_links.abs()

    road_freight = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Load",
        bus_carrier=[BusCarrier.H2, BusCarrier.AC, BusCarrier.CH4],
        carrier=[
            Carrier.road_freight_h2,
            Carrier.road_freight_ac,
            Carrier.road_freight_ch4,
        ],
    )

    bev_demand = collect_myopic_statistics(
        networks, statistic="bev_v2g", carrier=[Carrier.bev_demand]
    )
    bev_demand = bev_demand.clip(upper=0).abs()

    fed_transport = pd.concat(
        [transport_nodal, fed_transport_links, road_freight, bev_demand]
    )

    # level "carrier" holds the "sector" names in this evaluation
    return rename_aggregate(fed_transport, Group.transport)


def _fetch_fed_industry(networks: dict, result_path: Path) -> pd.Series | pd.DataFrame:
    """
    Calculate the final energy demand for the industry sector.

    Parameters
    ----------
    networks
        The loaded networks.
    result_path
        The path to the results directory, needed to locate resource
        files.

    Returns
    -------
    :
        A reshaped and renamed version of the industry statistics.
    """
    industry_statistics = fetch_industry_demand_statistics(networks, result_path)
    # This evaluation compares transport, heat, and industry sectors.
    # The granular industry sectors (Chemical, Iron and Steel, ...) are
    # not of interest, so we replace it by 'Industry' here.
    fed_industry = rename_aggregate(
        industry_statistics, Group.industry, level=DataModel.BUS_CARRIER
    )
    # we also rename the nodal industry file bus_carrier names here
    fed_industry = rename_aggregate(fed_industry, INDUSTRY_BUS_CARRIER)

    fed_industry = fed_industry.swaplevel(DataModel.CARRIER, DataModel.BUS_CARRIER)
    fed_industry.index.names = DataModel.YEAR_IDX_NAMES

    return fed_industry


def _fetch_fed_decentral_heat(networks: dict) -> pd.Series | pd.DataFrame:
    """
    Calculate the final energy demand for the decentral heat.

    Heat from Generator coomponents is assumend to be solar thermal
    technologies. This assumption is checked using an assert statement.

    Parameters
    ----------
    networks
        The loaded networks.

    Returns
    -------
    :
        Collective energy demand for different bus_carrier from
        decentral heat Links and Generators.
    """
    fed_heat_multiport = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_input",
        bus_carrier=[
            BusCarrier.HEAT_RURAL_SERVICES,
            BusCarrier.HEAT_RURAL_RESIDENTIAL,
        ],
        include_losses=True,  # return port "1" energies with efficiency = 1.0
    )
    fed_heat_multiport = drop_from_multtindex_by_regex(
        fed_heat_multiport, "DAC|water tank"
    )

    fed_heat_generation = collect_myopic_statistics(
        networks,
        "supply",
        comps="Generator",
        bus_carrier=[
            BusCarrier.HEAT_RURAL_SERVICES,
            BusCarrier.HEAT_RURAL_RESIDENTIAL,
        ],
    )
    assert all(
        fed_heat_generation.index.unique(DataModel.BUS_CARRIER)
        == ["residential rural heat", "services rural heat"]
    )
    fed_heat_generation = rename_aggregate(
        fed_heat_generation, Group.solar_thermal, level=DataModel.BUS_CARRIER
    )

    fed_heat = pd.concat([fed_heat_multiport, fed_heat_generation])
    # level "carrier" holds the "sector" names in this evaluation
    return rename_aggregate(fed_heat, Group.hh_and_services)


def _fetch_fed_district_heat(networks: dict) -> pd.Series | pd.DataFrame:
    """
    Calculate the final energy demand for central (district) heat.

    Grid losses and low temperature heat for industry parts are
    excluded.

    Parameters
    ----------
    networks
        The loaded networks.

    Returns
    -------
    :
        Collective energy demand for different bus_carrier from
        district heat Loads without grid losses and LTH for industry parts.
    """
    # fixme: why does the Toolbox use "p_set" instead of "p" for
    #  heat loads?
    fed_district_heat_loads = collect_myopic_statistics(
        networks,
        "withdrawal",
        comps="Load",
        bus_carrier=[
            BusCarrier.HEAT_URBAN_CENTRAL,
            BusCarrier.HEAT_URBAN_SERVICES,
            BusCarrier.HEAT_URBAN_RESIDENTIAL,
        ],
    ).drop(
        Carrier.low_temperature_heat_for_industry,
        level=DataModel.CARRIER,
        errors="ignore",
    )

    # some countries can never have district heat networks and the
    # evaluation should not show those.
    # Fixme: Why does the model have Loads for IT district heating? And
    #  how is it possible that those loads remain unmet?
    impossible_district_heat_nodes = ["IT0 0", "IT1 0"]
    fed_district_heat_loads = fed_district_heat_loads.drop(
        impossible_district_heat_nodes, level=DataModel.LOCATION, errors="ignore"
    )

    heat_loss_factor = get_heat_loss_factor(networks)
    fed_district_heat_loads = split_heat_hh_service_losses(
        fed_district_heat_loads, heat_loss_factor
    )
    fed_district_heat_loads = (
        fed_district_heat_loads.squeeze().abs()
    )  # plot demand as positives
    fed_district_heat_loads = filter_by(
        fed_district_heat_loads, carrier=Carrier.hh_and_services
    )

    # level "carrier" holds the "sector" names in this evaluation
    return rename_aggregate(fed_district_heat_loads, Group.hh_and_services)


def _fetch_fed_homes_and_trade(networks: dict) -> pd.Series | pd.DataFrame:
    """
    Calculate the final energy demand for domestic homes and trade.

    Parameters
    ----------
    networks
        The loaded networks.

    Returns
    -------
    :
        AC demand for heat from domestic homes and trade with the sector
        name "Households & Services".
    """
    return (
        collect_myopic_statistics(networks, statistic="ac_load_split")
        .pipe(filter_by, carrier=Carrier.domestic_homes_and_trade)
        .pipe(rename_aggregate, Group.hh_and_service)
        .mul(-1)
    )  # plot demand upwards
    # return rename_aggregate(fed_homes_and_trade, Group.hh_and_services)
