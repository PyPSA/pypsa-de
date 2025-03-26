"""Module for transport sector evaluations."""

from functools import partial
from pathlib import Path

from frozendict import frozendict

from constants import (
    TITLE_SUFFIX,
    BusCarrier,
    Carrier,
    DataModel,
    Group,
    Mapping,
)
from fileio import prepare_nodal_energy
from metric import Metric
from plots.barchart import ESMBarChart
from plots.facetbars import ESMGroupedBarChart
from statistic import (
    collect_myopic_statistics,
    get_location_and_carrier_and_bus_carrier,
)
from utils import (
    filter_by,
    make_evaluation_result_directories,
    rename_aggregate,
)

TRANSPORT_SECTORS = frozendict(
    {
        Carrier.bev_demand: "Passenger cars",
        Carrier.cng_long: "Passenger cars",
        Carrier.cng_short: "Passenger cars",
        Carrier.ft_rail: "Rail",
        Carrier.road_freight_ft: "Road freight",
        Carrier.electricity_rail: "Rail",
        Carrier.road_freight_ac: "Road freight",
        Carrier.fcev_long: "Passenger cars",
        Carrier.fcev_short: "Passenger cars",
        Carrier.ch4_navigation_domestic: "Navigation",
        Carrier.road_freight_ch4: "Road freight",
        "H2 domestic aviation": "Aviation",
        "H2 rail": "Rail",
        "H2 domestic navigation": "Navigation",
        Carrier.road_freight_h2: "Road freight",
        Carrier.hev_long: "Passenger cars",
        Carrier.hev_short: "Passenger cars",
        Carrier.ice_long: "Passenger cars",
        Carrier.ice_short: "Passenger cars",
        Carrier.ft_domestic_aviation: "Aviation",
        Carrier.ft_domestic_navigation: "Navigation",
        Carrier.phev_long: "Passenger cars",
        Carrier.phev_short: "Passenger cars",
    }
)


def eval_transport_total(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the transport demand per country.

    Writes 2 Excel files and 1 Barchart for the total energy demand
    for transport, and  2 Excel files and 1 GroupedBarchart for
    the sectoral transport energy demand

    Notes
    -----
    See eval docstring for parameter description.
    """
    statistics = _fetch_transport_statistics(networks, result_path)

    metric = Metric(
        "Energy Demand",
        is_unit="MWh",
        to_unit="TWh",
        statistics=statistics,
    )

    title = "Final Energy Demand Transport" + TITLE_SUFFIX
    metric.cfg.mapping = Mapping.bus_carrier

    metric.cfg.excel.title = title
    metric.cfg.excel.pivot_index = [
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]

    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.plot_category = DataModel.BUS_CARRIER
    metric.cfg.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    metric.cfg.plotly.file_name_template = "trans_sec_tot_{location}"
    metric.cfg.plotly.legend_header = "Energy Carrier"
    metric.cfg.plotly.footnotes = (
        "International aviation and navigation are not included.",
        "",
    )
    metric.cfg.plotly.xaxis_title = "<b>Years</b>"
    metric.cfg.plotly.cutoff = 0.0001

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_transport_sectoral(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate transport demand by sector.

    Writes 2 Excel files and 1 Barchart for the total energy demand
    for transport, and  2 Excel files and 1 GroupedBarchart for
    the sectoral transport energy demand

    Notes
    -----
    See eval docstring for parameter description.
    """
    statistics = _fetch_transport_statistics(networks, result_path)
    statistics_with_sector_as_carrier = [
        rename_aggregate(s, TRANSPORT_SECTORS) for s in statistics
    ]

    metric = Metric(
        "Energy Demand",
        is_unit="MWh",
        to_unit="TWh",
        statistics=statistics_with_sector_as_carrier,
    )

    title = "Final Energy Demand Transport by Carrier and Sector" + TITLE_SUFFIX
    metric.cfg.mapping = Mapping.bus_carrier

    metric.cfg.excel.chart_title = title

    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMGroupedBarChart
    metric.cfg.plotly.facet_column = DataModel.CARRIER  # sector subplots
    metric.cfg.plotly.plot_category = DataModel.BUS_CARRIER  # categories
    metric.cfg.plotly.file_name_template = "trans_sec_dem_{location}"
    metric.cfg.plotly.legend_header = "Energy Carrier"
    metric.cfg.plotly.xaxis_title = ""  # skip to use sector names
    metric.cfg.plotly.legend_font_size = 18
    metric.cfg.plotly.xaxis_font_size = 18
    metric.cfg.plotly.footnotes = (
        "<br>International aviation and navigation are not included.",
        "",
    )
    metric.cfg.plotly.category_orders = (
        Group.electrictiy,
        Group.h2,
        Group.ch4,
        Group.oil,
    )
    metric.cfg.plotly.cutoff = 0.00  # None in original Toolbox

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def _fetch_transport_statistics(networks: dict, result_path: Path) -> list:
    """Calculate the transport statistics.

    Parameters
    ----------
    networks
        The loaded postnetworks.
    result_path
        The path to the results directory, needed to locate resource
        files.

    Returns
    -------
    :
        A collection of transport sector statistics.
    """
    transport_buses = [BusCarrier.TRANSPORT_P, BusCarrier.TRANSPORT_P_LONG]

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
    nodal_transport = filter_by(nodal_energy, carrier=nodal_transport_carrier)

    # CNG demand has flipped signs in the network.
    cng_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Link",
        bus_carrier=[BusCarrier.CH4],
        carrier=[Carrier.cng_long, Carrier.cng_short],
    )

    p_hev_demand_without_ft = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Link",
        bus_carrier=[BusCarrier.AC] + transport_buses,
        carrier=[
            Carrier.phev_long,
            Carrier.phev_short,
            Carrier.hev_long,
            Carrier.hev_short,
        ],
    )

    ice_demand_without_ft = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Link",
        bus_carrier=transport_buses,
        carrier=[Carrier.ice_long, Carrier.ice_short],
    )

    fcev_demand = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Link",
        bus_carrier=[BusCarrier.H2] + transport_buses,
        carrier=[Carrier.fcev_long, Carrier.fcev_short],
    )

    # Fischer-Tropsch buses have location "EU" and not the country
    # of consumption. We calculate them using the bus
    # location from port 1 instead of bus 0 port location.
    vehicle_demand_ft = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Link",
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="1"),
        bus_carrier=BusCarrier.FT_1,
        carrier=[
            Carrier.ice_long,
            Carrier.ice_short,
            Carrier.phev_long,
            Carrier.phev_short,
            Carrier.hev_long,
            Carrier.hev_short,
        ],
    )

    road_freight = collect_myopic_statistics(
        networks,
        statistic="withdrawal",
        comps="Load",
        bus_carrier=[BusCarrier.H2, BusCarrier.AC, BusCarrier.CH4],
        carrier=[
            Carrier.road_freight_h2,
            Carrier.road_freight_ac,
            Carrier.road_freight_ch4,
            # "gas international navigation",
        ],
    )

    bev_demand = collect_myopic_statistics(
        networks, statistic="bev_v2g", carrier=[Carrier.bev_demand]
    )
    bev_demand = bev_demand.clip(upper=0).abs()

    return [
        nodal_transport,
        bev_demand,
        cng_demand,
        fcev_demand,
        ice_demand_without_ft,
        p_hev_demand_without_ft,
        vehicle_demand_ft,
        road_freight,
    ]
