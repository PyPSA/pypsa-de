"""Evaluate the nodal final energy demand by bus carrier."""

from pathlib import Path

from evals.constants import BusCarrier, Carrier, DataModel
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    calculate_input_share,
    filter_by,
    filter_for_carrier_connected_to,
    get_heat_loss_factor,
    split_urban_heat_losses_and_consumption,
)


def view_final_energy_demand(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the final energy demand per country.

    Parameters
    ----------
    result_path
    networks
    config
    subdir

    Returns
    -------
    :
    """
    #         _fetch_fed_transport(networks, result_path),
    #         _fetch_fed_industry(networks, result_path),
    # done    _fetch_fed_decentral_heat(networks),
    # done    _fetch_fed_district_heat(networks),
    # done  _fetch_fed_homes_and_trade(networks),  but broken

    link_supply_rural_heat = (
        collect_myopic_statistics(
            networks,
            comps="Link",
            statistic="energy_balance",
        )
        .pipe(filter_for_carrier_connected_to, BusCarrier.HEAT_RURAL, kind="supply")
        .pipe(calculate_input_share, BusCarrier.HEAT_RURAL)
    )

    generator_supply_rural_heat = collect_myopic_statistics(
        networks,
        comps="Generator",
        statistic="supply",
        bus_carrier=BusCarrier.HEAT_RURAL,
    )

    load_withdrawal_urban_heat = collect_myopic_statistics(
        networks,
        "withdrawal",
        comps="Load",
        bus_carrier=[BusCarrier.HEAT_URBAN_CENTRAL, BusCarrier.HEAT_URBAN_DECENTRAL],
    ).drop(
        Carrier.low_temperature_heat_for_industry,
        level=DataModel.CARRIER,
        # errors="ignore",
    )

    # # The Toolbox drops Italian urban heat technologies for unknown reasons.
    # load_withdrawal_urban_heat = load_withdrawal_urban_heat.drop(["IT0", "IT1", "IT2"], level=DataModel.LOCATION)
    # # todo: Is this correct? They probably had a good reason for that, but I just can't see it.

    loss_factor = get_heat_loss_factor(networks)
    load_split_urban_heat = split_urban_heat_losses_and_consumption(
        load_withdrawal_urban_heat, loss_factor
    )

    fed_homes_and_trade = collect_myopic_statistics(
        networks, statistic="ac_load_split"
    ).pipe(filter_by, carrier=Carrier.domestic_homes_and_trade)

    # todo: couldn't his be much easier? If we simply query for withdrawal at heat buses?
    from evals.utils import filter_by

    _ = (
        collect_myopic_statistics(
            networks,
            comps=("Generator", "Load"),
            statistic="energy_balance",
            bus_carrier=BusCarrier.heat_buses(),
            aggregate_time=False,
        )
        .pipe(
            filter_by,
            carrier="urban decentral biomass boiler",
            year="2020",
            location="CH",
        )
        .T
    )
    # fixme: "urban decentral biomass boiler" supplies to solid biomass and draws from heat!!!

    # todo: need to map carrier names to sector names in grouped barchart
    exporter = Exporter(
        statistics=[
            link_supply_rural_heat.mul(-1),
            generator_supply_rural_heat,
            load_split_urban_heat,
            fed_homes_and_trade,
        ],
        statistics_unit="MWh",
        view_config=config["view"],
    )

    # view specific static settings:
    exporter.defaults.plotly.chart = ESMBarChart
    exporter.defaults.excel.pivot_index = [
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    exporter.defaults.plotly.plot_category = DataModel.BUS_CARRIER
    exporter.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    # metric.defaults.plotly.file_name_template = "final_demand_{location}"
    # metric.defaults.plotly.legend_header = "Energy Carrier"
    # metric.cfg.plotly.footnotes = (
    #     "International aviation and navigation are not included.",
    #     "",
    # )
    # metric.defaults.plotly.category_orders = (
    #     Group.coal,
    #     Group.solar_thermal,
    #     Group.h2,
    #     Group.biomass,
    #     Group.oil,
    #     Group.ch4,
    #     Group.electrictiy,
    #     Group.heat_district,
    # )
    # metric.defaults.plotly.xaxis_title = "<b>Years</b>"
    # metric.defaults.plotly.cutoff = 0.04

    exporter.export(result_path, subdir=subdir)
