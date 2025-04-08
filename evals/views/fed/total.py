"""Evaluate the nodal final energy demand by bus carrier."""

from pathlib import Path

from evals.constants import DataModel
from evals.fileio import (
    Exporter,
    prepare_co2_emissions,
)
from evals.plots.barchart import ESMBarChart


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
    #         _fetch_fed_decentral_heat(networks),
    #         _fetch_fed_district_heat(networks),
    #         _fetch_fed_homes_and_trade(networks),
    # prepare_costs(result_path, 25)
    prepare_co2_emissions(result_path, networks["2020"], "THIA")
    # prepare_industry_demand()
    # prepare_nodal_energy(result_path, sub_directory="../resources")

    exporter = Exporter(
        statistics=[], statistics_unit="MWh", view_config=config["view"]
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
