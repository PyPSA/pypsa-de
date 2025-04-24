"""Module for transformation Sankey evaluations."""

from pathlib import Path

from statistic import collect_myopic_statistics
from utils import filter_by

# def eval_sankey_detail(
#     result_path: str | Path,
#     networks: dict,
#     subdir: str | Path = "esm_run/evaluation",
# ) -> None:
#     """"""
#     # set(carrier_target_dict)
#     carrier_ac = {
#         "H2 Electrolysis",
#         "H2 Fuel Cell",
#         "OCGT",
#         "PHS Dispatched Power from Stored",
#         "PHS Stored Power",
#         "V2G energy back to network",
#         "V2G energy demand",
#         "battery charger",
#         "battery discharger",
#         "coal power plant",
#         "coal power plant (CC)",
#         "lignite power plant",
#         "lignite power plant (CC)",
#         "oil power plant",
#         "residential rural CH4-powered PEMFC with internal SMR",
#         "residential rural H2-powered PEMFC",
#         "residential rural micro gas CHP",
#         "services rural CH4-powered PEMFC with internal SMR",
#         "services rural H2-powered PEMFC",
#         "services rural micro gas CHP",
#         "urban central air heat pump",
#         "urban central coal CHP CC electric",
#         "urban central coal CHP electric",
#         "urban central gas CHP CC electric",
#         "urban central gas CHP electric",
#         "urban central lignite CHP CC electric",
#         "urban central lignite CHP electric",
#         "urban central resistive heater",
#         "urban central solid biomass CHP",
#         "urban central solid biomass CHP CC",
#     }
#
#     ac_starts = [
#         "Battery",
#         "H2 Electrolysis",
#         "PHS",
#         "V2G",
#         "urban central air heat pump",
#         "urban central resistive heater",
#     ]
#     ac = collect_myopic_statistics(
#         networks,
#         "supply",
#         comps=["Link", "Store", "StorageUnit"],
#         bus_carrier=["AC", "Li ion", "Battery"],
#         carrier=ac_starts,
#     )
#
#     filter_by(ac, year="2020", location=[f"AT0 {i}" for i in range(10)]).sum() * 1e-6
#
#     # input: from eval results
#     #  - power balance
#     #    = electricity_production
#     #    = electricity_demand
#     #  - methane balance
#     #    = ch4 production
#     #    = ch4 demand
#     #  - h2 balance
#     #    = h2 production
#     #    = h2 demand
#     #  -
#     # input: from network
#     #  -  # numpydoc ignore=PR01
