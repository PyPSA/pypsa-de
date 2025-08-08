# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for Sankey diagram."""

import dataclasses
import re
from collections import defaultdict

import pandas as pd
import plotly
import pyam
from plotly.graph_objs import Figure, Sankey
from pyam.index import get_index_levels

from evals.constants import COLOUR, RUN_META_DATA
from evals.constants import DataModel as DM
from evals.plots._base import ESMChart
from evals.utils import (
    filter_by,
    prettify_number,
    rename_aggregate,
)

# Transformationsblöcke je ENergieträger
# Zusammenfassung Energieträger:
#  - AC (low voltage + AC) - uranium similar to primary but with additional step
#  - H2
#  - Gas
#  - Liquids (oil, methanol, NH3, electrobiofuels, naptha)
#  - Solids (waste, biomass, coal, lignite)
#  - Heat (central), connect decentral heat directly to FED
# Alle Losses in Grau und je Tranformationsblock (Energieträger)


BUS_CARRIER_COLORS = {
    "biogas": COLOUR.green_sage,
    "coal": COLOUR.grey_dark,
    "H2": COLOUR.green_mint,
    "NH3": COLOUR.yellow_canary,
    "lignite": COLOUR.brown_dark,
    "gas": COLOUR.brown_light,
    "municipal solid waste": COLOUR.grey_light,
    "AC": COLOUR.blue_celestial,
    "oil primary": COLOUR.red_deep,
    "rural heat": COLOUR.yellow_golden,
    "low voltage": COLOUR.blue_celestial,
    "solid biomass": COLOUR.green_sage,
    "uranium": COLOUR.orange_mellow,
    "urban central heat": COLOUR.yellow_golden,
    "urban decentral heat": COLOUR.yellow_golden,
    "EV battery": COLOUR.blue_celestial,
    "methanol": COLOUR.salmon,
    "oil": COLOUR.red_deep,
    "non-sequestered HVC": COLOUR.grey_light,
    "agriculture machinery oil": COLOUR.red_deep,
    "battery": COLOUR.blue_celestial,
    "ambient heat": COLOUR.yellow_golden,
    "home battery": COLOUR.blue_celestial,
    "industry methanol": COLOUR.salmon,
    "kerosene for aviation": COLOUR.red_deep,
    "shipping methanol": COLOUR.salmon,
    "gas for industry": COLOUR.brown_light,
    "naphtha for industry": COLOUR.red_deep,
    "solid biomass for industry": COLOUR.green_sage,
    "rural water tanks": COLOUR.yellow_golden,
    "urban central water pits": COLOUR.yellow_golden,
    "urban central water tanks": COLOUR.yellow_golden,
    "urban decentral water tanks": COLOUR.yellow_golden,
}

# idx = df_plot.sort_index(level="bus_carrier").droplevel("year").droplevel("location").index.drop_duplicates()
# idx = self._df.sort_index(level="carrier").index.drop_duplicates()
# dict.fromkeys(idx, ("", ""))
LINK_MAPPING = {
    # value - id (component, carrier, bus_carrier): (source, target) nodes
    ("Generator", "offwind-ac", "AC"): ("Wind Power", "AC Primary"),
    ("Generator", "onwind", "AC"): ("Wind Power", "AC Primary"),
    ("Generator", "ror", "AC"): ("Hydro Power", "AC Primary"),
    ("Generator", "solar", "AC"): ("Solar Power", "AC Primary"),
    ("Generator", "solar-hsat", "AC"): ("Solar Power", "AC Primary"),
    ("Generator", "offwind-dc", "AC"): ("Wind Power", "AC Primary"),
    ("Line", "Export Foreign", "AC"): ("AC Secondary Output", "Export"),
    ("Line", "Import Foreign", "AC"): ("Import", "AC Primary"),
    ("Line", "Export Domestic", "AC"): ("AC Secondary Output", "Export"),
    ("Line", "Import Domestic", "AC"): ("Import", "AC Primary"),
    ("Link", "CCGT methanol", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "CCGT methanol CC", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "DAC", "AC"): ("AC Secondary Output", "DAC"),
    ("Link", "H2 Fuel Cell", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "Haber-Bosch", "AC"): ("AC Secondary Output", "Haber-Bosch"),
    ("Link", "Haber-Bosch losses", "AC"): ("AC Secondary Output", "Losses"),
    ("Link", "OCGT", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "OCGT methanol", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "allam methanol", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "battery charger", "AC"): ("AC Secondary Output", "Battery"),
    ("Link", "battery charger losses", "AC"): ("AC Secondary Output", "Losses"),
    ("Link", "battery discharger", "AC"): ("Battery", "AC Secondary Input"),
    ("Link", "electricity distribution grid", "AC"): ("", ""),  # todo: drop
    ("Link", "electricity distribution grid losses", "AC"): (
        "AC Secondary Output",
        "Transmission Losses",
    ),
    ("Link", "methanolisation", "AC"): ("AC Secondary Output", "Methanolisation"),
    ("Link", "methanolisation losses", "AC"): ("AC Secondary Output", "Losses"),
    ("Link", "oil", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "waste CHP", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "waste CHP CC", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "CCGT", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "urban central coal CHP", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "urban central gas CHP", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "urban central gas CHP CC", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "urban central solid biomass CHP", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "urban central solid biomass CHP CC", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "coal", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "nuclear", "AC"): ("Nuclear Power", "AC Secondary Input"),
    ("Link", "solid biomass", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "urban central oil CHP", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "Export Foreign", "AC"): ("AC Secondary Output", "Export"),
    ("Link", "Import Foreign", "AC"): ("Import", "AC Primary"),
    ("Link", "lignite", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "urban central lignite CHP", "AC"): ("CHP", "AC Secondary Input"),
    ("Link", "Export Domestic", "AC"): ("AC Secondary Output", "Export"),
    ("Link", "Import Domestic", "AC"): ("Import", "AC Primary"),
    ("Link", "H2 Electrolysis", "AC"): ("AC Secondary Output", "Electrolysis"),
    ("Link", "H2 Electrolysis losses", "AC"): ("AC Secondary Output", "Losses"),
    ("Link", "H2 OCGT", "AC"): ("Power plant", "AC Secondary Input"),
    ("Link", "urban central H2 CHP", "AC"): ("CHP", "AC Secondary Input"),
    ("StorageUnit", "hydro", "AC"): ("", ""),  # todo: rename supply/demand
    ("StorageUnit", "PHS", "AC"): ("", ""),  # todo: rename supply/demand
    ("Link", "BEV charger", "EV battery"): ("", ""),
    ("Link", "V2G", "EV battery"): ("Car Battery", "AC Secondary Input"),
    ("Link", "V2G losses", "EV battery"): ("Car Battery", "Car Charger Losses"),
    ("Load", "land transport EV", "EV battery"): ("Transport", "Car Battery"),
    ("Store", "EV battery", "EV battery"): ("", ""),  # todo: drop
    ("Generator", "import H2", "H2"): ("Global Green Import", "H2 Primary"),
    ("Link", "Fischer-Tropsch", "H2"): ("H2 Secondary Output", "Oil Secondary Input"),
    ("Link", "Fischer-Tropsch losses", "H2"): (
        "H2 Secondary Output",
        "Fischer-Tropsch Losses",
    ),
    ("Link", "H2 Fuel Cell", "H2"): ("H2 Secondary Output", "CHP"),
    ("Link", "H2 Fuel Cell losses", "H2"): ("CHP", "CHP Losses"),
    ("Link", "Haber-Bosch", "H2"): ("H2 Secondary Output", "Haber-Bosch"),
    ("Link", "Haber-Bosch losses", "H2"): ("Haber-Bosch", "Haber-Bosch Losses"),
    ("Link", "Methanol steam reforming", "H2"): (
        "Methanol Steam Reforming",
        "H2 Secondary Input",
    ),
    ("Link", "Methanol steam reforming CC", "H2"): (
        "Methanol Steam Reforming",
        "H2 Secondary Input",
    ),
    ("Link", "SMR", "H2"): ("Steam Metehane Reforming", "H2 Secondary Input"),
    ("Link", "SMR CC", "H2"): ("Steam Metehane Reforming", "H2 Secondary Input"),
    ("Link", "Sabatier", "H2"): ("H2 Secondary Output", "Sabatier"),
    ("Link", "Sabatier losses", "H2"): ("Sabatier", "Sabatier Losses"),
    ("Link", "ammonia cracker", "H2"): ("Ammonia Cracker", "H2 Secondary Input"),
    ("Link", "electrobiofuels", "H2"): ("H2 Secondary Output", "Bio Fuels"),
    ("Link", "electrobiofuels losses", "H2"): ("Bio Fuels", "Bio Fuels Losses"),
    ("Link", "methanol-to-kerosene", "H2"): ("H2 Secondary Output", "Bio Fuels"),
    ("Link", "methanol-to-kerosene losses", "H2"): ("Bio Fuels", "Bio Fuels Losses"),
    ("Link", "methanolisation", "H2"): ("H2 Secondary Output", "Methanolisation"),
    ("Link", "methanolisation losses", "H2"): (
        "Methanolisation",
        "Methanolisation Losses",
    ),
    ("Link", "solid biomass to hydrogen", "H2"): (
        "Wood Hydrogen",
        "H2 Secondary Input",
    ),
    ("Link", "Export Foreign", "H2"): ("H2 Secondary Output", "Export"),
    ("Link", "H2 Electrolysis", "H2"): ("Electrolysis", "H2 Secondary Input"),
    ("Link", "Import Foreign", "H2"): ("H2 Import", "H2 Primary"),
    ("Link", "H2 OCGT", "H2"): ("H2 Secondary Output", "Power Plant"),
    ("Link", "H2 OCGT losses", "H2"): ("Power Plant", "Power Plant Losses"),
    ("Link", "urban central H2 CHP", "H2"): ("H2 Secondary Output", "CHP"),
    ("Link", "urban central H2 CHP losses", "H2"): ("CHP", "CHP Losses"),
    ("Link", "Export Domestic", "H2"): ("H2 Secondary Output", "Export"),
    ("Link", "Import Domestic", "H2"): ("H2 Import", "H2 Primary"),
    ("Load", "land transport fuel cell", "H2"): ("H2 Secondary Output", "Transport"),
    ("Load", "H2 for industry", "H2"): ("H2 Secondary Output", "Industry"),
    ("Store", "H2 Store", "H2"): ("", ""),  # todo: drop
    ("Generator", "import NH3", "NH3"): (
        "Global Green Import",
        "NH3 Primary",
    ),  # todo: regionalize imports?
    ("Link", "Haber-Bosch", "NH3"): ("Haber-Bosch", "NH3 Secondary Input"),
    ("Link", "ammonia cracker", "NH3"): ("NH3 Secondary Output", "Ammonia Cracker"),
    ("Link", "ammonia cracker losses", "NH3"): (
        "Ammonia Cracker",
        "Ammonia Cracker Losses",
    ),
    ("Load", "NH3", "NH3"): ("NH3 Secondary Output", "Agriculture"),
    ("Store", "ammonia store", "NH3"): ("", ""),  # todo: drop
    ("Link", "agriculture machinery oil", "agriculture machinery oil"): (
        "",
        "",
    ),  # todo: drop
    ("Load", "agriculture machinery oil", "agriculture machinery oil"): (
        "Oil Secondary Output",
        "Agriculture",
    ),
    ("Link", "rural air heat pump", "ambient heat"): ("Ambient Heat", "Decentral Heat"),
    ("Link", "rural ground heat pump", "ambient heat"): (
        "Ambient Heat",
        "Decentral Heat",
    ),
    ("Link", "urban decentral air heat pump", "ambient heat"): (
        "Ambient Heat",
        "Decentral Heat",
    ),
    ("Link", "electricity distribution grid", "ambient heat"): (
        "",
        "",
    ),  # fixme: should be losses
    ("Link", "urban central air heat pump", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),  # todo: jumper from Central heat to HH & Services
    ("Link", "urban central coal CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central gas boiler", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central ptes heat pump", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central solid biomass CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central solid biomass CHP CC", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central oil CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central gas CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "urban central lignite CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Central Heat",
    ),
    ("Link", "battery charger", "battery"): ("", ""),  # todo: already in supply. drop
    ("Link", "battery discharger", "battery"): (
        "",
        "",
    ),  # todo: already in supply. drop
    ("Link", "battery discharger losses", "battery"): (
        "",
        "",
    ),  # todo: already in supply. drop
    ("Store", "battery", "battery"): ("", ""),  # todo: drop
    ("Generator", "unsustainable biogas", "biogas"): ("Wet Biomass", "Biogas"),
    ("Generator", "biogas", "biogas"): ("Wet Biomass", "Biogas"),
    ("Link", "biogas to gas", "biogas"): ("", ""),  # todo: drop demand side
    ("Link", "biogas to gas CC", "biogas"): ("", ""),  # todo: drop demand side
    ("Generator", "coal", "coal"): ("", ""),  # todo: regionalize coal imports
    ("Link", "coal for industry", "coal"): ("", ""),  # todo: drop
    ("Link", "urban central coal CHP", "coal"): ("Coal Secondary Input", "CHP"),
    ("Link", "urban central coal CHP losses", "coal"): ("CHP", "CHP Losses"),
    ("Link", "coal", "coal"): ("Coal Secondary Input", "Power Plant"),
    ("Link", "coal losses", "coal"): ("Power Plant", "Power Plant Losses"),
    ("Store", "coal", "coal"): ("", ""),  # todo: drop
    ("Load", "coal for industry", "coal for industry"): (
        "Coal Secondary Output",
        "Industry",
    ),
    ("Generator", "production gas", "gas"): (
        "Domestic Fossil Production",
        "Methane Primary",
    ),
    ("Generator", "lng gas", "gas"): ("LNG Global Fossil Import", "Methane Primary"),
    ("Generator", "pipeline gas", "gas"): (
        "Pipeline Global Fossil Import",
        "Methane Primary",
    ),
    ("Link", "BioSNG", "gas"): ("Solid Biomass Primary", "Methane Secondary Input"),
    ("Link", "BioSNG CC", "gas"): ("Solid Biomass Primary", "Methane Secondary Input"),
    ("Link", "Export Foreign", "gas"): ("Methane Secondary Output", "Export"),
    ("Link", "Import Foreign", "gas"): ("Import", "Methane Primary"),
    ("Link", "OCGT", "gas"): ("Methane Secondary Output", "Power Plant"),
    ("Link", "OCGT losses", "gas"): ("Power Plant", "Power Plant Losses"),
    ("Link", "SMR", "gas"): ("Methane Secondary Output", "SMR"),
    ("Link", "SMR CC", "gas"): ("Methane Secondary Output", "SMR"),
    ("Link", "SMR CC losses", "gas"): ("SMR", "SMR Losses"),
    ("Link", "SMR losses", "gas"): ("SMR", "SMR Losses"),
    ("Link", "Sabatier", "gas"): ("Sabatier", "Methane Secondary Input"),
    ("Link", "gas for industry", "gas"): ("", ""),  # todo: drop
    ("Link", "gas for industry CC", "gas"): (
        "Methane Secondary Output",
        "Industry Losses",
    ),
    ("Link", "rural gas boiler", "gas"): ("Methane Secondary Output", "Decentral Heat"),
    ("Link", "rural gas boiler losses", "gas"): (
        "Decentral Heat",
        "Methane Losses",
    ),  # edge case
    ("Link", "urban decentral gas boiler", "gas"): (
        "Methane Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "urban decentral gas boiler losses", "gas"): (
        "Decentral Heat",
        "Methane Losses",
    ),
    ("Link", "CCGT", "gas"): ("Methane Secondary Output", "Power Plant"),
    ("Link", "CCGT losses", "gas"): ("Power Plant", "Power Plant Losses"),
    ("Link", "biogas to gas", "gas"): ("Biogas", "Methane Secondary Input"),
    ("Link", "biogas to gas CC", "gas"): ("Biogas", "Methane Secondary Input"),
    ("Link", "urban central gas CHP", "gas"): ("Methane Secondary Output", "CHP"),
    ("Link", "urban central gas CHP CC", "gas"): ("Methane Secondary Output", "CHP"),
    ("Link", "urban central gas CHP CC losses", "gas"): ("CHP", "CHP Losses"),
    ("Link", "urban central gas CHP losses", "gas"): ("CHP", "CHP Losses"),
    ("Link", "urban central gas boiler", "gas"): (
        "Methane Secondary Output",
        "CHP",
    ),  # fixme: not a CHP
    ("Link", "Export Domestic", "gas"): ("Methane Secondary Output", "Export"),
    ("Link", "Import Domestic", "gas"): ("Import", "Methane Primary"),
    ("Link", "import gas", "gas"): ("Green Global Import", "Methane Primary"),
    ("Store", "gas", "gas"): ("", ""),  # todo: drop
    ("Load", "gas for industry", "gas for industry"): (
        "Methane Secondary Output",
        "Industry",
    ),
    ("Link", "industry methanol", "industry methanol"): ("", ""),  # todo: drop
    ("Load", "industry methanol", "industry methanol"): (
        "Methanol Secondary Output",
        "Industry",
    ),
    ("Link", "kerosene for aviation", "kerosene for aviation"): (
        "Oil Secondary Output",
        "Transport",
    ),  # todo: drop
    ("Link", "methanol-to-kerosene", "kerosene for aviation"): (
        "Methanol Secondary Output",
        "Transport",
    ),
    ("Load", "kerosene for aviation", "kerosene for aviation"): ("", ""),  # todo: drop
    ("Link", "land transport oil", "land transport oil"): (
        "Oil Secondary Output",
        "Transport",
    ),
    ("Load", "land transport oil", "land transport oil"): ("", ""),  # todo: drop
    ("Generator", "lignite", "lignite"): ("Coal Primary", "Coal Secondary Input"),
    ("Link", "lignite", "lignite"): ("Coal Secondary Output", "Power Plant"),
    ("Link", "lignite losses", "lignite"): ("Power Plant", "Power Plant Losses"),
    ("Link", "urban central lignite CHP", "lignite"): ("Coal Secondary Output", "CHP"),
    ("Link", "urban central lignite CHP losses", "lignite"): ("CHP", "CHP Losses"),
    ("Store", "lignite", "lignite"): ("", ""),  # todo: drop
    ("Generator", "solar rooftop", "low voltage"): ("Solar Power", "AC Primary"),
    ("Link", "BEV charger", "low voltage"): ("AC Secondary Output", "Car Battery"),
    ("Link", "BEV charger losses", "low voltage"): (
        "Car Battery",
        "Car Battery Losses",
    ),
    ("Link", "electricity distribution grid", "low voltage"): ("", ""),  # todo: drop
    ("Link", "electricity distribution grid losses", "low voltage"): (
        "",
        "",
    ),  # fixme: wrong sign
    ("Link", "home battery charger", "home battery"): ("", ""),
    ("Link", "home battery discharger", "home battery"): ("", ""),
    ("Link", "home battery discharger losses", "home battery"): (
        "Battery",
        "Battery Losses",
    ),
    ("Store", "home battery", "home battery"): ("", ""),  # todo: drop
    ("Link", "home battery charger", "low voltage"): ("AC Secondary Output", "Battery"),
    ("Link", "home battery charger losses", "low voltage"): ("Battery", "Losses"),
    ("Link", "home battery discharger", "low voltage"): (
        "Battery",
        "AC Secondary Input",
    ),
    ("Link", "rural air heat pump", "low voltage"): (
        "AC Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "rural ground heat pump", "low voltage"): (
        "AC Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "rural resistive heater", "low voltage"): (
        "AC Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "rural resistive heater losses", "low voltage"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "urban decentral air heat pump", "low voltage"): (
        "AC Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "urban decentral resistive heater", "low voltage"): (
        "AC Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "urban decentral resistive heater losses", "low voltage"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "urban central air heat pump", "low voltage"): (
        "AC Secondary Output",
        "Central Heat",
    ),
    ("Link", "urban central ptes heat pump", "low voltage"): (
        "AC Secondary Output",
        "Central Heat",
    ),
    ("Link", "urban central resistive heater", "low voltage"): (
        "AC Secondary Output",
        "Central Heat",
    ),
    ("Link", "urban central resistive heater losses", "low voltage"): (
        "Central Heat",
        "Losses",
    ),
    ("Link", "urban central ptes heat pump losses", "low voltage"): (
        "Central Heat",
        "Central Heat Losses",
    ),
    ("Link", "V2G", "low voltage"): ("", ""),  # todo: drop demand
    ("Load", "agriculture electricity", "low voltage"): (
        "AC Secondary Output",
        "Agriculture",
    ),
    ("Load", "electricity", "low voltage"): ("AC Secondary Output", "AC Base Load"),
    ("Load", "industry electricity", "low voltage"): (
        "AC Secondary Output",
        "Industry",
    ),
    ("Link", "CCGT methanol", "methanol"): ("Methanol Secondary Output", "Power Plant"),
    ("Link", "CCGT methanol CC", "methanol"): (
        "Methanol Secondary Output",
        "Power Plant",
    ),
    ("Link", "CCGT methanol CC losses", "methanol"): (
        "Power Plant",
        "Power Plant Losses",
    ),
    ("Link", "CCGT methanol losses", "methanol"): ("Power Plant", "Power Plant Losses"),
    ("Link", "Methanol steam reforming", "methanol"): (
        "Methanol Secondary Output",
        "Methanol Steam Reforming",
    ),
    ("Link", "Methanol steam reforming CC", "methanol"): (
        "Methanol Secondary Output",
        "Methanol Steam Reforming",
    ),
    ("Link", "Methanol steam reforming CC losses", "methanol"): (
        "Methanol Steam Reforming",
        "Methanol Steam Reforming Losses",
    ),
    ("Link", "Methanol steam reforming losses", "methanol"): (
        "Methanol Steam Reforming",
        "Methanol Steam Reforming Losses",
    ),
    ("Link", "OCGT methanol", "methanol"): ("Methanol Secondary Output", "Power Plant"),
    ("Link", "OCGT methanol losses", "methanol"): ("Power Plant", "Power Plant Losses"),
    ("Link", "allam methanol", "methanol"): (
        "Methanol Secondary Output",
        "Power Plant",
    ),
    ("Link", "allam methanol losses", "methanol"): (
        "Power Plant",
        "Power Plant Losses",
    ),
    ("Link", "biomass-to-methanol", "methanol"): (
        "Methanol Secondary Input",
        "Solid Biomass Secondary Output",
    ),
    ("Link", "biomass-to-methanol CC", "methanol"): (
        "Methanol Secondary Input",
        "Solid Biomass Secondary Output",
    ),
    ("Link", "methanol-to-kerosene", "methanol"): (
        "Methanol Secondary Output",
        "Transport",
    ),
    ("Link", "methanol-to-kerosene losses", "methanol"): (
        "Methanol Secondary Output",
        "Methanol Losses",
    ),
    ("Link", "methanolisation", "methanol"): (
        "Methanol Secondary Input",
        "Methanolisation",
    ),
    ("Link", "industry methanol", "methanol"): (
        "Methanol Secondary Output",
        "Industry",
    ),
    ("Link", "import methanol", "methanol"): (
        "Methanol Green Global Import",
        "Methanol Primary",
    ),
    ("Link", "shipping methanol", "methanol"): (
        "Methanol Secondary Output",
        "Transport",
    ),
    ("Store", "methanol", "methanol"): ("", ""),  # todo: drop
    ("Generator", "municipal solid waste", "municipal solid waste"): (
        "Solid Waste",
        "Waste Primary",
    ),
    ("Link", "Export Foreign", "municipal solid waste"): (
        "Waste Secondary Output",
        "Export",
    ),
    ("Link", "Import Foreign", "municipal solid waste"): (
        "Solid Waste Import",
        "Waste Primary",
    ),
    ("Link", "municipal solid waste", "municipal solid waste"): ("", ""),  # todo: drop
    ("Link", "Export Domestic", "municipal solid waste"): (
        "Waste Secondary Output",
        "Export",
    ),
    ("Link", "Import Domestic", "municipal solid waste"): (
        "Solid Waste Import",
        "Waste Primary",
    ),
    ("Load", "naphtha for industry", "naphtha for industry"): ("", ""),  # todo: drop
    ("Link", "HVC to air", "non-sequestered HVC"): (
        "Waste Secondary Output",
        "HVC to Air",
    ),
    ("Link", "municipal solid waste", "non-sequestered HVC"): (
        "Waste Primary",
        "Waste Secondary Input",
    ),  # todo: drop
    ("Link", "waste CHP", "non-sequestered HVC"): ("Waste Secondary Output", "CHP"),
    ("Link", "waste CHP CC", "non-sequestered HVC"): ("Waste Secondary Output", "CHP"),
    ("Link", "waste CHP CC losses", "non-sequestered HVC"): ("CHP", "CHP Losses"),
    ("Link", "waste CHP losses", "non-sequestered HVC"): ("CHP", "CHP Losses"),
    ("Store", "non-sequestered HVC", "non-sequestered HVC"): (
        "",
        "",
    ),  # todo: drop if zero
    ("Link", "Fischer-Tropsch", "oil"): ("Fischer-Tropsch", "Oil Secondary Input"),
    ("Link", "agriculture machinery oil", "oil"): (
        "Oil Secondary Output",
        "Agriculture",
    ),
    ("Link", "biomass to liquid", "oil"): ("Biomass2Liquid", "Oil Secondary Input"),
    ("Link", "biomass to liquid CC", "oil"): ("Biomass2Liquid", "Oil Secondary Input"),
    ("Link", "electrobiofuels", "oil"): ("Bio Fuels", "Oil Secondary Input"),
    ("Link", "kerosene for aviation", "oil"): ("Oil Secondary Output", "Transport"),
    ("Link", "land transport oil", "oil"): ("Oil Secondary Output", "Transport"),
    ("Link", "naphtha for industry", "oil"): ("Oil Secondary Output", "Industry"),
    ("Link", "oil", "oil"): ("Oil Secondary Output", "Power Plant"),
    ("Link", "oil losses", "oil"): ("Power Plant", "Power Plant Losses"),
    ("Link", "shipping oil", "oil"): ("Oil Secondary Output", "Transport"),
    ("Link", "rural oil boiler", "oil"): ("Oil Secondary Output", "Decentral Heat"),
    ("Link", "rural oil boiler losses", "oil"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "unsustainable bioliquids", "oil"): (
        "Unsustainable Bio Liquids",
        "Oil Secondary Input",
    ),
    ("Link", "urban decentral oil boiler", "oil"): (
        "Oil Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "urban decentral oil boiler losses", "oil"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "urban central oil CHP", "oil"): ("Oil Secondary Output", "Central Heat"),
    ("Link", "urban central oil CHP losses", "oil"): (
        "Central Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "import oil", "oil"): ("Oil Green Global Import", "Oil Primary"),
    ("Link", "oil refining", "oil"): ("Oil Refining", "Oil Secondary Input"),
    ("Store", "oil", "oil"): ("", ""),  # todo: drop
    ("Generator", "oil primary", "oil primary"): ("", ""),  # todo: drop
    ("Link", "oil refining", "oil primary"): ("Raw Oil", "Oil Refining"),
    ("Link", "oil refining losses", "oil primary"): (
        "Oil Refining",
        "Oil Refining Losses",
    ),
    ("Generator", "rural heat vent", "rural heat"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Generator", "rural solar thermal", "rural heat"): (
        "Solar Thermal",
        "Decentral Heat",
    ),
    ("Link", "rural air heat pump", "rural heat"): ("", ""),  # todo: drop
    ("Link", "rural biomass boiler", "rural heat"): ("", ""),  # todo: drop
    ("Link", "rural gas boiler", "rural heat"): ("", ""),  # todo: drop
    ("Link", "rural ground heat pump", "rural heat"): ("", ""),  # todo: drop
    ("Link", "rural resistive heater", "rural heat"): ("", ""),  # todo: drop
    ("Link", "rural oil boiler", "rural heat"): ("", ""),  # todo: drop
    ("Load", "agriculture heat", "rural heat"): ("Decentral Heat", "Agriculture"),
    ("Load", "rural heat", "rural heat"): ("Decentral Heat", "HH & Services"),
    ("Store", "rural water tanks", "rural water tanks"): ("", ""),  # todo: drop
    ("Link", "shipping methanol", "shipping methanol"): (
        "Methanol Secondary Output",
        "Transport",
    ),
    ("Load", "shipping methanol", "shipping methanol"): ("", ""),  # todo: drop
    ("Link", "shipping oil", "shipping oil"): ("Oil Secondary Output", "Transport"),
    ("Load", "shipping oil", "shipping oil"): ("", ""),  # todo: drop
    ("Generator", "unsustainable solid biomass", "solid biomass"): (
        "Solid Biomass",
        "Solid Biomass Primary",
    ),
    ("Generator", "solid biomass", "solid biomass"): (
        "Solid Biomass",
        "Solid Biomass Primary",
    ),
    ("Link", "BioSNG", "solid biomass"): ("", ""),  # todo: drop
    ("Link", "BioSNG CC", "solid biomass"): ("", ""),  # todo: drop
    ("Link", "BioSNG CC losses", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Solid Biomass Secondary Output Losses",
    ),
    ("Link", "BioSNG losses", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Solid Biomass Secondary Output Losses",
    ),
    ("Link", "Export Foreign", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Export",
    ),
    ("Link", "Import Foreign", "solid biomass"): (
        "Import",
        "Solid Biomass Secondary Input",
    ),
    ("Link", "biomass to liquid", "solid biomass"): (
        "Solid Biomass Secondary Input",
        "Biomass2Liquid",
    ),
    ("Link", "biomass to liquid CC", "solid biomass"): (
        "Solid Biomass Secondary Input",
        "Biomass2Liquid",
    ),
    ("Link", "biomass to liquid CC losses", "solid biomass"): (
        "Biomass2Liquid",
        "Biomass2Liquid Losses",
    ),
    ("Link", "biomass to liquid losses", "solid biomass"): (
        "Biomass2Liquid",
        "Biomass2Liquid Losses",
    ),
    ("Link", "biomass-to-methanol", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Methanol Secondary Input",
    ),
    ("Link", "biomass-to-methanol CC", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Methanol Secondary Input",
    ),
    ("Link", "biomass-to-methanol CC losses", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Solid Biomass Secondary Output Losses",
    ),
    ("Link", "biomass-to-methanol losses", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Solid Biomass Secondary Output Losses",
    ),
    ("Link", "electrobiofuels", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Bio Fuels",
    ),
    ("Link", "electrobiofuels losses", "solid biomass"): (
        "Bio Fuels",
        "Bio Fuels Losses",
    ),
    ("Link", "rural biomass boiler", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "rural biomass boiler losses", "solid biomass"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "solid biomass for industry", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Industry",
    ),
    ("Link", "solid biomass for industry CC", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Industry",
    ),
    ("Link", "solid biomass to hydrogen", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "H2 Secondary Input",
    ),
    ("Link", "solid biomass to hydrogen losses", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Solid Biomass Secondary Output Losses",
    ),
    ("Link", "urban decentral biomass boiler", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Decentral Heat",
    ),
    ("Link", "urban decentral biomass boiler losses", "solid biomass"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Link", "urban central solid biomass CHP", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Central Heat",
    ),
    ("Link", "urban central solid biomass CHP CC", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Central Heat",
    ),
    ("Link", "Export Domestic", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Export",
    ),
    ("Link", "Import Domestic", "solid biomass"): (
        "Import",
        "Solid Biomass Secondary Input",
    ),
    ("Link", "solid biomass", "solid biomass"): (
        "Solid Biomass Secondary Output",
        "Power Plant",
    ),
    ("Link", "solid biomass losses", "solid biomass"): (
        "Power Plant",
        "Power Plant Losses",
    ),
    ("Link", "urban central solid biomass CHP losses", "solid biomass"): (
        "CHP",
        "CHP Losses",
    ),
    ("Load", "solid biomass for industry", "solid biomass for industry"): (
        "",
        "",
    ),  # todo: drop
    ("Generator", "unsustainable bioliquids", "unsustainable bioliquids"): (
        "",
        "",
    ),  # todo: drop
    ("Link", "unsustainable bioliquids", "unsustainable bioliquids"): (
        "",
        "",
    ),  # todo: drop
    ("Generator", "uranium", "uranium"): ("", ""),  # todo: drop
    ("Link", "nuclear", "uranium"): ("Uranium", "Nuclear Power"),
    ("Link", "nuclear losses", "uranium"): ("Nuclear Power", "Nuclear Power Losses"),
    ("Store", "uranium", "uranium"): ("", ""),  # todo: drop
    ("Generator", "urban central heat vent", "urban central heat"): (
        "Central Heat",
        "Central Heat Losses",
    ),
    ("Generator", "urban central solar thermal", "urban central heat"): (
        "Solar Thermal",
        "Central Heat",
    ),
    ("Link", "DAC", "urban central heat"): ("Central Heat", "DAC"),
    ("Link", "Fischer-Tropsch", "urban central heat"): (
        "Central Heat",
        "Fischer-Tropsch",
    ),
    ("Link", "H2 Fuel Cell", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "Haber-Bosch", "urban central heat"): ("Haber-Bosch", "Central Heat"),
    ("Link", "Sabatier", "urban central heat"): ("Sabatier", "Central Heat"),
    ("Link", "methanolisation", "urban central heat"): (
        "Methanolisation",
        "Central Heat",
    ),
    ("Link", "urban central air heat pump", "urban central heat"): (
        "",
        "",
    ),  # todo: drop
    ("Link", "urban central coal CHP", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "urban central gas CHP", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "urban central gas CHP CC", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "urban central gas boiler", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "urban central ptes heat pump", "urban central heat"): (
        "",
        "",
    ),  # todo: drop
    ("Link", "urban central resistive heater", "urban central heat"): ("", ""),
    ("Link", "urban central solid biomass CHP", "urban central heat"): (
        "CHP",
        "Central Heat",
    ),
    ("Link", "urban central solid biomass CHP CC", "urban central heat"): (
        "CHP",
        "Central Heat",
    ),
    ("Link", "waste CHP", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "waste CHP CC", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "urban central oil CHP", "urban central heat"): ("CHP", "Central Heat"),
    ("Link", "urban central lignite CHP", "urban central heat"): (
        "CHP",
        "Central Heat",
    ),
    ("Link", "H2 Electrolysis", "urban central heat"): ("Electrolysis", "Central Heat"),
    ("Link", "urban central H2 CHP", "urban central heat"): ("CHP", "Central Heat"),
    ("Load", "low-temperature heat for industry", "urban central heat"): (
        "Central Heat",
        "Industry",
    ),
    ("Load", "urban central heat", "urban central heat"): (
        "Central Heat",
        "HH & Services",
    ),
    ("Store", "urban central water pits", "urban central water pits"): (
        "",
        "",
    ),  # todo: drop
    ("Store", "urban central water tanks", "urban central water tanks"): (
        "",
        "",
    ),  # todo: drop
    ("Generator", "urban decentral heat vent", "urban decentral heat"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    ("Generator", "urban decentral solar thermal", "urban decentral heat"): (
        "Solar Thermal",
        "Decentral Heat",
    ),
    ("Link", "DAC", "urban decentral heat"): ("Decentral Heat", "DAC"),
    ("Link", "urban decentral air heat pump", "urban decentral heat"): (
        "",
        "",
    ),  # todo: reverse
    ("Link", "urban decentral biomass boiler", "urban decentral heat"): (
        "",
        "",
    ),  # todo: reverse
    ("Link", "urban decentral gas boiler", "urban decentral heat"): (
        "",
        "",
    ),  # todo: reverse
    ("Link", "urban decentral resistive heater", "urban decentral heat"): (
        "AC Secondary Output",
        "Decentral Heat",
    ),  # this one is correct, others need updates
    ("Link", "urban decentral oil boiler", "urban decentral heat"): (
        "Oil Secondary Output",
        "Decentral Heat",
    ),
    ("Load", "low-temperature heat for industry", "urban decentral heat"): (
        "",
        "",
    ),  # todo: drop
    ("Load", "urban decentral heat", "urban decentral heat"): (
        "Decentral Heat",
        "HH & Services",
    ),
    ("Store", "urban decentral water tanks", "urban decentral water tanks"): (
        "",
        "",
    ),  # todo: drop
}


LINK_MAPPING = {
    ("Link", "BEV charger", "EV battery"): ("", ""),
    ("Link", "BEV charger", "low voltage"): ("", ""),
    ("Link", "BEV charger", "low voltage losses"): ("", ""),
    ("Link", "BioSNG", "gas"): ("", ""),
    ("Link", "BioSNG", "solid biomass"): ("", ""),
    ("Link", "BioSNG", "solid biomass losses"): ("", ""),
    ("Link", "BioSNG CC", "gas"): ("", ""),
    ("Link", "BioSNG CC", "solid biomass"): ("", ""),
    ("Link", "BioSNG CC", "solid biomass losses"): ("", ""),
    ("Link", "CCGT", "AC"): ("", ""),
    ("Link", "CCGT", "gas"): ("", ""),
    ("Link", "CCGT", "gas losses"): ("", ""),
    ("Link", "CCGT methanol", "AC"): ("", ""),
    ("Link", "CCGT methanol", "methanol"): ("", ""),
    ("Link", "CCGT methanol", "methanol losses"): ("", ""),
    ("Link", "CCGT methanol CC", "AC"): ("", ""),
    ("Link", "CCGT methanol CC", "methanol"): ("", ""),
    ("Link", "CCGT methanol CC", "methanol losses"): ("", ""),
    ("Link", "DAC", "AC"): ("", ""),
    ("Link", "DAC", "urban central heat"): ("", ""),
    ("Link", "DAC", "urban decentral heat"): ("", ""),
    ("Store", "EV battery", "EV battery"): ("", ""),
    ("Link", "Fischer-Tropsch", "H2"): ("", ""),
    ("Link", "Fischer-Tropsch", "H2 losses"): ("", ""),
    ("Link", "Fischer-Tropsch", "oil"): ("", ""),
    ("Link", "Fischer-Tropsch", "urban central heat"): ("", ""),
    ("Link", "H2 Electrolysis", "AC losses"): (
        "Secondary AC In",
        "Secondary AC Losses",
    ),
    ("Link", "H2 Electrolysis", "H2"): ("Secondary AC In", "Secondary H2 Out"),
    ("Link", "H2 Electrolysis", "urban central heat"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    ("Link", "H2 Fuel Cell", "AC"): ("", ""),
    ("Link", "H2 Fuel Cell", "H2"): ("", ""),
    ("Link", "H2 Fuel Cell", "H2 losses"): ("", ""),
    ("Link", "H2 Fuel Cell", "urban central heat"): ("", ""),
    ("Link", "H2 OCGT", "AC"): ("", ""),
    ("Link", "H2 OCGT", "H2"): ("", ""),
    ("Link", "H2 OCGT", "H2 losses"): ("", ""),
    ("Store", "H2 Store", "H2"): ("", ""),
    ("Load", "H2 for industry", "H2"): ("", ""),
    ("Link", "HVC to air", "non-sequestered HVC"): ("", ""),
    ("Link", "Haber-Bosch", "AC"): ("", ""),
    ("Link", "Haber-Bosch", "AC losses"): ("", ""),
    ("Link", "Haber-Bosch", "H2"): ("", ""),
    ("Link", "Haber-Bosch", "H2 losses"): ("", ""),
    ("Link", "Haber-Bosch", "NH3"): ("", ""),
    ("Link", "Haber-Bosch", "urban central heat"): ("", ""),
    ("Link", "Methanol steam reforming", "H2"): ("", ""),
    ("Link", "Methanol steam reforming", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming", "methanol losses"): ("", ""),
    ("Link", "Methanol steam reforming CC", "H2"): ("", ""),
    ("Link", "Methanol steam reforming CC", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming CC", "methanol losses"): ("", ""),
    ("Load", "NH3", "NH3"): ("", ""),
    ("Link", "OCGT", "AC"): ("", ""),
    ("Link", "OCGT", "gas"): ("", ""),
    ("Link", "OCGT", "gas losses"): ("", ""),
    ("Link", "OCGT methanol", "AC"): ("", ""),
    ("Link", "OCGT methanol", "methanol"): ("", ""),
    ("Link", "OCGT methanol", "methanol losses"): ("", ""),
    ("StorageUnit", "PHS demand", "AC"): ("", ""),
    ("StorageUnit", "PHS supply", "AC"): ("", ""),
    ("Link", "SMR", "H2"): ("", ""),
    ("Link", "SMR", "gas"): ("", ""),
    ("Link", "SMR", "gas losses"): ("", ""),
    ("Link", "SMR CC", "H2"): ("", ""),
    ("Link", "SMR CC", "gas"): ("", ""),
    ("Link", "SMR CC", "gas losses"): ("", ""),
    ("Link", "Sabatier", "H2"): ("", ""),
    ("Link", "Sabatier", "H2 losses"): ("", ""),
    ("Link", "Sabatier", "gas"): ("", ""),
    ("Link", "Sabatier", "urban central heat"): ("", ""),
    ("Link", "V2G", "EV battery"): ("", ""),
    ("Link", "V2G", "EV battery losses"): ("", ""),
    ("Link", "V2G", "low voltage"): ("", ""),
    ("Load", "agriculture electricity", "low voltage"): ("", ""),
    ("Load", "agriculture heat", "rural heat"): ("", ""),
    ("Link", "agriculture machinery oil", "agriculture machinery oil"): ("", ""),
    ("Link", "agriculture machinery oil", "oil"): ("", ""),
    ("Load", "agriculture machinery oil", "agriculture machinery oil"): ("", ""),
    ("Link", "allam methanol", "AC"): ("", ""),
    ("Link", "allam methanol", "methanol"): ("", ""),
    ("Link", "allam methanol", "methanol losses"): ("", ""),
    ("Link", "ammonia cracker", "H2"): ("", ""),
    ("Link", "ammonia cracker", "NH3"): ("", ""),
    ("Link", "ammonia cracker", "NH3 losses"): ("", ""),
    ("Store", "ammonia store", "NH3"): ("", ""),
    ("Store", "battery", "battery"): ("", ""),
    ("Link", "battery charger", "AC"): ("", ""),
    ("Link", "battery charger", "AC losses"): ("", ""),
    ("Link", "battery charger", "battery"): ("", ""),
    ("Link", "battery discharger", "AC"): ("", ""),
    ("Link", "battery discharger", "battery"): ("", ""),
    ("Link", "battery discharger", "battery losses"): ("", ""),
    ("Generator", "biogas", "biogas"): ("", ""),
    ("Link", "biogas to gas", "biogas"): ("", ""),
    ("Link", "biogas to gas", "gas"): ("", ""),
    ("Link", "biogas to gas CC", "biogas"): ("", ""),
    ("Link", "biogas to gas CC", "gas"): ("", ""),
    ("Link", "biomass to liquid", "oil"): ("", ""),
    ("Link", "biomass to liquid", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid", "solid biomass losses"): ("", ""),
    ("Link", "biomass to liquid CC", "oil"): ("", ""),
    ("Link", "biomass to liquid CC", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid CC", "solid biomass losses"): ("", ""),
    ("Link", "biomass-to-methanol", "methanol"): ("", ""),
    ("Link", "biomass-to-methanol", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol", "solid biomass losses"): ("", ""),
    ("Link", "biomass-to-methanol CC", "methanol"): ("", ""),
    ("Link", "biomass-to-methanol CC", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol CC", "solid biomass losses"): ("", ""),
    ("Generator", "coal", "coal"): ("", ""),
    ("Link", "coal", "AC"): ("", ""),
    ("Link", "coal", "coal"): ("", ""),
    ("Link", "coal", "coal losses"): ("", ""),
    ("Store", "coal", "coal"): ("", ""),
    ("Load", "electricity", "low voltage"): ("Secondary AC Out", "Final AC"),
    ("Link", "electricity distribution grid", "losses"): ("", ""),
    ("Link", "electrobiofuels", "H2"): ("", ""),
    ("Link", "electrobiofuels", "H2 losses"): ("", ""),
    ("Link", "electrobiofuels", "oil"): ("", ""),
    ("Link", "electrobiofuels", "solid biomass"): ("", ""),
    ("Link", "electrobiofuels", "solid biomass losses"): ("", ""),
    ("Store", "gas", "gas"): ("", ""),
    ("Link", "gas for industry", "gas"): ("", ""),
    ("Link", "gas for industry", "gas for industry"): ("", ""),
    ("Load", "gas for industry", "gas for industry"): ("", ""),
    ("Link", "gas for industry CC", "gas"): ("", ""),
    ("Link", "gas for industry CC", "gas for industry"): ("", ""),
    ("Link", "gas for industry CC", "gas losses"): ("", ""),
    ("Store", "home battery", "home battery"): ("", ""),
    ("Link", "home battery charger", "home battery"): ("", ""),
    ("Link", "home battery charger", "low voltage"): ("", ""),
    ("Link", "home battery charger", "low voltage losses"): ("", ""),
    ("Link", "home battery discharger", "home battery"): ("", ""),
    ("Link", "home battery discharger", "home battery losses"): ("", ""),
    ("Link", "home battery discharger", "low voltage"): ("", ""),
    ("StorageUnit", "hydro supply", "AC"): ("", ""),  # is primary
    ("Generator", "import H2", "H2"): ("", ""),
    ("Generator", "import NH3", "NH3"): ("", ""),
    ("Link", "import gas", "gas"): ("", ""),
    ("Link", "import methanol", "methanol"): ("", ""),
    ("Link", "import oil", "oil"): ("", ""),
    ("Load", "industry electricity", "low voltage"): ("", ""),
    ("Link", "industry methanol", "industry methanol"): ("", ""),
    ("Link", "industry methanol", "methanol"): ("", ""),
    ("Load", "industry methanol", "industry methanol"): ("", ""),
    ("Link", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Link", "kerosene for aviation", "oil"): ("", ""),
    ("Load", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Load", "land transport EV", "EV battery"): ("", ""),
    ("Load", "land transport fuel cell", "H2"): ("", ""),
    ("Generator", "lignite", "lignite"): ("", ""),
    ("Link", "lignite", "AC"): ("", ""),
    ("Link", "lignite", "lignite"): ("", ""),
    ("Link", "lignite", "lignite losses"): ("", ""),
    ("Store", "lignite", "lignite"): ("", ""),
    ("Generator", "lng gas", "gas"): ("", ""),
    ("Load", "low-temperature heat for industry", "urban central heat"): ("", ""),
    ("Load", "low-temperature heat for industry", "urban decentral heat"): ("", ""),
    ("Store", "methanol", "methanol"): ("", ""),
    ("Link", "methanol-to-kerosene", "H2"): ("", ""),
    ("Link", "methanol-to-kerosene", "H2 losses"): ("", ""),
    ("Link", "methanol-to-kerosene", "kerosene for aviation"): ("", ""),
    ("Link", "methanol-to-kerosene", "methanol"): ("", ""),
    ("Link", "methanol-to-kerosene", "methanol losses"): ("", ""),
    ("Link", "methanolisation", "AC"): ("", ""),
    ("Link", "methanolisation", "AC losses"): ("", ""),
    ("Link", "methanolisation", "H2"): ("", ""),
    ("Link", "methanolisation", "H2 losses"): ("", ""),
    ("Link", "methanolisation", "methanol"): ("", ""),
    ("Link", "methanolisation", "urban central heat"): ("", ""),
    ("Generator", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Link", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Link", "municipal solid waste", "non-sequestered HVC"): ("", ""),
    ("Link", "naphtha for industry", "naphtha for industry"): ("", ""),
    ("Link", "naphtha for industry", "oil"): ("", ""),
    ("Load", "naphtha for industry", "naphtha for industry"): ("", ""),
    ("Store", "non-sequestered HVC", "non-sequestered HVC"): ("", ""),
    ("Link", "nuclear", "AC"): ("Primary Uranium", "Secondary AC In"),
    ("Link", "nuclear", "uranium losses"): ("Primary Uranium", "Secondary Losses"),
    ("Generator", "offwind-ac", "AC"): ("Wind Power", "Primary AC"),
    ("Generator", "offwind-dc", "AC"): ("Wind Power", "Primary AC"),
    ("Link", "oil", "AC"): ("", ""),
    ("Link", "oil", "oil"): ("", ""),
    ("Link", "oil", "oil losses"): ("", ""),
    ("Store", "oil", "oil"): ("", ""),
    ("Generator", "oil primary", "oil primary"): ("", ""),
    ("Link", "oil refining", "oil"): ("", ""),
    ("Link", "oil refining", "oil primary"): ("", ""),
    ("Link", "oil refining", "oil primary losses"): ("", ""),
    ("Generator", "onwind", "AC"): ("Wind Power", "Primary AC"),
    ("Generator", "pipeline gas", "gas"): ("", ""),
    ("Generator", "production gas", "gas"): ("", ""),
    ("Generator", "ror", "AC"): ("Hydro Power", "Primary AC"),
    ("Link", "rural air heat pump", "ambient heat"): ("", ""),
    ("Link", "rural air heat pump", "low voltage"): ("", ""),
    ("Link", "rural air heat pump", "rural heat"): ("", ""),
    ("Link", "rural biomass boiler", "rural heat"): ("", ""),
    ("Link", "rural biomass boiler", "solid biomass"): ("", ""),
    ("Link", "rural biomass boiler", "solid biomass losses"): ("", ""),
    ("Link", "rural gas boiler", "gas"): ("", ""),
    ("Link", "rural gas boiler", "gas losses"): ("", ""),
    ("Link", "rural gas boiler", "rural heat"): ("", ""),
    ("Link", "rural ground heat pump", "ambient heat"): ("", ""),
    ("Link", "rural ground heat pump", "low voltage"): ("", ""),
    ("Link", "rural ground heat pump", "rural heat"): ("", ""),
    ("Load", "rural heat", "rural heat"): ("", ""),
    ("Generator", "rural heat vent", "rural heat"): ("", ""),
    ("Link", "rural resistive heater", "low voltage"): ("", ""),
    ("Link", "rural resistive heater", "low voltage losses"): ("", ""),
    ("Link", "rural resistive heater", "rural heat"): ("", ""),
    ("Generator", "rural solar thermal", "rural heat"): (
        "Solar Heat",
        "HH & Services",
    ),  # to FED
    ("Store", "rural water tanks", "rural water tanks"): ("", ""),
    ("Link", "rural water tanks charger", "rural heat"): ("", ""),
    ("Link", "rural water tanks charger", "rural water tanks"): ("", ""),
    ("Link", "rural water tanks discharger", "rural heat"): ("", ""),
    ("Link", "rural water tanks discharger", "rural water tanks"): ("", ""),
    ("Link", "shipping methanol", "methanol"): ("", ""),
    ("Link", "shipping methanol", "shipping methanol"): ("", ""),
    ("Load", "shipping methanol", "shipping methanol"): ("", ""),
    ("Generator", "solar", "AC"): ("Solar Power", "Primary AC"),
    ("Generator", "solar rooftop", "low voltage"): ("Solar Power", "Primary AC"),
    ("Generator", "solar-hsat", "AC"): ("Solar Power", "Primary AC"),
    ("Generator", "solid biomass", "solid biomass"): (
        "Solid Biomass",
        "Primary Solids",
    ),
    ("Link", "solid biomass for industry", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry", "solid biomass for industry"): ("", ""),
    ("Load", "solid biomass for industry", "solid biomass for industry"): ("", ""),
    ("Link", "solid biomass for industry CC", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry CC", "solid biomass for industry"): ("", ""),
    ("Link", "solid biomass for industry CC", "solid biomass losses"): ("", ""),
    ("Link", "solid biomass to hydrogen", "H2"): ("", ""),
    ("Link", "solid biomass to hydrogen", "solid biomass"): ("", ""),
    ("Link", "solid biomass to hydrogen", "solid biomass losses"): ("", ""),
    ("Generator", "uranium", "uranium"): ("Uranium", "Primary Uranium"),
    ("Store", "uranium", "uranium"): ("", ""),
    ("Link", "urban central H2 CHP", "AC"): ("Secondary H2 In", "Secondary AC Out"),
    ("Link", "urban central H2 CHP", "H2 losses"): (
        "Secondary H2 In",
        "Secondary H2 Losses",
    ),
    ("Link", "urban central H2 CHP", "urban central heat"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central air heat pump", "ambient heat"): ("", ""),
    ("Link", "urban central air heat pump", "low voltage"): ("", ""),
    ("Link", "urban central air heat pump", "urban central heat"): ("", ""),
    ("Link", "urban central coal CHP", "AC"): ("Secondary coal In", "Secondary AC Out"),
    ("Link", "urban central coal CHP", "ambient heat"): (
        "Secondary Ambient Heat",
        "Secondary Heat Out",
    ),
    ("Link", "urban central coal CHP", "urban central heat"): (
        "Secondary coal In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas CHP", "AC"): ("Secondary gas In", "Secondary AC Out"),
    ("Link", "urban central gas CHP", "ambient heat"): (
        "Secondary Ambient Heat",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas CHP", "gas losses"): (
        "Secondary gas In",
        "Secondary gas Losses",
    ),
    ("Link", "urban central gas CHP", "urban central heat"): (
        "Secondary gas In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas CHP CC", "AC"): (
        "Secondary gas In",
        "Secondary AC Out",
    ),
    ("Link", "urban central gas CHP CC", "gas losses"): (
        "Secondary gas In",
        "Secondary gas Losses",
    ),
    ("Link", "urban central gas CHP CC", "urban central heat"): (
        "Secondary gas In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas boiler", "ambient heat"): ("", ""),
    ("Link", "urban central gas boiler", "gas"): ("", ""),
    ("Link", "urban central gas boiler", "urban central heat"): ("", ""),
    ("Load", "urban central heat", "urban central heat"): ("", ""),
    ("Generator", "urban central heat vent", "urban central heat"): ("", ""),
    ("Link", "urban central ptes heat pump", "ambient heat"): ("", ""),
    ("Link", "urban central ptes heat pump", "low voltage"): ("", ""),
    ("Link", "urban central ptes heat pump", "urban central heat"): ("", ""),
    ("Link", "urban central resistive heater", "low voltage"): ("", ""),
    ("Link", "urban central resistive heater", "low voltage losses"): ("", ""),
    ("Link", "urban central resistive heater", "urban central heat"): ("", ""),
    ("Generator", "urban central solar thermal", "urban central heat"): (
        "Solar Heat",
        "Primary Heat",
    ),
    ("Link", "urban central solid biomass CHP", "AC"): ("", ""),
    ("Link", "urban central solid biomass CHP", "ambient heat"): ("", ""),
    ("Link", "urban central solid biomass CHP", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP", "urban central heat"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "AC"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "ambient heat"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "urban central heat"): ("", ""),
    ("Store", "urban central water pits", "urban central water pits"): ("", ""),
    ("Link", "urban central water pits charger", "urban central heat"): ("", ""),
    ("Link", "urban central water pits charger", "urban central water pits"): ("", ""),
    ("Link", "urban central water pits discharger", "urban central heat"): ("", ""),
    ("Link", "urban central water pits discharger", "urban central water pits"): (
        "",
        "",
    ),
    ("Store", "urban central water tanks", "urban central water tanks"): ("", ""),
    ("Link", "urban central water tanks charger", "urban central heat"): ("", ""),
    ("Link", "urban central water tanks charger", "urban central water tanks"): (
        "",
        "",
    ),
    ("Link", "urban central water tanks discharger", "urban central heat"): ("", ""),
    ("Link", "urban central water tanks discharger", "urban central water tanks"): (
        "",
        "",
    ),
    ("Link", "urban decentral air heat pump", "ambient heat"): ("", ""),
    ("Link", "urban decentral air heat pump", "low voltage"): ("", ""),
    ("Link", "urban decentral air heat pump", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral biomass boiler", "solid biomass"): ("", ""),
    ("Link", "urban decentral biomass boiler", "solid biomass losses"): ("", ""),
    ("Link", "urban decentral biomass boiler", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral gas boiler", "gas"): ("", ""),
    ("Link", "urban decentral gas boiler", "gas losses"): ("", ""),
    ("Link", "urban decentral gas boiler", "urban decentral heat"): ("", ""),
    ("Load", "urban decentral heat", "urban decentral heat"): ("", ""),
    ("Generator", "urban decentral heat vent", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral resistive heater", "low voltage"): ("", ""),
    ("Link", "urban decentral resistive heater", "low voltage losses"): ("", ""),
    ("Link", "urban decentral resistive heater", "urban decentral heat"): ("", ""),
    ("Generator", "urban decentral solar thermal", "urban decentral heat"): ("", ""),
    ("Store", "urban decentral water tanks", "urban decentral water tanks"): ("", ""),
    ("Link", "urban decentral water tanks charger", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral water tanks charger", "urban decentral water tanks"): (
        "",
        "",
    ),
    ("Link", "urban decentral water tanks discharger", "urban decentral heat"): (
        "",
        "",
    ),
    ("Link", "urban decentral water tanks discharger", "urban decentral water tanks"): (
        "",
        "",
    ),
    ("Link", "waste CHP", "AC"): ("", ""),
    ("Link", "waste CHP", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP", "non-sequestered HVC losses"): ("", ""),
    ("Link", "waste CHP", "urban central heat"): ("", ""),
    ("Link", "waste CHP CC", "AC"): ("", ""),
    ("Link", "waste CHP CC", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP CC", "non-sequestered HVC losses"): ("", ""),
    ("Link", "waste CHP CC", "urban central heat"): ("", ""),
}


class SankeyChart(ESMChart):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.location = self._df.index.unique(DM.LOCATION).item()
        self.year = self._df.index.unique(DM.YEAR).item()
        self._df = self._df.droplevel(DM.YEAR).droplevel(DM.LOCATION)
        self._df.columns = ["value"]
        # self._df = self._df.abs()

    @staticmethod
    def _add_source_target_columns(idx: tuple) -> pd.Series:
        return pd.Series(LINK_MAPPING.get(idx, ""))

    def add_source_target_columns(self):
        _df = self._df.copy()
        _df["index"] = _df.index  # convert to tuples
        _df[["source", "target"]] = _df["index"].apply(self._add_source_target_columns)
        _df = _df.drop(columns=["index"])
        return _df.query("source != '' and target != ''")

    @staticmethod
    def get_label_mapping(df_agg):
        return {
            label: i
            for i, label in enumerate(
                set(pd.concat([df_agg["source"], df_agg["target"]]))
            )
        }

    def add_customdata(self, df_agg):
        to_concat = []
        for _, data in df_agg.groupby(["source", "target"]):
            data = data.reset_index()
            carrier_values = [
                f"{c}: {prettify_number(v)} {self.unit}"
                for c, v in zip(data["carrier"], data["value"])
            ]
            data["link_customdata"] = "<br>".join(carrier_values)
            to_concat.append(data)

        return pd.concat(to_concat)

    def add_id_source_target_columns(self, df_agg):
        label_mapping = self.get_label_mapping(df_agg)
        df_agg["source_id"] = df_agg["source"].map(label_mapping)
        df_agg["target_id"] = df_agg["target"].map(label_mapping)
        return df_agg

    @staticmethod
    def map_colors_from_bus_carrier(df_agg):
        df_agg["color"] = (
            df_agg["bus_carrier"].map(BUS_CARRIER_COLORS).fillna(COLOUR.grey_neutral)
        )
        return df_agg

    @staticmethod
    def combine_duplicates(df_agg):
        return (
            df_agg.groupby(["source", "target"])
            .agg(
                {
                    "value": "sum",
                    "source": "first",
                    "target": "first",
                    "source_id": "first",
                    "target_id": "first",
                    "link_customdata": "first",
                    "bus_carrier": "first",
                }
            )
            .reset_index(drop=True)
        )

    @staticmethod
    def get_label_group(lbl: str) -> str:
        if lbl in (
            "Hydro Power",
            "Solar Power",
            "Wind Power",
            "Solid Biomass",
            "Solar Heat",
            "Uranium",
        ):
            return "A10"
        elif lbl in ("Nuclear Power Plant",):
            return "A15"
        elif lbl.startswith("Primary"):
            return "A20"
        elif lbl.startswith("Primary") and lbl.endswith("Losses"):
            return "A25"
        elif lbl.startswith("Secondary") and lbl.endswith("In"):
            return "B10"
        elif lbl.startswith("Secondary") and lbl.endswith(("Out", "Losses")):
            return "B20"
        else:
            raise ValueError(f"Unknown label group: '{lbl}'")

    def plot(self):
        # Concatenate the data with source and target columns
        df_agg = self.add_source_target_columns()

        # add jumpers:
        for bus_carrier in ("AC",):
            primary_to_secondary = filter_by(df_agg, target=f"Primary {bus_carrier}")
            row_idx = ("Jumper", "primary to secondary", bus_carrier)
            df_agg.loc[row_idx, ["value", "source", "target"]] = [
                primary_to_secondary.value.sum(),
                f"Primary {bus_carrier}",
                f"Secondary {bus_carrier} In",
            ]

            secondary_demand = filter_by(df_agg, source=f"Secondary {bus_carrier} In")[
                "value"
            ].sum()
            primary_supply = filter_by(df_agg, target=f"Secondary {bus_carrier} In")[
                "value"
            ].sum()
            row_idx = ("Jumper", "secondary forwarding", bus_carrier)
            df_agg.loc[row_idx, ["value", "source", "target"]] = [
                primary_supply - secondary_demand,
                f"Secondary {bus_carrier} In",
                f"Secondary {bus_carrier} Out",
            ]

        df_agg = self.add_id_source_target_columns(df_agg)
        df_agg = self.add_customdata(df_agg)
        df_agg = self.combine_duplicates(df_agg)
        df_agg = self.map_colors_from_bus_carrier(df_agg)

        labels = pd.Series(list(self.get_label_mapping(df_agg)))
        # groups = [self.get_label_group(lbl) for lbl in labels.values]

        # # x pos is wrong
        # x_pos, y_pos = self.sankey_positions_auto_order_optimized(
        #     labels,
        #     groups,
        #     df_agg.source_id,
        #     df_agg.target_id,
        #     x_spacing=0.4,
        #     y_padding=0.15,
        # )
        # x_pos = [x / max(x_pos) for x in x_pos]

        self.fig = Figure(
            data=[
                Sankey(
                    valuesuffix=self.unit,
                    node=dict(
                        line=dict(color="black", width=0.5),
                        label=labels,
                        hovertemplate="%{label}: %{value}<extra></extra>",
                        # x=x_pos,
                        # y=y_pos,
                        pad=20,
                        thickness=20,
                        # color=df_agg.color,
                        # groups=[[1, 2], [3, 4]],
                    ),
                    link=dict(
                        # arrowlen=15,
                        source=df_agg.source_id,
                        target=df_agg.target_id,
                        value=df_agg.value.abs(),
                        color=df_agg.color,
                        customdata=df_agg.link_customdata,
                        hovertemplate="%{customdata} <extra></extra>",
                    ),
                )
            ]
        )

        self._set_base_layout()
        # self._style_title_and_legend_and_xaxis_label()
        # self._append_footnotes()

        # plotly.io.show(self.fig)  # todo: remove debugging

    def _set_base_layout(self):
        """Set various figure properties."""
        self.fig.update_layout(
            height=800,
            font_family="Calibri",
            plot_bgcolor="#ffffff",
            legend_title_text=self.cfg.legend_header,
        )
        # update axes
        self.fig.update_yaxes(
            showgrid=self.cfg.yaxes_showgrid, visible=self.cfg.yaxes_visible
        )
        self.fig.update_layout(
            xaxis={"categoryorder": "category ascending"},
            hovermode="x",  # all categories are shown by mouse-over
        )
        # trace order always needs to be reversed to show correct order
        # of legend entries for relative bar charts
        self.fig.update_layout(legend={"traceorder": "reversed"})

        # export the metadata directly in the Layout property for JSON
        self.fig.update_layout(meta=[RUN_META_DATA])

    @staticmethod
    def sankey_positions_auto_order_optimized(
        labels, groups, sources, targets, x_spacing=0.3, y_padding=0.05, iterations=10
    ):
        """
        Generate x/y positions for Sankey nodes grouped vertically,
        ordered left-to-right based on flows, and vertically arranged to reduce link crossings.

        Parameters
        ----------
        labels : list of str
            Node labels
        groups : list of str
            Group name for each node
        sources : list of int
            Source node indices
        targets : list of int
            Target node indices
        x_spacing : float
            Horizontal distance between groups
        y_padding : float
            Vertical distance between nodes in the same group
        iterations : int
            Number of vertical reordering passes to reduce link crossings

        Returns
        -------
        x_positions, y_positions : lists of floats
            Coordinates for Plotly Sankey `node.x` and `node.y`
        """
        # Map node index -> group
        node_to_group = {i: groups[i] for i in range(len(labels))}

        # Build group adjacency
        group_graph = defaultdict(set)
        for s, t in zip(sources, targets):
            g_s = node_to_group[s]
            g_t = node_to_group[t]
            if g_s != g_t:
                group_graph[g_s].add(g_t)

        # Topological sort of groups
        indegree = {g: 0 for g in set(groups)}
        for g in group_graph:
            for neigh in group_graph[g]:
                indegree[neigh] += 1

        queue = [g for g in indegree if indegree[g] == 0]
        ordered_groups = []
        while queue:
            g = queue.pop(0)
            ordered_groups.append(g)
            for neigh in sorted(group_graph[g]):
                indegree[neigh] -= 1
                if indegree[neigh] == 0:
                    queue.append(neigh)
        for g in set(groups):
            if g not in ordered_groups:
                ordered_groups.append(g)

        # Assign x positions
        group_x_map = {g: i * x_spacing for i, g in enumerate(ordered_groups)}

        # Start with alphabetical order within each group
        group_nodes = {
            g: sorted(
                [i for i, grp in enumerate(groups) if grp == g], key=lambda i: labels[i]
            )
            for g in ordered_groups
        }

        # Optimization iterations
        for _ in range(iterations):
            changed = False
            for g_idx, g in enumerate(ordered_groups):
                nodes = group_nodes[g]
                if g_idx > 0:
                    # Sort by avg y of connected nodes in the previous group
                    prev_g = ordered_groups[g_idx - 1]
                    neighbor_map = defaultdict(list)
                    for s, t in zip(sources, targets):
                        if node_to_group[t] == g and node_to_group[s] == prev_g:
                            neighbor_map[t].append(s)
                    if neighbor_map:
                        order = sorted(
                            nodes,
                            key=lambda n: sum(
                                group_nodes[prev_g].index(nb)
                                for nb in neighbor_map.get(n, [])
                            )
                            / (len(neighbor_map.get(n, [])) or 1),
                        )
                        if order != nodes:
                            group_nodes[g] = order
                            changed = True

                if g_idx < len(ordered_groups) - 1:
                    # Sort by avg y of connected nodes in next group
                    next_g = ordered_groups[g_idx + 1]
                    neighbor_map = defaultdict(list)
                    for s, t in zip(sources, targets):
                        if node_to_group[s] == g and node_to_group[t] == next_g:
                            neighbor_map[s].append(t)
                    if neighbor_map:
                        order = sorted(
                            group_nodes[g],
                            key=lambda n: sum(
                                group_nodes[next_g].index(nb)
                                for nb in neighbor_map.get(n, [])
                            )
                            / (len(neighbor_map.get(n, [])) or 1),
                        )
                        if order != group_nodes[g]:
                            group_nodes[g] = order
                            changed = True
            if not changed:
                break

        # Assign final positions
        x_positions = [None] * len(labels)
        y_positions = [None] * len(labels)

        for g in ordered_groups:
            indices = group_nodes[g]
            n = len(indices)
            if n == 1:
                y_coords = [0.5]
            else:
                total_height = (n - 1) * y_padding
                start_y = 0.5 - total_height / 2
                y_coords = [start_y + j * y_padding for j in range(n)]
            for idx, y in zip(indices, y_coords):
                x_positions[idx] = group_x_map[g]
                y_positions[idx] = y

        return x_positions, y_positions


@dataclasses.dataclass
class SankeyNode:
    name: str
    x: float
    y: float
    color: str
    # def __init__(self, bus_carrier, label, variables):
    #     self.bus_carrier: str = bus_carrier
    #     self.label = label
    #     self.variables = variables
    #     self.color = COLOUR_SCHEME[bus_carrier]


def sankey(df: pyam.IamDataFrame, mapping: dict) -> Figure:
    """
    Plot a sankey diagram.

    It is currently only possible to create this diagram for single years.

    Parameters
    ----------
    df
        Data to be plotted
    mapping
        Assigns the source and target component of a variable

    Returns
    -------
    :
        The generated plotly figure.
    """

    # Check for duplicates
    for col in [name for name in df.dimensions if name != "variable"]:
        levels = get_index_levels(df._data, col)
        if len(levels) > 1:
            raise ValueError(f"Non-unique values in column {col}: {levels}")

    # Concatenate the data with source and target columns
    _df = pd.DataFrame.from_dict(
        mapping, orient="index", columns=["source", "target"]
    ).merge(df._data, how="left", left_index=True, right_on="variable")
    label_mapping = {
        label: i
        for i, label in enumerate(set(pd.concat([_df["source"], _df["target"]])))
    }
    _df = _df.replace(label_mapping)

    def get_carrier_color(s) -> str:
        carrier = re.findall(r"\|AC", s)[0].strip("|")
        color_map = {"AC": COLOUR.red}
        return color_map[carrier]

    _df["color"] = _df.index.get_level_values("variable").map(get_carrier_color)

    region = get_index_levels(_df, "region")[0]
    unit = " " + get_index_levels(_df, "unit")[0]
    year = get_index_levels(_df, "year")[0]
    fig = Figure(
        data=[
            Sankey(
                valuesuffix=unit,
                node=dict(
                    # pad=15,
                    # thickness=10,
                    line=dict(color="black", width=0.5),
                    label=pd.Series(list(label_mapping)),
                    hovertemplate="%{label}: %{value}<extra></extra>",
                    color=_df.color,
                ),
                link=dict(
                    # arrowlen=15,
                    source=_df.source,
                    target=_df.target,
                    value=_df.value,
                    color=_df.color,
                    hovertemplate='"%{source.label}" to "%{target.label}": %{value}<extra></extra> ',
                ),
            )
        ]
    )
    fig.update_layout(title_text=f"region: {region}, year: {year}", font_size=10)
    return fig


def read_iamc_data_frame(filepath):
    xls = pd.read_excel(
        filepath,
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    return xls.stack()


def get_mapping(df) -> (dict, set):
    mapping = {}
    nodes = set()
    for v in df.index.unique("Variable"):
        # skip aggregations
        if v.count("|") < 2:
            continue

        if v.startswith("Primary"):
            _, bus_carrier, tech = v.split("|")
            mapping[v] = (tech, bus_carrier)
            nodes.add(tech)
            nodes.add(bus_carrier)
        elif v.startswith("Secondary"):
            _, bc_output, bc_input, tech = v.split("|")
            nodes.add(tech)
            nodes.add(bc_input)
            nodes.add(bc_output)
            if bc_output == "Demand":
                mapping[v] = (bc_input, tech)
            elif bc_output == "Losses":
                mapping[v] = (tech, bc_output)
            else:  # Link supply
                mapping[v] = (tech, bc_input)
        elif v.startswith("Final"):
            _, bus_carrier, tech = v.split("|")
            mapping[v] = (bus_carrier, tech)
            nodes.add(tech)
            nodes.add(bus_carrier)
        else:
            raise ValueError(f"Unexpected variable '{v}'")

    return mapping, nodes


def sort_mapping(k):
    if k.startswith("Primary"):
        return 0
    elif k.startswith("Secondary"):
        return 1
    elif k.startswith("Final"):
        return 2
    else:
        raise ValueError(f"Unexpected key '{k}'")


def remove_missing_variables(m: dict) -> dict:
    clean_mapping = {}
    variables = df.index.unique("Variable")
    for k, v in m.items():
        if k in variables:
            clean_mapping[k] = v
        else:
            print(f"Skipping '{k}' because it does not exist in AT {year}.")
    return clean_mapping


def get_xmap(nodes) -> dict:
    # dict.fromkeys(sorted(nodes), "")
    return {
        "AC": 0.5,
        "Agriculture": 1,
        "Air Heat Pump": 0.5,
        "Ambient Heat": "",
        "BEV charger": 0.5,
        "Base Load": "",
        "Biogas CC": 0.5,
        "Biomass": 0.5,
        "Boiler": 0.5,
        "CHP": 0.5,
        "Distribution Grid": 0.5,
        "Electrolysis": 0.5,
        "Export": 0.0,
        "Export Domestic": 1.0,
        "Export Foreign": 1.0,
        "Fischer-Tropsch": 0.5,
        "Gas": 0.5,
        "Gas Compressing": 0.5,
        "Ground Heat Pump": 0.5,
        "H2": 0.5,
        "H2 Compressing": 0.5,
        "HH & Services": 1.0,
        "HVC from naphtha": 0.5,
        "HVC to air": 0.5,
        "Heat": 0.5,
        "Import Domestic": 0.0,
        "Import Foreign": 0.0,
        "Import Global": 0.0,
        "Industry": 1.0,
        "Industry CC": 1.0,
        "Losses": 1.0,
        "Methanol": 0.5,
        "Methanolisation": 0.5,
        "Oil": 0.5,
        "Powerplant": 0.5,
        "Resistive Heater": 0.5,
        "Run-of-River": 0.0,
        "Sabatier": 0.5,
        "Solar Rooftop": 0.0,
        "Solar Utility": 0.0,
        "Solid": 0.5,
        "Transport": 1.0,
        "Waste": 0.5,
        "Water Pits": 0.5,
        "Water Tank": 0.5,
        "Wind Onshore": 0.0,
    }


if __name__ == "__main__":
    df = read_iamc_data_frame(
        filepath="/IdeaProjects/pypsa-at/results/v2025.03/AT10_KN2040/evaluation/exported_iamc_variables.xlsx"
    )
    mapping, nodes = get_mapping(df)
    mapping_sorted = {k: mapping[k] for k in sorted(mapping, key=sort_mapping)}
    xmap = get_xmap(nodes)
    df = rename_aggregate(df, "TWh", level="Unit").div(1e6)
    year = "2050"

    at_regions = [s for s in df.index.unique("Region") if s.startswith("AT")]
    df = filter_by(df, Year=year, Region=at_regions)
    df = rename_aggregate(df, "AT", level="Region")

    # AC
    variable_mapper_ac = {
        # "Primary Energy|AC|Import Domestic": "Primary Energy|AC|Total Import",
        "Primary Energy|AC|Import Foreign": "Primary Energy|AC|Total Import",
        "Primary Energy|AC|Reservoir": "Primary Energy|AC|Hydro Power",
        "Primary Energy|AC|Run-of-River": "Primary Energy|AC|Hydro Power",
        "Primary Energy|AC|Solar HSAT": "Primary Energy|AC|Solar Power",
        "Primary Energy|AC|Solar Rooftop": "Primary Energy|AC|Solar Power",
        "Primary Energy|AC|Solar Utility": "Primary Energy|AC|Solar Power",
        "Primary Energy|AC|Wind Onshore": "Primary Energy|AC|Wind Power",
        "Primary Energy|AC|Wind Offshore": "Primary Energy|AC|Wind Power",
        "Secondary Energy|AC|Biomass|CHP": "Secondary Energy|AC|CHP",
        "Secondary Energy|AC|Biomass|CHP CC": "Secondary Energy|AC|CHP",
        "Secondary Energy|AC|Gas|CHP": "Secondary Energy|AC|CHP",
        "Secondary Energy|AC|Gas|CHP CC": "Secondary Energy|AC|CHP",
        "Secondary Energy|AC|Waste|CHP": "Secondary Energy|AC|CHP",
        "Secondary Energy|AC|Waste|CHP CC": "Secondary Energy|AC|CHP",
        "Secondary Energy|AC|Gas|Powerplant": "Secondary Energy|AC|Power Plant",
        "Secondary Energy|AC|H2|Powerplant": "Secondary Energy|AC|Power Plant",
        "Secondary Energy|AC|Methanol|Powerplant": "Secondary Energy|AC|Power Plant",
        "Secondary Energy|Demand|AC|Air Heat Pump": "Secondary Energy|Demand|AC|Heat Pump",
        "Secondary Energy|Demand|AC|Ground Heat Pump": "Secondary Energy|Demand|AC|Heat Pump",
        # "Final Energy|AC|Export Domestic": "Final Energy|AC|Total Export",
        "Final Energy|AC|Export Foreign": "Final Energy|AC|Total Export",
        "Secondary Energy|Losses|AC|BEV charger": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Battery storage": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Distribution Grid": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Electrolysis": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Haber-Bosch": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Home Battery storage": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Methanolisation": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|Resistive Heater": "Secondary Energy|Losses|AC",
        "Secondary Energy|Losses|AC|V2G": "Secondary Energy|Losses|AC",
    }
    variable_mapper_gas = {
        "Primary Energy|Gas|Biogas": "Primary Energy|Gas|Biogas",
        "Primary Energy|Gas|Biogas CC": "Primary Energy|Gas|Biogas",
        "Primary Energy|Gas|Domestic Production": "Primary Energy|Gas|Biogas",
        "Primary Energy|Gas|Global Import LNG": "",
        "Primary Energy|Gas|Global Import Pipeline": "",
        "Primary Energy|Gas|Green Global Import": "",
        "Primary Energy|Gas|Import Domestic": "",
        "Primary Energy|Gas|Import Foreign": "",
    }
    df = rename_aggregate(df, variable_mapper_ac, level="Variable")

    # _idx = list(df.index[0])
    # _idx[3] = "Primary Energy|AC"
    # df.loc[_idx] = df.query("Variable.str.startswith('Primary Energy|AC|') & 'Domestic' not in Variable")

    mapping = {
        "Primary Energy|AC|Total Import": ("Import", "AC"),
        "Primary Energy|AC|Hydro Power": ("Hydro Power", "AC"),
        "Primary Energy|AC|Solar Power": ("Solar Power", "AC"),
        "Primary Energy|AC|Wind Power": ("Wind Power", "AC"),
        # "Primary Energy|AC": ("AC", "AC"),
        # supply
        "Secondary Energy|AC|CHP": ("CHP", "AC"),
        "Secondary Energy|AC|Power Plant": ("Power Plant", "AC"),
        # demand
        "Secondary Energy|Demand|AC|Heat Pump": ("AC", "Heat Pump"),
        "Secondary Energy|Demand|AC|DAC": ("AC", "DAC"),
        "Secondary Energy|Demand|AC|Electrolysis": ("AC", "Electrolysis"),
        "Secondary Energy|Demand|AC|Gas Compressing": (
            "AC",
            "Gas Compressing",
        ),
        "Secondary Energy|Demand|AC|H2 Compressing": ("AC", "H2 Compressing"),
        "Secondary Energy|Demand|AC|Haber-Bosch": ("AC", "Haber-Bosch"),
        "Secondary Energy|Demand|AC|Methanolisation": (
            "AC",
            "Methanolisation",
        ),
        "Secondary Energy|Demand|AC|Resistive Heater": (
            "AC",
            "Heat Secondary",
        ),
        "Secondary Energy|Losses|AC": ("AC", "Losses"),
        # Load
        # "Final Energy|AC": ("AC", "AC Final"),
        "Final Energy|AC|Agriculture": ("AC", "Agriculture"),
        "Final Energy|AC|Base Load": ("AC", "Base Load"),
        "Final Energy|AC|Total Export": ("AC", "Export"),
        "Final Energy|AC|Industry": ("AC", "Industry"),
        "Final Energy|AC|Transport": ("AC", "Transport"),
        # Secondary Energy|Losses|AC|BEV charger
        # Secondary Energy|Losses|AC|Battery storage
        # Secondary Energy|Losses|AC|Distribution Grid
        # Secondary Energy|Losses|AC|Electrolysis
        # Secondary Energy|Losses|AC|Haber-Bosch
        # Secondary Energy|Losses|AC|Home Battery storage
        # Secondary Energy|Losses|AC|Methanolisation
        # Secondary Energy|Losses|AC|Resistive Heater
        # Secondary Energy|Losses|AC|V2G
    }

    # clean_mapping = remove_missing_variables(mapping_sorted)

    iamc = pyam.IamDataFrame(df)

    iamc_fig = sankey(iamc, mapping=mapping)
    iamc_fig.update_layout(height=800)
    # node = iamc_fig.data[0].node.to_plotly_json()
    # link = iamc_fig.data[0].link.to_plotly_json()
    #
    # node["x"] = [xmap.get(label, 0.2) for label in node["label"]]
    # node["y"] = [xmap.get(label, 0.4) for label in node["label"]]
    #
    # new_sankey = Sankey(
    #     node=node,
    #     link=link,
    #     arrangement="fixed",  # necessary for x/y positions
    # )
    #
    # fig = Figure(data=[new_sankey])
    #
    # fig.update_layout(height=800)

    plotly.io.show(iamc_fig)
