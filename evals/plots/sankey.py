# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for Sankey diagram."""

import re

import pandas as pd
import plotly
import pyam
from plotly.graph_objs import Figure, Sankey
from pyam.index import get_index_levels

from evals.constants import COLOUR, COLOUR_SCHEME, RUN_META_DATA
from evals.constants import DataModel as DM
from evals.plots._base import ESMChart
from evals.utils import (
    filter_by,
    prettify_number,
    rename_aggregate,
)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)

# idx = df_plot.sort_index(level="bus_carrier").droplevel("year").droplevel("location").index.drop_duplicates()
# dict.fromkeys(idx, ("", ""))
mapping = {
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
    ("Link", "DAC", "AC"): ("AC Secondary Input", "DAC"),
    ("Link", "H2 Fuel Cell", "AC"): ("Power plant", "AC Secondary Input"),
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
    ("Link", "BEV charger", "EV battery"): ("AC Secondary Output", "Car Battery"),
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
    ("Link", "H2 Fuel Cell", "H2"): ("H2 Secondary Output", "Power Plant"),
    ("Link", "H2 Fuel Cell losses", "H2"): ("Power Plant", "Power Plant Losses"),
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
    ("Generator", "import NH3", "NH3"): ("", ""),
    ("Link", "Haber-Bosch", "NH3"): ("", ""),
    ("Link", "ammonia cracker", "NH3"): ("", ""),
    ("Link", "ammonia cracker losses", "NH3"): ("", ""),
    ("Load", "NH3", "NH3"): ("", ""),
    ("Store", "ammonia store", "NH3"): ("", ""),
    ("Link", "agriculture machinery oil", "agriculture machinery oil"): ("", ""),
    ("Load", "agriculture machinery oil", "agriculture machinery oil"): ("", ""),
    ("Link", "rural air heat pump", "ambient heat"): ("", ""),
    ("Link", "rural ground heat pump", "ambient heat"): ("", ""),
    ("Link", "urban decentral air heat pump", "ambient heat"): ("", ""),
    ("Link", "electricity distribution grid", "ambient heat"): ("", ""),
    ("Link", "urban central air heat pump", "ambient heat"): ("", ""),
    ("Link", "urban central coal CHP", "ambient heat"): ("", ""),
    ("Link", "urban central gas boiler", "ambient heat"): ("", ""),
    ("Link", "urban central ptes heat pump", "ambient heat"): ("", ""),
    ("Link", "urban central solid biomass CHP", "ambient heat"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "ambient heat"): ("", ""),
    ("Link", "urban central oil CHP", "ambient heat"): ("", ""),
    ("Link", "urban central gas CHP", "ambient heat"): ("", ""),
    ("Link", "urban central lignite CHP", "ambient heat"): ("", ""),
    ("Link", "battery charger", "battery"): ("", ""),
    ("Link", "battery discharger", "battery"): ("", ""),
    ("Link", "battery discharger losses", "battery"): ("", ""),
    ("Store", "battery", "battery"): ("", ""),
    ("Generator", "unsustainable biogas", "biogas"): ("", ""),
    ("Generator", "biogas", "biogas"): ("", ""),
    ("Link", "biogas to gas", "biogas"): ("", ""),
    ("Link", "biogas to gas CC", "biogas"): ("", ""),
    ("Generator", "coal", "coal"): ("", ""),
    ("Link", "coal for industry", "coal"): ("", ""),
    ("Link", "urban central coal CHP", "coal"): ("", ""),
    ("Link", "coal", "coal"): ("", ""),
    ("Link", "coal losses", "coal"): ("", ""),
    ("Link", "urban central coal CHP losses", "coal"): ("", ""),
    ("Store", "coal", "coal"): ("", ""),
    ("Load", "coal for industry", "coal for industry"): ("", ""),
    ("Generator", "production gas", "gas"): ("", ""),
    ("Generator", "lng gas", "gas"): ("", ""),
    ("Generator", "pipeline gas", "gas"): ("", ""),
    ("Link", "BioSNG", "gas"): ("", ""),
    ("Link", "BioSNG CC", "gas"): ("", ""),
    ("Link", "Export Foreign", "gas"): ("", ""),
    ("Link", "Import Foreign", "gas"): ("", ""),
    ("Link", "OCGT", "gas"): ("", ""),
    ("Link", "OCGT losses", "gas"): ("", ""),
    ("Link", "SMR", "gas"): ("", ""),
    ("Link", "SMR CC", "gas"): ("", ""),
    ("Link", "SMR CC losses", "gas"): ("", ""),
    ("Link", "SMR losses", "gas"): ("", ""),
    ("Link", "Sabatier", "gas"): ("", ""),
    ("Link", "gas for industry", "gas"): ("", ""),
    ("Link", "gas for industry CC", "gas"): ("", ""),
    ("Link", "rural gas boiler", "gas"): ("", ""),
    ("Link", "rural gas boiler losses", "gas"): ("", ""),
    ("Link", "urban decentral gas boiler", "gas"): ("", ""),
    ("Link", "urban decentral gas boiler losses", "gas"): ("", ""),
    ("Link", "CCGT", "gas"): ("", ""),
    ("Link", "CCGT losses", "gas"): ("", ""),
    ("Link", "biogas to gas", "gas"): ("", ""),
    ("Link", "biogas to gas CC", "gas"): ("", ""),
    ("Link", "urban central gas CHP", "gas"): ("", ""),
    ("Link", "urban central gas CHP CC", "gas"): ("", ""),
    ("Link", "urban central gas CHP CC losses", "gas"): ("", ""),
    ("Link", "urban central gas CHP losses", "gas"): ("", ""),
    ("Link", "urban central gas boiler", "gas"): ("", ""),
    ("Link", "Export Domestic", "gas"): ("", ""),
    ("Link", "Import Domestic", "gas"): ("", ""),
    ("Link", "import gas", "gas"): ("", ""),
    ("Store", "gas", "gas"): ("", ""),
    ("Load", "gas for industry", "gas for industry"): ("", ""),
    ("Link", "home battery charger", "home battery"): ("", ""),
    ("Link", "home battery discharger", "home battery"): ("", ""),
    ("Link", "home battery discharger losses", "home battery"): ("", ""),
    ("Store", "home battery", "home battery"): ("", ""),
    ("Link", "industry methanol", "industry methanol"): ("", ""),
    ("Load", "industry methanol", "industry methanol"): ("", ""),
    ("Link", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Link", "methanol-to-kerosene", "kerosene for aviation"): ("", ""),
    ("Load", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Link", "land transport oil", "land transport oil"): ("", ""),
    ("Load", "land transport oil", "land transport oil"): ("", ""),
    ("Generator", "lignite", "lignite"): ("", ""),
    ("Link", "lignite", "lignite"): ("", ""),
    ("Link", "lignite losses", "lignite"): ("", ""),
    ("Link", "urban central lignite CHP", "lignite"): ("", ""),
    ("Link", "urban central lignite CHP losses", "lignite"): ("", ""),
    ("Store", "lignite", "lignite"): ("", ""),
    ("Generator", "solar rooftop", "low voltage"): ("", ""),
    ("Link", "BEV charger", "low voltage"): ("", ""),
    ("Link", "BEV charger losses", "low voltage"): ("", ""),
    ("Link", "electricity distribution grid", "low voltage"): ("", ""),
    ("Link", "electricity distribution grid losses", "low voltage"): ("", ""),
    ("Link", "home battery charger", "low voltage"): ("", ""),
    ("Link", "home battery charger losses", "low voltage"): ("", ""),
    ("Link", "home battery discharger", "low voltage"): ("", ""),
    ("Link", "rural air heat pump", "low voltage"): ("", ""),
    ("Link", "rural ground heat pump", "low voltage"): ("", ""),
    ("Link", "rural resistive heater", "low voltage"): ("", ""),
    ("Link", "rural resistive heater losses", "low voltage"): ("", ""),
    ("Link", "urban decentral air heat pump", "low voltage"): ("", ""),
    ("Link", "urban decentral resistive heater", "low voltage"): ("", ""),
    ("Link", "urban decentral resistive heater losses", "low voltage"): ("", ""),
    ("Link", "urban central air heat pump", "low voltage"): ("", ""),
    ("Link", "urban central ptes heat pump", "low voltage"): ("", ""),
    ("Link", "urban central resistive heater", "low voltage"): ("", ""),
    ("Link", "urban central resistive heater losses", "low voltage"): ("", ""),
    ("Link", "urban central ptes heat pump losses", "low voltage"): ("", ""),
    ("Link", "V2G", "low voltage"): ("", ""),
    ("Load", "agriculture electricity", "low voltage"): ("", ""),
    ("Load", "electricity", "low voltage"): ("", ""),
    ("Load", "industry electricity", "low voltage"): ("", ""),
    ("Link", "CCGT methanol", "methanol"): ("", ""),
    ("Link", "CCGT methanol CC", "methanol"): ("", ""),
    ("Link", "CCGT methanol CC losses", "methanol"): ("", ""),
    ("Link", "CCGT methanol losses", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming CC", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming CC losses", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming losses", "methanol"): ("", ""),
    ("Link", "OCGT methanol", "methanol"): ("", ""),
    ("Link", "OCGT methanol losses", "methanol"): ("", ""),
    ("Link", "allam methanol", "methanol"): ("", ""),
    ("Link", "allam methanol losses", "methanol"): ("", ""),
    ("Link", "biomass-to-methanol", "methanol"): ("", ""),
    ("Link", "biomass-to-methanol CC", "methanol"): ("", ""),
    ("Link", "methanol-to-kerosene", "methanol"): ("", ""),
    ("Link", "methanol-to-kerosene losses", "methanol"): ("", ""),
    ("Link", "methanolisation", "methanol"): ("", ""),
    ("Link", "industry methanol", "methanol"): ("", ""),
    ("Link", "import methanol", "methanol"): ("", ""),
    ("Link", "shipping methanol", "methanol"): ("", ""),
    ("Store", "methanol", "methanol"): ("", ""),
    ("Generator", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Link", "Export Foreign", "municipal solid waste"): ("", ""),
    ("Link", "Import Foreign", "municipal solid waste"): ("", ""),
    ("Link", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Link", "Export Domestic", "municipal solid waste"): ("", ""),
    ("Link", "Import Domestic", "municipal solid waste"): ("", ""),
    ("Load", "naphtha for industry", "naphtha for industry"): ("", ""),
    ("Link", "HVC to air", "non-sequestered HVC"): ("", ""),
    ("Link", "municipal solid waste", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP CC", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP CC losses", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP losses", "non-sequestered HVC"): ("", ""),
    ("Store", "non-sequestered HVC", "non-sequestered HVC"): ("", ""),
    ("Link", "Fischer-Tropsch", "oil"): ("", ""),
    ("Link", "agriculture machinery oil", "oil"): ("", ""),
    ("Link", "biomass to liquid", "oil"): ("", ""),
    ("Link", "biomass to liquid CC", "oil"): ("", ""),
    ("Link", "electrobiofuels", "oil"): ("", ""),
    ("Link", "kerosene for aviation", "oil"): ("", ""),
    ("Link", "land transport oil", "oil"): ("", ""),
    ("Link", "naphtha for industry", "oil"): ("", ""),
    ("Link", "oil", "oil"): ("", ""),
    ("Link", "oil losses", "oil"): ("", ""),
    ("Link", "shipping oil", "oil"): ("", ""),
    ("Link", "rural oil boiler", "oil"): ("", ""),
    ("Link", "rural oil boiler losses", "oil"): ("", ""),
    ("Link", "unsustainable bioliquids", "oil"): ("", ""),
    ("Link", "urban decentral oil boiler", "oil"): ("", ""),
    ("Link", "urban decentral oil boiler losses", "oil"): ("", ""),
    ("Link", "urban central oil CHP", "oil"): ("", ""),
    ("Link", "urban central oil CHP losses", "oil"): ("", ""),
    ("Link", "import oil", "oil"): ("", ""),
    ("Link", "oil refining", "oil"): ("", ""),
    ("Store", "oil", "oil"): ("", ""),
    ("Generator", "oil primary", "oil primary"): ("", ""),
    ("Link", "oil refining", "oil primary"): ("", ""),
    ("Link", "oil refining losses", "oil primary"): ("", ""),
    ("Generator", "rural heat vent", "rural heat"): ("", ""),
    ("Generator", "rural solar thermal", "rural heat"): ("", ""),
    ("Link", "rural air heat pump", "rural heat"): ("", ""),
    ("Link", "rural biomass boiler", "rural heat"): ("", ""),
    ("Link", "rural gas boiler", "rural heat"): ("", ""),
    ("Link", "rural ground heat pump", "rural heat"): ("", ""),
    ("Link", "rural resistive heater", "rural heat"): ("", ""),
    ("Link", "rural oil boiler", "rural heat"): ("", ""),
    ("Load", "agriculture heat", "rural heat"): ("", ""),
    ("Load", "rural heat", "rural heat"): ("", ""),
    ("Store", "rural water tanks", "rural water tanks"): ("", ""),
    ("Link", "shipping methanol", "shipping methanol"): ("", ""),
    ("Load", "shipping methanol", "shipping methanol"): ("", ""),
    ("Link", "shipping oil", "shipping oil"): ("", ""),
    ("Load", "shipping oil", "shipping oil"): ("", ""),
    ("Generator", "unsustainable solid biomass", "solid biomass"): ("", ""),
    ("Generator", "solid biomass", "solid biomass"): ("", ""),
    ("Link", "BioSNG", "solid biomass"): ("", ""),
    ("Link", "BioSNG CC", "solid biomass"): ("", ""),
    ("Link", "BioSNG CC losses", "solid biomass"): ("", ""),
    ("Link", "BioSNG losses", "solid biomass"): ("", ""),
    ("Link", "Export Foreign", "solid biomass"): ("", ""),
    ("Link", "Import Foreign", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid CC", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid CC losses", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid losses", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol CC", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol CC losses", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol losses", "solid biomass"): ("", ""),
    ("Link", "electrobiofuels", "solid biomass"): ("", ""),
    ("Link", "electrobiofuels losses", "solid biomass"): ("", ""),
    ("Link", "rural biomass boiler", "solid biomass"): ("", ""),
    ("Link", "rural biomass boiler losses", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry CC", "solid biomass"): ("", ""),
    ("Link", "solid biomass to hydrogen", "solid biomass"): ("", ""),
    ("Link", "solid biomass to hydrogen losses", "solid biomass"): ("", ""),
    ("Link", "urban decentral biomass boiler", "solid biomass"): ("", ""),
    ("Link", "urban decentral biomass boiler losses", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "solid biomass"): ("", ""),
    ("Link", "Export Domestic", "solid biomass"): ("", ""),
    ("Link", "Import Domestic", "solid biomass"): ("", ""),
    ("Link", "solid biomass", "solid biomass"): ("", ""),
    ("Link", "solid biomass losses", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP losses", "solid biomass"): ("", ""),
    ("Load", "solid biomass for industry", "solid biomass for industry"): ("", ""),
    ("Generator", "unsustainable bioliquids", "unsustainable bioliquids"): ("", ""),
    ("Link", "unsustainable bioliquids", "unsustainable bioliquids"): ("", ""),
    ("Generator", "uranium", "uranium"): ("", ""),
    ("Link", "nuclear", "uranium"): ("", ""),
    ("Link", "nuclear losses", "uranium"): ("", ""),
    ("Store", "uranium", "uranium"): ("", ""),
    ("Generator", "urban central heat vent", "urban central heat"): ("", ""),
    ("Generator", "urban central solar thermal", "urban central heat"): ("", ""),
    ("Link", "DAC", "urban central heat"): ("", ""),
    ("Link", "Fischer-Tropsch", "urban central heat"): ("", ""),
    ("Link", "H2 Fuel Cell", "urban central heat"): ("", ""),
    ("Link", "Haber-Bosch", "urban central heat"): ("", ""),
    ("Link", "Sabatier", "urban central heat"): ("", ""),
    ("Link", "methanolisation", "urban central heat"): ("", ""),
    ("Link", "urban central air heat pump", "urban central heat"): ("", ""),
    ("Link", "urban central coal CHP", "urban central heat"): ("", ""),
    ("Link", "urban central gas CHP", "urban central heat"): ("", ""),
    ("Link", "urban central gas CHP CC", "urban central heat"): ("", ""),
    ("Link", "urban central gas boiler", "urban central heat"): ("", ""),
    ("Link", "urban central ptes heat pump", "urban central heat"): ("", ""),
    ("Link", "urban central resistive heater", "urban central heat"): ("", ""),
    ("Link", "urban central solid biomass CHP", "urban central heat"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "urban central heat"): ("", ""),
    ("Link", "waste CHP", "urban central heat"): ("", ""),
    ("Link", "waste CHP CC", "urban central heat"): ("", ""),
    ("Link", "urban central oil CHP", "urban central heat"): ("", ""),
    ("Link", "urban central lignite CHP", "urban central heat"): ("", ""),
    ("Link", "H2 Electrolysis", "urban central heat"): ("", ""),
    ("Link", "urban central H2 CHP", "urban central heat"): ("", ""),
    ("Load", "low-temperature heat for industry", "urban central heat"): ("", ""),
    ("Load", "urban central heat", "urban central heat"): ("", ""),
    ("Store", "urban central water pits", "urban central water pits"): ("", ""),
    ("Store", "urban central water tanks", "urban central water tanks"): ("", ""),
    ("Generator", "urban decentral heat vent", "urban decentral heat"): ("", ""),
    ("Generator", "urban decentral solar thermal", "urban decentral heat"): ("", ""),
    ("Link", "DAC", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral air heat pump", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral biomass boiler", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral gas boiler", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral resistive heater", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral oil boiler", "urban decentral heat"): ("", ""),
    ("Load", "low-temperature heat for industry", "urban decentral heat"): ("", ""),
    ("Load", "urban decentral heat", "urban decentral heat"): ("", ""),
    ("Store", "urban decentral water tanks", "urban decentral water tanks"): ("", ""),
}


class SankeyChart(ESMChart):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.location = self._df.index.unique(DM.LOCATION).item()
        self.year = self._df.index.unique(DM.YEAR).item()
        self._df = self._df.droplevel(DM.YEAR).droplevel(DM.LOCATION)
        self._df.columns = ["value"]

    @staticmethod
    def _add_source_target_columns(idx: tuple) -> pd.Series:
        return pd.Series(mapping.get(idx, ""))

    def plot(self):
        # Concatenate the data with source and target columns

        _df = self._df.copy()  # preserve original

        # [*_df.groupby(level=_df.index.names)][0]

        # df_agg = _df.groupby(level=0).apply(self.custom_aggregate).reset_index()

        _df["index"] = _df.index  # convert to tuples
        _df[["source", "target"]] = _df["index"].apply(self._add_source_target_columns)
        _df = _df.drop(columns=["index"])
        df_agg = _df

        # add jumpers:
        for bus_carrier in ("AC", "H2"):
            primary_to_secondary = filter_by(df_agg, target=f"{bus_carrier} Primary")
            row_idx = ("Jumper", "primary to secondary", bus_carrier)
            df_agg.loc[row_idx, ["value", "source", "target"]] = [
                primary_to_secondary.value.sum(),
                f"{bus_carrier} Primary",
                f"{bus_carrier} Secondary Input",
            ]

            secondary_forwarding = filter_by(
                df_agg, target=f"{bus_carrier} Secondary Input"
            )
            row_idx = ("Jumper", "secondary forwaring", bus_carrier)
            df_agg.loc[row_idx, ["value", "source", "target"]] = [
                secondary_forwarding.value.sum(),
                f"{bus_carrier} Secondary Input",
                f"{bus_carrier} Secondary Output",
            ]

        label_mapping = {
            label: i
            for i, label in enumerate(set(pd.concat([_df["source"], _df["target"]])))
        }
        df_agg["source_id"] = df_agg["source"].map(label_mapping)
        df_agg["target_id"] = df_agg["target"].map(label_mapping)

        # _df = df_mapping.merge(df._data, how="left", left_index=True, right_on="variable")
        # _df = _df.replace(label_mapping)

        # df_agg["customdata"] = df_agg.groupby(["source", "target"]).apply(self._map_customdata)  # , include_groups=False
        to_concat = []
        for _, data in df_agg.groupby(["source", "target"]):
            source, target = _

            if source == "" and target == "":
                continue

            # assert data.index.unique("bus_carrier").item()
            # assert data.index.unique("component").item()
            data = data.reset_index()

            carrier_values = [
                f"{c}: {prettify_number(v)} {self.unit}"
                for c, v in zip(data["carrier"], data["value"])
            ]
            data["link_customdata"] = "<br>".join(carrier_values)

            to_concat.append(data)

        df_agg = pd.concat(to_concat)

        df_agg = (
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

        # def get_carrier_color(s) -> str:
        #     carrier = re.findall(r"\|AC", s)[0].strip("|")
        #     color_map = {"AC": COLOUR.red}
        #     return color_map[carrier]

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

        df_agg["color"] = df_agg["bus_carrier"].map(BUS_CARRIER_COLORS)

        label = pd.Series(list(label_mapping))

        self.fig = Figure(
            data=[
                Sankey(
                    valuesuffix=self.unit,
                    node=dict(
                        # pad=15,
                        # thickness=10,
                        line=dict(color="black", width=0.5),
                        label=label,
                        hovertemplate="%{label}: %{value}<extra></extra>",
                        color=df_agg.color,
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

        plotly.io.show(self.fig)  # todo: remove debugging

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


class SankeyNode:
    def __init__(self, bus_carrier, label, variables):
        self.bus_carrier: str = bus_carrier
        self.label = label
        self.variables = variables
        self.color = COLOUR_SCHEME[bus_carrier]


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
