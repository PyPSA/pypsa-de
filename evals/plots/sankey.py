# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for Sankey diagram."""

import re

import numpy as np
import pandas as pd
import plotly
import pyam
from plotly.graph_objs import Figure, Sankey
from pyam.index import get_index_levels

from evals.constants import COLOUR, COLOUR_SCHEME
from evals.constants import DataModel as DM
from evals.plots._base import ESMChart
from evals.utils import filter_by, rename_aggregate

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)

mapping = {
    ("Generator", "biogas", "biogas"): ("src", "dst"),
    ("Generator", "coal", "coal"): ("src", "dst"),
    ("Generator", "import H2", "H2"): ("src", "dst"),
    ("Generator", "import NH3", "NH3"): ("src", "dst"),
    ("Generator", "lignite", "lignite"): ("src", "dst"),
    ("Generator", "lng gas", "gas"): ("src", "dst"),
    ("Generator", "municipal solid waste", "municipal solid waste"): ("src", "dst"),
    ("Generator", "offwind-ac", "AC"): ("src", "dst"),
    ("Generator", "offwind-dc", "AC"): ("src", "dst"),
    ("Generator", "oil primary", "oil primary"): ("src", "dst"),
    ("Generator", "onwind", "AC"): ("src", "dst"),
    ("Generator", "pipeline gas", "gas"): ("src", "dst"),
    ("Generator", "production gas", "gas"): ("src", "dst"),
    ("Generator", "ror", "AC"): ("src", "dst"),
    ("Generator", "rural heat vent", "rural heat"): ("src", "dst"),
    ("Generator", "rural solar thermal", "rural heat"): ("src", "dst"),
    ("Generator", "solar", "AC"): ("src", "dst"),
    ("Generator", "solar rooftop", "low voltage"): ("src", "dst"),
    ("Generator", "solar-hsat", "AC"): ("src", "dst"),
    ("Generator", "solid biomass", "solid biomass"): ("src", "dst"),
    ("Generator", "unsustainable biogas", "biogas"): ("src", "dst"),
    ("Generator", "unsustainable bioliquids", "unsustainable bioliquids"): (
        "src",
        "dst",
    ),
    ("Generator", "unsustainable solid biomass", "solid biomass"): ("src", "dst"),
    ("Generator", "uranium", "uranium"): ("src", "dst"),
    ("Generator", "urban central heat vent", "urban central heat"): ("src", "dst"),
    ("Generator", "urban central solar thermal", "urban central heat"): ("src", "dst"),
    ("Generator", "urban decentral heat vent", "urban decentral heat"): ("src", "dst"),
    ("Generator", "urban decentral solar thermal", "urban decentral heat"): (
        "src",
        "dst",
    ),
    ("Line", "AC", "AC"): ("src", "dst"),
    ("Link", "BEV charger", "EV battery"): ("src", "dst"),
    ("Link", "BEV charger", "low voltage"): ("src", "dst"),
    ("Link", "BioSNG", "co2"): ("src", "dst"),
    ("Link", "BioSNG", "gas"): ("src", "dst"),
    ("Link", "BioSNG", "solid biomass"): ("src", "dst"),
    ("Link", "BioSNG CC", "co2"): ("src", "dst"),
    ("Link", "BioSNG CC", "co2 stored"): ("src", "dst"),
    ("Link", "BioSNG CC", "gas"): ("src", "dst"),
    ("Link", "BioSNG CC", "solid biomass"): ("src", "dst"),
    ("Link", "CCGT", "AC"): ("src", "dst"),
    ("Link", "CCGT", "co2"): ("src", "dst"),
    ("Link", "CCGT", "gas"): ("src", "dst"),
    ("Link", "CCGT methanol", "AC"): ("src", "dst"),
    ("Link", "CCGT methanol", "co2"): ("src", "dst"),
    ("Link", "CCGT methanol", "methanol"): ("src", "dst"),
    ("Link", "CCGT methanol CC", "AC"): ("src", "dst"),
    ("Link", "CCGT methanol CC", "co2"): ("src", "dst"),
    ("Link", "CCGT methanol CC", "co2 stored"): ("src", "dst"),
    ("Link", "CCGT methanol CC", "methanol"): ("src", "dst"),
    ("Link", "CO2 pipeline", "co2 stored"): ("src", "dst"),
    ("Link", "DAC", "AC"): ("src", "dst"),
    ("Link", "DAC", "co2"): ("src", "dst"),
    ("Link", "DAC", "co2 stored"): ("src", "dst"),
    ("Link", "DAC", "urban central heat"): ("src", "dst"),
    ("Link", "DAC", "urban decentral heat"): ("src", "dst"),
    ("Link", "DC", "AC"): ("src", "dst"),
    ("Link", "Fischer-Tropsch", "H2"): ("src", "dst"),
    ("Link", "Fischer-Tropsch", "co2 stored"): ("src", "dst"),
    ("Link", "Fischer-Tropsch", "oil"): ("src", "dst"),
    ("Link", "Fischer-Tropsch", "urban central heat"): ("src", "dst"),
    ("Link", "H2 Fuel Cell", "AC"): ("src", "dst"),
    ("Link", "H2 Fuel Cell", "H2"): ("src", "dst"),
    ("Link", "H2 Fuel Cell", "urban central heat"): ("src", "dst"),
    ("Link", "HVC to air", "co2"): ("src", "dst"),
    ("Link", "HVC to air", "non-sequestered HVC"): ("src", "dst"),
    ("Link", "Haber-Bosch", "AC"): ("src", "dst"),
    ("Link", "Haber-Bosch", "H2"): ("src", "dst"),
    ("Link", "Haber-Bosch", "NH3"): ("src", "dst"),
    ("Link", "Haber-Bosch", "urban central heat"): ("src", "dst"),
    ("Link", "Methanol steam reforming", "H2"): ("src", "dst"),
    ("Link", "Methanol steam reforming", "co2"): ("src", "dst"),
    ("Link", "Methanol steam reforming", "methanol"): ("src", "dst"),
    ("Link", "Methanol steam reforming CC", "H2"): ("src", "dst"),
    ("Link", "Methanol steam reforming CC", "co2"): ("src", "dst"),
    ("Link", "Methanol steam reforming CC", "co2 stored"): ("src", "dst"),
    ("Link", "Methanol steam reforming CC", "methanol"): ("src", "dst"),
    ("Link", "OCGT", "AC"): ("src", "dst"),
    ("Link", "OCGT", "co2"): ("src", "dst"),
    ("Link", "OCGT", "gas"): ("src", "dst"),
    ("Link", "OCGT methanol", "AC"): ("src", "dst"),
    ("Link", "OCGT methanol", "co2"): ("src", "dst"),
    ("Link", "OCGT methanol", "methanol"): ("src", "dst"),
    ("Link", "SMR", "H2"): ("src", "dst"),
    ("Link", "SMR", "co2"): ("src", "dst"),
    ("Link", "SMR", "gas"): ("src", "dst"),
    ("Link", "SMR CC", "H2"): ("src", "dst"),
    ("Link", "SMR CC", "co2"): ("src", "dst"),
    ("Link", "SMR CC", "co2 stored"): ("src", "dst"),
    ("Link", "SMR CC", "gas"): ("src", "dst"),
    ("Link", "Sabatier", "H2"): ("src", "dst"),
    ("Link", "Sabatier", "co2 stored"): ("src", "dst"),
    ("Link", "Sabatier", "gas"): ("src", "dst"),
    ("Link", "Sabatier", "urban central heat"): ("src", "dst"),
    ("Link", "agriculture machinery oil", "agriculture machinery oil"): ("src", "dst"),
    ("Link", "agriculture machinery oil", "co2"): ("src", "dst"),
    ("Link", "agriculture machinery oil", "oil"): ("src", "dst"),
    ("Link", "allam methanol", "AC"): ("src", "dst"),
    ("Link", "allam methanol", "co2"): ("src", "dst"),
    ("Link", "allam methanol", "co2 stored"): ("src", "dst"),
    ("Link", "allam methanol", "methanol"): ("src", "dst"),
    ("Link", "ammonia cracker", "H2"): ("src", "dst"),
    ("Link", "ammonia cracker", "NH3"): ("src", "dst"),
    ("Link", "battery charger", "AC"): ("src", "dst"),
    ("Link", "battery charger", "battery"): ("src", "dst"),
    ("Link", "battery discharger", "AC"): ("src", "dst"),
    ("Link", "battery discharger", "battery"): ("src", "dst"),
    ("Link", "biogas to gas", "biogas"): ("src", "dst"),
    ("Link", "biogas to gas", "co2"): ("src", "dst"),
    ("Link", "biogas to gas", "gas"): ("src", "dst"),
    ("Link", "biogas to gas CC", "biogas"): ("src", "dst"),
    ("Link", "biogas to gas CC", "co2"): ("src", "dst"),
    ("Link", "biogas to gas CC", "co2 stored"): ("src", "dst"),
    ("Link", "biogas to gas CC", "gas"): ("src", "dst"),
    ("Link", "biomass to liquid", "co2"): ("src", "dst"),
    ("Link", "biomass to liquid", "oil"): ("src", "dst"),
    ("Link", "biomass to liquid", "solid biomass"): ("src", "dst"),
    ("Link", "biomass to liquid CC", "co2"): ("src", "dst"),
    ("Link", "biomass to liquid CC", "co2 stored"): ("src", "dst"),
    ("Link", "biomass to liquid CC", "oil"): ("src", "dst"),
    ("Link", "biomass to liquid CC", "solid biomass"): ("src", "dst"),
    ("Link", "biomass-to-methanol", "co2"): ("src", "dst"),
    ("Link", "biomass-to-methanol", "methanol"): ("src", "dst"),
    ("Link", "biomass-to-methanol", "solid biomass"): ("src", "dst"),
    ("Link", "biomass-to-methanol CC", "co2"): ("src", "dst"),
    ("Link", "biomass-to-methanol CC", "co2 stored"): ("src", "dst"),
    ("Link", "biomass-to-methanol CC", "methanol"): ("src", "dst"),
    ("Link", "biomass-to-methanol CC", "solid biomass"): ("src", "dst"),
    ("Link", "coal", "AC"): ("src", "dst"),
    ("Link", "coal", "co2"): ("src", "dst"),
    ("Link", "coal", "coal"): ("src", "dst"),
    ("Link", "coal for industry", "co2"): ("src", "dst"),
    ("Link", "coal for industry", "coal"): ("src", "dst"),
    ("Link", "coal for industry", "coal for industry"): ("src", "dst"),
    ("Link", "electricity distribution grid", "AC"): ("src", "dst"),
    ("Link", "electricity distribution grid", "low voltage"): ("src", "dst"),
    ("Link", "electrobiofuels", "H2"): ("src", "dst"),
    ("Link", "electrobiofuels", "co2"): ("src", "dst"),
    ("Link", "electrobiofuels", "oil"): ("src", "dst"),
    ("Link", "electrobiofuels", "solid biomass"): ("src", "dst"),
    ("Link", "gas for industry", "co2"): ("src", "dst"),
    ("Link", "gas for industry", "gas"): ("src", "dst"),
    ("Link", "gas for industry", "gas for industry"): ("src", "dst"),
    ("Link", "gas for industry CC", "co2"): ("src", "dst"),
    ("Link", "gas for industry CC", "co2 stored"): ("src", "dst"),
    ("Link", "gas for industry CC", "gas"): ("src", "dst"),
    ("Link", "gas for industry CC", "gas for industry"): ("src", "dst"),
    ("Link", "gas pipeline", "AC"): ("src", "dst"),
    ("Link", "gas pipeline", "gas"): ("src", "dst"),
    ("Link", "gas pipeline new", "gas"): ("src", "dst"),
    ("Link", "home battery charger", "home battery"): ("src", "dst"),
    ("Link", "home battery charger", "low voltage"): ("src", "dst"),
    ("Link", "home battery discharger", "home battery"): ("src", "dst"),
    ("Link", "home battery discharger", "low voltage"): ("src", "dst"),
    ("Link", "import gas", "co2"): ("src", "dst"),
    ("Link", "import gas", "gas"): ("src", "dst"),
    ("Link", "import methanol", "co2"): ("src", "dst"),
    ("Link", "import methanol", "methanol"): ("src", "dst"),
    ("Link", "import oil", "co2"): ("src", "dst"),
    ("Link", "import oil", "oil"): ("src", "dst"),
    ("Link", "industry methanol", "co2"): ("src", "dst"),
    ("Link", "industry methanol", "industry methanol"): ("src", "dst"),
    ("Link", "industry methanol", "methanol"): ("src", "dst"),
    ("Link", "kerosene for aviation", "co2"): ("src", "dst"),
    ("Link", "kerosene for aviation", "kerosene for aviation"): ("src", "dst"),
    ("Link", "kerosene for aviation", "oil"): ("src", "dst"),
    ("Link", "land transport oil", "co2"): ("src", "dst"),
    ("Link", "land transport oil", "land transport oil"): ("src", "dst"),
    ("Link", "land transport oil", "oil"): ("src", "dst"),
    ("Link", "lignite", "AC"): ("src", "dst"),
    ("Link", "lignite", "co2"): ("src", "dst"),
    ("Link", "lignite", "lignite"): ("src", "dst"),
    ("Link", "methanol-to-kerosene", "H2"): ("src", "dst"),
    ("Link", "methanol-to-kerosene", "co2"): ("src", "dst"),
    ("Link", "methanol-to-kerosene", "kerosene for aviation"): ("src", "dst"),
    ("Link", "methanol-to-kerosene", "methanol"): ("src", "dst"),
    ("Link", "methanolisation", "AC"): ("src", "dst"),
    ("Link", "methanolisation", "H2"): ("src", "dst"),
    ("Link", "methanolisation", "co2 stored"): ("src", "dst"),
    ("Link", "methanolisation", "methanol"): ("src", "dst"),
    ("Link", "methanolisation", "urban central heat"): ("src", "dst"),
    ("Link", "municipal solid waste", "co2"): ("src", "dst"),
    ("Link", "municipal solid waste", "municipal solid waste"): ("src", "dst"),
    ("Link", "municipal solid waste", "non-sequestered HVC"): ("src", "dst"),
    ("Link", "municipal solid waste transport", "municipal solid waste"): (
        "src",
        "dst",
    ),
    ("Link", "naphtha for industry", "naphtha for industry"): ("src", "dst"),
    ("Link", "naphtha for industry", "oil"): ("src", "dst"),
    ("Link", "naphtha for industry", "process emissions"): ("src", "dst"),
    ("Link", "nuclear", "AC"): ("src", "dst"),
    ("Link", "nuclear", "uranium"): ("src", "dst"),
    ("Link", "oil", "AC"): ("src", "dst"),
    ("Link", "oil", "co2"): ("src", "dst"),
    ("Link", "oil", "oil"): ("src", "dst"),
    ("Link", "oil refining", "co2"): ("src", "dst"),
    ("Link", "oil refining", "oil"): ("src", "dst"),
    ("Link", "oil refining", "oil primary"): ("src", "dst"),
    ("Link", "process emissions", "co2"): ("src", "dst"),
    ("Link", "process emissions", "process emissions"): ("src", "dst"),
    ("Link", "process emissions CC", "co2"): ("src", "dst"),
    ("Link", "process emissions CC", "co2 stored"): ("src", "dst"),
    ("Link", "process emissions CC", "process emissions"): ("src", "dst"),
    ("Link", "rural air heat pump", "low voltage"): ("src", "dst"),
    ("Link", "rural air heat pump", "rural heat"): ("src", "dst"),
    ("Link", "rural biomass boiler", "rural heat"): ("src", "dst"),
    ("Link", "rural biomass boiler", "solid biomass"): ("src", "dst"),
    ("Link", "rural gas boiler", "co2"): ("src", "dst"),
    ("Link", "rural gas boiler", "gas"): ("src", "dst"),
    ("Link", "rural gas boiler", "rural heat"): ("src", "dst"),
    ("Link", "rural ground heat pump", "low voltage"): ("src", "dst"),
    ("Link", "rural ground heat pump", "rural heat"): ("src", "dst"),
    ("Link", "rural oil boiler", "co2"): ("src", "dst"),
    ("Link", "rural oil boiler", "oil"): ("src", "dst"),
    ("Link", "rural oil boiler", "rural heat"): ("src", "dst"),
    ("Link", "rural resistive heater", "low voltage"): ("src", "dst"),
    ("Link", "rural resistive heater", "rural heat"): ("src", "dst"),
    ("Link", "rural water tanks charger", "rural heat"): ("src", "dst"),
    ("Link", "rural water tanks charger", "rural water tanks"): ("src", "dst"),
    ("Link", "rural water tanks discharger", "rural heat"): ("src", "dst"),
    ("Link", "rural water tanks discharger", "rural water tanks"): ("src", "dst"),
    ("Link", "shipping oil", "co2"): ("src", "dst"),
    ("Link", "shipping oil", "oil"): ("src", "dst"),
    ("Link", "shipping oil", "shipping oil"): ("src", "dst"),
    ("Link", "solid biomass", "AC"): ("src", "dst"),
    ("Link", "solid biomass", "solid biomass"): ("src", "dst"),
    ("Link", "solid biomass for industry", "solid biomass"): ("src", "dst"),
    ("Link", "solid biomass for industry", "solid biomass for industry"): (
        "src",
        "dst",
    ),
    ("Link", "solid biomass for industry CC", "co2"): ("src", "dst"),
    ("Link", "solid biomass for industry CC", "co2 stored"): ("src", "dst"),
    ("Link", "solid biomass for industry CC", "solid biomass"): ("src", "dst"),
    ("Link", "solid biomass for industry CC", "solid biomass for industry"): (
        "src",
        "dst",
    ),
    ("Link", "solid biomass to hydrogen", "H2"): ("src", "dst"),
    ("Link", "solid biomass to hydrogen", "co2"): ("src", "dst"),
    ("Link", "solid biomass to hydrogen", "co2 stored"): ("src", "dst"),
    ("Link", "solid biomass to hydrogen", "solid biomass"): ("src", "dst"),
    ("Link", "solid biomass transport", "solid biomass"): ("src", "dst"),
    ("Link", "unsustainable bioliquids", "co2"): ("src", "dst"),
    ("Link", "unsustainable bioliquids", "oil"): ("src", "dst"),
    ("Link", "unsustainable bioliquids", "unsustainable bioliquids"): ("src", "dst"),
    ("Link", "urban central air heat pump", "low voltage"): ("src", "dst"),
    ("Link", "urban central air heat pump", "urban central heat"): ("src", "dst"),
    ("Link", "urban central coal CHP", "AC"): ("src", "dst"),
    ("Link", "urban central coal CHP", "co2"): ("src", "dst"),
    ("Link", "urban central coal CHP", "coal"): ("src", "dst"),
    ("Link", "urban central coal CHP", "urban central heat"): ("src", "dst"),
    ("Link", "urban central gas CHP", "AC"): ("src", "dst"),
    ("Link", "urban central gas CHP", "co2"): ("src", "dst"),
    ("Link", "urban central gas CHP", "gas"): ("src", "dst"),
    ("Link", "urban central gas CHP", "urban central heat"): ("src", "dst"),
    ("Link", "urban central gas CHP CC", "AC"): ("src", "dst"),
    ("Link", "urban central gas CHP CC", "co2"): ("src", "dst"),
    ("Link", "urban central gas CHP CC", "co2 stored"): ("src", "dst"),
    ("Link", "urban central gas CHP CC", "gas"): ("src", "dst"),
    ("Link", "urban central gas CHP CC", "urban central heat"): ("src", "dst"),
    ("Link", "urban central gas boiler", "co2"): ("src", "dst"),
    ("Link", "urban central gas boiler", "gas"): ("src", "dst"),
    ("Link", "urban central gas boiler", "urban central heat"): ("src", "dst"),
    ("Link", "urban central lignite CHP", "AC"): ("src", "dst"),
    ("Link", "urban central lignite CHP", "co2"): ("src", "dst"),
    ("Link", "urban central lignite CHP", "lignite"): ("src", "dst"),
    ("Link", "urban central lignite CHP", "urban central heat"): ("src", "dst"),
    ("Link", "urban central oil CHP", "AC"): ("src", "dst"),
    ("Link", "urban central oil CHP", "co2"): ("src", "dst"),
    ("Link", "urban central oil CHP", "oil"): ("src", "dst"),
    ("Link", "urban central oil CHP", "urban central heat"): ("src", "dst"),
    ("Link", "urban central ptes heat pump", "low voltage"): ("src", "dst"),
    ("Link", "urban central ptes heat pump", "urban central heat"): ("src", "dst"),
    ("Link", "urban central resistive heater", "low voltage"): ("src", "dst"),
    ("Link", "urban central resistive heater", "urban central heat"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP", "AC"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP", "solid biomass"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP", "urban central heat"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP CC", "AC"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP CC", "co2"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP CC", "co2 stored"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP CC", "solid biomass"): ("src", "dst"),
    ("Link", "urban central solid biomass CHP CC", "urban central heat"): (
        "src",
        "dst",
    ),
    ("Link", "urban central water pits charger", "urban central heat"): ("src", "dst"),
    ("Link", "urban central water pits charger", "urban central water pits"): (
        "src",
        "dst",
    ),
    ("Link", "urban central water pits discharger", "urban central heat"): (
        "src",
        "dst",
    ),
    ("Link", "urban central water pits discharger", "urban central water pits"): (
        "src",
        "dst",
    ),
    ("Link", "urban central water tanks charger", "urban central heat"): ("src", "dst"),
    ("Link", "urban central water tanks charger", "urban central water tanks"): (
        "src",
        "dst",
    ),
    ("Link", "urban central water tanks discharger", "urban central heat"): (
        "src",
        "dst",
    ),
    ("Link", "urban central water tanks discharger", "urban central water tanks"): (
        "src",
        "dst",
    ),
    ("Link", "urban decentral air heat pump", "low voltage"): ("src", "dst"),
    ("Link", "urban decentral air heat pump", "urban decentral heat"): ("src", "dst"),
    ("Link", "urban decentral biomass boiler", "solid biomass"): ("src", "dst"),
    ("Link", "urban decentral biomass boiler", "urban decentral heat"): ("src", "dst"),
    ("Link", "urban decentral gas boiler", "co2"): ("src", "dst"),
    ("Link", "urban decentral gas boiler", "gas"): ("src", "dst"),
    ("Link", "urban decentral gas boiler", "urban decentral heat"): ("src", "dst"),
    ("Link", "urban decentral oil boiler", "co2"): ("src", "dst"),
    ("Link", "urban decentral oil boiler", "oil"): ("src", "dst"),
    ("Link", "urban decentral oil boiler", "urban decentral heat"): ("src", "dst"),
    ("Link", "urban decentral resistive heater", "low voltage"): ("src", "dst"),
    ("Link", "urban decentral resistive heater", "urban decentral heat"): (
        "src",
        "dst",
    ),
    ("Link", "urban decentral water tanks charger", "urban decentral heat"): (
        "src",
        "dst",
    ),
    ("Link", "urban decentral water tanks charger", "urban decentral water tanks"): (
        "src",
        "dst",
    ),
    ("Link", "urban decentral water tanks discharger", "urban decentral heat"): (
        "src",
        "dst",
    ),
    ("Link", "urban decentral water tanks discharger", "urban decentral water tanks"): (
        "src",
        "dst",
    ),
    ("Link", "waste CHP", "AC"): ("src", "dst"),
    ("Link", "waste CHP", "co2"): ("src", "dst"),
    ("Link", "waste CHP", "non-sequestered HVC"): ("src", "dst"),
    ("Link", "waste CHP", "urban central heat"): ("src", "dst"),
    ("Link", "waste CHP CC", "AC"): ("src", "dst"),
    ("Link", "waste CHP CC", "co2"): ("src", "dst"),
    ("Link", "waste CHP CC", "co2 stored"): ("src", "dst"),
    ("Link", "waste CHP CC", "non-sequestered HVC"): ("src", "dst"),
    ("Link", "waste CHP CC", "urban central heat"): ("src", "dst"),
    ("Link", "H2 Electrolysis", "AC"): ("src", "dst"),
    ("Link", "H2 Electrolysis", "H2"): ("src", "dst"),
    ("Link", "H2 Electrolysis", "urban central heat"): ("src", "dst"),
    ("Link", "H2 OCGT", "AC"): ("src", "dst"),
    ("Link", "H2 OCGT", "H2"): ("src", "dst"),
    ("Link", "H2 pipeline", "AC"): ("src", "dst"),
    ("Link", "H2 pipeline", "H2"): ("src", "dst"),
    ("Link", "H2 pipeline (Kernnetz)", "AC"): ("src", "dst"),
    ("Link", "H2 pipeline (Kernnetz)", "H2"): ("src", "dst"),
    ("Link", "V2G", "EV battery"): ("src", "dst"),
    ("Link", "V2G", "low voltage"): ("src", "dst"),
    ("Link", "co2 sequestered", "co2 sequestered"): ("src", "dst"),
    ("Link", "co2 sequestered", "co2 stored"): ("src", "dst"),
    ("Link", "shipping methanol", "co2"): ("src", "dst"),
    ("Link", "shipping methanol", "methanol"): ("src", "dst"),
    ("Link", "shipping methanol", "shipping methanol"): ("src", "dst"),
    ("Link", "urban central H2 CHP", "AC"): ("src", "dst"),
    ("Link", "urban central H2 CHP", "H2"): ("src", "dst"),
    ("Link", "urban central H2 CHP", "urban central heat"): ("src", "dst"),
    ("Link", "H2 pipeline retrofitted", "H2"): ("src", "dst"),
    ("Load", "H2 for industry", "H2"): ("src", "dst"),
    ("Load", "NH3", "NH3"): ("src", "dst"),
    ("Load", "agriculture electricity", "low voltage"): ("src", "dst"),
    ("Load", "agriculture heat", "rural heat"): ("src", "dst"),
    ("Load", "agriculture machinery oil", "agriculture machinery oil"): ("src", "dst"),
    ("Load", "coal for industry", "coal for industry"): ("src", "dst"),
    ("Load", "electricity", "low voltage"): ("src", "dst"),
    ("Load", "gas for industry", "gas for industry"): ("src", "dst"),
    ("Load", "industry electricity", "low voltage"): ("src", "dst"),
    ("Load", "industry methanol", "industry methanol"): ("src", "dst"),
    ("Load", "kerosene for aviation", "kerosene for aviation"): ("src", "dst"),
    ("Load", "land transport EV", "EV battery"): ("src", "dst"),
    ("Load", "land transport fuel cell", "H2"): ("src", "dst"),
    ("Load", "land transport oil", "land transport oil"): ("src", "dst"),
    ("Load", "low-temperature heat for industry", "urban central heat"): ("src", "dst"),
    ("Load", "low-temperature heat for industry", "urban decentral heat"): (
        "src",
        "dst",
    ),
    ("Load", "naphtha for industry", "naphtha for industry"): ("src", "dst"),
    ("Load", "process emissions", "process emissions"): ("src", "dst"),
    ("Load", "rural heat", "rural heat"): ("src", "dst"),
    ("Load", "shipping oil", "shipping oil"): ("src", "dst"),
    ("Load", "solid biomass for industry", "solid biomass for industry"): (
        "src",
        "dst",
    ),
    ("Load", "urban central heat", "urban central heat"): ("src", "dst"),
    ("Load", "urban decentral heat", "urban decentral heat"): ("src", "dst"),
    ("Load", "shipping methanol", "shipping methanol"): ("src", "dst"),
    ("StorageUnit", "PHS", "AC"): ("src", "dst"),
    ("StorageUnit", "hydro", "AC"): ("src", "dst"),
    ("Store", "H2 Store", "H2"): ("src", "dst"),
    ("Store", "ammonia store", "NH3"): ("src", "dst"),
    ("Store", "battery", "battery"): ("src", "dst"),
    ("Store", "co2", "co2"): ("src", "dst"),
    ("Store", "co2 stored", "co2 stored"): ("src", "dst"),
    ("Store", "coal", "coal"): ("src", "dst"),
    ("Store", "gas", "gas"): ("src", "dst"),
    ("Store", "home battery", "home battery"): ("src", "dst"),
    ("Store", "lignite", "lignite"): ("src", "dst"),
    ("Store", "methanol", "methanol"): ("src", "dst"),
    ("Store", "non-sequestered HVC", "non-sequestered HVC"): ("src", "dst"),
    ("Store", "oil", "oil"): ("src", "dst"),
    ("Store", "rural water tanks", "rural water tanks"): ("src", "dst"),
    ("Store", "uranium", "uranium"): ("src", "dst"),
    ("Store", "urban central water pits", "urban central water pits"): ("src", "dst"),
    ("Store", "urban central water tanks", "urban central water tanks"): ("src", "dst"),
    ("Store", "urban decentral water tanks", "urban decentral water tanks"): (
        "src",
        "dst",
    ),
    ("Store", "EV battery", "EV battery"): ("src", "dst"),
    ("Store", "co2 sequestered", "co2 sequestered"): ("src", "dst"),
}


class SankeyChart(ESMChart):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.location = self._df.index.unique(DM.LOCATION).item()
        self.year = self._df.index.unique(DM.YEAR).item()
        self._df = self._df.droplevel(DM.YEAR).droplevel(DM.LOCATION)
        self._df.columns = ["value"]

    @staticmethod
    def _add_source_target_columns(idx: tuple):
        return pd.Series(mapping.get(idx))

    def plot(self):
        # Concatenate the data with source and target columns

        _df = self._df.copy()  # preserve original
        _df["index"] = _df.index  # convert to tuples
        _df[["source", "target"]] = _df["index"].apply(self._add_source_target_columns)
        label_mapping = {
            label: i
            for i, label in enumerate(set(pd.concat([_df["source"], _df["target"]])))
        }
        _df["source_id"] = _df["source"].map(label_mapping)
        _df["target_id"] = _df["target"].map(label_mapping)

        # df_agg = _df.groupby(["source", "target"]).sum()

        _df = _df.drop(columns=["index"])

        # _df = df_mapping.merge(df._data, how="left", left_index=True, right_on="variable")
        label_mapping = {
            label: i
            for i, label in enumerate(set(pd.concat([_df["source"], _df["target"]])))
        }
        _df = _df.replace(label_mapping)

        # def get_carrier_color(s) -> str:
        #     carrier = re.findall(r"\|AC", s)[0].strip("|")
        #     color_map = {"AC": COLOUR.red}
        #     return color_map[carrier]

        _df["color"] = _df.index.get_level_values("bus_carrier").map(COLOUR_SCHEME)
        # _df["label"] =
        _df["color"] = _df["color"].replace(
            np.nan, COLOUR.red
        )  # todo: map all colors and remove me

        self.fig = Figure(
            data=[
                Sankey(
                    valuesuffix=self.unit,
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
                        source=_df.source_id,
                        target=_df.target_id,
                        value=_df.value,
                        color=_df.color,
                        hovertemplate='"%{source.label}" to "%{target.label}": %{value}<extra></extra> ',
                    ),
                )
            ]
        )

        self._set_base_layout()
        self._style_title_and_legend_and_xaxis_label()
        self._append_footnotes()

        plotly.io.show(self.fig)  # todo: remove debugging


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
