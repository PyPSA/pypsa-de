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

# idx = df_plot.droplevel("year").droplevel("location").index.drop_duplicates()
# dict.fromkeys(idx, ("", ""))
mapping = {
    ("Generator", "offwind-ac", "AC"): ("", ""),
    ("Generator", "onwind", "AC"): ("", ""),
    ("Generator", "ror", "AC"): ("", ""),
    ("Generator", "rural heat vent", "rural heat"): ("", ""),
    ("Generator", "rural solar thermal", "rural heat"): ("", ""),
    ("Generator", "solar", "AC"): ("", ""),
    ("Generator", "solar rooftop", "low voltage"): ("", ""),
    ("Generator", "solar-hsat", "AC"): ("", ""),
    ("Generator", "unsustainable solid biomass", "solid biomass"): ("", ""),
    ("Generator", "urban decentral heat vent", "urban decentral heat"): ("", ""),
    ("Generator", "urban decentral solar thermal", "urban decentral heat"): ("", ""),
    ("Generator", "production gas", "gas"): ("", ""),
    ("Generator", "unsustainable biogas", "biogas"): ("", ""),
    ("Generator", "unsustainable bioliquids", "unsustainable bioliquids"): ("", ""),
    ("Generator", "urban central heat vent", "urban central heat"): ("", ""),
    ("Generator", "urban central solar thermal", "urban central heat"): ("", ""),
    ("Generator", "import H2", "H2"): ("", ""),
    ("Generator", "lng gas", "gas"): ("", ""),
    ("Generator", "offwind-dc", "AC"): ("", ""),
    ("Generator", "pipeline gas", "gas"): ("", ""),
    ("Generator", "biogas", "biogas"): ("", ""),
    ("Generator", "coal", "coal"): ("", ""),
    ("Generator", "import NH3", "NH3"): ("", ""),
    ("Generator", "lignite", "lignite"): ("", ""),
    ("Generator", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Generator", "oil primary", "oil primary"): ("", ""),
    ("Generator", "solid biomass", "solid biomass"): ("", ""),
    ("Generator", "uranium", "uranium"): ("", ""),
    ("Line", "Export Foreign", "AC"): ("", ""),
    ("Line", "Import Foreign", "AC"): ("", ""),
    ("Line", "Export Domestic", "AC"): ("", ""),
    ("Line", "Import Domestic", "AC"): ("", ""),
    ("Link", "BEV charger", "EV battery"): ("", ""),
    ("Link", "BEV charger", "low voltage"): ("", ""),
    ("Link", "BioSNG", "co2"): ("", ""),
    ("Link", "BioSNG", "gas"): ("", ""),
    ("Link", "BioSNG", "solid biomass"): ("", ""),
    ("Link", "BioSNG CC", "co2"): ("", ""),
    ("Link", "BioSNG CC", "co2 stored"): ("", ""),
    ("Link", "BioSNG CC", "gas"): ("", ""),
    ("Link", "BioSNG CC", "solid biomass"): ("", ""),
    ("Link", "CCGT methanol", "AC"): ("", ""),
    ("Link", "CCGT methanol", "co2"): ("", ""),
    ("Link", "CCGT methanol", "methanol"): ("", ""),
    ("Link", "CCGT methanol CC", "AC"): ("", ""),
    ("Link", "CCGT methanol CC", "co2"): ("", ""),
    ("Link", "CCGT methanol CC", "co2 stored"): ("", ""),
    ("Link", "CCGT methanol CC", "methanol"): ("", ""),
    ("Link", "DAC", "AC"): ("", ""),
    ("Link", "DAC", "co2"): ("", ""),
    ("Link", "DAC", "co2 stored"): ("", ""),
    ("Link", "DAC", "urban decentral heat"): ("", ""),
    ("Link", "Export Foreign", "co2 stored"): ("", ""),
    ("Link", "Export Foreign", "gas"): ("", ""),
    ("Link", "Export Foreign", "municipal solid waste"): ("", ""),
    ("Link", "Export Foreign", "solid biomass"): ("", ""),
    ("Link", "Fischer-Tropsch", "H2"): ("", ""),
    ("Link", "Fischer-Tropsch", "co2 stored"): ("", ""),
    ("Link", "Fischer-Tropsch", "oil"): ("", ""),
    ("Link", "H2 Fuel Cell", "AC"): ("", ""),
    ("Link", "H2 Fuel Cell", "H2"): ("", ""),
    ("Link", "HVC to air", "co2"): ("", ""),
    ("Link", "HVC to air", "non-sequestered HVC"): ("", ""),
    ("Link", "Haber-Bosch", "AC"): ("", ""),
    ("Link", "Haber-Bosch", "H2"): ("", ""),
    ("Link", "Haber-Bosch", "NH3"): ("", ""),
    ("Link", "Import Foreign", "co2 stored"): ("", ""),
    ("Link", "Import Foreign", "gas"): ("", ""),
    ("Link", "Import Foreign", "municipal solid waste"): ("", ""),
    ("Link", "Import Foreign", "solid biomass"): ("", ""),
    ("Link", "Methanol steam reforming", "H2"): ("", ""),
    ("Link", "Methanol steam reforming", "co2"): ("", ""),
    ("Link", "Methanol steam reforming", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming CC", "H2"): ("", ""),
    ("Link", "Methanol steam reforming CC", "co2"): ("", ""),
    ("Link", "Methanol steam reforming CC", "co2 stored"): ("", ""),
    ("Link", "Methanol steam reforming CC", "methanol"): ("", ""),
    ("Link", "OCGT", "AC"): ("", ""),
    ("Link", "OCGT", "co2"): ("", ""),
    ("Link", "OCGT", "gas"): ("", ""),
    ("Link", "OCGT methanol", "AC"): ("", ""),
    ("Link", "OCGT methanol", "co2"): ("", ""),
    ("Link", "OCGT methanol", "methanol"): ("", ""),
    ("Link", "SMR", "H2"): ("", ""),
    ("Link", "SMR", "co2"): ("", ""),
    ("Link", "SMR", "gas"): ("", ""),
    ("Link", "SMR CC", "H2"): ("", ""),
    ("Link", "SMR CC", "co2"): ("", ""),
    ("Link", "SMR CC", "co2 stored"): ("", ""),
    ("Link", "SMR CC", "gas"): ("", ""),
    ("Link", "Sabatier", "H2"): ("", ""),
    ("Link", "Sabatier", "co2 stored"): ("", ""),
    ("Link", "Sabatier", "gas"): ("", ""),
    ("Link", "agriculture machinery oil", "agriculture machinery oil"): ("", ""),
    ("Link", "agriculture machinery oil", "co2"): ("", ""),
    ("Link", "agriculture machinery oil", "oil"): ("", ""),
    ("Link", "allam methanol", "AC"): ("", ""),
    ("Link", "allam methanol", "co2"): ("", ""),
    ("Link", "allam methanol", "co2 stored"): ("", ""),
    ("Link", "allam methanol", "methanol"): ("", ""),
    ("Link", "ammonia cracker", "H2"): ("", ""),
    ("Link", "ammonia cracker", "NH3"): ("", ""),
    ("Link", "battery charger", "AC"): ("", ""),
    ("Link", "battery charger", "battery"): ("", ""),
    ("Link", "battery discharger", "AC"): ("", ""),
    ("Link", "battery discharger", "battery"): ("", ""),
    ("Link", "biomass to liquid", "co2"): ("", ""),
    ("Link", "biomass to liquid", "oil"): ("", ""),
    ("Link", "biomass to liquid", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid CC", "co2"): ("", ""),
    ("Link", "biomass to liquid CC", "co2 stored"): ("", ""),
    ("Link", "biomass to liquid CC", "oil"): ("", ""),
    ("Link", "biomass to liquid CC", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol", "co2"): ("", ""),
    ("Link", "biomass-to-methanol", "methanol"): ("", ""),
    ("Link", "biomass-to-methanol", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol CC", "co2"): ("", ""),
    ("Link", "biomass-to-methanol CC", "co2 stored"): ("", ""),
    ("Link", "biomass-to-methanol CC", "methanol"): ("", ""),
    ("Link", "biomass-to-methanol CC", "solid biomass"): ("", ""),
    ("Link", "coal for industry", "co2"): ("", ""),
    ("Link", "coal for industry", "coal"): ("", ""),
    ("Link", "coal for industry", "coal for industry"): ("", ""),
    ("Link", "electricity distribution grid", "AC"): ("", ""),
    ("Link", "electricity distribution grid", "low voltage"): ("", ""),
    ("Link", "electrobiofuels", "H2"): ("", ""),
    ("Link", "electrobiofuels", "co2"): ("", ""),
    ("Link", "electrobiofuels", "oil"): ("", ""),
    ("Link", "electrobiofuels", "solid biomass"): ("", ""),
    ("Link", "gas for industry", "co2"): ("", ""),
    ("Link", "gas for industry", "gas"): ("", ""),
    ("Link", "gas for industry", "gas for industry"): ("", ""),
    ("Link", "gas for industry CC", "co2"): ("", ""),
    ("Link", "gas for industry CC", "co2 stored"): ("", ""),
    ("Link", "gas for industry CC", "gas"): ("", ""),
    ("Link", "gas for industry CC", "gas for industry"): ("", ""),
    ("Link", "home battery charger", "home battery"): ("", ""),
    ("Link", "home battery charger", "low voltage"): ("", ""),
    ("Link", "home battery discharger", "home battery"): ("", ""),
    ("Link", "home battery discharger", "low voltage"): ("", ""),
    ("Link", "kerosene for aviation", "co2"): ("", ""),
    ("Link", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Link", "kerosene for aviation", "oil"): ("", ""),
    ("Link", "land transport oil", "co2"): ("", ""),
    ("Link", "land transport oil", "land transport oil"): ("", ""),
    ("Link", "land transport oil", "oil"): ("", ""),
    ("Link", "methanol-to-kerosene", "H2"): ("", ""),
    ("Link", "methanol-to-kerosene", "co2"): ("", ""),
    ("Link", "methanol-to-kerosene", "kerosene for aviation"): ("", ""),
    ("Link", "methanol-to-kerosene", "methanol"): ("", ""),
    ("Link", "methanolisation", "AC"): ("", ""),
    ("Link", "methanolisation", "H2"): ("", ""),
    ("Link", "methanolisation", "co2 stored"): ("", ""),
    ("Link", "methanolisation", "methanol"): ("", ""),
    ("Link", "municipal solid waste", "co2"): ("", ""),
    ("Link", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Link", "municipal solid waste", "non-sequestered HVC"): ("", ""),
    ("Link", "naphtha for industry", "naphtha for industry"): ("", ""),
    ("Link", "naphtha for industry", "oil"): ("", ""),
    ("Link", "naphtha for industry", "process emissions"): ("", ""),
    ("Link", "oil", "AC"): ("", ""),
    ("Link", "oil", "co2"): ("", ""),
    ("Link", "oil", "oil"): ("", ""),
    ("Link", "process emissions", "co2"): ("", ""),
    ("Link", "process emissions", "process emissions"): ("", ""),
    ("Link", "process emissions CC", "co2"): ("", ""),
    ("Link", "process emissions CC", "co2 stored"): ("", ""),
    ("Link", "process emissions CC", "process emissions"): ("", ""),
    ("Link", "rural air heat pump", "low voltage"): ("", ""),
    ("Link", "rural air heat pump", "rural heat"): ("", ""),
    ("Link", "rural biomass boiler", "rural heat"): ("", ""),
    ("Link", "rural biomass boiler", "solid biomass"): ("", ""),
    ("Link", "rural gas boiler", "co2"): ("", ""),
    ("Link", "rural gas boiler", "gas"): ("", ""),
    ("Link", "rural gas boiler", "rural heat"): ("", ""),
    ("Link", "rural ground heat pump", "low voltage"): ("", ""),
    ("Link", "rural ground heat pump", "rural heat"): ("", ""),
    ("Link", "rural resistive heater", "low voltage"): ("", ""),
    ("Link", "rural resistive heater", "rural heat"): ("", ""),
    ("Link", "rural water tanks charger", "rural heat"): ("", ""),
    ("Link", "rural water tanks charger", "rural water tanks"): ("", ""),
    ("Link", "rural water tanks discharger", "rural heat"): ("", ""),
    ("Link", "rural water tanks discharger", "rural water tanks"): ("", ""),
    ("Link", "shipping oil", "co2"): ("", ""),
    ("Link", "shipping oil", "oil"): ("", ""),
    ("Link", "shipping oil", "shipping oil"): ("", ""),
    ("Link", "solid biomass for industry", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry", "solid biomass for industry"): ("", ""),
    ("Link", "solid biomass for industry CC", "co2"): ("", ""),
    ("Link", "solid biomass for industry CC", "co2 stored"): ("", ""),
    ("Link", "solid biomass for industry CC", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry CC", "solid biomass for industry"): ("", ""),
    ("Link", "solid biomass to hydrogen", "H2"): ("", ""),
    ("Link", "solid biomass to hydrogen", "co2"): ("", ""),
    ("Link", "solid biomass to hydrogen", "co2 stored"): ("", ""),
    ("Link", "solid biomass to hydrogen", "solid biomass"): ("", ""),
    ("Link", "urban decentral air heat pump", "low voltage"): ("", ""),
    ("Link", "urban decentral air heat pump", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral biomass boiler", "solid biomass"): ("", ""),
    ("Link", "urban decentral biomass boiler", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral gas boiler", "co2"): ("", ""),
    ("Link", "urban decentral gas boiler", "gas"): ("", ""),
    ("Link", "urban decentral gas boiler", "urban decentral heat"): ("", ""),
    ("Link", "urban decentral resistive heater", "low voltage"): ("", ""),
    ("Link", "urban decentral resistive heater", "urban decentral heat"): ("", ""),
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
    ("Link", "waste CHP", "co2"): ("", ""),
    ("Link", "waste CHP", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP CC", "AC"): ("", ""),
    ("Link", "waste CHP CC", "co2"): ("", ""),
    ("Link", "waste CHP CC", "co2 stored"): ("", ""),
    ("Link", "waste CHP CC", "non-sequestered HVC"): ("", ""),
    ("Link", "CCGT", "AC"): ("", ""),
    ("Link", "CCGT", "co2"): ("", ""),
    ("Link", "CCGT", "gas"): ("", ""),
    ("Link", "DAC", "urban central heat"): ("", ""),
    ("Link", "Fischer-Tropsch", "urban central heat"): ("", ""),
    ("Link", "H2 Fuel Cell", "urban central heat"): ("", ""),
    ("Link", "Haber-Bosch", "urban central heat"): ("", ""),
    ("Link", "Sabatier", "urban central heat"): ("", ""),
    ("Link", "biogas to gas", "biogas"): ("", ""),
    ("Link", "biogas to gas", "co2"): ("", ""),
    ("Link", "biogas to gas", "gas"): ("", ""),
    ("Link", "biogas to gas CC", "biogas"): ("", ""),
    ("Link", "biogas to gas CC", "co2"): ("", ""),
    ("Link", "biogas to gas CC", "co2 stored"): ("", ""),
    ("Link", "biogas to gas CC", "gas"): ("", ""),
    ("Link", "industry methanol", "co2"): ("", ""),
    ("Link", "industry methanol", "industry methanol"): ("", ""),
    ("Link", "industry methanol", "methanol"): ("", ""),
    ("Link", "methanolisation", "urban central heat"): ("", ""),
    ("Link", "rural oil boiler", "co2"): ("", ""),
    ("Link", "rural oil boiler", "oil"): ("", ""),
    ("Link", "rural oil boiler", "rural heat"): ("", ""),
    ("Link", "unsustainable bioliquids", "co2"): ("", ""),
    ("Link", "unsustainable bioliquids", "oil"): ("", ""),
    ("Link", "unsustainable bioliquids", "unsustainable bioliquids"): ("", ""),
    ("Link", "urban central air heat pump", "low voltage"): ("", ""),
    ("Link", "urban central air heat pump", "urban central heat"): ("", ""),
    ("Link", "urban central coal CHP", "AC"): ("", ""),
    ("Link", "urban central coal CHP", "co2"): ("", ""),
    ("Link", "urban central coal CHP", "coal"): ("", ""),
    ("Link", "urban central coal CHP", "urban central heat"): ("", ""),
    ("Link", "urban central gas CHP", "AC"): ("", ""),
    ("Link", "urban central gas CHP", "co2"): ("", ""),
    ("Link", "urban central gas CHP", "gas"): ("", ""),
    ("Link", "urban central gas CHP", "urban central heat"): ("", ""),
    ("Link", "urban central gas CHP CC", "AC"): ("", ""),
    ("Link", "urban central gas CHP CC", "co2"): ("", ""),
    ("Link", "urban central gas CHP CC", "co2 stored"): ("", ""),
    ("Link", "urban central gas CHP CC", "gas"): ("", ""),
    ("Link", "urban central gas CHP CC", "urban central heat"): ("", ""),
    ("Link", "urban central gas boiler", "co2"): ("", ""),
    ("Link", "urban central gas boiler", "gas"): ("", ""),
    ("Link", "urban central gas boiler", "urban central heat"): ("", ""),
    ("Link", "urban central ptes heat pump", "low voltage"): ("", ""),
    ("Link", "urban central ptes heat pump", "urban central heat"): ("", ""),
    ("Link", "urban central resistive heater", "low voltage"): ("", ""),
    ("Link", "urban central resistive heater", "urban central heat"): ("", ""),
    ("Link", "urban central solid biomass CHP", "AC"): ("", ""),
    ("Link", "urban central solid biomass CHP", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP", "urban central heat"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "AC"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "co2"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "co2 stored"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "urban central heat"): ("", ""),
    ("Link", "urban central water pits charger", "urban central heat"): ("", ""),
    ("Link", "urban central water pits charger", "urban central water pits"): ("", ""),
    ("Link", "urban central water pits discharger", "urban central heat"): ("", ""),
    ("Link", "urban central water pits discharger", "urban central water pits"): (
        "",
        "",
    ),
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
    ("Link", "urban decentral oil boiler", "co2"): ("", ""),
    ("Link", "urban decentral oil boiler", "oil"): ("", ""),
    ("Link", "urban decentral oil boiler", "urban decentral heat"): ("", ""),
    ("Link", "waste CHP", "urban central heat"): ("", ""),
    ("Link", "waste CHP CC", "urban central heat"): ("", ""),
    ("Link", "Export Domestic", "co2 stored"): ("", ""),
    ("Link", "Export Domestic", "gas"): ("", ""),
    ("Link", "Export Domestic", "municipal solid waste"): ("", ""),
    ("Link", "Export Domestic", "solid biomass"): ("", ""),
    ("Link", "Import Domestic", "co2 stored"): ("", ""),
    ("Link", "Import Domestic", "gas"): ("", ""),
    ("Link", "Import Domestic", "municipal solid waste"): ("", ""),
    ("Link", "Import Domestic", "solid biomass"): ("", ""),
    ("Link", "coal", "AC"): ("", ""),
    ("Link", "coal", "co2"): ("", ""),
    ("Link", "coal", "coal"): ("", ""),
    ("Link", "nuclear", "AC"): ("", ""),
    ("Link", "nuclear", "uranium"): ("", ""),
    ("Link", "solid biomass", "AC"): ("", ""),
    ("Link", "solid biomass", "solid biomass"): ("", ""),
    ("Link", "urban central oil CHP", "AC"): ("", ""),
    ("Link", "urban central oil CHP", "co2"): ("", ""),
    ("Link", "urban central oil CHP", "oil"): ("", ""),
    ("Link", "urban central oil CHP", "urban central heat"): ("", ""),
    ("Link", "Export Foreign", "AC"): ("", ""),
    ("Link", "Import Foreign", "AC"): ("", ""),
    ("Link", "import gas", "co2"): ("", ""),
    ("Link", "import gas", "gas"): ("", ""),
    ("Link", "lignite", "AC"): ("", ""),
    ("Link", "lignite", "co2"): ("", ""),
    ("Link", "lignite", "lignite"): ("", ""),
    ("Link", "urban central lignite CHP", "AC"): ("", ""),
    ("Link", "urban central lignite CHP", "co2"): ("", ""),
    ("Link", "urban central lignite CHP", "lignite"): ("", ""),
    ("Link", "urban central lignite CHP", "urban central heat"): ("", ""),
    ("Link", "Export Domestic", "AC"): ("", ""),
    ("Link", "Import Domestic", "AC"): ("", ""),
    ("Link", "import methanol", "co2"): ("", ""),
    ("Link", "import methanol", "methanol"): ("", ""),
    ("Link", "import oil", "co2"): ("", ""),
    ("Link", "import oil", "oil"): ("", ""),
    ("Link", "oil refining", "co2"): ("", ""),
    ("Link", "oil refining", "oil"): ("", ""),
    ("Link", "oil refining", "oil primary"): ("", ""),
    ("Link", "Export Foreign", "H2"): ("", ""),
    ("Link", "H2 Electrolysis", "AC"): ("", ""),
    ("Link", "H2 Electrolysis", "H2"): ("", ""),
    ("Link", "H2 Electrolysis", "urban central heat"): ("", ""),
    ("Link", "Import Foreign", "H2"): ("", ""),
    ("Link", "V2G", "EV battery"): ("", ""),
    ("Link", "V2G", "low voltage"): ("", ""),
    ("Link", "shipping methanol", "co2"): ("", ""),
    ("Link", "shipping methanol", "methanol"): ("", ""),
    ("Link", "shipping methanol", "shipping methanol"): ("", ""),
    ("Link", "H2 OCGT", "AC"): ("", ""),
    ("Link", "H2 OCGT", "H2"): ("", ""),
    ("Link", "urban central H2 CHP", "AC"): ("", ""),
    ("Link", "urban central H2 CHP", "H2"): ("", ""),
    ("Link", "urban central H2 CHP", "urban central heat"): ("", ""),
    ("Link", "Export Domestic", "H2"): ("", ""),
    ("Link", "Import Domestic", "H2"): ("", ""),
    ("Link", "co2 sequestered", "co2 sequestered"): ("", ""),
    ("Link", "co2 sequestered", "co2 stored"): ("", ""),
    ("Load", "agriculture electricity", "low voltage"): ("", ""),
    ("Load", "agriculture heat", "rural heat"): ("", ""),
    ("Load", "agriculture machinery oil", "agriculture machinery oil"): ("", ""),
    ("Load", "coal for industry", "coal for industry"): ("", ""),
    ("Load", "electricity", "low voltage"): ("", ""),
    ("Load", "gas for industry", "gas for industry"): ("", ""),
    ("Load", "industry electricity", "low voltage"): ("", ""),
    ("Load", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Load", "land transport EV", "EV battery"): ("", ""),
    ("Load", "land transport fuel cell", "H2"): ("", ""),
    ("Load", "land transport oil", "land transport oil"): ("", ""),
    ("Load", "low-temperature heat for industry", "urban decentral heat"): ("", ""),
    ("Load", "naphtha for industry", "naphtha for industry"): ("", ""),
    ("Load", "process emissions", "process emissions"): ("", ""),
    ("Load", "rural heat", "rural heat"): ("", ""),
    ("Load", "shipping oil", "shipping oil"): ("", ""),
    ("Load", "solid biomass for industry", "solid biomass for industry"): ("", ""),
    ("Load", "urban decentral heat", "urban decentral heat"): ("", ""),
    ("Load", "H2 for industry", "H2"): ("", ""),
    ("Load", "industry methanol", "industry methanol"): ("", ""),
    ("Load", "low-temperature heat for industry", "urban central heat"): ("", ""),
    ("Load", "urban central heat", "urban central heat"): ("", ""),
    ("Load", "NH3", "NH3"): ("", ""),
    ("Load", "shipping methanol", "shipping methanol"): ("", ""),
    ("StorageUnit", "hydro", "AC"): ("", ""),
    ("StorageUnit", "PHS", "AC"): ("", ""),
    ("Store", "H2 Store", "H2"): ("", ""),
    ("Store", "battery", "battery"): ("", ""),
    ("Store", "co2 stored", "co2 stored"): ("", ""),
    ("Store", "gas", "gas"): ("", ""),
    ("Store", "home battery", "home battery"): ("", ""),
    ("Store", "non-sequestered HVC", "non-sequestered HVC"): ("", ""),
    ("Store", "rural water tanks", "rural water tanks"): ("", ""),
    ("Store", "urban decentral water tanks", "urban decentral water tanks"): ("", ""),
    ("Store", "urban central water pits", "urban central water pits"): ("", ""),
    ("Store", "urban central water tanks", "urban central water tanks"): ("", ""),
    ("Store", "ammonia store", "NH3"): ("", ""),
    ("Store", "co2", "co2"): ("", ""),
    ("Store", "coal", "coal"): ("", ""),
    ("Store", "lignite", "lignite"): ("", ""),
    ("Store", "methanol", "methanol"): ("", ""),
    ("Store", "oil", "oil"): ("", ""),
    ("Store", "uranium", "uranium"): ("", ""),
    ("Store", "EV battery", "EV battery"): ("", ""),
    ("Store", "co2 sequestered", "co2 sequestered"): ("", ""),
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
