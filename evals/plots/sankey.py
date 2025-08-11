# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for Sankey diagram."""

import pandas as pd
from plotly.graph_objs import Figure, Sankey

from evals.constants import COLOUR, RUN_META_DATA
from evals.constants import DataModel as DM
from evals.plots._base import ESMChart
from evals.utils import (
    filter_by,
    prettify_number,
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


LINK_MAPPING = {
    # ("Link", "BEV charger", "EV battery"): ("Secondary AC Out", "Transport"),
    ("Link", "BEV charger", "low voltage losses"): ("Transport", "Losses Transport"),
    # Wood gasification treated as primary energy to simplify energy flows
    ("Link", "BioSNG", "gas"): ("Solid Biomass", "Primary Gas"),
    ("Link", "BioSNG CC", "gas"): ("Solid Biomass", "Primary Gas"),
    # ("Link", "BioSNG", "solid biomass"): ("", ""),
    ("Link", "BioSNG", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),  # todo: need to subtract from Biomass primary
    # ("Link", "BioSNG CC", "solid biomass"): ("", ""),
    ("Link", "BioSNG CC", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "CCGT", "AC"): ("Secondary Gas In", "Secondary AC Out"),
    # ("Link", "CCGT", "gas"): ("", ""),
    ("Link", "CCGT", "gas losses"): ("Secondary Gas In", "Transformation Losses"),
    ("Link", "CCGT methanol", "AC"): ("Secondary Liquids In", "Secondary AC Out"),
    # ("Link", "CCGT methanol", "methanol"): ("", ""),
    ("Link", "CCGT methanol", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    ("Link", "CCGT methanol CC", "AC"): ("Secondary Liquids In", "Secondary AC Out"),
    # ("Link", "CCGT methanol CC", "methanol"): ("", ""),
    ("Link", "CCGT methanol CC", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    ("Link", "DAC", "AC"): ("", ""),
    ("Link", "DAC", "urban central heat"): ("", ""),
    ("Link", "DAC", "urban decentral heat"): ("", ""),
    # ("Store", "EV battery", "EV battery"): ("", ""),
    # ("Link", "Fischer-Tropsch", "H2"): ("", ""),
    ("Link", "Fischer-Tropsch", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
    ),
    ("Link", "Fischer-Tropsch", "oil"): ("Secondary H2 In", "Secondary Liquids Out"),
    ("Link", "Fischer-Tropsch", "urban central heat"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    ("Link", "H2 Electrolysis", "AC losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    ("Link", "H2 Electrolysis", "H2"): ("Secondary AC In", "Secondary H2 Out"),
    ("Link", "H2 Electrolysis", "urban central heat"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    ("Link", "H2 Fuel Cell", "AC"): ("Secondary H2 In", "Secondary AC Out"),
    # ("Link", "H2 Fuel Cell", "H2"): ("", ""),
    ("Link", "H2 Fuel Cell", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
    ),
    ("Link", "H2 Fuel Cell", "urban central heat"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    ("Link", "H2 OCGT", "AC"): ("Secondary H2 In", "Secondary AC Out"),
    # ("Link", "H2 OCGT", "H2"): ("", ""),
    ("Link", "H2 OCGT", "H2 losses"): ("Secondary H2 In", "Transformation Losses"),
    # ("Store", "H2 Store", "H2"): ("", ""),
    ("Store", "H2 Store supply", "H2"): ("Secondary H2 Out", "Secondary H2 In"),
    # ("Store", "H2 Store demand", "H2"): ("Secondary H2 In", "Secondary H2 Out"),
    ("Load", "H2 for industry", "H2"): ("Secondary H2 Out", "Industry"),
    ("Link", "HVC to air", "non-sequestered HVC"): ("", ""),
    # ("Link", "Haber-Bosch", "AC"): ("", ""),
    ("Link", "Haber-Bosch", "AC losses"): ("Secondary AC In", "Transformation Losses"),
    # ("Link", "Haber-Bosch", "H2"): ("", ""),
    ("Link", "Haber-Bosch", "H2 losses"): ("Secondary H2 In", "Transformation Losses"),
    # ("Link", "Haber-Bosch", "NH3"): ("", ""),
    ("Link", "Haber-Bosch", "NH3 from AC"): (
        "Secondary AC In",
        "Secondary Liquids Out",
    ),
    ("Link", "Haber-Bosch", "NH3 from H2"): (
        "Secondary H2 In",
        "Secondary Liquids Out",
    ),
    # ("Link", "Haber-Bosch", "urban central heat"): ("", ""),
    ("Link", "Haber-Bosch", "urban central heat from AC"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    ("Link", "Haber-Bosch", "urban central heat from H2"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    ("Link", "Methanol steam reforming", "H2"): (
        "Secondary Liquids In",
        "Secondary H2 Out",
    ),
    # ("Link", "Methanol steam reforming", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    ("Link", "Methanol steam reforming CC", "H2"): (
        "Secondary Liquids In",
        "Secondary H2 Out",
    ),
    # ("Link", "Methanol steam reforming CC", "methanol"): ("", ""),
    ("Link", "Methanol steam reforming CC", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    ("Load", "NH3", "NH3"): ("", ""),
    ("Link", "OCGT", "AC"): ("Secondary Gas In", "Secondary AC Out"),
    # ("Link", "OCGT", "gas"): ("", ""),
    ("Link", "OCGT", "gas losses"): ("Secondary Gas In", "Transformation Losses"),
    ("Link", "OCGT methanol", "AC"): ("Secondary Liquids In", "Secondary AC Out"),
    ("Link", "OCGT methanol", "methanol"): ("", ""),
    ("Link", "OCGT methanol", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    ("StorageUnit", "PHS demand", "AC"): ("Secondary AC In", "AC Storage"),
    ("StorageUnit", "PHS supply", "AC"): ("AC Storage", "Secondary AC Out"),
    # todo: PHS losses
    ("Link", "SMR", "H2"): ("Secondary Gas In", "Secondary H2 Out"),
    # ("Link", "SMR", "gas"): ("", ""),
    ("Link", "SMR", "gas losses"): ("Secondary Gas In", "Transformation Losses"),
    ("Link", "SMR CC", "H2"): ("Secondary Gas In", "Secondary H2 Out"),
    # ("Link", "SMR CC", "gas"): ("", ""),
    ("Link", "SMR CC", "gas losses"): ("Secondary Gas In", "Transformation Losses"),
    # ("Link", "Sabatier", "H2"): ("", ""),
    ("Link", "Sabatier", "H2 losses"): ("Secondary H2 In", "Transformation Losses"),
    ("Link", "Sabatier", "gas"): ("Secondary H2 In", "Secondary Gas Out"),
    ("Link", "Sabatier", "urban central heat"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    # todo: treat V2G as a storage technology same as battery
    ("Link", "V2G", "EV battery"): ("", ""),
    # ("Link", "V2G", "EV battery losses"): ("Transport", "Losses"),
    # ("Link", "V2G", "low voltage"): ("Transport", "Secondary AC In"),
    ("Load", "agriculture electricity", "low voltage"): (
        "Secondary AC Out",
        "Agriculture",
    ),
    ("Load", "agriculture heat", "rural heat"): ("", ""),
    # prefer Load for tests to find unexpected losses in Links
    # ("Link", "agriculture machinery oil", "agriculture machinery oil"): ("Secondary Liquids Out", "Agriculture"),
    # ("Link", "agriculture machinery oil", "oil"): ("", ""),
    ("Load", "agriculture machinery oil", "agriculture machinery oil"): (
        "Secondary Liquids Out",
        "Agriculture",
    ),
    ("Link", "allam methanol", "AC"): ("Secondary Liquids In", "Secondary AC Out"),
    # ("Link", "allam methanol", "methanol"): ("", ""),
    ("Link", "allam methanol", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    ("Link", "ammonia cracker", "H2"): ("Secondary Liquids In", "Secondary H2 Out"),
    # ("Link", "ammonia cracker", "NH3"): ("", ""),
    ("Link", "ammonia cracker", "NH3 losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    # ("Store", "ammonia store", "NH3"): ("", ""),
    # ("Store", "battery", "battery"): ("", ""),
    # ("Store", "home battery", "home battery"): ("", ""),
    ("Link", "battery charger", "AC"): (
        "Secondary AC In",
        "AC Storage",
    ),  # including losses
    ("Link", "home battery charger", "low voltage"): ("Secondary AC In", "AC Storage"),
    ("Link", "battery charger", "AC losses"): ("AC Storage", "Storage Losses"),
    ("Link", "home battery charger", "low voltage losses"): (
        "AC Storage",
        "Storage Losses",
    ),
    # ("Link", "battery charger", "battery"): ("Secondary AC In", "AC Storage"),
    # ("Link", "home battery charger", "home battery"): ("", ""),
    ("Link", "battery discharger", "AC"): ("AC Storage", "Secondary AC Out"),
    ("Link", "home battery discharger", "low voltage"): (
        "AC Storage",
        "Secondary AC Out",
    ),
    # ("Link", "battery discharger", "battery"): ("", ""),
    # ("Link", "home battery discharger", "home battery"): ("", ""),
    ("Link", "battery discharger", "battery losses"): ("AC Storage", "Storage Losses"),
    ("Link", "home battery discharger", "home battery losses"): (
        "AC Storage",
        "Storage Losses",
    ),
    ("Link", "biogas to gas", "gas"): ("Biogas", "Primary Gas"),
    ("Link", "biogas to gas CC", "gas"): ("Biogas", "Primary Gas"),
    ("Link", "biomass to liquid", "oil"): (
        "Secondary Solids In",
        "Secondary Liquids Out",
    ),
    # ("Link", "biomass to liquid", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "biomass to liquid CC", "oil"): (
        "Secondary Solids In",
        "Secondary Liquids Out",
    ),
    # ("Link", "biomass to liquid CC", "solid biomass"): ("", ""),
    ("Link", "biomass to liquid CC", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "biomass-to-methanol", "methanol"): (
        "Secondary Solids In",
        "Secondary Liquids Out",
    ),
    # ("Link", "biomass-to-methanol", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "biomass-to-methanol CC", "methanol"): (
        "Secondary Solids In",
        "Secondary Liquids Out",
    ),
    # ("Link", "biomass-to-methanol CC", "solid biomass"): ("", ""),
    ("Link", "biomass-to-methanol CC", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Generator", "coal", "coal"): ("Coal", "Primary Solids"),
    ("Link", "coal", "AC"): ("Secondary Solids In", "Secondary AC Out"),
    # ("Link", "coal", "coal"): ("", ""),
    ("Link", "coal", "coal losses"): ("Secondary Solids In", "Transformation Losses"),
    # ("Store", "coal", "coal"): ("", ""),
    ("Load", "electricity", "low voltage"): ("Secondary AC Out", "Final AC"),
    ("Link", "electricity distribution grid", "losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    # ("Link", "electrobiofuels", "H2"): ("", ""),
    ("Link", "electrobiofuels", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
    ),
    # ("Link", "electrobiofuels", "oil"): ("", ""),
    ("Link", "electrobiofuels", "oil from solid biomass"): (
        "Secondary Solids In",
        "Secondary Liquids Out",
    ),
    ("Link", "electrobiofuels", "oil from H2"): (
        "Secondary H2 In",
        "Secondary Liquids Out",
    ),
    # ("Link", "electrobiofuels", "solid biomass"): ("", ""),
    ("Link", "electrobiofuels", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Store", "gas", "gas"): ("", ""),  # todo: gas supply/demand
    ("Store", "gas supply", "gas"): ("Secondary Gas Out", "Secondary Gas In"),
    # ("Store", "gas demand", "gas"): ("Secondary Gas In", "Secondary Gas Out"),
    # ("Link", "gas for industry", "gas"): ("", ""),
    ("Link", "gas for industry", "gas for industry"): ("Secondary Gas Out", "Industry"),
    # ("Load", "gas for industry", "gas for industry"): ("", ""),
    # ("Link", "gas for industry CC", "gas"): ("", ""),
    ("Link", "gas for industry CC", "gas for industry"): (
        "Secondary Gas Out",
        "Industry",
    ),
    ("Link", "gas for industry CC", "gas losses"): ("Industry", "Losses Industry"),
    ("StorageUnit", "hydro supply", "AC"): ("Hydro Power", "Primary AC"),  # is primary
    ("Generator", "import H2", "H2"): ("Green Hydrogen", "Primary H2"),
    ("Generator", "import NH3", "NH3"): ("Green Liquids", "Primary Liquids"),
    ("Link", "import gas", "gas"): ("Green Gas", "Primary Gas"),
    ("Link", "import methanol", "methanol"): ("Green Liquids", "Primary Liquids"),
    ("Link", "import oil", "oil"): ("Green Liquids", "Primary Liquids"),
    ("Load", "industry electricity", "low voltage"): ("Secondary AC Out", "Industry"),
    # ("Link", "industry methanol", "industry methanol"): ("", ""),
    # ("Link", "industry methanol", "methanol"): ("", ""),
    ("Load", "industry methanol", "industry methanol"): (
        "Secondary Liquids Out",
        "Industry",
    ),
    # ("Link", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    # ("Link", "kerosene for aviation", "oil"): ("", ""),
    ("Load", "kerosene for aviation", "kerosene for aviation"): (
        "Secondary Liquids Out",
        "Transport",
    ),
    ("Load", "land transport EV", "EV battery"): ("Secondary AC Out", "Transport"),
    ("Load", "land transport fuel cell", "H2"): ("Secondary H2 Out", "Transport"),
    ("Generator", "lignite", "lignite"): ("Coal", "Primary Solids"),
    ("Link", "lignite", "AC"): ("Secondary Solids In", "Secondary AC Out"),
    # ("Link", "lignite", "lignite"): ("", ""),
    ("Link", "lignite", "lignite losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    # ("Store", "lignite", "lignite"): ("", ""),
    ("Generator", "lng gas", "gas"): ("LNG", "Primary Gas"),
    ("Load", "low-temperature heat for industry", "urban central heat"): (
        "Secondary Heat Out",
        "Industry",
    ),
    ("Load", "low-temperature heat for industry", "urban decentral heat"): (
        "Decentral Heat",
        "Industry",
    ),
    # ("Store", "methanol", "methanol"): ("", ""),
    # ("Link", "methanol-to-kerosene", "H2"): ("", ""),
    ("Link", "methanol-to-kerosene", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
    ),
    # ("Link", "methanol-to-kerosene", "kerosene for aviation"): ("", ""),
    ("Link", "methanol-to-kerosene", "kerosene for aviation from H2"): (
        "Secondary H2 In",
        "Secondary Liquids Out",
    ),
    # ("Link", "methanol-to-kerosene", "kerosene for aviation from methanol"): ("", ""),
    # ("Link", "methanol-to-kerosene", "methanol"): ("", ""),
    # must include Liquids to Liquids conversion losses
    ("Link", "methanol-to-kerosene", "methanol losses"): (
        "Secondary Liquids In",
        "Transformation Losses",
    ),
    # ("Link", "methanolisation", "AC"): ("", ""),
    ("Link", "methanolisation", "AC losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    # ("Link", "methanolisation", "H2"): ("", ""),
    ("Link", "methanolisation", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
    ),
    # ("Link", "methanolisation", "methanol"): ("", ""),
    ("Link", "methanolisation", "methanol from AC"): (
        "Secondary AC In",
        "Secondary Liquids Out",
    ),
    ("Link", "methanolisation", "methanol from H2"): (
        "Secondary H2 In",
        "Secondary Liquids Out",
    ),
    # ("Link", "methanolisation", "urban central heat"): ("", ""),
    ("Link", "methanolisation", "urban central heat from AC"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    ("Link", "methanolisation", "urban central heat from H2"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    # ("Generator", "municipal solid waste", "municipal solid waste"): ("Waste", "Primary Solids"),
    # ("Link", "municipal solid waste", "municipal solid waste"): ("", ""),
    ("Link", "municipal solid waste", "non-sequestered HVC"): (
        "Waste",
        "Primary Solids",
    ),
    # ("Link", "naphtha for industry", "naphtha for industry"): ("", ""),
    # ("Link", "naphtha for industry", "oil"): ("", ""),
    ("Load", "naphtha for industry", "naphtha for industry"): (
        "Secondary Liquids Out",
        "Industry",
    ),
    ("Store", "non-sequestered HVC", "non-sequestered HVC"): (
        "Waste",
        "Primary Solids",
    ),  # todo: why is this store not balanced?
    ("Link", "nuclear", "AC"): ("Nuclear Power", "Primary AC"),
    ("Generator", "offwind-ac", "AC"): ("Wind Power", "Primary AC"),
    ("Generator", "offwind-dc", "AC"): ("Wind Power", "Primary AC"),
    ("Link", "oil", "AC"): ("Secondary Liquids In", "Secondary AC Out"),
    # ("Link", "oil", "oil"): ("", ""),
    ("Link", "oil", "oil losses"): ("Secondary Liquids In", "Transformation Losses"),
    # ("Store", "oil", "oil"): ("", ""),
    # ("Generator", "oil primary", "oil primary"): ("", ""),
    # skipping oil primary and refining losses
    ("Link", "oil refining", "oil"): ("Oil", "Primary Liquids"),
    # ("Link", "oil refining", "oil primary"): ("", ""),
    # ("Link", "oil refining", "oil primary losses"): ("", ""),
    ("Generator", "onwind", "AC"): ("Wind Power", "Primary AC"),
    ("Generator", "pipeline gas", "gas"): ("Pipeline", "Primary Gas"),
    ("Generator", "production gas", "gas"): ("Production", "Primary Gas"),
    ("Generator", "ror", "AC"): ("Hydro Power", "Primary AC"),
    ("Link", "rural air heat pump", "ambient heat"): (
        "Ambient Heat Decentral",
        "Decentral Heat",
    ),
    # ("Link", "rural air heat pump", "low voltage"): ("", ""),
    ("Link", "rural air heat pump", "rural heat"): (
        "Secondary AC In",
        "Decentral Heat",
    ),
    ("Link", "rural biomass boiler", "rural heat"): (
        "Secondary Solids In",
        "Decentral Heat",
    ),
    # ("Link", "rural biomass boiler", "solid biomass"): ("", ""),
    ("Link", "rural biomass boiler", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    # ("Link", "rural gas boiler", "gas"): ("", ""),
    ("Link", "rural gas boiler", "gas losses"): (
        "Secondary Gas In",
        "Transformation Losses",
    ),
    ("Link", "rural gas boiler", "rural heat"): ("Secondary Gas In", "Decentral Heat"),
    ("Link", "rural ground heat pump", "ambient heat"): (
        "Ambient Heat Decentral",
        "Decentral Heat",
    ),
    # ("Link", "rural ground heat pump", "low voltage"): ("", ""),
    ("Link", "rural ground heat pump", "rural heat"): (
        "Secondary AC In",
        "Decentral Heat",
    ),
    ("Load", "rural heat", "rural heat"): ("Decentral Heat", "HH & Services"),
    ("Generator", "rural heat vent", "rural heat"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    # ("Link", "rural resistive heater", "low voltage"): ("", ""),
    ("Link", "rural resistive heater", "low voltage losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    ("Link", "rural resistive heater", "rural heat"): (
        "Secondary AC In",
        "Decentral Heat",
    ),
    ("Generator", "rural solar thermal", "rural heat"): (
        "Solar Heat",
        "Decentral Heat",
    ),
    ("Store", "rural water tanks", "rural water tanks"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    # ("Link", "rural water tanks charger", "rural heat"): ("", ""),
    # ("Link", "rural water tanks charger", "rural water tanks"): ("", ""),
    # ("Link", "rural water tanks discharger", "rural heat"): ("", ""),
    # ("Link", "rural water tanks discharger", "rural water tanks"): ("", ""),
    # ("Link", "shipping methanol", "methanol"): ("", ""),
    # ("Link", "shipping methanol", "shipping methanol"): ("", ""),
    ("Load", "shipping methanol", "shipping methanol"): (
        "Secondary Liquids Out",
        "Transport",
    ),
    ("Generator", "solar", "AC"): ("Solar Power", "Primary AC"),
    ("Generator", "solar rooftop", "low voltage"): ("Solar Power", "Primary AC"),
    ("Generator", "solar-hsat", "AC"): ("Solar Power", "Primary AC"),
    ("Generator", "solid biomass", "solid biomass"): (
        "Solid Biomass",
        "Primary Solids",
    ),
    # ("Link", "solid biomass for industry", "solid biomass"): ("", ""),
    ("Link", "solid biomass for industry", "solid biomass for industry"): (
        "Secondary Solids Out",
        "Industry",
    ),
    # ("Load", "solid biomass for industry", "solid biomass for industry"): ("Secondary Solids Out", "Industry"),
    ("Link", "solid biomass for industry CC", "solid biomass"): (
        "Secondary Solids Out",
        "Industry",
    ),
    # ("Link", "solid biomass for industry CC", "solid biomass for industry"): ("Secondary Solids Out", "Industry"),
    ("Link", "solid biomass for industry CC", "solid biomass losses"): (
        "Industry",
        "Losses Industry",
    ),
    ("Link", "solid biomass to hydrogen", "H2"): (
        "Secondary Solids In",
        "Secondary H2 Out",
    ),
    # ("Link", "solid biomass to hydrogen", "solid biomass"): ("", ""),
    ("Link", "solid biomass to hydrogen", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "urban central H2 CHP", "AC"): ("Secondary H2 In", "Secondary AC Out"),
    ("Link", "urban central H2 CHP", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
    ),
    ("Link", "urban central H2 CHP", "urban central heat"): (
        "Secondary H2 In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central air heat pump", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    # ("Link", "urban central air heat pump", "low voltage"): ("", ""),
    ("Link", "urban central air heat pump", "urban central heat"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central coal CHP", "AC"): (
        "Secondary Solids In",
        "Secondary AC Out",
    ),
    ("Link", "urban central coal CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    ("Link", "urban central coal CHP", "urban central heat"): (
        "Secondary Solids In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas CHP", "AC"): ("Secondary Gas In", "Secondary AC Out"),
    ("Link", "urban central gas CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas CHP", "gas losses"): (
        "Secondary Gas In",
        "Transformation Losses",
    ),
    ("Link", "urban central gas CHP", "urban central heat"): (
        "Secondary Gas In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas CHP CC", "AC"): (
        "Secondary Gas In",
        "Secondary AC Out",
    ),
    ("Link", "urban central gas CHP CC", "gas losses"): (
        "Secondary Gas In",
        "Transformation Losses",
    ),
    ("Link", "urban central gas CHP CC", "urban central heat"): (
        "Secondary Gas In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central gas boiler", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    # ("Link", "urban central gas boiler", "gas"): ("", ""),
    ("Link", "urban central gas boiler", "urban central heat"): (
        "Secondary Gas In",
        "Secondary Heat Out",
    ),
    ("Load", "urban central heat", "urban central heat"): (
        "Secondary Heat Out",
        "HH & Services",
    ),
    ("Generator", "urban central heat vent", "urban central heat"): (
        "Secondary Heat Out",
        "Ventilation Losses",
    ),
    ("Link", "urban central ptes heat pump", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    ("Link", "urban central ptes heat pump", "low voltage"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    # ("Link", "urban central ptes heat pump", "urban central heat"): ("", ""),
    # ("Link", "urban central resistive heater", "low voltage"): ("", ""),,
    ("Link", "urban central resistive heater", "low voltage losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    # ("Link", "urban central resistive heater", "urban central heat"): ("Secondary AC In", "Secondary Heat Out"),
    ("Generator", "urban central solar thermal", "urban central heat"): (
        "Solar Heat",
        "Primary Heat",
    ),
    ("Link", "urban central solid biomass CHP", "AC"): (
        "Secondary Solids In",
        "Secondary AC Out",
    ),
    ("Link", "urban central solid biomass CHP", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    # ("Link", "urban central solid biomass CHP", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP", "urban central heat"): (
        "Secondary Solids In",
        "Secondary Heat Out",
    ),
    ("Link", "urban central solid biomass CHP CC", "AC"): (
        "Secondary Solids In",
        "Secondary AC Out",
    ),
    ("Link", "urban central solid biomass CHP CC", "ambient heat"): (
        "Ambient Heat Central",
        "Secondary Heat Out",
    ),
    # ("Link", "urban central solid biomass CHP CC", "solid biomass"): ("", ""),
    ("Link", "urban central solid biomass CHP CC", "urban central heat"): (
        "Secondary Solids In",
        "Secondary Heat Out",
    ),
    ("Store", "urban central water pits", "urban central water pits"): (
        "Secondary Heat In",
        "Transformation Losses",
    ),
    # ("Link", "urban central water pits charger", "urban central heat"): ("", ""),
    # ("Link", "urban central water pits charger", "urban central water pits"): ("", ""),
    # ("Link", "urban central water pits discharger", "urban central heat"): ("", ""),
    # ("Link", "urban central water pits discharger", "urban central water pits"): ("", ""),
    ("Store", "urban central water tanks", "urban central water tanks"): (
        "Secondary Heat In",
        "Transformation Losses",
    ),
    # ("Link", "urban central water tanks charger", "urban central heat"): ("", ""),
    # ("Link", "urban central water tanks charger", "urban central water tanks"): ("", ""),
    # ("Link", "urban central water tanks discharger", "urban central heat"): ("", ""),
    # ("Link", "urban central water tanks discharger", "urban central water tanks"): ("", ""),
    ("Link", "urban decentral air heat pump", "ambient heat"): (
        "Ambient Heat Decentral",
        "Decentral Heat",
    ),
    # ("Link", "urban decentral air heat pump", "low voltage"): ("", ""),
    ("Link", "urban decentral air heat pump", "urban decentral heat"): (
        "Secondary AC In",
        "Decentral Heat",
    ),
    # ("Link", "urban decentral biomass boiler", "solid biomass"): ("", ""),
    ("Link", "urban decentral biomass boiler", "solid biomass losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "urban decentral biomass boiler", "urban decentral heat"): (
        "Secondary Solids In",
        "Secondary Heat Out",
    ),
    # ("Link", "urban decentral gas boiler", "gas"): ("", ""),
    ("Link", "urban decentral gas boiler", "gas losses"): (
        "Secondary Gas In",
        "Transformation Losses",
    ),
    ("Link", "urban decentral gas boiler", "urban decentral heat"): (
        "Secondary Gas In",
        "Decentral Heat",
    ),
    ("Load", "urban decentral heat", "urban decentral heat"): ("", ""),
    ("Generator", "urban decentral heat vent", "urban decentral heat"): (
        "Decentral Heat",
        "HH & Services",
    ),
    # ("Link", "urban decentral resistive heater", "low voltage"): ("", ""),
    ("Link", "urban decentral resistive heater", "low voltage losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    ("Link", "urban decentral resistive heater", "urban decentral heat"): (
        "Secondary AC In",
        "Secondary Heat Out",
    ),
    ("Generator", "urban decentral solar thermal", "urban decentral heat"): (
        "Solar Heat",
        "Decentral Heat",
    ),
    ("Store", "urban decentral water tanks", "urban decentral water tanks"): (
        "Decentral Heat",
        "Decentral Heat Losses",
    ),
    # ("Link", "urban decentral water tanks charger", "urban decentral heat"): ("", ""),
    # ("Link", "urban decentral water tanks charger", "urban decentral water tanks"): ("", ""),
    # ("Link", "urban decentral water tanks discharger", "urban decentral heat"): ("", ""),
    # ("Link", "urban decentral water tanks discharger", "urban decentral water tanks"): ("", ""),
    ("Link", "waste CHP", "AC"): ("Secondary Solids In", "Secondary AC Out"),
    # ("Link", "waste CHP", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP", "non-sequestered HVC losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "waste CHP", "urban central heat"): (
        "Secondary Solids In",
        "Secondary Heat Out",
    ),
    ("Link", "waste CHP CC", "AC"): ("Secondary Solids In", "Secondary AC Out"),
    # ("Link", "waste CHP CC", "non-sequestered HVC"): ("", ""),
    ("Link", "waste CHP CC", "non-sequestered HVC losses"): (
        "Secondary Solids In",
        "Transformation Losses",
    ),
    ("Link", "waste CHP CC", "urban central heat"): (
        "Secondary Solids In",
        "Secondary Heat Out",
    ),
}


class SankeyChart(ESMChart):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.location = self._df.index.unique(DM.LOCATION).item()
        self.year = self._df.index.unique(DM.YEAR).item()
        self._df = self._df.droplevel(DM.YEAR).droplevel(DM.LOCATION)
        self._df.columns = ["value"]
        # self._df = self._df.abs()

    def plot(self):
        # Concatenate the data with source and target columns
        df_agg = self.add_source_target_columns()
        df_agg = self.add_jumpers(df_agg)
        label_mapping = self.get_label_mapping(df_agg)
        df_agg = self.add_id_source_target_columns(df_agg, label_mapping)
        df_agg = self.add_customdata(df_agg, self.unit)
        df_agg = self.combine_duplicates(df_agg)
        df_agg = self.map_colors_from_bus_carrier(df_agg)

        labels = list(label_mapping)
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
            # hovermode="x",  # show all categories on mouse-over
        )
        # trace order always needs to be reversed to show correct order
        # of legend entries for relative bar charts
        self.fig.update_layout(legend={"traceorder": "reversed"})

        # export the metadata directly in the Layout property for JSON
        self.fig.update_layout(meta=[RUN_META_DATA])

    @staticmethod
    def add_jumpers(df_agg):
        for bus_carrier in ("AC", "H2", "Gas", "Liquids", "Solids", "Heat", "Biomass"):
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
        return df_agg

    @staticmethod
    def _source_target_mapping(idx: tuple) -> pd.Series:
        return pd.Series(LINK_MAPPING.get(idx, ""))

    def add_source_target_columns(self):
        _df = self._df.copy()
        _df["index"] = _df.index  # convert to tuples
        _df[["source", "target"]] = _df["index"].apply(self._source_target_mapping)
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

    @staticmethod
    def add_customdata(df_agg, unit):
        to_concat = []
        for _, data in df_agg.groupby(["source", "target"]):
            data = data.reset_index()
            carrier_values = [
                f"{c}: {prettify_number(v)} {unit}"
                for c, v in zip(data["carrier"], data["value"])
            ]
            data["link_customdata"] = "<br>".join(carrier_values)
            to_concat.append(data)

        return pd.concat(to_concat)

    @staticmethod
    def add_id_source_target_columns(df_agg, label_mapping):
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
