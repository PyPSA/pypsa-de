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
    ("Link", "BEV charger", "EV battery"): ("Secondary AC Out", "Transport"),
    ("Link", "BEV charger", "low voltage losses"): (
        "Secondary AC In",
        "Transformation Losses",
    ),
    # Wood gasification treated as primary energy to simplify energy flows
    ("Link", "BioSNG", "gas"): ("Solid Biomass", "Primary Gas"),
    ("Link", "BioSNG CC", "gas"): ("Solid Biomass", "Primary Gas"),
    # ("Link", "BioSNG", "solid biomass"): ("", ""),
    ("Link", "BioSNG", "solid biomass losses"): (
        "Solids Secondary In",
        "Transformation Losses",
    ),  # todo: need to subtract from Biomass primary?
    # ("Link", "BioSNG CC", "solid biomass"): ("", ""),
    ("Link", "BioSNG CC", "solid biomass losses"): (
        "Solids Secondary In",
        "Transformation Losses",
    ),
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
        "Transformation Losses",
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
    # todo: treat V2G as a storage technology same as battery
    ("Link", "V2G", "EV battery"): ("", ""),
    ("Link", "V2G", "EV battery losses"): ("Transport", "Losses"),
    ("Link", "V2G", "low voltage"): ("Transport", "Secondary AC In"),
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
    ("Link", "biogas to gas", "gas"): ("Biogas", "Primary gas"),
    ("Link", "biogas to gas CC", "gas"): ("Biogas", "Primary gas"),
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
    ("Generator", "coal", "coal"): ("Coal", "Primary Solids"),
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
    ("StorageUnit", "hydro supply", "AC"): ("Hydro Power", "Primary AC"),  # is primary
    ("Generator", "import H2", "H2"): ("Green Hydrogen", "Primary H2"),
    ("Generator", "import NH3", "NH3"): ("Green Liquids", "Primary Liquids"),
    ("Link", "import gas", "gas"): ("Green Gas", "Primary Gas"),
    ("Link", "import methanol", "methanol"): ("Green Liquids", "Primary Liquids"),
    ("Link", "import oil", "oil"): ("Green Liquids", "Primary Liquids"),
    ("Load", "industry electricity", "low voltage"): ("", ""),
    ("Link", "industry methanol", "industry methanol"): ("", ""),
    ("Link", "industry methanol", "methanol"): ("", ""),
    ("Load", "industry methanol", "industry methanol"): ("", ""),
    ("Link", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Link", "kerosene for aviation", "oil"): ("", ""),
    ("Load", "kerosene for aviation", "kerosene for aviation"): ("", ""),
    ("Load", "land transport EV", "EV battery"): ("", ""),
    ("Load", "land transport fuel cell", "H2"): ("", ""),
    ("Generator", "lignite", "lignite"): ("Coal", "Primary Solids"),
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
    ("Link", "nuclear", "AC"): ("Nuclear Power", "Primary AC"),
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
    ("Link", "urban central H2 CHP", "AC"): ("Secondary H2 In", "Secondary AC Out"),
    ("Link", "urban central H2 CHP", "H2 losses"): (
        "Secondary H2 In",
        "Transformation Losses",
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
        "Transformation Losses",
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
        "Transformation Losses",
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
            hovermode="x",  # all categories are shown by mouse-over
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
