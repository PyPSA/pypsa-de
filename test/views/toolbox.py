"""Collects Toolbox functions used to verify migrations."""


def get_filter_dictionaries(
    component: str = "capacities",
    language_parameter: str = "DE",
    int_ext: str = "int",
) -> dict:
    """
    Return filter and mapping dictionary for evaluations.

    Parameters
    ----------
    component
          which dictionary is choosen
    language_parameter
       "DE": German, "EN": English
    int_ext
        "int" mapping for internal use, f.ex. excel,
        "ext" mapping for external use, f.ex., plotly

    Returns
    -------
    :
        Dictionary in the chosen language
    """
    mappings_internal = {
        "capacities": {
            "urban central lignite CHP CC electric": ["Kohle", "Coal"],
            "urban central lignite CHP electric": ["Kohle", "Coal"],
            "lignite power plant (CC)": ["Kohle", "Coal"],
            "lignite power plant": ["Kohle", "Coal"],
            "urban central coal CHP CC electric": ["Kohle", "Coal"],
            "urban central coal CHP electric": ["Kohle", "Coal"],
            "coal power plant": ["Kohle", "Coal"],
            "coal power plant (CC)": ["Kohle", "Coal"],
            "urban central gas CHP CC electric": ["Methan", "Methane"],
            "OCGT": ["Methan", "Methane"],
            "urban central gas CHP electric": ["Methan", "Methane"],
            "residential urban decentral micro gas CHP": ["Methan", "Methane"],
            "residential rural micro gas CHP": ["Methan", "Methane"],
            "services rural micro gas CHP": ["Methan", "Methane"],
            "services urban decentral micro gas CHP": ["Methan", "Methane"],
            "oil power plant": ["Mineralöl", "Oil"],
            "nuclear": ["Kernenergie", "Nuclear Power"],
            "onwind-1": ["Onshore", "Onshore"],
            "onwind-2": ["Onshore", "Onshore"],
            "onwind-3": ["Onshore", "Onshore"],
            "onwind-4": ["Onshore", "Onshore"],
            "offwind-ac": ["Offshore", "Offshore"],
            "offwind-dc": ["Offshore", "Offshore"],
            "ror": ["Laufwasserkraftwerke", "Run-of-River"],
            "PHS": ["Pumpspeicher", "Pumped Hydro Storage"],
            "hydro": ["Reservoirs", "Reservoir"],
            "urban central solid biomass CHP CC": ["Biomasse", "Biomass"],
            "urban central solid biomass CHP": ["Biomasse", "Biomass"],
            "solar-rooftop": ["PV-Dach", "PV-Rooftop"],
            "solar-utility": ["PV-Freifläche", "PV-Utility"],
            "H2 Fuel Cell": ["Wasserstoff", "Hydrogen"],
            "value of lost load": ["Schlupfvariable Last", "Power Disconnect"],
            "battery discharger": ["Batterie Speicher", "Battery Storage"],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "services urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "residential urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
        },
        "fuel_capacities": {
            "Fischer-Tropsch 1": ["Fischer-Tropsch 1", "Fischer-Tropsch 1"],
            "Fischer-Tropsch 2": ["Fischer-Tropsch 2", "Fischer-Tropsch 2"],
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["HT Elektrolyse", "Electrolysis HT"],
            "SMR": ["Dampfreformierung", "SMR"],
            "SMR CC": ["Dampfreformierung CC", "SMR"],
            "Sabatier": ["Methanisierung", "Methanation"],
            "biogas approximation": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            # "biogas to gas": ["Biomethanaufbereitung","Bio Methane Processing"],
            "helmeth": ["Helmeth", "Helmeth"],
        },
        "gas_store_capacities": {
            "gas": ["Methan Speicher", "Methane Store"],
            "H2 underground": [
                "Wasserstoff-Untertagespeicher",
                "Hydrogen Underground Storage",
            ],
            "H2 tube": ["Wasserstoff-Roehrenspeicher", "Hydrogen Tube Storage"],
        },
        "fuel_net_capacities": {
            "Fischer-Tropsch import link 1": [
                "Oel Importkapazitaet",
                "Oil import capacity",
            ],
            "Fischer-Tropsch import link 2": [
                "Oel Importkapazitaet",
                "Oil import capacity",
            ],
            "import capacity H2 foreign": [
                "Netzkapazitaet H2 Ausland",
                "Net Capacity H2 Foreign",
            ],
            "import capacity H2 domestic": [
                "Netzkapazitaet H2 Inland",
                "Net Capacity H2 Domestic",
            ],
            "import capacity gas foreign": [
                "Netzkapazitaet Gas Ausland",
                "Net Capacity Gas Foreign",
            ],
            "import capacity gas domestic": [
                "Netzkapazitaet Gas Inland",
                "Net Capacity Gas Domestic",
            ],
        },
        "elec_demand": {
            "domestic homes and trade": [
                "Haushalte & Gewerbe",
                "Households & Services",
            ],
            "electricity road freight": ["Güterverkehr", "Road Freight"],
            "industry": ["Industrie", "Industry"],
            "industry new electricity": [
                "Elektrif. Industrie",
                "Electrif. Industry",
            ],
            "grid losses": ["Netzverluste", "Grid Losses"],
            "electricity rail": ["Strom Bahn", "Electricity Rail"],
            "PHEV short": ["PKW Verkehr", "Passenger Transport PHEV"],
            "PHEV long": ["PKW Verkehr", "Passenger Transport PHEV"],
            "BEV to passenger demand": [
                "PKW Verkehr",
                "Passenger Transport BEV",
            ],
            # "V2G energy demand": ["V2G Speicher", "V2G Charger"],
            # "battery charger" : ["Batteriespeicher","Battery Storage"],
            "DAC": ["Direct Air Capture", "Direct Air Capture"],
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["Elektrolyse", "Electrolysis HT"],
            "residential rural ground heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services rural ground heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services rural resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential rural resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral air heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services urban decentral air heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services urban decentral resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "urban central air heat pump": ["Fernwärme", "District Heat"],
            "urban central resistive heater": ["Fernwärme", "District Heat"],
            "helmeth": ["P2G", "P2G"],
            # "PHS Stored Power" : ["Pumpspeicherverbrauch","Pumped Hydro Storage"],
            "foreign export": ["Export Ausland", "Export European"],
            "domestic export": ["Export Inland", "Export Domestic"],
        },
        "elec_prod": {
            "urban central lignite CHP CC electric": ["Kohle", "Coal"],
            "urban central lignite CHP electric": ["Kohle", "Coal"],
            "lignite power plant (CC)": ["Kohle", "Coal"],
            "lignite power plant": ["Kohle", "Coal"],
            "urban central coal CHP CC electric": ["Kohle", "Coal"],
            "urban central coal CHP electric": ["Kohle", "Coal"],
            "coal power plant": ["Kohle", "Coal"],
            "coal power plant (CC)": ["Kohle", "Coal"],
            "urban central gas CHP CC electric": ["Methan", "Methane"],
            "OCGT": ["Methan", "Methane"],
            "urban central gas CHP electric": ["Methan", "Methane"],
            "residential urban decentral micro gas CHP": ["Methan", "Methane"],
            "residential rural micro gas CHP": ["Methan", "Methane"],
            "services rural micro gas CHP": ["Methan", "Methane"],
            "services urban decentral micro gas CHP": ["Methan", "Methane"],
            "oil power plant": ["Mineralöl", "Oil"],
            "nuclear": ["Kernenergie", "Nuclear Power"],
            "onwind-1": ["Onshore", "Onshore"],
            "onwind-2": ["Onshore", "Onshore"],
            "onwind-3": ["Onshore", "Onshore"],
            "onwind-4": ["Onshore", "Onshore"],
            "offwind-ac": ["Offshore", "Offshore"],
            "offwind-dc": ["Offshore", "Offshore"],
            # "V2G energy back to network": ["V2G", "V2G Discharge"],
            "ror": ["Laufwasserkraftwerke", "Run-of-River"],
            "PHS Dispatched Power from Inflow": [
                "Pumpspeicher",
                "Pumped Hydro Storage",
            ],
            "hydro Dispatched Power": ["Reservoirs", "Reservoir"],
            "urban central solid biomass CHP CC": [
                "Biomasse CHP",
                "Biomass CHP",
            ],
            "urban central solid biomass CHP": ["Biomasse CHP", "Biomass CHP"],
            "solar-rooftop": ["PV-Dach", "PV-Rooftop"],
            "solar-utility": ["PV-Freifläche", "PV-Utility"],
            "H2 Fuel Cell": [
                "Wasserstoff Brennstoffzelle",
                "Hydrogen Fuel Cell",
            ],
            "value of lost load": ["Schlupfvariable Last", "Power Disconnect"],
            # "battery discharger": ["Sonstige","Miscellaneous"],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "services urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "residential urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
            "domestic import": ["Import Inland", "Import Domestic"],
            "foreign import": ["Import Ausland", "Import European"],
        },
        "gas_production_mapping_EU": {
            "Gas from Sabatier": ["Methanisierung", "Methanation"],
            # "gas store draw" : ["Methan Speichernutzung", "Methane Store Discharge"],
            # "biogas approximation": ["Biomethanaufbereitung","Bio Methane Processing"],
            "biogas to gas": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            "helmeth": ["Helmeth", "Helmeth"],
            "gas generator": ["Global Market*", "Global Market*"],
        },
        "gas_production_mapping": {
            "Gas from Sabatier": ["Methanisierung", "Methanation"],
            # "gas store draw" : ["Methan Speichernutzung", "Methane Store Discharge"],
            # "biogas approximation": ["Biomethanaufbereitung","Bio Methane Processing"],
            "biogas to gas": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            "helmeth": ["Helmeth", "Helmeth"],
            "gas foreign import": ["Import International", "Import European"],
            "gas generator": ["Global Market*", "Global Market*"],
        },
        "gas_demand_time_series": {
            "CNG long": ["Transport", "Transport"],
            "CNG short": ["Transport", "Transport"],
            "OCGT": ["Strom OCGT", "Electricity OCGT"],
            "SMR": ["Dampfreformierung", "SMR"],
            "SMR CC": ["Dampfreformierung", "SMR"],
            "gas Store": ["Speicher", "Storage In"],
            "gas domestic navigation": ["PKW", "Transport"],
            "gas feedstock": ["Industrie", "Industry"],
            "gas for industry": ["Industrie", "Industry"],
            "gas for industry CC": ["Industrie", "Industry"],
            "gas international navigation": ["PKW", "Transport"],
            "gas road freight": ["PKW", "Transport"],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential rural micro gas CHP": ["Sonstige", "Miscellaneous"],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural gas boiler": ["Dezentrale Wärme", "Decentral Heat"],
            "services rural micro gas CHP": ["Sonstige", "Miscellaneous"],
            "urban central gas CHP CC electric": [
                "Strom (KWK)",
                "Electricity CHP",
            ],
            "urban central gas CHP CC heat": ["Fernwärme", "District Heat"],
            "urban central gas CHP electric": [
                "Strom (KWK)",
                "Electricity CHP",
            ],
            "urban central gas CHP heat": ["Fernwärme", "District Heat"],
            "urban central gas boiler": ["Fernwärme", "District Heat"],
            "Net Export": ["Netto Export", "Net Export"],
        },
        "gas_demand_mapping_EU": {
            "OCGT": ["Strom (OCGT)", "Electricity (OCGT)"],
            "Gas for SMR CC": ["Dampfreformierung", "SMR"],
            "Gas for SMR": ["Dampfreformierung", "SMR"],
            "CNG short": ["Transport", "Transport"],
            "CNG long": ["Transport", "Transport"],
            "residential rural gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services rural gas boiler": ["Dezentrale Wärme", "Decentral Heat"],
            "urban central gas boiler": ["Fernwärme", "District Heat"],
            "residential rural micro gas CHP": ["Sonstige", "Miscellaneous"],
            "services rural micro gas CHP": ["Sonstige", "Miscellaneous"],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "urban central gas CHP electric": [
                "Strom (KWK)",
                "Electricity CHP",
            ],
            "urban central gas CHP CC electric": [
                "Strom (KWK)",
                "Electricity CHP",
            ],
            "urban central gas CHP heat": ["Fernwärme", "District Heat"],
            "urban central gas CHP CC heat": ["Fernwärme", "District Heat"],
            "services urban decentral gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral micro gas CHP": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services urban decentral micro gas CHP": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "gas for industry": ["Industrie", "Industry"],
            "gas for industry CC": ["Industrie", "Industry"],
            "gas feedstock": ["Industrie", "Industry"],
            "gas road freight": ["PKW", "Transport"],
            "gas domestic navigation": ["PKW", "Transport"],
            "gas international navigation": ["PKW", "Transport"],
        },
        "gas_demand_mapping": {
            "OCGT": ["Strom (OCGT)", "Electricity (OCGT)"],
            "Gas for SMR CC": ["Dampfreformierung", "SMR"],
            "Gas for SMR": ["Dampfreformierung", "SMR"],
            "CNG short": ["Transport", "Transport"],
            "CNG long": ["Transport", "Transport"],
            "residential rural gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services rural gas boiler": ["Dezentrale Wärme", "Decentral Heat"],
            "urban central gas boiler": ["Fernwärme", "District Heat"],
            "residential rural micro gas CHP": ["Sonstige", "Miscellaneous"],
            "services rural micro gas CHP": ["Sonstige", "Miscellaneous"],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "urban central gas CHP electric": [
                "Strom (KWK)",
                "Electricity CHP",
            ],
            "urban central gas CHP CC electric": [
                "Strom (KWK)",
                "Electricity CHP",
            ],
            "urban central gas CHP heat": ["Fernwärme", "District Heat"],
            "urban central gas CHP CC heat": ["Fernwärme", "District Heat"],
            "residential urban decentral gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services urban decentral gas boiler": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral micro gas CHP": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services urban decentral micro gas CHP": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "gas for industry": ["Industrie", "Industry"],
            "gas for industry CC": ["Industrie", "Industry"],
            "gas feedstock": ["Industrie", "Industry"],
            "gas road freight": ["PKW", "Transport"],
            "gas domestic navigation": ["PKW", "Transport"],
            "gas international navigation": ["PKW", "Transport"],
            # "gas store fill": ["Speicher","Storage"],
            "gas foreign export": ["Ausland Export", "Export European"],
            "gas domestic export": ["Inland Export", "Export Domestic"],
        },
        "gas_prod_time_series": {
            "Sabatier": ["Methanisierung", "Methanation"],
            "biogas to gas": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            # "gas": ["Netto Import","Net Import"],
            "gas Store": ["Speicher", "Storage Out"],
            "Net Import": ["Netto Import", "Net Import"],
        },
        "h2_production_mapping_EU": {
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["Elektrolyse", "Electrolysis"],
            # "H2 store draw" : ["Speichernutzung", "Store Discharge"],
            "H2 from SMR": ["Dampfreformierung", "SMR"],
            "H2 from SMR CC": ["Dampfreformierung CC", "SMR"],
            "H2 Import RU": ["Nicht-EU Import", "Non-EU Import"],
            "H2 Import NAF": ["Nicht-EU Import", "Non-EU Import"],
        },
        "h2_production_mapping": {
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["Elektrolyse", "Electrolysis"],
            # "H2 store draw" : ["Speichernutzung", "Store Discharge"],
            "H2 from SMR": ["Dampfreformierung", "SMR"],
            "H2 from SMR CC": ["Dampfreformierung CC", "SMR"],
            "H2 foreign import": ["Import International", "Import European"],
            "H2 domestic import": ["Import National", "Import Domestic"],
            "H2 retro foreign import": [
                "Import International",
                "Import European",
            ],
            "H2 retro domestic import": ["Import National", "Import Domestic"],
            "H2 Import RU": ["Import International", "Import Global"],
            "H2 Import NAF": ["Import International", "Import Global"],
        },
        "h2_demand_mapping_EU": {
            "FCEV long": ["Transport", "Transport"],
            "FCEV short": ["Transport", "Transport"],
            "Fischer-Tropsch 1": ["Synth. Kraftst.", "Synth. Fuels"],
            "Fischer-Tropsch 2": ["Synth. Kraftst.", "Synth. Fuels"],
            "services rural H2-powered PEMFC": ["Sonstiges", "Miscellaneous"],
            "residential rural H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "services urban decentral H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "residential urban decentral H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "H2 for Sabatier": ["Synth. Kraftst.", "Synth. Fuels"],
            "H2 Fuel Cell": ["Brennstoffzelle", "Fuel Cell"],
            "H2 road freight": ["Transport", "Transport"],
            "H2 for industry": ["Industrie", "Industry"],
            "H2 for shipping": ["Transport", "Transport"],
            "H2 for rail": ["Transport", "Transport"],
            "H2 for aviation": ["Transport", "Transport"],
            # "H2 foreign export": ["Ausland Export", "Export European"],
            # "H2 retro foreign export": ["Ausland Export", "Export European"],
            # "H2 domestic export": ["Inland Export", "Export Domestic"],
            # "H2 retro domestic export": ["Inland Export", "Export Domestic"],
        },
        "h2_demand_mapping": {
            "FCEV long": ["Transport", "Transport"],
            "FCEV short": ["Transport", "Transport"],
            "Fischer-Tropsch 1": ["Synth. Kraftst.", "Synth. Fuels"],
            "Fischer-Tropsch 2": ["Synth. Kraftst.", "Synth. Fuels"],
            "services rural H2-powered PEMFC": ["Sonstiges", "Miscellaneous"],
            "residential rural H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "services urban decentral H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "residential urban decentral H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "H2 for Sabatier": ["Methanisierung", "Methanation"],
            "H2 Fuel Cell": ["Brennstoffzelle", "Fuel Cell"],
            "H2 road freight": ["Transport", "Transport"],
            "H2 for industry": ["Industrie", "Industry"],
            "H2 for shipping": ["Transport", "Transport"],
            "H2 for rail": ["Transport", "Transport"],
            "H2 for aviation": ["Transport", "Transport"],
            # "H2 store fill": ["Speicher","Storage"],
            "H2 foreign export": ["Ausland Export", "Export European"],
            "H2 retro foreign export": ["Ausland Export", "Export European"],
            "H2 domestic export": ["Inland Export", "Export Domestic"],
            "H2 retro domestic export": ["Inland Export", "Export Domestic"],
        },
        "co2_emissions_mapping": {
            "Fischer-Tropsch road freight": [
                "Fischer-Tropsch Straßenverkehr Fracht",
                "Fischer-Tropsch road freight",
            ],
            "Fischer-Tropsch rail": [
                "Fischer-Tropsch Schienenverkehr",
                "Fischer-Tropsch rail",
            ],
            "Fischer-Tropsch domestic navigation": [
                "Fischer-Tropsch Schiffverkehr inland",
                "Fischer-Tropsch domestic navigation",
            ],
            "Fischer-Tropsch domestic aviation": [
                "Fischer-Tropsch Flugverkehr inland",
                "Fischer-Tropsch domestic aviation",
            ],
            "Fischer-Tropsch industry": [
                "Fischer-Tropsch Industrie",
                "Fischer-Tropsch industry",
            ],
            "hard coal industry": ["Kohle Industrie", "hard coal industry"],
            "process emissions": ["Proezssemissionen", "process emissions"],
            "process emissions CC": [
                "Prozessemissionen CC",
                "process emissions CC",
            ],
            "gas domestic navigation": [
                "Gas Schiffverkehr inland",
                "gas domestic navigation",
            ],
            "gas road freight": [
                "Gas Straßenverkehr Fracht",
                "gas road freight",
            ],
            "OCGT": ["OCGT", "OCGT"],
            "SMR": ["Dampfreformierung", "SMR"],
            "coal power plant": ["Kohlekraftwerke", "Coal PP"],
            "lignite power plant": ["Kohlekraftwerke", "Coal PP"],
            "ICE short": ["ICE Kurzstrecke", "ICE short"],
            "ICE long": ["ICE Langstrecke", "ICE long"],
            "CNG short": ["CNG Kurzstrecke", "CNG short"],
            "CNG long": ["CNG Langstrecke", "CNG long"],
            "HEV short": ["HEV Kurzstrecke", "HEV short"],
            "HEV long": ["HEV Langstrecke", "HEV long"],
            "residential rural gas boiler": [
                "Wärme dezentral gas boiler",
                "residential rural gas boiler",
            ],
            "services rural gas boiler": [
                "Wärme dezentral gas boiler",
                "services rural gas boiler",
            ],
            "urban central gas boiler": [
                "Wärme zentral gas boiler",
                "urban central gas boiler",
            ],
            "urban central gas CHP electric": [
                "Gas BHKW",
                "urban central gas CHP electric",
            ],
            "urban central gas CHP heat": [
                "Gas BHKW",
                "urban central gas CHP heat",
            ],
            "urban central coal CHP electric": ["Kohle BHKW", "Coal CHP"],
            "urban central coal CHP heat": ["Kohle BHKW", "Coal CHP"],
            "urban central lignite CHP electric": ["Kohle BHKW", "Coal CHP"],
            "urban central lignite CHP heat": ["Kohle BHKW", "Coal CHP"],
            "biogas to gas": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            "gas for industry": ["Gas Industrie", "gas for industry"],
            "services rural oil boiler": [
                "Wärme dezentral Öl boiler",
                "services rural oil boiler",
            ],
            "residential rural oil boiler": [
                "Wärme dezentral Öl boiler",
                "residential rural oil boiler",
            ],
            "oil power plant": ["Ölkraftwerke", "Oil PP"],
            "PHEV short": ["PHEV Kurzstrecke", "PHEV short"],
            "PHEV long": ["PHEV Langstrecke", "PHEV long"],
            "residential rural micro gas CHP": [
                "Wärme dezentral micro gas CHP",
                "residential rural micro gas CHP",
            ],
            "services rural micro gas CHP": [
                "Wärme dezentral micro gas CHP",
                "services rural micro gas CHP",
            ],
            "DAC": ["Direct Air Capture", "Direct Air Capture"],
            "Fischer-Tropsch import 1": ["Bioölimporte", "Import Biofuels"],
            "SMR CC": ["Dampfreformierung CC", "SMR CC"],
            "coal power plant (CC)": ["Kohlekraftwerke CC", "Coal PP CC"],
            "lignite power plant (CC)": ["Kohlekraftwerke CC", "Coal PP CC"],
            "urban central gas CHP CC electric": ["Gas BHKW CC", "Gas CHP CC"],
            "urban central gas CHP CC heat": ["Gas BHKW CC", "Gas CHP CC"],
            "urban central coal CHP CC electric": [
                "Kohle BHKW CC",
                "Coal CHP CC",
            ],
            "urban central coal CHP CC heat": ["Kohle BHKW CC", "Coal CHP CC"],
            "urban central lignite CHP CC electric": [
                "Kohle BHKW CC",
                "Coal CHP CC",
            ],
            "urban central lignite CHP CC heat": [
                "Kohle BHKW CC",
                "Coal CHP CC",
            ],
            "gas for industry CC": ["Gas Industrie CC", "gas for industry CC"],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Wärme dezentral CH4-powered PEMFC with internal SMR",
                "residential rural CH4-powered PEMFC with internal SMR",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Wärme dezentral CH4-powered PEMFC with internal SMR",
                "services rural CH4-powered PEMFC with internal SMR",
            ],
            "co2 vent": ["co2 vent", "co2 vent"],
        },
        "h2_demand_time_series": {
            "H2 Store": ["Speicher", "Storage In"],
            "FCEV long": ["Transport", "Transport"],
            "FCEV short": ["Transport", "Transport"],
            "Fischer-Tropsch 1": ["Fischer-Tropsch", "Fischer-Tropsch"],
            "Fischer-Tropsch 2": ["Fischer-Tropsch", "Fischer-Tropsch"],
            "H2 Fuel Cell": ["Brennstoffzelle", "Fuel Cell"],
            "H2 for aviation": ["Transport", "Transport"],
            "H2 for industry": ["Industrie", "Industry"],
            "H2 for rail": ["Transport", "Transport"],
            "H2 for shipping": ["Transport", "Transport"],
            "H2 road freight": ["Transport", "Transport"],
            "Sabatier": ["Methanisierung", "Methanation"],
            "residential rural H2-powered PEMFC": [
                "Sonstiges",
                "Miscellaneous",
            ],
            "services rural H2-powered PEMFC": ["Sonstiges", "Miscellaneous"],
            "Net Export": ["Netto Export", "Net Export"],
        },
        "h2_prod_time_series": {
            "H2 Store": ["Speicher", "Storage Out"],
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "SMR": ["Dampfreformierung", "SMR"],
            "SMR CC": ["Dampfreformierung", "SMR"],
            "Net Import": ["Netto Import", "Net Import"],
        },
        "district_heat_prod_time_series": {
            "urban central air heat pump": ["", "Heat Pump"],
            "urban central coal CHP CC heat": ["", "Coal CHP"],
            "urban central coal CHP heat": ["", "Coal CHP"],
            "urban central gas CHP CC heat": ["", "Gas CHP"],
            "urban central gas CHP heat": ["", "Gas CHP"],
            "urban central gas boiler": ["", "Gas Boiler"],
            "urban central lignite CHP CC heat": ["", "Coal CHP"],
            "urban central lignite CHP heat": ["", "Coal CHP"],
            "urban central resistive heater": ["", "Resistive Heater"],
            "H2 Fuel Cell": ["", "Fuel Cell (Heat)"],
            "urban central solid biomass CHP": ["", "Biomass CHP"],
            "urban central solid biomass CHP CC": ["", "Biomass CHP"],
            "urban central solid biomass boiler": ["", "Solid Biomass Boiler"],
            "urban central solid biomass boiler CC": [
                "",
                "Solid Biomass Boiler",
            ],
            "Fischer-Tropsch 1": ["Fischer-Tropsch", "Fischer-Tropsch"],
            "Fischer-Tropsch 2": ["Fischer-Tropsch", "Fischer-Tropsch"],
            "urban central solar thermal collector": ["", "Solar Thermal"],
            "urban central water tanks discharger": ["Speicher", "Storage Out"],
        },
        "district_heat_dem_time_series": {
            "low-temperature heat for industry": [
                "Industrie (Niedrigtemp.)",
                "Industry",
            ],
            "hh and services": ["HH und Gewerbe", "HH and Services (Heat)"],
            "grid losses": ["Wärmenetzverlust", "Grid Losses"],
            "urban central water tanks charger": ["Speicher", "Storage In"],
            "DAC": ["DAC", "Direct Air Capture"],
        },
    }
    mappings_external = {
        "capacities": {
            "urban central lignite CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central lignite CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "lignite power plant (CC)": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "lignite power plant": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central coal CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central coal CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "coal power plant": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "coal power plant (CC)": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central gas CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "OCGT": ["Thermische Kraftwerke", "Thermal Powerplants"],
            "urban central gas CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential urban decentral micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential rural micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services rural micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services urban decentral micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "oil power plant": ["Thermische Kraftwerke", "Thermal Powerplants"],
            "nuclear": ["Kernenergie", "Nuclear Power"],
            "onwind-1": ["Windkraft", "Wind Power"],
            "onwind-2": ["Windkraft", "Wind Power"],
            "onwind-3": ["Windkraft", "Wind Power"],
            "onwind-4": ["Windkraft", "Wind Power"],
            "offwind-ac": ["Windkraft", "Wind Power"],
            "offwind-dc": ["Windkraft", "Wind Power"],
            "ror": ["Laufwasserkraftwerke", "Run-of-River"],
            "PHS": ["Pumpspeicher", "Pumped Hydro Storage"],
            "hydro": ["Reservoirs", "Reservoir"],
            "urban central solid biomass CHP CC": ["Biomasse", "Biomass"],
            "urban central solid biomass CHP": ["Biomasse", "Biomass"],
            "solar-rooftop": ["Photovoltaik", "Photovoltaics"],
            "solar-utility": ["Photovoltaik", "Photovoltaics"],
            "H2 Fuel Cell": ["Wasserstoff", "Hydrogen"],
            # "value of lost load": ["Schlupfvariable Last","Power Disconnect"],
            "battery discharger": ["Batterie Speicher", "Battery Storage"],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "services urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "residential urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
        },
        "fuel_capacities": {
            "Fischer-Tropsch 1": ["Fischer-Tropsch", "Fischer-Tropsch"],
            "Fischer-Tropsch 2": ["Fischer-Tropsch", "Fischer-Tropsch"],
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["HT Elektrolyse", "Electrolysis HT"],
            "SMR": ["Dampfreformierung", "SMR"],
            "SMR CC": ["Dampfreformierung", "SMR"],
            "Sabatier": ["Methanisierung", "Methanation"],
            "biogas approximation": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            # "biogas to gas": ["Biomethanaufbereitung","Bio Methane Processing"],
            "helmeth": ["Helmeth", "Helmeth"],
        },
        "fuel_net_capacities": {
            "Fischer-Tropsch import link 1": [
                "Oel Importkapazität",
                "Oil import capacity",
            ],
            "Fischer-Tropsch import link 2": [
                "Oel Importkapazität",
                "Oil import capacity",
            ],
            "import capacity H2 foreign": [
                "Netzkapazität H2 Ausland",
                "Net Capacity H2 Foreign",
            ],
            "import capacity H2 domestic": [
                "Netzkapazität H2 Inland",
                "Net Capacity H2 Domestic",
            ],
            "import capacity gas foreign": [
                "Netzkapazität Gas Ausland",
                "Net Capacity Gas Foreign",
            ],
            "import capacity gas domestic": [
                "Netzkapazität Gas Inland",
                "Net Capacity Gas Domestic",
            ],
        },
        "elec_demand": {
            "domestic homes and trade": [
                "Haushalte & Gewerbe",
                "Households & Services",
            ],
            "electricity road freight": ["Verkehr", "Transport"],
            "industry": ["Industrie", "Industry"],
            "industry new electricity": ["Elektrif. Industrie", "Industry"],
            # "grid losses": ["Netzverluste", "Grid Losses"],
            "electricity rail": ["Verkehr", "Transport"],
            "PHEV short": ["Verkehr", "Transport"],
            "PHEV long": ["Verkehr", "Transport"],
            "BEV to passenger demand": ["Verkehr", "Transport"],
            # "V2G energy demand": ["V2G Speicher", "V2G Charger"],
            # "battery charger" : ["Batteriespeicher","Battery Storage"],
            "DAC": ["Direct Air Capture", "Direct Air Capture"],
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["Elektrolyse", "Electrolysis"],
            "residential rural ground heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services rural ground heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services rural resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential rural resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral air heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "residential urban decentral resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services urban decentral air heat pump": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "services urban decentral resistive heater": [
                "Dezentrale Wärme",
                "Decentral Heat",
            ],
            "urban central air heat pump": ["Fernwärme", "District Heat"],
            "urban central resistive heater": ["Fernwärme", "District Heat"],
            "helmeth": ["P2G", "P2G"],
            # "PHS Stored Power" : ["Pumpspeicherverbrauch","Pumped Hydro Stored"],
            "foreign export": ["Export Ausland", "Export European"],
            "domestic export": ["Export Inland", "Export Domestic"],
        },
        "elec_prod": {
            "urban central lignite CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central lignite CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "lignite power plant (CC)": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "lignite power plant": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central coal CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central coal CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "coal power plant": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "coal power plant (CC)": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central gas CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "OCGT": ["Thermische Kraftwerke", "Thermal Powerplants"],
            "urban central gas CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential urban decentral micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential rural micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services rural micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services urban decentral micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "oil power plant": ["Thermische Kraftwerke", "Thermal Powerplants"],
            "nuclear": ["Kernenergie", "Nuclear Power"],
            "onwind-1": ["Windkraft", "Wind Power"],
            "onwind-2": ["Windkraft", "Wind Power"],
            "onwind-3": ["Windkraft", "Wind Power"],
            "onwind-4": ["Windkraft", "Wind Power"],
            "offwind-ac": ["Windkraft", "Wind Power"],
            "offwind-dc": ["Windkraft", "Wind Power"],
            # "V2G energy back to network": ["V2G", "V2G Discharge"],
            "PHS Dispatched Power from Inflow": [
                "Wasserkraft",
                "Inflow Hydro Storage",
            ],
            "ror": ["Wasserkraft ror", "Run-of-River"],
            "hydro Dispatched Power": ["Wasserkraft", "Inflow Hydro Storage"],
            "urban central solid biomass CHP CC": [
                "Biomasse CHP",
                "Biomass CHP",
            ],
            "urban central solid biomass CHP": ["Biomasse CHP", "Biomass CHP"],
            "solar-rooftop": ["Photovoltaik", "Photovoltaics"],
            "solar-utility": ["Photovoltaik", "Photovoltaics"],
            "H2 Fuel Cell": [
                "Wasserstoff Brennstoffzelle",
                "Hydrogen Fuel Cell",
            ],
            # "battery discharger": ["Batterie Entladung","Battery Discharge"],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Sonstige",
                "Miscellaneous",
            ],
            "services rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "services urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
            "residential rural H2-powered PEMFC": ["Sonstige", "Miscellaneous"],
            "residential urban decentral H2-powered PEMFC": [
                "Sonstige",
                "Miscellaneous",
            ],
            "domestic import": ["Import Inland", "Import Domestic"],
            "foreign import": ["Import Ausland", "Import European"],
        },
        "elec_times_series_dem": {
            # Base Load
            "Base Load": ["Basis Last", "Base Load"],
            "grid losses": ["Basis Last", "Base Load"],
            # Industry
            "industry": ["Industrie", "Industry"],  #
            "industry new electricity": ["Industrie", "Industry"],
            "industry current electricity": ["Industrie", "Industry"],
            "DAC": ["Industrie", "Industry"],
            # transport
            "electricity road freight": ["Verkehr", "Transport"],  #
            "electricity rail": ["Verkehr", "Transport"],  #
            "PHEV short": ["Verkehr", "Transport"],
            "PHEV long": ["Verkehr", "Transport"],
            "BEV charger": ["Verkehr", "Transport"],
            # storage in
            "battery charger": ["Einspeichern", "Storage In"],
            "PHS": ["Einspeichern", "Storage In"],
            # electrolysis
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["Elektrolyse", "Electrolysis"],
            # Heat
            "residential rural ground heat pump": ["Wärme", "Heat"],
            "services rural ground heat pump": ["Wärme", "Heat"],
            "services rural resistive heater": ["Wärme", "Heat"],
            "residential rural resistive heater": ["Wärme", "Heat"],
            "residential urban decentral air heat pump": ["Wärme", "Heat"],
            "residential urban decentral resistive heater": ["Wärme", "Heat"],
            "services urban decentral air heat pump": ["Wärme", "Heat"],
            "services urban decentral resistive heater": ["Wärme", "Heat"],
            "urban central air heat pump": ["Wärme", "Heat"],
            "urban central resistive heater": ["Wärme", "Heat"],
            # Net Export
            "Net Export": ["Net Export", "Net Export"],
        },
        "elec_times_series_prod": {
            # Thermal Powerplants
            "urban central lignite CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central lignite CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "lignite power plant (CC)": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "lignite power plant": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central coal CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central coal CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "coal power plant": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "coal power plant (CC)": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central gas CHP CC electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "OCGT": ["Thermische Kraftwerke", "Thermal Powerplants"],
            "urban central gas CHP electric": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential urban decentral micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential rural micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services rural micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services urban decentral micro gas CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "oil power plant": ["Thermische Kraftwerke", "Thermal Powerplants"],
            "services urban decentral CH4-powered PEMFC with internal SMR": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential urban decentral CH4-powered PEMFC with internal SMR": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services rural H2-powered PEMFC": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "services urban decentral H2-powered PEMFC": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential rural H2-powered PEMFC": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "residential urban decentral H2-powered PEMFC": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central solid biomass CHP CC": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "urban central solid biomass CHP": [
                "Thermische Kraftwerke",
                "Thermal Powerplants",
            ],
            "H2 Fuel Cell": ["Thermische Kraftwerke", "Thermal Powerplants"],
            # Nuclear
            "nuclear": ["Kernenergie", "Nuclear Power"],
            # Wind power
            "onwind-1": ["Windkraft", "Wind Power"],
            "onwind-2": ["Windkraft", "Wind Power"],
            "onwind-3": ["Windkraft", "Wind Power"],
            "onwind-4": ["Windkraft", "Wind Power"],
            "offwind-ac": ["Windkraft", "Wind Power"],
            "offwind-dc": ["Windkraft", "Wind Power"],
            # Storage out
            "V2G": ["Speicher Aus", "Storage Out"],
            "hydro": ["Speicher Aus", "Storage Out"],
            "battery discharger": ["Speicher Aus", "Storage Out"],
            # changed PHS from "Storage Out" to "Storage In" because
            # the new mapping combines production/demand time series
            # and the new evaluation explicitly renames positive and
            # negative amounts of the PHS time series.
            "PHS": ["Speicher Aus", "Storage In"],
            # run-of-river
            "ror": ["Wasserkraft ror", "Run-of-River"],
            # Photovoltaics
            "solar-rooftop": ["Photovoltaik", "Photovoltaics"],
            "solar-utility": ["Photovoltaik", "Photovoltaics"],
            # Net Import
            "Net Import": ["Net Import", "Net Import"],
        },
        "gas_production_mapping": {
            "Gas from Sabatier": ["Methanisierung", "Methanation"],
            # "gas store draw" : ["Methan Speichernutzung", "Methane Store Discharge"],
            # "biogas approximation": ["Biomethanaufbereitung","Bio Methane Processing"],
            "biogas to gas": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            "helmeth": ["Helmeth", "Helmeth"],
            "gas foreign import": ["Import International", "Import European"],
            "gas domestic import": ["Import National", "Import Domestic"],
            "gas generator": ["Globaler Markt", "Global Market*"],
        },
        "h2_production_mapping": {
            "H2 Electrolysis": ["Elektrolyse", "Electrolysis"],
            "H2 HT Electrolysis": ["Elektrolyse", "Electrolysis"],
            # "H2 store draw" : ["Wasserstoff Speichernutzung", "Hydrogen Store Discharge"],
            "H2 from SMR": ["Dampfreformierung", "SMR"],
            "H2 from SMR CC": ["Dampfreformierung CC", "SMR"],
            "H2 foreign import": ["Import International", "Import European"],
            "H2 domestic import": ["Import National", "Import Domestic"],
            "H2 retro foreign import": [
                "Import International",
                "Import European",
            ],
            "H2 retro domestic import": ["Import National", "Import Domestic"],
            "H2 Import RU": ["Import International", "Import Global"],
            "H2 Import NAF": ["Import International", "Import Global"],
        },
        "co2_emissions_mapping": {
            "Fischer-Tropsch road freight": ["Transport", "Transport"],
            "Fischer-Tropsch rail": ["Transport", "Transport"],
            "Fischer-Tropsch domestic navigation": ["Transport", "Transport"],
            "Fischer-Tropsch domestic aviation": ["Transport", "Transport"],
            "Fischer-Tropsch industry": ["Industrie", "Industry"],
            "hard coal industry": ["Industrie", "Industry"],
            "process emissions": ["Industrie", "Industry"],
            "process emissions CC": ["Industrie CC", "Industry CC"],
            "gas domestic navigation": ["Transport", "Transport"],
            "gas road freight": ["Transport", "Transport"],
            "OCGT": ["OCGT", "OCGT"],
            "SMR": ["Dampfreformierung", "SMR"],
            "coal power plant": ["Kohlekraftwerke", "Coal PP"],
            "lignite power plant": ["Kohlekraftwerke", "Coal PP"],
            "ICE short": ["Transport", "Transport"],
            "ICE long": ["Transport", "Transport"],
            "CNG short": ["Transport", "Transport"],
            "CNG long": ["Transport", "Transport"],
            "HEV short": ["Transport", "Transport"],
            "HEV long": ["Transport", "Transport"],
            "residential rural gas boiler": [
                "Wärme dezentral",
                "Decentral Heat",
            ],
            "services rural gas boiler": ["Wärme dezentral", "Decentral Heat"],
            "urban central gas boiler": ["Wärme zentral", "District Heat"],
            "urban central gas CHP electric": ["Gas BHKW", "Gas CHP"],
            "urban central gas CHP heat": ["Gas BHKW", "Gas CHP"],
            "urban central coal CHP electric": ["Kohle BHKW", "Coal CHP"],
            "urban central coal CHP heat": ["Kohle BHKW", "Coal CHP"],
            "urban central lignite CHP electric": ["Kohle BHKW", "Coal CHP"],
            "urban central lignite CHP heat": ["Kohle BHKW", "Coal CHP"],
            "biogas to gas": [
                "Biomethanaufbereitung",
                "Bio Methane Processing",
            ],
            "gas for industry": ["Industrie", "Industry"],
            "services rural oil boiler": ["Wärme dezentral", "Decentral Heat"],
            "residential rural oil boiler": [
                "Wärme dezentral",
                "Decentral Heat",
            ],
            "oil power plant": ["Ölkraftwerke", "Oil PP"],
            "PHEV short": ["Transport", "Transport"],
            "PHEV long": ["Transport", "Transport"],
            "residential rural micro gas CHP": [
                "Wärme dezentral",
                "Decentral Heat",
            ],
            "services rural micro gas CHP": [
                "Wärme dezentral",
                "Decentral Heat",
            ],
            "DAC": ["Direct Air Capture", "Direct Air Capture"],
            "Fischer-Tropsch import 1": ["Bioölimporte", "Import Biofuels"],
            "SMR CC": ["Dampfreformierung CC", "SMR CC"],
            "coal power plant (CC)": ["Kohlekraftwerke CC", "Coal PP CC"],
            "lignite power plant (CC)": ["Kohlekraftwerke CC", "Coal PP CC"],
            "urban central gas CHP CC electric": ["Gas BHKW CC", "Gas CHP CC"],
            "urban central gas CHP CC heat": ["Gas BHKW CC", "Gas CHP CC"],
            "urban central coal CHP CC electric": [
                "Kohle BHKW CC",
                "Coal CHP CC",
            ],
            "urban central coal CHP CC heat": ["Kohle BHKW CC", "Coal CHP CC"],
            "urban central lignite CHP CC electric": [
                "Kohle BHKW CC",
                "Coal CHP CC",
            ],
            "urban central lignite CHP CC heat": [
                "Kohle BHKW CC",
                "Coal CHP CC",
            ],
            "gas for industry CC": ["Industrie CC", "Industry CC"],
            "residential rural CH4-powered PEMFC with internal SMR": [
                "Wärme dezentral",
                "Decentral Heat",
            ],
            "services rural CH4-powered PEMFC with internal SMR": [
                "Wärme dezentral",
                "Decentral Heat",
            ],
            "co2 vent": ["co2 vent", "co2 vent"],
        },
    }
    lang_index = 0
    if language_parameter == "DE":
        lang_index = 0
    elif language_parameter == "EN":
        lang_index = 1

    if int_ext == "int":
        return_dict = {
            key: mappings_internal[component][key][lang_index]
            for key in mappings_internal[component].keys()
        }
    elif int_ext == "ext":
        if component in mappings_external:  # if there exists no external dict
            return_dict = {
                key: mappings_external[component][key][lang_index]
                for key in mappings_external[component].keys()
            }
        else:
            return_dict = {
                key: mappings_internal[component][key][lang_index]
                for key in mappings_internal[component].keys()
            }
    else:
        return_dict = {}

    return return_dict
