# -*- coding: utf-8 -*-
"""
Collect constant values and identifiers used for evaluations.

Values in this module do not need to be changed during runtime.
"""

import importlib
import re
from datetime import datetime as dt
from importlib.metadata import PackageNotFoundError
from subprocess import CalledProcessError, check_output

from frozendict import frozendict

# represents constants import time
NOW: str = dt.now().strftime("%Y%m%d%H%M%S")


class DataModel:
    """Metric data model constants."""

    LOCATION: str = "location"
    COMPONENT: str = "component"
    CARRIER: str = "carrier"
    BUS_CARRIER: str = "bus_carrier"
    METRIC: str = "metric"
    YEAR: str = "year"
    SNAPSHOTS: str = "snapshots"
    IDX_NAMES: list = [LOCATION, CARRIER, BUS_CARRIER]
    YEAR_IDX_NAMES: list = [YEAR, LOCATION, CARRIER, BUS_CARRIER]


class BusCarrier:
    """Container to collect all bus carrier names."""

    AC: str = "AC"
    DC: str = "DC"
    CH4: str = "gas"
    H2: str = "H2"
    TRANSPORT_P: str = "passenger transport"
    TRANSPORT_P_LONG: str = "passenger transport long"
    FT: str = "Fischer-Tropsch"
    FT_1: str = "Fischer-Tropsch 1"
    FT_2: str = "Fischer-Tropsch 2"
    HEAT_URBAN_CENTRAL: str = "urban central heat"
    HEAT_URBAN_DECENTRAL: str = "urban decentral heat"
    HEAT_RURAL: str = "rural heat"
    # ESM heat buses:
    # HEAT_URBAN_SERVICES: str = "services urban decentral heat"
    # HEAT_URBAN_RESIDENTIAL: str = "residential urban decentral heat"
    # HEAT_RURAL_SERVICES: str = "services rural heat"
    # HEAT_RURAL_RESIDENTIAL: str = "residential rural heat"
    LI_ION: str = "Li ion"
    BATTERY: str = "battery"
    HOME_BATTERY: str = "home battery"
    EV_BATTERY: str = "EV battery"
    SOLID_BIOMASS: str = "solid biomass"

    @classmethod
    def ac_stores(cls) -> list:
        return [
            cls.AC,
            cls.DC,
            cls.LI_ION,
            cls.BATTERY,
            cls.HOME_BATTERY,
            cls.EV_BATTERY,
        ]

    @classmethod
    def heat_buses(cls) -> list:
        return [cls.HEAT_URBAN_CENTRAL, cls.HEAT_URBAN_DECENTRAL, cls.HEAT_RURAL]


class Carrier:
    """Container to collect all carrier names."""

    chp_urban_central_lignite_cc: str = "urban central lignite CHP CC electric"
    chp_urban_central_lignite: str = "urban central lignite CHP electric"
    chp_urban_central_coal_cc: str = "urban central coal CHP CC electric"
    chp_urban_central_coal: str = "urban central coal CHP electric"
    chp_urban_central_ch4_cc: str = "urban central gas CHP CC electric"
    chp_urban_central_ch4: str = "urban central gas CHP electric"
    chp_urban_central_solid_biomass_cc: str = "urban central solid biomass CHP CC"
    chp_urban_central_solid_biomass: str = "urban central solid biomass CHP"

    chp_urban_decentral_micro_ch4: str = "residential urban decentral micro gas CHP"
    chp_urban_decentral_services_micro_ch4: str = (
        "services urban decentral micro gas CHP"
    )

    chp_rural_residential_micro_ch4: str = "residential rural micro gas CHP"
    chp_rural_services_micro_ch4: str = "services rural micro gas CHP"

    pemfc_urban_services_decentral_ch4_smr: str = (
        "services urban decentral CH4-powered PEMFC with internal SMR"
    )
    pemfc_rural_services_decentral_ch4_smr: str = (
        "residential rural CH4-powered PEMFC with internal SMR"
    )
    pemfc_rural_services_ch4_smr: str = (
        "services rural CH4-powered PEMFC with internal SMR"
    )
    pemfc_urban_residential_decentral_ch4_smr: str = (
        "residential urban decentral CH4-powered PEMFC with internal SMR"
    )
    pemfc_rural_services_h2_smr: str = "services rural H2-powered PEMFC"

    pemfc_urban_services_decentral_h2: str = "services urban decentral H2-powered PEMFC"
    pemfc_rural_residential_h2: str = "residential rural H2-powered PEMFC"
    pemfc_urban_residential_decentral_h2_smr: str = (
        "residential urban decentral H2-powered PEMFC"
    )

    pp_lignite_cc: str = "lignite power plant (CC)"
    pp_lignite: str = "lignite power plant"
    pp_coal: str = "coal power plant"
    pp_coal_cc: str = "coal power plant (CC)"
    pp_oil: str = "oil power plant"

    ocgt: str = "OCGT"

    nuclear: str = "nuclear"
    onwind_1: str = "onwind-1"
    onwind_2: str = "onwind-2"
    onwind_3: str = "onwind-3"
    onwind_4: str = "onwind-4"
    offwind_ac: str = "offwind-ac"
    offwind_dc: str = "offwind-dc"
    ror: str = "ror"
    phs: str = "PHS"
    hydro: str = "hydro"

    solar_rooftop: str = "solar-rooftop"
    solar_utility: str = "solar-utility"
    h2_fuel_cell: str = "H2 Fuel Cell"
    lost_load: str = "value of lost load"
    battery_discharger: str = "battery discharger"
    battery: str = "battery"

    ft_1: str = "Fischer-Tropsch 1"
    ft_2: str = "Fischer-Tropsch 2"
    h2_electrolysis: str = "H2 Electrolysis"
    h2_electrolysis_ht: str = "H2 HT Electrolysis"
    smr: str = "SMR"
    smr_cc: str = "SMR CC"
    sabatier: str = "Sabatier"
    biogas_approximation: str = "biogas approximation"
    helmeth: str = "helmeth"
    ch4: str = "gas"
    h2_cavern: str = "H2 cavern"
    h2_tube: str = "H2 tube"
    ft_import_link_1: str = "Fischer-Tropsch import link 1"
    ft_import_link_2: str = "Fischer-Tropsch import link 2"
    h2_import_capacity_foreign: str = "import capacity H2 foreign"
    h2_import_capacity_domestic: str = "import capacity H2 domestic"
    ch4_import_capacity_foreign: str = "import capacity gas foreign"
    ch4_import_capacity_domestic: str = "import capacity gas domestic"
    domestic_homes_and_trade: str = "domestic homes and trade"
    road_freight_ac: str = "electricity road freight"
    industry: str = "industry"
    industry_new_electricity: str = "industry new electricity"
    grid_losses: str = "urban central heat losses"
    electricity_rail: str = "electricity rail"
    phev_short: str = "PHEV short"
    phev_long: str = "PHEV long"

    v2g: str = "V2G"
    bev: str = "BEV"
    bev_charger: str = "BEV charger"
    bev_passenger_withdrawal: str = "BEV to passenger used"
    bev_charger_supply: str = "BEV charger out"
    bev_charger_draw: str = "BEV charger in"
    bev_charger_losses: str = "BEV charger losses"
    v2g_supply: str = "V2G energy back to network"
    v2g_withdrawal: str = "V2G energy draw"

    dac: str = "DAC"
    heat_pump_residential_rural_ground: str = "residential rural ground heat pump"
    heat_pump_ground_services_rural: str = "services rural ground heat pump"
    resistive_heater_rural_services: str = "services rural resistive heater"
    resistive_heater_rural_residential: str = "residential rural resistive heater"
    heat_pump_air_urban_residential_decentral: str = (
        "residential urban decentral air heat pump"
    )
    resistive_heater_urban_decentral_residential: str = (
        "residential urban decentral resistive heater"
    )
    heat_pump_air_services_urban_decentral: str = (
        "services urban decentral air heat pump"
    )
    resistive_heater_services_urban_decentral: str = (
        "services urban decentral resistive heater"
    )

    heat_pump_air_urban_central: str = "urban central air heat pump"
    resistive_heater_urban_central: str = "urban central resistive heater"
    export_foreign: str = "foreign export"
    export_domestic: str = "domestic export"
    phs_dispatched_power_inflow: str = "PHS Dispatched Power from Inflow"
    hydro_dispatched_power: str = "hydro Dispatched Power"
    import_domestic: str = "domestic import"
    import_foreign: str = "foreign import"
    ch4_from_sabatier: str = "Gas from Sabatier"
    biogas_to_ch4: str = "biogas to gas"
    gas_pipepline: str = "gas pipeline"
    gas_pipepline_new: str = "gas pipeline new"
    ch4_generator: str = "gas generator"
    ch4_import_foreign: str = "gas foreign import"
    cng_long: str = "CNG long"
    cng_short: str = "CNG short"
    ch4_store: str = "gas Store"
    ch4_navigation_domestic: str = "gas domestic navigation"
    ch4_feedstock: str = "gas feedstock"
    ch4_industry: str = "gas for industry"
    ch4_industry_cc: str = "gas for industry CC"
    ch4_navigation_international: str = "gas international navigation"
    road_freight_ch4: str = "gas road freight"
    ch4_boiler_residential_rural: str = "residential rural gas boiler"
    ch4_boiler_services_rural: str = "services rural gas boiler"

    chp_urban_central_ch4_heat_cc: str = "urban central gas CHP CC heat"
    chp_urban_central_ch4_heat: str = "urban central gas CHP heat"
    ch4_boiler_urban_central: str = "urban central gas boiler"
    export_net: str = "Net Export"
    ch4_for_smr_cc: str = "Gas for SMR CC"
    ch4_for_smr: str = "Gas for SMR"
    ch4_boiler_urban_decentral_services: str = "services urban decentral gas boiler"
    ch4_boiler_urban_decentral_residential: str = (
        "residential urban decentral gas boiler"
    )
    export_ch4_foreign: str = "gas foreign export"
    export_ch4_domestic: str = "gas domestic export"
    import_net: str = "Net Import"

    h2_from_smr: str = "H2 from SMR"
    h2_from_smr_cc: str = "H2 from SMR CC"
    h2_import_russia: str = "H2 Import RU"
    h2_import_naf: str = "H2 Import NAF"
    h2_import_foreign: str = "H2 foreign import"
    h2_import_domestic: str = "H2 domestic import"
    h2_import_foreign_retro: str = "H2 retro foreign import"
    h2_import_domestic_retro: str = "H2 retro domestic import"
    fcev_long: str = "FCEV long"
    fcev_short: str = "FCEV short"
    h2_sabatier: str = "H2 for Sabatier"
    h2_pipeline: str = "H2 pipeline"
    h2_pipeline_retro: str = "H2 pipeline retrofitted"
    h2_pipeline_kernnetz: str = "H2 pipeline (Kernnetz)"
    road_freight_h2: str = "H2 road freight"
    h2_industry: str = "H2 for industry"
    h2_shipping: str = "H2 for shipping"
    h2_rail: str = "H2 for rail"
    h2_aviation: str = "H2 for aviation"
    h2_export_foreign: str = "H2 foreign export"
    h2_export_foreign_retro: str = "H2 retro foreign export"
    h2_export_domestic: str = "H2 domestic export"
    h2_export_domestic_retro: str = "H2 retro domestic export"
    road_freight_ft: str = "Fischer-Tropsch road freight"
    ft_rail: str = "Fischer-Tropsch rail"
    ft_domestic_navigation: str = "Fischer-Tropsch domestic navigation"
    ft_domestic_aviation: str = "Fischer-Tropsch domestic aviation"
    ft_industry: str = "Fischer-Tropsch industry"
    hard_coal_industry: str = "hard coal industry"
    process_emissions: str = "process emissions"
    process_emissions_cc: str = "process emissions CC"
    ice_short: str = "ICE short"
    ice_long: str = "ICE long"
    hev_short: str = "HEV short"
    hev_long: str = "HEV long"

    chp_urban_central_coal_heat: str = "urban central coal CHP heat"
    chp_urban_central_lignite_heat: str = "urban central lignite CHP heat"
    chp_urban_central_coal_heat_cc: str = "urban central coal CHP CC heat"
    chp_urban_central_lignite_heat_cc: str = "urban central lignite CHP CC heat"
    oil: str = "oil"
    oil_boiler_rural_services: str = "services rural oil boiler"
    oil_boiler_rural_residential: str = "residential rural oil boiler"
    oil_boiler_urban_residential: str = "residential urban decentral oil boiler"
    ft_import_1: str = "Fischer-Tropsch import 1"
    co2_vent: str = "co2 vent"
    h2_store: str = "H2 Store"
    solid_biomass_boiler_urban_central: str = "urban central solid biomass boiler"
    solid_biomass_boiler_urban_central_cc: str = "urban central solid biomass boiler CC"
    solar_thermal_collector_urban_central: str = "urban central solar thermal collector"
    water_tanks_discharger_urban_central: str = "urban central water tanks discharger"
    water_tanks_charger_urban_central: str = "urban central water tanks charger"
    low_temperature_heat_for_industry: str = "low-temperature heat for industry"
    hh_and_services: str = "hh and services"
    value_lost_load: str = "value of lost load"

    # derivative metric names
    bev_demand: str = "BEV to passenger demand"
    bev_losses: str = "BEV to passenger losses"
    v2g_demand: str = "V2G energy demand"
    v2g_losses: str = "V2G energy total losses"


class Group:
    """Container to collect all carrier nice names."""

    phs_inflow: str = "Inflow Hydro Storage"
    base_load: str = "Base Load"
    battery_storage: str = "Battery Storage"
    biomass: str = "Biomass"
    ch4_bio_processing: str = "Bio Methane Processing"
    chp_biomass: str = "Biomass CHP"
    cng_long: str = "CNG long"
    cng_short: str = "CNG short"
    coal: str = "Coal"
    chp_coal: str = "Coal CHP"
    chp_coal_cc: str = "Coal CHP CC"
    pp_coal: str = "Coal PP"
    pp_coal_cc: str = "Coal PP CC"
    heat_decentral: str = "Decentral Heat"
    dac: str = "Direct Air Capture"
    heat_district: str = "District Heat"
    heat: str = "Heat"
    electrictiy: str = "Electricity"
    rail_elecricity: str = "Electricity Rail"
    ocgt_electricity: str = "Electricity OCGT"
    chp_electricity: str = "Electricity CHP"
    industry_electrification: str = "Electrif. Industry"
    electrolysis: str = "Electrolysis"
    electrolysis_ht: str = "Electrolysis HT"
    export_domestic: str = "Export Domestic"
    export_foreign: str = "Export Foreign"
    ft: str = "Fischer-Tropsch"
    ft_1: str = "Fischer-Tropsch 1"
    ft_2: str = "Fischer-Tropsch 2"
    ft_domestic_aviation: str = "Fischer-Tropsch domestic aviation"
    ft_domestic_navigation: str = "Fischer-Tropsch domestic navigation"
    ft_industry: str = "Fischer-Tropsch industry"
    ft_rail: str = "Fischer-Tropsch rail"
    ft_road_freight: str = "Fischer-Tropsch road freight"
    fuel_cell: str = "Fuel Cell"
    fuel_cell_heat: str = "Fuel Cell (Heat)"
    ch4_boiler: str = "Gas Boiler"
    chp_ch4: str = "Gas CHP"
    chp_ch4_cc: str = "Gas CHP CC"
    global_market: str = "Global Market*"
    grid_losses: str = "Grid Losses"
    hev_long: str = "HEV long"
    hev_short: str = "HEV short"
    hh_and_services_heat: str = "HH and Services (Heat)"
    # HT_Electrolysis: str = "HT Electrolysis"
    heat_pump: str = "Heat Pump"
    helmeth: str = "Helmeth"
    hh_and_services: str = "Households & Services"
    h2: str = "Hydrogen"
    h2_fuel_cell: str = "Hydrogen Fuel Cell"
    h2_tube_storage: str = "Hydrogen Tube Storage"
    h2_underground_storage: str = "Hydrogen Underground Storage"
    ice_long: str = "ICE long"
    ice_short: str = "ICE short"
    import_biofuels: str = "Import Biofuels"
    import_domestic: str = "Import Domestic"
    import_foreign: str = "Import Foreign"
    import_global: str = "Import Global"
    industry: str = "Industry"
    industry_cc: str = "Industry CC"
    methanation: str = "Methanation"
    ch4: str = "Methane"
    ch4_store: str = "Methane Store"
    misc: str = "Miscellaneous"
    ch4_capacity_domestic_net: str = "Net Capacity Gas Domestic"
    ch4_capacity_foreign_net: str = "Net Capacity Gas Foreign"
    h2_capacity_domestic_net: str = "Net Capacity H2 Domestic"
    h2_capacity_foreign_net: str = "Net Capacity H2 Foreign"
    export_net: str = "Net Export"
    import_net: str = "Net Import"
    import_non_eu: str = "Non-EU Import"
    nuclear_power: str = "Nuclear Power"
    ocgt: str = "OCGT"
    wind_offshore: str = "Offshore"
    wind: str = "Wind Power"
    oil: str = "Oil"
    oil_boiler: str = "Oil Boiler"
    pp_oil: str = "Oil PP"
    pp_thermal: str = "Thermal Powerplants"
    import_capacity_oil: str = "Oil import capacity"
    wind_onshore: str = "Onshore"
    p2g: str = "P2G"
    phev_long: str = "PHEV long"
    phev_short: str = "PHEV short"
    pv: str = "Photovoltaics"
    pv_rooftop: str = "PV-Rooftop"
    pv_utility: str = "PV-Utility"
    bev_passenger_transport: str = "Passenger Transport BEV"
    phev: str = "Passenger Transport PHEV"
    power_disconnect: str = "Power Disconnect"
    phs: str = "Pumped Hydro Storage"
    reservoir: str = "Reservoir"
    resistive_heater: str = "Resistive Heater"
    road_freight: str = "Road Freight"
    ror: str = "Run-of-River"
    smr: str = "SMR"
    smr_cc: str = "SMR CC"
    solar_thermal: str = "Solar Thermal"
    solid_biomass_boiler: str = "Solid Biomass Boiler"
    storage_in: str = "Storage In"
    storage_out: str = "Storage Out"
    synth_fuels: str = "Synth. Fuels"
    transport: str = "Transport"
    co2_vent: str = "co2 vent"
    ch4_domestic_navigation: str = "gas domestic navigation"
    ch4_industry: str = "gas for industry"
    ch4_industry_cc: str = "gas for industry CC"
    ch4_road_freight: str = "gas road freight"
    coal_industry: str = "hard coal industry"
    process_emissions: str = "process emissions"
    process_emissions_cc: str = "process emissions CC"
    ch4_rural_residential_pemfc_smr: str = (
        "residential rural CH4-powered PEMFC with internal SMR"
    )
    ch4_boiler_rural_residential: str = "residential rural gas boiler"
    ch4_rural_residential_chp: str = "residential rural micro gas CHP"
    ch4_rural_residential_oil_boiler: str = "residential rural oil boiler"
    ch4_rural_services_pemfc_smr: str = (
        "services rural CH4-powered PEMFC with internal SMR"
    )
    ch4_rural_services_boiler: str = "services rural gas boiler"
    ch4_rural_services_chp: str = "services rural micro gas CHP"
    oil_rural_services_boiler: str = "services rural oil boiler"
    ch4_urban_central_chp: str = "urban central gas CHP electric"
    ch4_urban_central_chp_heat: str = "urban central gas CHP heat"
    ch4_urban_central_boiler: str = "urban central gas boiler"
    soc: str = "State of Charge"
    soc_max: str = "Max State of Charge"
    turbine_cum: str = "Accumulated Turbining"
    pumping_cum: str = "Accumulated Pumping"
    spill_cum: str = "Accumulated Outflow Spill"
    inflow_cum: str = "Accumulated Natural Inflow"


class Regex:
    """A collection of regular expression patterns."""

    # ends with 4 digits
    year: re.Pattern = re.compile(r"\d{4}$")

    # matches: startswith 2 capital letters, followed by up to 3 digits,
    # 1 space, and any number of digits for optional subnets.
    region: re.Pattern = re.compile(r"^[A-Z]{2}[\d,A-G]{0,3}\s*\d*")

    # matches: startswith 2 capital letters, followed by up to 3 digits,
    # groups: only the first 2 letters that are the country code
    country: re.Pattern = re.compile(r"^([A-Z]{2})[\d,A-G]{0,3}\s*")

    # match anything inside parenthesis.
    unit: re.Pattern = re.compile(r"\([^()]*\)")


class Mapping:
    """Collects dictionaries used to rename technologies and carrier."""

    carrier = frozendict(
        {
            Carrier.bev_charger: Group.transport,
            Carrier.bev_charger_draw: Group.transport,
            Carrier.bev_charger_supply: Group.transport,
            Carrier.bev_demand: {
                "internal": "Passenger Transport BEV",
                "external": Group.transport,
            },
            Carrier.bev_losses: Group.transport,
            Carrier.bev_passenger_withdrawal: Group.transport,
            "Base Load": Group.base_load,
            Carrier.cng_long: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.cng_short: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.dac: {
                "internal": {
                    "default": Group.dac,
                    "electricity": Group.industry,
                },
                "external": {
                    "default": Group.dac,
                    "electricity": Group.industry,
                },
            },
            Carrier.fcev_long: Group.transport,
            Carrier.fcev_short: Group.transport,
            Carrier.ft_1: {
                "internal": {
                    "default": Group.ft_1,
                    "h2_demand": Group.synth_fuels,
                    "h2_balance": Group.synth_fuels,
                    "h2_demand_ts": Group.ft,
                    "h2_demand_eu": Group.synth_fuels,
                    "district_heat": Group.ft,
                },
                "external": {
                    "default": Group.ft,
                    "h2_demand": Group.synth_fuels,
                    "h2_balance": Group.synth_fuels,
                    "h2_demand_ts": Group.ft,
                    "h2_demand_eu": Group.synth_fuels,
                },
            },
            Carrier.ft_2: {
                "internal": {
                    "default": Group.ft_2,
                    "h2_demand": Group.synth_fuels,
                    "h2_balance": Group.synth_fuels,
                    "h2_demand_ts": Group.ft,
                    "h2_demand_eu": Group.synth_fuels,
                    "district_heat": Group.ft,
                },
                "external": {
                    "default": Group.ft,
                    "h2_demand": Group.synth_fuels,
                    "h2_balance": Group.synth_fuels,
                    "h2_demand_ts": Group.ft,
                    "h2_demand_eu": Group.synth_fuels,
                },
            },
            Carrier.ft_domestic_aviation: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.ft_domestic_navigation: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.ft_import_1: Group.import_biofuels,
            Carrier.ft_import_link_1: Group.import_capacity_oil,
            Carrier.ft_import_link_2: Group.import_capacity_oil,
            Carrier.ft_industry: {
                "internal": {"default": Group.industry, "co2": None},
                "external": Group.industry,
            },
            "Fischer-Tropsch international aviation": Group.transport,
            "Fischer-Tropsch international navigation": Group.transport,
            Carrier.ft_rail: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            "Fischer-Tropsch road freight": {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            "Gas for SMR": Group.smr,
            "Gas for SMR CC": Group.smr,
            "Gas from Sabatier": Group.methanation,
            Carrier.h2_electrolysis: Group.electrolysis,
            Carrier.h2_fuel_cell: {
                "internal": {
                    "default": Group.h2,
                    "capacity": Group.fuel_cell,
                    "h2_demand": Group.fuel_cell,
                    "h2_demand_ts": Group.fuel_cell,
                    "h2_demand_eu": Group.fuel_cell,
                    "e_production": Group.h2_fuel_cell,
                    "electricity": Group.pp_thermal,
                    "district_heat": Group.fuel_cell_heat,
                },
                "external": {
                    "default": Group.h2,
                    "capacity": Group.fuel_cell,
                    "h2_demand": Group.fuel_cell,
                    "h2_demand_ts": Group.fuel_cell,
                    "h2_demand_eu": Group.fuel_cell,
                    "e_production": Group.h2_fuel_cell,
                    "electricity": Group.pp_thermal,
                    "district_heat": Group.fuel_cell_heat,
                },
            },
            Carrier.h2_electrolysis_ht: {
                "internal": {
                    "default": Group.electrolysis,
                    "capacity": Group.electrolysis_ht,
                    "e_demand": Group.electrolysis_ht,
                },
                "external": {
                    "default": Group.electrolysis,
                    "capacity": Group.electrolysis_ht,
                },
            },
            Carrier.h2_store: {
                "internal": {"default": Group.storage_in, "h2_prod": Group.storage_out},
                "external": {"default": Group.storage_in, "h2_prod": Group.storage_out},
            },
            Carrier.h2_cavern: "Hydrogen Underground Storage",
            "H2 domestic aviation": Group.transport,
            "H2 domestic navigation": Group.transport,
            Carrier.h2_sabatier: {
                "internal": {
                    "default": Group.methanation,
                    "h2_demand_eu": Group.synth_fuels,
                },
                "external": {
                    "default": Group.methanation,
                    "h2_demand_eu": Group.synth_fuels,
                },
            },
            Carrier.h2_aviation: Group.transport,
            Carrier.h2_industry: Group.industry,
            Carrier.h2_rail: Group.transport,
            Carrier.h2_shipping: Group.transport,
            "H2 international aviation": Group.transport,
            "H2 international navigation": Group.transport,
            "H2 rail": "H2 for rail",
            "H2 road freight": Group.transport,
            Carrier.h2_tube: "Hydrogen Tube Storage",
            Carrier.hev_long: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.hev_short: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.ice_long: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.ice_short: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.export_net: Group.export_net,
            Carrier.import_net: Group.import_net,
            Carrier.ocgt: {
                "internal": {
                    "default": Group.ch4,
                    "co2": Group.ocgt,
                    "capacity": Group.ocgt_electricity,
                },
                "external": {
                    "default": Group.pp_thermal,
                    "co2": Group.ocgt,
                    "capacity": Group.ocgt_electricity,
                },
            },
            Carrier.phev_long: {
                "internal": {"default": "Passenger Transport PHEV", "co2": None},
                "external": Group.transport,
            },
            Carrier.phev_short: {
                "internal": {"default": "Passenger Transport PHEV", "co2": None},
                "external": Group.transport,
            },
            Carrier.phs: {
                "internal": {
                    "default": Group.phs,
                    "electricity": Group.storage_in,
                },
                "external": {
                    "default": Group.phs,
                    "electricity": Group.storage_in,
                },
            },
            Carrier.phs_dispatched_power_inflow: {
                "internal": {
                    "default": Group.phs_inflow,
                    "e_production": Group.phs,
                },
                "external": Group.phs_inflow,
            },
            Carrier.smr: Group.smr,
            Carrier.smr_cc: {
                "internal": {
                    "default": Group.smr_cc,
                    "capacity": Group.smr,
                    "h2_prod": Group.smr,
                },
                "external": {
                    "default": Group.smr_cc,
                    "capacity": Group.smr,
                    "h2_prod": Group.smr,
                    "h2_balance": Group.smr,
                },
            },
            Carrier.sabatier: Group.methanation,
            Carrier.v2g: Group.storage_out,
            Carrier.v2g_supply: Group.storage_in,
            Carrier.v2g_demand: Group.storage_in,
            Carrier.v2g_withdrawal: Group.storage_out,
            Carrier.v2g_losses: Group.storage_out,
            "battery charger": Group.storage_in,
            Carrier.battery_discharger: {
                "internal": {
                    "default": Group.battery_storage,
                    "electricity": Group.storage_in,
                },
                "external": {
                    "default": Group.battery_storage,
                    "electricity": Group.storage_in,
                },
            },
            Carrier.battery: Group.battery_storage,
            Carrier.biogas_approximation: Group.ch4_bio_processing,
            Carrier.biogas_to_ch4: Group.ch4_bio_processing,
            Carrier.co2_vent: Group.co2_vent,
            Carrier.pp_coal: {
                "internal": {"default": Group.coal, "co2": Group.pp_coal},
                "external": {"default": Group.pp_thermal, "co2": Group.pp_coal},
            },
            Carrier.pp_coal_cc: {
                "internal": {"default": Group.coal, "co2": Group.pp_coal_cc},
                "external": {"default": Group.pp_thermal, "co2": Group.pp_coal_cc},
            },
            Carrier.export_domestic: Group.export_domestic,
            Carrier.domestic_homes_and_trade: Group.hh_and_services,
            Carrier.import_domestic: Group.import_domestic,
            "electricity": Group.base_load,
            Carrier.electricity_rail: {
                "internal": Group.rail_elecricity,
                "external": Group.transport,
            },
            "electricity road freight": {
                "internal": Group.road_freight,
                "external": Group.transport,
            },
            Carrier.export_foreign: Group.export_foreign,
            Carrier.import_foreign: Group.import_foreign,
            Carrier.ch4_store: {
                "internal": {
                    "default": Group.storage_in,
                    "gas_prod": Group.storage_out,
                },
                "external": {
                    "default": Group.storage_in,
                    "gas_prod": Group.storage_out,
                },
            },
            "gas domestic import": "Import Domestic",
            Carrier.ch4_navigation_domestic: {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.ch4_feedstock: Group.industry,
            Carrier.ch4_industry: {
                "internal": {"default": Group.industry, "co2": None},
                "external": Group.industry,
            },
            Carrier.ch4_industry_cc: {
                "internal": {
                    "default": Group.industry_cc,
                    "capacity": Group.industry,
                    "co2": None,
                },
                "external": {"default": Group.industry, "co2": Group.industry_cc},
            },
            "gas foreign import": "Import European",
            "gas generator": "Global Market*",
            "gas international navigation": Group.transport,
            "gas road freight": {
                "internal": {"default": Group.transport, "co2": None},
                "external": Group.transport,
            },
            Carrier.grid_losses: {
                "internal": {
                    "default": Group.base_load,
                    "district_heat": Group.grid_losses,
                    "e_demand": Group.grid_losses,
                },
                "external": {
                    "default": Group.base_load,
                    "district_heat": Group.grid_losses,
                },
            },
            Carrier.hard_coal_industry: {
                "internal": {"default": Group.industry, "co2": None},
                "external": Group.industry,
            },
            Carrier.helmeth: {
                "internal": {"e_demand": Group.p2g, "default": Group.helmeth},
                "external": {"e_demand": Group.p2g, "default": Group.helmeth},
            },
            Carrier.hh_and_services: Group.hh_and_services_heat,
            Carrier.hydro: {
                "internal": {
                    "default": Group.reservoir,
                    "electricity": Group.storage_out,
                },
                "external": {
                    "default": Group.reservoir,
                    "electricity": Group.storage_out,
                },
            },
            Carrier.hydro_dispatched_power: {
                "internal": {
                    "default": Group.phs_inflow,
                    "e_production": Group.reservoir,
                },
                "external": Group.phs_inflow,
            },
            Carrier.h2_import_capacity_domestic: "Net Capacity H2 Domestic",
            Carrier.h2_import_capacity_foreign: "Net Capacity H2 Foreign",
            Carrier.ch4_import_capacity_domestic: "Net Capacity Gas Domestic",
            Carrier.ch4_import_capacity_foreign: "Net Capacity Gas Foreign",
            Carrier.industry: Group.industry,
            "industry current electricity": Group.industry,
            Carrier.industry_new_electricity: {
                "internal": "Electrif. Industry",
                "external": Group.industry,
            },
            Carrier.pp_lignite: {
                "internal": {"default": Group.coal, "co2": Group.pp_coal},
                "external": {"default": Group.pp_thermal, "co2": Group.pp_coal},
            },
            Carrier.pp_lignite_cc: {
                "internal": {"default": Group.coal, "co2": Group.pp_coal_cc},
                "external": {"default": Group.pp_thermal, "co2": Group.pp_coal_cc},
            },
            Carrier.low_temperature_heat_for_industry: Group.industry,
            Carrier.nuclear: Group.nuclear_power,
            Carrier.offwind_ac: {
                "internal": Group.wind_offshore,
                "external": Group.wind,
            },
            Carrier.offwind_dc: {
                "internal": Group.wind_offshore,
                "external": Group.wind,
            },
            Carrier.pp_oil: {
                "internal": {"default": Group.oil, "co2": Group.pp_oil},
                "external": {"default": Group.pp_thermal, "co2": Group.pp_oil},
            },
            Carrier.onwind_1: {"internal": Group.wind_onshore, "external": Group.wind},
            Carrier.onwind_2: {"internal": Group.wind_onshore, "external": Group.wind},
            Carrier.onwind_3: {"internal": Group.wind_onshore, "external": Group.wind},
            Carrier.onwind_4: {"internal": Group.wind_onshore, "external": Group.wind},
            Carrier.process_emissions: {
                "internal": {"default": Group.industry, "co2": None},
                "external": Group.industry,
            },
            Carrier.process_emissions_cc: {
                "internal": {"default": Group.industry_cc, "co2": None},
                "external": Group.industry_cc,
            },
            Carrier.pemfc_rural_services_decentral_ch4_smr: {
                "internal": {
                    "default": Group.misc,
                    "co2": None,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "co2": Group.heat_decentral,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.pemfc_rural_residential_h2: {
                "internal": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.ch4_boiler_residential_rural: {
                "internal": {"default": Group.heat_decentral, "co2": None},
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.heat_pump_residential_rural_ground: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.chp_rural_residential_micro_ch4: {
                "internal": {
                    "default": Group.ch4,
                    "capacity": Group.misc,
                    "co2": None,
                },
                "external": {
                    "default": Group.pp_thermal,
                    "capacity": Group.misc,
                    "co2": Group.heat_decentral,
                },
            },
            Carrier.oil_boiler_rural_residential: {
                "internal": {"default": Group.heat_decentral, "co2": None},
                "external": Group.heat_decentral,
            },
            Carrier.oil_boiler_urban_residential: {
                "internal": Group.heat_decentral,
                "external": Group.heat_decentral,
            },
            Carrier.resistive_heater_rural_residential: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            "residential rural solid biomass boiler": Group.heat_decentral,
            Carrier.pemfc_urban_residential_decentral_ch4_smr: {
                "internal": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                    "co2": None,
                },
                "external": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.pemfc_urban_residential_decentral_h2_smr: {
                "internal": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.heat_pump_air_urban_residential_decentral: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.ch4_boiler_urban_decentral_residential: Group.heat_decentral,
            Carrier.chp_urban_decentral_micro_ch4: {
                "internal": {"default": Group.ch4, "capacity": Group.misc},
                "external": {
                    "default": Group.pp_thermal,
                    "capacity": Group.misc,
                },
            },
            Carrier.resistive_heater_urban_decentral_residential: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.ror: Group.ror,
            Carrier.pemfc_rural_services_ch4_smr: {
                "internal": {
                    "default": Group.misc,
                    "co2": None,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "co2": Group.heat_decentral,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.pemfc_rural_services_h2_smr: {
                "internal": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.ch4_boiler_services_rural: {
                "internal": {"default": Group.heat_decentral, "co2": None},
                "external": Group.heat_decentral,
            },
            Carrier.heat_pump_ground_services_rural: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.chp_rural_services_micro_ch4: {
                "internal": {
                    "default": Group.ch4,
                    "capacity": Group.misc,
                    "co2": None,
                },
                "external": {
                    "default": Group.pp_thermal,
                    "capacity": Group.misc,
                    "co2": Group.heat_decentral,
                },
            },
            Carrier.oil_boiler_rural_services: {
                "internal": {"default": Group.heat_decentral, "co2": None},
                "external": Group.heat_decentral,
            },
            Carrier.resistive_heater_rural_services: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            "services rural solid biomass boiler": Group.heat_decentral,
            Carrier.pemfc_urban_services_decentral_ch4_smr: {
                "internal": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.pemfc_urban_services_decentral_h2: {
                "internal": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
                "external": {
                    "default": Group.misc,
                    "electricity": Group.pp_thermal,
                },
            },
            Carrier.heat_pump_air_services_urban_decentral: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.ch4_boiler_urban_decentral_services: Group.heat_decentral,
            Carrier.chp_urban_decentral_services_micro_ch4: {
                "internal": {"default": Group.ch4, "capacity": Group.misc},
                "external": {
                    "default": Group.pp_thermal,
                    "capacity": Group.misc,
                },
            },
            Carrier.resistive_heater_services_urban_decentral: {
                "internal": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_decentral,
                    "electricity": Group.heat,
                },
            },
            Carrier.solar_rooftop: {"internal": Group.pv_rooftop, "external": Group.pv},
            Carrier.solar_utility: {"internal": Group.pv_utility, "external": Group.pv},
            Carrier.heat_pump_air_urban_central: {
                "internal": {
                    "default": Group.heat_district,
                    "district_heat": Group.heat_pump,
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_district,
                    "district_heat": Group.heat_pump,
                    "electricity": Group.heat,
                },
            },
            Carrier.chp_urban_central_coal_cc: {
                "internal": {"default": Group.coal, "co2": Group.chp_coal_cc},
                "external": {"default": Group.pp_thermal, "co2": Group.chp_coal_cc},
            },
            Carrier.chp_urban_central_coal_heat_cc: {
                "internal": {
                    "default": Group.chp_coal_cc,
                    "district_heat": Group.chp_coal,
                },
                "external": {
                    "default": Group.chp_coal_cc,
                    "district_heat": Group.chp_coal,
                },
            },
            Carrier.chp_urban_central_coal: {
                "internal": {"default": Group.coal, "co2": Group.chp_coal},
                "external": {"default": Group.pp_thermal, "co2": Group.chp_coal},
            },
            Carrier.chp_urban_central_coal_heat: Group.chp_coal,
            Carrier.chp_urban_central_ch4_cc: {
                "internal": {
                    "default": Group.ch4,
                    "capacity": Group.chp_electricity,
                    "co2": Group.chp_ch4_cc,
                },
                "external": {
                    "default": Group.pp_thermal,
                    "capacity": Group.chp_electricity,
                    "co2": Group.chp_ch4_cc,
                },
            },
            Carrier.chp_urban_central_ch4_heat_cc: {
                "internal": {
                    "default": Group.chp_ch4_cc,
                    "capacity": Group.heat_district,
                    "district_heat": Group.chp_ch4,
                },
                "external": {
                    "default": Group.chp_ch4_cc,
                    "capacity": Group.heat_district,
                    "district_heat": Group.chp_ch4,
                },
            },
            Carrier.chp_urban_central_ch4: {  # fixme: removable?
                "internal": {
                    "default": Group.ch4,
                    "capacity": Group.chp_electricity,
                    "co2": None,
                },
                "external": {
                    "default": Group.pp_thermal,
                    "capacity": Group.chp_electricity,
                    "co2": Group.chp_ch4,
                },
            },
            "urban central gas CHP heat": {  # fixme: removable?
                "internal": {
                    "default": Group.chp_ch4,
                    "capacity": Group.heat_district,
                    "co2": None,
                },
                "external": {"default": Group.chp_ch4, "capacity": Group.heat_district},
            },
            Carrier.ch4_boiler_urban_central: {
                "internal": {
                    "default": Group.heat_district,
                    "district_heat": "Gas Boiler",
                    "co2": None,
                },
                "external": {
                    "default": Group.heat_district,
                    "district_heat": "Gas Boiler",
                },
            },
            "urban central lignite CHP CC electric": {  # fixme: removable?
                "internal": {"default": Group.coal, "co2": Group.chp_coal_cc},
                "external": {"default": Group.pp_thermal, "co2": Group.chp_coal_cc},
            },
            "urban central lignite CHP CC heat": {  # fixme: removable?
                "internal": {
                    "default": Group.chp_coal_cc,
                    "district_heat": Group.chp_coal,
                },
                "external": {
                    "default": Group.chp_coal_cc,
                    "district_heat": Group.chp_coal,
                },
            },
            "urban central lignite CHP electric": {  # fixme: removable?
                "internal": {"default": Group.coal, "co2": Group.chp_coal},
                "external": {"default": Group.pp_thermal, "co2": Group.chp_coal},
            },
            "urban central lignite CHP heat": Group.chp_coal,  # fixme: removable?
            Carrier.resistive_heater_urban_central: {
                "internal": {
                    "default": Group.heat_district,
                    "district_heat": "Resistive Heater",
                    "electricity": Group.heat,
                },
                "external": {
                    "default": Group.heat_district,
                    "district_heat": "Resistive Heater",
                    "electricity": Group.heat,
                },
            },
            "urban central solar thermal": "Solar Thermal",
            Carrier.solar_thermal_collector_urban_central: "Solar Thermal",
            Carrier.chp_urban_central_solid_biomass: {
                "internal": {
                    "default": Group.biomass,
                    "e_production": Group.chp_biomass,
                    "e_demand": Group.chp_biomass,
                    "electricity": Group.pp_thermal,
                    "district_heat": Group.chp_biomass,
                    "biomass": Group.electrictiy,
                },
                "external": {
                    "default": Group.biomass,
                    "e_production": Group.chp_biomass,
                    "e_demand": Group.chp_biomass,
                    "electricity": Group.pp_thermal,
                    "district_heat": Group.chp_biomass,
                    "biomass": Group.electrictiy,
                },
            },
            Carrier.chp_urban_central_solid_biomass_cc: {
                "internal": {
                    "default": Group.biomass,
                    "e_production": Group.chp_biomass,
                    "e_demand": Group.chp_biomass,
                    "electricity": Group.pp_thermal,
                    "district_heat": Group.chp_biomass,
                    "biomass": Group.electrictiy,
                },
                "external": {
                    "default": Group.biomass,
                    "e_production": Group.chp_biomass,
                    "e_demand": Group.chp_biomass,
                    "electricity": Group.pp_thermal,
                    "district_heat": Group.chp_biomass,
                    "biomass": Group.electrictiy,
                },
            },
            Carrier.solid_biomass_boiler_urban_central: Group.solid_biomass_boiler,
            Carrier.solid_biomass_boiler_urban_central_cc: Group.solid_biomass_boiler,
            Carrier.water_tanks_charger_urban_central: Group.storage_in,
            Carrier.water_tanks_discharger_urban_central: Group.storage_out,
            Carrier.value_lost_load: Group.power_disconnect,
        }
    )

    bus_carrier = frozendict(
        {
            BusCarrier.AC: Group.electrictiy,
            BusCarrier.CH4: Group.ch4,
            BusCarrier.H2: Group.h2,
            BusCarrier.TRANSPORT_P: Group.oil,
            BusCarrier.TRANSPORT_P_LONG: Group.oil,
            BusCarrier.FT: Group.oil,
            BusCarrier.FT_1: Group.oil,
            BusCarrier.FT_2: Group.oil,
            BusCarrier.SOLID_BIOMASS: Group.biomass,
            # BusCarrier.HEAT_RURAL_SERVICES: Group.heat,
            # BusCarrier.HEAT_RURAL_RESIDENTIAL: Group.industry,
            # BusCarrier.HEAT_URBAN_CENTRAL: Group.heat_district,
            # BusCarrier.HEAT_URBAN_SERVICES: Group.heat_district,
            # BusCarrier.HEAT_URBAN_RESIDENTIAL: Group.heat_district,
        }
    )


TITLE_SUFFIX: str = " {location} in {unit}"

UNITS: frozendict = frozendict(
    {
        "W": 1e-6,
        "Wh": 1e-6,
        "KW": 1e-3,
        "kW": 1e-3,  # alias
        "KWh": 1e-3,
        "kWh": 1e-3,  # alias
        "MW": 1,  # model base unit
        "MWh": 1,  # model base unit
        "GW": 1e3,
        "GWh": 1e3,
        "TW": 1e6,
        "TWh": 1e6,
        "PW": 1e9,
        "PWh": 1e9,
        "currency": 1,
        "EUR": 1,  # base currency
        "t": 1,
        "kt": 1e3,
        "Mt": 1e6,
    }
)

# transmission technologies
TRANSMISSION_CARRIER: tuple = (
    "",
    "DC",
    Carrier.gas_pipepline,
    Carrier.gas_pipepline_new,
    Carrier.h2_pipeline,
    Carrier.h2_pipeline_retro,
    Carrier.h2_pipeline_kernnetz,
)


class TradeTypes:
    """Collect trade type names."""

    LOCAL: str = "local"  # same node
    DOMESTIC: str = "domestic"  # same country, but different node
    FOREIGN: str = "foreign"  # different country


class COLOUR:
    """Container to collect colour codes in hex format."""

    coral: str = "#E8B5B1"
    raspberry: str = "#961454"
    salmon: str = "#E19990"
    rose: str = "#D5A1BB"
    peach: str = "#EBBFBA"

    red: str = "#CA0638"
    red_chestnut: str = "#96332C"
    red_bright: str = "#E53212"
    red_deep: str = "#B20633"
    red_fire: str = "#E63313"

    green: str = "#3C703E"
    green_light: str = "#509554"
    green_ocean: str = "#3DCCBF"
    green_mint: str = "#B0D4B2"
    green_sage: str = "#82B973"
    turquoise: str = "#e8e8e8"

    grey_light: str = "#ECECEC"
    grey_dark: str = "#535353"
    grey_charcoal: str = "#485055"
    grey_deep: str = "#3C3C3C"
    grey_cool: str = "#919699"
    grey_silver: str = "#D0D0D0"
    grey_neutral: str = "#9F9F9F"

    black: str = "#000000"

    brown: str = "#C58000"
    brown_dark: str = "#b37400"
    brown_sallow: str = "#bf9c5c"
    brown_light: str = "#e8cc99"
    brown_deep: str = "#4d3200"

    blue_pastel: str = "#B5C9D5"
    blue_moonstone: str = "#3DACBF"
    blue_dark: str = "#5F5F5F"
    blue_persian: str = "#0064A2"
    blue_celestial: str = "#4F8FCD"
    blue_cerulean: str = "#005082"
    blue_sky: str = "#99C1DA"
    blue_lavender: str = "#636EFA"

    orange: str = "#FF6600"
    orange_mellow: str = "#FECB52"

    yellow_bright: str = "#FED500"
    yellow_vivid: str = "#FEC500"
    yellow_canary: str = "#FFDE53"
    yellow_golden: str = "#FFB200"


# cannot freeze, because plotly manipulates the dictionary
COLOUR_SCHEME_BMK: dict = {
    # dark blue - coal
    Group.coal: COLOUR.blue_dark,
    Group.pp_coal: COLOUR.blue_dark,
    Group.chp_coal: COLOUR.blue_dark,
    Group.chp_coal_cc: COLOUR.turquoise,
    "Coal Import": COLOUR.blue_dark,
    # red - oil
    "Oil": COLOUR.red,
    "Oil PP": COLOUR.red,
    "Fischer-Tropsch": COLOUR.red,
    "Oil Import": COLOUR.red,
    # dark green - biogas
    "Biogas": COLOUR.green,
    "Bio Methane Processing": COLOUR.green,
    # light green - biomass
    Group.biomass: COLOUR.green_light,
    "Wet Biomass": COLOUR.green_light,
    "Solid Biomass": COLOUR.green_light,
    Group.chp_biomass: COLOUR.green_light,
    Group.solid_biomass_boiler: COLOUR.green,
    # brown - methane
    "Methane": COLOUR.brown,
    "Gas PP": COLOUR.brown,
    Group.chp_ch4: COLOUR.brown_dark,
    "Methanation": COLOUR.brown,
    "Gas Boiler": COLOUR.brown_sallow,
    Group.chp_ch4_cc: COLOUR.brown_light,
    "Methane Import": COLOUR.brown,
    "Thermal Powerplants": COLOUR.brown,
    "OCGT": COLOUR.brown,
    # light grey - hydrogen
    "Hydrogen": COLOUR.blue_pastel,
    "Electrolysis": COLOUR.blue_pastel,
    "SMR": COLOUR.yellow_bright,
    "Hydrogen Tube Storage": COLOUR.blue_pastel,
    "Hydrogen Underground Storage": COLOUR.grey_charcoal,
    "SMR CC": COLOUR.grey_cool,
    "Hydrogen Import": COLOUR.blue_pastel,
    # teal - wind power
    "Wind Power": COLOUR.blue_moonstone,
    "Onshore": COLOUR.blue_moonstone,
    "Offshore": COLOUR.green_ocean,
    # blue - hydro
    "Hydro Power": COLOUR.blue_persian,
    "Run-of-River": COLOUR.blue_persian,
    "Reservoir": COLOUR.blue_cerulean,
    "Pumped Hydro Storage": COLOUR.red_chestnut,
    "Inflow Hydro Storage": COLOUR.blue_cerulean,
    # blue - heat
    "Resistive Heater": COLOUR.blue_persian,
    "Heat Pump": COLOUR.blue_celestial,
    "Fuel Cell (Heat)": COLOUR.blue_pastel,
    "Demand": COLOUR.grey_neutral,
    # yellow - solar
    "Solar Power": COLOUR.yellow_bright,
    "Photovoltaics": COLOUR.yellow_bright,
    "PV-Utility": COLOUR.yellow_bright,
    "PV-Rooftop": COLOUR.yellow_vivid,
    "Solar Thermal": COLOUR.yellow_canary,
    # red - nuclear
    "Nuclear": COLOUR.orange,
    "Nuclear Power": COLOUR.orange,
    # light blue - electricity
    "Electricity": COLOUR.blue_celestial,
    "Electricity CHP": COLOUR.blue_celestial,
    "Battery Storage": COLOUR.coral,
    "Electricity Import": COLOUR.blue_celestial,
    "Electricity OCGT": COLOUR.blue_celestial,
    # purple - heat supply
    "District Heat": COLOUR.raspberry,
    "Decentral Heat": COLOUR.salmon,
    "Heat": COLOUR.rose,
    # light pink
    "HH and Services (Heat)": COLOUR.salmon,
    # orange - ambient heat
    "Ambient Heat": COLOUR.red_bright,
    # light green - DAC, Fuel cell
    "Heat for DAC": COLOUR.green_mint,
    "Direct Air Capture": COLOUR.green_mint,
    "Fuel Cell": COLOUR.blue_celestial,
    "Hydrogen Fuel Cell": COLOUR.blue_celestial,  # fixme: yellow in old Toolbox?!
    # grey - losses, misc
    "Transformation Losses": COLOUR.grey_silver,
    "Miscellaneous": COLOUR.grey_dark,
    "Losses": COLOUR.grey_silver,
    "Storage": COLOUR.grey_light,
    "co2 vent": COLOUR.grey_silver,
    Group.import_foreign: COLOUR.grey_silver,
    Group.export_foreign: COLOUR.grey_silver,
    Group.import_domestic: COLOUR.blue_lavender,
    Group.export_domestic: COLOUR.orange_mellow,
    Group.power_disconnect: COLOUR.grey_dark,
    "Grid Losses": COLOUR.grey_silver,
    # Sectors
    "Industry": COLOUR.red,
    "Households & Services": COLOUR.grey_neutral,
    "Transport": COLOUR.grey_deep,
    "Industry CC": COLOUR.red_deep,
    # Time Series
    "Inflexible Demand": COLOUR.black,
    "Base Load": COLOUR.yellow_golden,
    "Storage In": COLOUR.green_sage,
    "Storage Out": COLOUR.green_sage,
    "Net Import": COLOUR.grey_silver,
    "Net Export": COLOUR.grey_silver,
    Group.import_global: COLOUR.grey_silver,
    Group.global_market: COLOUR.blue_lavender,
    "State of Charge": COLOUR.blue_sky,
    "Max State of Charge": COLOUR.grey_silver,
    "Accumulated Turbining": COLOUR.blue_celestial,
    "Accumulated Pumping": COLOUR.peach,
    "Accumulated Outflow Spill": COLOUR.grey_silver,
    "Accumulated Natural Inflow": COLOUR.blue_cerulean,
    "Residualload": COLOUR.red_fire,
    "Waste": COLOUR.raspberry,
    "Methanolisation": COLOUR.salmon,
    "Methane Compression": COLOUR.brown,
    "Hydrogen Compression": COLOUR.blue_pastel,
    "Haber-Bosch": COLOUR.red,
    "Agriculture": COLOUR.green_light,
    "Distribution Grid": COLOUR.grey_silver,
    "Ammonia Cracking": COLOUR.red_chestnut,
    "Sabatier": COLOUR.yellow_canary,
    "Synth. Fuels": COLOUR.red,
    "Methanol Steam Reforming": COLOUR.salmon,
    "H2 from Solid Biomass": COLOUR.green_mint,
}

ALIAS_COUNTRY: frozendict = frozendict(
    {
        "EU": "Europe",
        "AL": "Albania",
        "AT": "Austria",
        "BA": "Bosnia and Herzegovina",
        "BE": "Belgium",
        "BG": "Bulgaria",
        "CH": "Switzerland",
        "CZ": "Czech Republic",
        "DE": "Germany",
        "DK": "Denmark",
        "EE": "Estonia",
        "ES": "Spain",
        "FI": "Finland",
        "FR": "France",
        "GB": "Great Britain",
        "GR": "Greece",
        "HR": "Croatia",
        "HU": "Hungary",
        "IE": "Ireland",
        "IT": "Italy",
        "LT": "Lithuania",
        "LU": "Luxembourg",
        "LV": "Latvia",
        "ME": "Montenegro",
        "MK": "North Macedonia",
        "NL": "Netherlands",
        "NO": "Norway",
        "PL": "Poland",
        "PT": "Portugal",
        "RO": "Romania",
        "RS": "Serbia",
        "SE": "Sweden",
        "SI": "Slovenia",
        "SK": "Slovakia",
    }
)
ALIAS_COUNTRY_REV: frozendict = frozendict({v: k for k, v in ALIAS_COUNTRY.items()})

ALIAS_REGION: frozendict = frozendict(
    {
        "AT11": "Burgenland (AT)",
        "AT12": "Lower Austria (AT)",
        "AT13": "Vienna (AT)",
        "AT21": "Carinthia (AT)",
        "AT22": "Styria (AT)",
        "AT31": "Upper Austria (AT)",
        "AT32": "Salzburg (AT)",
        "AT33": "Tyrol (AT)",
        "AT333": "East Tyrol (AT)",
        "AT34": "Vorarlberg (AT)",
        # German NUTS1
        "DE1": "Baden-Württemberg",
        "DE2": "Bavaria",
        "DE3": "Berlin",
        "DE4": "Brandenburg",
        "DE5": "Bremen",
        "DE6": "Hamburg",
        "DE7": "Hesse",
        "DE8": "Mecklenburg-Western Pomerania",
        "DE9": "Lower Saxony",
        "DEA": "North Rhine-Westphalia",
        "DEB": "Rhineland-Palatinate",
        "DEC": "Saarland",
        "DED": "Saxony",
        "DEE": "Saxony-Anhalt",
        "DEF": "Schleswig-Holstein",
        "DEG": "Thuringia",
    }
)
ALIAS_REGION_REV: frozendict = frozendict({v: k for k, v in ALIAS_REGION.items()})

ALIAS_LOCATION: frozendict = ALIAS_COUNTRY | ALIAS_REGION
ALIAS_LOCATION_REV: frozendict = frozendict({v: k for k, v in ALIAS_LOCATION.items()})


try:
    esmtools_version = importlib.metadata.version("esmtools")
except PackageNotFoundError:
    esmtools_version = "esmtools not installed."

try:
    git_hash = check_output("git rev-parse HEAD").strip().decode("utf-8")
except (CalledProcessError, FileNotFoundError):
    git_hash = "Not a git repo."

RUN_META_DATA = {
    "esmtools": f"{esmtools_version} ({git_hash})",
    # todo: include model version and git hash
    # todo: include model input data version and git hash
}
