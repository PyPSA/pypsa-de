c_e_battery = ['battery', 'home battery']
c_e_bev = ['EV battery']
c_e_DSM = ['DSM']
c_e_PHS = ['PHS']
c_e_H2 = ['H2 Store']
c_e_tank = ['rural water tanks', 'urban decentral water tanks', 'urban central water tanks']

c_power = ['AC', 'DC']
c_battery = ['battery', 'home battery']

c_pg_natgas = ['urban central gas CHP CC', 'urban central gas CHP', 'CCGT', 'OCGT']
c_pg_coal_oil = ['urban central lignite CHP', 'urban central oil CHP', 'coal']
c_pg_h2 = ['H2 Fuel Cell', 'H2 OCGT', 'urban central H2 CHP']
c_pg_retrofit_h2 = ['H2 retrofit OCGT', 'urban central H2 retrofit CHP']

c_pg_import = c_power
c_pg_inner_import = c_power

c_pg_battery = ['battery discharger', 'home battery discharger']
c_pg_phs = ['PHS']
c_pg_dsm = ['DSM dispatch']

c_pg_biomass = ['urban central solid biomass CHP', 'urban central solid biomass CHP CC']
c_pg_waste = ['waste CHP', 'waste CHP CC']
c_pg_water = ['ror', 'hydro']
c_pg_biogas = ['biogas']


c_pg_onwind = ['onwind']
c_pg_offwind = ['offwind-dc', 'offwind-ac']
c_pg_pv = ['solar', 'solar-hsat', 'solar rooftop']
c_pg_wind = ['onwind', 'offwind-ac', 'offwind-dc']



c_pc_ghd = ['electricity']
c_pc_industry = ['industry electricity']
c_pc_agriculture = ['agriculture electricity']

c_pc_bev_charger = ['BEV charger']
c_pc_electrolysis = ['H2 Electrolysis']

c_pc_heat_pump = ['urban central air heat pump', 'rural air heat pump', 'rural ground heat pump', 'urban decentral air heat pump']
c_pc_dac = ['DAC']
c_pc_resistive = ['urban central resistive heater', 'rural resistive heater', 'urban decentral resistive heater']

c_pc_battery = ['battery charger', 'home battery charge']
c_pc_phs = ['PHS']
c_pc_dsm = ['DSM store']
c_pc_export = c_power
c_pc_inner_export = c_power

map_name = {
  'c_pg_natgas': 'natgas generation',
  'c_pg_h2': 'hydrogen generation',
  'c_pg_retrofit_h2': 'retrofit H2 generation',
  'c_pg_import': 'import power',
  'c_pg_inner_import': 'transmission(receive)',
  'c_pg_battery': 'battery dispatch',
  'c_pg_phs': 'PHS dispatch',
  'c_pg_dsm': 'DSM dispatch',
  'c_pg_biomass': 'biomass generation',
  'c_pg_waste': 'waste generation',
  'c_pg_water': 'water',
  'c_pg_onwind': 'onwind',
  'c_pg_offwind': 'offwind',
  'c_pg_pv': 'solar',
  'c_pg_biogas': 'biogas',
  'c_pg_coal_oil': 'coal/oil',

  'c_pc_ghd': 'electricity GHD/Private Haushalte',
  'c_pc_industry': 'electricity industry',
  'c_pc_agriculture': 'electricity agriculture',
  'c_pc_bev_charger': 'BEV charger',
  'c_pc_electrolysis': 'electrolysis',
  'c_pc_heat_pump': 'heat pump',
  'c_pc_dac': 'DAC',
  'c_pc_resistive': 'resistive heat',
  'c_pc_battery': 'battery charger',
  'c_pc_phs': 'PHS charger',
  'c_pc_dsm': 'DSM charger',
  'c_pc_export': 'export power',
  'c_pc_inner_export': 'transmission(transmit)',

  'c_e_battery': 'battery',
  'c_e_bev': 'EV battery',
  'c_e_DSM': 'DSM',
  'c_e_PHS': 'PHS',
  'c_e_H2': 'H2',
  'c_e_tank': 'water tank'

}

map_color = {
  'c_pg_natgas': '#4A4A4A',
  'c_pg_h2': '#6F6F6F',
  'c_pg_retrofit_h2': '#939393',
  'c_pg_coal_oil': '#B7B7B7',
  'c_pg_inner_import': '#765FB4',
  'c_pg_import': '#5438A1',
  'c_pg_battery': '#3777B4',
  'c_pg_phs': '#5F92C3',
  'c_pg_dsm': '#87AED2',
  'c_pg_biomass': '#48A299',
  'c_pg_biogas': '#B6DAD6',
  'c_pg_waste': '#6DB4AD',
  'c_pg_water': '#91C7C2',
  'c_pg_onwind': '#4F9C59',
  'c_pg_offwind': '#72B07A',
  'c_pg_pv': '#95C49B',

  'c_pc_ghd': '#CB9B47',
  'c_pc_industry': '#D6AF6B',
  'c_pc_agriculture': '#E0C390',
  'c_pc_bev_charger': '#B9703C',
  'c_pc_electrolysis': '#C78D63',
  'c_pc_heat_pump': '#9C372A',
  'c_pc_dac': '#C3877F',
  'c_pc_resistive': '#B05F54',
  'c_pc_battery': '#A3397E',
  'c_pc_phs': '#B66198',
  'c_pc_dsm': '#C888B2',
  'c_pc_inner_export': '#BBAFD9',
  'c_pc_export': '#9887C7',

  'c_e_battery': '#3777B4',
  'c_e_bev': '#B9703C',
  'c_e_DSM': '#C888B2',
  'c_e_PHS': '#5F92C3',
  'c_e_H2': '#C78D63',
  'c_e_tank':'#4F9C59',
}

map_name_color = {map_name[key]: value for key, value in map_color.items() if key in map_name}









## power supply
onwind = ['onwind']
offwind = ['offwind-ac', 'offwind-dc', 'offwind-float']
solar = ['solar-hsat', 'solar']
solarRooftop = ['solar rooftop']
phs = ['PHS']
hydro = ['hydro']
batteryDischarger = ['battery discharger', 'home battery discharger']
gasCHP = ['urban central CHP', 'urban central CHP CC']
biomassCHP = ['urban central solid biomass CHP', 'urban central solid biomass CHP CC']
geothermalORC = ['geothermal organic rankine cycle']
h2FC =['H2 Fuel Cell']
v2g = ['V2G']
gasPower = ['OCGT', 'CCGT']+ gasCHP
h2Power = ['H2 turbine']
gasH2Power= gasPower  + h2Power
power = ['AC', 'DC']
distribution = ['electricity distribution grid'] # there is two

## power use
h2Electrolysis = ['H2 Electrolysis'] # load
batteryCharger = ['battery charger', 'home battery charger']
power2warm = ['DAC', 'urban central air heat pump', 'urban central resistive heater', 'rural air heat pump', 'rural ground heat pump', 'rural resistive heater', 'urban decentral air heat pump', 'urban decentral resistive heater']
bevCharger = ['BEV charger']
powerUse = ['electricity'] #lo a d
powerUseIndustry =['industry electricity', 'agriculture electricity', ]
landTransportEV = ['land transport EV'] # load

allCharge = bevCharger + batteryCharger

## battery
battery = ['battery', 'home battery']
batteryCharger = batteryCharger
batteryDischarger = batteryDischarger

#3 biomass
biomassCHP = biomassCHP
biomassBoiler = ['rural biomass boiler', 'urban decentral biomass boiler']

## geothermal
geothermalORC = geothermalORC
geothermalWarm = ['geothermal district heat']

## ror
ror = ['ror'] # generator

## warm
tank = ['rural water tanks', 'urban decentral water tanks', 'urban decentral water tanks', 'urban central water tanks']
tankDischarger = ['rural water tanks discharger', 'urban decentral water tanks discharger', 'urban central water tanks discharge'],
tankCharger = ['rural water tanks charger', 'urban decentral water tanks charger', 'urban central water tanks charger']

heatPump = ['rural ground heat pump','rural air heat pump','urban decentral air heat pump', 'urban central air heat pump']
resistiveHeat = ['rural resistive heater','urban decentral resistive heater','urban central resistive heate']
dac = ['DAC']
gasBoiler= ['rural gas boiler','urban central gas boiler']
geothermalWarm = geothermalWarm
biomassBoiler = biomassBoiler

warmLoad = ['urban decentral heat', 'rural heat', 'agriculture heat', 'urban central heat', 'low-temperature heat for industry']

## h2
h2FC =['H2 Fuel Cell']
h2Turbine = ['H2 turbine']
h2Store=['H2 Store']


## filter
def filter2045 (df):
  return df.index.str.endswith('-2045')

def filter2040 (df):
  return df.index.str.endswith('-2040')

def filter2035 (df):
  return df.index.str.endswith('-2035')

def filter2030 (df):
  return df.index.str.endswith('-2030')

mapName = {
  'onwind': 'Offshore-Wind',
  'offwind': 'Onshore-Wind',
  'solar': 'PV-Freifläche',
  'solarRooftop': 'PV-Aufdach',
  'phs': 'PHS',
  'phsCharge': 'PHS-Aufladung',
  'hydro': 'Wasserkraft',
  'batteryDischarger': 'Batterie',
  'gasH2Power':'Gas/Wasserstoff',
  'biomassCHP': 'Biomasse',
  'geothermalORC': 'Geothermie',
  'v2g': 'V2G',
  'importPower': 'Stromimport',
  'otherSupply': 'sonstige Versorgung',

  'h2Electrolysis': 'Elektrolyse',
  'batteryCharger': 'Batterie aufladen',
  'power2warm': 'Wärme',
  'bevCharger': 'BEV',
  'powerUse': 'Elektrizität',
  'exportPower': 'Stromexport',
  'otherUse': 'sonstiger Stromverbrauch',
  'Smart': 'Smart',
  'Gas': 'Gas'
}

## embellishment
colors= [ '#002e4f', '#11337A', '#2274A5', '#3f8d9e', '#44AF69',
         '#5e9a7a', '#F1C40F', '#d4a017', '#FE7F2D', '#b86f33'
         '#9370DB'
         ]

lightColors = ['#28A5FF' '#68A0D2', '#7EBEE4', '#97CAC5', '#9ECAB3', 
               '#ADCEDC', '#F8E185', '#F1D384', '#FFBF96', '#E2B792',
               '#C9B7ED']
mapColor= {
  'onwind': '#11337A',
  'offwind': '#2274A5',
  'solar': '#F1C40F',
  'solarRooftop': '#d4a017',
  # 'phs': 'PHS',
  # 'hydro': 'Wasserkraft',
  'batteryDischarger': '#44AF69',
  'gasH2Power':'#b86f33',
  'biomassCHP': '#5e9a7a',
  # 'geothermalORC': 'Geothermie',
  'v2g': '#FE7F2D',
  'importPower': '#002e4f',
  'otherSupply': '#9C9C9C',

  'h2Electrolysis': '#E2B792',
  'batteryCharger': '#9ECAB3',
  'power2warm': '#F8E185',
  'bevCharger': '#FFBF96',
  'powerUse': '#ADCEDC',
  'exportPower': '#68A0D2',
  'otherUse': '#D6D6D6',
  'Gas': '#333333',
  'Smart': '#333333'
}

mapNameColor = {mapName[key]: value for key, value in mapColor.items() if key in mapName}



