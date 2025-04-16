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



