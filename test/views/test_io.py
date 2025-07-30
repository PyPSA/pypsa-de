import pandas as pd
import pytest
import yaml
from esmtools.fileio import prepare_co2_emissions, prepare_costs
from numpy import nan


@pytest.mark.unit
def test_prepare_costs(result_path, with_mock_config, expected_costs):
    """Check if the cost calculations return as expected."""
    costs = prepare_costs(result_path, 25)
    costs = costs.loc[["2015"]].droplevel("year")
    pd.testing.assert_frame_equal(costs, expected_costs)


@pytest.mark.unit
def test_read_historical_country_emissions(
    result_path, with_mock_config, expected_co2_share
):
    """Check if the historical country emissions return as expected."""
    p = 15  # encoding precision
    data_path = result_path / ".." / "data"
    co2_share = prepare_co2_emissions(data_path, "THI")
    co2_share.columns = expected_co2_share.columns  # they are not the same
    assert co2_share.round(p).equals(expected_co2_share.round(p))


@pytest.fixture(scope="module")
def expected_costs() -> pd.DataFrame:
    """
    Return the expected costs.

    Expected data is the result from
    Toolbox.Import_functions.prepare_costs and for 2030 only.
    """
    expected_data = {
        "CO2 intensity": {
            "BEV vehicle": 0.0,
            "CCGT": 0.0,
            "CH4 pipeline": 0.0,
            "CH4 pipeline DIN1000": 0.0,
            "CH4 pipeline DIN1400": 0.0,
            "CH4 pipeline DIN600": 0.0,
            "CH4 pipeline retrofit": 0.0,
            "CH4 pipeline retrofit DIN1000": 0.0,
            "CH4 pipeline retrofit DIN1400": 0.0,
            "CH4 pipeline retrofit DIN600": 0.0,
            "CH4-powered PEMFC with internal SMR": 0.0,
            "CNG vehicle": 0.0,
            "DSM-household: charger": 0.0,
            "DSM-household: discharger": 0.0,
            "DSM-household: storage": 0.0,
            "DSM_Aluminum": 0.0,
            "DSM_Cement": 0.0,
            "DSM_Chlorine": 0.0,
            "DSM_ChlorineM": 0.0,
            "DSM_ChlorineQ": 0.0,
            "DSM_Household (turn-off)": 0.0,
            "DSM_Paper": 0.0,
            "DSM_PaperH": 0.0,
            "DSM_PaperS": 0.0,
            "DSM_Steel": 0.0,
            "FCEV vehicle": 0.0,
            "Fischer-Tropsch": 0.0,
            "H2 pipeline": 0.0,
            "H2-powered PEMFC with internal SMR": 0.0,
            "HEV vehicle": 0.0,
            "HVAC overhead": 0.0,
            "HVDC inverter pair": 0.0,
            "HVDC overhead": 0.0,
            "HVDC submarine": 0.0,
            "ICE vehicle": 0.0,
            "OCGT": 0.0,
            "PHEV vehicle": 0.0,
            "PHS": 0.0,
            "SMR": 0.0,
            "SMR CC": 0.0,
            "battery inverter": 0.0,
            "battery storage": 0.0,
            "biogas": 0.0,
            "biogas import": 0.0,
            "biomass": 0.0,
            "biomass CHP capture": 0.0,
            "central CHP": 0.0,
            "central air-sourced heat pump": 0.0,
            "central gas CHP": 0.0,
            "central gas CHP CC": 0.0,
            "central gas boiler": 0.0,
            "central resistive heater": 0.0,
            "central solar thermal": 0.0,
            "central solid biomass CHP": 0.0,
            "central solid biomass CHP CC": 0.0,
            "central water tank storage": 0.0,
            "coal": 0.3379999999999999,
            "coal CHP power plant": 0.0,
            "coal CHP power plant (CC)": 0.0,
            "coal power plant": 0.0,
            "coal power plant (CC)": 0.0,
            "decentral CHP": 0.0,
            "decentral air-sourced heat pump": 0.0,
            "decentral coal boiler": 0.0,
            "decentral gas boiler": 0.0,
            "decentral ground-sourced heat pump": 0.0,
            "decentral oil boiler": 0.0,
            "decentral resistive heater": 0.0,
            "decentral solar thermal": 0.0,
            "decentral solid biomass boiler": 0.0,
            "decentral water tank storage": 0.0,
            "direct air capture": 0.0,
            "electrolysis": 0.0,
            "electrolysis_soec": 0.0,
            "fischer tropsch import": 0.0,
            "fuel cell": 0.0,
            "gas": 0.2044,
            "geothermal": 0.026,
            "helmeth": 0.0,
            "hydro": 0.0,
            "hydrogen": 0.0,
            "hydrogen pore storage underground": 0.0,
            "hydrogen storage": 0.0,
            "hydrogen storage underground": 0.0,
            "industry CC": 0.0,
            "lignite": 0.4062,
            "lignite CHP power plant": 0.0,
            "lignite CHP power plant (CC)": 0.0,
            "lignite power plant": 0.0,
            "lignite power plant (CC)": 0.0,
            "methanation": 0.0,
            "micro CHP": 0.0,
            "nuclear": 0.0,
            "offwind": 0.0,
            "offwind-ac-connection-submarine": 0.0,
            "offwind-ac-connection-underground": 0.0,
            "offwind-ac-station": 0.0,
            "offwind-dc-connection-submarine": 0.0,
            "offwind-dc-connection-underground": 0.0,
            "offwind-dc-station": 0.0,
            "oil": 0.2847,
            "oil power plant": 0.0,
            "onwind-1": 0.0,
            "onwind-2": 0.0,
            "onwind-3": 0.0,
            "onwind-4": 0.0,
            "onwind-landcosts": 0.0,
            "retrofitting I": 0.0,
            "retrofitting II": 0.0,
            "ror": 0.0,
            "solar": 0.0,
            "solar-rooftop": 0.0,
            "solar-utility": 0.0,
            "solid biomass": 0.0,
            "value of lost load": 0.0,
            "water tank charger": 0.0,
            "water tank discharger": 0.0,
        },
        "FOM": {
            "BEV vehicle": 0.0,
            "CCGT": 2.54,
            "CH4 pipeline": 0.8,
            "CH4 pipeline DIN1000": 0.0,
            "CH4 pipeline DIN1400": 0.0,
            "CH4 pipeline DIN600": 0.0,
            "CH4 pipeline retrofit": 0.8,
            "CH4 pipeline retrofit DIN1000": 0.0,
            "CH4 pipeline retrofit DIN1400": 0.0,
            "CH4 pipeline retrofit DIN600": 0.0,
            "CH4-powered PEMFC with internal SMR": 6.091911217,
            "CNG vehicle": 0.0,
            "DSM-household: charger": 0.0,
            "DSM-household: discharger": 0.0,
            "DSM-household: storage": 0.0,
            "DSM_Aluminum": 0.0,
            "DSM_Cement": 0.0,
            "DSM_Chlorine": 0.0,
            "DSM_ChlorineM": 0.0,
            "DSM_ChlorineQ": 0.0,
            "DSM_Household (turn-off)": 0.0,
            "DSM_Paper": 0.0,
            "DSM_PaperH": 0.0,
            "DSM_PaperS": 0.0,
            "DSM_Steel": 0.0,
            "FCEV vehicle": 0.0,
            "Fischer-Tropsch": 4.428571428571429,
            "H2 pipeline": 0.8,
            "H2-powered PEMFC with internal SMR": 6.275946276142857,
            "HEV vehicle": 0.0,
            "HVAC overhead": 2.0,
            "HVDC inverter pair": 2.0,
            "HVDC overhead": 2.0,
            "HVDC submarine": 2.0,
            "ICE vehicle": 0.0,
            "OCGT": 3.25,
            "PHEV vehicle": 0.0,
            "PHS": 0.985221675,
            "SMR": 7.030470914,
            "SMR CC": 5.152860411857143,
            "battery inverter": 0.4942857142857143,
            "battery storage": 0.0,
            "biogas": 0.0,
            "biogas import": 0.0,
            "biomass": 5.009624684142857,
            "biomass CHP capture": 3.0,
            "central CHP": 4.019607843,
            "central air-sourced heat pump": 0.2428571428571428,
            "central gas CHP": 3.13,
            "central gas CHP CC": 2.963418377714286,
            "central gas boiler": 3.36,
            "central resistive heater": 0.98,
            "central solar thermal": 0.61,
            "central solid biomass CHP": 2.84150798,
            "central solid biomass CHP CC": 2.7190364125714286,
            "central water tank storage": 0.3,
            "coal": 0.0,
            "coal CHP power plant": 1.838956454,
            "coal CHP power plant (CC)": 2.751448417,
            "coal power plant": 2.1,
            "coal power plant (CC)": 2.716221271,
            "decentral CHP": 3.023166023428572,
            "decentral air-sourced heat pump": 2.9857142857142858,
            "decentral coal boiler": 3.0403477960000003,
            "decentral gas boiler": 4.697590702857143,
            "decentral ground-sourced heat pump": 1.896428571428572,
            "decentral oil boiler": 3.366530612,
            "decentral resistive heater": 0.0471788967142857,
            "decentral solar thermal": 2.2846638657142857,
            "decentral solid biomass boiler": 4.717074829714285,
            "decentral water tank storage": 4.06504065,
            "direct air capture": 4.95,
            "electrolysis": 2.005714285714286,
            "electrolysis_soec": 3.8314285714285714,
            "fischer tropsch import": 0.0,
            "fuel cell": 5.0,
            "gas": 0.0,
            "geothermal": 3.815132915,
            "helmeth": 3.0,
            "hydro": 0.51369863,
            "hydrogen": 0.0,
            "hydrogen pore storage underground": 0.0,
            "hydrogen storage": 0.0,
            "hydrogen storage underground": 0.0,
            "industry CC": 2.0,
            "lignite": 0.0,
            "lignite CHP power plant": 2.3146331890000003,
            "lignite CHP power plant (CC)": 1.822590859,
            "lignite power plant": 2.24,
            "lignite power plant (CC)": 2.815190835,
            "methanation": 3.0342857142857143,
            "micro CHP": 5.6,
            "nuclear": 3.069515498,
            "offwind": 4.806273755714286,
            "offwind-ac-connection-submarine": 0.0,
            "offwind-ac-connection-underground": 0.0,
            "offwind-ac-station": 0.0,
            "offwind-dc-connection-submarine": 0.0,
            "offwind-dc-connection-underground": 0.0,
            "offwind-dc-station": 0.0,
            "oil": 0.0,
            "oil power plant": 2.4693877551428574,
            "onwind-1": 2.3814285714285712,
            "onwind-2": 2.3814285714285712,
            "onwind-3": 2.3814285714285712,
            "onwind-4": 2.3814285714285712,
            "onwind-landcosts": 0.0,
            "retrofitting I": 1.0,
            "retrofitting II": 1.0,
            "ror": 0.24,
            "solar": 4.166667,
            "solar-rooftop": 2.258535492,
            "solar-utility": 2.317247969428572,
            "solid biomass": 0.0,
            "value of lost load": 0.0,
            "water tank charger": 0.0,
            "water tank discharger": 0.0,
        },
        "VOM": {
            "BEV vehicle": 0.0,
            "CCGT": 2.91,
            "CH4 pipeline": 0.0,
            "CH4 pipeline DIN1000": 0.0,
            "CH4 pipeline DIN1400": 0.0,
            "CH4 pipeline DIN600": 0.0,
            "CH4 pipeline retrofit": 0.0,
            "CH4 pipeline retrofit DIN1000": 0.0,
            "CH4 pipeline retrofit DIN1400": 0.0,
            "CH4 pipeline retrofit DIN600": 0.0,
            "CH4-powered PEMFC with internal SMR": 0.0,
            "CNG vehicle": 0.0,
            "DSM-household: charger": 0.0,
            "DSM-household: discharger": 0.0,
            "DSM-household: storage": 0.0,
            "DSM_Aluminum": 200.0,
            "DSM_Cement": 100.0,
            "DSM_Chlorine": 200.0,
            "DSM_ChlorineM": 200.0,
            "DSM_ChlorineQ": 200.0,
            "DSM_Household (turn-off)": 20000.0,
            "DSM_Paper": 300.0,
            "DSM_PaperH": 300.0,
            "DSM_PaperS": 300.0,
            "DSM_Steel": 125.0,
            "FCEV vehicle": 0.0,
            "Fischer-Tropsch": 0.0,
            "H2 pipeline": 0.0,
            "H2-powered PEMFC with internal SMR": 0.0,
            "HEV vehicle": 0.0,
            "HVAC overhead": 0.0,
            "HVDC inverter pair": 0.0,
            "HVDC overhead": 0.0,
            "HVDC submarine": 0.0,
            "ICE vehicle": 0.0,
            "OCGT": 3.0,
            "PHEV vehicle": 0.0,
            "PHS": 0.0,
            "SMR": 0.0,
            "SMR CC": 0.0,
            "battery inverter": 0.0,
            "battery storage": 0.0,
            "biogas": 0.0,
            "biogas import": 0.0,
            "biomass": 0.0,
            "biomass CHP capture": 0.0,
            "central CHP": 0.0,
            "central air-sourced heat pump": 0.0,
            "central gas CHP": 2.2,
            "central gas CHP CC": 4.57,
            "central gas boiler": 0.0,
            "central resistive heater": 0.0,
            "central solar thermal": 0.0,
            "central solid biomass CHP": 3.142857142857143,
            "central solid biomass CHP CC": 6.5285714285714285,
            "central water tank storage": 0.0,
            "coal": 0.0,
            "coal CHP power plant": 6.0,
            "coal CHP power plant (CC)": 4.826467726000001,
            "coal power plant": 6.0,
            "coal power plant (CC)": 9.75,
            "decentral CHP": 0.0,
            "decentral air-sourced heat pump": 0.0,
            "decentral coal boiler": 0.0,
            "decentral gas boiler": 0.0,
            "decentral ground-sourced heat pump": 0.0,
            "decentral oil boiler": 0.0,
            "decentral resistive heater": 0.0,
            "decentral solar thermal": 0.0,
            "decentral solid biomass boiler": 0.0,
            "decentral water tank storage": 0.0,
            "direct air capture": 0.0,
            "electrolysis": 0.0,
            "electrolysis_soec": 0.0,
            "fischer tropsch import": 0.0,
            "fuel cell": 0.0,
            "gas": 0.0,
            "geothermal": 0.0,
            "helmeth": 0.0,
            "hydro": 0.0,
            "hydrogen": 0.0,
            "hydrogen pore storage underground": 0.0,
            "hydrogen storage": 0.0,
            "hydrogen storage underground": 0.0,
            "industry CC": 0.0,
            "lignite": 0.0,
            "lignite CHP power plant": 8.5,
            "lignite CHP power plant (CC)": 5.215698993999999,
            "lignite power plant": 8.5,
            "lignite power plant (CC)": 12.25,
            "methanation": 0.0,
            "micro CHP": 0.0,
            "nuclear": 2.1,
            "offwind": 0.01,
            "offwind-ac-connection-submarine": 0.0,
            "offwind-ac-connection-underground": 0.0,
            "offwind-ac-station": 0.0,
            "offwind-dc-connection-submarine": 0.0,
            "offwind-dc-connection-underground": 0.0,
            "offwind-dc-station": 0.0,
            "oil": 0.0,
            "oil power plant": 6.0,
            "onwind-1": 0.01,
            "onwind-2": 0.01,
            "onwind-3": 0.01,
            "onwind-4": 0.01,
            "onwind-landcosts": 0.0,
            "retrofitting I": 0.0,
            "retrofitting II": 0.0,
            "ror": 0.0,
            "solar": 0.01,
            "solar-rooftop": 0.0,
            "solar-utility": 0.0,
            "solid biomass": 0.0,
            "value of lost load": 0.0,
            "water tank charger": 0.0,
            "water tank discharger": 0.0,
        },
        "c_b": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": nan,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": 0.7,
            "central gas CHP CC": 0.7,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": 1.01,
            "central solid biomass CHP CC": 1.01,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": 1.01,
            "coal CHP power plant (CC)": 1.01,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": 1.01,
            "lignite CHP power plant (CC)": 1.01,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "c_v": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": nan,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": 0.17,
            "central gas CHP CC": 0.17,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": 0.15,
            "central solid biomass CHP CC": 0.15,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": 0.15,
            "coal CHP power plant (CC)": 0.15,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": 0.15,
            "lignite CHP power plant (CC)": 0.15,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "capture_rate": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": 0.9115763548571428,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "compression-electricity-input": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": 0.0916995074285714,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": 0.15,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "compression-heat-output": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": 0.1497044337142857,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": 0.2,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "discount rate": {
            "BEV vehicle": 0.07,
            "CCGT": 0.07,
            "CH4 pipeline": 0.07,
            "CH4 pipeline DIN1000": 0.07,
            "CH4 pipeline DIN1400": 0.07,
            "CH4 pipeline DIN600": 0.07,
            "CH4 pipeline retrofit": 0.07,
            "CH4 pipeline retrofit DIN1000": 0.07,
            "CH4 pipeline retrofit DIN1400": 0.07,
            "CH4 pipeline retrofit DIN600": 0.07,
            "CH4-powered PEMFC with internal SMR": 0.07,
            "CNG vehicle": 0.07,
            "DSM-household: charger": 0.07,
            "DSM-household: discharger": 0.07,
            "DSM-household: storage": 0.07,
            "DSM_Aluminum": 0.07,
            "DSM_Cement": 0.07,
            "DSM_Chlorine": 0.07,
            "DSM_ChlorineM": 0.07,
            "DSM_ChlorineQ": 0.07,
            "DSM_Household (turn-off)": 0.07,
            "DSM_Paper": 0.07,
            "DSM_PaperH": 0.07,
            "DSM_PaperS": 0.07,
            "DSM_Steel": 0.07,
            "FCEV vehicle": 0.07,
            "Fischer-Tropsch": 0.07,
            "H2 pipeline": 0.07,
            "H2-powered PEMFC with internal SMR": 0.07,
            "HEV vehicle": 0.07,
            "HVAC overhead": 0.07,
            "HVDC inverter pair": 0.07,
            "HVDC overhead": 0.07,
            "HVDC submarine": 0.07,
            "ICE vehicle": 0.07,
            "OCGT": 0.07,
            "PHEV vehicle": 0.07,
            "PHS": 0.07,
            "SMR": 0.07,
            "SMR CC": 0.07,
            "battery inverter": 0.07,
            "battery storage": 0.07,
            "biogas": 0.07,
            "biogas import": 0.07,
            "biomass": 0.07,
            "biomass CHP capture": 0.07,
            "central CHP": 0.07,
            "central air-sourced heat pump": 0.07,
            "central gas CHP": 0.07,
            "central gas CHP CC": 0.07,
            "central gas boiler": 0.07,
            "central resistive heater": 0.07,
            "central solar thermal": 0.07,
            "central solid biomass CHP": 0.07,
            "central solid biomass CHP CC": 0.07,
            "central water tank storage": 0.07,
            "coal": 0.07,
            "coal CHP power plant": 0.1,
            "coal CHP power plant (CC)": 0.1,
            "coal power plant": 0.1,
            "coal power plant (CC)": 0.1,
            "decentral CHP": 0.07,
            "decentral air-sourced heat pump": 0.07,
            "decentral coal boiler": 0.07,
            "decentral gas boiler": 0.07,
            "decentral ground-sourced heat pump": 0.07,
            "decentral oil boiler": 0.07,
            "decentral resistive heater": 0.07,
            "decentral solar thermal": 0.07,
            "decentral solid biomass boiler": 0.07,
            "decentral water tank storage": 0.07,
            "direct air capture": 0.07,
            "electrolysis": 0.07,
            "electrolysis_soec": 0.07,
            "fischer tropsch import": 0.07,
            "fuel cell": 0.07,
            "gas": 0.07,
            "geothermal": 0.07,
            "helmeth": 0.07,
            "hydro": 0.07,
            "hydrogen": 0.07,
            "hydrogen pore storage underground": 0.07,
            "hydrogen storage": 0.07,
            "hydrogen storage underground": 0.07,
            "industry CC": 0.07,
            "lignite": 0.07,
            "lignite CHP power plant": 0.1,
            "lignite CHP power plant (CC)": 0.1,
            "lignite power plant": 0.1,
            "lignite power plant (CC)": 0.1,
            "methanation": 0.07,
            "micro CHP": 0.07,
            "nuclear": 0.1,
            "offwind": 0.07,
            "offwind-ac-connection-submarine": 0.07,
            "offwind-ac-connection-underground": 0.07,
            "offwind-ac-station": 0.07,
            "offwind-dc-connection-submarine": 0.07,
            "offwind-dc-connection-underground": 0.07,
            "offwind-dc-station": 0.07,
            "oil": 0.07,
            "oil power plant": 0.07,
            "onwind-1": 0.07,
            "onwind-2": 0.07,
            "onwind-3": 0.07,
            "onwind-4": 0.07,
            "onwind-landcosts": 0.07,
            "retrofitting I": 0.07,
            "retrofitting II": 0.07,
            "ror": 0.07,
            "solar": 0.07,
            "solar-rooftop": 0.07,
            "solar-utility": 0.07,
            "solid biomass": 0.07,
            "value of lost load": 0.07,
            "water tank charger": 0.07,
            "water tank discharger": 0.07,
        },
        "efficiency": {
            "BEV vehicle": 1.0,
            "CCGT": 0.61,
            "CH4 pipeline": 1.0,
            "CH4 pipeline DIN1000": 1.0,
            "CH4 pipeline DIN1400": 1.0,
            "CH4 pipeline DIN600": 1.0,
            "CH4 pipeline retrofit": 1.0,
            "CH4 pipeline retrofit DIN1000": 1.0,
            "CH4 pipeline retrofit DIN1400": 1.0,
            "CH4 pipeline retrofit DIN600": 1.0,
            "CH4-powered PEMFC with internal SMR": 0.3458571428571428,
            "CNG vehicle": 0.3139567945714286,
            "DSM-household: charger": 1.0,
            "DSM-household: discharger": 1.0,
            "DSM-household: storage": 1.0,
            "DSM_Aluminum": 1.0,
            "DSM_Cement": 1.0,
            "DSM_Chlorine": 1.0,
            "DSM_ChlorineM": 1.0,
            "DSM_ChlorineQ": 1.0,
            "DSM_Household (turn-off)": 1.0,
            "DSM_Paper": 1.0,
            "DSM_PaperH": 1.0,
            "DSM_PaperS": 1.0,
            "DSM_Steel": 1.0,
            "FCEV vehicle": 0.590288061,
            "Fischer-Tropsch": 0.6071428571428572,
            "H2 pipeline": 1.0,
            "H2-powered PEMFC with internal SMR": 0.4714285714285714,
            "HEV vehicle": 0.4319483494285714,
            "HVAC overhead": 1.0,
            "HVDC inverter pair": 1.0,
            "HVDC overhead": 1.0,
            "HVDC submarine": 1.0,
            "ICE vehicle": 0.3047020382857143,
            "OCGT": 0.4,
            "PHEV vehicle": 0.613950273,
            "PHS": 0.75,
            "SMR": 0.74,
            "SMR CC": 0.67,
            "battery inverter": 0.9142857142857144,
            "battery storage": 1.0,
            "biogas": 1.0,
            "biogas import": 1.0,
            "biomass": 0.3,
            "biomass CHP capture": 1.0,
            "central CHP": 1.0,
            "central air-sourced heat pump": 3.8428571428571434,
            "central gas CHP": 0.45,
            "central gas CHP CC": 0.36,
            "central gas boiler": 0.8442857142857143,
            "central resistive heater": 0.98,
            "central solar thermal": 1.0,
            "central solid biomass CHP": 0.2921428571428571,
            "central solid biomass CHP CC": 0.3234285714285714,
            "central water tank storage": 1.0,
            "coal": 1.0,
            "coal CHP power plant": 0.38,
            "coal CHP power plant (CC)": 0.341652142,
            "coal power plant": 0.46,
            "coal power plant (CC)": 0.41,
            "decentral CHP": 1.0,
            "decentral air-sourced heat pump": 3.8928571428571432,
            "decentral coal boiler": 0.77,
            "decentral gas boiler": 0.83,
            "decentral ground-sourced heat pump": 4.221428571428572,
            "decentral oil boiler": 0.87,
            "decentral resistive heater": 0.99,
            "decentral solar thermal": 1.0,
            "decentral solid biomass boiler": 0.77,
            "decentral water tank storage": 1.0,
            "direct air capture": 1.0,
            "electrolysis": 0.7928571428571429,
            "electrolysis_soec": 1.0,
            "fischer tropsch import": 1.0,
            "fuel cell": 0.4928571428571429,
            "gas": 1.0,
            "geothermal": 0.2437142857142857,
            "helmeth": 0.7928571428571428,
            "hydro": 0.9,
            "hydrogen": 1.0,
            "hydrogen pore storage underground": 1.0,
            "hydrogen storage": 1.0,
            "hydrogen storage underground": 1.0,
            "industry CC": 0.9,
            "lignite": 1.0,
            "lignite CHP power plant": 0.38,
            "lignite CHP power plant (CC)": 0.341652142,
            "lignite power plant": 0.46,
            "lignite power plant (CC)": 0.4,
            "methanation": 0.77,
            "micro CHP": 0.36,
            "nuclear": 0.33,
            "offwind": 1.0,
            "offwind-ac-connection-submarine": 1.0,
            "offwind-ac-connection-underground": 1.0,
            "offwind-ac-station": 1.0,
            "offwind-dc-connection-submarine": 1.0,
            "offwind-dc-connection-underground": 1.0,
            "offwind-dc-station": 1.0,
            "oil": 1.0,
            "oil power plant": 0.37,
            "onwind-1": 1.0,
            "onwind-2": 1.0,
            "onwind-3": 1.0,
            "onwind-4": 1.0,
            "onwind-landcosts": 1.0,
            "retrofitting I": 1.0,
            "retrofitting II": 1.0,
            "ror": 0.9,
            "solar": 1.0,
            "solar-rooftop": 1.0,
            "solar-utility": 1.0,
            "solid biomass": 1.0,
            "value of lost load": 1.0,
            "water tank charger": 0.9,
            "water tank discharger": 0.9,
        },
        "efficiency-heat": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": 0.597,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": 0.4914285714285714,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": nan,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": 0.45,
            "central gas CHP CC": 0.36,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": 0.483,
            "central solid biomass CHP CC": 0.483,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": 0.15,
            "coal CHP power plant (CC)": 0.135283019,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": 0.15,
            "lignite CHP power plant (CC)": 0.135283019,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": 0.44,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "electric_efficiency": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": nan,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": 0.8257142857142857,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "electricity-input": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": 0.0268177337142857,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": 0.32,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "fuel": {
            "BEV vehicle": 0.0,
            "CCGT": 0.0,
            "CH4 pipeline": 0.0,
            "CH4 pipeline DIN1000": 0.0,
            "CH4 pipeline DIN1400": 0.0,
            "CH4 pipeline DIN600": 0.0,
            "CH4 pipeline retrofit": 0.0,
            "CH4 pipeline retrofit DIN1000": 0.0,
            "CH4 pipeline retrofit DIN1400": 0.0,
            "CH4 pipeline retrofit DIN600": 0.0,
            "CH4-powered PEMFC with internal SMR": 0.0,
            "CNG vehicle": 0.0,
            "DSM-household: charger": 0.0,
            "DSM-household: discharger": 0.0,
            "DSM-household: storage": 0.0,
            "DSM_Aluminum": 0.0,
            "DSM_Cement": 0.0,
            "DSM_Chlorine": 0.0,
            "DSM_ChlorineM": 0.0,
            "DSM_ChlorineQ": 0.0,
            "DSM_Household (turn-off)": 0.0,
            "DSM_Paper": 0.0,
            "DSM_PaperH": 0.0,
            "DSM_PaperS": 0.0,
            "DSM_Steel": 0.0,
            "FCEV vehicle": 0.0,
            "Fischer-Tropsch": 0.0,
            "H2 pipeline": 0.0,
            "H2-powered PEMFC with internal SMR": 0.0,
            "HEV vehicle": 0.0,
            "HVAC overhead": 0.0,
            "HVDC inverter pair": 0.0,
            "HVDC overhead": 0.0,
            "HVDC submarine": 0.0,
            "ICE vehicle": 0.0,
            "OCGT": 0.0,
            "PHEV vehicle": 0.0,
            "PHS": 0.0,
            "SMR": 0.0,
            "SMR CC": 0.0,
            "battery inverter": 0.0,
            "battery storage": 0.0,
            "biogas": 60.0,
            "biogas import": 116.85714285714286,
            "biomass": 0.0,
            "biomass CHP capture": 0.0,
            "central CHP": 0.0,
            "central air-sourced heat pump": 0.0,
            "central gas CHP": 0.0,
            "central gas CHP CC": 0.0,
            "central gas boiler": 0.0,
            "central resistive heater": 0.0,
            "central solar thermal": 0.0,
            "central solid biomass CHP": 0.0,
            "central solid biomass CHP CC": 0.0,
            "central water tank storage": 0.0,
            "coal": 8.168662582285714,
            "coal CHP power plant": 0.0,
            "coal CHP power plant (CC)": 0.0,
            "coal power plant": 0.0,
            "coal power plant (CC)": 0.0,
            "decentral CHP": 0.0,
            "decentral air-sourced heat pump": 0.0,
            "decentral coal boiler": 0.0,
            "decentral gas boiler": 0.0,
            "decentral ground-sourced heat pump": 0.0,
            "decentral oil boiler": 0.0,
            "decentral resistive heater": 0.0,
            "decentral solar thermal": 0.0,
            "decentral solid biomass boiler": 0.0,
            "decentral water tank storage": 0.0,
            "direct air capture": 0.0,
            "electrolysis": 0.0,
            "electrolysis_soec": 0.0,
            "fischer tropsch import": 124.0,
            "fuel cell": 0.0,
            "gas": 21.799515062857143,
            "geothermal": 0.0,
            "helmeth": 0.0,
            "hydro": 0.0,
            "hydrogen": 101.0,
            "hydrogen pore storage underground": 0.0,
            "hydrogen storage": 0.0,
            "hydrogen storage underground": 0.0,
            "industry CC": 0.0,
            "lignite": 2.535986628,
            "lignite CHP power plant": 0.0,
            "lignite CHP power plant (CC)": 0.0,
            "lignite power plant": 0.0,
            "lignite power plant (CC)": 0.0,
            "methanation": 0.0,
            "micro CHP": 0.0,
            "nuclear": 3.146303094,
            "offwind": 0.0,
            "offwind-ac-connection-submarine": 0.0,
            "offwind-ac-connection-underground": 0.0,
            "offwind-ac-station": 0.0,
            "offwind-dc-connection-submarine": 0.0,
            "offwind-dc-connection-underground": 0.0,
            "offwind-dc-station": 0.0,
            "oil": 34.49733088142857,
            "oil power plant": 0.0,
            "onwind-1": 0.0,
            "onwind-2": 0.0,
            "onwind-3": 0.0,
            "onwind-4": 0.0,
            "onwind-landcosts": 0.0,
            "retrofitting I": 0.0,
            "retrofitting II": 0.0,
            "ror": 0.0,
            "solar": 0.0,
            "solar-rooftop": 0.0,
            "solar-utility": 0.0,
            "solid biomass": 75.0,
            "value of lost load": 1000.0,
            "water tank charger": 0.0,
            "water tank discharger": 0.0,
        },
        "heat-input": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": 0.7740197045714285,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": 2.071428571428572,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "heat-output": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": 0.7740197045714285,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": 1.0357142857142858,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "heat_efficiency": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": nan,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": nan,
            "central gas CHP CC": nan,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": nan,
            "central solid biomass CHP CC": nan,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": nan,
            "coal CHP power plant (CC)": nan,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": 0.8257142857142857,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": nan,
            "lignite CHP power plant (CC)": nan,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "investment": {
            "BEV vehicle": 105972669.43714283,
            "CCGT": 747310.0,
            "CH4 pipeline": 0.0,
            "CH4 pipeline DIN1000": 362.22,
            "CH4 pipeline DIN1400": 207.23,
            "CH4 pipeline DIN600": 852.08,
            "CH4 pipeline retrofit": 0.0,
            "CH4 pipeline retrofit DIN1000": 140.05,
            "CH4 pipeline retrofit DIN1400": 106.39,
            "CH4 pipeline retrofit DIN600": 291.68,
            "CH4-powered PEMFC with internal SMR": 9003148.368714286,
            "CNG vehicle": 30253636.032857142,
            "DSM-household: charger": 0.0,
            "DSM-household: discharger": 0.0,
            "DSM-household: storage": 13333.333331428572,
            "DSM_Aluminum": 0.0,
            "DSM_Cement": 0.0,
            "DSM_Chlorine": 0.0,
            "DSM_ChlorineM": 0.0,
            "DSM_ChlorineQ": 0.0,
            "DSM_Household (turn-off)": 0.0,
            "DSM_Paper": 0.0,
            "DSM_PaperH": 0.0,
            "DSM_PaperS": 0.0,
            "DSM_Steel": 0.0,
            "FCEV vehicle": 102925297.08857144,
            "Fischer-Tropsch": 548476.4285714285,
            "H2 pipeline": 658.57,
            "H2-powered PEMFC with internal SMR": 12244897.959,
            "HEV vehicle": 47023996.71428572,
            "HVAC overhead": 500.0,
            "HVDC inverter pair": 600000.0,
            "HVDC overhead": 2000.0,
            "HVDC submarine": 2000.0,
            "ICE vehicle": 28052254.285714284,
            "OCGT": 392310.7101000001,
            "PHEV vehicle": 60813013.28571428,
            "PHS": 1218000.0,
            "SMR": 416740.0,
            "SMR CC": 596500.0,
            "battery inverter": 191428.57142857145,
            "battery storage": 271000.0,
            "biogas": 0.0,
            "biogas import": 0.0,
            "biomass": 3230351.7720000003,
            "biomass CHP capture": 2323152.709142857,
            "central CHP": 1000392.311,
            "central air-sourced heat pump": 827372.4436,
            "central gas CHP": 964680.0,
            "central gas CHP CC": 2806528.5714285714,
            "central gas boiler": 55714.285714285725,
            "central resistive heater": 112244.898,
            "central solar thermal": 224642.58285714284,
            "central solid biomass CHP": 2759428.571428572,
            "central solid biomass CHP CC": 5805357.142857143,
            "central water tank storage": 175.0,
            "coal": 0.0,
            "coal CHP power plant": 2824637.1130000004,
            "coal CHP power plant (CC)": 5870351.398714287,
            "coal power plant": 1480800.0,
            "coal power plant (CC)": 4515114.2857142845,
            "decentral CHP": 1627007.984285714,
            "decentral air-sourced heat pump": 897142.8571428572,
            "decentral coal boiler": 196448.21428571426,
            "decentral gas boiler": 178600.0,
            "decentral ground-sourced heat pump": 1428571.4285714284,
            "decentral oil boiler": 226465.7142857143,
            "decentral resistive heater": 623196.2482142857,
            "decentral solar thermal": 495238.0952571428,
            "decentral solid biomass boiler": 398095.2381142857,
            "decentral water tank storage": 13666.66667,
            "direct air capture": 5714285.714285715,
            "electrolysis": 999860.0,
            "electrolysis_soec": 1445464.2857142857,
            "fischer tropsch import": 0.0,
            "fuel cell": 1428571.4285714289,
            "gas": 0.0,
            "geothermal": 9818883.222714284,
            "helmeth": 3121205.1947142854,
            "hydro": 2336000.0,
            "hydrogen": 0.0,
            "hydrogen pore storage underground": 210.0,
            "hydrogen storage": 41571.42857142857,
            "hydrogen storage underground": 210.0,
            "industry CC": 311.38130275714286,
            "lignite": 0.0,
            "lignite CHP power plant": 2451941.938,
            "lignite CHP power plant (CC)": 6606946.621142857,
            "lignite power plant": 1632518.5714285716,
            "lignite power plant (CC)": 5473408.571428571,
            "methanation": 553421.4285714285,
            "micro CHP": 6726710.203999999,
            "nuclear": 3259121.224,
            "offwind": 2113573.950285714,
            "offwind-ac-connection-submarine": 2685.0,
            "offwind-ac-connection-underground": 1342.0,
            "offwind-ac-station": 250000.0,
            "offwind-dc-connection-submarine": 2000.0,
            "offwind-dc-connection-underground": 1000.0,
            "offwind-dc-station": 400000.0,
            "oil": 0.0,
            "oil power plant": 344000.0,
            "onwind-1": 1330714.2857142857,
            "onwind-2": 1330714.2857142857,
            "onwind-3": 1330714.2857142857,
            "onwind-4": 1330714.2857142857,
            "onwind-landcosts": 0.0,
            "retrofitting I": 102.6673377,
            "retrofitting II": 181.17765480000003,
            "ror": 5000000.0,
            "solar": 600000.0,
            "solar-rooftop": 784948.3457142857,
            "solar-utility": 667908.9838857142,
            "solid biomass": 0.0,
            "value of lost load": 0.0,
            "water tank charger": 0.0,
            "water tank discharger": 0.0,
        },
        "lifetime": {
            "BEV vehicle": 13.0,
            "CCGT": 33.0,
            "CH4 pipeline": 55.0,
            "CH4 pipeline DIN1000": 25.0,
            "CH4 pipeline DIN1400": 25.0,
            "CH4 pipeline DIN600": 25.0,
            "CH4 pipeline retrofit": 55.0,
            "CH4 pipeline retrofit DIN1000": 25.0,
            "CH4 pipeline retrofit DIN1400": 25.0,
            "CH4 pipeline retrofit DIN600": 25.0,
            "CH4-powered PEMFC with internal SMR": 20.0,
            "CNG vehicle": 13.0,
            "DSM-household: charger": 25.0,
            "DSM-household: discharger": 25.0,
            "DSM-household: storage": 10.0,
            "DSM_Aluminum": 25.0,
            "DSM_Cement": 25.0,
            "DSM_Chlorine": 25.0,
            "DSM_ChlorineM": 25.0,
            "DSM_ChlorineQ": 25.0,
            "DSM_Household (turn-off)": 25.0,
            "DSM_Paper": 25.0,
            "DSM_PaperH": 25.0,
            "DSM_PaperS": 25.0,
            "DSM_Steel": 25.0,
            "FCEV vehicle": 13.0,
            "Fischer-Tropsch": 25.142857142857142,
            "H2 pipeline": 55.0,
            "H2-powered PEMFC with internal SMR": 12.571428571428571,
            "HEV vehicle": 13.0,
            "HVAC overhead": 40.0,
            "HVDC inverter pair": 40.0,
            "HVDC overhead": 40.0,
            "HVDC submarine": 40.0,
            "ICE vehicle": 13.0,
            "OCGT": 25.0,
            "PHEV vehicle": 13.0,
            "PHS": 45.0,
            "SMR": 25.0,
            "SMR CC": 25.0,
            "battery inverter": 21.428571428571427,
            "battery storage": 21.428571428571427,
            "biogas": 25.0,
            "biogas import": 25.0,
            "biomass": 30.0,
            "biomass CHP capture": 25.0,
            "central CHP": 30.0,
            "central air-sourced heat pump": 25.0,
            "central gas CHP": 30.0,
            "central gas CHP CC": 30.0,
            "central gas boiler": 25.0,
            "central resistive heater": 20.0,
            "central solar thermal": 30.0,
            "central solid biomass CHP": 25.0,
            "central solid biomass CHP CC": 25.0,
            "central water tank storage": 40.0,
            "coal": 25.0,
            "coal CHP power plant": 45.0,
            "coal CHP power plant (CC)": 40.0,
            "coal power plant": 45.0,
            "coal power plant (CC)": 45.0,
            "decentral CHP": 20.0,
            "decentral air-sourced heat pump": 25.0,
            "decentral coal boiler": 25.0,
            "decentral gas boiler": 25.0,
            "decentral ground-sourced heat pump": 25.0,
            "decentral oil boiler": 25.0,
            "decentral resistive heater": 25.0,
            "decentral solar thermal": 25.0,
            "decentral solid biomass boiler": 25.0,
            "decentral water tank storage": 30.0,
            "direct air capture": 20.0,
            "electrolysis": 21.285714285714285,
            "electrolysis_soec": 19.285714285714285,
            "fischer tropsch import": 25.0,
            "fuel cell": 10.0,
            "gas": 25.0,
            "geothermal": 30.0,
            "helmeth": 25.0,
            "hydro": 75.0,
            "hydrogen": 25.0,
            "hydrogen pore storage underground": 50.0,
            "hydrogen storage": 30.0,
            "hydrogen storage underground": 50.0,
            "industry CC": 25.0,
            "lignite": 25.0,
            "lignite CHP power plant": 45.0,
            "lignite CHP power plant (CC)": 40.0,
            "lignite power plant": 45.0,
            "lignite power plant (CC)": 45.0,
            "methanation": 25.142857142857142,
            "micro CHP": 22.142857142857142,
            "nuclear": 60.0,
            "offwind": 25.0,
            "offwind-ac-connection-submarine": 25.0,
            "offwind-ac-connection-underground": 25.0,
            "offwind-ac-station": 25.0,
            "offwind-dc-connection-submarine": 25.0,
            "offwind-dc-connection-underground": 25.0,
            "offwind-dc-station": 25.0,
            "oil": 25.0,
            "oil power plant": 25.0,
            "onwind-1": 25.428571428571427,
            "onwind-2": 25.428571428571427,
            "onwind-3": 25.428571428571427,
            "onwind-4": 25.428571428571427,
            "onwind-landcosts": 25.0,
            "retrofitting I": 50.0,
            "retrofitting II": 50.0,
            "ror": 100.0,
            "solar": 25.0,
            "solar-rooftop": 25.0,
            "solar-utility": 25.0,
            "solid biomass": 25.0,
            "value of lost load": 25.0,
            "water tank charger": 25.0,
            "water tank discharger": 25.0,
        },
        "p_nom_ratio": {
            "BEV vehicle": nan,
            "CCGT": nan,
            "CH4 pipeline": nan,
            "CH4 pipeline DIN1000": nan,
            "CH4 pipeline DIN1400": nan,
            "CH4 pipeline DIN600": nan,
            "CH4 pipeline retrofit": nan,
            "CH4 pipeline retrofit DIN1000": nan,
            "CH4 pipeline retrofit DIN1400": nan,
            "CH4 pipeline retrofit DIN600": nan,
            "CH4-powered PEMFC with internal SMR": nan,
            "CNG vehicle": nan,
            "DSM-household: charger": nan,
            "DSM-household: discharger": nan,
            "DSM-household: storage": nan,
            "DSM_Aluminum": nan,
            "DSM_Cement": nan,
            "DSM_Chlorine": nan,
            "DSM_ChlorineM": nan,
            "DSM_ChlorineQ": nan,
            "DSM_Household (turn-off)": nan,
            "DSM_Paper": nan,
            "DSM_PaperH": nan,
            "DSM_PaperS": nan,
            "DSM_Steel": nan,
            "FCEV vehicle": nan,
            "Fischer-Tropsch": nan,
            "H2 pipeline": nan,
            "H2-powered PEMFC with internal SMR": nan,
            "HEV vehicle": nan,
            "HVAC overhead": nan,
            "HVDC inverter pair": nan,
            "HVDC overhead": nan,
            "HVDC submarine": nan,
            "ICE vehicle": nan,
            "OCGT": nan,
            "PHEV vehicle": nan,
            "PHS": nan,
            "SMR": nan,
            "SMR CC": nan,
            "battery inverter": nan,
            "battery storage": nan,
            "biogas": nan,
            "biogas import": nan,
            "biomass": nan,
            "biomass CHP capture": nan,
            "central CHP": nan,
            "central air-sourced heat pump": nan,
            "central gas CHP": 1.046511628,
            "central gas CHP CC": 1.046511628,
            "central gas boiler": nan,
            "central resistive heater": nan,
            "central solar thermal": nan,
            "central solid biomass CHP": 1.0,
            "central solid biomass CHP CC": 1.0,
            "central water tank storage": nan,
            "coal": nan,
            "coal CHP power plant": 1.0,
            "coal CHP power plant (CC)": 1.0,
            "coal power plant": nan,
            "coal power plant (CC)": nan,
            "decentral CHP": nan,
            "decentral air-sourced heat pump": nan,
            "decentral coal boiler": nan,
            "decentral gas boiler": nan,
            "decentral ground-sourced heat pump": nan,
            "decentral oil boiler": nan,
            "decentral resistive heater": nan,
            "decentral solar thermal": nan,
            "decentral solid biomass boiler": nan,
            "decentral water tank storage": nan,
            "direct air capture": nan,
            "electrolysis": nan,
            "electrolysis_soec": nan,
            "fischer tropsch import": nan,
            "fuel cell": nan,
            "gas": nan,
            "geothermal": nan,
            "helmeth": nan,
            "hydro": nan,
            "hydrogen": nan,
            "hydrogen pore storage underground": nan,
            "hydrogen storage": nan,
            "hydrogen storage underground": nan,
            "industry CC": nan,
            "lignite": nan,
            "lignite CHP power plant": 1.0,
            "lignite CHP power plant (CC)": 1.0,
            "lignite power plant": nan,
            "lignite power plant (CC)": nan,
            "methanation": nan,
            "micro CHP": nan,
            "nuclear": nan,
            "offwind": nan,
            "offwind-ac-connection-submarine": nan,
            "offwind-ac-connection-underground": nan,
            "offwind-ac-station": nan,
            "offwind-dc-connection-submarine": nan,
            "offwind-dc-connection-underground": nan,
            "offwind-dc-station": nan,
            "oil": nan,
            "oil power plant": nan,
            "onwind-1": nan,
            "onwind-2": nan,
            "onwind-3": nan,
            "onwind-4": nan,
            "onwind-landcosts": nan,
            "retrofitting I": nan,
            "retrofitting II": nan,
            "ror": nan,
            "solar": nan,
            "solar-rooftop": nan,
            "solar-utility": nan,
            "solid biomass": nan,
            "value of lost load": nan,
            "water tank charger": nan,
            "water tank discharger": nan,
        },
        "fixed": {
            "BEV vehicle": 316992994.43543404,
            "CCGT": 1939420.1313143142,
            "CH4 pipeline": 0.0,
            "CH4 pipeline DIN1000": 777.0571386917374,
            "CH4 pipeline DIN1400": 444.5628370909633,
            "CH4 pipeline DIN600": 1827.9356378346185,
            "CH4 pipeline retrofit": 0.0,
            "CH4 pipeline retrofit DIN1000": 300.44407341885545,
            "CH4 pipeline retrofit DIN1400": 228.23452317766532,
            "CH4 pipeline retrofit DIN600": 625.7302915730936,
            "CH4-powered PEMFC with internal SMR": 34957433.01951048,
            "CNG vehicle": 90496830.26342472,
            "DSM-household: charger": 0.0,
            "DSM-household: discharger": 0.0,
            "DSM-household: storage": 47459.167569008336,
            "DSM_Aluminum": 0.0,
            "DSM_Cement": 0.0,
            "DSM_Chlorine": 0.0,
            "DSM_ChlorineM": 0.0,
            "DSM_ChlorineQ": 0.0,
            "DSM_Household (turn-off)": 0.0,
            "DSM_Paper": 0.0,
            "DSM_PaperH": 0.0,
            "DSM_PaperS": 0.0,
            "DSM_Steel": 0.0,
            "FCEV vehicle": 307877477.2830954,
            "Fischer-Tropsch": 1781317.1357237825,
            "H2 pipeline": 1312.7988113627976,
            "H2-powered PEMFC with internal SMR": 56620541.78119119,
            "HEV vehicle": 140661527.24052146,
            "HVAC overhead": 1187.614235920129,
            "HVDC inverter pair": 1425137.0831041548,
            "HVDC overhead": 4750.456943680516,
            "HVDC submarine": 4750.456943680516,
            "ICE vehicle": 83911900.43549176,
            "OCGT": 1160362.0755784402,
            "PHEV vehicle": 181908215.4339298,
            "PHS": 2538061.936829858,
            "SMR": 1626486.4858385946,
            "SMR CC": 2048069.6469713722,
            "battery inverter": 461340.805337028,
            "battery storage": 619619.238173816,
            "biogas": 0.0,
            "biogas import": 0.0,
            "biomass": 10553773.278399581,
            "biomass CHP capture": 6726137.9207106205,
            "central CHP": 3020746.6559345997,
            "central air-sourced heat pump": 1825164.7599189065,
            "central gas CHP": 2698364.3934774683,
            "central gas CHP CC": 7733430.684603044,
            "central gas boiler": 166321.79184306998,
            "central resistive heater": 292378.1080593327,
            "central solar thermal": 486836.4395833128,
            "central solid biomass CHP": 7879934.400182774,
            "central solid biomass CHP CC": 16400261.841833973,
            "central water tank storage": 341.2899825720451,
            "coal": 0.0,
            "coal CHP power plant": 8458416.03161026,
            "coal CHP power plant (CC)": 19045459.67678399,
            "coal power plant": 4530914.9936380405,
            "coal power plant (CC)": 14510811.708975274,
            "decentral CHP": 5069129.910594032,
            "decentral air-sourced heat pump": 2594260.376030847,
            "decentral coal boiler": 570750.7957062591,
            "decentral gas boiler": 592891.3842728434,
            "decentral ground-sourced heat pump": 3741957.2476768326,
            "decentral oil boiler": 676429.4418179413,
            "decentral resistive heater": 1344270.2375879618,
            "decentral solar thermal": 1345279.072824608,
            "decentral solid biomass boiler": 1323480.214026275,
            "decentral water tank storage": 41422.576763899386,
            "direct air capture": 20556132.24903653,
            "electrolysis": 2794282.3516328316,
            "electrolysis_soec": 4855476.84388241,
            "fischer tropsch import": 0.0,
            "fuel cell": 6870625.097405882,
            "gas": 0.0,
            "geothermal": 29146798.278232962,
            "helmeth": 9036709.698792242,
            "hydro": 4413730.755371954,
            "hydrogen": 0.0,
            "hydrogen pore storage underground": 380.41421008294026,
            "hydrogen storage": 83752.29793476197,
            "hydrogen storage underground": 380.41421008294026,
            "industry CC": 823.6854174394489,
            "lignite": 0.0,
            "lignite CHP power plant": 7633957.960206666,
            "lignite CHP power plant (CC)": 19901002.443025615,
            "lignite power plant": 5052277.8533213595,
            "lignite power plant (CC)": 17726032.265076686,
            "methanation": 1604470.3868927888,
            "micro CHP": 24578218.762564607,
            "nuclear": 10675631.599490786,
            "offwind": 7073775.598458461,
            "offwind-ac-connection-submarine": 5760.030968437179,
            "offwind-ac-connection-underground": 2878.942852753331,
            "offwind-ac-station": 536315.73262916,
            "offwind-dc-connection-submarine": 4290.52586103328,
            "offwind-dc-connection-underground": 2145.26293051664,
            "offwind-dc-station": 858105.1722066561,
            "oil": 0.0,
            "oil power plant": 950337.79504001,
            "onwind-1": 3628672.4606899726,
            "onwind-2": 3628672.4606899726,
            "onwind-3": 3628672.4606899726,
            "onwind-4": 3628672.4606899726,
            "onwind-landcosts": 0.0,
            "retrofitting I": 211.64833048435224,
            "retrofitting II": 373.49705386867447,
            "ror": 9060095.575382503,
            "solar": 1912157.8083099837,
            "solar-rooftop": 2127129.0129768173,
            "solar-utility": 1819768.0682571267,
            "solid biomass": 0.0,
            "value of lost load": 0.0,
            "water tank charger": 0.0,
            "water tank discharger": 0.0,
        },
    }
    return pd.DataFrame(expected_data)


@pytest.fixture(scope="module")
def expected_co2_share() -> pd.DataFrame:
    """
    Return the expected CO2 share.

    Expected data is the result from
    Toolbox.Import_functions.get_historical_emissions_per_country
    with all countries, except 'EU28', all options 'THI' and for
    cluster AT10.
    """
    expected_data = {
        "electricity": {
            "NO": 0.007554612594072355,
            "CH": 0.009263012800390899,
            "BA": 0.004616271350009418,
            "RS": 0.008070030900771885,
            "AL": 0.0008989582704655664,
            "ME": 0.000509592572557623,
            "MK": 0.0016722791676916139,
            "FR": 0.08553429072403554,
            "DE": 0.2334729845980263,
            "GB": 0.12936769012653634,
            "IT": 0.09738494374574719,
            "ES": 0.04862866406876785,
            "PL": 0.0839751065477898,
            "SE": 0.012160194267656712,
            "NL": 0.027784485024652417,
            "BE": 0.023442976689618852,
            "FI": 0.01210541014340007,
            "CZ": 0.036756937759289225,
            "DK": 0.010914575064984751,
            "PT": 0.00964540728294732,
            "RO": 0.038538898030006416,
            "BG": 0.0175903779301616,
            "EE": 0.008238046001386896,
            "GR": 0.016463711384387653,
            "LV": 0.004007342658643644,
            "HU": 0.01640009002046733,
            "IE": 0.007076177344366184,
            "SK": 0.0139744906409257,
            "LT": 0.007928551728504044,
            "HR": 0.005149010504087956,
            "LU": 0.00258617956325792,
            "SI": 0.0033887268687538695,
            "CY": 0.0008466457303740906,
            "MT": 0.00019264541023779257,
            "AT0 0": 0.00045938826262773984,
            "AT0 1": 0.00262243966784562,
            "AT0 2": 0.0029818260933870927,
            "AT0 3": 0.00087140422819407,
            "AT0 4": 0.0019340433703506491,
            "AT0 5": 0.002323497644667211,
            "AT0 6": 0.0008685471651508107,
            "AT0 7": 0.0006198900689433906,
            "AT0 8": 7.534558368106797e-05,
            "AT0 9": 0.0011043004001795992,
        }
    }
    return pd.DataFrame(expected_data)


@pytest.fixture
def with_mock_config(monkeypatch):
    """"""
    import esmtools.fileio as mod

    def mock(*args, **kwargs):
        s = """version: 0.1
tutorial: false
logging:
  level: INFO
workdir_root: /opt/esm/runs/manual__2023-11-13T16:08:29.443991+01:00/
repo_root: /opt/esm/runs/manual__2023-11-13T16:08:29.443991+01:00/
pypsa_eur_scripts_dir: ./pypsa-eur/scripts/
results_dir: results/
summary_dir: results
run: esm_run
foresight: myopic
scenario:
  simpl:
  - ''
  ll:
  - ''
  lv:
  - 1.0
  clusters:
  - AT10
  type: shape
  opts:
  - ''
  sector_opts:
  - Co2L0p79-120H-T-H-B-I
  planning_horizons:
  - '2015'
  - '2020'
  - '2030'
  - '2040'
  - '2050'
  name: szenario_AT10
countries:
- AL
- AT
- BA
- BE
- BG
- CH
- CZ
- DE
- DK
- EE
- ES
- FI
- FR
- GB
- GR
- HR
- HU
- IE
- IT
- LT
- LU
- LV
- ME
- MK
- NL
- 'NO'
- PL
- PT
- RO
- RS
- SE
- SI
- SK
constraints:
  co2_limit_per_country: true
  capacity_limit_per_country_carrier: true
  capacity_limit_carrier_global: true
  dispatch_limit_global: true
  dispatch_limit_country: true
  renewable_feedstock_industry: false
  fischer-tropsch_import_limit: false
  ft_import_share: 0.0
  h2_import_limit: false
  h2_import_share: 0.0
  h2_country_production: true
  h2_country_production_share: 0.5
  optimize_transmission_lines: true
  capacity_limit_ac_dc: true
  biogas_import_limit: false
  gas_import_share: 0.0
  german_renewable_share: true
  e_mobility: true
  bev_constraint_years:
  - '2030'
  - '2040'
  ice_constraint_years:
  - '2015'
  - '2020'
  cng_constraint_years:
  - '2015'
  - '2020'
  - '2025'
  - '2030'
  - '2035'
  - '2040'
  - '2045'
  - '2050'
  co2_transport_limit: true
  co2_transport_limit_years:
  - '2025'
  - '2030'
  - '2035'
  - '2040'
  - '2050'
  h2_import_limit_for_industry: false
  h2_import_limit_share_per_cluster: 0.0
  ft_sabatier_limit: false
  national_renewable_electricity:
    activate: true
    countries:
      AT:
      - '2030'
      - '2040'
  PV_rooftop_utility_ratio:
    activate: true
    non_countries:
    - AT
    - DE
    - IT1 0
    - ES
    - NL
    - EE
    years:
    - '2030'
    share_rooftop: 0.5
  decentral_heating_dispatch_profiles:
    profiles_and_dismantling: true
    allowed_dismantled_technologies:
    - gas boiler
    - oil boiler
    - coal boiler
    - biomass boiler
    - resistive heater
    allowed_dismantle_ratio_per_year: 0.05
snapshots:
  start: '2015-01-01'
  end: '2016-01-01'
  closed: left
enable:
  prepare_links_p_nom: false
  retrieve_databundle: false
  build_cutout: true
  build_natura_raster: true
  powerplantmatching: false
  build_renewable_profiles: true
customize_data:
  hydro_cap: true
  biomass:
    scenario: Med
    custom_mapping: true
    regional_biomass_share: true
    custom_data_source: false
  powerplants: true
  load: true
  use_hotmaps_data: true
  nuts_region_mapping_for_hotmaps: true
  heat_profile: true
  consider_national_statistics:
    heat_volumes: true
    total_road_demand: true
    electricity_rail: true
  industry_demand_per_cluster: true
  delete_p_nom_max:
    active: true
    countries:
      AT:
      - solar-utility
      - solar-rooftop
      - onwind-1
      - onwind-2
      - onwind-3
      - onwind-4
      DE:
      - solar-rooftop
  update_hydro_data: true
  custom_industry_input: true
  gas_generator_update: true
electricity:
  voltages:
  - 220.0
  - 300.0
  - 380.0
  co2limit: 77500000.0
  co2base: 3100000000.0
  extendable_carriers:
    Generator: []
    StorageUnit: []
    Store: []
  max_hours:
    battery: 6
    H2: 168
  powerplants_filter: false
  custom_powerplants: false
  conventional_carriers:
  - nuclear
  - oil
  - coal
  - lignite
  - OCGT
  - CCGT
  - ror
customize_network:
  switch: true
atlite:
  nprocesses: 1
  cutouts:
    europe-2015-era5:
      module: era5
      xs:
      - -12.0
      - 35.0
      ys:
      - 72.0
      - 33.0
      years:
      - 2015
      - 2015
      dx: 0.1
      dy: 0.1
  cutout_dir: pypsa-eur/cutouts
  cutout_name: europe-2015-era5
cutout_tech: onwind-1
renewable:
  onwind-1:
    cutout: europe-2015-era5
    resource:
      method: wind
      turbine: Suzlon_S82_1.5_MW
    capacity_per_sqkm: 3
    corine:
      grid_codes:
      - 23
      - 24
      - 25
      - 29
      distance: 1000
      distance_grid_codes:
      - 1
      - 2
      - 3
      - 4
      - 5
      - 6
    natura: true
    potential: simple
    clip_p_max_pu: 0.01
  onwind-2:
    cutout: europe-2015-era5
    resource:
      method: wind
      turbine: Vestas_V80_2MW_gridstreamer
    capacity_per_sqkm: 3
    corine:
      grid_codes:
      - 15
      - 16
      - 17
      - 19
      - 20
      - 21
      - 22
      distance: 1000
      distance_grid_codes:
      - 1
      - 2
      - 3
      - 4
      - 5
      - 6
    natura: true
    potential: simple
    clip_p_max_pu: 0.01
  onwind-3:
    cutout: europe-2015-era5
    resource:
      method: wind
      turbine: Enercon_E82_3000kW
    capacity_per_sqkm: 3
    corine:
      grid_codes:
      - 12
      - 13
      - 14
      distance: 1000
      distance_grid_codes:
      - 1
      - 2
      - 3
      - 4
      - 5
      - 6
    natura: true
    potential: simple
    clip_p_max_pu: 0.01
  onwind-4:
    cutout: europe-2015-era5
    resource:
      method: wind
      turbine: Enercon_E126_7500kW
    capacity_per_sqkm: 3
    corine:
      grid_codes:
      - 18
      - 26
      - 27
      - 28
      - 30
      - 31
      - 32
      distance: 1000
      distance_grid_codes:
      - 1
      - 2
      - 3
      - 4
      - 5
      - 6
    natura: true
    potential: simple
    clip_p_max_pu: 0.01
  offwind-ac:
    cutout: europe-2015-era5
    resource:
      method: wind
      turbine: NREL_ReferenceTurbine_5MW_offshore
    capacity_per_sqkm: 3
    corine:
    - 44
    - 255
    natura: true
    max_depth: 50
    max_shore_distance: 50000
    potential: simple
    clip_p_max_pu: 0.01
  offwind-dc:
    cutout: europe-2015-era5
    resource:
      method: wind
      turbine: NREL_ReferenceTurbine_5MW_offshore
    capacity_per_sqkm: 3
    corine:
    - 44
    - 255
    natura: true
    max_depth: 50
    min_shore_distance: 50000
    potential: simple
    clip_p_max_pu: 0.01
  solar-utility:
    cutout: europe-2015-era5
    resource:
      method: pv
      panel: CSi
      orientation:
        slope: 35.0
        azimuth: 180.0
    capacity_per_sqkm: 1.7
    correction_factor: 0.854337
    corine:
    - 4
    - 7
    - 8
    - 12
    - 13
    - 14
    - 15
    - 16
    - 17
    - 18
    - 19
    - 20
    - 26
    - 31
    - 32
    natura: true
    potential: simple
    clip_p_max_pu: 0.01
  solar-rooftop:
    cutout: europe-2015-era5
    resource:
      method: pv
      panel: CSi
      orientation:
        slope: 35.0
        azimuth: 180.0
    capacity_per_sqkm: 1.7
    correction_factor: 0.854337
    corine:
    - 1
    - 2
    - 3
    - 5
    - 6
    - 11
    natura: true
    potential: simple
    clip_p_max_pu: 0.01
  hydro:
    cutout: europe-2015-era5
    carriers:
    - ror
    - PHS
    - hydro
    PHS_max_hours: 6
    hydro_max_hours: energy_capacity_totals_by_country
    clip_min_inflow: 1.0
lines:
  types:
    220.0: Al/St 240/40 2-bundle 220.0
    300.0: Al/St 240/40 3-bundle 300.0
    380.0: Al/St 240/40 4-bundle 380.0
  s_max_pu: 0.7
  length_factor: 1.25
  under_construction: keep
links:
  p_max_pu: 1.0
  include_tyndp: true
  under_construction: keep
transformers:
  x: 0.1
  s_nom: 2000.0
  type: ''
pipelines:
  update_AT: true
  small_retrofit_pipelines: 100.0
  small_gas_pipelines: 500.0
  gas_h2_capacity_ratio: 0.85
PECD_wind_PV: true
sector:
  central: true
  urban_decentral_countries: []
  dsm_restriction_value: 0.75
  dsm_restriction_time: 7
  transport_heating_deadband_upper: 20.0
  transport_heating_deadband_lower: 15.0
  ICE_lower_degree_factor: 0.375
  ICE_upper_degree_factor: 1.6
  EV_lower_degree_factor: 0.98
  EV_upper_degree_factor: 0.63
  district_heating_loss: 0.15
  bev: true
  bev_availability: 0.5
  bev_power_per_year: 2.598
  bev_dsm_battery: 0.05
  v2g: false
  time_dep_hp_cop: true
  space_heating_fraction: 1.0
  tes: true
  tes_tau: 3.0
  boilers: true
  solid biomass boiler: true
  coal boiler: false
  coal chp: true
  oil_power_plant: true
  biogas_import: true
  oil_boilers: true
  chp: true
  chp_parameters:
    c_v: 0.15
    p_nom_ratio: 1.0
  solar_thermal: true
  solar_cf_correction: 0.788457
  marginal_cost_storage: 0.0
  methanation: true
  helmeth: false
  dac: true
  co2_vent: true
  SMR: true
  cc_fraction: 0.9
  h2_non_eu_import: true
  h2_NAF_import_countries:
  - IT
  - ES
  - GR
  h2_RU_import_countries:
  - DE
  - SK
  - PL
  h2_non_eu_import_p_min_pu: 0.125
  split_electrolysis: false
  hydrogen_underground_storage: true
  gas_network: true
  gas_pipe_p_max_pu: 0.6
  h2_retrofit: true
  h2_retro_pipe_p_max_pu: 0.6
  h2_retro_pipeline_utilisation: 0.85
  h2_new_pipe_p_max_pu: 0.6
  use_fischer_tropsch_waste_heat: true
  use_fuel_cell_waste_heat: true
  electricity_distribution_grid: false
  electricity_distribution_grid_cost_factor: 1.0
  value_of_lost_load: true
  grid_losses: true
  grid_losses_factor: 0.055
  ft_import: true
  help_store_restriction_time: 23
  help_store_restriction_value: 0.05
  pemfc_h2: true
  pemfc_ch4: true
  ht_electrolysis_share_heat: 0.2
  usable_industrial_waste_heat: 0.1
costs:
  year: 2015
  lifetime: 25
  discountrate: 0.07
  USD2013_to_EUR2013: 0.7532
  emission_prices:
    co2: 0.0
  lines:
    length_factor: 1.25
solving:
  options:
    formulation: kirchhoff
    load_shedding: false
    noisy_costs: true
    min_iterations: 1
    max_iterations: 1
    clip_p_max_pu: 0.01
    use_linmod: true
  solver:
    name: gurobi
    method: 2
    crossover: 0
    BarConvTol: 0.0001
    FeasibilityTol: 0.0001
    AggFill: 0
    PreDual: 0
    BarHomogeneous: 1
    Aggregate: 0
    NumericFocus: 1
    IISMethod: 1
    threads: 8
    InfeasHandling: 2
industry:
  DRI_ratio: 1.0
  H2_DRI: 2.3
  elec_DRI: 0.23
  clinker_factor: 0.81
  paper_heat: 0.8
  paper_hp: 2
  H2_for_NH3: 1.06878951066018
  H2_from_SMR_2015: 0.7032
  SMR_Efficiency_2015: 0.74
  process_emissions_basic_chemicals: 0.51
  efficiency_gain: 0.7
  elec_gas_split: 0.7
  steel_2050: 80
  ammonium_2050: 22
plotting:
  map:
    figsize:
    - 7
    - 7
    boundaries:
    - -10.2
    - 29
    - 35
    - 72
    p_nom:
      bus_size_factor: 50000.0
      linewidth_factor: 3000.0
  costs_max: 1200
  costs_threshold: 1
  energy_max: 20000.0
  energy_min: -15000.0
  energy_threshold: 50.0
  vre_techs:
  - onwind-1
  - onwind-2
  - onwind-3
  - onwind-4
  - offwind-ac
  - offwind-dc
  - solar-rooftop
  - solar-utility
  - ror
  renewable_storage_techs:
  - PHS
  - hydro
  conv_techs:
  - OCGT
  - CCGT
  - Nuclear
  - Coal
  storage_techs:
  - hydro+PHS
  - battery
  - H2
  load_carriers:
  - AC load
  AC_carriers:
  - AC line
  - AC transformer
  link_carriers:
  - DC line
  - Converter AC-DC
  heat_links:
  - heat pump
  - resistive heater
  - CHP heat
  - CHP electric
  - gas boiler
  - central heat pump
  - central resistive heater
  - central CHP heat
  - central CHP electric
  - central gas boiler
  heat_generators:
  - gas boiler
  - central gas boiler
  - solar thermal collector
  - central solar thermal collector
mobility:
  add_bev: true
  add_fcev: true
  add_ice: true
  add_cng: true
  add_hev: true
  add_phev: true
  phev_fischer_tropsch_split:
    long: 0.9
    short: 0.3
  efficiency_bev_charger: 0.9
  fraction_short_distance: 0.75
"""
        return yaml.safe_load(s)

    monkeypatch.setattr(mod, "read_config", mock)
    # mod.read_config = mock_read_config
