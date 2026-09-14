import logging
import os
import sys

import pandas as pd
import pypsa

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + "/../.."))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from export_ariadne_variables import _get_fuel_fractions, get_export_import

from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)

pypsa.options.params.statistics.round = 10
logger = logging.getLogger(__name__)

# Defining global variables

TWh2MWh = 1e6
MWh2TWh = 1e-6
MW2GW = 1e-3
MW2TW = 1e-6
t2Mt = 1e-6
toe_to_MWh = 11.630  # GWh/ktoe OR MWh/toe


EUR20TOEUR23 = 1.1076


def get_electricity_balance(
    n,
    region="DE",
    return_intermediates=False,
):

    carrier_to_label = {
        # Secondary Energy|Electricity|Carrier
        "hydro": "Secondary Energy|Electricity|Hydro",
        "ror": "Secondary Energy|Electricity|Hydro",
        "onwind": "Secondary Energy|Electricity|Wind|Onshore",
        "offwind-ac": "Secondary Energy|Electricity|Wind|Offshore",
        "offwind-dc": "Secondary Energy|Electricity|Wind|Offshore",
        "solar": "Secondary Energy|Electricity|Solar",
        "solar rooftop": "Secondary Energy|Electricity|Solar",
        "solar-hsat": "Secondary Energy|Electricity|Solar",
        "nuclear": "Secondary Energy|Electricity|Nuclear",
        "solid biomass": "Secondary Energy|Electricity|Biomass",
        "biogas CHP": "Secondary Energy|Electricity|Biomass",
        "urban central solid biomass CHP": "Secondary Energy|Electricity|Biomass",
        "urban central solid biomass CHP CC": "Secondary Energy|Electricity|Biomass",
        "H2 CCGT": "Secondary Energy|Electricity|Hydrogen",
        "H2 Fuel Cell": "Secondary Energy|Electricity|Hydrogen",
        "H2 OCGT": "Secondary Energy|Electricity|Hydrogen",
        "H2 retrofit OCGT": "Secondary Energy|Electricity|Hydrogen",
        "H2 retrofit CCGT": "Secondary Energy|Electricity|Hydrogen",
        "oil": "Secondary Energy|Electricity|Oil",
        "urban central H2 CHP": "Secondary Energy|Electricity|Hydrogen",
        "urban central H2 retrofit CHP": "Secondary Energy|Electricity|Hydrogen",
        "urban central oil CHP": "Secondary Energy|Electricity|Oil",
        "waste CHP": "Secondary Energy|Electricity|Waste",
        "waste CHP CC": "Secondary Energy|Electricity|Waste",
        "CCGT": "Secondary Energy|Electricity|Gas",
        "OCGT": "Secondary Energy|Electricity|Gas",
        "urban central gas CHP": "Secondary Energy|Electricity|Gas",
        "urban central gas CHP CC": "Secondary Energy|Electricity|Gas",
        "coal": "Secondary Energy|Electricity|Coal|Hard Coal",
        "lignite": "Secondary Energy|Electricity|Coal|Lignite",
        "urban central coal CHP": "Secondary Energy|Electricity|Coal|Hard Coal",
        "urban central lignite CHP": "Secondary Energy|Electricity|Coal|Lignite",
        "electricity distribution grid": "Secondary Energy|Electricity|Transmission Losses",  # in the joint balance at the AC + low voltage bus, the distribution grid should only show up for losses
        # Secondary Energy Input|Electricity|Carrier
        "methanolisation": "Secondary Energy Input|Electricity|Liquids",
        "urban central air heat pump": "Secondary Energy Input|Electricity|Heat",
        "urban central resistive heater": "Secondary Energy Input|Electricity|Heat",
        "H2 Electrolysis": "Secondary Energy Input|Electricity|Hydrogen",
        "H2 pipeline": "Secondary Energy Input|Electricity|Hydrogen",
        "H2 pipeline retrofitted": "Secondary Energy Input|Electricity|Hydrogen",
        "H2 pipeline (Kernnetz)": "Secondary Energy Input|Electricity|Hydrogen",
        # Final Energy|Sector|Electricity
        "DAC": "Final Energy|Carbon Dioxide Removal|Electricity",
        "electricity": "Final Energy|Residential and Commercial|Electricity",
        "rural air heat pump": "Final Energy|Residential and Commercial|Electricity",  # "Stromnachfrage (fix)",
        "rural ground heat pump": "Final Energy|Residential and Commercial|Electricity",  # "Stromnachfrage (fix)",
        "rural resistive heater": "Final Energy|Residential and Commercial|Electricity",  # "Stromnachfrage (fix)",
        "urban decentral air heat pump": "Final Energy|Residential and Commercial|Electricity",  # "Stromnachfrage (fix)",
        "urban decentral resistive heater": "Final Energy|Residential and Commercial|Electricity",  # "Stromnachfrage (fix)",
        "industry electricity": "Final Energy|Industry|Electricity",
        "BEV charger": "Final Energy|Transportation|Electricity",
        "agriculture electricity": "Final Energy|Agriculture|Electricity",
        "agriculture machinery electric": "Final Energy|Agriculture|Electricity",
        # Other
        "AC": "Trade|Electricity",
        "DC": "Trade|Electricity",
        "PHS": "Storage|Electricity",
        "battery charger": "Storage|Electricity",
        "battery discharger": "Storage|Electricity",
        "home battery charger": "Storage|Electricity",
        "home battery discharger": "Storage|Electricity",
    }

    electricity_balance = n.statistics.energy_balance(
        bus_carrier=["AC", "low voltage"],
        groupby=["bus", "name", "carrier"],
        nice_names=False,
        aggregate_time=False,
    )
    electricity_balance_by_carrier = (
        electricity_balance[
            electricity_balance.index.get_level_values("bus").str.startswith(region)
        ]
        .groupby("carrier")
        .sum()
    )

    missing_carriers = sorted(
        set(electricity_balance_by_carrier.index) - set(carrier_to_label)
    )
    assert not missing_carriers, (
        f"Missing carriers in electricity balance: {missing_carriers}"
    )
    electricity_balance_by_variable = (
        electricity_balance_by_carrier.rename(index=carrier_to_label)
        .groupby(level=0)
        .sum()
    )
    electricity_variables = electricity_balance_by_variable.copy()
    include_gross_generation = False
    if include_gross_generation:
        gross_generation_factors = {
            "Secondary Energy|Electricity|Gas": 1.03,
            "Secondary Energy|Electricity|Coal": 1.10,
            "Secondary Energy|Electricity|Wind": 1.02,
            "Secondary Energy|Electricity|Solar": 1.02,
            "Secondary Energy|Electricity|Biomass|Gaseous and Liquid": 1.05,
            "Secondary Energy|Electricity|Biomass|Solid": 1.05,
            "Secondary Energy|Electricity|Waste": 1.27,
        }
        total_generation_losses = electricity_variables.iloc[0] * 0
        for carrier, gross_generation_factor in gross_generation_factors.items():
            generation_losses = electricity_balance_by_variable.loc[carrier] * (
                gross_generation_factor - 1
            )
            electricity_variables.loc[carrier] += generation_losses
            total_generation_losses += generation_losses

        electricity_variables.loc[
            "Secondary Energy|Electricity|Plant Losses"
        ] = -total_generation_losses

    ac_exports, ac_imports = get_export_import(n, region, ["AC"], aggregate=False)
    dc_exports, dc_imports = get_export_import(n, region, ["DC"], aggregate=False)
    hourly_exports = ac_exports.add(dc_exports, fill_value=0).div(
        n.snapshot_weightings.generators, axis=0
    )
    hourly_imports = ac_imports.add(dc_imports, fill_value=0).div(
        n.snapshot_weightings.generators, axis=0
    )

    hourly_transmission_losses = (
        hourly_imports.sum(axis=1)
        - hourly_exports.sum(axis=1)
        - electricity_variables.loc["Trade|Electricity"]
    )
    assert all(hourly_transmission_losses >= 0), (
        "Net imports cannot be negative after accounting for Stromimporte"
    )
    electricity_variables.loc[
        "Trade|Secondary Energy|Electricity|Gross Import|Volume"
    ] = hourly_imports.sum(axis=1)
    # TODO this is a slight misnomer, as the template asks for the gross export, not the net export )
    electricity_variables.loc[
        "Trade|Secondary Energy|Electricity|Volume"
    ] = -hourly_exports.sum(axis=1)
    electricity_variables.loc["Secondary Energy|Electricity|Transmission Losses"] -= (
        hourly_transmission_losses
    )
    electricity_variables.drop("Trade|Electricity", inplace=True)
    # TODO This could be improved by getting the fuel fraction at every hour
    biomethane_share = _get_fuel_fractions(n, region, "gas").get("Biomass", 0)
    biomethane_generation = (
        electricity_variables.loc["Secondary Energy|Electricity|Gas"] * biomethane_share
    )
    electricity_variables.loc["Secondary Energy|Electricity|Gas"] -= (
        biomethane_generation
    )
    electricity_variables.loc["Secondary Energy|Electricity|Gas|Natural Gas"] = (
        electricity_variables.loc["Secondary Energy|Electricity|Gas"]
    )
    electricity_variables.drop("Secondary Energy|Electricity|Gas", inplace=True)
    # TODO this is a misnomer, it should be Gas|Biomethane, but the template is not consistent with this
    electricity_variables.loc["Secondary Energy|Electricity|Gas|Efuel"] = (
        biomethane_generation
    )

    electricity_variables.loc["Secondary Energy|Electricity|Storage"] = (
        electricity_variables.loc["Storage|Electricity"].clip(upper=0)
    )
    electricity_variables.loc["Secondary Energy|Electricity|Storage Discharge"] = (
        electricity_variables.loc["Storage|Electricity"].clip(lower=0)
    )
    electricity_variables.drop("Storage|Electricity", inplace=True)

    if return_intermediates:
        return {
            "electricity_balance_by_variable": electricity_balance_by_variable,
            "hourly_imports": hourly_imports,
            "hourly_exports": hourly_exports,
            "hourly_transmission_losses": hourly_transmission_losses,
            "biomethane_generation": biomethane_generation,
            "electricity_variables": electricity_variables,
        }

    return electricity_variables


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_dynamic_variables",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Bal_v5",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    ariadne_template = pd.read_excel(snakemake.input.template, sheet_name=None)
    var2unit = ariadne_template["variable_definitions"].set_index("Variable")["Unit"]
    var2unit = var2unit.str.replace("GJ", "MWh")
    var2unit = var2unit.str.replace("PJ", "TWh")

    # Load data
    networks = [pypsa.Network(fn) for fn in snakemake.input.networks]

    modelyears = [fn[-7:-3] for fn in snakemake.input.networks]

    if "debug" == "debug":  # For debugging
        var = pd.Series()
        idx = 1
        n = networks[idx]
        region = "DE"
        kwargs = {
            "groupby": ["bus", "carrier"],
            "at_port": True,
            "nice_names": False,
        }

    yearly_dfs = []
    for n, year in zip(networks, modelyears):
        logger.info(f"Getting data for year {year}...")
        eleba = get_electricity_balance(n, region="DE")
        _df = (
            eleba.rename_axis(index="Variable", columns="Time")
            .stack()
            .rename("Value")
            .reset_index()
        )
        _df["Time"] = _df["Time"].map(
            lambda timestamp: timestamp.replace(year=int(year))
        )
        _df["Model"] = "PyPSA-Eur v0.10"
        _df["Scenario"] = snakemake.wildcards.run
        _df["Region"] = "Deutschland"
        _df["Unit"] = _df["Variable"].map(var2unit).fillna("NA")
        yearly_dfs.append(_df)

    df = pd.concat(yearly_dfs, ignore_index=True)
    not_in_template = df.loc[df["Unit"] == "NA"].index
    logger.info(
        "Dropping variables which are not in the template: %s", not_in_template.tolist()
    )
    ariadne_df = df.drop(not_in_template)

    meta = pd.Series(
        {
            "Model": "PyPSA-Eur v0.10",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "preliminary",
            "Internal usage within Kopernikus AG Szenarien": "yes",
            "Release for publication": "no",
        }
    )

    # For export to the Ariadne-internal DB, convert most Wh-based entries to J
    ariadne_df.loc[ariadne_df["Unit"] == "TWh/yr", "Value"] *= 3.6
    ariadne_df.loc[ariadne_df["Unit"] == "TWh/yr", "Unit"] = "PJ/yr"
    ariadne_df.loc[ariadne_df["Unit"] == "EUR2020/MWh", "Value"] /= 3.6
    ariadne_df.loc[ariadne_df["Unit"] == "EUR2020/MWh", "Unit"] = "EUR2020/GJ"

    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        ariadne_df.round(5).to_excel(writer, sheet_name="data", index=False)
        meta.to_frame().T.to_excel(writer, sheet_name="meta", index=False)
