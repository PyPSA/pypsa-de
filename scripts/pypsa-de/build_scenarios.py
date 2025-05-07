# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# This script reads in data from the IIASA database to create the scenario.yaml file
import logging
from pathlib import Path

import pandas as pd
import ruamel.yaml

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def get_transport_growth(df, planning_horizons):
    aviation = df.loc["Final Energy|Bunkers|Aviation", "TWh/yr"]

    aviation[2020] = 111.25  # Ariadne2-internal DB, Aladin model
    aviation_growth_factor = aviation / aviation[2020]

    return aviation_growth_factor[planning_horizons]


def get_primary_steel_share(df, planning_horizons):
    # Get share of primary steel production
    model = snakemake.params.leitmodelle["industry"]
    total_steel = df.loc[model, "Production|Steel"]
    primary_steel = df.loc[model, "Production|Steel|Primary"]

    total_steel[2020] = 40.621  # Ariadne2-internal DB, FORECAST, 2021
    primary_steel[2020] = 28.53  # Ariadne2-internal DB, FORECAST, 2021

    primary_steel_share = primary_steel / total_steel
    primary_steel_share = primary_steel_share[planning_horizons]

    return primary_steel_share.set_index(pd.Index(["Primary_Steel_Share"]))


def get_DRI_share(df, planning_horizons):
    # Get share of DRI steel production
    model = "FORECAST v1.0"
    total_steel = df.loc[model, "Production|Steel|Primary"]
    # Assuming that only hydrogen DRI steel is sustainable and DRI using natural gas is phased out
    DRI_steel = df.loc[model, "Production|Steel|Primary|Direct Reduction Hydrogen"]

    total_steel[2020] = 40.621  # Ariadne2-internal DB, FORECAST, 2021
    DRI_steel[2020] = 0  # Ariadne2-internal DB, FORECAST, 2021

    DRI_steel_share = DRI_steel / total_steel
    DRI_steel_share = DRI_steel_share[planning_horizons]

    return DRI_steel_share.set_index(pd.Index(["DRI_Steel_Share"]))


def get_co2_budget(df, source):
    # relative to the DE emissions in 1990 *including bunkers*; also
    # account for non-CO2 GHG and allow extra room for international
    # bunkers which are excluded from the national targets

    # Baseline emission in DE in 1990 in Mt as understood by the KSG and by PyPSA
    baseline_co2 = 1251
    baseline_pypsa = 1052
    if source == "KSG":
        ## GHG targets according to KSG
        initial_years_co2 = pd.Series(
            index=[2020, 2025, 2030],
            data=[813, 643, 438],
        )

        later_years_co2 = pd.Series(
            index=[2035, 2040, 2045, 2050],
            data=[0.77, 0.88, 1.0, 1.0],
        )

        targets_co2 = pd.concat(
            [initial_years_co2, (1 - later_years_co2) * baseline_co2],
        )
    elif source == "UBA":
        ## For Zielverfehlungsszenarien use UBA Projektionsbericht
        targets_co2 = pd.Series(
            index=[2020, 2025, 2030, 2035, 2040, 2045, 2050],
            data=[813, 655, 455, 309, 210, 169, 157],
        )
    else:
        raise ValueError("Invalid source for CO2 budget.")
    ## Compute nonco2 from Ariadne-Leitmodell (REMIND)

    # co2 = (
    #     df.loc["Emissions|CO2 incl Bunkers","Mt CO2/yr"]
    #     - df.loc["Emissions|CO2|Land-Use Change","Mt CO2-equiv/yr"]
    #     - df.loc["Emissions|CO2|Energy|Demand|Bunkers","Mt CO2/yr"]
    # )
    # ghg = (
    #     df.loc["Emissions|Kyoto Gases","Mt CO2-equiv/yr"]
    #     - df.loc["Emissions|Kyoto Gases|Land-Use Change","Mt CO2-equiv/yr"]
    #     # No Kyoto Gas emissions for Bunkers recorded in Ariadne DB
    # )

    try:
        co2_land_use_change = df.loc["Emissions|CO2|Land-Use Change", "Mt CO2-equiv/yr"]
    except KeyError:  # Key not in Ariadne public database
        co2_land_use_change = df.loc["Emissions|CO2|AFOLU", "Mt CO2/yr"]

    co2 = df.loc["Emissions|CO2", "Mt CO2/yr"] - co2_land_use_change

    try:
        kyoto_land_use_change = df.loc[
            "Emissions|Kyoto Gases|Land-Use Change", "Mt CO2-equiv/yr"
        ]
    except KeyError:  # Key not in Ariadne public database
        # Guesstimate of difference from Ariadne 2 data
        kyoto_land_use_change = co2_land_use_change + 4.5

    ghg = df.loc["Emissions|Kyoto Gases", "Mt CO2-equiv/yr"] - kyoto_land_use_change

    nonco2 = ghg - co2

    ## PyPSA disregards nonco2 GHG emissions, but includes bunkers

    targets_pypsa = targets_co2 - nonco2

    target_fractions_pypsa = targets_pypsa.loc[targets_co2.index] / baseline_pypsa
    target_fractions_pypsa[2020] = (
        0.671  # Hard-coded based on REMIND data from ariadne2-internal DB
    )

    return target_fractions_pypsa.round(3)


def write_to_scenario_yaml(input, output, scenarios, df):
    # read in yaml file
    yaml = ruamel.yaml.YAML()
    file_path = Path(input)
    config = yaml.load(file_path)
    for scenario in scenarios:
        reference_scenario = (
            config[scenario]
            .get("iiasa_database", {})
            .get("reference_scenario", "KN2045_Mix")
        )

        planning_horizons = [
            2020,
            2025,
            2030,
            2035,
            2040,
            2045,
            2050,
        ]
        logger.info(
            "Using hard-coded values for the year 2020 for aviation demand, steel shares and non-co2 emissions. Source: Model results in the Ariadne2-internal database"
        )

        aviation_demand_factor = get_transport_growth(
            df.loc[snakemake.params.leitmodelle["transport"], reference_scenario, :],
            planning_horizons,
        )

        co2_budget_source = config[scenario]["co2_budget_DE_source"]

        co2_budget_fractions = get_co2_budget(
            df.loc[snakemake.params.leitmodelle["general"], reference_scenario],
            co2_budget_source,
        )

        if not config[scenario].get("sector"):
            config[scenario]["sector"] = {}

        if config[scenario]["sector"].get("aviation_demand_factor") is not None:
            logger.warning(f"Overwriting aviation_demand_factor in {scenario} scenario")
        else:
            config[scenario]["sector"]["aviation_demand_factor"] = {}

        for year in planning_horizons:
            config[scenario]["sector"]["aviation_demand_factor"][year] = round(
                aviation_demand_factor.loc[year].item(), 4
            )

        st_primary_fraction = get_primary_steel_share(
            df.loc[:, reference_scenario, :], planning_horizons
        )

        dri_fraction = get_DRI_share(
            df.loc[:, reference_scenario, :], planning_horizons
        )
        if not config[scenario].get("industry"):
            config[scenario]["industry"] = {}

        if config[scenario]["industry"].get("St_primary_fraction") is not None:
            logger.warning(f"Overwriting St_primary_fraction in {scenario} scenario")
        else:
            config[scenario]["industry"]["St_primary_fraction"] = {}

        for year in st_primary_fraction.columns:
            config[scenario]["industry"]["St_primary_fraction"][year] = round(
                st_primary_fraction.loc["Primary_Steel_Share", year].item(), 4
            )

        if config[scenario]["industry"].get("DRI_fraction") is not None:
            logger.warning(f"Overwriting DRI_fraction in {scenario} scenario")
        else:
            config[scenario]["industry"]["DRI_fraction"] = {}

        for year in dri_fraction.columns:
            config[scenario]["industry"]["DRI_fraction"][year] = round(
                dri_fraction.loc["DRI_Steel_Share", year].item(), 4
            )
        if not config[scenario].get("solving"):
            config[scenario]["solving"] = {}
        if not config[scenario]["solving"].get("constraints"):
            config[scenario]["solving"]["constraints"] = {}
        if (
            config[scenario]["solving"]["constraints"].get("co2_budget_national")
            is not None
        ):
            logger.warning(f"Overwriting co2_budget_national in {scenario} scenario")
        else:
            config[scenario]["solving"]["constraints"]["co2_budget_national"] = {}

        for year, target in co2_budget_fractions.items():
            config[scenario]["solving"]["constraints"]["co2_budget_national"][year] = {}
            config[scenario]["solving"]["constraints"]["co2_budget_national"][year][
                "DE"
            ] = target

    # write back to yaml file
    yaml.dump(config, Path(output))


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("build_scenarios")

    configure_logging(snakemake)
    # Set USERNAME and PASSWORD for the Ariadne DB
    ariadne_db = pd.read_csv(
        snakemake.input.ariadne_database,
        index_col=["model", "scenario", "region", "variable", "unit"],
    )
    ariadne_db.columns = ariadne_db.columns.astype(int)

    df = ariadne_db.loc[:, :, "Deutschland"]

    scenarios = snakemake.params.scenarios

    input = snakemake.input.scenario_yaml
    output = snakemake.output.scenario_yaml

    # for scenario in scenarios:
    write_to_scenario_yaml(input, output, scenarios, df)
