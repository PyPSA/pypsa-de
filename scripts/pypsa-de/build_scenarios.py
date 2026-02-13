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
    aviation = df.xs(
        ("Final Energy|Bunkers|Aviation", "TWh/yr"), level=("variable", "unit")
    ).squeeze()

    aviation[2020] = 111.25  # Ariadne2-internal DB, Aladin model
    aviation_growth_factor = aviation / aviation[2020]

    return aviation_growth_factor[planning_horizons]


def get_primary_steel_share(df, planning_horizons):
    # Get share of primary steel production
    model = snakemake.params.leitmodelle["industry"]
    model_df = df.xs(model, level="model")
    total_steel = model_df.xs("Production|Steel", level="variable").squeeze()
    primary_steel = model_df.xs("Production|Steel|Primary", level="variable").squeeze()

    total_steel[2020] = 40.621  # Ariadne2-internal DB, FORECAST, 2021
    primary_steel[2020] = 28.53  # Ariadne2-internal DB, FORECAST, 2021

    primary_steel_share = primary_steel / total_steel
    primary_steel_share = primary_steel_share[planning_horizons]

    return pd.DataFrame([primary_steel_share], index=["Primary_Steel_Share"])


def get_DRI_share(df, planning_horizons):
    # Get share of DRI steel production
    model = "FORECAST v1.0"
    model_df = df.xs(model, level="model")
    total_steel = model_df.xs("Production|Steel|Primary", level="variable").squeeze()
    # Assuming that only hydrogen DRI steel is sustainable and DRI using natural gas is phased out
    DRI_steel = model_df.xs(
        "Production|Steel|Primary|Direct Reduction Hydrogen", level="variable"
    ).squeeze()

    total_steel[2020] = 40.621  # Ariadne2-internal DB, FORECAST, 2021
    DRI_steel[2020] = 0  # Ariadne2-internal DB, FORECAST, 2021

    DRI_steel_share = DRI_steel / total_steel
    DRI_steel_share = DRI_steel_share[planning_horizons]

    return pd.DataFrame([DRI_steel_share], index=["DRI_Steel_Share"])


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
            data=[
                813,  # historisch
                # 643,  # Ziel laut KSG
                640,  # geschätzt, https://www.agora-energiewende.de/fileadmin/Projekte/2025/2025-28_DE_JAW25/A-EW_391_Die_Energiewende_in_Deutschland_Stand_der_Dinge_2025_WEB.pdf,
                # 438,  # Ziel laut KSG
                456.8,  # Durch UBA anhand der Emissionsdaten 2024 vom 15.3.2025 rechnerisch ermittelte Vorschläge für die Anpassung der Jahresemissionsgesamtmengen nach § 4 Abs. 2 KSG
            ],
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

    try:
        co2_land_use_change = df.xs(
            ("Emissions|CO2|Land-Use Change", "Mt CO2-equiv/yr"),
            level=("variable", "unit"),
        ).squeeze()
    except KeyError:  # Key not in Ariadne public database
        co2_land_use_change = df.xs(
            ("Emissions|CO2|AFOLU", "Mt CO2/yr"), level=("variable", "unit")
        ).squeeze()

    co2 = (
        df.xs(("Emissions|CO2", "Mt CO2/yr"), level=("variable", "unit")).squeeze()
        - co2_land_use_change
    )

    try:
        kyoto_land_use_change = df.xs(
            ("Emissions|Kyoto Gases|Land-Use Change", "Mt CO2-equiv/yr"),
            level=("variable", "unit"),
        ).squeeze()
    except KeyError:  # Key not in Ariadne public database
        # Guesstimate of difference from Ariadne 2 data
        kyoto_land_use_change = co2_land_use_change + 4.5

    ghg = (
        df.xs(
            ("Emissions|Kyoto Gases", "Mt CO2-equiv/yr"),
            level=("variable", "unit"),
        ).squeeze()
        - kyoto_land_use_change
    )

    nonco2 = ghg - co2

    # Hard-code values for 2020 and 2025 from UBA reports / projections
    # Table 1, https://reportnet.europa.eu/public/dataflow/1478, GHG - CO2
    nonco2[2020] = 733.7 - 647.9
    nonco2[2025] = 628.8 - 554.6

    ## PyPSA disregards nonco2 GHG emissions, but includes bunkers

    targets_pypsa = targets_co2 - nonco2

    logger.info("Non-CO2 GHG emissions assumed (in Mt CO2-equiv/yr):")
    for year in nonco2.index:
        logger.info(f"{year}: {nonco2.loc[year]:.1f}")

    target_fractions_pypsa = targets_pypsa.loc[targets_co2.index] / baseline_pypsa

    return target_fractions_pypsa.round(3)


def write_weather_dependent_config(config, scenario, weather_year):
    # Insert weather-dependent configuration dynamically
    cutout_name = f"europe-{weather_year}-sarah3-era5"

    # atlite section
    config[scenario]["atlite"] = {
        "default_cutout": cutout_name,
        "cutouts": {
            cutout_name: {
                "module": ["sarah", "era5"],
                "x": [-12.0, 42.0],
                "y": [33.0, 72.0],
                "dx": 0.3,
                "dy": 0.3,
                "time": [str(weather_year), str(weather_year)],
            }
        },
    }

    # snapshots section
    config[scenario]["snapshots"] = {
        "start": f"{weather_year}-01-01",
        "end": f"{int(weather_year) + 1}-01-01",
        "inclusive": "left",
    }

    # renewable section
    config[scenario]["renewable"] = {
        "onwind": {"cutout": cutout_name},
        "offwind-ac": {"cutout": cutout_name},
        "offwind-dc": {"cutout": cutout_name},
        "offwind-float": {"cutout": cutout_name},
        "solar": {"cutout": cutout_name},
        "solar-hsat": {"cutout": cutout_name},
        "hydro": {"cutout": cutout_name},
    }

    # lines section
    config[scenario]["lines"] = {"dynamic_line_rating": {"cutout": cutout_name}}


def write_to_scenario_yaml(input, output, scenarios, df):
    # read in yaml file
    yaml = ruamel.yaml.YAML()
    file_path = Path(input)
    config = yaml.load(file_path)
    for scenario in scenarios:
        if config.get(scenario) is None:
            logger.warning(
                f"Found an empty scenario config for {scenario}. Using default config `pypsa.de.yaml`."
            )
            config[scenario] = {}
        if config[scenario].get("weather_year", False):
            weather_year = config[scenario]["weather_year"]
            default_weather_year = int(snakemake.config["snapshots"]["start"][:4])
            if (
                snakemake.config["run"]["shared_resources"]["policy"] != False
                and weather_year != default_weather_year
            ):
                raise ValueError(
                    f"The run uses shared resources, but weather year {weather_year} in scenario {scenario} does not match the start year of the snapshots {default_weather_year}. If you are running scenarios with multiple weather years, make sure to deactivate shared_resources!"
                )
            write_weather_dependent_config(config, scenario, weather_year)

        reference_scenario = (
            config[scenario]
            .get("pypsa-de", {})
            .get(
                "reference_scenario",
                snakemake.config["pypsa-de"]["reference_scenario"],
            )  # Using the default reference scenario from pypsa.de.yaml
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

        transport_df = df.xs(snakemake.params.leitmodelle["transport"], level="model")
        transport_df = transport_df.xs(reference_scenario, level="scenario")
        aviation_demand_factor = get_transport_growth(transport_df, planning_horizons)
        if not config[scenario].get("co2_budget_DE_source"):
            logger.info(
                f"No CO2 budget source for DE specified in the scenario config. Using KSG targets and REMIND emissions from {reference_scenario} for the {scenario} scenario."
            )
            co2_budget_source = "KSG"
        else:
            co2_budget_source = config[scenario]["co2_budget_DE_source"]

        general_df = df.xs(snakemake.params.leitmodelle["general"], level="model")
        general_df = general_df.xs(reference_scenario, level="scenario")
        co2_budget_fractions = get_co2_budget(general_df, co2_budget_source)

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

        scenario_df = df.xs(reference_scenario, level="scenario")
        st_primary_fraction = get_primary_steel_share(scenario_df, planning_horizons)

        dri_fraction = get_DRI_share(scenario_df, planning_horizons)
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
        if not config[scenario]["solving"]["constraints"].get("co2_budget_national"):
            config[scenario]["solving"]["constraints"]["co2_budget_national"] = {}
        if (
            config[scenario]["solving"]["constraints"]["co2_budget_national"].get("DE")
            is not None
        ):
            logger.warning(
                f"Overwriting co2_budget_national for DE in {scenario} scenario"
            )
        else:
            config[scenario]["solving"]["constraints"]["co2_budget_national"]["DE"] = {}

        for year, target in co2_budget_fractions.items():
            config[scenario]["solving"]["constraints"]["co2_budget_national"]["DE"][
                year
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

    df = ariadne_db.xs("Deutschland", level="region")

    scenarios = snakemake.params.scenarios

    input = snakemake.input.scenario_yaml
    output = snakemake.output.scenario_yaml

    # for scenario in scenarios:
    write_to_scenario_yaml(input, output, scenarios, df)
