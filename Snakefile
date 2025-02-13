# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from shutil import unpack_archive

from snakemake.utils import min_version


include: "submodules/pypsa-eur/rules/common.smk"


import yaml, sys
from shutil import move

sys.path.append("workflow/submodules/pypsa-eur/scripts")

from _helpers import path_provider, get_scenarios, get_rdir

min_version("8.5")


configfile: "workflow/submodules/pypsa-eur/config/config.default.yaml"
configfile: "config/config.yaml"
configfile: "config/config.personal.yaml"


run = config["run"]
scenarios = get_scenarios(run)  # global variable
RDIR = get_rdir(run)
policy = run["shared_resources"]["policy"]
exclude = run["shared_resources"]["exclude"]
logs = path_provider("logs/", RDIR, policy, exclude)
benchmarks = path_provider("benchmarks/", RDIR, policy, exclude)
resources = path_provider("resources/", RDIR, policy, exclude)

RESULTS = "results/" + RDIR


wildcard_constraints:
    simpl="[a-zA-Z0-9]*",
    clusters="[0-9]+(m|c)?|all",
    ll="(v|c)([0-9\.]+|opt)",
    opts="[-+a-zA-Z0-9\.]*",
    sector_opts="[-+a-zA-Z0-9\.\s]*",
    planning_horizons="[0-9]*",
    year="[0-9]*",


module pypsaeur:
    snakefile:
        "submodules/pypsa-eur/Snakefile"
    config:
        config


use rule * from pypsaeur


from pathlib import Path

data_dir = Path("workflow/submodules/pypsa-eur/data")


rule get_data:
    output:
        [
            str(Path("data") / p.relative_to(data_dir))
            for p in data_dir.rglob("*")
            if p.is_file()
        ],
    shell:
        """
        mkdir -p data
        cp -nR {data_dir}/. data/
        """


rule clean:
    message:
        "Remove all build results but keep downloaded data."
    run:
        import shutil

        shutil.rmtree("resources")
        shutil.rmtree("results")
        print("Data downloaded to data/ has not been cleaned.")


rule retrieve_egon_data:
    output:
        spatial="data/egon/demandregio_spatial_2018.json",
        mapping="data/egon/mapping_technologies.json",
    shell:
        """
        wget -O {output.spatial} "https://api.opendata.ffe.de/demandregio/demandregio_spatial?id_spatial=5&year=2018"
        wget -O {output.mapping} "https://api.opendata.ffe.de/demandregio/demandregio_spatial_description?id_spatial=5"
        """


rule retrieve_ariadne_database:
    params:
        db_name=config_provider("iiasa_database", "db_name"),
        leitmodelle=config_provider("iiasa_database", "leitmodelle"),
        scenarios=config_provider("iiasa_database", "scenarios"),
    output:
        data=resources("ariadne_database.csv"),
    log:
        "logs/retrieve_ariadne_database.log",
    resources:
        mem_mb=1000,
    script:
        "scripts/retrieve_ariadne_database.py"


def input_profile_offwind(w):
    return {
        f"profile_{tech}": resources("profile_{clusters}_" + tech + ".nc")
        for tech in ["offwind-ac", "offwind-dc", "offwind-float"]
        if (tech in config_provider("electricity", "renewable_carriers")(w))
    }


use rule prepare_sector_network from pypsaeur with:
    input:
        unpack(input_profile_offwind),
        **{
            k: v
            for k, v in rules.prepare_sector_network.input.items()
            if k != "district_heat_share"
        },
        district_heat_share=resources(
            "district_heat_share_base_s_{clusters}_{planning_horizons}-modified.csv"
        ),


rule modify_cost_data:
    params:
        file_path="ariadne-data/costs/",
        file_name="costs_{planning_horizons}.csv",
        cost_horizon=config_provider("costs", "horizon"),
        NEP=config_provider("costs", "NEP"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_price_add_on_fossils=config_provider("co2_price_add_on_fossils"),
    input:
        modifications=lambda w: (
            "ariadne-data/costs_2019-modifications.csv"
            if w.planning_horizons == "2020"
            and config_provider("energy", "energy_totals_year") == 2019
            else "ariadne-data/costs_{planning_horizons}-modifications.csv"
        ),
    output:
        resources("costs_{planning_horizons}.csv"),
    resources:
        mem_mb=1000,
    script:
        "scripts/modify_cost_data.py"


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cost_data", True):

    ruleorder: modify_cost_data > retrieve_cost_data


rule build_mobility_demand:
    params:
        db_name=config_provider("iiasa_database", "db_name"),
        reference_scenario=config_provider("iiasa_database", "reference_scenario"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        leitmodelle=config_provider("iiasa_database", "leitmodelle"),
    input:
        ariadne=resources("ariadne_database.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
    output:
        mobility_demand=resources(
            "mobility_demand_aladin_{clusters}_{planning_horizons}.csv"
        ),
    resources:
        mem_mb=1000,
    log:
        logs("build_mobility_demand_{clusters}_{planning_horizons}.log"),
    script:
        "scripts/build_mobility_demand.py"


rule build_egon_data:
    input:
        demandregio_spatial="data/egon/demandregio_spatial_2018.json",
        mapping_38_to_4=storage(
            "https://ffeopendatastorage.blob.core.windows.net/opendata/mapping_from_4_to_38.json",
            keep_local=True,
        ),
        mapping_technologies="data/egon/mapping_technologies.json",
        nuts3=resources("nuts3_shapes.geojson"),
    output:
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
    script:
        "scripts/build_egon_data.py"


baseyear_value = config["scenario"]["planning_horizons"][0]


rule add_district_heating_subnodes:
    params:
        district_heating=config_provider("sector", "district_heating"),
        baseyear=config_provider("scenario", "planning_horizons", 0),
    input:
        network=RESULTS
        + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
        nuts3=resources("nuts3_shapes.geojson"),
        regions_onshore=resources(
            "regions_onshore_base_s_{clusters}.geojson"
        ),
        fernwaermeatlas="data/fernwaermeatlas/fernwaermeatlas.xlsx",
        cities="data/fernwaermeatlas/cities_geolocations.geojson",
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        existing_heating_distribution=resources(
            f"existing_heating_distribution_base_s_{{clusters}}_{baseyear_value}.csv"
        ),
        lau=storage(
            "https://gisco-services.ec.europa.eu/distribution/v2/lau/download/ref-lau-2021-01m.geojson.zip",
            keep_local=True,
        ),
    output:
        network=RESULTS
        + "prenetworks/base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        district_heating_subnodes=resources(
            "district_heating_subnodes_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.geojson"
        ),
        cop_profiles_extended=resources(
            "cop_profiles_base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        existing_heating_distribution_extended=(
            resources(
                "existing_heating_distribution_base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.csv"
            )
            if baseyear_value != "{planning_horizons}"
            else []
        ),
    resources:
        mem_mb=1000,
    script:
        "scripts/add_district_heating_subnodes.py"


ruleorder: modify_district_heat_share > build_district_heat_share


rule modify_district_heat_share:
    params:
        district_heating=config_provider("sector", "district_heating"),
    input:
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        district_heat_share=resources(
            "district_heat_share_base_s_{clusters}_{planning_horizons}.csv"
        ),
    output:
        district_heat_share=resources(
            "district_heat_share_base_s_{clusters}_{planning_horizons}-modified.csv"
        ),
    resources:
        mem_mb=1000,
    script:
        "scripts/modify_district_heat_share.py"


rule modify_prenetwork:
    params:
        enable_kernnetz=config_provider("wasserstoff_kernnetz", "enable"),
        costs=config_provider("costs"),
        max_hours=config_provider("electricity", "max_hours"),
        length_factor=config_provider("lines", "length_factor"),
        technology_occurrence=config_provider("first_technology_occurrence"),
        fossil_boiler_ban=config_provider("new_decentral_fossil_boiler_ban"),
        coal_ban=config_provider("coal_generation_ban"),
        nuclear_ban=config_provider("nuclear_generation_ban"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        H2_transmission_efficiency=config_provider(
            "sector", "transmission_efficiency", "H2 pipeline"
        ),
        H2_retrofit=config_provider("sector", "H2_retrofit"),
        H2_retrofit_capacity_per_CH4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        transmission_costs=config_provider("costs", "transmission"),
        biogas_must_run=config_provider("must_run_biogas"),
        clustering=config_provider("clustering", "temporal", "resolution_sector"),
        H2_plants=config_provider("electricity", "H2_plants_DE"),
        land_transport_electric_share=config_provider(
            "sector", "land_transport_electric_share"
        ),
        onshore_nep_force=config_provider("onshore_nep_force"),
        offshore_nep_force=config_provider("offshore_nep_force"),
        shipping_methanol_efficiency=config_provider(
            "sector", "shipping_methanol_efficiency"
        ),
        shipping_oil_efficiency=config_provider("sector", "shipping_oil_efficiency"),
        shipping_methanol_share=config_provider("sector", "shipping_methanol_share"),
        mwh_meoh_per_tco2=config_provider("sector", "MWh_MeOH_per_tCO2"),
    input:
        costs_modifications="ariadne-data/costs_{planning_horizons}-modifications.csv",
        network=RESULTS
        + "prenetworks-brownfield/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        wkn=(
            resources("wasserstoff_kernnetz_base_s_{clusters}.csv")
            if config_provider("wasserstoff_kernnetz", "enable")
            else []
        ),
        costs=resources("costs_{planning_horizons}.csv"),
        aladin_demand=resources(
            "mobility_demand_aladin_{clusters}_{planning_horizons}.csv"
        ),
        transport_data=resources("transport_data_s_{clusters}.csv"),
        biomass_potentials=resources(
            "biomass_potentials_s_{clusters}_{planning_horizons}.csv"
        ),
        industrial_demand=resources(
            "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
        ),
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        shipping_demand=resources("shipping_demand_s_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        offshore_connection_points="ariadne-data/offshore_connection_points.csv",
    output:
        network=RESULTS
        + "prenetworks-final/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    resources:
        mem_mb=1000,
    log:
        RESULTS
        + "logs/modify_prenetwork_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    script:
        "scripts/modify_prenetwork.py"


ruleorder: modify_industry_demand > build_industrial_production_per_country_tomorrow


use rule solve_sector_network_myopic from pypsaeur with:
    params:
        **{
            k: v
            for k, v in rules.solve_sector_network_myopic.params.items()
            if k != "custom_extra_functionality"
        },
        custom_extra_functionality=os.path.join(
            os.path.dirname(workflow.snakefile), "scripts/additional_functionality.py"
        ),
    input:
        **{
            k: v
            for k, v in rules.solve_sector_network_myopic.input.items()
            if k != "network"
        },
        network=RESULTS
        + "prenetworks-final/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        co2_totals_name=resources("co2_totals.csv"),


rule modify_existing_heating:
    params:
        iiasa_reference_scenario=config_provider("iiasa_database", "reference_scenario"),
        fallback_reference_scenario=config_provider(
            "iiasa_database", "fallback_reference_scenario"
        ),
    input:
        ariadne=resources("ariadne_database.csv"),
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
    output:
        existing_heating=resources("existing_heating.csv"),
    resources:
        mem_mb=1000,
    script:
        "scripts/modify_existing_heating.py"


rule retrieve_mastr:
    input:
        storage(
            "https://zenodo.org/records/8225106/files/bnetza_open_mastr_2023-08-08_B.zip",
            keep_local=True,
        ),
    params:
        "data/mastr",
    output:
        "data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
        "data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
    run:
        unpack_archive(input[0], params[0])


rule build_existing_chp_de:
    params:
        add_district_heating_subnodes=config_provider(
            "sector", "district_heating", "add_subnodes"
        ),
    input:
        mastr_biomass="data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
        mastr_combustion="data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
        plz_mapping=storage(
            "https://raw.githubusercontent.com/WZBSocialScienceCenter/plz_geocoord/master/plz_geocoord.csv",
            keep_local=True,
        ),
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        district_heating_subnodes=(
            resources(
                "district_heating_subnodes_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.geojson"
            )
            if config["sector"]["district_heating"].get("add_subnodes", True)
            else []
        ),
    output:
        german_chp=resources(
            "german_chp_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.csv"
        ),
    script:
        "scripts/build_existing_chp_de.py"


use rule add_existing_baseyear from pypsaeur with:
    params:
        **{k: v for k, v in rules.add_existing_baseyear.params.items()},
        add_district_heating_subnodes=config_provider(
            "sector", "district_heating", "add_subnodes"
        ),
    input:
        **{
            k: v
            for k, v in rules.add_existing_baseyear.input.items()
            if k != "network"
            and k != "cop_profiles"
            and k != "existing_heating_distribution"
        },
        network=(
            RESULTS
            + "prenetworks/base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
            if config["sector"]["district_heating"].get("add_subnodes", True)
            else RESULTS
            + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        cop_profiles=(
            resources(
                "cop_profiles_base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
            )
            if config["sector"]["district_heating"].get("add_subnodes", True)
            else resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc")
        ),
        custom_powerplants=resources(
            "german_chp_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.csv"
        ),
        existing_heating_distribution=(
            resources(
                "existing_heating_distribution_base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.csv"
            )
            if config["sector"]["district_heating"].get("add_subnodes", True)
            else resources(
                "existing_heating_distribution_base_s_{clusters}_{planning_horizons}.csv"
            )
        ),


def input_profile_tech_brownfield(w):
    return {
        f"profile_{tech}": resources("profile_{clusters}_" + tech + ".nc")
        for tech in config_provider("electricity", "renewable_carriers")(w)
        if tech != "hydro"
    }


use rule add_brownfield from pypsaeur with:
    input:
        unpack(input_profile_tech_brownfield),
        **{k: v for k, v in rules.add_brownfield.input.items() if k != "network"},
        network=(
            RESULTS
            + "prenetworks/base-extended_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
            if config["sector"]["district_heating"].get("add_subnodes", True)
            else RESULTS
            + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),


use rule build_existing_heating_distribution from pypsaeur with:
    input:
        **{
            k: v
            for k, v in rules.build_existing_heating_distribution.input.items()
            if k != "existing_heating"
        },
        existing_heating=resources("existing_heating.csv"),


rule modify_industry_demand:
    params:
        db_name=config_provider("iiasa_database", "db_name"),
    input:
        ariadne=resources("ariadne_database.csv"),
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}.csv"
        ),
    output:
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}-modified.csv"
        ),
    resources:
        mem_mb=1000,
    script:
        "scripts/modify_industry_demand.py"


use rule build_industrial_production_per_node from pypsaeur with:
    input:
        **{
            k: v
            for k, v in rules.build_industrial_production_per_node.input.items()
            if k != "industrial_production_per_country_tomorrow"
        },
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}-modified.csv"
        ),


rule modify_energy_totals:
    input:
        ariadne=resources("ariadne_database.csv"),
        energy_totals=resources("energy_totals.csv"),
    output:
        energy_totals=resources("energy_totals-modified.csv"),
    resources:
        mem_mb=1000,
    script:
        "scripts/modify_energy_totals.py"


use rule build_population_weighted_energy_totals from pypsaeur with:
    input:
        **{
            k: v
            for k, v in rules.build_population_weighted_energy_totals.input.items()
            if k != "energy_totals"
        },
        energy_totals=resources("energy_totals-modified.csv"),


rule build_wasserstoff_kernnetz:
    params:
        reload_locations=config_provider("wasserstoff_kernnetz", "reload_locations"),
    input:
        wasserstoff_kernnetz_1=storage(
            "https://fnb-gas.de/wp-content/uploads/2023/11/2023_11_15_Anlage2_Leitungsmeldungen_weiterer_potenzieller_Wasserstoffnetzbetreiber_Veroeffentlichung_final.xlsx",
            keep_local=True,
        ),
        wasserstoff_kernnetz_2=storage(
            "https://fnb-gas.de/wp-content/uploads/2023/11/2023_11_15_Anlage3_FNB_Massnahmenliste_Veroeffentlichung_final.xlsx",
            keep_local=True,
        ),
        gadm=storage(
            "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_DEU_1.json.zip",
            keep_local=True,
        ),
        locations="ariadne-data/wasserstoff_kernnetz/locations_wasserstoff_kernnetz.csv",
    output:
        cleaned_wasserstoff_kernnetz=resources("wasserstoff_kernnetz.csv"),
    script:
        "scripts/build_wasserstoff_kernnetz.py"


rule cluster_wasserstoff_kernnetz:
    params:
        kernnetz=config_provider("wasserstoff_kernnetz"),
    input:
        cleaned_h2_network=resources("wasserstoff_kernnetz.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
    output:
        clustered_h2_network=resources(
            "wasserstoff_kernnetz_base_s_{clusters}.csv"
        ),
    script:
        "scripts/cluster_wasserstoff_kernnetz.py"


rule download_ariadne_template:
    input:
        storage(
            "https://github.com/iiasa/ariadne-intern-workflow/raw/main/attachments/2024-05-29_template_Ariadne.xlsx",
            keep_local=True,
        ),
    output:
        resources("template_ariadne_database.xlsx"),
    run:
        move(input[0], output[0])


rule export_ariadne_variables:
    params:
        length_factor=config_provider("lines", "length_factor"),
        dg_cost_factor=config_provider(
            "sector", "electricity_distribution_grid_cost_factor"
        ),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        hours=config_provider("clustering", "temporal", "resolution_sector"),
        costs=config_provider("costs"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        co2_price_add_on_fossils=config_provider("co2_price_add_on_fossils"),
        co2_sequestration_cost=config_provider("sector", "co2_sequestration_cost"),
    input:
        template=resources("template_ariadne_database.xlsx"),
        industry_demands=expand(
            resources(
                "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
            ),
            **config["scenario"],
            allow_missing=True,
        ),
        networks=expand(
            RESULTS
            + "postnetworks/base_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            allow_missing=True,
        ),
        costs=expand(
            resources("costs_{planning_horizons}.csv"),
            **config["scenario"],
            allow_missing=True,
        ),
        industrial_production_per_country_tomorrow=expand(
            resources(
                "industrial_production_per_country_tomorrow_{planning_horizons}-modified.csv"
            ),
            **config["scenario"],
            allow_missing=True,
        ),
        industry_sector_ratios=expand(
            resources("industry_sector_ratios_{planning_horizons}.csv"),
            **config["scenario"],
            allow_missing=True,
        ),
        industrial_production=resources("industrial_production_per_country.csv"),
        energy_totals=resources("energy_totals.csv"),
    output:
        exported_variables=RESULTS + "ariadne/exported_variables.xlsx",
        exported_variables_full=RESULTS + "ariadne/exported_variables_full.xlsx",
    log:
        RESULTS + "logs/export_ariadne_variables.log",
    script:
        "scripts/export_ariadne_variables.py"


rule plot_ariadne_variables:
    params:
        iiasa_scenario=config_provider("iiasa_database", "reference_scenario"),
        fallback_reference_scenario=config_provider(
            "iiasa_database", "fallback_reference_scenario"
        ),
    input:
        exported_variables_full=RESULTS + "ariadne/exported_variables_full.xlsx",
        ariadne_database=resources("ariadne_database.csv"),
    output:
        primary_energy=RESULTS + "ariadne/primary_energy.png",
        primary_energy_detailed=RESULTS + "ariadne/primary_energy_detailed.png",
        secondary_energy=RESULTS + "ariadne/secondary_energy.png",
        secondary_energy_detailed=RESULTS + "ariadne/secondary_energy_detailed.png",
        final_energy=RESULTS + "ariadne/final_energy.png",
        final_energy_detailed=RESULTS + "ariadne/final_energy_detailed.png",
        capacity=RESULTS + "ariadne/capacity.png",
        capacity_detailed=RESULTS + "ariadne/capacity_detailed.png",
        energy_demand_emissions=RESULTS + "ariadne/energy_demand_emissions.png",
        energy_supply_emissions=RESULTS + "ariadne/energy_supply_emissions.png",
        co2_emissions=RESULTS + "ariadne/co2_emissions.png",
        primary_energy_price=RESULTS + "ariadne/primary_energy_price.png",
        secondary_energy_price=RESULTS + "ariadne/secondary_energy_price.png",
        #final_energy_residential_price = RESULTS + "ariadne/final_energy_residential_price.png",
        final_energy_industry_price=RESULTS + "ariadne/final_energy_industry_price.png",
        final_energy_transportation_price=RESULTS
        + "ariadne/final_energy_transportation_price.png",
        final_energy_residential_commercial_price=RESULTS
        + "ariadne/final_energy_residential_commercial_price.png",
        all_prices=RESULTS + "ariadne/all_prices.png",
        policy_carbon=RESULTS + "ariadne/policy_carbon.png",
        investment_energy_supply=RESULTS + "ariadne/investment_energy_supply.png",
        elec_val_2020=RESULTS + "ariadne/elec_val_2020.png",
        trade=RESULTS + "ariadne/trade.png",
    script:
        "scripts/plot_ariadne_variables.py"


rule ariadne_all:
    input:
        expand(RESULTS + "graphs/costs.svg", run=config_provider("run", "name")),
        expand(
            RESULTS + "ariadne/capacity_detailed.png",
            run=config_provider("run", "name"),
        ),
        exported_variables=expand(
            RESULTS + "ariadne/exported_variables_full.xlsx",
            run=config_provider("run", "name"),
        ),
    script:
        "scripts/plot_ariadne_scenario_comparison.py"


rule build_scenarios:
    params:
        scenarios=config_provider("run", "name"),
        db_name=config_provider("iiasa_database", "db_name"),
        leitmodelle=config_provider("iiasa_database", "leitmodelle"),
    input:
        ariadne_database=resources("ariadne_database.csv"),
        scenario_yaml=config["run"]["scenarios"]["manual_file"],
    output:
        scenario_yaml=config["run"]["scenarios"]["file"],
    log:
        "logs/build_scenarios.log",
    script:
        "scripts/build_scenarios.py"


rule check_sector_ratios:
    input:
        network=RESULTS
        + "postnetworks/base_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    log:
        "logs/check_sector_ratios.log",
    script:
        "scripts/check_sector_ratios.py"
