# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from pathlib import Path
import yaml
import sys
from os.path import normpath, exists
from shutil import copyfile, move, rmtree, unpack_archive
from snakemake.utils import min_version

min_version("8.11")

from scripts._helpers import path_provider, copy_default_files, get_scenarios, get_rdir


copy_default_files(workflow)


configfile: "config/config.default.yaml"
configfile: "config/config.yaml"
configfile: "config/config.personal.yaml"


run = config["run"]
scenarios = get_scenarios(run)
RDIR = get_rdir(run)

policy = run["shared_resources"]["policy"]
exclude = run["shared_resources"]["exclude"]

shared_resources = run["shared_resources"]["policy"]
exclude_from_shared = run["shared_resources"]["exclude"]
logs = path_provider("logs/", RDIR, shared_resources, exclude_from_shared)
benchmarks = path_provider("benchmarks/", RDIR, shared_resources, exclude_from_shared)
resources = path_provider("resources/", RDIR, shared_resources, exclude_from_shared)

# logs = path_provider("logs/", RDIR, policy, exclude)
# benchmarks = path_provider("benc€€hmarks/", RDIR, policy, exclude)
# resources = path_provider("resources/", RDIR, policy, exclude)


CDIR = "" if run["shared_cutouts"] else RDIR
RESULTS = "results/" + RDIR


localrules:
    purge,


wildcard_constraints:
    clusters="[0-9]+(m|c)?|all",
    ll=r"(v|c)([0-9\.]+|opt)",
    opts=r"[-+a-zA-Z0-9\.]*",
    sector_opts=r"[-+a-zA-Z0-9\.\s]*",


include: "rules/common.smk"
include: "rules/collect.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/solve_electricity.smk"
include: "rules/postprocess.smk"
include: "rules/validate.smk"
include: "rules/development.smk"


if config["foresight"] == "overnight":

    include: "rules/solve_overnight.smk"


if config["foresight"] == "myopic":

    include: "rules/solve_myopic.smk"


if config["foresight"] == "perfect":

    include: "rules/solve_perfect.smk"


rule all:
    input:
        expand(RESULTS + "graphs/costs.svg", run=config["run"]["name"]),
    default_target: True


rule create_scenarios:
    output:
        config["run"]["scenarios"]["file"],
    conda:
        "envs/retrieve.yaml"
    script:
        "config/create_scenarios.py"


rule purge:
    run:
        import builtins

        do_purge = builtins.input(
            "Do you really want to delete all generated resources, \nresults and docs (downloads are kept)? [y/N] "
        )
        if do_purge == "y":
            rmtree("resources/", ignore_errors=True)
            rmtree("results/", ignore_errors=True)
            rmtree("doc/_build", ignore_errors=True)
            print("Purging generated resources, results and docs. Downloads are kept.")
        else:
            raise Exception(f"Input {do_purge}. Aborting purge.")


rule dag:
    message:
        "Creating DAG of workflow."
    output:
        dot=resources("dag.dot"),
        pdf=resources("dag.pdf"),
        png=resources("dag.png"),
    conda:
        "envs/environment.yaml"
    shell:
        r"""
        snakemake --rulegraph all | sed -n "/digraph/,\$p" > {output.dot}
        dot -Tpdf -o {output.pdf} {output.dot}
        dot -Tpng -o {output.png} {output.dot}
        """


rule doc:
    message:
        "Build documentation."
    output:
        directory("doc/_build"),
    shell:
        "make -C doc html"


rule sync:
    params:
        cluster=f"{config['remote']['ssh']}:{config['remote']['path']}",
    shell:
        """
        rsync -uvarh --ignore-missing-args --files-from=.sync-send . {params.cluster}
        rsync -uvarh --no-g {params.cluster}/resources . || echo "No resources directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/results . || echo "No results directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/logs . || echo "No logs directory, skipping rsync"
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
        "logs/pypsa-de/retrieve_ariadne_database.log",
    resources:
        mem_mb=1000,
    script:
        "scripts/pypsa-de/retrieve_ariadne_database.py"


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
    log:
        logs("modify_cost_data_{planning_horizons}.log"),
    script:
        "scripts/pypsa-de/modify_cost_data.py"


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
        "scripts/pypsa-de/build_mobility_demand.py"


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
    log:
        logs("build_egon_data.log"),
    script:
        "scripts/pypsa-de/build_egon_data.py"


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
    log:
        logs("modify_district_heat_share_{clusters}_{planning_horizons}.log"),
    script:
        "scripts/pypsa-de/modify_district_heat_share.py"


rule modify_prenetwork:
    params:
        efuel_export_ban=config_provider("solving", "constraints", "efuel_export_ban"),
        enable_kernnetz=config_provider("wasserstoff_kernnetz", "enable"),
        costs=config_provider("costs"),
        max_hours=config_provider("electricity", "max_hours"),
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
        must_run=config_provider("must_run"),
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
        scale_capacity=config_provider("scale_capacity"),
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
        mem_mb=4000,
    log:
        RESULTS
        + "logs/modify_prenetwork_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    script:
        "scripts/pypsa-de/modify_prenetwork.py"


ruleorder: modify_industry_demand > build_industrial_production_per_country_tomorrow


rule modify_existing_heating:
    params:
        iiasa_reference_scenario=config_provider("iiasa_database", "reference_scenario"),
        leitmodelle=config_provider("iiasa_database", "leitmodelle"),
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
    log:
        logs("modify_existing_heating.log"),
    script:
        "scripts/pypsa-de/modify_existing_heating.py"


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
    input:
        mastr_biomass="data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
        mastr_combustion="data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
        plz_mapping=storage(
            "https://raw.githubusercontent.com/WZBSocialScienceCenter/plz_geocoord/master/plz_geocoord.csv",
            keep_local=True,
        ),
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        german_chp=resources("german_chp_{clusters}.csv"),
    log:
        logs("build_existing_chp_de_{clusters}.log"),
    script:
        "scripts/pypsa-de/build_existing_chp_de.py"


rule modify_industry_demand:
    params:
        db_name=config_provider("iiasa_database", "db_name"),
        reference_scenario=config_provider("iiasa_database", "reference_scenario"),
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
    log:
        logs("modify_industry_demand_{planning_horizons}.log"),
    script:
        "scripts/pypsa-de/modify_industry_demand.py"


rule build_wasserstoff_kernnetz:
    params:
        kernnetz=config_provider("wasserstoff_kernnetz"),
    input:
        wasserstoff_kernnetz_1=storage(
            "https://fnb-gas.de/wp-content/uploads/2024/07/2024_07_22_Anlage2_Leitungsmeldungen_weiterer_potenzieller_Wasserstoffnetzbetreiber.xlsx",
            keep_local=True,
        ),
        wasserstoff_kernnetz_2=storage(
            "https://fnb-gas.de/wp-content/uploads/2024/07/2024_07_22_Anlage3_FNB_Massnahmenliste_Neubau.xlsx",
            keep_local=True,
        ),
        wasserstoff_kernnetz_3=storage(
            "https://fnb-gas.de/wp-content/uploads/2024/07/2024_07_22_Anlage4_FNB_Massnahmenliste_Umstellung.xlsx",
            keep_local=True,
        ),
        gadm=storage(
            "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_DEU_1.json.zip",
            keep_local=True,
        ),
        locations="ariadne-data/wasserstoff_kernnetz/locations_wasserstoff_kernnetz.csv",
        regions_onshore=resources("regions_onshore_base_s.geojson"),
        regions_offshore=resources("regions_offshore_base_s.geojson"),
    output:
        cleaned_wasserstoff_kernnetz=resources("wasserstoff_kernnetz.csv"),
    log:
        logs("build_wasserstoff_kernnetz.log"),
    script:
        "scripts/pypsa-de/build_wasserstoff_kernnetz.py"


rule cluster_wasserstoff_kernnetz:
    params:
        kernnetz=config_provider("wasserstoff_kernnetz"),
    input:
        cleaned_h2_network=resources("wasserstoff_kernnetz.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
    output:
        clustered_h2_network=resources("wasserstoff_kernnetz_base_s_{clusters}.csv"),
    log:
        logs("cluster_wasserstoff_kernnetz_{clusters}.log"),
    script:
        "scripts/pypsa-de/cluster_wasserstoff_kernnetz.py"


rule download_ariadne_template:
    input:
        storage(
            "https://github.com/iiasa/ariadne-intern-workflow/raw/main/attachments/2025-01-27_template_Ariadne.xlsx",
            keep_local=True,
        ),
    output:
        resources("template_ariadne_database.xlsx"),
    run:
        move(input[0], output[0])


rule export_ariadne_variables:
    params:
        planning_horizons=config_provider("scenario", "planning_horizons"),
        hours=config_provider("clustering", "temporal", "resolution_sector"),
        costs=config_provider("costs"),
        config_industry=config_provider("industry"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        co2_price_add_on_fossils=config_provider("co2_price_add_on_fossils"),
        co2_sequestration_cost=config_provider("sector", "co2_sequestration_cost"),
        post_discretization=config_provider("solving", "options", "post_discretization"),
        NEP_year=config_provider("costs", "NEP"),
        NEP_transmission=config_provider("costs", "transmission"),
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
            + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
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
    resources:
        # A heuristic to assign more memory for larger networks. Should probably be quadratic??
        mem_mb=(lambda w: 400 * int(config_provider("scenario", "clusters")(w)[0])),
    log:
        RESULTS + "logs/export_ariadne_variables.log",
    script:
        "scripts/pypsa-de/export_ariadne_variables.py"


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
        NEP_plot=RESULTS + "ariadne/NEP_plot.png",
        NEP_Trassen_plot=RESULTS + "ariadne/NEP_Trassen_plot.png",
        transmission_investment_csv=RESULTS + "ariadne/transmission_investment.csv",
        trassenlaenge_csv=RESULTS + "ariadne/trassenlaenge.csv",
        Kernnetz_Investment_plot=RESULTS + "ariadne/Kernnetz_Investment_plot.png",
        elec_trade=RESULTS + "ariadne/elec-trade-DE.pdf",
        h2_trade=RESULTS + "ariadne/h2-trade-DE.pdf",
        trade_balance=RESULTS + "ariadne/trade-balance-DE.pdf",
    log:
        RESULTS + "logs/plot_ariadne_variables.log",
    script:
        "scripts/pypsa-de/plot_ariadne_variables.py"


rule ariadne_all:
    input:
        expand(RESULTS + "graphs/costs.svg", run=config_provider("run", "name")),
        expand(
            RESULTS + "ariadne/capacity_detailed.png",
            run=config_provider("run", "name"),
        ),
        expand(
            RESULTS
            + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_incl_kernnetz_{planning_horizons}.pdf",
            run=config_provider("run", "name"),
            **config["scenario"],
            allow_missing=True,
        ),
        exported_variables=expand(
            RESULTS + "ariadne/exported_variables_full.xlsx",
            run=config_provider("run", "name"),
        ),
    script:
        "scripts/pypsa-de/plot_ariadne_scenario_comparison.py"


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
        "scripts/pypsa-de/build_scenarios.py"


rule plot_hydrogen_network_incl_kernnetz:
    params:
        plotting=config_provider("plotting"),
        foresight=config_provider("foresight"),
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        map=RESULTS
        + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_incl_kernnetz_{planning_horizons}.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        RESULTS
        + "logs/plot_hydrogen_network_incl_kernnetz/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_hydrogen_network_incl_kernnetz/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "scripts/pypsa-de/plot_hydrogen_network_incl_kernnetz.py"


rule plot_ariadne_report:
    params:
        planning_horizons=config_provider("scenario", "planning_horizons"),
        plotting=config_provider("plotting"),
        run=config_provider("run", "name"),
        foresight=config_provider("foresight"),
        costs=config_provider("costs"),
        post_discretization=config_provider("solving", "options", "post_discretization"),
        NEP_year=config_provider("costs", "NEP"),
        hours=config_provider("clustering", "temporal", "resolution_sector"),
        NEP_transmission=config_provider("costs", "transmission"),
    input:
        networks=expand(
            RESULTS
            + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            allow_missing=True,
        ),
        regions_onshore_clustered=expand(
            resources("regions_onshore_base_s_{clusters}.geojson"),
            clusters=config["scenario"]["clusters"],
            allow_missing=True,
        ),
        rc="matplotlibrc",
        costs=expand(
            resources("costs_{planning_horizons}.csv"),
            **config["scenario"],
            allow_missing=True,
        ),
    output:
        elec_price_duration_curve=RESULTS
        + "ariadne/report/elec_price_duration_curve.pdf",
        elec_price_duration_hist=RESULTS + "ariadne/report/elec_price_duration_hist.pdf",
        backup_capacity=RESULTS + "ariadne/report/backup_capacity.pdf",
        backup_generation=RESULTS + "ariadne/report/backup_generation.pdf",
        elec_prices_spatial_de=RESULTS + "ariadne/report/elec_prices_spatial_de.pdf",
        results=directory(RESULTS + "ariadne/report"),
        elec_transmission=directory(RESULTS + "ariadne/report/elec_transmission"),
        h2_transmission=directory(RESULTS + "ariadne/report/h2_transmission"),
        co2_transmission=directory(RESULTS + "ariadne/report/co2_transmission"),
        elec_balances=directory(RESULTS + "ariadne/report/elec_balance_timeseries"),
        heat_balances=directory(RESULTS + "ariadne/report/heat_balance_timeseries"),
        nodal_balances=directory(RESULTS + "ariadne/report/balance_timeseries_2045"),
    resources:
        mem_mb=30000,
    log:
        RESULTS + "logs/plot_ariadne_report.log",
    script:
        "scripts/pypsa-de/plot_ariadne_report.py"


rule ariadne_report_only:
    input:
        expand(
            RESULTS + "ariadne/report/elec_price_duration_curve.pdf",
            run=config_provider("run", "name"),
        ),
