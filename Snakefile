# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from pathlib import Path
import yaml
import sys
from os.path import normpath, exists, join
from shutil import copyfile, move, rmtree, unpack_archive
from snakemake.utils import min_version

min_version("8.11")

from scripts._helpers import (
    get_rdir,
    get_scenarios,
    get_shadow,
    path_provider,
)


configfile: "config/config.default.yaml"
configfile: "config/plotting.default.yaml"
configfile: "config/config.de.yaml"


run = config["run"]
scenarios = get_scenarios(run)
RDIR = get_rdir(run)
shadow_config = get_shadow(run)

policy = run["shared_resources"]["policy"]
exclude = run["shared_resources"]["exclude"]

shared_resources = run["shared_resources"]["policy"]
exclude_from_shared = run["shared_resources"]["exclude"]
logs = path_provider("logs/", RDIR, shared_resources, exclude_from_shared)
benchmarks = path_provider("benchmarks/", RDIR, shared_resources, exclude_from_shared)
resources = path_provider("resources/", RDIR, shared_resources, exclude_from_shared)

RESULTS = "results/" + RDIR


localrules:
    purge,


wildcard_constraints:
    clusters="[0-9]+(m|c)?|all|adm",
    opts=r"[-+a-zA-Z0-9\.]*",
    sector_opts=r"[-+a-zA-Z0-9\.\s]*",
    planning_horizons=r"[0-9]{4}",


include: "rules/common.smk"


# Data constants
OSM_DATASET = dataset_version("osm")


include: "rules/collect.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/solve_electricity.smk"
include: "rules/postprocess.smk"
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
        expand(resources("maps/power-network.pdf"), run=config["run"]["name"]),
        expand(
            resources("maps/power-network-s-{clusters}.pdf"),
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        # COP profiles plots
        expand(
            RESULTS + "graphs/cop_profiles_s_{clusters}_{planning_horizons}.html",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf"
                if config_provider("sector", "H2_network")(w)
                else []
            ),
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf"
                if config_provider("sector", "gas_network")(w)
                else []
            ),
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: expand(
            (
                RESULTS + "csvs/cumulative_costs.csv"
                if config_provider("foresight")(w) == "myopic"
                else []
            ),
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "graphics/balance_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "graphics/heatmap_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        # Explicitly list heat source types for temperature maps
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_temperature_map_river_water.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "river_water"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_temperature_map_sea_water.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "sea_water"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_temperature_map_ambient_air.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "air"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        # Only river_water has energy maps
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_energy_map_river_water.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "river_water"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "graphics/balance_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "graphics/heatmap_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "graphics/interactive_bus_balance/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: balance_map_paths("static", w),
        lambda w: balance_map_paths("interactive", w),
    default_target: True


rule create_scenarios:
    output:
        config["run"]["scenarios"]["file"],
    script:
        "config/create_scenarios.py"


rule purge:
    run:
        import builtins

        do_purge = builtins.input(
            "Do you really want to delete all generated files?\n"
            "\t* resources\n"
            "\t* results\n"
            "\t* docs\n"
            "Downloaded files are kept.\n"
            "Delete all files in the folders above? [y/N] "
        )
        if do_purge == "y":

            # Remove the directories and recreate them with .gitkeep
            for dir_path in ["resources/", "results/"]:
                rmtree(dir_path, ignore_errors=True)
                Path(dir_path).mkdir(parents=True, exist_ok=True)
                (Path(dir_path) / ".gitkeep").touch()

            rmtree("doc/_build", ignore_errors=True)
            print(
                "Purging all generated resources, results and docs. Downloads are kept."
            )
        else:
            raise Exception(f"Input {do_purge}. Aborting purge.")


rule dump_graph_config:
    """Dump the current Snakemake configuration to a YAML file for graph generation."""
    output:
        config_file=temp(resources("dag_final_config.yaml")),
    run:
        import yaml

        with open(output.config_file, "w") as f:
            yaml.dump(config, f)


rule rulegraph:
    """Generates Rule DAG in DOT, PDF, PNG, and SVG formats using the final configuration."""
    message:
        "Creating RULEGRAPH dag in multiple formats using the final configuration."
    input:
        config_file=rules.dump_graph_config.output.config_file,
    output:
        dot=resources("dag_rulegraph.dot"),
        pdf=resources("dag_rulegraph.pdf"),
        png=resources("dag_rulegraph.png"),
        svg=resources("dag_rulegraph.svg"),
    shell:
        r"""
        # Generate DOT file using nested snakemake with the dumped final config
        echo "[Rule rulegraph] Using final config file: {input.config_file}"
        snakemake --rulegraph --configfile {input.config_file} --quiet | sed -n "/digraph/,\$p" > {output.dot}

        # Generate visualizations from the DOT file
        if [ -s {output.dot} ]; then
            dot -c

            echo "[Rule rulegraph] Generating PDF from DOT"
            dot -Tpdf -o {output.pdf} {output.dot} || {{ echo "Error: Failed to generate PDF. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule rulegraph] Generating PNG from DOT"
            dot -Tpng -o {output.png} {output.dot} || {{ echo "Error: Failed to generate PNG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule rulegraph] Generating SVG from DOT"
            dot -Tsvg -o {output.svg} {output.dot} || {{ echo "Error: Failed to generate SVG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule rulegraph] Successfully generated all formats."
        else
            echo "[Rule rulegraph] Error: Failed to generate valid DOT content." >&2
            exit 1
        fi
        """


rule filegraph:
    """Generates File DAG in DOT, PDF, PNG, and SVG formats using the final configuration."""
    message:
        "Creating FILEGRAPH dag in multiple formats using the final configuration."
    input:
        config_file=rules.dump_graph_config.output.config_file,
    output:
        dot=resources("dag_filegraph.dot"),
        pdf=resources("dag_filegraph.pdf"),
        png=resources("dag_filegraph.png"),
        svg=resources("dag_filegraph.svg"),
    shell:
        r"""
        # Generate DOT file using nested snakemake with the dumped final config
        echo "[Rule filegraph] Using final config file: {input.config_file}"
        snakemake --filegraph all --configfile {input.config_file} --quiet | sed -n "/digraph/,\$p" > {output.dot}

        # Generate visualizations from the DOT file
        if [ -s {output.dot} ]; then
            echo "[Rule filegraph] Generating PDF from DOT"
            dot -Tpdf -o {output.pdf} {output.dot} || {{ echo "Error: Failed to generate PDF. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule filegraph] Generating PNG from DOT"
            dot -Tpng -o {output.png} {output.dot} || {{ echo "Error: Failed to generate PNG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule filegraph] Generating SVG from DOT"
            dot -Tsvg -o {output.svg} {output.dot} || {{ echo "Error: Failed to generate SVG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule filegraph] Successfully generated all formats."
        else
            echo "[Rule filegraph] Error: Failed to generate valid DOT content." >&2
            exit 1
        fi
        """


rule doc:
    message:
        "Build documentation."
    output:
        directory("doc/_build"),
    shell:
        "pixi run build-docs {output} html"


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


rule sync_dry:
    params:
        cluster=f"{config['remote']['ssh']}:{config['remote']['path']}",
    shell:
        """
        rsync -uvarh --ignore-missing-args --files-from=.sync-send . {params.cluster} -n
        rsync -uvarh --no-g {params.cluster}/resources . -n || echo "No resources directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/results . -n || echo "No results directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/logs . -n || echo "No logs directory, skipping rsync"
        """


rule clean:
    message:
        "Remove all build results but keep downloaded data."
    run:
        import shutil

        shutil.rmtree("resources")
        shutil.rmtree("results")
        print("Data downloaded to data/ has not been cleaned.")


if (ARIADNE_DATABASE := dataset_version("ariadne_database"))["source"] in ["primary"]:

    rule retrieve_ariadne_database:
        params:
            source="primary",
        output:
            data="data/ariadne_database.csv",
        log:
            "logs/retrieve_ariadne_database_primary.log",
        resources:
            mem_mb=1000,
        script:
            "scripts/pypsa-de/retrieve_ariadne_database.py"


if (ARIADNE_DATABASE := dataset_version("ariadne_database"))["source"] in ["archive"]:

    rule retrieve_ariadne_database:
        params:
            source="archive",
        input:
            raw_xlsx=storage(ARIADNE_DATABASE["url"]),
        output:
            data="data/ariadne_database.csv",
        log:
            "logs/retrieve_ariadne_database_archive.log",
        resources:
            mem_mb=1000,
        script:
            "scripts/pypsa-de/retrieve_ariadne_database.py"


if (ARIADNE_TEMPLATE := dataset_version("ariadne_template"))["source"] in [
    "primary",
]:

    rule retrieve_ariadne_template:
        input:
            storage(
                "https://github.com/iiasa/ariadne-intern-workflow/raw/main/attachments/2025-01-27_template_Ariadne.xlsx",
            ),
        output:
            "data/template_ariadne_database.xlsx",
        run:
            move(input[0], output[0])


if (OPEN_MASTR := dataset_version("open_mastr"))["source"] in ["primary"]:

    rule retrieve_open_mastr:
        input:
            storage(OPEN_MASTR["url"]),
        params:
            "data/mastr",
        output:
            "data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
            "data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
        run:
            unpack_archive(input[0], params[0])


if (EGON := dataset_version("egon"))["source"] in ["primary"]:

    rule retrieve_egon_data:
        input:
            spatial=storage(
                f"{EGON['url']}?id_spatial=5&year=2018",
            ),
            mapping=storage(
                f"{EGON['url']}_description?id_spatial=5",
            ),
        output:
            spatial="data/egon/demandregio_spatial_2018.json",
            mapping="data/egon/mapping_technologies.json",
        run:
            move(input.spatial, output.spatial)
            move(input.mapping, output.mapping)


rule build_exogenous_mobility_data:
    params:
        reference_scenario=config_provider("pypsa-de", "reference_scenario"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        leitmodelle=config_provider("pypsa-de", "leitmodelle"),
        ageb_for_mobility=config_provider("pypsa-de", "ageb_for_mobility"),
        uba_for_mobility=config_provider("pypsa-de", "uba_for_mobility"),
        shipping_oil_share=config_provider("sector", "shipping_oil_share"),
        aviation_demand_factor=config_provider("sector", "aviation_demand_factor"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        ariadne="data/ariadne_database.csv",
        energy_totals=resources("energy_totals.csv"),
    output:
        mobility_data=resources(
            "modified_mobility_data_{clusters}_{planning_horizons}.csv"
        ),
    resources:
        mem_mb=1000,
    log:
        logs("build_exogenous_mobility_data_{clusters}_{planning_horizons}.log"),
    script:
        "scripts/pypsa-de/build_exogenous_mobility_data.py"


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


rule prepare_district_heating_subnodes:
    params:
        district_heating=config_provider("sector", "district_heating"),
        baseyear=config_provider("scenario", "planning_horizons", 0),
    input:
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        fernwaermeatlas="data/fernwaermeatlas/fernwaermeatlas.xlsx",
        cities="data/fernwaermeatlas/cities_geolocations.geojson",
        lau_regions=rules.retrieve_lau_regions.output["zip"],
        census=storage(
            "https://www.destatis.de/static/DE/zensus/gitterdaten/Zensus2022_Heizungsart.zip",
            keep_local=True,
        ),
        osm_land_cover=storage(
            "https://heidata.uni-heidelberg.de/api/access/datafile/23053?format=original&gbrecs=true",
            keep_local=True,
        ),
        natura=ancient("data/bundle/natura/natura.tiff"),
        groundwater_depth=storage(
            "http://thredds-gfnl.usc.es/thredds/fileServer/GLOBALWTDFTP/annualmeans/EURASIA_WTD_annualmean.nc",
            keep_local=True,
        ),
    output:
        district_heating_subnodes=resources(
            "district_heating_subnodes_base_s_{clusters}.geojson"
        ),
        regions_onshore_extended=resources(
            "regions_onshore_base-extended_s_{clusters}.geojson"
        ),
        regions_onshore_restricted=resources(
            "regions_onshore_base-restricted_s_{clusters}.geojson"
        ),
    resources:
        mem_mb=20000,
    script:
        "scripts/pypsa-de/prepare_district_heating_subnodes.py"


def baseyear_value(wildcards):
    return config_provider("scenario", "planning_horizons", 0)(wildcards)


rule add_district_heating_subnodes:
    params:
        district_heating=config_provider("sector", "district_heating"),
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        heat_pump_sources=config_provider(
            "sector", "heat_pump_sources", "urban central"
        ),
        heat_utilisation_potentials=config_provider(
            "sector", "district_heating", "heat_utilisation_potentials"
        ),
        direct_utilisation_heat_sources=config_provider(
            "sector", "district_heating", "direct_utilisation_heat_sources"
        ),
        adjustments=config_provider("adjustments", "sector"),
    input:
        unpack(input_heat_source_power),
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        subnodes=resources("district_heating_subnodes_base_s_{clusters}.geojson"),
        nuts3=resources("nuts3_shapes.geojson"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        fernwaermeatlas="data/fernwaermeatlas/fernwaermeatlas.xlsx",
        cities="data/fernwaermeatlas/cities_geolocations.geojson",
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        direct_heat_source_utilisation_profiles=resources(
            "direct_heat_source_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        existing_heating_distribution=lambda w: resources(
            f"existing_heating_distribution_base_s_{{clusters}}_{baseyear_value(w)}.csv"
        ),
        lau_regions=rules.retrieve_lau_regions.output["zip"],
    output:
        network=resources(
            "networks/base-extended_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        district_heating_subnodes=resources(
            "district_heating_subnodes_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.geojson"
        ),
        existing_heating_distribution_extended=(
            resources(
                "existing_heating_distribution_base-extended_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv"
            )
            if baseyear_value != "{planning_horizons}"
            else []
        ),
    resources:
        mem_mb=10000,
    script:
        "scripts/pypsa-de/add_district_heating_subnodes.py"


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
        H2_plants=config_provider("electricity", "H2_plants"),
        onshore_nep_force=config_provider("onshore_nep_force"),
        offshore_nep_force=config_provider("offshore_nep_force"),
        shipping_methanol_efficiency=config_provider(
            "sector", "shipping_methanol_efficiency"
        ),
        shipping_oil_efficiency=config_provider("sector", "shipping_oil_efficiency"),
        shipping_methanol_share=config_provider("sector", "shipping_methanol_share"),
        mwh_meoh_per_tco2=config_provider("sector", "MWh_MeOH_per_tCO2"),
        scale_capacity=config_provider("scale_capacity"),
        bev_charge_rate=config_provider("sector", "bev_charge_rate"),
        bev_energy=config_provider("sector", "bev_energy"),
        bev_dsm_availability=config_provider("sector", "bev_dsm_availability"),
        uba_for_industry=config_provider("pypsa-de", "uba_for_industry", "enable"),
        scale_industry_non_energy=config_provider(
            "pypsa-de", "uba_for_industry", "scale_industry_non_energy"
        ),
        limit_cross_border_flows_ac=config_provider(
            "pypsa-de", "limit_cross_border_flows_ac"
        ),
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
        wkn=lambda w: (
            resources("wasserstoff_kernnetz_base_s_{clusters}.csv")
            if config_provider("wasserstoff_kernnetz", "enable")(w)
            else []
        ),
        costs=resources("costs_{planning_horizons}_processed.csv"),
        modified_mobility_data=resources(
            "modified_mobility_data_{clusters}_{planning_horizons}.csv"
        ),
        biomass_potentials=resources(
            "biomass_potentials_s_{clusters}_{planning_horizons}.csv"
        ),
        industrial_demand=resources(
            "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
        ),
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}-modified.csv"
        ),
        industry_sector_ratios=resources(
            "industry_sector_ratios_{planning_horizons}.csv"
        ),
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        shipping_demand=resources("shipping_demand_s_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        offshore_connection_points="data/pypsa-de/offshore_connection_points.csv",
        new_industrial_energy_demand="data/pypsa-de/UBA_Projektionsbericht2025_Abbildung31_MWMS.csv",
    output:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_final.nc"
        ),
    resources:
        mem_mb=4000,
    log:
        RESULTS
        + "logs/modify_prenetwork_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    script:
        "scripts/pypsa-de/modify_prenetwork.py"


ruleorder: modify_industry_production > build_industrial_production_per_country_tomorrow


rule modify_existing_heating:
    input:
        ariadne="data/ariadne_database.csv",
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
    output:
        existing_heating=resources("existing_heating.csv"),
    resources:
        mem_mb=1000,
    log:
        logs("modify_existing_heating.log"),
    script:
        "scripts/pypsa-de/modify_existing_heating.py"


rule build_existing_chp_de:
    params:
        district_heating_subnodes=config_provider(
            "sector", "district_heating", "subnodes"
        ),
    input:
        mastr_biomass="data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
        mastr_combustion="data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
        plz_mapping=storage(
            "https://raw.githubusercontent.com/WZBSocialScienceCenter/plz_geocoord/master/plz_geocoord.csv",
            keep_local=True,
        ),
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        district_heating_subnodes=lambda w: (
            resources("district_heating_subnodes_base_s_{clusters}.geojson")
            if config_provider("sector", "district_heating", "subnodes", "enable")(w)
            else []
        ),
    output:
        german_chp=resources("german_chp_base_s_{clusters}.csv"),
    resources:
        mem_mb=4000,
    log:
        logs("build_existing_chp_de_{clusters}.log"),
    script:
        "scripts/pypsa-de/build_existing_chp_de.py"


rule modify_industry_production:
    params:
        reference_scenario=config_provider("pypsa-de", "reference_scenario"),
    input:
        ariadne="data/ariadne_database.csv",
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
        logs("modify_industry_production_{planning_horizons}.log"),
    script:
        "scripts/pypsa-de/modify_industry_production.py"


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
        locations="data/pypsa-de/wasserstoff_kernnetz/locations_wasserstoff_kernnetz.csv",
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


rule export_ariadne_variables:
    params:
        planning_horizons=config_provider("scenario", "planning_horizons"),
        hours=config_provider("clustering", "temporal", "resolution_sector"),
        config_industry=config_provider("industry"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        co2_sequestration_cost=config_provider("sector", "co2_sequestration_cost"),
        post_discretization=config_provider("solving", "options", "post_discretization"),
        NEP_year=lambda w: config_provider("costs", "custom_cost_fn")(w)[-8:-4],
        NEP_transmission=config_provider("costs", "transmission"),
    input:
        template="data/template_ariadne_database.xlsx",
        industry_demands=expand(
            resources(
                "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
            ),
            **config["scenario"],
            allow_missing=True,
        ),
        networks=expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            allow_missing=True,
        ),
        costs=expand(
            resources("costs_{planning_horizons}_processed.csv"),
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
        mem_mb=16000,
    log:
        RESULTS + "logs/export_ariadne_variables.log",
    script:
        "scripts/pypsa-de/export_ariadne_variables.py"


rule plot_ariadne_variables:
    params:
        reference_scenario=config_provider("pypsa-de", "reference_scenario"),
    input:
        exported_variables_full=RESULTS + "ariadne/exported_variables_full.xlsx",
        ariadne_database="data/ariadne_database.csv",
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
        # expand(
        #     RESULTS + "ariadne/capacity_detailed.png",
        #     run=config_provider("run", "name"),
        # ),
        expand(
            RESULTS
            + "maps/base_s_{clusters}_{opts}_{sector_opts}-h2_network_incl_kernnetz_{planning_horizons}.pdf",
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
        scenarios=config["run"]["name"],
        leitmodelle=config["pypsa-de"]["leitmodelle"],
    input:
        ariadne_database="data/ariadne_database.csv",
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
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        map=RESULTS
        + "maps/base_s_{clusters}_{opts}_{sector_opts}-h2_network_incl_kernnetz_{planning_horizons}.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        RESULTS
        + "logs/plot_hydrogen_network_incl_kernnetz/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_hydrogen_network_incl_kernnetz/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "scripts/pypsa-de/plot_hydrogen_network_incl_kernnetz.py"


rule plot_ariadne_report:
    params:
        planning_horizons=config_provider("scenario", "planning_horizons"),
        plotting=config_provider("plotting"),
        run=config_provider("run", "name"),
        foresight=config_provider("foresight"),
        post_discretization=config_provider("solving", "options", "post_discretization"),
        NEP_year=lambda w: config_provider("costs", "custom_cost_fn")(w)[-8:-4],
        hours=config_provider("clustering", "temporal", "resolution_sector"),
        NEP_transmission=config_provider("costs", "transmission"),
    input:
        networks=expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
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
            resources("costs_{planning_horizons}_processed.csv"),
            **config["scenario"],
            allow_missing=True,
        ),
        exported_variables_full=RESULTS + "ariadne/exported_variables_full.xlsx",
    output:
        elec_price_duration_curve=RESULTS
        + "ariadne/report/elec_price_duration_curve.pdf",
        elec_price_duration_hist=RESULTS + "ariadne/report/elec_price_duration_hist.pdf",
        backup_capacity=RESULTS + "ariadne/report/backup_capacity.pdf",
        backup_generation=RESULTS + "ariadne/report/backup_generation.pdf",
        results=directory(RESULTS + "ariadne/report"),
        elec_transmission=directory(RESULTS + "ariadne/report/elec_transmission"),
        h2_transmission=directory(RESULTS + "ariadne/report/h2_transmission"),
        co2_transmission=directory(RESULTS + "ariadne/report/co2_transmission"),
        elec_balances=directory(RESULTS + "ariadne/report/elec_balance_timeseries"),
        heat_balances=directory(RESULTS + "ariadne/report/heat_balance_timeseries"),
        nodal_balances=directory(RESULTS + "ariadne/report/balance_timeseries_2045"),
    resources:
        mem_mb=32000,
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
