rule modify_nuts3_shapes:
    params:
        clustering=config_provider("clustering", "mode"),
        admin_levels=config_provider("clustering", "administrative"),
    input:
        nuts3_shapes=resources("nuts3_shapes-raw.geojson"),
    output:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
    log:
        logs("modify_nuts3_shapes.log"),
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/pypsa-at/modify_nuts3_shapes.py"


rule export_iamc_variables:
    params:
        planning_horizons=config_provider("scenario", "planning_horizons"),
        hours=config_provider("clustering", "temporal", "resolution_sector"),
        max_hours=config_provider("electricity", "max_hours"),
        costs=config_provider("costs"),
        config_industry=config_provider("industry"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        co2_price_add_on_fossils=config_provider("co2_price_add_on_fossils"),
        co2_sequestration_cost=config_provider("sector", "co2_sequestration_cost"),
        post_discretization=config_provider("solving", "options", "post_discretization"),
        NEP_year=config_provider("costs", "NEP"),
        NEP_transmission=config_provider("costs", "transmission"),
    input:
        results_path=RESULTS,
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
        exported_variables=RESULTS + "evaluation/exported_iamc_variables.xlsx",
        exported_variables_full=RESULTS + "evaluation/exported_iamc_variables_full.xlsx",
    resources:
        mem_mb=16000,
    log:
        RESULTS + "logs/export_iamc_variables.log",
    script:
        "scripts/pypsa-at/export_iamc_variables.py"


rule plot_iamc_variables:
    input:
        iamc_variables=RESULTS + "evaluation/exported_iamc_variables.xlsx",
    output:
        sankey=RESULTS + "evaluation/HTML/sankey_diagram_EU_2050.html",
    resources:
        mem_mb=16000,
    log:
        RESULTS + "logs/plot_iamc_variables.log",
    script:
        "scripts/pypsa-at/plot_iamc_variables.py"


rule validate_pypsa_at:
    params:
        clsutering=config_provider("clustering"),
    input:
        results_path=RESULTS,
    output:
        validity_report=RESULTS + "validity_report.html",
    resources:
        mem_mb=16000,
    log:
        RESULTS + "logs/validity_report.log",
    script:
        "scripts/pypsa-at/validity_report.py"
