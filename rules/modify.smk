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


rule modify_population_layouts:
    input:
        pop_layout_total=resources("pop_layout_total-raw.nc"),
        pop_layout_urban=resources("pop_layout_urban-raw.nc"),
        pop_layout_rural=resources("pop_layout_rural-raw.nc"),
    output:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
    resources:
        mem_mb=2000,
    log:
        logs("modify_population_layouts.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/pypsa-at/modify_population_layouts.py"


rule export_iamc_variables:
    input:
        networks=expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            allow_missing=True,
        ),
    output:
        exported_variables=RESULTS + "evaluation/exported_iamc_variables.xlsx",
    resources:
        mem_mb=16000,
    log:
        RESULTS + "logs/export_iamc_variables.log",
    script:
        "../scripts/pypsa-at/export_iamc_variables.py"


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
        "../scripts/pypsa-at/plot_iamc_variables.py"


rule validate_pypsa_at:  # rule to run integration tests on solved networks
    params:
        clustering=config_provider("clustering"),
        rdir=RESULTS,
    input:
        expand(
            RESULTS + "evaluation/HTML/sankey_diagram_EU_2050.html",
            run=config["run"]["name"],
        ),
    output:
        validity_report=RESULTS + "validity_report.html",
    resources:
        mem_mb=16000,
    shell:
        # setting RESULTS here does not work, but why?!
        # "pytest test/test_base_network.py --html {params.rdir}/validity_report.html --results_dir {params.rdir}"
        "pytest --html {params.rdir}/validity_report.html --results_dir {params.rdir}"
