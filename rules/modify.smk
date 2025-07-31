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
        exported_variables=RESULTS + "evaluation/exported_iamc_variables.xlsx",
    output:
        touch(
            RESULTS + "evaluation/HTML/sankey_diagram_EU_2050.html",
        ),


rule validate_pypsa_at:
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
        # fixme: remove unit mark once tests pass
        "pytest -m unit --html {params.rdir}/validity_report.html --result-path={params.rdir}"
