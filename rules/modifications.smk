rule custom_administrative_clustering:
    params:
        clustering=config_provider("clustering", "mode"),
        admin_levels=config_provider("clustering", "administrative"),
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
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


if (
    config_provider("custom_administrative_clustering", "enable")
    # and config_provider("clustering", "mode") == "administrative"
):

    ruleorder: custom_administrative_clustering > base_network
