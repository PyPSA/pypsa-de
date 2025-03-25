rule modify_nuts3_shapes:
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
    config["custom_administrative_clustering"].get("enable", False)
    and config_provider("clustering", "mode") == "administrative"
):

    ruleorder: modify_nuts3_shapes > base_network
