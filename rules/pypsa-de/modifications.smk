# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC BY 4.0


rule build_scenarios:
    input:
        ariadne_database="data/ariadne_database.csv",
        scenario_yaml=config["run"]["scenarios"]["manual_file"],
    output:
        scenario_yaml=config["run"]["scenarios"]["file"],
    log:
        "logs/build_scenarios.log",
    params:
        scenarios=config["run"]["name"],
        leitmodelle=config["pypsa-de"]["leitmodelle"],
    script:
        scripts("pypsa-de/build_scenarios.py")


rule build_exogenous_mobility_data:
    input:
        ariadne="data/ariadne_database.csv",
        energy_totals=resources("energy_totals.csv"),
    output:
        mobility_data=resources("modified_mobility_data_{horizon}.csv"),
    log:
        logs("build_exogenous_mobility_data_{horizon}.log"),
    resources:
        mem_mb=1000,
    params:
        reference_scenario=config_provider("pypsa-de", "reference_scenario"),
        leitmodelle=config_provider("pypsa-de", "leitmodelle"),
        uba_for_mobility=config_provider("pypsa-de", "uba_for_mobility"),
        shipping_oil_share=config_provider("sector", "shipping_oil_share"),
        aviation_demand_factor=config_provider("sector", "aviation_demand_factor"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    script:
        scripts("pypsa-de/build_exogenous_mobility_data.py")


rule build_egon_data:
    input:
        demandregio_spatial=f"{EGON['folder']}/demandregio_spatial_2018.json",
        mapping_38_to_4=storage(
            "https://ffeopendatastorage.blob.core.windows.net/opendata/mapping_from_4_to_38.json",
            keep_local=True,
        ),
        mapping_technologies=f"{EGON['folder']}/mapping_technologies.json",
        nuts3=resources("nuts3_shapes.geojson"),
    output:
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
    log:
        logs("build_egon_data.log"),
    script:
        scripts("pypsa-de/build_egon_data.py")


rule prepare_district_heating_subnodes:
    input:
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
        onshore_regions=resources("onshore_regions.geojson"),
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
        district_heating_subnodes=resources("district_heating_subnodes.geojson"),
        onshore_regions_extended=resources("onshore_regions_extended.geojson"),
        onshore_regions_restricted=resources("onshore_regions_restricted.geojson"),
    resources:
        mem_mb=20000,
    params:
        district_heating=config_provider("sector", "district_heating"),
        baseyear=config_provider("planning_horizons", 0),
    script:
        scripts("pypsa-de/prepare_district_heating_subnodes.py")


rule extend_existing_heating_distribution:
    input:
        existing_heating_distribution=resources(
            f"existing_heating_distribution_{config['planning_horizons'][0]}.csv"
        ),
        subnodes=resources("district_heating_subnodes.geojson"),
    output:
        existing_heating_distribution_extended=resources(
            f"existing_heating_distribution_extended_{config['planning_horizons'][0]}.csv"
        ),
        district_heating_subnodes_selected=resources(
            "district_heating_subnodes_selected.geojson"
        ),
    params:
        nlargest=config_provider("sector", "district_heating", "subnodes", "nlargest"),
    script:
        scripts("pypsa-de/extend_existing_heating_distribution.py")


rule add_district_heating_subnodes:
    input:
        unpack(input_heat_source_power),
        network=resources("networks/composed_{horizon}.nc"),
        subnodes=resources("district_heating_subnodes.geojson"),
        nuts3=resources("nuts3_shapes.geojson"),
        onshore_regions=resources("onshore_regions.geojson"),
        fernwaermeatlas="data/fernwaermeatlas/fernwaermeatlas.xlsx",
        cities="data/fernwaermeatlas/cities_geolocations.geojson",
        cop_profiles=resources("cop_profiles_{horizon}.nc"),
        direct_heat_source_utilisation_profiles=resources(
            "direct_heat_source_utilisation_profiles_{horizon}.nc"
        ),
        lau_regions=rules.retrieve_lau_regions.output["zip"],
    output:
        network=resources("networks/composed_with_subnodes_{horizon}.nc"),
    resources:
        mem_mb=10000,
    params:
        district_heating=config_provider("sector", "district_heating"),
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
    script:
        scripts("pypsa-de/add_district_heating_subnodes.py")


ruleorder: modify_district_heat_share > build_district_heat_share


rule modify_district_heat_share:
    input:
        heating_technologies_nuts3=resources("heating_technologies_nuts3.geojson"),
        onshore_regions=resources("onshore_regions.geojson"),
        district_heat_share=resources("district_heat_share_{horizon}.csv"),
    output:
        district_heat_share=resources("district_heat_share_{horizon}-modified.csv"),
    log:
        logs("modify_district_heat_share_{horizon}.log"),
    resources:
        mem_mb=1000,
    params:
        district_heating=config_provider("sector", "district_heating"),
    script:
        scripts("pypsa-de/modify_district_heat_share.py")


ruleorder: modify_industry_production > build_industrial_production_per_country_tomorrow


rule modify_existing_heating:
    input:
        ariadne="data/ariadne_database.csv",
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
    output:
        existing_heating=resources("existing_heating.csv"),
    log:
        logs("modify_existing_heating.log"),
    resources:
        mem_mb=1000,
    script:
        scripts("pypsa-de/modify_existing_heating.py")


rule build_existing_chp_de:
    input:
        mastr_biomass="data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
        mastr_combustion="data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
        plz_mapping=storage(
            "https://raw.githubusercontent.com/WZBSocialScienceCenter/plz_geocoord/master/plz_geocoord.csv",
            keep_local=True,
        ),
        regions=resources("onshore_regions.geojson"),
        district_heating_subnodes=lambda w: (
            resources("district_heating_subnodes.geojson")
            if config_provider("sector", "district_heating", "subnodes", "enable")(w)
            else []
        ),
    output:
        german_chp=resources("german_chp.csv"),
    log:
        logs("build_existing_chp_de.log"),
    resources:
        mem_mb=4000,
    params:
        district_heating_subnodes=config_provider(
            "sector", "district_heating", "subnodes"
        ),
    script:
        scripts("pypsa-de/build_existing_chp_de.py")


rule modify_industry_production:
    input:
        ariadne="data/ariadne_database.csv",
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{horizon}.csv"
        ),
    output:
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{horizon}-modified.csv"
        ),
    log:
        logs("modify_industry_production_{horizon}.log"),
    resources:
        mem_mb=1000,
    params:
        reference_scenario=config_provider("pypsa-de", "reference_scenario"),
    script:
        scripts("pypsa-de/modify_industry_production.py")


rule build_wasserstoff_kernnetz:
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
        onshore_regions=resources("onshore_regions_base.geojson"),
        regions_offshore=resources("offshore_regions_base.geojson"),
    output:
        cleaned_wasserstoff_kernnetz=resources("wasserstoff_kernnetz.csv"),
    log:
        logs("build_wasserstoff_kernnetz.log"),
    params:
        kernnetz=config_provider("wasserstoff_kernnetz"),
    script:
        scripts("pypsa-de/build_wasserstoff_kernnetz.py")


rule cluster_wasserstoff_kernnetz:
    input:
        cleaned_h2_network=resources("wasserstoff_kernnetz.csv"),
        onshore_regions=resources("onshore_regions.geojson"),
        regions_offshore=resources("offshore_regions.geojson"),
    output:
        clustered_h2_network=resources("wasserstoff_kernnetz_clustered.csv"),
    log:
        logs("cluster_wasserstoff_kernnetz.log"),
    params:
        kernnetz=config_provider("wasserstoff_kernnetz"),
    script:
        scripts("pypsa-de/cluster_wasserstoff_kernnetz.py")
