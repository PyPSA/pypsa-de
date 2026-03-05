# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC BY 4.0


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
        scripts("pypsa-de/build_scenarios.py")


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
        scripts("pypsa-de/prepare_district_heating_subnodes.py")


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
        scripts("pypsa-de/add_district_heating_subnodes.py")


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
        scripts("pypsa-de/modify_district_heat_share.py")


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
        scripts("pypsa-de/modify_prenetwork.py")


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
        scripts("pypsa-de/modify_existing_heating.py")


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
        scripts("pypsa-de/build_existing_chp_de.py")


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
        scripts("pypsa-de/modify_industry_production.py")


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
        scripts("pypsa-de/build_wasserstoff_kernnetz.py")


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
        scripts("pypsa-de/cluster_wasserstoff_kernnetz.py")
