# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

"""
Production implementation of streamlined PyPSA-EUR workflow compose rules.

This file implements the streamlined workflow structure:
base → simplified → clustered → composed → solved

All configuration is now driven by config sections rather than wildcards.
"""


def get_compose_inputs(w):
    """Determine inputs for compose rule based on foresight and horizon."""
    cfg = get_config(w)
    foresight = cfg["foresight"]
    horizon = int(w.horizon)
    sector_enabled = cfg["sector"]["enabled"]
    horizons = cfg["planning_horizons"]
    elec = cfg["electricity"]
    sector = cfg["sector"]
    import_carriers = sector["imports"]["price"] if sector["imports"]["enable"] else {}

    # Electricity-only inputs (always included)
    inputs = dict(
        **input_profile_tech(w),
        **input_class_regions(w),
        **input_conventional(w),
        tech_costs=resources(f"costs_{horizon}_processed.csv"),
        powerplants=(
            resources("powerplants.csv")
            if elec["conventional_carriers"]
            or elec["extendable_carriers"]["Generator"]
            or "hydro" in elec["renewable_carriers"]
            or elec["estimate_battery_capacities"]
            or (
                elec["estimate_renewable_capacities"]["enable"]
                and elec["estimate_renewable_capacities"]["from_powerplantmatching"]
            )
            or (foresight != "overnight" and horizon == horizons[0])
            else []
        ),
        hydro_capacities=ancient("data/hydro_capacities.csv"),
        unit_commitment="data/unit_commitment.csv",
        fuel_price=(
            resources("monthly_fuel_price.csv")
            if cfg["conventional"]["dynamic_fuel_price"]
            else []
        ),
        co2_price=resources("co2_price.csv"),
        eurostat=(
            resources("eurostat_energy_balances.csv")
            if cfg["co2_budget"]["relative"]
            else []
        ),
        co2=(
            rules.retrieve_ghg_emissions.output["csv"]
            if cfg["co2_budget"]["relative"]
            else []
        ),
        load=resources("electricity_demand.nc"),
        snapshot_weightings=resources("snapshot_weightings.csv"),
        network=(
            resources("networks/clustered.nc")
            if not cfg["sector"]["district_heating"]["subnodes"]["enable"]
            else resources("networks/clustered-extended.nc")
        ),
        solar_rooftop_potentials=(
            resources("solar_rooftop_potentials.csv")
            if "solar" in cfg["electricity"]["renewable_carriers"]
            else []
        ),
    )

    # Sector-specific inputs (only when sector coupling is enabled)
    if sector_enabled:
        sector_inputs = dict(
            **input_heat_source_power(w),
            clustered_gas_network=(
                rules.cluster_gas_network.output.clustered_gas_network
                if sector["gas_network"] or sector["H2_retrofit"]
                else []
            ),
            gas_input_nodes_simplified=(
                rules.build_gas_input_locations.output.gas_input_nodes_simplified
                if sector["gas_network"] or {"gas", "H2"} & set(import_carriers)
                else []
            ),
            pop_weighted_energy_totals=resources("pop_weighted_energy_totals.csv"),
            pop_weighted_heat_totals=(
                resources("pop_weighted_heat_totals.csv") if sector["heating"] else []
            ),
            shipping_demand=resources("shipping_demand.csv"),
            transport_demand=resources("transport_demand.csv"),
            transport_data=resources("transport_data.csv"),
            avail_profile=resources("avail_profile.csv"),
            dsm_profile=resources("dsm_profile.csv"),
            heat_dsm_profile=resources("residential_heat_dsm_profile.csv"),
            biomass_potentials=resources("biomass_potentials_{horizon}.csv"),
            h2_cavern=(
                resources("salt_cavern_potentials.csv")
                if sector["hydrogen_underground_storage"]
                else []
            ),
            clustered_pop_layout=resources("pop_layout.csv"),
            industrial_demand=resources("industrial_energy_demand_{horizon}.csv"),
            hourly_heat_demand_total=resources("hourly_heat_demand_total.nc"),
            industrial_production=resources("industrial_production_{horizon}.csv"),
            district_heat_share=resources("district_heat_share_{horizon}-modified.csv"),
            heating_efficiencies=resources("heating_efficiencies.csv"),
            existing_heating_distribution=(
                resources("existing_heating_distribution_{horizon}.csv")
                if not cfg["sector"]["district_heating"]["subnodes"]["enable"]
                else resources(
                    f"existing_heating_distribution_extended_{horizons[0]}.csv"
                )
            ),
            german_chps=resources("german_chp.csv"),
            temp_soil_total=resources("temp_soil_total.nc"),
            temp_air_total=resources("temp_air_total.nc"),
            cop_profiles=resources("cop_profiles_{horizon}.nc"),
            direct_heat_source_utilisation_profiles=resources(
                "direct_heat_source_utilisation_profiles_{horizon}.nc"
            ),
            retro_cost=(
                resources("retro_cost.csv")
                if cfg["sector"]["retrofitting"]["retro_endogen"]
                else []
            ),
            floor_area=(
                resources("floor_area.csv")
                if cfg["sector"]["retrofitting"]["retro_endogen"]
                else []
            ),
            biomass_transport_costs=(
                resources("biomass_transport_costs.csv")
                if cfg["sector"]["biomass_transport"]
                or cfg["sector"]["biomass_spatial"]
                else []
            ),
            sequestration_potential=(
                resources("co2_sequestration_potential.csv")
                if cfg["sector"]["regional_co2_sequestration_potential"]["enable"]
                else []
            ),
            ptes_e_max_pu_profiles=(
                resources("ptes_e_max_pu_profiles_{horizon}.nc")
                if cfg["sector"]["district_heating"]["ptes"]["dynamic_capacity"]
                else []
            ),
            ptes_direct_utilisation_profiles=(
                resources("ptes_direct_utilisation_profiles_{horizon}.nc")
                if cfg["sector"]["district_heating"]["ptes"]["supplemental_heating"][
                    "enable"
                ]
                else []
            ),
            solar_thermal_total=(
                resources("solar_thermal_total.nc")
                if cfg["sector"]["solar_thermal"]
                else []
            ),
            egs_potentials=(
                resources("egs_potentials.csv")
                if cfg["sector"]["enhanced_geothermal"]["enable"]
                else []
            ),
            egs_overlap=(
                resources("egs_overlap.csv")
                if cfg["sector"]["enhanced_geothermal"]["enable"]
                else []
            ),
            egs_capacity_factors=(
                resources("egs_capacity_factors.csv")
                if cfg["sector"]["enhanced_geothermal"]["enable"]
                else []
            ),
            ates_potentials=(
                resources("ates_potentials_{horizon}.csv")
                if cfg["sector"]["district_heating"]["ates"]["enable"]
                else []
            ),
        )
        inputs.update(sector_inputs)
        # pypsa-de specific inputs
        uba_industry_enabled = horizon in cfg["pypsa-de"]["uba_for_industry"]["enable"]
        inputs.update(
            modified_mobility_data=(
                resources(f"modified_mobility_data_{horizon}.csv")
                if sector["transport"]
                else []
            ),
            industrial_demand_2025=(
                resources("industrial_energy_demand_2025.csv")
                if sector["industry"]
                else []
            ),
            industrial_production_per_country_tomorrow=(
                resources(
                    f"industrial_production_per_country_tomorrow_{horizon}-modified.csv"
                )
                if sector["industry"] and uba_industry_enabled
                else []
            ),
            industry_sector_ratios=(
                resources(f"industry_sector_ratios_{horizon}.csv")
                if sector["industry"] and uba_industry_enabled
                else []
            ),
            new_industrial_energy_demand=(
                "data/pypsa-de/UBA_Projektionsbericht2025_Abbildung31_MWMS.csv"
                if sector["industry"] and uba_industry_enabled
                else []
            ),
            onshore_regions=resources("onshore_regions.geojson"),
            regions_offshore=resources("offshore_regions.geojson"),
            offshore_connection_points="data/pypsa-de/offshore_connection_points.csv",
            wkn=(
                rules.cluster_wasserstoff_kernnetz.output.clustered_h2_network
                if config_provider("wasserstoff_kernnetz", "enable")(w)
                else []
            ),
        )

    # Add brownfield inputs for non-first horizons
    if foresight == "overnight" and len(horizons) > 1:
        raise ValueError(
            "Overnight optimization can only be run for a single planning horizon."
        )

    if horizon != horizons[0]:
        # Not first horizon - need previous network
        prev_horizon = horizons[horizons.index(horizon) - 1]

        if foresight == "myopic":
            # Myopic uses solved network from previous horizon
            inputs["network_previous"] = RESULTS + f"networks/solved_{prev_horizon}.nc"
        elif foresight == "perfect":
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = resources(
                f"networks/composed_{prev_horizon}.nc"
            )
        else:
            raise ValueError(f"Invalid foresight type: {foresight}")

    # imported by compose_network.py; listed so code changes trigger reruns
    inputs["code_dependencies"] = [
        "scripts/add_electricity.py",
        "scripts/add_existing_baseyear.py",
        "scripts/add_brownfield.py",
        "scripts/prepare_network.py",
        "scripts/prepare_perfect_foresight.py",
        "scripts/prepare_sector_network.py",
        "scripts/pypsa-de/modify_prenetwork.py",
        "scripts/_helpers.py",
    ]

    return inputs


# Main composition rule - combines all network building steps
rule compose_network:
    input:
        unpack(get_compose_inputs),
    output:
        resources("networks/composed_{horizon}.nc"),
    log:
        logs("compose_network_{horizon}.log"),
    benchmark:
        benchmarks("compose_network_{horizon}")
    threads: 1
    resources:
        mem_mb=10000,
    params:
        foresight=config_provider("foresight"),
        electricity=config_provider("electricity"),
        sector=config_provider("sector"),
        clustering=config_provider("clustering"),
        clustering_temporal=config_provider("clustering", "temporal"),
        existing_capacities=config_provider("existing_capacities"),
        pypsa_eur=config_provider("pypsa_eur"),
        renewable=config_provider("renewable"),
        conventional=config_provider("conventional"),
        costs=config_provider("costs"),
        emission_prices=config_provider("costs", "emission_prices"),
        load=config_provider("load"),
        lines=config_provider("lines"),
        links=config_provider("links"),
        transmission_losses=config_provider("solving", "options", "transmission_losses"),
        industry=config_provider("industry"),
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        add_district_heating_subnodes=config_provider(
            "sector", "district_heating", "subnodes", "enable"
        ),
        horizons=config_provider("planning_horizons"),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        conventional_carriers=config_provider("electricity", "conventional_carriers"),
        fuel_carriers=config_provider("existing_capacities", "conventional_carriers"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        h2_retrofit=config_provider("sector", "H2_retrofit"),
        h2_retrofit_capacity_per_ch4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        capacity_threshold=config_provider("existing_capacities", "threshold_capacity"),
        temperature_limited_stores=config_provider(
            "sector", "district_heating", "temperature_limited_stores"
        ),
        tes=config_provider("sector", "tes"),
        dynamic_ptes_capacity=config_provider(
            "sector", "district_heating", "ptes", "dynamic_capacity"
        ),
        direct_utilisation_heat_sources=config_provider(
            "sector", "district_heating", "direct_utilisation_heat_sources"
        ),
        co2_budget=config_provider("co2_budget"),
        adjustments=config_provider("adjustments"),
        # pypsa-de specific
        planning_horizons=config_provider("planning_horizons"),
        efuel_export_ban=config_provider("solving", "constraints", "efuel_export_ban"),
        enable_kernnetz=config_provider("wasserstoff_kernnetz", "enable"),
        pypsa_de_enabled=config_provider("pypsa-de", "enable"),
        technology_occurrence=config_provider("first_technology_occurrence"),
        fossil_boiler_ban=config_provider("new_decentral_fossil_boiler_ban"),
        coal_ban=config_provider("coal_generation_ban"),
        nuclear_ban=config_provider("nuclear_generation_ban"),
        H2_transmission_efficiency=config_provider(
            "sector", "transmission_efficiency", "H2 pipeline"
        ),
        H2_retrofit=config_provider("sector", "H2_retrofit"),
        transmission_costs=config_provider("costs", "transmission"),
        must_run=config_provider("must_run"),
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
            "pypsa-de", "uba_for_industry", "scale_non_energy"
        ),
        limit_cross_border_flows_ac=config_provider(
            "pypsa-de", "limit_cross_border_flows_ac"
        ),
        space_heat_DE_factor=config_provider("pypsa-de", "reduce_space_heat_DE_factor"),
        space_heat_EU_factor=config_provider(
            "sector", "reduce_space_heat_exogenously_factor"
        ),
        deactivate_early_transmission_expansion=config_provider(
            "pypsa-de", "deactivate_early_transmission_expansion"
        ),
    message:
        "Composing network for horizon {wildcards.horizon}"
    script:
        scripts("compose_network.py")
