# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC BY 4.0


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
        scripts("pypsa-de/export_ariadne_variables.py")


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
        scripts("pypsa-de/plot_ariadne_variables.py")


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
        scripts("pypsa-de/plot_hydrogen_network_incl_kernnetz.py")


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
        scripts("pypsa-de/plot_ariadne_report.py")


rule ariadne_report_only:
    input:
        expand(
            RESULTS + "ariadne/report/elec_price_duration_curve.pdf",
            run=config_provider("run", "name"),
        ),


rule plot_scenario_comparison:
    input:
        exported_variables=expand(
            RESULTS + "ariadne/exported_variables_full.xlsx",
            run=config_provider("run", "name"),
        ),
    output:
        price_carbon="results/"
        + config["run"]["prefix"]
        + "/scenario_comparison/Price-Carbon.png",
    script:
        scripts("pypsa-de/plot_scenario_comparison.py")


rule compare_scenarios:
    input:
        price_carbon="results/"
        + config["run"]["prefix"]
        + "/scenario_comparison/Price-Carbon.png",
        # expand(
        #     RESULTS + "ariadne/capacity_detailed.png",
        #     run=config_provider("run", "name"),
        # ),


rule ariadne_all:
    input:
        expand(
            RESULTS + "ariadne/report/elec_price_duration_curve.pdf",
            run=config_provider("run", "name"),
        ),
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
        price_carbon="results/"
        + config["run"]["prefix"]
        + "/scenario_comparison/Price-Carbon.png",
