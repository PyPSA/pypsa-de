import logging
import sys

import pandas as pd
from xarray import DataArray

from scripts.prepare_sector_network import determine_emission_sectors

logger = logging.getLogger(__name__)


def add_capacity_limits(
    n, investment_year, limits_capacity, snakemake, sense="maximum"
):
    for c in n.iterate_components(limits_capacity):
        logger.info(f"Adding {sense} constraints for {c.list_name}")

        attr = "e" if c.name == "Store" else "p"
        units = "MWh or tCO2" if c.name == "Store" else "MW"

        for carrier in limits_capacity[c.name]:
            for ct in limits_capacity[c.name][carrier]:
                if investment_year not in limits_capacity[c.name][carrier][ct].keys():
                    continue

                limit = 1e3 * limits_capacity[c.name][carrier][ct][investment_year]

                logger.info(
                    f"Adding constraint on {c.name} {carrier} capacity in {ct} to be {sense} {limit} {units}"
                )

                valid_components = (
                    (c.df.index.str[:2] == ct)
                    & (c.df.carrier.str[: len(carrier)] == carrier)
                    & ~c.df.carrier.str.contains("thermal")
                )  # exclude solar thermal

                existing_index = c.df.index[
                    valid_components & ~c.df[attr + "_nom_extendable"]
                ]
                extendable_index = c.df.index[
                    valid_components & c.df[attr + "_nom_extendable"]
                ]

                existing_capacity = c.df.loc[existing_index, attr + "_nom"].sum()

                logger.info(
                    f"Existing {c.name} {carrier} capacity in {ct}: {existing_capacity} {units}"
                )
                if extendable_index.empty:
                    logger.warning(
                        f"No extendable {c.name} {carrier} capacities found in {ct}. Skipping."
                    )
                    continue
                nom = n.model[c.name + "-" + attr + "_nom"].loc[extendable_index]

                lhs = nom.sum()

                cname = f"capacity_{sense}-{ct}-{c.name}-{carrier.replace(' ', '-')}"

                if cname in n.global_constraints.index:
                    logger.warning(
                        f"Global constraint {cname} already exists. Dropping and adding it again."
                    )
                    n.global_constraints.drop(cname, inplace=True)

                rhs = limit - existing_capacity

                if sense == "maximum":
                    if rhs <= 0:
                        logger.warning(
                            f"Existing capacity in {ct} for carrier {carrier} already exceeds the limit of {limit} MW. Limiting capacity expansion for this investment period to 0."
                        )
                        rhs = 0

                    n.model.add_constraints(
                        lhs <= rhs,
                        name=f"GlobalConstraint-{cname}",
                    )
                    n.add(
                        "GlobalConstraint",
                        cname,
                        constant=rhs,
                        sense="<=",
                        type="",
                        carrier_attribute="",
                    )

                elif sense == "minimum":
                    n.model.add_constraints(
                        lhs >= rhs,
                        name=f"GlobalConstraint-{cname}",
                    )
                    n.add(
                        "GlobalConstraint",
                        cname,
                        constant=rhs,
                        sense=">=",
                        type="",
                        carrier_attribute="",
                    )
                else:
                    logger.error("sense {sense} not recognised")
                    sys.exit()


def add_power_limits(n, investment_year, limits_power_max):
    """
    " Restricts the maximum inflow/outflow of electricity from/to a country.
    """
    for ct in limits_power_max:
        if investment_year not in limits_power_max[ct].keys():
            continue

        lim = 1e3 * limits_power_max[ct][investment_year]  # in MW

        logger.info(
            f"Adding constraint on electricity import/export from/to {ct} to be < {lim} MW"
        )
        # identify interconnectors

        incoming_lines = n.lines[
            (n.lines.carrier == "AC")
            & (n.lines.bus0.str[:2] != ct)
            & (n.lines.bus1.str[:2] == ct)
            & n.lines.active
        ]
        outgoing_lines = n.lines[
            (n.lines.carrier == "AC")
            & (n.lines.bus0.str[:2] == ct)
            & (n.lines.bus1.str[:2] != ct)
            & n.lines.active
        ]
        incoming_links = n.links[
            (n.links.carrier == "DC")
            & (n.links.bus0.str[:2] != ct)
            & (n.links.bus1.str[:2] == ct)
            & n.links.active
        ]
        outgoing_links = n.links[
            (n.links.carrier == "DC")
            & (n.links.bus0.str[:2] == ct)
            & (n.links.bus1.str[:2] != ct)
            & n.links.active
        ]

        for t in n.snapshots:
            # For incoming flows s > 0 means imports, s < 0 exports
            # For outgoing flows s > 0 means exports, s < 0 imports
            # to get the positive and negative parts separately, we use auxiliary variables
            incoming_lines_var = n.model["Line-s"].loc[t, incoming_lines.index]
            n.model.add_variables(
                coords=[incoming_lines.index],
                name=f"Line-s-incoming-{ct}-aux-pos-{t}",
                lower=0,
                upper=incoming_lines.s_nom_max,
            )
            n.model.add_variables(
                coords=[incoming_lines.index],
                name=f"Line-s-incoming-{ct}-aux-neg-{t}",
                lower=-incoming_lines.s_nom_max,
                upper=0,
            )
            n.model.add_constraints(
                n.model[f"Line-s-incoming-{ct}-aux-pos-{t}"] >= incoming_lines_var,
                name=f"Line-s-incoming-{ct}-aux-pos-constr-{t}",
            )
            n.model.add_constraints(
                n.model[f"Line-s-incoming-{ct}-aux-neg-{t}"] <= incoming_lines_var,
                name=f"Line-s-incoming-{ct}-aux-neg-constr-{t}",
            )

            outgoing_lines_var = n.model["Line-s"].loc[t, outgoing_lines.index]
            n.model.add_variables(
                coords=[outgoing_lines.index],
                name=f"Line-s-outgoing-{ct}-aux-pos-{t}",
                lower=0,
                upper=outgoing_lines.s_nom_max,
            )
            n.model.add_variables(
                coords=[outgoing_lines.index],
                name=f"Line-s-outgoing-{ct}-aux-neg-{t}",
                lower=-outgoing_lines.s_nom_max,
                upper=0,
            )
            n.model.add_constraints(
                n.model[f"Line-s-outgoing-{ct}-aux-pos-{t}"] >= outgoing_lines_var,
                name=f"Line-s-outgoing-{ct}-aux-pos-constr-{t}",
            )
            n.model.add_constraints(
                n.model[f"Line-s-outgoing-{ct}-aux-neg-{t}"] <= outgoing_lines_var,
                name=f"Line-s-outgoing-{ct}-aux-neg-constr-{t}",
            )

            incoming_links_var = n.model["Link-p"].loc[t, incoming_links.index]
            n.model.add_variables(
                coords=[incoming_links.index],
                name=f"Link-p-incoming-{ct}-aux-pos-{t}",
                lower=0,
                upper=incoming_links.p_nom_max,
            )
            n.model.add_variables(
                coords=[incoming_links.index],
                name=f"Link-p-incoming-{ct}-aux-neg-{t}",
                lower=-incoming_links.p_nom_max,
                upper=0,
            )
            n.model.add_constraints(
                n.model[f"Link-p-incoming-{ct}-aux-pos-{t}"] >= incoming_links_var,
                name=f"Link-p-incoming-{ct}-aux-pos-constr-{t}",
            )
            n.model.add_constraints(
                n.model[f"Link-p-incoming-{ct}-aux-neg-{t}"] <= incoming_links_var,
                name=f"Link-p-incoming-{ct}-aux-neg-constr-{t}",
            )

            outgoing_links_var = n.model["Link-p"].loc[t, outgoing_links.index]
            n.model.add_variables(
                coords=[outgoing_links.index],
                name=f"Link-p-outgoing-{ct}-aux-pos-{t}",
                lower=0,
                upper=outgoing_links.p_nom_max,
            )
            n.model.add_variables(
                coords=[outgoing_links.index],
                name=f"Link-p-outgoing-{ct}-aux-neg-{t}",
                lower=-outgoing_links.p_nom_max,
                upper=0,
            )
            n.model.add_constraints(
                n.model[f"Link-p-outgoing-{ct}-aux-pos-{t}"] >= outgoing_links_var,
                name=f"Link-p-outgoing-{ct}-aux-pos-constr-{t}",
            )
            n.model.add_constraints(
                n.model[f"Link-p-outgoing-{ct}-aux-neg-{t}"] <= outgoing_links_var,
                name=f"Link-p-outgoing-{ct}-aux-neg-constr-{t}",
            )
            # To constraint the absolute values of imports and exports, we have to sum the
            # corresponding positive and negative flows separately, using auxiliary variables
            import_lhs = (
                n.model[f"Link-p-incoming-{ct}-aux-pos-{t}"].sum()
                + n.model[f"Line-s-incoming-{ct}-aux-pos-{t}"].sum()
                - n.model[f"Link-p-outgoing-{ct}-aux-neg-{t}"].sum()
                - n.model[f"Line-s-outgoing-{ct}-aux-neg-{t}"].sum()
            ) / 10  # divide by 10 to improve numerical stability
            export_lhs = (
                n.model[f"Link-p-outgoing-{ct}-aux-pos-{t}"].sum()
                + n.model[f"Line-s-outgoing-{ct}-aux-pos-{t}"].sum()
                - n.model[f"Link-p-incoming-{ct}-aux-neg-{t}"].sum()
                - n.model[f"Line-s-incoming-{ct}-aux-neg-{t}"].sum()
            ) / 10

            n.model.add_constraints(
                import_lhs <= lim / 10, name=f"Power-import-limit-{ct}-{t}"
            )
            n.model.add_constraints(
                export_lhs <= lim / 10, name=f"Power-export-limit-{ct}-{t}"
            )


def h2_import_limits(n, investment_year, limits_volume_max):
    for ct in limits_volume_max["h2_import"]:
        limit = limits_volume_max["h2_import"][ct][investment_year] * 1e6

        logger.info(f"limiting H2 imports in {ct} to {limit / 1e6} TWh/a")
        pipeline_carrier = [
            "H2 pipeline",
            "H2 pipeline (Kernnetz)",
            "H2 pipeline retrofitted",
        ]
        incoming = n.links.index[
            (n.links.carrier.isin(pipeline_carrier))
            & (n.links.bus0.str[:2] != ct)
            & (n.links.bus1.str[:2] == ct)
        ]
        outgoing = n.links.index[
            (n.links.carrier.isin(pipeline_carrier))
            & (n.links.bus0.str[:2] == ct)
            & (n.links.bus1.str[:2] != ct)
        ]

        incoming_p = (
            n.model["Link-p"].loc[:, incoming] * n.snapshot_weightings.generators
        ).sum()
        outgoing_p = (
            n.model["Link-p"].loc[:, outgoing] * n.snapshot_weightings.generators
        ).sum()

        lhs = incoming_p - outgoing_p

        cname = f"H2_import_limit-{ct}"

        n.model.add_constraints(lhs <= limit, name=f"GlobalConstraint-{cname}")

        if cname in n.global_constraints.index:
            logger.warning(
                f"Global constraint {cname} already exists. Dropping and adding it again."
            )
            n.global_constraints.drop(cname, inplace=True)

        n.add(
            "GlobalConstraint",
            cname,
            constant=limit,
            sense="<=",
            type="",
            carrier_attribute="",
        )

        logger.info("Adding H2 export ban")

        cname = f"H2_export_ban-{ct}"

        n.model.add_constraints(lhs >= 0, name=f"GlobalConstraint-{cname}")

        if cname in n.global_constraints.index:
            logger.warning(
                f"Global constraint {cname} already exists. Dropping and adding it again."
            )
            n.global_constraints.drop(cname, inplace=True)

        n.add(
            "GlobalConstraint",
            cname,
            constant=0,
            sense=">=",
            type="",
            carrier_attribute="",
        )


def h2_production_limits(n, investment_year, limits_volume_min, limits_volume_max):
    for ct in limits_volume_max["electrolysis"]:
        if ct not in limits_volume_min["electrolysis"]:
            logger.warning(
                f"no lower limit for H2 electrolysis in {ct} assuming 0 TWh/a"
            )
            limit_lower = 0
        else:
            limit_lower = limits_volume_min["electrolysis"][ct][investment_year] * 1e6

        limit_upper = limits_volume_max["electrolysis"][ct][investment_year] * 1e6

        logger.info(
            f"limiting H2 electrolysis in DE between {limit_lower / 1e6} and {limit_upper / 1e6} TWh/a"
        )

        production = n.links[
            (n.links.carrier == "H2 Electrolysis") & (n.links.bus0.str.contains(ct))
        ].index
        efficiency = n.links.loc[production, "efficiency"]

        lhs = (
            n.model["Link-p"].loc[:, production]
            * n.snapshot_weightings.generators
            * efficiency
        ).sum()

        cname_upper = f"H2_production_limit_upper-{ct}"
        cname_lower = f"H2_production_limit_lower-{ct}"

        n.model.add_constraints(
            lhs <= limit_upper, name=f"GlobalConstraint-{cname_upper}"
        )

        n.model.add_constraints(
            lhs >= limit_lower, name=f"GlobalConstraint-{cname_lower}"
        )

        if cname_upper not in n.global_constraints.index:
            n.add(
                "GlobalConstraint",
                cname_upper,
                constant=limit_upper,
                sense="<=",
                type="",
                carrier_attribute="",
            )
        if cname_lower not in n.global_constraints.index:
            n.add(
                "GlobalConstraint",
                cname_lower,
                constant=limit_lower,
                sense=">=",
                type="",
                carrier_attribute="",
            )


def electricity_import_limits(n, investment_year, limits_volume_max):
    for ct in limits_volume_max["electricity_import"]:
        limit = limits_volume_max["electricity_import"][ct][investment_year] * 1e6

        if limit < 0:
            limit *= n.snapshot_weightings.generators.sum() / 8760

        logger.info(f"limiting electricity imports in {ct} to {limit / 1e6} TWh/a")

        incoming_line = n.lines.index[
            (n.lines.carrier == "AC")
            & (n.lines.bus0.str[:2] != ct)
            & (n.lines.bus1.str[:2] == ct)
        ]
        outgoing_line = n.lines.index[
            (n.lines.carrier == "AC")
            & (n.lines.bus0.str[:2] == ct)
            & (n.lines.bus1.str[:2] != ct)
        ]

        incoming_link = n.links.index[
            (n.links.carrier == "DC")
            & (n.links.bus0.str[:2] != ct)
            & (n.links.bus1.str[:2] == ct)
        ]
        outgoing_link = n.links.index[
            (n.links.carrier == "DC")
            & (n.links.bus0.str[:2] == ct)
            & (n.links.bus1.str[:2] != ct)
        ]

        incoming_line_p = (
            n.model["Line-s"].loc[:, incoming_line] * n.snapshot_weightings.generators
        ).sum()
        outgoing_line_p = (
            n.model["Line-s"].loc[:, outgoing_line] * n.snapshot_weightings.generators
        ).sum()

        incoming_link_p = (
            n.model["Link-p"].loc[:, incoming_link] * n.snapshot_weightings.generators
        ).sum()
        outgoing_link_p = (
            n.model["Link-p"].loc[:, outgoing_link] * n.snapshot_weightings.generators
        ).sum()

        lhs = (incoming_link_p - outgoing_link_p) + (incoming_line_p - outgoing_line_p)

        cname = f"Electricity_import_limit-{ct}"

        n.model.add_constraints(lhs <= limit, name=f"GlobalConstraint-{cname}")

        if cname in n.global_constraints.index:
            logger.warning(
                f"Global constraint {cname} already exists. Dropping and adding it again."
            )
            n.global_constraints.drop(cname, inplace=True)

        n.add(
            "GlobalConstraint",
            cname,
            constant=limit,
            sense="<=",
            type="",
            carrier_attribute="",
        )


def add_national_co2_budgets(n, snakemake, national_co2_budgets, investment_year):
    """
    Add a set of emissions limit constraints for specified countries.

    The countries and emissions limits are specified in the config file entry 'co2_budget_national'.

    Parameters
    ----------
    n : pypsa.Network
    snakemake : snakemake.io.Snakemake
    national_co2_budgets : dict
    investment_year : int

    """
    logger.info("Adding national CO2 budgets")
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

    sectors = determine_emission_sectors(n.config["sector"])

    # convert MtCO2 to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_total_totals = co2_totals[sectors].sum(axis=1) * nyears

    for ct in national_co2_budgets:
        if ct != "DE":
            logger.error(
                f"CO2 budget for countries other than `DE` is not yet supported. Found country {ct}. Please check the config file."
            )

        limit = co2_total_totals[ct] * national_co2_budgets[ct][investment_year]
        logger.info(
            f"Limiting emissions in country {ct} to {national_co2_budgets[ct][investment_year]:.1%} of "
            f"1990 levels, i.e. {limit:,.2f} tCO2/a",
        )

        lhs = []

        for port in [col[3:] for col in n.links if col.startswith("bus")]:
            links = n.links.index[
                (n.links.index.str[:2] == ct)
                & (n.links[f"bus{port}"] == "co2 atmosphere")
                & (
                    n.links.carrier != "kerosene for aviation"
                )  # first exclude aviation to multiply it with a domestic factor later
            ]

            logger.info(
                f"For {ct} adding following link carriers to port {port} CO2 constraint: {n.links.loc[links, 'carrier'].unique()}"
            )

            if port == "0":
                efficiency = -1.0
            elif port == "1":
                efficiency = n.links.loc[links, "efficiency"]
            else:
                efficiency = n.links.loc[links, f"efficiency{port}"]

            lhs.append(
                (
                    n.model["Link-p"].loc[:, links]
                    * efficiency
                    * n.snapshot_weightings.generators
                ).sum()
            )

        # Aviation demand
        energy_totals = pd.read_csv(snakemake.input.energy_totals, index_col=[0, 1])
        domestic_aviation = energy_totals.loc[
            (ct, snakemake.params.energy_year), "total domestic aviation"
        ]
        international_aviation = energy_totals.loc[
            (ct, snakemake.params.energy_year), "total international aviation"
        ]
        domestic_factor = domestic_aviation / (
            domestic_aviation + international_aviation
        )
        aviation_links = n.links[
            (n.links.index.str[:2] == ct) & (n.links.carrier == "kerosene for aviation")
        ]
        lhs.append
        (
            n.model["Link-p"].loc[:, aviation_links.index]
            * aviation_links.efficiency2
            * n.snapshot_weightings.generators
        ).sum() * domestic_factor
        logger.info(
            f"Adding domestic aviation emissions for {ct} with a factor of {domestic_factor}"
        )

        # Adding Efuel imports and exports to constraint
        incoming_oil = n.links.index[n.links.index == f"EU renewable oil -> {ct} oil"]
        outgoing_oil = n.links.index[n.links.index == f"{ct} renewable oil -> EU oil"]

        lhs.append(
            (
                -1
                * n.model["Link-p"].loc[:, incoming_oil]
                * 0.2571
                * n.snapshot_weightings.generators
            ).sum()
        )
        lhs.append(
            (
                n.model["Link-p"].loc[:, outgoing_oil]
                * 0.2571
                * n.snapshot_weightings.generators
            ).sum()
        )

        incoming_methanol = n.links.index[
            n.links.index == f"EU methanol -> {ct} methanol"
        ]
        outgoing_methanol = n.links.index[
            n.links.index == f"{ct} methanol -> EU methanol"
        ]

        lhs.append(
            (
                -1
                * n.model["Link-p"].loc[:, incoming_methanol]
                / snakemake.config["sector"]["MWh_MeOH_per_tCO2"]
                * n.snapshot_weightings.generators
            ).sum()
        )

        lhs.append(
            (
                n.model["Link-p"].loc[:, outgoing_methanol]
                / snakemake.config["sector"]["MWh_MeOH_per_tCO2"]
                * n.snapshot_weightings.generators
            ).sum()
        )

        # Methane
        incoming_CH4 = n.links.index[n.links.index == f"EU renewable gas -> {ct} gas"]
        outgoing_CH4 = n.links.index[n.links.index == f"{ct} renewable gas -> EU gas"]

        lhs.append(
            (
                -1
                * n.model["Link-p"].loc[:, incoming_CH4]
                * 0.198
                * n.snapshot_weightings.generators
            ).sum()
        )

        lhs.append(
            (
                n.model["Link-p"].loc[:, outgoing_CH4]
                * 0.198
                * n.snapshot_weightings.generators
            ).sum()
        )

        lhs = sum(lhs)

        cname = f"co2_limit-{ct}"

        n.model.add_constraints(
            lhs <= limit,
            name=f"GlobalConstraint-{cname}",
        )

        if cname in n.global_constraints.index:
            logger.warning(
                f"Global constraint {cname} already exists. Dropping and adding it again."
            )
            n.global_constraints.drop(cname, inplace=True)

        n.add(
            "GlobalConstraint",
            cname,
            constant=limit,
            sense="<=",
            type="",
            carrier_attribute="",
        )


def add_decentral_heat_pump_budgets(n, decentral_heat_pump_budgets, investment_year):
    carriers = [
        "rural air heat pump",
        "rural ground heat pump",
        "urban decentral air heat pump",
        "rural resistive heater",
        "urban decentral resistive heater",
    ]

    heat_pumps = n.links.index[n.links.carrier.isin(carriers)]

    if heat_pumps.empty:
        logger.warning(
            "No heat pumps found in the network. Skipping decentral heat pump budgets."
        )
        return

    if investment_year not in decentral_heat_pump_budgets["DE"].keys():
        logger.warning(
            f"No decentral heat pump budget for {investment_year} found in the config file. Skipping."
        )
        return

    logger.info("Adding decentral heat pump budgets")

    for ct in decentral_heat_pump_budgets:
        if ct != "DE":
            logger.error(
                f"Heat pump budget for countries other than `DE` is not yet supported. Found country {ct}. Please check the config file."
            )

        limit = decentral_heat_pump_budgets[ct][investment_year] * 1e6

        logger.info(
            f"Limiting decentral heat pump electricity consumption in country {ct} to {decentral_heat_pump_budgets[ct][investment_year]:.1%} MWh.",
        )
        heat_pumps = heat_pumps[heat_pumps.str.startswith(ct)]

        lhs = []

        lhs.append(
            (
                n.model["Link-p"].loc[:, heat_pumps] * n.snapshot_weightings.generators
            ).sum()
        )

        lhs = sum(lhs)

        cname = f"decentral_heat_pump_limit-{ct}"
        if cname in n.global_constraints.index:
            logger.warning(
                f"Global constraint {cname} already exists. Dropping and adding it again."
            )
            n.global_constraints.drop(cname, inplace=True)

        n.model.add_constraints(
            lhs <= limit,
            name=f"GlobalConstraint-{cname}",
        )
        n.add(
            "GlobalConstraint",
            cname,
            constant=limit,
            sense="<=",
            type="",
            carrier_attribute="",
        )


def force_boiler_profiles_existing_per_boiler(n):
    """
    This scales each boiler dispatch to be proportional to the load profile.
    """

    logger.info(
        "Forcing each existing boiler dispatch to be proportional to the load profile"
    )

    decentral_boilers = n.links.index[
        n.links.carrier.str.contains("boiler")
        & ~n.links.carrier.str.contains("urban central")
        & ~n.links.p_nom_extendable
    ]

    if decentral_boilers.empty:
        return

    boiler_loads = n.links.loc[decentral_boilers, "bus1"]
    boiler_loads = boiler_loads[boiler_loads.isin(n.loads_t.p_set.columns)]
    decentral_boilers = boiler_loads.index
    boiler_profiles_pu = n.loads_t.p_set[boiler_loads].div(
        n.loads_t.p_set[boiler_loads].max(), axis=1
    )
    boiler_profiles_pu.columns = decentral_boilers
    boiler_profiles = DataArray(
        boiler_profiles_pu.multiply(n.links.loc[decentral_boilers, "p_nom"], axis=1)
    )

    # will be per unit
    n.model.add_variables(coords=[decentral_boilers], name="Link-fixed_profile_scaling")

    lhs = (
        (1, n.model["Link-p"].loc[:, decentral_boilers]),
        (
            -boiler_profiles,
            n.model["Link-fixed_profile_scaling"],
        ),
    )

    n.model.add_constraints(lhs, "=", 0, "Link-fixed_profile_scaling")

    # hack so that PyPSA doesn't complain there is nowhere to store the variable
    n.links["fixed_profile_scaling_opt"] = 0.0


def add_h2_derivate_limit(n, investment_year, limits_volume_max):
    for ct in limits_volume_max["h2_derivate_import"]:
        limit = limits_volume_max["h2_derivate_import"][ct][investment_year] * 1e6

        logger.info(f"limiting H2 derivate imports in {ct} to {limit / 1e6} TWh/a")

        incoming = n.links.loc[
            [
                "EU renewable oil -> DE oil",
                "EU methanol -> DE methanol",
                "EU renewable gas -> DE gas",
            ]
        ].index
        outgoing = n.links.loc[
            [
                "DE renewable oil -> EU oil",
                "DE methanol -> EU methanol",
                "DE renewable gas -> EU gas",
            ]
        ].index

        carrier_idx_dict = {
            # Every carrier should respect the limit individually
            "renewable_oil": 0,
            "methanol": 1,
            "renewable_gas": 2,
            # Exports of one carrier should not compensate for imports of another carrier
            "H2_derivate_oil_meoh": [0, 1],
            "H2_derivate_oil_gas": [0, 2],
            "H2_derivate_meoh_gas": [1, 2],
            # The sum of all carriers should respect the limit
            "H2_derivate_oil_meoh_gas": [0, 1, 2],
        }
        for carrier, idx in carrier_idx_dict.items():
            cname = f"{carrier}_import_limit-{ct}"

            incoming_p = (
                n.model["Link-p"].loc[:, incoming[idx]]
                * n.snapshot_weightings.generators
            ).sum()
            outgoing_p = (
                n.model["Link-p"].loc[:, outgoing[idx]]
                * n.snapshot_weightings.generators
            ).sum()

            lhs = incoming_p - outgoing_p

            n.model.add_constraints(lhs <= limit, name=f"GlobalConstraint-{cname}")

            if cname in n.global_constraints.index:
                logger.warning(
                    f"Global constraint {cname} already exists. Dropping and adding it again."
                )
                n.global_constraints.drop(cname, inplace=True)

            n.add(
                "GlobalConstraint",
                cname,
                constant=limit,
                sense="<=",
                type="",
                carrier_attribute="",
            )

    # Export bans on efuels are implemented in modify_prenetwork by restricting p_max_pu of the DE -> EU links


def adapt_nuclear_output(n):
    logger.info(
        "limiting german electricity generation from nuclear to 2020 value of 61 TWh"
    )
    limit = 61e6

    nuclear_de_index = n.links.index[
        (n.links.carrier == "nuclear") & (n.links.index.str[:2] == "DE")
    ]

    nuclear_gen = (
        n.model["Link-p"].loc[:, nuclear_de_index]
        * n.links.loc[nuclear_de_index, "efficiency"]
        * n.snapshot_weightings.generators
    ).sum()

    lhs = nuclear_gen

    cname = "Nuclear_generation_limit-DE"

    n.model.add_constraints(lhs <= limit, name=f"GlobalConstraint-{cname}")

    if cname in n.global_constraints.index:
        logger.warning(
            f"Global constraint {cname} already exists. Dropping and adding it again."
        )
        n.global_constraints.drop(cname, inplace=True)

    n.add(
        "GlobalConstraint",
        cname,
        constant=limit,
        sense="<=",
        type="",
        carrier_attribute="",
    )


def add_empty_co2_atmosphere_store_constraint(n):
    """
    Ensures that the CO2 atmosphere store at the last snapshot is empty.
    """
    logger.info(
        "Adding constraint for empty CO2 atmosphere store at the last snapshot."
    )
    cname = "empty_co2_atmosphere_store"

    last_snapshot = n.snapshots.values[-1]
    lhs = n.model["Store-e"].loc[last_snapshot, "co2 atmosphere"]

    n.model.add_constraints(lhs == 0, name=cname)


def additional_functionality(n, snapshots, snakemake):
    logger.info("Adding Ariadne-specific functionality")

    investment_year = int(snakemake.wildcards.planning_horizons[-4:])
    constraints = snakemake.params.solving["constraints"]

    if not snakemake.params.get("regret_run"):
        add_capacity_limits(
            n, investment_year, constraints["limits_capacity_min"], snakemake, "minimum"
        )

        add_capacity_limits(
            n, investment_year, constraints["limits_capacity_max"], snakemake, "maximum"
        )

        if snakemake.wildcards.clusters != "1":
            h2_import_limits(n, investment_year, constraints["limits_volume_max"])

            electricity_import_limits(
                n, investment_year, constraints["limits_volume_max"]
            )

        if investment_year >= 2025:
            h2_production_limits(
                n,
                investment_year,
                constraints["limits_volume_min"],
                constraints["limits_volume_max"],
            )
        add_h2_derivate_limit(n, investment_year, constraints["limits_volume_max"])

        if isinstance(constraints["co2_budget_national"], dict):
            add_national_co2_budgets(
                n,
                snakemake,
                constraints["co2_budget_national"],
                investment_year,
            )
        else:
            logger.warning("No national CO2 budget specified!")

    add_power_limits(n, investment_year, constraints["limits_power_max"])

    force_boiler_profiles_existing_per_boiler(n)

    if isinstance(constraints.get("decentral_heat_pump_budgets"), dict):
        add_decentral_heat_pump_budgets(
            n,
            constraints["decentral_heat_pump_budgets"],
            investment_year,
        )

    if investment_year == 2020:
        adapt_nuclear_output(n)

    if "co2 atmosphere" in n.generators.index:
        logger.warning(
            "CO2 atmosphere generator found. Adding constraint for empty CO2 atmosphere store at the last snapshot."
        )
        add_empty_co2_atmosphere_store_constraint(n)
