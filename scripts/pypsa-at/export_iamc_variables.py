"""
Export variables in IAMC data model for all regions.

IAMC variable naming convention:
Category|Subcategory|Specification

Notes
-----
https://docs.ece.iiasa.ac.at/iamc.html
https://pyam-iamc.readthedocs.io/en/stable/
"""

import pandas as pd
import pypsa
from pyam import IamDataFrame
from pypsa.statistics import port_efficiency

from evals.constants import TradeTypes
from evals.fileio import read_networks
from evals.utils import filter_by, insert_index_level
from scripts._helpers import configure_logging, mock_snakemake


def get_nodal_system_cost(n: pypsa.Network, year: str) -> pd.DataFrame:
    """
    Extract total energy system cost for all regions.

    Parameters
    ----------
    n
        The pypsa network instance.
    year
        The planning horizon for the network.

    Returns
    -------
    :
        Total CAPEX plus OPEX per model region in billions EUR (2020).
    """

    # Nodal OPEX plus nodal Capex in billion EUR2020
    system_cost = (
        n.statistics.capex(groupby="location", aggregate_across_components="sum")
        .add(n.statistics.opex(groupby="location", aggregate_across_components="sum"))
        .div(1e9)
        # centralize all below
        .round(4)
        .pipe(insert_index_level, "billion EUR2020", "Unit")
        .pipe(insert_index_level, "Cost|Total Energy System Cost", "Variable")
        .pipe(insert_index_level, snakemake.wildcards.run, "Scenario")
        .pipe(insert_index_level, "PyPSA-AT", "Model")
        .pipe(insert_index_level, year, "Year", pos=-1)
    )
    system_cost.index = system_cost.index.rename({"location": "Region"})

    return system_cost


def get_primary_energy(n, year) -> pd.DataFrame:
    """
    Extract all primary energy variables from the networks.

    In general, primary energy is the supply side of
    network components. If a component has both, supply and demand,
    the balance is returned and a warning is raised.

    Variables for primary energy follow the naming scheme:
    Primary Energy|BusCarrierGroup|Technology

    Parameters
    ----------
    n
        The pypsa network instance.
    year
        The planning horizon for the network.

    Returns
    -------
    :
        The primary energy in MWh for all regions and
    """

    var = {}

    # "Primary Energy|Oil"  (= sum of all oil variables)
    # "Primary Energy|Oil|Heat"
    # "Primary Energy|Oil|Electricity"
    # "Primary Energy|Oil|Liquids" (= Oil minus rest)
    # what about Fischer-Tropsch?

    # need to reduce primary energy by renewable oil production in the same region
    # renewable_oil_production = n.statistics.supply(
    #     groupby=["location", "carrier", "bus_carrier"], comps="Link", bus_carrier="oil"
    # ).drop(
    #     ["unsustainable bioliquids", "oil refining", "import oil"],
    #     axis=0,
    #     level="carrier",
    #     errors="ignore",
    # )
    # oil_usage = (
    #     n.statistics.withdrawal(
    #         groupby=["location", "carrier", "bus_carrier"], bus_carrier="oil"
    #     )
    #     .drop("Store", axis=0, level="component", errors="ignore")
    #     .droplevel("component")
    #     .mul(-1)  # withdrawal is negative
    # )
    #
    # # assuming that regional oil production is consumed in the same region
    # oil_saldo = (
    #     oil_usage.groupby("location")
    #     .sum()
    #     .add(renewable_oil_production.groupby("location").sum(), fill_value=0)
    # )
    # oil_import = oil_saldo[oil_saldo.le(0)]
    # oil_export = oil_saldo[oil_saldo.gt(0)]
    #
    # var["Primary Energy|Oil"] = oil_import.divide(oil_refining_efficiency)
    # var["Final Energy|Oil"] = oil_export

    oil_balance = n.statistics.energy_balance(
        groupby="location", bus_carrier="oil", comps="Link"
    ).drop("EU", errors="ignore")
    liquids_all = oil_balance[oil_balance.le(0)].abs()

    # increase fossil oil demand by this factor to account for refining losses
    oil_refining_efficiency = (
        port_efficiency(n, "Link", "1").filter(like="oil refining").item()
    )
    # calculate the split of fossil oil generation to liquids production and assume
    # the same split for all regions
    oil_fossil_eu = n.statistics.supply(bus_carrier="oil primary").item()
    renewable_liquids = n.statistics.supply(
        groupby="bus_carrier", bus_carrier="oil", comps="Link"
    ).item()
    fossil_fraction = liquids_all * oil_fossil_eu / renewable_liquids
    # increase fossil oil share by refining losses
    liquids_w_losses = (
        liquids_all - fossil_fraction + fossil_fraction / oil_refining_efficiency
    )
    # assuming that regional oil production is consumed in the same region, the netted
    # energy balance becomes the primary energy (import) and final energy (export) of
    # the region. It is important to note that this is not the same as fossil oil imports,
    # but all oil imports including renewable oil production from other regions.
    # Lacking an oil network that connects regions, we assume that all oil is consumed
    # locally and only oil production surplus becomes exported as final energy.
    var["Primary Energy|Liquids"] = liquids_w_losses
    var["Primary Energy|Liquids|oil refining losses"] = liquids_w_losses - liquids_all

    # TEST: oil import must be same as total sum of oil primary
    # pd.testing.assert_series_equal(
    #     oil_import, oil_energy_balance[oil_energy_balance.le(0)], check_names=False
    # )  # PASS -> its the same

    # "Primary Energy|Gas"
    # "Primary Energy|Gas|Heat"
    # "Primary Energy|Gas|Electricity"
    # "Primary Energy|Gas|Hydrogen"
    # "Primary Energy|Gas|Gases" (?) Sabatier?

    for scope in (TradeTypes.DOMESTIC, TradeTypes.FOREIGN):
        gas_trade = n.statistics.trade_energy(
            bus_carrier="gas", direction="import", scope=scope
        )
        gas_trade = gas_trade[gas_trade.gt(0)].groupby("location").sum()
        var[f"Primary Energy|Gas|Import {scope.title()}"] = gas_trade

    # gas_trade_foreign = n.statistics.trade_energy(
    #     bus_carrier="gas", direction="import", scope="foreign"
    # )
    # gas_trade_domestic = n.statistics.trade_energy(
    #     bus_carrier="gas", direction="import", scope="domestic"
    # )
    #
    # gas_trade_foreign = (
    #     gas_trade_foreign[gas_trade_foreign.gt(0)].groupby("location").sum()
    # )
    # gas_trade_domestic = (
    #     gas_trade_domestic[gas_trade_domestic.gt(0)].groupby("location").sum()
    # )

    gas_generation = n.statistics.supply(
        groupby=["location", "carrier"], bus_carrier="gas", comps="Generator"
    )
    # var["Primary Energy|Gas|Import Foreign"] = gas_trade_foreign
    # var["Primary Energy|Gas|Import Domestic"] = gas_trade_domestic
    var["Primary Energy|Gas|Production"] = (
        filter_by(gas_generation, carrier="gas").groupby("location").sum()
    )
    var["Primary Energy|Gas|Import Global"] = (
        filter_by(gas_generation, carrier="import gas").groupby("location").sum()
    )
    # /IdeaProjects/pypsa-at/resources/v2025.02/KN2045_Mix/gas_input_locations_s_adm_simplified.csv
    var["Primary Energy|Gas"] = (
        pd.concat(
            [
                var[f"Primary Energy|Gas|Import {TradeTypes.FOREIGN}"],
                var[f"Primary Energy|Gas|Import {TradeTypes.DOMESTIC}"],
                gas_generation,
            ]
        )
        .groupby("location")
        .sum()
    )

    # non-sequestered HVC + municipal solid waste
    # municipal solid waste is generated regionally
    # municipal solid waste supplies to HVC bus + draws from CO2 atmosphere
    # naphtha for industry withdraws from oil and supplies naphtha + waste to HVC bus + process emissions
    # waste CHP only draws from HVC bus
    # Primary Energy|Waste is only municipal solid waste Generators, the rest is secondary energy
    # plus 'municipal solid waste transport' import amounts
    # n.statistics.supply(groupby=["location", "carrier", "bus_carrier"], bus_carrier=["non-sequestered HVC", "municipal solid waste"])
    waste_generation = n.statistics.supply(
        groupby="location", bus_carrier="municipal solid waste", comps="Generator"
    )
    waste_import = n.statistics.trade_energy(
        bus_carrier="municipal solid waste",
        direction="import",
        scope=("domestic", "foreign"),
    )
    waste_import = waste_import[waste_import.gt(0)].groupby("location").sum()

    var["Primary Energy|Waste"] = waste_generation.add(waste_import, fill_value=0)
    var["Primary Energy|Waste|Import"] = waste_import
    # "Primary Energy|Waste"  (= non-sequestered HVC)
    # "Primary Energy|Waste|Heat"
    # "Primary Energy|Waste|Electricity"

    # similar to oil, but without refining losses
    # Primary Energy|Coal
    # Primary Energy|Coal|Hard Coal
    # Primary Energy|Coal|Lignite
    # Primary Energy|Coal|Heat
    # Primary Energy|Coal|Electricity

    # Primary Energy|Fossil (= Coal + Gas + Oil + non-renewable HVC)

    # Primary Energy|Biomass
    # Primary Energy|Biomass|Gases  (BioMethane, BioSNG)
    # Primary Energy|Biomass|Liquids
    # Primary Energy|Biomass|Electricity
    # Primary Energy|Biomass|Heat
    # Primary Energy|Biomass|Solids
    # Primary Energy|Biomass|Hydrogen
    # Primary Energy|Biomass|Methanol

    # Primary Energy|Ammonia

    # Primary Energy|Nuclear

    # Primary Energy|Hydro
    # Primary Energy|Hydro|Pumped Storage
    # Primary Energy|Hydro|Run-of-River
    # Primary Energy|Solar
    # Primary Energy|Solar|PV-Rooftop
    # Primary Energy|Solar|PV-Utility
    # Primary Energy|Solar|PV-HSAT
    # Primary Energy|Wind
    # Primary Energy|Wind|Onshore
    # Primary Energy|Wind|Offshore

    # Primary Energy|Heat|Solar-Thermal
    # Primary Energy|Heat|Ambient Heat  (heat pumps)
    # Primary Energy|Heat|Latent Heat  (CHPs with efficiency > 1)
    # Primary Energy|Heat|Geothermal  (Geothermic heat)

    # Primary Energy|Renewable (= Biomass + Hydro + Solar + Wind + renewable HVC)

    # EU and IEA statistics tend to report imported electricity and fuels as
    # contributing to the region’s primary energy supply.
    # Primary Energy|Import|Electricity
    # Primary Energy|Import|Oil
    # Primary Energy|Import|Gas
    # Primary Energy|Import|Coal  (Hard + Lignite)
    # Primary Energy|Import|Waste
    # Primary Energy|Import|Biomass
    # Primary Energy|Import|Hydrogen


def get_secondary_energy(n, year) -> pd.DataFrame:
    """Extract all secondary energy variables from the networks."""


def get_final_energy(n, year) -> pd.DataFrame:
    """Extract all final energy variables from the networks."""


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_iamc_variables",
            # simpl="",
            # clusters="adm",
            # opts="",
            # ll="vopt",
            # sector_opts="None",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

    template = pd.read_excel(snakemake.input.template, sheet_name="data")

    # want to use 2050 first, to catch all technologies
    # during development and debugging sessions
    networks = read_networks(sorted(snakemake.input.networks, reverse=True))

    iamc_variables = []
    for year, n in networks.items():
        iamc_variables.append(get_nodal_system_cost(n, year))
        iamc_variables.append(get_primary_energy(n, year))
        iamc_variables.append(get_secondary_energy(n, year))
        iamc_variables.append(get_final_energy(n, year))

    df = pd.concat(iamc_variables)

    # drop all values for dummy template
    dummy_data = {
        "Model": ["PyPSA-AT"],
        "Scenario": [snakemake.wildcards.run],
        "Region": ["Europe"],
        "Variable": ["Category|Subcategory|Specification"],
        "Unit": ["Snuckels"],
        2020: [1],
        2025: [1],
        2030: [1],
        2035: [1],
        2040: [1],
        2045: [1],
        2050: [1],
    }
    new_row = pd.DataFrame(dummy_data)
    df = pd.concat([template, new_row])
    iamc = IamDataFrame(df)

    with pd.ExcelWriter(snakemake.output.exported_variables_full) as writer:
        iamc.to_excel(writer, sheet_name="data", index=False)

    meta = pd.Series(
        {
            "Model": "PyPSA-AT",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "draft",
            "Release for publication": "no",
        }
    ).to_frame("value")
    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        iamc.to_excel(writer, sheet_name="data", index=False)
        meta.to_excel(writer, sheet_name="meta", index=False)
