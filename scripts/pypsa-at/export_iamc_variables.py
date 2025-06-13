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

from evals.fileio import read_networks
from evals.utils import insert_index_level
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

    # increase oil demand by this factor to account for refining losses
    oil_refining_efficiency = (
        port_efficiency(n, "Link", "1").filter(like="oil refining").item()
    )
    # need to reduce primary energy by renewable oil production in this region by Fischer-Tropsch
    renewable_oil_production = n.statistics.supply(
        groupby=["location", "carrier", "bus_carrier"], comps="Link", bus_carrier="oil"
    ).drop(
        ["unsustainable bioliquids", "oil refining", "import oil"],
        axis=0,
        level="carrier",
        errors="ignore",
    )
    oil_usage = (
        n.statistics.withdrawal(
            groupby=["location", "carrier", "bus_carrier"], bus_carrier="oil"
        )
        .drop("Store", axis=0, level="component", errors="ignore")
        .droplevel("component")
        .mul(-1)  # withdrawal is negative
    )

    # assuming that regional oil production is consumed in the same region
    oil_saldo = (
        oil_usage.groupby("location")
        .sum()
        .add(renewable_oil_production.groupby("location").sum(), fill_value=0)
    )
    oil_import = oil_saldo[oil_saldo.gt(0)]
    oil_export = oil_saldo[oil_saldo.le(0)]

    var["Primary Energy|Oil"] = oil_import.divide(oil_refining_efficiency)
    var["Final Energy|Oil"] = oil_export

    oil_energy_balance = n.statistics.energy_balance(
        groupby="location", bus_carrier="oil", comps="Link"
    ).drop("EU", errors="ignore")
    var["Primary Energy|Oil"] = oil_energy_balance[oil_energy_balance.gt(0)]

    # TEST: oil import must be same as total sum of oil primary

    # "Primary Energy|Gas"
    # "Primary Energy|Gas|Heat"
    # "Primary Energy|Gas|Electricity"
    # "Primary Energy|Gas|Hydrogen"
    # "Primary Energy|Gas|Gases" (?) Sabatier?

    # "Primary Energy|Waste"  (= non-sequestered HVC)
    # "Primary Energy|Waste|Heat"
    # "Primary Energy|Waste|Electricity"

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
            simpl="",
            clusters="adm",
            opts="",
            ll="vopt",
            sector_opts="None",
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
