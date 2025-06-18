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

from evals.fileio import read_networks
from evals.iamcvars.primary import (
    primary_ammonia,
    primary_biomass,
    primary_coal,
    primary_fossil_oil,
    primary_gas,
    primary_heat,
    primary_hydro,
    primary_hydrogen,
    primary_nuclear,
    primary_solar,
    primary_waste,
    primary_wind,
)
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


def collect_primary_energy(n, year) -> pd.DataFrame:
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

    var = primary_fossil_oil(n, {})
    var = primary_gas(n, var)
    var = primary_hydrogen(n, var)
    var = primary_waste(n, var)
    var = primary_coal(n, var)
    var = primary_biomass(n, var)
    var = primary_hydro(n, var)
    var = primary_solar(n, var)
    var = primary_nuclear(n, var)
    var = primary_ammonia(n, var)
    var = primary_wind(n, var)
    var = primary_heat(n, var)

    # "Primary Energy|Gas"
    # "Primary Energy|Gas|Heat"
    # "Primary Energy|Gas|Electricity"
    # "Primary Energy|Gas|Hydrogen"
    # "Primary Energy|Gas|Gases" (?) Sabatier?

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

    return var


def collect_secondary_energy(n, year) -> pd.DataFrame:
    """Extract all secondary energy variables from the networks."""


def collect_final_energy(n, year) -> pd.DataFrame:
    """Extract all final energy variables from the networks."""


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_iamc_variables",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

    # template = pd.read_excel(snakemake.input.template, sheet_name="data")

    # want to use 2050 first, to catch all technologies
    # during development and debugging sessions
    networks = read_networks(sorted(snakemake.input.networks, reverse=True))

    iamc_variables = []
    for year, n in networks.items():
        iamc_variables.append(get_nodal_system_cost(n, year))
        iamc_variables.append(collect_primary_energy(n, year))
        iamc_variables.append(collect_secondary_energy(n, year))
        iamc_variables.append(collect_final_energy(n, year))

    df = pd.concat(iamc_variables)

    # # drop all values for dummy template
    iamc = IamDataFrame(df)

    meta = pd.Series(
        {
            "Model": "PyPSA-AT",
            "Commit SHA": "",
            "Repository": "https://gitlab.aggm.at/philip.worschischek/pypsa-at",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "draft",
            "Release for publication": "no",
        }
    ).to_frame("value")
    meta.index.name = "key"
    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        iamc.to_excel(writer, sheet_name="data", index=False)
        meta.to_excel(writer, sheet_name="meta", index=False)
