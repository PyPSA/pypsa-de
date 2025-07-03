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
    primary_gas,
    primary_heat,
    primary_hydro,
    primary_hydrogen,
    primary_nuclear,
    primary_oil,
    primary_solar,
    primary_waste,
    primary_wind,
)
from evals.statistic import collect_myopic_statistics
from evals.utils import insert_index_level
from scripts._helpers import configure_logging, mock_snakemake


def combine_variables(var: dict, unit: str, year: str) -> pd.Series:
    """
    Combine variables into a single dataframe.

    Parameters
    ----------
    var
    unit
    year

    Returns
    -------
    :
    """
    # for name, ds in var.items():
    #     ds = insert_index_level(ds, unit, "Unit")
    #     ds = insert_index_level(ds, year, "Year")
    #     var[name] = ds
    ds = (
        pd.concat({k: v for k, v in var.items() if not v.empty})
        .pipe(insert_index_level, unit, "Unit")
        .pipe(insert_index_level, year, "Year")
    )
    ds.index = ds.index.rename({None: "Variable"})

    return ds


def collect_system_cost(n: pypsa.Network) -> pd.DataFrame:
    """
    Extract total energy system cost per region.

    Parameters
    ----------
    n
        The pypsa network instance.

    Returns
    -------
    :
        Total CAPEX plus OPEX per model region in billions EUR (2020).
    """
    # Nodal OPEX and nodal CAPEX in billion EUR2020
    var = {}
    var["System Costs|CAPEX"] = n.statistics.capex(
        groupby="location", aggregate_across_components="sum"
    ).div(1e9)
    var["System Costs|OPEX"] = n.statistics.opex(
        groupby="location", aggregate_across_components="sum"
    ).div(1e9)
    var["System Costs"] = var["System Costs|CAPEX"] + var["System Costs|OPEX"]

    # todo: enable, or test later
    # assert vars["System Costs"].mul(1e9).sum() == n.objective, "Total system costs do not match the optimization result."

    # system_cost = (
    #     n.statistics.capex(groupby="location", aggregate_across_components="sum")
    #     .add(n.statistics.opex(groupby="location", aggregate_across_components="sum"))
    #     .div(1e9)  # to Billion (Milliarden)
    #     # # centralize all below
    #     # .round(4)
    #     # .pipe(insert_index_level, "billion EUR2020", "Unit")
    #     # .pipe(insert_index_level, "Cost|Total Energy System Cost", "Variable")
    #     # .pipe(insert_index_level, snakemake.wildcards.run, "Scenario")
    #     # .pipe(insert_index_level, "PyPSA-AT", "Model")
    #     # .pipe(insert_index_level, year, "Year", pos=-1)
    # )
    # # system_cost.index = system_cost.index.rename({"location": "Region"})

    return combine_variables(var, "billion EUR2020", n.year)


def collect_primary_energy(n: pypsa.Network) -> pd.Series:
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

    Returns
    -------
    :
        The primary energy in MWh for all regions and
    """
    var = {}
    primary_oil(n, var)
    primary_gas(n, var)
    primary_hydrogen(n, var)
    primary_waste(n, var)
    primary_coal(n, var)
    primary_biomass(n, var)
    primary_hydro(n, var)
    primary_solar(n, var)
    primary_nuclear(n, var)
    primary_ammonia(n, var)
    primary_wind(n, var)
    primary_heat(n, var)

    return combine_variables(var, "MWh", n.year)


def collect_secondary_energy(n) -> pd.DataFrame:
    """Extract all secondary energy variables from the networks."""


def collect_final_energy(n) -> pd.DataFrame:
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

    kwargs = {
        "groupby": ["location", "carrier", "bus_carrier", "unit"],
        "aggregate_components": False,
        "drop_zeros": False,
        "drop_unit": False,
    }
    # myopic_energy_balance = collect_myopic_statistics(networks, "energy_balance", **kwargs)
    # calculate all statistics and process them to IAMC data model.
    myopic_supply = collect_myopic_statistics(networks, "supply", **kwargs)
    myopic_withdrawal = collect_myopic_statistics(networks, "withdrawal", **kwargs)
    myopic_opex = collect_myopic_statistics(networks, "opex", **kwargs)  # wrong unit
    myopic_capex = collect_myopic_statistics(networks, "capex", **kwargs)  # wrong unit
    # Idea: calculate all once and extract from global mutable series
    # to ensure nothing is forgotten. The global series however is a risk.

    iamc_variables = []
    for year, n in networks.items():
        iamc_variables.append(collect_system_cost(n))
        iamc_variables.append(collect_primary_energy(n))
        iamc_variables.append(collect_secondary_energy(n))
        iamc_variables.append(collect_final_energy(n))

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
