"""
Export variables in IAMC data model for all regions.

Notes
-----
https://docs.ece.iiasa.ac.at/iamc.html
https://pyam-iamc.readthedocs.io/en/stable/
"""

import pandas as pd
from pyam import IamDataFrame

from evals.fileio import read_networks
from evals.utils import insert_index_level
from scripts._helpers import configure_logging, mock_snakemake


def get_economy(n, year):
    """Extract total energy system cost for all regions."""

    # Nodal OPEX plus nodal Capex in billion EUR2020
    system_cost = (
        n.statistics.capex(groupby="location", aggregate_across_components="sum")
        .add(n.statistics.opex(groupby="location", aggregate_across_components="sum"))
        .div(1e9)
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
    """Extract all primary energy variables from the networks."""


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
    # networks = [pypsa.Network(n) for n in snakemake.input.networks]
    networks = read_networks(snakemake.input.results_path)

    iamc_variables = []
    for year, n in networks.items():
        iamc_variables.append(get_economy(n, year))
        iamc_variables.append(get_primary_energy(n, year))

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
