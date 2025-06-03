"""Export pyam variables for all regions."""

import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

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

    with pd.ExcelWriter(snakemake.output.exported_variables_full) as writer:
        df.round(5).to_excel(writer, sheet_name="data", index=False)

    meta = pd.Series(
        {
            "Model": "PyPSA-AT",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "draft",
            "Release for publication": "no",
        }
    )
    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        df.round(5).to_excel(writer, sheet_name="data", index=False)
        meta.to_frame().T.to_excel(writer, sheet_name="meta", index=False)
