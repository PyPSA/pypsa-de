"""Plot pyam IAMC graphs for all regions."""

from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_iamc_variables",
            simpl="",
            clusters="adm",
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

    vars = pd.read_excel(snakemake.input.iamc_variables, sheet_name="data")
    with Path(snakemake.output.sankey).open("w+") as fh:
        fh.write("touched")
