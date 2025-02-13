# -*- coding: utf-8 -*-
import pyam

if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_ariadne_database")

    db = pyam.read_iiasa(
        snakemake.params.db_name,
        model=snakemake.params.leitmodelle.values(),
        scenario=snakemake.params.scenarios,
        # Download only the most recent iterations of scenarios
    )

    db.timeseries().to_csv(snakemake.output.data)
