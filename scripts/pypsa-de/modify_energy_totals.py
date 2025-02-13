# -*- coding: utf-8 -*-


import pandas as pd

energy_totals = pd.read_csv(snakemake.input.energy_totals, index_col=0)

ariadne = pd.read_csv(
    snakemake.input.ariadne,
    index_col=["model", "scenario", "region", "variable", "unit"],
)

energy_totals.to_csv(snakemake.output.energy_totals)
