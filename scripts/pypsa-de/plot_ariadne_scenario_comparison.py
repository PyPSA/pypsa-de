# -*- coding: utf-8 -*-

import os

import matplotlib.pyplot as plt
import pandas as pd

from scripts._helpers import mock_snakemake


def scenario_plot(df, var):
    unit = df._get_label_or_level_values("Unit")[0]
    if var.startswith("Investment"):
        unit = "billion EUR2020/yr"
    df = df.droplevel("Unit")
    ax = df.T.plot(xlabel="years", ylabel=str(unit), title=str(var))
    plt.close()
    prefix = snakemake.config["run"]["prefix"]
    var = var.replace("|", "-").replace("\\", "-").replace(" ", "-").replace("/", "-")
    ax.figure.savefig(f"results/{prefix}/ariadne_comparison/{var}", bbox_inches="tight")


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "ariadne_all",
            # simpl="",
            # clusters=22,
            # opts="",
            # ll="vopt",
            # sector_opts="None",
            # planning_horizons="2050",
            # run="KN2045_Bal_v4"
        )

    dfs = []
    for file in snakemake.input.exported_variables:
        _df = pd.read_excel(
            file, index_col=list(range(5)), sheet_name="data"
        ).droplevel(["Model", "Region"])
        dfs.append(_df)

    df = pd.concat(dfs, axis=0)

    prefix = snakemake.config["run"]["prefix"]
    if not os.path.exists(f"results/{prefix}/ariadne_comparison/"):
        os.mkdir(f"results/{prefix}/ariadne_comparison/")

    for var in df._get_label_or_level_values("Variable"):
        scenario_plot(df.xs(var, level="Variable"), var)

    var = "Price|Carbon"
