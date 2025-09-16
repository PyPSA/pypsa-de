import os

import matplotlib

matplotlib.use("Agg")  # Use a non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd

from scripts._helpers import mock_snakemake


def scenario_plot(df, output_dir, var):
    unit = df._get_label_or_level_values("Unit")[0]
    if var.startswith("Investment"):
        unit = "billion EUR2020/yr"
    df = df.droplevel("Unit")
    ax = df.T.plot(xlabel="years", ylabel=str(unit), title=str(var))
    var = var.replace("|", "-").replace("\\", "-").replace(" ", "-").replace("/", "-")
    ax.figure.savefig(f"{output_dir}/{var}.png", bbox_inches="tight", dpi=100)
    plt.close(ax.figure)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_scenario_comparison_regrets",
            regret_dir="regret_networks",
            # simpl="",
            # clusters=22,
            # opts="",
            # ll="vopt",
            # sector_opts="None",
            # planning_horizons="2050",
            # run="KN2045_Mix"
        )

    dfs = []
    fns = snakemake.input.exported_variables
    if "regret_variables" in fns[0] and len(fns) == 4:
        # reorder indices of fns as 0312
        fns = [fns[i] for i in [0, 3, 2, 1] if i < len(fns)]
    if "regret_variables" in fns[0] and len(fns) in [9, 16]:
        fns = [
            fn for fn in fns if "NoFlex/" not in fn
        ]  # !!! CAVEAT DISPATCHING ON FILENAME
    for file in fns:
        _df = pd.read_excel(
            file, index_col=list(range(5)), sheet_name="data"
        ).droplevel(["Model", "Region"])
        dfs.append(_df)

    df = pd.concat(dfs, axis=0)

    output_dir = snakemake.params.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for var in df._get_label_or_level_values("Variable"):
        scenario_plot(df.xs(var, level="Variable"), output_dir, var)
