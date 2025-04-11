import marimo

__generated_with = "0.12.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys

    sys.path.insert(0, "/IdeaProjects/pypsa-at")

    import marimo as mo

    from evals.fileio import prepare_co2_emissions, read_networks

    return mo, prepare_co2_emissions, read_networks, sys


@app.cell
def _(read_networks):
    networks = read_networks("/IdeaProjects/pypsa-at/results", "evals-dev/networks")
    return (networks,)


@app.cell
def _(networks):
    n = networks["2030"]
    return (n,)


@app.cell
def _():
    from scripts._helpers import (
        copy_default_files,
        get_rdir,
        get_scenarios,
        get_shadow,
        path_provider,
    )

    return (
        copy_default_files,
        get_rdir,
        get_scenarios,
        get_shadow,
        path_provider,
    )


@app.cell
def _(get_rdir, n):
    get_rdir(n.meta["run"])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
