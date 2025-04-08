import marimo

__generated_with = "0.11.28"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import pypsa

    return mo, pd, pypsa


@app.cell
def _(pypsa):
    n_2020 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2020.nc"
    )
    n_2030 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2030.nc"
    )
    n_2040 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2040.nc"
    )
    n_2045 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2045.nc"
    )
    return n_2020, n_2030, n_2040, n_2045


@app.cell
def _(n_2020, n_2030, n_2040, n_2045, pd):
    # create a quick DataFrame
    df = pd.DataFrame()
    df["inst_cap_2020"] = n_2020.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2020"] = n_2020.statistics.optimal_capacity(comps="Store").round(2)

    df["inst_cap_2030"] = n_2030.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2030"] = n_2030.statistics.optimal_capacity(comps="Store").round(2)

    df["inst_cap_2040"] = n_2040.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2040"] = n_2040.statistics.optimal_capacity(comps="Store").round(2)
    df["inst_cap_2045"] = n_2045.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2045"] = n_2045.statistics.optimal_capacity(comps="Store").round(2)
    return (df,)


@app.cell
def _(df):
    df
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        for some stores the transition from optimal capacity of year a is not used as installed capacity in year a + 1.
        Compare "H2 Store" between 2020 and 2045 and "gas" in the same time span.
        Supply and demand are calculated using the optimal capacity as a maximum.
        """
    )
    return


@app.cell
def _(n_2045, pd):
    def get_location(n: n_2045, c: str, port: str = "") -> pd.Series:
        """"""
        bus = f"bus{port}"
        return n.static(c)[bus].map(n.buses.location).rename("location")

    from pypsa.statistics import groupers

    groupers.add_grouper("location", get_location)
    return get_location, groupers


@app.cell
def _(n_2045):
    n_2045.statistics.withdrawal(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    )
    return


@app.cell
def _(n_2045):
    n_2045.statistics.supply(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""Capax, tho, are calculated using difference in installed and optimized capacity. (opt_cap - installed_cap)*investment cost per MW."""
    )
    return


@app.cell
def _(n_2020, n_2030, pd):
    df_capex = pd.DataFrame()
    df_capex["instal_cap_GR_2020"] = (
        n_2020.statistics.installed_capacity(
            groupby=["location", "carrier", "bus_carrier"], comps="Store"
        )
        .loc["GR"]
        .loc["gas"]
    )
    df_capex["opt_cap_GR_2020"] = (
        n_2020.statistics.optimal_capacity(
            groupby=["location", "carrier", "bus_carrier"], comps="Store"
        )
        .loc["GR"]
        .loc["gas"]
    )
    df_capex["capex_GR_2020"] = (
        n_2020.statistics.capex(
            groupby=["location", "carrier", "bus_carrier"], comps="Store"
        )
        .loc["GR"]
        .loc["gas"]
    )
    df_capex["instal_cap_GR_2030"] = (
        n_2030.statistics.installed_capacity(
            groupby=["location", "carrier", "bus_carrier"], comps="Store"
        )
        .loc["GR"]
        .loc["gas"]
    )
    df_capex["opt_cap_GR_2030"] = (
        n_2030.statistics.optimal_capacity(
            groupby=["location", "carrier", "bus_carrier"], comps="Store"
        )
        .loc["GR"]
        .loc["gas"]
    )
    df_capex["capex_GR_2030"] = (
        n_2030.statistics.capex(
            groupby=["location", "carrier", "bus_carrier"], comps="Store"
        )
        .loc["GR"]
        .loc["gas"]
    )

    return (df_capex,)


@app.cell
def _(df_capex):
    df_capex
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
