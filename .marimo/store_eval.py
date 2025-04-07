import marimo

__generated_with = "0.11.28"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import pypsa

    return pd, pypsa


@app.cell
def _(pypsa):
    n_2020 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2020.nc"
    )
    return (n_2020,)


@app.cell
def _(pypsa):
    n_20245 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2045.nc"
    )
    return (n_20245,)


@app.cell
def _(n_20245):
    n_20245.stores[
        n_20245.stores.bus.str.contains("AT")
        & n_20245.stores.carrier.str.contains("H2")
    ]
    return


@app.cell
def _(n_2020):
    n_2020.stores[
        n_2020.stores.bus.str.contains("AT") & n_2020.stores.carrier.str.contains("gas")
    ]
    return


@app.cell
def _(n_2020):
    n_2020.statistics.installed_capacity(comps="Store")
    return


@app.cell
def _(n_20245):
    n_20245.statistics.installed_capacity(comps="Store")
    return


@app.cell
def _(n_20245):
    n_20245.statistics.optimal_capacity(comps="Store")
    return


@app.cell
def _(pypsa):
    n_2030 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2030.nc"
    )
    n_2040 = pypsa.Network(
        "/mnt/storage/pypsa-at-AT10-365H/networks/base_s_adm__none_2040.nc"
    )
    return n_2030, n_2040


@app.cell
def _(n_2020, n_20245, n_2030, n_2040, pd):
    df = pd.DataFrame()
    df["inst_cap_2020"] = n_2020.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2020"] = n_2020.statistics.optimal_capacity(comps="Store").round(2)

    df["inst_cap_2030"] = n_2030.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2030"] = n_2030.statistics.optimal_capacity(comps="Store").round(2)

    df["inst_cap_2040"] = n_2040.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2040"] = n_2040.statistics.optimal_capacity(comps="Store").round(2)
    df["inst_cap_2045"] = n_20245.statistics.installed_capacity(comps="Store").round(2)
    df["Opt_cap_2045"] = n_20245.statistics.optimal_capacity(comps="Store").round(2)

    return (df,)


@app.cell
def _(df):
    df
    return


@app.cell
def _(n_20245):
    n_20245.stores_t
    return


@app.cell
def _(n_20245, pd):
    def get_location(n: n_20245, c: str, port: str = "") -> pd.Series:
        """"""
        bus = f"bus{port}"
        return n.static(c)[bus].map(n.buses.location).rename("location")

    from pypsa.statistics import groupers

    groupers.add_grouper("location", get_location)

    return get_location, groupers


@app.cell
def _(n_20245):
    n_20245.statistics.withdrawal(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    )
    return


@app.cell
def _(n_20245):
    n_20245.statistics.supply(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    )
    return


@app.cell
def _(n_2030):
    n_2030.static("Store").loc["GR gas Store"].filter(like="e_nom").T
    return


@app.cell
def _(n_2020):
    n_2020.static("Store").loc["GR gas Store"].filter(like="e_nom").T
    return


@app.cell
def _(n_2020):
    n_2020.statistics.capex(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    ).loc["GR"].loc["gas"]
    return


@app.cell
def _(n_2020):
    n_2020.statistics.installed_capacity(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    ).loc["GR"].loc["gas"]
    return


@app.cell
def _(n_2020):
    n_2020.statistics.optimal_capacity(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    ).loc["GR"].loc["gas"]
    return


@app.cell
def _(n_2030):
    n_2030.statistics.capex(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    ).loc["GR"].loc["gas"]
    return


@app.cell
def _(n_2030):
    n_2030.statistics.optimal_capacity(
        groupby=["location", "carrier", "bus_carrier"], comps="Store"
    ).loc["GR"].loc["gas"]
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
