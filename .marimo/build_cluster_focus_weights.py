import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd

    return mo, pd


@app.cell
def _():
    countries = [
        "AL",
        "AT",
        "BA",
        "BE",
        "BG",
        "CH",
        "CZ",
        "DE",
        "DK",
        "EE",
        "ES",
        "FI",
        "FR",
        "GB",
        "GR",
        "HR",
        "HU",
        "IE",
        "IT",
        "LT",
        "LU",
        "LV",
        "ME",
        "MK",
        "NL",
        "NO",
        "PL",
        "PT",
        "RO",
        "RS",
        "SE",
        "SI",
        "SK",
        "XK",
    ]
    return (countries,)


@app.cell
def _(countries):
    nodes_per_country = {
        "DE": 16,
        "AT": 10,
        "IT": 3,
        "DK": 2,
        # "UK": 2,
        # "ES": 2,
        # "GR": 2
    }
    country_nodes = {c: nodes_per_country.get(c, 1) for c in countries}
    country_nodes
    return country_nodes, nodes_per_country


@app.cell
def _(country_nodes):
    n_cluster = sum(country_nodes.values())
    n_cluster
    return (n_cluster,)


@app.cell
def _(country_nodes, n_cluster):
    focus_weights = {
        c: round(w / n_cluster - 5e-5, 4) for c, w in country_nodes.items()
    }
    focus_weights
    return (focus_weights,)


@app.cell
def _(focus_weights):
    focus_weights["AT"] += 1 - sum(focus_weights.values())
    focus_weights["AT"] = round(focus_weights["AT"], 4)
    focus_weights
    return


@app.cell
def _(focus_weights):
    assert sum(focus_weights.values()) == 1.0, (
        f"Sum of focus weights is not 1.0 but {sum(focus_weights.values())}"
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
