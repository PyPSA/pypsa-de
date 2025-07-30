import marimo

__generated_with = "0.14.9"
app = marimo.App(width="medium")


@app.cell
def _(mod):
    mod.md()
    return


@app.cell
def _(mo, radio):
    mo.vstack(["Year:", radio])
    return


@app.cell(hide_code=True)
def _(fig):
    fig.show()
    return


@app.cell
def _(mo, networks):
    options = sorted(networks.keys())
    radio = mo.ui.radio(options=options, inline=True)
    return (radio,)


@app.cell
def _():
    import sys

    sys.path.insert(0, "..")

    import marimo as mo
    import plotly.graph_objects as go

    return go, mo


@app.cell
def _():
    from pathlib import Path

    results_path = Path("/IdeaProjects/pypsa-at/results/v2025.03/AT10_KN2040")
    return (results_path,)


@app.cell
def _(networks, radio):
    n = networks[radio.value]
    return (n,)


@app.cell
def _(results_path):
    from evals.fileio import read_networks

    networks = read_networks(result_path=results_path)
    return (networks,)


@app.cell
def _(n):
    df = n.statistics.optimal_capacity(
        groupby=["location", "carrier"], comps=["Generator", "Link"]
    )
    return (df,)


@app.cell
def _(df):
    from evals.utils import rename_aggregate

    # df_clean = df.drop(ignore_carrier, level="carrier", errors="ignore")
    df_clean = df.filter(regex="wind|solar|H2 Electrolysis|nuclear|ror|hydro", axis=0)
    solar_thermals = df_clean.filter(like="thermal", axis=0)
    solar_thermals
    df_clean.drop(solar_thermals.index, inplace=True, axis=0)
    print(dict.fromkeys(df_clean.index.unique("carrier"), ""))
    mapper = {
        "solar": "solar",
        "solar rooftop": "solar",
        "onwind": "wind",
        "solar-hsat": "solar",
        "offwind-ac": "wind",
        "offwind-dc": "wind",
        "H2 Electrolysis": "electrolysis",
        "nuclear": "nuclear",
        "ror": "hydro",
        "hydro": "hydro",
    }
    df_clean = rename_aggregate(df_clean, mapper)
    # df_clean
    return (df_clean,)


@app.cell
def _(df_clean, go):
    from plotly.subplots import make_subplots

    plot_df = df_clean.to_frame("values").pivot_table(
        values="values", index="location", columns="carrier", aggfunc="sum"
    )
    # nplots = len(plot_df.columns)-1
    fig = make_subplots(
        rows=2,
        cols=2,
        shared_yaxes=True,
    )
    x = plot_df["electrolysis"].to_numpy()
    plot_df.drop("electrolysis", inplace=True, axis=1)

    for i, column in enumerate(plot_df, start=1):
        y = plot_df[column]

        _idx = {
            1: (1, 1),
            2: (2, 1),
            3: (1, 2),
            4: (2, 2),
        }

        row, col = _idx[i]

        # model = LinearRegression()
        # model.fit(x, y)
        # y_pred = model.predict(x)

        fig.add_trace(
            go.Scatter(x=x, y=y, mode="markers+text", name=column, text=y.index),
            row=row,
            col=col,
        )
        fig.update_xaxes(title_text=column, row=row, col=col)

    #     fig = px.scatter(x=y, y=y_pred, labels={'x': 'ground truth', 'y': 'prediction'})

    # color_discrete_map = {
    #     "H2 Electrolysis": "grey",
    #     "nuclear": "orange",
    #     "offwind-ac": "blue",
    #     "offwind-dc": "blue",
    #     "onwind": "blue",
    #     "solar": "yellow",
    #     "solar rooftop": "yellow",
    #     "solar-hsat": "yellow",
    # }

    # fig = px.scatter_matrix(
    #     plot_df,
    #     dimensions=plot_df.columns,
    #     color_discrete_map=color_discrete_map
    # )
    fig.update_yaxes(title_text="Electrolysis", row=1, col=1)
    fig.update_yaxes(title_text="Electrolysis", row=2, col=1)
    fig.update_layout(height=1000, width=1000, title_text="Electrolysis vs Renewables")
    return (fig,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
