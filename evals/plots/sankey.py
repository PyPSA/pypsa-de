# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for Sankey diagram."""

import pandas as pd
import plotly
import pyam
from plotly.graph_objs import Figure, Sankey
from pyam.index import get_index_levels

from evals.utils import filter_by, rename_aggregate

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)


def sankey(df: pyam.IamDataFrame, mapping: dict) -> Figure:
    """
    Plot a sankey diagram.

    It is currently only possible to create this diagram for single years.

    Parameters
    ----------
    df
        Data to be plotted
    mapping
        Assigns the source and target component of a variable

        .. code-block:: python

            {
                variable: (source, target),
            }

    Returns
    -------
    :
        The generated plotly figure.
    """

    # Check for duplicates
    for col in [name for name in df.dimensions if name != "variable"]:
        levels = get_index_levels(df._data, col)
        if len(levels) > 1:
            raise ValueError(f"Non-unique values in column {col}: {levels}")

    # Concatenate the data with source and target columns
    _df = pd.DataFrame.from_dict(
        mapping, orient="index", columns=["source", "target"]
    ).merge(df._data, how="left", left_index=True, right_on="variable")
    label_mapping = {
        label: i
        for i, label in enumerate(set(pd.concat([_df["source"], _df["target"]])))
    }
    _df.replace(label_mapping, inplace=True)
    region = get_index_levels(_df, "region")[0]
    unit = get_index_levels(_df, "unit")[0]
    year = get_index_levels(_df, "year")[0]
    fig = Figure(
        data=[
            Sankey(
                valuesuffix=unit,
                node=dict(
                    pad=15,
                    thickness=10,
                    line=dict(color="black", width=0.5),
                    label=pd.Series(list(label_mapping)),
                    hovertemplate="%{label}: %{value}<extra></extra>",
                    color="blue",
                ),
                link=dict(
                    arrowlen=15,
                    source=_df.source,
                    target=_df.target,
                    value=_df.value,
                    hovertemplate='"%{source.label}" to "%{target.label}": \
                %{value}<extra></extra>',
                ),
            )
        ]
    )
    fig.update_layout(title_text=f"region: {region}, year: {year}", font_size=10)
    return fig


def read_iamc_data_frame(filepath):
    xls = pd.read_excel(
        filepath,
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    return xls.stack()


def get_mapping(df) -> (dict, set):
    mapping = {}
    nodes = set()
    for v in df.index.unique("Variable"):
        # skip aggregations
        if v.count("|") < 2:
            continue

        if v.startswith("Primary"):
            _, bus_carrier, tech = v.split("|")
            mapping[v] = (tech, bus_carrier)
            nodes.add(tech)
            nodes.add(bus_carrier)
        elif v.startswith("Secondary"):
            _, bc_output, bc_input, tech = v.split("|")
            nodes.add(tech)
            nodes.add(bc_input)
            nodes.add(bc_output)
            if bc_output == "Demand":
                mapping[v] = (bc_input, tech)
            elif bc_output == "Losses":
                mapping[v] = (tech, bc_output)
            else:  # Link supply
                mapping[v] = (tech, bc_input)
        elif v.startswith("Final"):
            _, bus_carrier, tech = v.split("|")
            mapping[v] = (bus_carrier, tech)
            nodes.add(tech)
            nodes.add(bus_carrier)
        else:
            raise ValueError(f"Unexpected variable '{v}'")

    return mapping, nodes


def sort_mapping(k):
    if k.startswith("Primary"):
        return 0
    elif k.startswith("Secondary"):
        return 1
    elif k.startswith("Final"):
        return 2
    else:
        raise ValueError(f"Unexpected key '{k}'")


def remove_missing_variables(m: dict) -> dict:
    clean_mapping = {}
    variables = df.index.unique("Variable")
    for k, v in m.items():
        if k in variables:
            clean_mapping[k] = v
        else:
            print(f"Skipping '{k}' because it does not exist in AT {year}.")
    return clean_mapping


def get_xmap(nodes) -> dict:
    # dict.fromkeys(sorted(nodes), "")
    return {
        "AC": 0.5,
        "Agriculture": 1,
        "Air Heat Pump": 0.5,
        "Ambient Heat": "",
        "BEV charger": 0.5,
        "Base Load": "",
        "Biogas CC": 0.5,
        "Biomass": 0.5,
        "Boiler": 0.5,
        "CHP": 0.5,
        "Distribution Grid": 0.5,
        "Electrolysis": 0.5,
        "Export": 0.0,
        "Export Domestic": 1.0,
        "Export Foreign": 1.0,
        "Fischer-Tropsch": 0.5,
        "Gas": 0.5,
        "Gas Compressing": 0.5,
        "Ground Heat Pump": 0.5,
        "H2": 0.5,
        "H2 Compressing": 0.5,
        "HH & Services": 1.0,
        "HVC from naphtha": 0.5,
        "HVC to air": 0.5,
        "Heat": 0.5,
        "Import Domestic": 0.0,
        "Import Foreign": 0.0,
        "Import Global": 0.0,
        "Industry": 1.0,
        "Industry CC": 1.0,
        "Losses": 1.0,
        "Methanol": 0.5,
        "Methanolisation": 0.5,
        "Oil": 0.5,
        "Powerplant": 0.5,
        "Resistive Heater": 0.5,
        "Run-of-River": 0.0,
        "Sabatier": 0.5,
        "Solar Rooftop": 0.0,
        "Solar Utility": 0.0,
        "Solid": 0.5,
        "Transport": 1.0,
        "Waste": 0.5,
        "Water Pits": 0.5,
        "Water Tank": 0.5,
        "Wind Onshore": 0.0,
    }


if __name__ == "__main__":
    df = read_iamc_data_frame(
        filepath="/IdeaProjects/pypsa-at/results/v2025.03/AT10_KN2040/evaluation/exported_iamc_variables.xlsx"
    )
    mapping, nodes = get_mapping(df)
    mapping_sorted = {k: mapping[k] for k in sorted(mapping, key=sort_mapping)}
    xmap = get_xmap(nodes)
    df = rename_aggregate(df, "TWh", level="Unit").div(1e6)
    year = "2050"

    at_regions = [s for s in df.index.unique("Region") if s.startswith("AT")]
    df = filter_by(df, Year=year, Region=at_regions)
    df = rename_aggregate(df, "AT", level="Region")

    # AC
    variable_mapper_ac = {
        "Primary Energy|AC|Import Domestic": "Primary Energy|AC|Total Import",
        "Primary Energy|AC|Import Foreign": "Primary Energy|AC|Total Import",
        "Primary Energy|AC|Reservoir": "Primary Energy|AC|Hydro Power",
        "Primary Energy|AC|Run-of-River": "Primary Energy|AC|Hydro Power",
        "Primary Energy|AC|Solar HSAT": "Primary Energy|AC|Solar Power",
        "Primary Energy|AC|Solar Rooftop": "Primary Energy|AC|Solar Power",
        "Primary Energy|AC|Solar Utility": "Primary Energy|AC|Solar Power",
        "Primary Energy|AC|Wind Onshore": "Primary Energy|AC|Wind Power",
        "Primary Energy|AC|Wind Offshore": "Primary Energy|AC|Wind Power",
    }
    df = rename_aggregate(df, variable_mapper_ac, level="Variable")

    mapping = {
        "Primary Energy|AC|Total Import": ("Import", "AC"),
        "Primary Energy|AC|Hydro Power": ("Hydro Power", "AC"),
        "Primary Energy|AC|Solar Power": ("Solar Power", "AC"),
        "Primary Energy|AC|Wind Power": ("Wind Power", "AC"),
    }

    # clean_mapping = remove_missing_variables(mapping_sorted)

    iamc = pyam.IamDataFrame(df)

    iamc_fig = sankey(iamc, mapping=mapping)
    iamc_fig.update_layout(height=800)
    # node = iamc_fig.data[0].node.to_plotly_json()
    # link = iamc_fig.data[0].link.to_plotly_json()
    #
    # node["x"] = [xmap.get(label, 0.2) for label in node["label"]]
    # node["y"] = [xmap.get(label, 0.4) for label in node["label"]]
    #
    # new_sankey = Sankey(
    #     node=node,
    #     link=link,
    #     arrangement="fixed",  # necessary for x/y positions
    # )
    #
    # fig = Figure(data=[new_sankey])
    #
    # fig.update_layout(height=800)

    plotly.io.show(iamc_fig)
