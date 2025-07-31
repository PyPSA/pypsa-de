# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for Sankey diagram."""

import pandas as pd
import plotly
import pyam
from plotly.graph_objs import Figure, Sankey

from evals.utils import filter_by, rename_aggregate

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)


def read_iamc_data_frame():
    xls = pd.read_excel(
        "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx",
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    return xls.stack()


def main():
    year = "2050"
    region = "GB0"

    df = read_iamc_data_frame()
    df = rename_aggregate(df, "TWh", level="Unit").div(1e6)
    df = filter_by(df, Year=year, Region=region)

    raw_mapping = get_mapping(df)
    variables = df.index.unique("Variable")
    for k, v in raw_mapping.items():
        if k in variables:
            clean_mapping[k] = v
        else:
            print(f"Skipping '{k}' because it does not exist in {region} {year}.")

    iamc = pyam.IamDataFrame(df)
    fig = iamc.plot.sankey(mapping=clean_mapping)

    fig.update_layout(height=600)

    plotly.io.show(fig)


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
    FILEPATH = "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx"
    df = read_iamc_data_frame()
    mapping, nodes = get_mapping(df)
    mapping_sorted = {k: mapping[k] for k in sorted(mapping, key=sort_mapping)}
    xmap = get_xmap(nodes)
    df = rename_aggregate(df, "TWh", level="Unit").div(1e6)
    year = "2050"
    region = "GB0"
    df = filter_by(df, Year=year, Region=region)

    clean_mapping = {}
    variables = df.index.unique("Variable")
    for k, v in mapping_sorted.items():
        if k in variables:
            clean_mapping[k] = v
        else:
            print(f"Skipping '{k}' because it does not exist in {region} {year}.")

    iamc = pyam.IamDataFrame(df)

    iamc_fig = iamc.plot.sankey(mapping=clean_mapping)
    node = iamc_fig.data[0].node.to_plotly_json()
    link = iamc_fig.data[0].link.to_plotly_json()

    node["x"] = [xmap.get(label, 0.2) for label in node["label"]]
    node["y"] = [xmap.get(label, 0.4) for label in node["label"]]
    # node["y"] = [ymap[label] for label in node["label"]]
    # node["color"] = [colormap[label] for label in node["label"]]

    new_sankey = Sankey(
        node=node,
        link=link,
        arrangement="fixed",  # necessary for x/y positions
    )

    fig = Figure(data=[new_sankey])

    fig.update_layout(height=600)

    plotly.io.show(fig)
