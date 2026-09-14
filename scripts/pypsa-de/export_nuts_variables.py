import logging
import os
import sys

import pandas as pd
import pypsa

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + "/../.."))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from functools import reduce

from export_ariadne_variables import get_capacities

from scripts._helpers import (
    configure_logging,
    mock_snakemake,
    set_scenario_config,
    update_config_from_wildcards,
)

pypsa.options.params.statistics.round = 10
logger = logging.getLogger(__name__)

# Defining global variables

TWh2MWh = 1e6
MWh2TWh = 1e-6
MW2GW = 1e-3
MW2TW = 1e-6
t2Mt = 1e-6
toe_to_MWh = 11.630  # GWh/ktoe OR MWh/toe


EUR20TOEUR23 = 1.1076


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_nuts_variables",
            simpl="",
            clusters=27,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Bal_v5",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    ariadne_template = pd.read_excel(snakemake.input.template, sheet_name=None)
    var2unit = ariadne_template["variable_definitions"].set_index("Variable")["Unit"]
    var2unit = var2unit.str.replace("GJ", "MWh")
    var2unit = var2unit.str.replace("PJ", "TWh")

    # Load data
    networks = [pypsa.Network(fn) for fn in snakemake.input.networks]

    modelyears = [fn[-7:-3] for fn in snakemake.input.networks]

    if "debug" == "debug":  # For debugging
        var = pd.Series()
        idx = 1
        n = networks[idx]
        region = "DE"
        kwargs = {
            "groupby": ["bus", "carrier"],
            "at_port": True,
            "nice_names": False,
        }

    nuts_regions = n.buses.query("carrier == 'AC' and name.str.startswith('DE')").index
    assert "DEA" in nuts_regions, (
        "DEA not in NUTS regions, probably the run did not use the right NUTS level."
    )
    yearly_dfs = []
    for n, year in zip(networks, modelyears):
        regional_dfs = []
        for region in nuts_regions:
            logger.info(f"Getting data for year {year} and region {region}...")
            var = get_capacities(n, region=region)
            _df = var.rename_axis("Variable").to_frame(name=year).reset_index()
            _df["Model"] = "PyPSA-DE v0.1.1"
            _df["Scenario"] = snakemake.wildcards.run
            _df["Region"] = region
            _df["Unit"] = _df["Variable"].map(var2unit).fillna("NA")
            _df = _df[["Model", "Scenario", "Region", "Variable", "Unit", year]]
            regional_dfs.append(_df)
        yearly_dfs.append(pd.concat(regional_dfs, ignore_index=True))

    df = reduce(
        lambda left, right: pd.merge(
            left, right, on=["Model", "Scenario", "Region", "Variable", "Unit"]
        ),
        yearly_dfs,
    )

    not_in_template = df.loc[df["Unit"] == "NA"].index
    logger.info(
        "Dropping variables which are not in the template: %s", not_in_template.tolist()
    )
    ariadne_df = df.drop(not_in_template)

    meta = pd.Series(
        {
            "Model": "PyPSA-DE v0.1.1",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "preliminary",
            "Internal usage within Kopernikus AG Szenarien": "yes",
            "Release for publication": "no",
        }
    )

    # For export to the Ariadne-internal DB, convert most Wh-based entries to J
    ariadne_df.loc[ariadne_df["Unit"] == "TWh/yr", modelyears] *= 3.6
    ariadne_df.loc[ariadne_df["Unit"] == "TWh/yr", "Unit"] = "PJ/yr"
    ariadne_df.loc[ariadne_df["Unit"] == "EUR2020/MWh", modelyears] /= 3.6
    ariadne_df.loc[ariadne_df["Unit"] == "EUR2020/MWh", "Unit"] = "EUR2020/GJ"

    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        ariadne_df.round(5).to_excel(writer, sheet_name="data", index=False)
        meta.to_frame().T.to_excel(writer, sheet_name="meta", index=False)
