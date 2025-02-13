# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Load and prepares the data of the eGo^N DemandRegio project <https://opendata.ffe.de/project/demandregio/>
about district heating in Germany on NUTS3 level.

Inputs:
    - resources/nuts3.geojson: Path to the GeoJSON file containing NUTS3 regions data.
    - resources/mapping_technologies.json: Path to the JSON file containing the mapping of technologies.
    - resources/demandregio_spatial.json: Path to the JSON file containing spatially resolved heating structure data from the DemandRegio project.
    - resources/mapping_38_to_4.json: Path to the JSON file containing the mapping of regions.

Outputs:
    - resources/heating_technologies_nuts3.geojson: Path to the GeoJSON file where the processed heating technologies data will be saved.
"""

import logging

logger = logging.getLogger(__name__)
import json
import re

import geopandas as gpd
import pandas as pd

if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_egon_data",
            run="KN2045_Bal_v4",
        )

logger.info("Retrieving and cleaning egon data")

nuts3 = gpd.read_file(snakemake.input.nuts3)[
    ["index", "pop", "geometry"]
]  # Keep only necessary columns

#  Parse dictionary from string
description = pd.read_json(snakemake.input.mapping_technologies).loc["resources"][
    "oep_metadata"
][0]["schema"]["fields"][5]["description"]
# Extract the part after the colon
key_value_part = description.split(":", 1)[1].strip()

# Split into individual key-value pairs
pairs = re.split(r";\s*", key_value_part)

# Initialize dictionary
internal_id = {}

# Process each pair
for pair in pairs:
    key, value = pair.split(":", 1)
    internal_id[int(key.strip())] = value.strip()

with open(snakemake.input.demandregio_spatial) as datafile:
    data = json.load(datafile)["data"]
df = pd.DataFrame(data)

id_region = pd.read_json(snakemake.input.mapping_38_to_4)

df["internal_id"] = df["internal_id"].apply(lambda x: x[0])

df["nuts3"] = df.id_region.map(id_region.set_index(id_region.id_region_from).kuerzel_to)

heat_tech_per_region = df.groupby([df.nuts3, df.internal_id]).sum().value.unstack()
heat_tech_per_region.rename(columns=internal_id, inplace=True)

egon_df = heat_tech_per_region.merge(nuts3, left_on="nuts3", right_on="index")
egon_gdf = gpd.GeoDataFrame(egon_df)  # Convert merged DataFrame to GeoDataFrame
egon_gdf = egon_gdf.to_crs("EPSG:4326")

egon_gdf.to_file(snakemake.output.heating_technologies_nuts3)
