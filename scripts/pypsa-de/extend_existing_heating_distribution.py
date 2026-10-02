# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/PyPSA/pypsa-eur>
# SPDX-License-Identifier: MIT
"""Extend existing heating shares to district heating subnodes."""

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging


def extend_heating_distribution(
    existing_heating_distribution: pd.DataFrame, subnodes: gpd.GeoDataFrame
) -> pd.DataFrame:
    """
    Extend heating distribution by subnodes mirroring the distribution of the
    corresponding mother node.

    Parameters
    ----------
    existing_heating_distribution : pd.DataFrame
        DataFrame containing the existing heating distribution.
    subnodes : gpd.GeoDataFrame
        GeoDataFrame containing information about district heating subnodes.

    Returns
    -------
    pd.DataFrame
        Extended DataFrame with heating distribution for subnodes.
    """
    # Merge the existing heating distribution with subnodes on the cluster name
    mother_nodes = (
        existing_heating_distribution.loc[subnodes.cluster.unique()]
        .unstack(-1)
        .to_frame()
    )
    cities_within_cluster = subnodes.groupby("cluster")["Stadt"].apply(list)
    mother_nodes["cities"] = mother_nodes.apply(
        lambda i: cities_within_cluster[i.name[2]], axis=1
    )
    # Explode the list of cities
    mother_nodes = mother_nodes.explode("cities")

    # Reset index to temporarily flatten it
    mother_nodes_reset = mother_nodes.reset_index()

    # Append city name to the third level of the index
    mother_nodes_reset["name"] = (
        mother_nodes_reset["name"] + " " + mother_nodes_reset["cities"]
    )

    # Set the index back
    mother_nodes = mother_nodes_reset.set_index(["heat name", "technology", "name"])

    # Drop the temporary 'cities' column
    mother_nodes.drop("cities", axis=1, inplace=True)

    # Reformat to match the existing heating distribution
    mother_nodes = mother_nodes.squeeze().unstack(-1).T

    # Combine the exploded data with the existing heating distribution
    existing_heating_distribution_extended = pd.concat(
        [existing_heating_distribution, mother_nodes]
    )
    return existing_heating_distribution_extended


if "snakemake" not in globals():
    from scripts._helpers import mock_snakemake

    snakemake = mock_snakemake("extend_existing_heating_distribution")

configure_logging(snakemake)
existing_heating_distribution = pd.read_csv(
    snakemake.input.existing_heating_distribution,
    header=[0, 1],
    index_col=0,
)
subnodes = gpd.read_file(snakemake.input.subnodes)
extend_heating_distribution(existing_heating_distribution, subnodes).to_csv(
    snakemake.output.existing_heating_distribution_extended
)
subnodes.sort_values(by="Wärmeeinspeisung in GWh/a", ascending=False).head(
    snakemake.params.nlargest
).to_file(snakemake.output.district_heating_subnodes_selected, driver="GeoJSON")
