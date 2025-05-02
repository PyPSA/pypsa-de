"""Modify the NUTS3 shapefile for custom administrative clustering."""

import sys
import logging

import geopandas as gpd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def override_nuts(nuts_code: str | tuple, override: str, level: str = "level1") -> None:
    """
    Update the NUTS codes.

    Parameters
    ----------
    nuts_code
        The NUTS codes substrings in the index used to identify
        regions that should be updated.
    override
        The value to set for the specified regions.
    level
        The level to set the override value for.

    Returns
    -------
    :
    """
    mask = nuts3_regions.index.str.startswith(nuts_code)
    nuts3_regions.loc[mask, level] = override


def assert_expected_number_of_entries(nuts_code: str, expected: int, lvl: int = 1):
    """
    Ensure that a specific number of entries are present for a NUTS code.

    Parameters
    ----------
    nuts_code
        The NUTS code to check for.
    expected
        The expected number of entries.
    lvl
        The level to check the `nuts_code` at.

    Raises
    ------
    AssertionError
        If the number of entries does not match the expected value.
    """
    regions_at_level = nuts3_regions.query(f"level{lvl}.str.startswith(@nuts_code)")
    entries = regions_at_level[f"level{lvl}"].unique()
    assert len(entries) == expected


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("modify_nuts3_shapes")

    configure_logging(snakemake)
    config = snakemake.config

    admin_levels = snakemake.params.get("admin_levels")
    nuts3_regions = gpd.read_file(snakemake.input.nuts3_shapes).set_index("index")

    if not (
        config.get("modify_nuts3_shapes", {}).get("enable")
        and config["clustering"]["mode"] == "administrative"
    ):
        logger.info("Skipping NUTS3 shapefile modification.")
        nuts3_regions.to_file(snakemake.output.nuts3_shapes)
        sys.exit(0)

    assert admin_levels.get("level") == 0
    logger.info("Applying custom administrative clustering.")

    # AT: 10
    assert admin_levels.get("AT") == 2
    override_nuts("AT333", "AT333", "level2")
    assert_expected_number_of_entries("AT", expected=10, lvl=2)
    # IT: 3
    assert admin_levels.get("IT") == 1
    override_nuts("IT", "IT0")  # mainland
    override_nuts("ITG1", "IT1")  # Sicily
    override_nuts("ITG2", "IT2")  # Sardinia
    assert_expected_number_of_entries("IT", expected=3)
    # DK: 2
    assert admin_levels.get("DK") == 1
    override_nuts("DK", "DK0")
    override_nuts(("DK01", "DK02"), "DK1")  # Sjaelland
    assert_expected_number_of_entries("DK", expected=2)
    # UK: 2
    assert admin_levels.get("GB") == 1
    override_nuts("GB", "GB0")
    override_nuts("GBN", "GB1")  # North Ireland
    assert_expected_number_of_entries("GB", expected=2)
    # FR: 2
    assert admin_levels.get("FR") == 1
    override_nuts("FR", "FR0")
    override_nuts("FRM0", "FR1")  # Corsica
    assert_expected_number_of_entries("FR", expected=2)
    # ES: 2
    assert admin_levels.get("ES") == 1
    override_nuts("ES", "ES0")
    override_nuts("ES53", "ES1")  # Balearic Islands
    assert_expected_number_of_entries("ES", expected=2)

    nuts3_regions.to_file(snakemake.output.nuts3_shapes)
