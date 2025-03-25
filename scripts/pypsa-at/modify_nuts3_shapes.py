# -*- coding: utf-8 -*-
"""Modify the NUTS3 shapefile to lift NUTS regions."""

import logging

import geopandas as gpd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def override_nuts(nuts_code: str | tuple, override: str, level: str = "level1") -> None:
    """Update the NUTS codes.

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


def assert_expected_number_of_entries(
    nuts_code: str, expected: int, level: str = "level1"
):
    """Ensure that a specific number of entries are present for a NUTS code.

    Parameters
    ----------
    nuts_code
        The NUTS code to check for.
    expected
        The expected number of entries.
    level
        The level to check the `nuts_code` at.

    Raises
    ------
    AssertionError
        If the number of entries does not match the expected value.
    """
    entries = nuts3_regions.query("@level.str.startswith(@nuts_code)")[level].unique()
    assert len(entries) == expected


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("modify_nuts3_shapes")

    configure_logging(snakemake)

    logger.info("Applying custom administrative clustering.")
    assert admin_levels.get("level") == 0

    config = snakemake.config

    nuts3_shapes = snakemake.input.nuts3_shapes
    admin_levels = snakemake.params.get("admin_levels")

    nuts3_regions = gpd.read_file(nuts3_shapes).set_index("index")

    # AT: 10
    assert admin_levels.get("AT") == 2
    override_nuts("AT333", "AT333", "level2")
    assert_expected_number_of_entries("AT", expected=10, level="level2")
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
    # UK: 2a
    assert admin_levels.get("GB") == 1
    override_nuts("GB", "UK0")
    override_nuts("GBN", "UK1")  # North Ireland
    assert_expected_number_of_entries("UK", expected=2)
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
