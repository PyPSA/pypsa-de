"""Update population layouts for urban, rural, or total."""

import logging

import xarray as xr

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def main():
    logger.info("Modify Population Layouts.")
    for fp_input, fp_output in zip(snakemake.input, snakemake.output):
        # dummy placeholder until update data is available
        xr.open_dataset(fp_input).to_netcdf(fp_output)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "modify_population_layouts",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

    if snakemake.config["mods"].get("modify_population_layouts"):
        main()
