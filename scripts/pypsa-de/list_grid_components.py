# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Authoring helper for :mod:`stochastic_grid_scenarios` (see
notebooks/stochastic_grid_uncertainty.ipynb): dumps the current AC line and
DC link names/capacities of a network to CSV, giving a ready-to-edit
template for a new grid-topology scenario at that network's `clusters`
resolution. Not part of the main pipeline DAG - run standalone, e.g.:

    python scripts/pypsa-de/list_grid_components.py \
        resources/.../base_s_27_..._2030_final.nc \
        data/pypsa-de/grid_scenarios/lines_template_de+inter_27cl.csv \
        data/pypsa-de/grid_scenarios/links_template_de+inter_27cl.csv
"""

import logging
import sys

import pypsa

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def write_templates(network_path, lines_path, links_path):
    n = pypsa.Network(network_path)

    n.lines[["bus0", "bus1", "s_nom"]].to_csv(lines_path)
    logger.info(f"Wrote {len(n.lines)} lines to {lines_path}")

    dc = n.links.query("carrier == 'DC'")
    dc[["bus0", "bus1", "p_nom"]].to_csv(links_path)
    logger.info(f"Wrote {len(dc)} DC links to {links_path}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        if len(sys.argv) == 4:
            network_path, lines_path, links_path = sys.argv[1:]
            write_templates(network_path, lines_path, links_path)
            sys.exit(0)

        snakemake = mock_snakemake(
            "list_grid_components",
            clusters=27,
            opts="",
            sector_opts="none",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    write_templates(
        snakemake.input.network,
        snakemake.output.lines_template,
        snakemake.output.links_template,
    )
