# -*- coding: utf-8 -*-
import logging

import pandas as pd
import pypsa

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "check_sector_ratios",
            simpl="",
            clusters=22,
            opts="",
            ll="v1.2",
            sector_opts="None",
            planning_horizons="2030",
            run="KN2045_Bal_v4",
        )
    logger.info("Check sector ratios")

    n = pypsa.Network(snakemake.input.network)

    # check the heating sector

    # check the industry sector
