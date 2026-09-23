# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Precompute real ``base_year`` EUR daily fuel prices and store them in the CSV.

The daily fuel price CSV keeps the raw nominal prices (one ``<fuel>`` column)
for provenance and a precomputed deflated ``<fuel>_real<base_year>`` column that
``modify_prenetwork.add_daily_fuel_prices`` consumes at model runtime. Doing the
deflation here, once, removes the fragile ``pydeflate`` network dependency from
every model run.

Deflation uses the German GDP deflator via ``pydeflate``, the same method as
``build_monthly_prices``. It is a one-off maintenance script, run manually when
the nominal series changes::

    pixi run python scripts/pypsa-de/deflate_daily_fuel_prices.py
"""

import logging

import pandas as pd

logger = logging.getLogger(__name__)

FN = "data/pypsa-de/daily_fuel_prices.csv"
BASE_YEAR = 2020
FUELS = ["gas"]


def _patch_imf_session() -> None:
    """
    Route ``imf_reader`` through its uncached session.

    Works around a ``requests_cache`` x ``cattrs`` incompatibility in the pinned
    environment that crashes when caching the IMF WEO download. pydeflate still
    caches the parsed dataset as parquet, so the download runs only once.
    """
    import imf_reader.cache.config as cfg
    import imf_reader.utils as utils

    utils.get_session = cfg.get_uncached_session


def deflate(prices: pd.DataFrame, base_year: int, iso: str = "DEU") -> pd.DataFrame:
    """
    Deflate nominal EUR prices to real ``base_year`` EUR (German GDP deflator).

    Parameters
    ----------
    prices : pd.DataFrame
        Nominal prices with a ``DatetimeIndex`` and one column per fuel.
    base_year : int
        Reference year to deflate to.
    iso : str
        ISO3 country code selecting the deflator series.

    Returns
    -------
    pd.DataFrame
        Real prices in ``base_year`` EUR, same shape as ``prices``.
    """
    from pydeflate import imf_gdp_deflate, set_pydeflate_path

    set_pydeflate_path("../data/pydeflate/")
    df = prices.assign(year=prices.index.year, iso_code=iso)
    real = pd.DataFrame(index=prices.index)
    for col in prices.columns:
        real[col] = imf_gdp_deflate(
            df,
            value_column=col,
            source_currency="EUR",
            target_currency="EUR",
            base_year=base_year,
        )["value"].values
    return real


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    _patch_imf_session()

    df = pd.read_csv(FN, index_col=0, parse_dates=True)
    real = deflate(df[FUELS], BASE_YEAR)
    for fuel in FUELS:
        df[f"{fuel}_real{BASE_YEAR}"] = real[fuel]
        annual = real[fuel].groupby(real.index.year).mean().round(2)
        logger.info(f"{fuel} real{BASE_YEAR} annual mean:\n{annual.to_string()}")

    df.to_csv(FN)
    logger.info(f"Wrote {FN} with columns {list(df.columns)}.")
