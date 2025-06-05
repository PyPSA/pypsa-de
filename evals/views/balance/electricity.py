"""Module for electricity evaluations."""

from pathlib import Path

from evals.plotting import simple_bus_balance


def view_balance_electricity(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """
    Evaluate the electricity production & demand by country and year.

    Returns
    -------
    :

    Notes
    -----
    Balances do nat add up to zero, because of transmission losses and
    storage cycling (probably).
    """
    simple_bus_balance(networks, config, result_path)
