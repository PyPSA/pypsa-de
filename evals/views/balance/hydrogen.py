"""Module for hydrogen nodal balances."""

from pathlib import Path

from evals.views.common import simple_bus_balance


def view_balance_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """
    Evaluate the Hydrogen balance.

    Returns
    -------
    :

    Notes
    -----
    See eval module docstring for parameter description.
    """
    simple_bus_balance(networks, config, result_path)
