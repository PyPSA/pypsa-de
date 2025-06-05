"""Module for methane nodal balances."""

from pathlib import Path

from evals.plotting import simple_bus_balance


def view_balance_methane(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """
    Evaluate the methane balance.

    Returns
    -------
    :
    """
    simple_bus_balance(networks, config, result_path)
