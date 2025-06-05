"""Module for carbon dioxide nodal balances."""

from pathlib import Path

from evals.plotting import simple_bus_balance


def view_balance_carbon(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the carbon balance.

    Returns
    -------
    :
    """
    simple_bus_balance(networks, config, result_path)
