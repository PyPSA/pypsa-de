from pathlib import Path

from evals.views.common import simple_timeseries


def view_timeseries_electricity(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """Evaluate the electricity balance time series."""
    simple_timeseries(networks, config, result_path)


def view_timeseries_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """Evaluate the Hydrogen balance time series."""
    simple_timeseries(networks, config, result_path)


def view_timeseries_methane(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """Evaluate the Methane balance time series."""
    simple_timeseries(networks, config, result_path)


def view_timeseries_carbon(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """Evaluate the Carbon balance time series."""
    simple_timeseries(networks, config, result_path)
