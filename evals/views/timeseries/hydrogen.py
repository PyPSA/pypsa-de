"""Export time series views for hydrogen."""

from pathlib import Path

from evals.views.common import simple_timeseries


def view_timeseries_hydrogen(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """Evaluate the Hydrogen balance time series."""
    simple_timeseries(networks, config, result_path)
