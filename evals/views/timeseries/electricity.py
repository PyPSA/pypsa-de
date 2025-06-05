"""Export time series views for electricity."""

from pathlib import Path

from evals.views.timeseries.common import simple_timeseries


def view_timeseries_electricity(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """Evaluate the Methane balance time series."""
    simple_timeseries(networks, config, result_path)
