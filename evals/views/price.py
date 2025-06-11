"""Evaluate nodal prices per energy bus carrier."""

from pathlib import Path

from evals.fileio import Exporter
from evals.statistic import collect_myopic_statistics


def view_price_map(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """
    Export nodal prices to file using Folium.

    Parameters
    ----------
    result_path : str | Path
        The path to the results directory.
    networks : dict
        A dictionary of networks.
    config : dict
        Configuration dictionary.
    """
    statistics = []

    collect_myopic_statistics(
        networks,
    )

    exporter = Exporter(statistics=statistics, view_config=config["view"])
    exporter.export(result_path, subdir=config["view"]["subdir"])
