from pathlib import Path

from evals.fileio import Exporter
from evals.plots import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.views.common import _parse_view_config_items, simple_optimal_capacity


def view_capacity_gas_storage(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate optimal storage capacities for CH4 and H2.

    Returns
    -------
    :

    Notes
    -----
    FixMe: No Hydrogen Storage with current config?
    """
    (
        bus_carrier,
        transmission_comps,
        transmission_carrier,
        storage_carrier,
    ) = _parse_view_config_items(networks, config)

    gas_stores = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        storage=True,
        bus_carrier=bus_carrier,
    )

    metric = Exporter(
        statistics=[gas_stores],
        view_config=config["view"],
    )

    metric.defaults.plotly.chart = ESMBarChart
    # prevent dropping empty years
    metric.defaults.plotly.cutoff_drop = False

    metric.export(result_path, subdir)


def view_capacity_electricity_production(
    result_path: str | Path,
    networks: dict,
    config: dict,
) -> None:
    """
    Evaluate the optimal capacity for AC technologies.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country.
    """
    simple_optimal_capacity(networks, config, result_path, kind="production")
