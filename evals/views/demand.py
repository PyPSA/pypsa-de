from pathlib import Path

from evals.constants import BusCarrier, DataModel
from evals.fileio import Exporter
from evals.plots import ESMBarChart
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    calculate_input_share,
    drop_from_multtindex_by_regex,
    filter_for_carrier_connected_to,
)


def view_demand_heat(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:
    """
    Evaluate the energy required for heat production and generation.

    Results are grouped by bus_carrier and not by carrier
    as usual to show the input energy carrier mix.

    Returns
    -------
    :

    Notes
    -----
    See eval docstring for parameter description.
    """
    energy_for_heat = (
        collect_myopic_statistics(networks, comps="Link", statistic="energy_balance")
        # todo: is dropping CO2 really justified? Discussions needed, or disclaimer in graph.
        # .drop(["co2", "co2 stored"], level=DataModel.BUS_CARRIER)
        .pipe(drop_from_multtindex_by_regex, "water tanks")
        .pipe(filter_for_carrier_connected_to, BusCarrier.heat_buses(), kind="supply")
        .pipe(calculate_input_share, BusCarrier.HEAT_RURAL)
    )

    generator_supply = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=BusCarrier.heat_buses(),
    )

    exporter = Exporter(
        statistics=[energy_for_heat.mul(-1), generator_supply],
        view_config=config["view"],
    )

    # view specific static settings:
    exporter.defaults.plotly.chart = ESMBarChart
    exporter.defaults.excel.pivot_index = [DataModel.LOCATION, DataModel.BUS_CARRIER]
    exporter.defaults.plotly.plot_category = DataModel.BUS_CARRIER
    exporter.defaults.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]

    exporter.export(result_path, subdir=subdir)
