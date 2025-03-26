"""Module for primary energy evaluations."""

from pathlib import Path

from esmtools.constants import TITLE_SUFFIX, BusCarrier, DataModel
from esmtools.metric import Metric
from esmtools.plots.barchart import ESMBarChart
from esmtools.statistic import collect_myopic_statistics
from esmtools.utils import (
    drop_from_multtindex_by_regex,
    make_evaluation_result_directories,
)


def eval_heat_primary_energy(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate the energy required for heat production and generation.

    Results are grouped by bus_carrier and not by carrier
    as usual to show the input energy carrier mix.

    TODO: add missing tables in Excel

    Writes 2 Excel files and 1 BarChart per country.

    Notes
    -----
    See esmtools.eval docstring for parameter description.
    """
    heat_bus_carrier = [
        BusCarrier.HEAT_URBAN_CENTRAL,
        BusCarrier.HEAT_RURAL_SERVICES,
        BusCarrier.HEAT_RURAL_RESIDENTIAL,
        # fixme: why are not all heat buses considered?
        # BusCarrier.HEAT_URBAN_RESIDENTIAL,
        # BusCarrier.HEAT_URBAN_SERVICES
    ]

    link_input = collect_myopic_statistics(
        networks,
        comps="Link",
        statistic="energy_input",
        bus_carrier=heat_bus_carrier,
    ).drop("", level=DataModel.CARRIER, errors="ignore")
    link_input = drop_from_multtindex_by_regex(link_input, "DAC|water tank")

    generator_input = collect_myopic_statistics(
        networks,
        statistic="energy_input",
        comps="Generator",
        bus_carrier=heat_bus_carrier,
    )
    generator_input = drop_from_multtindex_by_regex(generator_input, "DAC|water tank")

    metric = Metric(
        "Primary Energy Demand for Heat Production",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[link_input, generator_input],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX
    metric.cfg.mapping = "district_heat"
    metric.cfg.excel.chart_title = title
    metric.cfg.excel.pivot_index = [DataModel.LOCATION, DataModel.BUS_CARRIER]
    metric.cfg.plotly.title = title
    metric.cfg.plotly.chart = ESMBarChart
    metric.cfg.plotly.file_name_template = "heat_mix_{location}"
    metric.cfg.plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.CARRIER,
    ]

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
