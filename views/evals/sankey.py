"""Module for Sankey evaluations."""

from functools import partial
from pathlib import Path

import pandas as pd

from esmtools.constants import BusCarrier, Carrier, DataModel, TradeTypes
from esmtools.fileio import prepare_industry_demand
from esmtools.metric import Metric
from esmtools.plots.sankey import OverviewSankey
from esmtools.statistic import (
    collect_myopic_statistics,
    get_location_and_carrier_and_bus_carrier,
)
from esmtools.utils import (
    Mapping,
    filter_by,
    get_mapping,
    insert_index_level,
    rename_aggregate,
)


def eval_sankey(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Produce Sankey diagrams.

    # https://plotly.com/blog/sankey-diagrams/

    Returns
    -------
    :
        Produces one HTML file per node and year.
    """

    ambient_heat = ""

    generation = build_metric_generation(networks, result_path)
    transform_in = build_metric_transformation_in(networks)
    # transform_out = build_metric_transformation_out(networks)
    # demand = build_metric_demand(networks, result_path)

    generation = format_metric_for_sankey(generation)
    transform_in = format_metric_for_sankey(transform_in)
    # transform_out = format_metric_for_sankey(transform_out)
    # demand = format_metric_for_sankey(demand)

    to_concat = [generation.squeeze(), transform_in.squeeze()]
    sankey_data = pd.concat(to_concat)

    for (year, location), df in sankey_data.groupby(
        [DataModel.YEAR, DataModel.LOCATION]
    ):
        overview_sankey = OverviewSankey(df, year, location)
        overview_sankey.plot()
        file_name = f"sankey_overview_{location}_{year}.html"
        file_path = Path(result_path, file_name)
        overview_sankey.fig.write_html(file_path, include_plotlyjs="cdn")
        break


def build_metric_generation(networks: dict, result_path: Path) -> Metric:
    """Build the left side input data for the sankey chart.

    Parameters
    ----------
    networks
    result_path

    Returns
    -------
    :
        The built metric object.
    """
    generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="1"),
        comps=["Generator", "StorageUnit"],
    ).drop("PHS", level="carrier")  # included in phs statistic

    trade_import = collect_myopic_statistics(
        networks,
        statistic="trade_energy",
        scope=(TradeTypes.FOREIGN, TradeTypes.DOMESTIC),
        direction="import",
    )

    biogas_generation = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Link",
        groupby=partial(get_location_and_carrier_and_bus_carrier, location_port="1"),
        bus_carrier="gas",
        carrier=["biogas to gas"],
    )

    _industry = prepare_industry_demand(result_path, networks)
    biomass_for_industry = filter_by(_industry, carrier=BusCarrier.SOLID_BIOMASS)
    biomass_for_industry = rename_aggregate(
        biomass_for_industry, "solid biomass", level="bus_carrier"
    )
    biomass_for_industry = rename_aggregate(
        biomass_for_industry, "solid biomass for industry", level="carrier"
    )

    hydro_supply = collect_myopic_statistics(
        networks,
        statistic="phs_split",
        carrier=[Carrier.phs_dispatched_power_inflow, Carrier.hydro_dispatched_power],
        # + Carrier.ror ?
    )

    # todo: ambient heat statistic

    # optionally configure the metric
    return Metric(
        metric_name="Generation",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[
            generation,
            trade_import,
            biogas_generation,
            biomass_for_industry,
            hydro_supply,
        ],
    )


def build_metric_transformation_in(networks: dict) -> Metric:
    """Build the left side of the transformation and storage node."""
    transformation_input = collect_myopic_statistics(
        networks,
        "energy_input",
        comps="Link",
        include_losses=True,
    )
    return Metric(
        metric_name="Transformation Input",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[transformation_input],
    )


def build_metric_transformation_out(networks: dict) -> Metric:
    """Build the right side of the transformation and storage node."""
    return Metric(
        "Transformation Output", is_unit="", to_unit="", statistics=[pd.DataFrame()]
    )


def build_metric_demand(networks: dict, result_path: Path) -> Metric:
    """Build the right side output data for the sankey chart."""
    return Metric(
        "Transformation Output", is_unit="", to_unit="", statistics=[pd.DataFrame()]
    )


def format_metric_for_sankey(metric: Metric) -> pd.DataFrame:
    """Apply aggregations and return a data frame.

    Parameters
    ----------
    metric

    Returns
    -------
    :
    """
    mapper = get_mapping(metric.cfg.mapping)
    df = rename_aggregate(metric.df, mapper=mapper)
    df = rename_aggregate(df, mapper=Mapping.bus_carrier, level=DataModel.BUS_CARRIER)

    return insert_index_level(df, df.attrs["name"], "level")


if __name__ == "__main__":
    from esmtools.fileio import read_networks

    _result_path = Path().resolve().parents[2] / "pypsa-eur-sec" / "results"
    _networks = read_networks(_result_path)
    eval_sankey(result_path=_result_path, networks=_networks)
