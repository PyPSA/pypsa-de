# -*- coding: utf-8 -*-
"""Module for gas storage evaluations."""

from pathlib import Path

import pandas as pd
from configs import PlotConfig
from constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel, Group
from metric import Metric
from plots.barchart import ESMBarChart
from plots.timeseries import ESMTimeSeriesChart
from statistic import collect_myopic_statistics
from utils import make_evaluation_result_directories


def eval_gas_storage_capacities(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate optimal storage capacities for CH4 and H2.

    Returns
    -------
    :
        Writes 2 Excel files and 1 BarChart per country. The BarChart
        will only contain H2 store capacities, because gas capacities
        are not yet optimized.

    Notes
    -----
    See eval docstring for parameter description.
    """
    gas_stores = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        storage=True,
        # include CH4 once CH4 stores are subject to optimization
        bus_carrier=[BusCarrier.H2],
    )

    # The pypsa statistic wrongly returns MW as a unit in
    # df.attrs. We override it with MWh here. Plus, we need
    # a data frame for the function below to work.
    gas_stores = gas_stores.to_frame(f"{gas_stores.attrs['name']} (MWh)")

    def update_carrier_with_storage_type_suffix(df: pd.DataFrame) -> pd.DataFrame:
        """Append the storage type to the carrier name.

        Every H2 bus has one store attached. The capital costs per
        energy unit of tube stores are much higher than the capital
        cost for underground caverns. The 'e_nom_opt' value serves
        as a decision criteria.

        Parameters
        ----------
        df
            The input data frame with the carrier index level to update.

        Returns
        -------
        :
            The data frame with updated carrier index level values.
        """
        _year = df.index.unique("year")[0]
        stores = networks[_year].df("Store")

        index_list = []
        for (year, location, carrier, bus_carrier), _ in df.iterrows():
            if carrier != "H2":
                index_list.append((year, location, carrier, bus_carrier))
                continue
            e_nom_opt = stores.loc[f"{location} H2 Store", "e_nom_opt"]
            storage_type = "cavern" if e_nom_opt > 1 else "tube"
            index_list.append(
                (year, location, f"{carrier} {storage_type}", bus_carrier)
            )

        df.index = pd.MultiIndex.from_tuples(index_list, names=DataModel.YEAR_IDX_NAMES)

        return df

    gas_stores = gas_stores.groupby("year", group_keys=False).apply(
        update_carrier_with_storage_type_suffix
    )

    metric = Metric(
        metric_name="Hydrogen Store Volume",  # change once CH4 is added
        is_unit="MWh",
        to_unit="TWh",
        statistics=[gas_stores.squeeze()],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX

    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "StoreVol_{location}"

    metric.defaults.plotly.cutoff = 0.005  # TWh
    metric.defaults.plotly.cutoff_drop = False

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)


def eval_phs_hydro_operation(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate storage operation over time for PHS and hydro stores.

    The evaluation shows accumulated time series for natural inflow,
    pumped inflow, power generation outflow and natural losses (spill)
    along the storage fill level (state of charge) and the available
    volume.

    Returns
    -------
    :
        For both carrier (PHS and hydro), 2 Excel files and 1 TimeSeries
        Chart are written on a per country and per year basis.

    Notes
    -----
    See eval docstring for parameter description.
    """
    output_path = make_evaluation_result_directories(result_path, subdir)

    evaluation_types = [
        (
            Carrier.phs,
            "hydro_phs_storage_operation_{year}_{location}",
            "Pumped Hydro Storage Operation",
        ),
        (
            Carrier.hydro,
            "hydro_res_storage_operation_{year}_{location}",
            "Reservoir Hydro Storage Operation",
        ),
    ]

    phs_hydro_operation = collect_myopic_statistics(
        networks, statistic="phs_hydro_operation"
    )

    cfg_plotly = PlotConfig()
    cfg_plotly.chart = ESMTimeSeriesChart
    cfg_plotly.plotby = [DataModel.YEAR, DataModel.LOCATION]
    cfg_plotly.plot_category = DataModel.BUS_CARRIER
    cfg_plotly.pivot_index = [
        DataModel.YEAR,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    cfg_plotly.unit = "TWh"
    cfg_plotly.stacked = False
    cfg_plotly.yaxis_color = "black"
    cfg_plotly.line_shape = "spline"
    cfg_plotly.line_dash = {"Max State of Charge": "dash"}
    cfg_plotly.legend_header = "Types of Storage"
    cfg_plotly.footnotes = (
        "Storages are located according to their market connection.",
        "",
    )
    cfg_plotly.category_orders = (
        Group.inflow_cum,
        Group.pumping_cum,
        Group.turbine_cum,
        Group.soc_max,
        Group.soc,
    )

    # evaluation for PHS and hydro are almost identical. Combining them
    # in one evaluation function seems justified.
    for carrier, file_name_template, metric_name in evaluation_types:
        statistic = phs_hydro_operation.query(f"{DataModel.CARRIER} == '{carrier}'")
        metric = Metric(
            metric_name=metric_name, is_unit="", to_unit=1e-6, statistics=[statistic]
        )
        metric.defaults.plotly = cfg_plotly  # use defaults

        ts_names = metric.df.index.unique(DataModel.BUS_CARRIER)
        metric.defaults.excel.chart_title = (
            metric.df.attrs["name"] + " in" + TITLE_SUFFIX
        )
        metric.defaults.excel.axis_labels = [metric.df.attrs["name"], "TWh"]
        metric.defaults.plotly.title = (
            metric.df.attrs["name"] + " in {location} in {year} in {unit}"
        )
        metric.defaults.plotly.file_name_template = file_name_template
        metric.defaults.plotly.line_width = {
            c: 5 if c in (Group.soc, Group.soc_max) else 3 for c in ts_names
        }
        metric.defaults.plotly.fill = {
            c: "tozeroy" if c == Group.soc else None for c in ts_names
        }
        # metric.export_excel(output_path)
        # metric.export_csv(output_path)
        metric.export_plotly(output_path)


# def eval_frequency_analysis(
#     result_path: str | Path,
#     networks: dict,
#     subdir: str | Path = "esm_run/evaluation",
# ) -> None:  # numpydoc ignore=PR01
#     """"""


def eval_electricity_storage(
    result_path: str | Path,
    networks: dict,
    subdir: str | Path = "esm_run/evaluation",
) -> None:  # numpydoc ignore=PR01
    """Evaluate electricity Stores."""
    ac_storage = collect_myopic_statistics(
        networks,
        statistic="optimal_capacity",
        storage=True,
        bus_carrier=[BusCarrier.AC, BusCarrier.BATTERY],
    )

    metric = Metric(
        metric_name="Power Storage Volumes",
        is_unit="MWh",
        to_unit="TWh",
        statistics=[ac_storage],
    )

    title = metric.df.attrs["name"] + TITLE_SUFFIX

    metric.defaults.excel.chart_title = title
    metric.defaults.plotly.title = title
    metric.defaults.plotly.chart = ESMBarChart
    metric.defaults.plotly.file_name_template = "elec_storage_volumes_{location}"

    metric.defaults.plotly.cutoff = 1e-6  # 1 MWh  # 0.005
    metric.defaults.plotly.cutoff_drop = False

    metric.defaults.plotly.category_orders = (
        Group.reservoir,
        Group.battery_storage,
        Group.phs,
    )

    metric.defaults.plotly.footnotes = (
        " Further potential power storage volumes resulting from vehicle-to-grid "
        "options are not included here. <br> Storage volumes smaller than 1 TWh "
        "are not displayed. <br> Storages are located according to their market "
        "connection.",
        "",
    )

    output_path = make_evaluation_result_directories(result_path, subdir)
    metric.export_excel(output_path)
    metric.export_csv(output_path)
    metric.export_plotly(output_path)
