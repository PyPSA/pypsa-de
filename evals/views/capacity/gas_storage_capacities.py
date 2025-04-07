# -*- coding: utf-8 -*-
"""Module for gas storage evaluations."""

from pathlib import Path

import pandas as pd

from evals.constants import TITLE_SUFFIX, BusCarrier, Carrier, DataModel, Group
from evals.fileio import Exporter
from evals.plots.barchart import ESMBarChart
from evals.plots.timeseries import ESMTimeSeriesChart
from evals.statistic import collect_myopic_statistics


def view_gas_storage_capacities(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:  # numpydoc ignore=PR01
    """
    Evaluate optimal storage capacities for CH4 and H2.

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
        bus_carrier=[BusCarrier.H2, BusCarrier.CH4],
    )

    # The pypsa statistic wrongly returns MW as a unit in
    # df.attrs. We override it with MWh here. Plus, we need
    # a data frame for the function below to work.
    gas_stores = gas_stores.to_frame(f"{gas_stores.attrs['name']} (MWh)")

    def update_carrier_with_storage_type_suffix(df: pd.DataFrame) -> pd.DataFrame:
        """
        Append the storage type to the carrier name.

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

    metric = Exporter(
        statistics=[gas_stores.squeeze()],
        statistics_unit="MWh",
        view_config=config["view"],
    )

    metric.defaults.plotly.chart = ESMBarChart
    # prevent dropping empty years
    metric.defaults.plotly.cutoff_drop = False

    metric.export(result_path, subdir)
