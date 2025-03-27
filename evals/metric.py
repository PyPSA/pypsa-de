# -*- coding: utf-8 -*-
"""Extends the StatisticsAccessor with additional metrics."""

import logging
from functools import cached_property
from pathlib import Path

import pandas as pd

from evals.configs import ViewDefaults
from evals.constants import NOW, TITLE_SUFFIX, DataModel, Group
from evals.fileio import (
    export_excel_countries,
    export_excel_regions_at,
    export_vamos_jsons,
)
from evals.utils import (
    add_dummy_rows,
    aggregate_locations,
    rename_aggregate,
    scale,
    verify_metric_format,
)

logger = logging.getLogger(__name__)


def combine_statistics(
    statistics: list,
    metric_name: str,
    is_unit: str,
    to_unit: str,
    keep_regions: tuple = ("AT",),
    region_nice_names: bool = True,
) -> pd.DataFrame:
    """Build the metric data frame from statistics.

    Parameters
    ----------
    statistics
        The statistics to combine.
    metric_name
        The metric name used in plot titles and column labels.
    is_unit
        The common unit of input statistics.
    to_unit
        The desired unit of the output metric.
    keep_regions
        A collection of country codes for which original input
        cluster codes will be included in the metric locations.
    region_nice_names
        Whether to replace location country codes with country/region
        names.

    Returns
    -------
    :
        The formatted metric in the desired unit and locations.
    """
    df = pd.concat(statistics)

    if was_series := isinstance(df, pd.Series):
        df = df.to_frame(f"{metric_name} ({is_unit})")

    df = aggregate_locations(df, keep_regions, region_nice_names)

    df.attrs["name"] = metric_name
    df.attrs["unit"] = to_unit

    df.columns.name = DataModel.METRIC if was_series else DataModel.SNAPSHOTS
    if df.columns.name == DataModel.SNAPSHOTS:
        df.columns = pd.to_datetime(df.columns, errors="raise")

    if to_unit and (is_unit != to_unit):
        df = scale(df, to_unit=to_unit)

    df = _split_trade_saldo_to_netted_import_export(df)

    verify_metric_format(df)

    return df


class Metric:
    """A class to build metrics from statistics.

    The metric data frame consists of multiple joined statistics,
    aggregated to countries and scaled to a specified unit. The
    data frame format is verified and expected by export functions.

    Parameters
    ----------
    statistics
        A list of Series for time aggregated statistics or list of
        data frames for statistics with snapshots as columns.
    statistics_unit
        The input statistics unit.
    keep_regions
        A tuple of location prefixes that are used to match
        locations to keep during aggregation.
    region_nice_names
        Whether, or not to rename country codes after aggregation
        to show the full country name.
    """

    def __init__(  # noqa: PLR0913, PLR0917
        self,
        statistics: list,
        statistics_unit: str,
        view_config: dict,
        keep_regions: tuple = ("AT",),
        region_nice_names: bool = True,
    ) -> None:
        self.statistics = statistics
        self.is_unit = statistics_unit
        self.metric_name = view_config["name"]
        self.to_unit = view_config["unit"]
        self.keep_regions = keep_regions
        self.region_nice_names = region_nice_names
        self.view_config = view_config
        self.defaults = ViewDefaults()

        # update defaults from config for this view
        self.defaults.excel.title = view_config["name"] + TITLE_SUFFIX
        self.defaults.plotly.title = view_config["name"] + TITLE_SUFFIX
        self.defaults.plotly.file_name_template = view_config["file_name"]
        self.defaults.plotly.cutoff = view_config["cutoff"]
        self.defaults.plotly.category_orders = view_config["legend_order"]

    @cached_property
    def df(self) -> pd.DataFrame:
        """Build the metric and store it as a cached property.

        (This is useful, because users do not need to remember
        building the metric data frame. It will be built once if needed)

        Returns
        -------
        :
            The cached metric data frame.
        """
        return combine_statistics(
            self.statistics,
            self.metric_name,
            self.is_unit,
            self.to_unit,
            self.keep_regions,
            self.region_nice_names,
        )

    def export_plotly(self, output_path: Path) -> None:
        """Create the plotly figure and export it as HTML and JSON.

        Parameters
        ----------
        output_path
            The path to the HTML folder with all the html files are
            stored.
        """
        cfg = self.defaults.plotly
        df = rename_aggregate(
            self.df, level=cfg.plot_category, mapper=self.view_config["categories"]
        )

        df_plot = df.pivot_table(
            index=cfg.pivot_index, columns=cfg.pivot_columns, aggfunc="sum"
        )

        df_plot = add_dummy_rows(df_plot, self.keep_regions)
        df_plot = df_plot.drop(cfg.drop_years, level=DataModel.YEAR, errors="ignore")

        json_file_paths = []
        for idx, data in df_plot.groupby(cfg.plotby):
            chart = cfg.chart(data, cfg)
            chart.plot()
            chart.to_html(output_path, cfg.plotby, idx)
            json_file_path = chart.to_json(output_path, cfg.plotby, idx)
            json_file_paths.append(json_file_path)

        if cfg.export_vamos_jsons:
            export_vamos_jsons(json_file_paths, cfg.file_name_template)

    def export_excel(self, output_path: Path) -> None:
        """Export metrics to Excel files for countries and regions.

        Parameters
        ----------
        output_path
            The path where the Excel files will be saved.
        """
        file_name_stem = self.view_config["file_name"].split("_{")[0]
        file_path = output_path / "XLSX" / f"{file_name_stem}_{NOW}.xlsx"
        with pd.ExcelWriter(file_path, engine="openpyxl") as writer:
            export_excel_countries(
                self.df, writer, self.defaults.excel, self.view_config
            )

        if self.df.columns.name == DataModel.SNAPSHOTS:
            return  # skips region sheets for time series

        file_path_at = output_path / f"{file_name_stem}_AT_{NOW}.xlsx"
        with pd.ExcelWriter(file_path_at, engine="openpyxl") as writer:
            export_excel_regions_at(
                self.df, writer, self.defaults.excel, self.view_config
            )

    def export_csv(self, output_path: Path) -> None:
        """Encode the metric da frame to a CSV file.

        Parameters
        ----------
        output_path
            The path to the CSV folder with all the csv files are
            stored.

        Returns
        -------
        :
            Writes the metric to a CSV file.
        """
        file_name = self.defaults.plotly.file_name_template.split("_{", maxsplit=1)[0]
        file_path = output_path / "CSV" / f"{file_name}_{NOW}.csv"
        self.df.to_csv(file_path, encoding="utf-8")

    def export(self, output_path: Path, export_config: dict) -> None:
        """Export the metric to formats specified in the config.

        Parameters
        ----------
        output_path
            The path to the CSV folder with all the csv files are
            stored.
        export_config
            The export configuration from the TOML file for this view.

        Returns
        -------
        :
        """
        if export_config["plotly"]:
            self.export_plotly(output_path)
        if export_config["excel"]:
            self.export_excel(output_path)
        if export_config["csv"]:
            self.export_csv(output_path)

    def consistency_checks(self, config_checks: dict) -> None:
        """Assert all categories are assigned to a group.

        The method typically is called after exporting the metric.
        Unmapped categories do not cause evaluations to fail, but
        the evaluation function should return in error state to obviate
        missing entries in the mapping.

        Parameter
        ---------
        config_checks
            A dictionary with flags for every test.

        Returns
        -------
        :

        Raises
        ------
        AssertionError
            If not all technologies or bus_carrier are
            assigned to a group.
        """
        # todo: refactor mapping to categories in a new multiindex level used in plots package
        category = self.defaults.plotly.plot_category
        categories = self.view_config["categories"]

        if config_checks["all_categories_mapped"]:
            assert self.df.index.unique(category).isin(categories.keys()).all(), (
                f"Incomplete categories detected. There are technologies in the metric "
                f"data frame, that are not assigned to a group (nice name)."
                f"\nMissing items: "
                f"{self.df.index.unique(category).difference(categories.keys())}"
            )

        if config_checks["no_superfluous_categories"]:
            superfluous_categories = self.df.index.unique(category).difference(
                categories.keys()
            )
            assert (
                len(superfluous_categories) == 0
            ), f"Superfluous categories found: {superfluous_categories}"

        if config_checks["legend_entry_order"]:
            a = set(self.view_config["legend_order"])
            b = set(categories.values())
            additional = a.difference(b)
            assert (
                not additional
            ), f"Superfluous categories defined in legend order: {additional}"
            missing = b.difference(a)
            assert (
                not missing
            ), f"Some categories are not defined in legend order: {missing}"


def _split_trade_saldo_to_netted_import_export(df: pd.DataFrame) -> pd.DataFrame:
    """Split the trade saldo carrier into import and export.

    The splitting needs to happen after the location aggregation.
    Otherwise, resulting netted import/export values are incorrect
    for countries with multiple regions, if the regions become
    aggregated, e.g. Germany.

    Parameters
    ----------
    df
        The input data frame with the foreign saldo carrier.

    Returns
    -------
    :
        The output data frame with positive trade values
        as import and negative values as export.
    """
    saldo = df.query("carrier.str.contains('saldo')")

    if saldo.empty:
        return df

    net_import = rename_aggregate(saldo.mul(saldo.gt(0)), Group.import_net)
    net_export = rename_aggregate(saldo.mul(saldo.le(0)), Group.export_net)

    saldo_carrier = saldo.index.unique("carrier")
    df_without_saldo = df.drop(saldo_carrier, level=DataModel.CARRIER)

    return pd.concat([df_without_saldo, net_import, net_export]).sort_index()
