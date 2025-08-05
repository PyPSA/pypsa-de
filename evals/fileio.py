# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Input - Output related functions."""

import logging
import re
from functools import cached_property
from importlib import resources
from pathlib import Path
from typing import Callable

import pandas as pd
import pypsa
import tomllib
from pydantic.v1.utils import deep_update

from evals.configs import ViewDefaults
from evals.constants import (
    ALIAS_COUNTRY,
    ALIAS_REGION,
    COLOUR_SCHEME,
    NOW,
    TITLE_SUFFIX,
    DataModel,
    Regex,
)
from evals.excel import export_excel_countries, export_excel_regions_at
from evals.utils import (
    combine_statistics,
    insert_index_level,
    rename_aggregate,
)
from scripts._helpers import get_rdir, path_provider


def read_networks(
    result_path: str | Path | list, sub_directory: str = "networks"
) -> dict:
    """
    Read network results from NetCDF (.nc) files.

    The function returns a dictionary of data frames. The planning
    horizon (year) is used as dictionary key and added to the network
    as an attribute to associate the year with it. Network snapshots
    are equal for all networks, although the year changes. This is
    required to align timestamp columns in a data frame. Snapshots
    will become fixed late in the evaluation process (just before
    export to file).

    In addition, the function patches the statistics accessor attached
    to loaded networks and adds the configuration under n.meta if it is
    missing.

    Parameters
    ----------
    result_path
        Absolute or relative path to the run results folder that
        contains all model results (typically ends with "results",
        or is a time-stamp).
    sub_directory
        The subdirectory name to read files from relative to the
        result folder.

    Returns
    -------
    :
        A Dictionary that contains pypsa.Network objects as values the
        year from the end of the file name as keys.
    """
    # delayed import to prevent circular dependency error
    from evals.statistic import ESMStatistics

    if isinstance(result_path, list):
        file_paths = [Path(p) for p in result_path]  # assuming snakemake.input.networks
    else:
        input_path = Path(result_path) / sub_directory
        file_paths = input_path.glob(r"*[0-9].nc")

    networks = {}
    for file_path in file_paths:
        year = re.search(Regex.year, file_path.stem).group()
        n = pypsa.Network(file_path)
        n.statistics = ESMStatistics(n, result_path)
        n.name = f"PyPSA-AT Network {year}"
        n.year = year
        networks[year] = n

    assert networks, f"No networks found in {file_paths}."

    return networks


def read_views_config(
    func: Callable, config_override: str = "config.override.toml"
) -> dict:
    """
    Return the configuration for a view function.

    The function reads the default configuration from the
    TOML file and optionally updates it using the config
    file from the override file. The configuration returned
    is stripped down to the relevant parts that matter for the
    called view function.

    Parameters
    ----------
    func
        The view function to be called by the CLI module.
    config_override
        A file name as a string as passed to the CLI module.

    Returns
    -------
    :
        The default configuration with optional overrides from
        a second configuration file.
    """
    default_fp = resources.files("evals") / "config.default.toml"
    default = tomllib.load(default_fp.open("rb"))
    default_global = default["global"]
    default_view = default[func.__name__]

    if config_override:
        override_fp = Path(resources.files("evals")) / config_override
        override = tomllib.load(override_fp.open("rb"))
        default_global = deep_update(default_global, override["global"])

        if override_view := override.get(func.__name__, {}):
            default_view = deep_update(default_view, override_view)

    config = {"global": default_global, "view": default_view}

    logger = logging.getLogger()
    logger.debug(f"Configuration items: {config}")

    return config


def read_csv_files(
    result_path: str | Path, glob: str, sub_directory: str
) -> pd.DataFrame:
    """
    Read CSV files from disk.

    Assumes, that if the file name ends with an underscore and 4
    digits, the 4 digits represent the year. The year is prepended
    to the result dataframe index. Otherwise, the first column in the
    CSV file will be the index.

    The function caches result with the same input arguments.

    Parameters
    ----------
    result_path
        Absolute or relative path to the run results folder that
        contains all model results (typically ends with "results",
        or is a time-stamp).
    glob
        The search pattern to filter file names. The asterix can
        be used as wildcard character that matches anything.
    sub_directory
        The subdirectory name to read files from relative to the
        result folder.

    Returns
    -------
    :
        All CSV files concatenated into one DataFrame along the
        index axis.
    """
    input_path = Path(result_path) / sub_directory
    assert input_path.is_dir(), f"Input path does not exist: {input_path.resolve()}"
    file_paths = input_path.glob(glob)

    df_list = []
    for file_path in file_paths:
        _df = pd.read_csv(file_path, index_col=0)
        if year := re.search(Regex.year, file_path.stem):
            _df = insert_index_level(_df, year.group(), DataModel.YEAR)
        df_list.append(_df)

    # must assert after the loop, because file_paths is a generator
    assert df_list, f"No files named like '{glob}' in {input_path.resolve()}."

    return pd.concat(df_list, sort=True)


def get_resources_directory(n: pypsa.Network) -> Callable:
    """Return a path provider to the resources directory for a network."""
    run = n.meta["run"]
    return path_provider(
        "../resources/",  # assuming CWD is evals/cli.py
        get_rdir(run),
        run["shared_resources"]["policy"],
        run["shared_resources"]["exclude"],
    )


def _add_dummy_rows(df: pd.DataFrame, keep_regions: tuple) -> pd.DataFrame:
    """
    Add rows for missing year - country combinations.

    This is required to export empty figures. Empty figures
    are used in the VAMOS interface to show that a metric has
    no data for a country. For example, Italy has no district
    heat network and, as a result, no data in the respective
    district heat production capacities evaluation chart.

    Parameters
    ----------
    df
        The data frame with a locations index level.
    keep_regions
        The regions to add empty rows for.

    Returns
    -------
    :
        The input data frame one with additional emtpy row
        per missing country.
    """
    attrs = df.attrs
    years = df.index.unique(DataModel.YEAR)  # assuming all required years are present
    countries = list(ALIAS_COUNTRY.values())
    regions = [loc for k, loc in ALIAS_REGION.items() if k.startswith(keep_regions)]
    locations = countries + regions

    idx_names_required = DataModel.YEAR_IDX_NAMES[:2]  # year, location
    n_levels_to_add = df.index.nlevels - len(idx_names_required)
    idx_required = pd.MultiIndex.from_product(
        [years, locations], names=idx_names_required
    )

    idx_present = df.reset_index().set_index(idx_names_required).index.unique()
    idx_missing_year_loc = idx_required.difference(idx_present)

    if idx_missing_year_loc.empty:
        return df

    missing_items = [idx + ("",) * n_levels_to_add for idx in idx_missing_year_loc]
    idx_missing = pd.MultiIndex.from_tuples(missing_items, names=df.index.names)
    rows_missing = pd.DataFrame(index=idx_missing, columns=df.columns, data=pd.NA)
    result = pd.concat([rows_missing, df])
    result.attrs = attrs

    return result


class Exporter:
    """
    A class to export statistics.

    The exporter data frame consists of multiple joined statistics,
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

    def __init__(
        self,
        statistics: list,
        view_config: dict,
        keep_regions: tuple = (
            "AT",
            "GB",
            "ES",
            "FR",
            "DE",
            "IT",
        ),  # todo: move to global config
        region_nice_names: bool = True,
    ) -> None:
        self.statistics = statistics
        units = {stat.attrs["unit"] for stat in statistics}
        assert len(units) == 1, f"Mixed units cannot be exported: {units}."
        self.is_unit = units.pop()
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
        """
        Build the metric and store it as a cached property.

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
        """
        Create the plotly figure and export it as HTML and JSON.

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

        df_plot = _add_dummy_rows(df_plot, self.keep_regions)

        for idx, data in df_plot.groupby(cfg.plotby):
            chart = cfg.chart(data, cfg)
            chart.plot()
            chart.to_html(output_path, cfg.plotby, idx)
            chart.to_json(output_path, cfg.plotby, idx)

    def export_excel(self, output_path: Path) -> None:
        """
        Export metrics to Excel files for countries and regions.

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
        """
        Encode the metric da frame to a CSV file.

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

    def export(self, result_path: Path, subdir: str) -> None:
        """
        Export the metric to formats specified in the config.

        Parameters
        ----------
        result_path
            The path to the results folder.
        subdir
            The subdirectory inside the results folder to store evaluation results under.

        Returns
        -------
        :
        """
        output_path = self.make_evaluation_result_directories(result_path, subdir)

        self.export_plotly(output_path)

        if "excel" in self.view_config.get("exports", []):
            self.export_excel(output_path)
        if "csv" in self.view_config.get("exports", []):
            self.export_csv(output_path)

        # always run tests after the export
        self.consistency_checks()

    def consistency_checks(self) -> None:
        """
        Run plausibility and consistency checks on a metric.

        The method typically is called after exporting the metric.
        Unmapped categories do not cause evaluations to fail, but
        the evaluation function should return in error state to obviate
        missing entries in the mapping.

        Parameter
        ---------
        config_checks
            A dictionary with flags for every test to run.

        Returns
        -------
        :

        Raises
        ------
        AssertionError
            In case one of the checks fails.
        """
        self.default_checks()

        if "balances_almost_zero" in self.view_config.get("checks", []):
            groups = [DataModel.YEAR, DataModel.LOCATION]
            yearly_sum = self.df.groupby(groups).sum().abs()
            balanced = yearly_sum < self.view_config["cutoff"]
            if isinstance(balanced, pd.DataFrame):
                assert balanced.all().all(), (
                    f"Imbalances detected: {yearly_sum[balanced == False].dropna(how='all').sort_values(by=balanced.columns[0], na_position='first').tail()}"
                )
            else:  # Series
                assert balanced.all().item(), (
                    f"Imbalances detected: {yearly_sum[balanced.squeeze() == False].squeeze().sort_values().tail()}"
                )

    def default_checks(self) -> None:
        """Perform integrity checks for views."""
        category = self.defaults.plotly.plot_category
        categories = self.view_config["categories"]

        assert self.df.index.unique(category).isin(categories.keys()).all(), (
            f"Incomplete categories detected. There are technologies in the metric "
            f"data frame, that are not assigned to a group (nice name)."
            f"\nMissing items: "
            f"{self.df.index.unique(category).difference(categories.keys())}"
        )

        superfluous_categories = self.df.index.unique(category).difference(
            categories.keys()
        )
        assert len(superfluous_categories) == 0, (
            f"Superfluous categories found: {superfluous_categories}"
        )

        a = set(self.view_config["legend_order"])
        b = set(categories.values())
        additional = a.difference(b)
        assert not additional, (
            f"Superfluous categories defined in legend order: {additional}"
        )
        missing = b.difference(a)
        assert not missing, (
            f"Some categories are not defined in legend order: {missing}"
        )

        no_color = [c for c in categories.values() if c not in COLOUR_SCHEME]
        assert len(no_color) == 0, (
            f"Some categories used in the view do not have a color assigned: {no_color}"
        )

    def make_evaluation_result_directories(
        self, result_path: Path, subdir: Path | str
    ) -> Path:
        """
        Create all directories needed to store evaluations results.

        Parameters
        ----------
        result_path
            The path of the result folder.
        subdir
            A relative path inside the result folder.

        Returns
        -------
        :
            The joined path: result_dir / subdir.
        """
        output_path = self.make_directory(result_path, subdir)
        self.make_directory(output_path, "HTML")
        self.make_directory(output_path, "JSON")
        self.make_directory(output_path, "CSV")
        self.make_directory(output_path, "XLSX")

        return output_path

    @staticmethod
    def make_directory(base: Path, subdir: Path | str) -> Path:
        """
        Create a directory and return its path.

        Parameters
        ----------
        base
            The path to base of the new folder.
        subdir
            A relative path inside the base folder.

        Returns
        -------
        :
            The joined path: result_dir / subdir / now.
        """
        base = Path(base).resolve()
        assert base.is_dir(), f"Base path does not exist: {base}."
        directory_path = base / subdir
        directory_path.mkdir(parents=True, exist_ok=True)

        return directory_path
