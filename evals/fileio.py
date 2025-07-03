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
from evals.constants import COLOUR_SCHEME_BMK, NOW, TITLE_SUFFIX, DataModel, Regex
from evals.excel import export_excel_countries, export_excel_regions_at
from evals.utils import (
    add_dummy_rows,
    calculate_cost_annuity,
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


def prepare_costs(
    result_path: str | Path,
    n_years: int,
    n: pypsa.Network,
    raw: bool = False,
    sub_directory: str = "esm_run/interpolated_data",
) -> pd.DataFrame:
    """
    Read cost files and calculate fixed costs.

    Parameters
    ----------
    result_path
        Absolute or relative path to the run results folder that
        contains all model results (typically ends with "results",
        or is a time-stamp).
    n_years
        The number of years used to calculate the annual costs.
    n
        The Notwork with the meta attribute to obtain the config.
    raw
        Whether to return the CSV data as is, i.e. without added cost
        calculations.
    sub_directory
        The location of the costs files relative to the results folder
        (Only required during testing).

    Returns
    -------
    :
        Costs by year, technology and parameter.
    """
    costs = read_csv_files(result_path, "costs_*.csv", sub_directory)
    costs = costs.set_index("parameter", append=True, drop=True)
    if raw:
        return costs

    default_costs = {
        "CO2 intensity": 0,
        "FOM": 0,
        "VOM": 0,
        "discount rate": n.meta["costs"]["discountrate"],
        "efficiency": 1,
        "fuel": 0,
        "investment": 0,
        "lifetime": n.meta["costs"]["lifetime"],
    }
    costs.loc[costs["unit"].str.contains("/kW"), "value"] *= 1e3  # to 1/MW
    usd_to_eur = n.meta["costs"]["USD2013_to_EUR2013"]
    costs.loc[costs["unit"].str.contains("USD"), "value"] *= usd_to_eur

    costs = costs["value"]
    costs = costs.unstack(level="parameter")
    costs = costs.groupby(costs.index.names).sum(min_count=1)  # keep NaNs
    costs = costs.fillna(default_costs)

    def _calculate_fixed_costs(row: pd.Series) -> pd.Series:
        """
        Calculate the fixed costs per row.

        Parameters
        ----------
        row
            A row in the costs data frame.

        Returns
        -------
        :
            The fixed costs for the row.
        """
        annuity = calculate_cost_annuity(row["lifetime"], row["discount rate"])
        annuity += row["FOM"] / 100.0
        return annuity * row["investment"] * n_years

    costs["fixed"] = costs.apply(_calculate_fixed_costs, axis=1)

    # align index names with data model
    costs.index.names = [DataModel.YEAR, DataModel.CARRIER]

    return costs


def prepare_co2_emissions(
    n: pypsa.Network,
) -> pd.DataFrame:
    """
    Read emissions from file and calculate the country share.

    The function reads the emissions from file. It also uses the run
    configuration to determine the cluster settings. If the cluster is
    named 'AT10', the region population shares are read from file.
    The region population shares split the country level emissions into
    region emissions.

    Parameters
    ----------
    n
        The Notwork with the meta attribute to obtain the config.

    Returns
    -------
    :
        The CO2 emission share by country.

    Notes
    -----
    This function needs to be updated for Agriculture emissions and
    for Austrian NUTS2 regions shares, that are not accurate otherwise.
    """
    run = n.meta["run"]
    res = get_resources_directory(n)
    co2 = read_csv_files(res(""), "co2_totals.csv", "")

    options = n.meta["wildcards"]["sector_opts"]
    cols = ["electricity"]
    if "T" in options or options == "none":  # Transport
        cols.extend(["rail non-elec", "road non-elec"])
    if "H" in options or options == "none":  # Households
        cols.extend(["residential non-elec", "services non-elec"])
    if "I" in options or options == "none":  # Industry
        cols.extend(
            [
                "industrial non-elec",
                "industrial processes",
                "domestic aviation",
                "domestic navigation",
            ]
        )
    if "A" in options or options == "none":  # Agriculture
        cols.extend(["agriculture"])

    co2 = co2[cols].sum(axis=1)
    co2_dist = co2 / co2.sum()
    cluster = n.meta["wildcards"]["clusters"]
    # disaggregate country emissions by nodal population share
    pop = read_csv_files(
        res(run["prefix"]), f"pop_layout_base_s_{cluster}.csv", run["name"][0]
    )

    def population_fraction(ct):
        return ct["fraction"] * co2_dist[ct.name]

    co2_dist_pop = pop.groupby("ct", group_keys=False).apply(
        population_fraction, include_groups=False
    )  # .droplevel("ct")
    co2_dist_pop.index = co2_dist_pop.index.set_names(DataModel.LOCATION)

    return co2_dist_pop


def prepare_industry_demand(
    result_path: str | Path,
    networks: dict,
    n: pypsa.Network,
    sub_directory: str = "interpolated_data",
    rename: bool = True,
) -> pd.Series:
    """
    Read industry demand from resource file and correct methane.

    The industry demand for methane must be increased by the carbon
    capture amounts. The 'gas for industry CC' Links have efficiencies
    smaller than 1. Therefore, resulting losses must be added to the
    methane values in the resource file.

    The industry demand is read from resource files, and not
    extracted from the network, because the network contains
    demands in an aggregated form and the CSV files are on a per-sector
    granularity.

    Parameters
    ----------
    result_path
        Absolute or relative path to the run results folder that
        contains all model results. (typically ends with "results",
        or is a time-stamp).
    networks
        The pypsa networks, used to extract the methane amounts
        including CC.
    n
        The Notwork with the meta attribute to obtain the config.
    sub_directory : optional
        The location of the CSV files relative to the results folder.
    rename : optional
        Whether, or not to rename sectors and apply aggregate groups.

    Returns
    -------
    :
        The industry demand by year, country, and sector group in MWh.
    """
    # delayed import to prevent circular import error
    from statistic import collect_myopic_statistics

    simpl = n.meta["scenario"]["simpl"][0]
    cluster = n.meta["scenario"]["clusters"][0]

    file_name_pattern = f"industrial_demand_by_sector_elec_s{simpl}_{cluster}_*.csv"
    industry_demand = read_csv_files(result_path, file_name_pattern, sub_directory)
    industry_demand = industry_demand.set_index("sector", append=True)

    if rename:
        industry_sectors = {
            "Integrated steelworks": "Iron and Steel",
            "Electric arc": "Iron and Steel",
            "Cement": "Non-metallic Minerals",
            "Ceramics & other NMM": "Non-metallic Minerals",
            "Glass production": "Non-metallic Minerals",
            "Basic chemicals": "Chemical",
            "Other chemicals": "Chemical",
            "Pharmaceutical products etc.": "Chemical",
            "Pulp production": "Paper, Pulp and Print",
            "Paper production": "Paper, Pulp and Print",
            "Printing and media reproduction": "Paper, Pulp and Print",
            "Alumina production": "Non-ferrous Metals",
            "Aluminium - primary production": "Non-ferrous Metals",
            "Aluminium - secondary production": "Non-ferrous Metals",
            "Other non-ferrous metals": "Non-ferrous Metals",
            "Food, beverages and tobacco": "Food and Tabacco",
            "Transport Equipment": "Transport Equipment",
            "Machinery Equipment": "Machinery",
            "Textiles and leather": "Textile and Leather",
            "Wood and wood products": "Wood and Wood Products",
            "Other Industrial Sectors": "Non-specified Industry",
        }
        industry_demand = rename_aggregate(
            industry_demand, industry_sectors, level="sector"
        )

    industry_demand = industry_demand * 1e6  # to MWh base unit

    gas_for_industry = collect_myopic_statistics(
        networks,
        "withdrawal",
        comps="Link",
        bus_carrier="gas",
        carrier=["gas for industry", "gas for industry CC"],
        drop_zero_rows=False,
        at_port=False,  # only port 0
    )

    # combine gas demand (CC and non-CC) by naming them the same
    mapper = {"gas for industry CC": "gas for industry"}
    gas_for_industry = rename_aggregate(
        gas_for_industry, mapper, level=DataModel.CARRIER
    )

    def increase_gas_demand_by_cc_share(ser: pd.Series) -> pd.Series:
        """
        Correct the gas demand for industry by CC fraction.

        The value in the resource file does not include energy for
        carbon capture. The sector values must be increased by
        the amount of energy used by carbon capture (if any).

        Parameters
        ----------
        ser
            The input Series with gas demand values.

        Returns
        -------
        :
            Gas demands increased by the energy amounts needed for
            carbon capture.
        """
        demand_cc = gas_for_industry[pd.IndexSlice[ser.name]].to_numpy()[0]
        demand = ser.sum()
        scaling_factor = demand_cc / demand
        return ser if demand < 1 else ser * scaling_factor

    industry_demand["methane"] = (
        industry_demand["methane"]
        .groupby(["year", "country"], group_keys=False)
        .apply(increase_gas_demand_by_cc_share)
    )

    # rename index levels for consistency with the statistics data model
    # The 'sector' in the industry demand interpolated data file becomes
    # the DataModel 'bus_carrier'. The unnamed column (technologies?)
    # becomes the DataModel 'carrier'
    industry_demand = industry_demand.stack().swaplevel(-1, -2)
    industry_demand.index.names = [DataModel.YEAR] + DataModel.IDX_NAMES

    return industry_demand


def prepare_nodal_energy(
    result_path: str | Path, sub_directory: str = "../resources"
) -> pd.DataFrame:
    """
    Prepare nodal energy data for analysis.

    Since bus_carrier information is not included in the nodal energy
    CSV files, it is explicitly added here through a mapping of CSV
    file column labels to bus_carrier names.

    Note that all values returned are negative, although demand usually
    has negative sign in the evaluation routines.

    Parameters
    ----------
    result_path
        Absolute or relative path to the results folder in the project
        root directory.
    sub_directory
        The subdirectory relative to the result path where the CSV
        files are located.

    Returns
    -------
    :
        Processed nodal energy data with aggregated locations and
        scaled to MWh.

    Notes
    -----
    The result does not contain bus_carrier information. This should
    be added here
    """
    try:
        nodal_energy = read_csv_files(
            result_path, "nodal_energy_totals_*.csv", sub_directory=sub_directory
        )  # in TWh
    except AssertionError:
        nodal_energy = read_csv_files(
            result_path,
            "nodal_energy_totals_*.csv",
            sub_directory=sub_directory.replace("../", "./"),
        )  # in TWh

    nodal_energy = nodal_energy.stack()
    nodal_energy.index.names = DataModel.YEAR_IDX_NAMES[:3]
    nodal_energy = nodal_energy * 1e6  # hard cast to MWh

    carrier_to_bus_carrier = {
        "BEV road freight LNF": "AC",
        "BEV road freight Lkw<12t": "AC",
        "BEV road freight Lkw>12t": "AC",
        "BEV road passenger": "AC",
        "Benzin HEV road passenger": "Fischer-Tropsch",
        "Benzin PHEV road passenger": "Fischer-Tropsch",
        "Benzin road passenger": "Fischer-Tropsch",
        "CNG HEV road passenger": "gas",
        "CNG PHEV road passenger": "gas",
        "CNG road freight LNF": "gas",
        "CNG road freight Lkw<12t": "gas",
        "CNG road freight Lkw>12t": "gas",
        "CNG road passenger": "gas",
        "Diesel HEV road passenger": "Fischer-Tropsch",
        "Diesel PHEV road freight LNF": "Fischer-Tropsch",
        "Diesel PHEV road passenger": "Fischer-Tropsch",
        "Diesel road freight LNF": "Fischer-Tropsch",
        "Diesel road freight Lkw<12t": "Fischer-Tropsch",
        "Diesel road freight Lkw>12t": "Fischer-Tropsch",
        "Diesel road passenger": "Fischer-Tropsch",
        "Diesel-electric hybrid road freight Lkw>12t": "Fischer-Tropsch",
        "Fischer-Tropsch Extra-EU aviation passenger": "Fischer-Tropsch",
        "Fischer-Tropsch Intra-EU aviation passenger": "Fischer-Tropsch",
        "Fischer-Tropsch domestic aviation": "Fischer-Tropsch",
        "Fischer-Tropsch domestic aviation freight": "Fischer-Tropsch",
        "Fischer-Tropsch domestic aviation passenger": "Fischer-Tropsch",
        "Fischer-Tropsch domestic navigation": "Fischer-Tropsch",
        "Fischer-Tropsch domestic navigation freight": "Fischer-Tropsch",
        "Fischer-Tropsch domestic navigation passenger": "Fischer-Tropsch",
        "Fischer-Tropsch heavy duty road freight": "Fischer-Tropsch",
        "Fischer-Tropsch international aviation": "Fischer-Tropsch",
        "Fischer-Tropsch international aviation freight": "Fischer-Tropsch",
        "Fischer-Tropsch international aviation passenger": "Fischer-Tropsch",
        "Fischer-Tropsch international navigation": "Fischer-Tropsch",
        "Fischer-Tropsch international navigation freight": "Fischer-Tropsch",
        "Fischer-Tropsch international navigation passenger": "Fischer-Tropsch",
        "Fischer-Tropsch light duty road freight": "Fischer-Tropsch",
        "Fischer-Tropsch other road passenger": "Fischer-Tropsch",
        "Fischer-Tropsch passenger cars": "Fischer-Tropsch",
        "Fischer-Tropsch rail": "Fischer-Tropsch",
        "Fischer-Tropsch rail freight": "Fischer-Tropsch",
        "Fischer-Tropsch rail passenger": "Fischer-Tropsch",
        "Fischer-Tropsch road": "Fischer-Tropsch",
        "Fischer-Tropsch road freight": "Fischer-Tropsch",
        "Fischer-Tropsch road passenger": "Fischer-Tropsch",
        "Fischer-Tropsch two-wheel": "Fischer-Tropsch",
        "H2 FCV road freight LNF": "H2",
        "H2 FCV road freight Lkw<12t": "H2",
        "H2 FCV road freight Lkw>12t": "H2",
        "H2 FCV road passenger": "H2",
        "H2 domestic aviation": "H2",
        "H2 domestic navigation": "H2",
        "H2 domestic navigation freight": "H2",
        "H2 domestic navigation passenger": "H2",
        "H2 international aviation": "H2",
        "H2 international navigation": "H2",
        "H2 international navigation freight": "H2",
        "H2 international navigation passenger": "H2",
        "H2 rail": "H2",
        "H2 rail passenger": "H2",
        "H2 road": "H2",
        "H2 road freight": "H2",
        "H2 road freight LNF": "H2",
        "H2 road freight Lkw<12t": "H2",
        "H2 road freight Lkw>12t": "H2",
        "H2 road passenger": "H2",
        "LH2 Extra-EU aviation passenger": "LH2",
        "LH2 Intra-EU aviation passenger": "LH2",
        "LH2 domestic aviation freight": "LH2",
        "LH2 domestic aviation passenger": "LH2",
        "LH2 domestic navigation freight": "LH2",
        "LH2 domestic navigation passenger": "LH2",
        "LH2 international aviation freight": "LH2",
        "LH2 international aviation passenger": "LH2",
        "LH2 international navigation freight": "LH2",
        "LH2 international navigation passenger": "LH2",
        "LH2 rail freight": "LH2",
        "LNG domestic navigation freight": "gas",
        "LNG international navigation freight": "gas",
        "LNG road freight Lkw<12t": "gas",
        "LNG road freight Lkw>12t": "gas",
        "biogas light duty road freight": "gas",
        "biogas other road passenger": "gas",
        "biogas passenger cars": "gas",
        "biogas road": "gas",
        "electricity light duty road freight": "AC",
        "electricity other road passenger": "AC",
        "electricity passenger cars": "AC",
        "electricity rail": "AC",
        "electricity rail freight": "AC",
        "electricity rail passenger": "AC",
        "electricity residential": "AC",
        "electricity residential cooking": "AC",
        "electricity residential space": "AC",
        "electricity residential water": "AC",
        "electricity road": "AC",
        "electricity road freight": "AC",
        "electricity road passenger": "AC",
        "electricity services": "AC",
        "electricity services cooking": "AC",
        "electricity services space": "AC",
        "electricity services water": "AC",
        "gas domestic navigation": "gas",
        "gas domestic navigation passenger": "gas",
        "gas international navigation": "gas",
        "gas international navigation passenger": "gas",
        "gas light duty road freight": "gas",
        "gas other road passenger": "gas",
        "gas passenger cars": "gas",
        "gas road": "gas",
        "gas road freight": "gas",
        "gas road passenger": "gas",
        "oil rail": "Fischer-Tropsch",
        "oil road": "Fischer-Tropsch",
        "renewables rail": "AC",
        "renewables road": "AC",
        "total aviation freight": "mixed",
        "total aviation passenger": "mixed",
        "total domestic aviation": "mixed",
        "total domestic aviation freight": "mixed",
        "total domestic aviation passenger": "mixed",
        "total domestic navigation": "mixed",
        "total electricity equivalent passenger": "mixed",
        "total heavy duty road freight": "mixed",
        "total international aviation": "mixed",
        "total international aviation freight": "mixed",
        "total international aviation passenger": "mixed",
        "total international navigation": "mixed",
        "total light duty road freight": "mixed",
        "total other road passenger": "mixed",
        "total passenger cars": "mixed",
        "total rail": "mixed",
        "total rail freight": "mixed",
        "total rail passenger": "mixed",
        "total residential": "mixed",
        "total residential cooking": "mixed",
        "total residential space": "mixed",
        "total residential water": "mixed",
        "total road": "mixed",
        "total road freight": "mixed",
        "total road passenger": "mixed",
        "total services": "mixed",
        "total services cooking": "mixed",
        "total services space": "mixed",
        "total services water": "mixed",
        "total two-wheel": "mixed",
    }

    carrier = nodal_energy.index.get_level_values("carrier")
    nodal_energy = nodal_energy.to_frame()
    nodal_energy["bus_carrier"] = carrier.map(carrier_to_bus_carrier)
    nodal_energy = nodal_energy.set_index("bus_carrier", append=True, drop=True)
    nodal_energy = nodal_energy.squeeze()

    nodal_energy.attrs["name"] = "Energy Totals"
    nodal_energy.attrs["unit"] = "MWh"

    return nodal_energy


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
        keep_regions: tuple = ("AT",),
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

        df_plot = add_dummy_rows(df_plot, self.keep_regions)

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

        no_color = [c for c in categories.values() if c not in COLOUR_SCHEME_BMK]
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


def get_resources_directory(n: pypsa.Network) -> Callable:
    """Return a path provider to the resources directory for a network."""
    run = n.meta["run"]
    return path_provider(
        "../resources/",  # assuming CWD is evals/cli.py
        get_rdir(run),
        run["shared_resources"]["policy"],
        run["shared_resources"]["exclude"],
    )


def read_pypsa_csv(
    result_path: str | Path, kind: str, index_cols: int, sub_directory: str = "csvs"
) -> pd.Series:
    """
    Read nodal energy totals from CSV files.

    Parameters
    ----------
    result_path
        Absolute or relative path to the results folder in the project
        root directory.
    kind
        The kind of nodal energy totals to read, e.g. "energy_balance".
    index_cols
        The number of index columns in the CSV file.
    sub_directory
        The subdirectory relative to the result path where the CSV
        files are located.

    Returns
    -------
    :
        Nodal energy totals as a Series.
    """
    df = pd.read_csv(
        Path(result_path) / sub_directory / f"{kind}.csv",
        index_col=list(range(index_cols)),
        header=0,
        skiprows=[0, 1, 3],
    )
    df.index.names = [
        DataModel.COMPONENT,
        DataModel.CARRIER,
        DataModel.LOCATION,
        DataModel.BUS_CARRIER,
    ]
    df.columns.name = DataModel.YEAR
    ser = df.stack().reorder_levels(  # .pivot_table(index=DataModel.IDX_NAMES, aggfunc="sum")
        [
            DataModel.YEAR,
            DataModel.COMPONENT,
            DataModel.LOCATION,
            DataModel.CARRIER,
            DataModel.BUS_CARRIER,
        ]
    )

    name = kind
    unit = "undefined"

    ser.name = f"{name} ({unit})"
    ser.attrs["name"] = name
    ser.attrs["unit"] = unit

    return ser
