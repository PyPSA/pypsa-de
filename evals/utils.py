"""Collect package helper functions."""

import logging
import re
from contextlib import contextmanager
from decimal import ROUND_HALF_UP, Decimal, localcontext
from pathlib import Path

import numpy as np
import pandas as pd

from constants import (
    ALIAS_COUNTRY,
    ALIAS_LOCATION,
    ALIAS_REGION,
    UNITS,
    Carrier,
    DataModel,
    Group,
    Mapping,
    Regex,
    TradeTypes,
)


def verify_metric_format(metric: pd.DataFrame) -> None:
    """Ensure correct metric format.

    Parameters
    ----------
    metric
        The metric data frame. This format is supported by export
        functions.

    Raises
    ------
    AssertionError
        If the metric does not comply with the data model.
    """
    assert isinstance(metric, pd.DataFrame), (
        f"Metric must be a DataFrame, " f"but {type(metric)} was passed."
    )
    assert set(metric.index.names).issubset(set(DataModel.YEAR_IDX_NAMES)), (
        f"Metric index levels must contain {DataModel.YEAR_IDX_NAMES}, "
        f"but {metric.index.names} is set."
    )
    assert metric.columns.names in ([DataModel.METRIC], [DataModel.SNAPSHOTS]), (
        f"Metric column level names must be [{DataModel.METRIC}] or "
        f"[{DataModel.SNAPSHOTS}], but {metric.columns.names} is set."
    )

    assert metric.attrs.get("name"), "Must set the metric name in 'metric.attrs'."
    assert metric.attrs.get("unit"), "Must set the metric unit in 'metric.attrs'."

    if metric.columns.names == [DataModel.METRIC]:
        assert all(
            "(" in c and ")" in c for c in metric.columns
        ), f"All columns must have a unit in braces: {metric.columns}"

        assert len(metric.columns) == 1, "Multiple aggregated metrics are not allowed."

    elif metric.columns.name == DataModel.SNAPSHOTS:
        assert isinstance(
            metric.columns, pd.DatetimeIndex
        ), "Snapshot columns must be of type DatetimeIndex."


def insert_index_level(
    df: pd.DataFrame | pd.Series,
    value: str,
    index_name: str,
    axis: int = 0,
    pos: int = 0,
) -> pd.DataFrame | pd.Series:
    """Add an index level to the data frame.

    Parameters
    ----------
    df
        The data frame that will receive the new outer level index.
    value
        The new index values.
    index_name
        The new index level name.
    axis : optional
        The index axis. Pass 0 for row index and 1 for column index.
    pos : optional
        Move the new index name to this position. 0 is outer left,
        1 is the second, and so on.

    Returns
    -------
    :
        The data frame with the new index level.
    """
    result = pd.concat({value: df}, names=[index_name], axis=axis)
    if pos == 0:  # no need to reorder levels. We are done inserting.
        return result
    idx = df.index if axis == 0 else df.columns
    idx_names = list(idx.names)
    idx_names.insert(pos, index_name)
    if isinstance(result, pd.DataFrame):
        return result.reorder_levels(idx_names, axis=axis)
    return result.reorder_levels(idx_names)


def calculate_cost_annuity(n: float, r: float | pd.Series = 0.07) -> float | pd.Series:
    """Calculate the annuity factor for an asset.

    Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20,0.05)*20 = 1.6

    Parameters
    ----------
    n
        The lifetime of the asset in years.
    r
        The discount rate of the asset.

    Returns
    -------
    :
        The calculated annuity factors.

    Notes
    -----
    This function was adopted from the abandoned package "vresutils".
    """
    if isinstance(r, pd.Series):
        ser = pd.Series(1 / n, index=r.index)
        return ser.where(r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n))
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def get_unit(s: str) -> str:
    """Parse the unit from a string.

    The unit must be inside round parentheses. If multiple
    parenthesis are found in the input string, returns the last one.

    Parameters
    ----------
    s
        The input string that should contain a unit.

    Returns
    -------
    :
        All characters inside the last pair of parenthesis without
        the enclosing parenthesis, or an empty string.
    """
    if matches := re.findall(Regex.unit, s):
        return matches[-1].strip("()")
    return ""


def get_trade_type(bus_a: str, bus_b: str) -> str:
    """Determine the trade type between two buses.

    Parameters
    ----------
    bus_a
        1st string that should start with a region substring.
    bus_b
        2nd string that should start with a region substring.

    Returns
    -------
    :
        The trade type. One of constants.TRADE_TYPES.
    """
    loc_a = re.findall(Regex.region, bus_a)[:1]
    loc_b = re.findall(Regex.region, bus_b)[:1]
    if not loc_a or not loc_b:  # no region(s) found
        return ""
    elif loc_a[0] == loc_b[0]:
        # transformation link in same region, e.g. heat
        return TradeTypes.LOCAL
    elif loc_a[0][:2] == loc_b[0][:2]:  # only country codes match
        return TradeTypes.DOMESTIC
    else:
        return TradeTypes.FOREIGN


def trade_mask(
    comp: pd.DataFrame, scopes: str | tuple, buses: tuple = ("bus0", "bus1")
) -> pd.Series:
    """Get the mask for a given trade type.

    The logic only compares bus0 and bus1 in a given component.

    Parameters
    ----------
    comp
        The component data frame. Should be one a branch_component,
        i.e. 'Line', 'Link', or 'Transformer'.
    scopes
        The trade scope(s) to match. One or multiple of 'local',
        'domestic', 'foreign'.
    buses
        Two buses to determine the trade type from. The trade type will
        be 'local', 'domestic', or 'foreign', for same location, same
        country code, or different country code, respectively.

    Returns
    -------
    :
        A pandas Series with the same index as component index and 1
        or 0 as values for match or differ, respectively.

    Raises
    ------
    ValueError
        In case the passed trade type is not supported and to prevent
        unintended string matches.
    """
    scopes = (scopes,) if isinstance(scopes, str) else scopes
    if unknown_scopes := set(scopes).difference(
        {TradeTypes.LOCAL, TradeTypes.DOMESTIC, TradeTypes.FOREIGN}
    ):
        raise ValueError(f"Invalid trade scopes detected: {unknown_scopes}.")
    df = comp[[*buses]]
    trade = df.apply(lambda row: get_trade_type(row[buses[0]], row[buses[1]]), axis=1)
    return trade.isin(scopes)


def filter_by(
    df: pd.DataFrame | pd.Series, exclude: bool = False, **kwargs: object
) -> pd.DataFrame | pd.Series:
    """Filter a data frame by key value pairs.

    Constructs a pandas query using the pandas.Index.isin() method.
    Since the pandas query API is only available for data frames,
    any passed pandas Series is converted to frame and reset to
    series.

    Parameters
    ----------
    df
        The data frame or Series to filter.
    exclude
        Set to True to exclude the filter result from the original
        data set, and return the difference.
    **kwargs
        Key=value pairs, used in the filter expression. Valid keys are
        index level names or column labels.

    Returns
    -------
    :
        The filtered data frame in the same format as the input
        dataframe.
    """
    if was_series := isinstance(df, pd.Series):
        df = df.to_frame()

    where_clauses = []
    for key, vals in kwargs.items():
        vals = [vals] if np.isscalar(vals) else vals
        where_clauses.append(f"{key} in {vals}")

    expression = " & ".join(where_clauses)
    result = df.query(expression)

    if exclude:
        result = df.drop(result.index)

    # squeeze(axis=1) to preserve index even for single rows
    return result.squeeze(axis=1) if was_series else result


def expand_to_time_series(
    df: pd.DataFrame | pd.Series, snapshots: pd.Index, nhours: int = 8760
) -> pd.DataFrame:
    """Convert time aggregated values into a time series.

    Any column label will be dropped and replaced by the given
    snapshots. It is assumed, that the metric holds yearly values, as
    produced by time aggregation methods. The data frame index and
    attrs are preserved. Time series value will become the yearly value
    divided by the number hours per year, i.e. the hourly values.

    Parameters
    ----------
    df
        A data frame input data frame with one column.
    snapshots
        The columns labels to use in the result (snapshot time stamps).
    nhours
        Divide values in the input by this number..

    Returns
    -------
    :
        The time series data frame with values average values
        representing the time interval between snapshots.

    Raises
    ------
    NotImplementedError
        If a data frame with more than one column is passed.
    """
    if isinstance(df, pd.DataFrame):
        if df.shape[1] > 1:
            raise NotImplementedError(
                f"Broadcasting multiple columns is not supported. "
                f"Only single column data frames may be passed as "
                f"input, but found {df.shape[1]} columns."
            )
        df = df.squeeze(axis=1)

    hourly = df / nhours
    values = np.tile(hourly.to_numpy(), (len(snapshots), 1)).T
    result = pd.DataFrame(index=df.index, columns=snapshots, data=values)
    result.attrs = df.attrs
    return result


def split_location_carrier(index: pd.MultiIndex, names: list) -> pd.MultiIndex:
    r"""Split location and carrier in the index.

    The location must be encoded in the string and match the regex
    '^[A-Z]{2}\\d\\s\\d'. Subsequent characters become the carrier
    name. The location defaults to an emtpy string if the regex
    does not match.

    Parameters
    ----------
    index
        A pandas Multiindex with the innermost level to split.
    names
        The list of output Multiindex names.

    Returns
    -------
    :
        The resulting Multiindex with one additional
        level due to the splitting.
    """
    idx_split = []
    for *prefix, loc_category in index:
        matches = re.match(Regex.region, loc_category)
        location = matches.group() if matches else ""
        technology = loc_category.removeprefix(location).strip()
        idx_split.append((*prefix, location, technology))

    return pd.MultiIndex.from_tuples(idx_split, names=names)


def rename_aggregate(
    df: pd.DataFrame | pd.Series,
    mapper: dict | str,
    level: str = DataModel.CARRIER,
    agg: str = "sum",
) -> pd.Series | pd.DataFrame:
    """Rename index values and aggregate duplicates.

    In case the supplied mapper is a string, all values in the
    supplied level are replaced by this string.

    Parameters
    ----------
    df
        The input data frame.
    mapper
        A Dictionary with key-value pairs to rename index values, or
        a string used to replace all values in the given level.
    level
        The index level name.
    agg
        The aggregation method for duplicated index values after
        renaming.

    Returns
    -------
    :
        A data frame with renamed index values and aggregated values.

    Notes
    -----
    Support for column axis mapping was removed, because the groupby
    operation along axis=1 removes column level names and does not
    work correctly.
    """
    if isinstance(mapper, str):
        mapper = dict.fromkeys(df.index.unique(level=level), mapper)
    renamed = df.rename(mapper, level=level)
    return renamed.groupby(df.index.names).agg(agg)


def apply_cutoff(df: pd.DataFrame, limit: float, drop: bool = True) -> pd.DataFrame:
    """Replace small absolute values with NaN.

    The limit boundary is not inclusive, i.e. the limit value itself
    will not be replaced by NaN.

    Parameters
    ----------
    df
        The data frame to remove values from.
    limit
        Absolute values smaller than the limit will be dropped.
    drop
        Whether to drop all NaN rows from the returned data frame.

    Returns
    -------
    :
        A data frame without values that are smaller than the limit.
    """
    result = df.mask(cond=df.abs() <= abs(limit), other=pd.NA)
    if drop:
        result = result.dropna(how="all", axis=0)
    return result


def get_mapping(map_name: str, map_type: str = "external") -> dict:
    """Extract a dictionary from the nested mapping definition.

    Parameters
    ----------
    map_name
        The name of the mapping component to extract. Different
        evaluation contexts require different grouping of carrier.

    map_type
        The mapping type to extract. May be "internal" or "external",
        defaults to "external". The internal mapping returns groups
        of finer granularity than external.

    Returns
    -------
    :
        The combined mapping for the requested components.
    """
    # some mappings are provided directly in the evals functions.
    if isinstance(map_name, dict):
        return map_name

    result = {}
    for carrier, group_names in Mapping.carrier.items():
        if isinstance(group_names, str):  # single entry
            result[carrier] = group_names
        elif isinstance(group_names, dict):
            group = group_names.get(map_type)
            if isinstance(group, str):
                result[carrier] = group
            elif isinstance(group, dict):
                value = group.get(map_name, group["default"])
                # skips renaming altogether if value is empty.
                result[carrier] = value or carrier

    return result


def aggregate_eu(df: pd.DataFrame, agg: str = "sum") -> pd.DataFrame:
    """Calculate the EU region as the sum of all country regions.

    The carrier 'import net', 'export net', 'Import European' and '
    Export European' need to be removed from the EU data set.
    The total import and export over all countries evens out and
    is not required for EU location. The non-EU imports
    are named differently, e.g. 'global import'.

    Parameters
    ----------
    df
        The data frame with one MultiIndex level named 'location'.
    agg
        The aggregation function.

    Returns
    -------
    :
        Summed metric with one location named 'EU'.
    """
    df = df.query(f"{DataModel.LOCATION} not in ['EU', '']")  # valid countries only
    totals = rename_aggregate(df, "EU", level=DataModel.LOCATION, agg=agg)
    excluded = [
        Group.import_net,  # required for CH4 and H2!
        Group.export_net,
        Group.import_european,
        Group.export_european,
        # exclude domestic trade for EU region
        Group.import_domestic,
        Group.export_domestic,
        Carrier.import_domestic,
        Carrier.export_domestic,
    ]
    return totals.drop(excluded, level=DataModel.CARRIER, errors="ignore")


def aggregate_locations(
    df: pd.DataFrame,
    keep_regions: tuple = ("AT",),
    nice_names: bool = True,
) -> pd.DataFrame:
    """Aggregate to countries, including EU and keeping certain regions.

    The input data frame is expected to contain locations as regions,
    e.g. "AT0 1", "FR0 0", etc.

    Parameters
    ----------
    df
        The input data frame with a locations index level.
    keep_regions
        A tuple of regions, that should be preserved in the output,
        i.e. they are added to the result as before the aggregation.
    nice_names
        Whether, or not to use the nice country names instead of the
        country codes.

    Returns
    -------
    :
        A data frame with aggregated countries, plus any region in
        'keep_regions' and Europe/EU.
    """
    country_code_map = {loc: loc[:2] for loc in df.index.unique(DataModel.LOCATION)}
    if "EU" in country_code_map.values():
        logger = logging.getLogger(__name__)
        logger.warning(
            "Values for 'EU' node found in input data frame. "
            "This can lead to value doubling during location aggregation.",
        )
    countries = rename_aggregate(df, country_code_map, level=DataModel.LOCATION)
    # domestic trade only makes sense between regions. Aggregated
    # countries could have domestic trade, but import and export nets
    # to zero.
    countries = countries.drop(
        [
            Carrier.export_domestic,
            Carrier.import_domestic,
            Group.import_domestic,
            Group.export_domestic,
        ],
        level=DataModel.CARRIER,
        errors="ignore",
    )
    europe = aggregate_eu(df)
    mask = df.index.get_level_values(DataModel.LOCATION).str.startswith(keep_regions)
    regions = df.loc[mask, :]
    result = pd.concat([countries, regions, europe]).sort_index(axis=0)
    if nice_names:
        result = result.rename(index=ALIAS_LOCATION, level=DataModel.LOCATION)
    return result


def add_dummy_rows(df: pd.DataFrame, keep_regions: tuple) -> pd.DataFrame:
    """Add rows for missing year - country combinations.

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


def scale(df: pd.DataFrame, to_unit: str) -> pd.DataFrame:
    """Scale metric values to the specified target unit.

    Multiplies all columns in the metric by a scaling factor.
    The scaling factor is calculated from the unit in the data frame
    columns and the given target unit. Also updates the unit
    names encoded in the data frame columns for time aggregated
    metrics.

    Parameters
    ----------
    df
        The input data frame with valid units in the column labels.
    to_unit
        The target unit. See constants.UNITS for possible
        units.

    Returns
    -------
    :
        The scaled data frame with replaced units in column labels.

    Raises
    ------
    raises KeyError
        If the 'to_unit' is not found in UNITS, or if the attrs
        dictionary has no unit field.
    raises ValueError
        If input units are inconsistent, i.e. mixed power and energy
        columns.
    """
    if df.columns.name == DataModel.SNAPSHOTS:
        is_unit = df.attrs["unit"]
        scaling_factor = is_unit / to_unit
        result = df.mul(scaling_factor)
    else:
        scale_to = to_unit if isinstance(to_unit, float) else UNITS[to_unit]
        units_in = list(map(get_unit, df.columns))
        if to_unit.endswith("h") and not all(u.endswith("h") for u in units_in):
            raise ValueError("Denying to convert units from power to energy.")
        if to_unit.endswith("W") and not all(u.endswith("W") for u in units_in):
            raise ValueError("Denying to convert unit from energy to power.")
        scale_in = [UNITS[s] for s in units_in]
        scaling_factors = [x / scale_to for x in scale_in]

        result = df.mul(scaling_factors, axis=1)
        result.columns = result.columns.str.replace(
            "|".join(units_in), to_unit, regex=True
        )

    result.attrs["unit"] = to_unit

    return result


def drop_from_multtindex_by_regex(
    df: pd.DataFrame, pattern: str, level: str = DataModel.CARRIER
) -> pd.DataFrame | pd.Series:
    """Drop all rows that match the regex in the index level.

    This function is needed, because pandas.DataFrame.filter cannot
    be applied to MultiIndexes.

    Parameters
    ----------
    df
        The input data frame with a multi index.
    pattern
        The regular expression pattern as a raw string.
    level
        The multi index level to match the regex to.

    Returns
    -------
    :
        The input data where the regular expression does not match.
    """
    mask = df.index.get_level_values(level).str.contains(pattern, regex=True)
    return df[~mask]


@contextmanager
def operations_override(networks: dict, component: str, operation: str) -> None:
    """Patch the used operations time series.

    Note, that monkeypatching does not work anymore since PyPSA
    >0.30, because the 'get_operations' function is not used
    by the pypsa.statistics anymore.

    Parameters
    ----------
    networks
        The PyPSA network dictionary.
    component
        The component to patch, e.g. Link, Store, etc.
    operation
        The desired operations time series to use instead of 'p'.

    Yields
    ------
    None, passes to the with statement block.
    """
    _temp_key = "_tmp"

    for n in networks.values():
        c = n.pnl(component)
        c[_temp_key] = c["p"]  # save a copy
        c["p"] = c[operation]  # overwrite

    yield  # run anything in the with statement

    for n in networks.values():
        c = n.pnl(component)
        c["p"] = c.pop(_temp_key)  # restore original


def prettify_number(x: float) -> str:
    """Format a float for display on trace hover actions.

    Parameters
    ----------
    x
        The imprecise value to format.

    Returns
    -------
    :
        The formatted number as a string with 1 or 0 decimal places,
        depending on the magnitude of the input value.
    """
    if abs(round(x, 0)) >= 10:
        with localcontext():
            return str(round(round(Decimal(x), 1), 0))
    else:
        with localcontext() as ctx:
            ctx.rounding = ROUND_HALF_UP
            return str(round(round(Decimal(x), 2), 1))


def make_evaluation_result_directories(result_path: Path, subdir: Path | str) -> Path:
    """Create all directories needed to store evaluations results.

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
    output_path = make_directory(result_path, subdir)
    make_directory(output_path, "HTML")
    make_directory(output_path, "JSON")
    make_directory(output_path, "CSV")

    return output_path


def make_directory(base: Path, subdir: Path | str) -> Path:
    """Create a directory and return its path.

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


def swap_multiindex_values(
    df: pd.DataFrame, lvl0: str | int, lvl1: str | int
) -> pd.DataFrame:
    """Swap values in two multiindex index levels.

    Parameters
    ----------
    df
        The inout data frame with the two index level in its MultiIndex.
    lvl0
        The first multiindex level name.
    lvl1
        The second multiindex level name.

    Returns
    -------
    :
        A data frame with the same MultiIndex levels as the input,
        but the contained MultiIndex values are replaced by each other.
    """
    idx_names = list(df.index.names)
    pos0, pos1 = idx_names.index(lvl0), idx_names.index(lvl1)
    # idx_names[pos1], idx_names[pos0] = idx_names[pos0], idx_names[pos1]
    # For the most curious and trust scattering reasons, this line
    # needs to run twice.
    # idx_names[pos1], idx_names[pos0] = idx_names[pos0], idx_names[pos1]
    df_reversed_values = df.swaplevel(lvl0, lvl1)
    df_reversed_values.index.names = idx_names

    return df_reversed_values
