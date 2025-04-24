import numpy as np
import pandas as pd
import pytest
from esmtools.constants import DataModel
from esmtools.plots._base import ESMChart
from esmtools.plots.timeseries import ESMTimeSeriesChart
from esmtools.utils import (
    aggregate_locations,
    apply_cutoff,
    calculate_cost_annuity,
    expand_to_time_series,
    filter_by,
    get_trade_type,
    get_unit,
    insert_index_level,
    prettify_number,
    rename_aggregate,
    scale,
    split_location_carrier,
    trade_mask,
    verify_metric_format,
)


@pytest.fixture(scope="module")
def simple_data_frame():
    """Produce a simple data frame."""
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]}, index=["1", "2"])
    df.columns.name = "columns"
    df.index.name = "index"
    return df


@pytest.fixture(scope="module")
def df_buses():
    """Produce a valid multiport data frame."""
    return pd.DataFrame(
        {
            "bus0": ["DE0 0", "GB0 0", "FR0 0"],
            "bus1": ["DE0 0", "GB1 0", "GB0 0"],
        }
    )


@pytest.fixture(scope="module")
def df_multi_index():
    """Produce a simple data frame with a multiindex."""
    midx = pd.MultiIndex.from_product(
        [["A", "B"], ["1", "2", "3"]], names=["idx1", "idx2"]
    )
    return pd.DataFrame(data={"col": [*range(6)]}, index=midx)


@pytest.fixture(scope="module")
def ser_multi_index():
    """Produce a simple series with a multiindex."""
    midx = pd.MultiIndex.from_product(
        [["A", "B"], ["1", "2", "3"]], names=["idx1", "idx2"]
    )
    return pd.Series(data=[*range(6)], index=midx, name="col")


@pytest.fixture(scope="module")
def df_sort():
    """Produce a data frame for sorting."""
    return pd.DataFrame({"A": ["a", "b", "c", "d"]})


def df_locations():
    """Produce a data frame to test location aggregation."""
    idx = [
        ("FR0 0", "A"),
        ("AT0 1", "A"),
        ("AT0 2", "A"),
        ("AT0 2", "B"),
        ("DE0 1", "A"),
        ("DE0 2", "A"),
        ("", "A"),
    ]
    return pd.DataFrame(
        {"a": [*range(len(idx))]},
        index=pd.MultiIndex.from_tuples(idx, names=["location", "carrier"]),
    )


def df_metric(name="foo", unit="(bar)"):
    """Return a data frame with correct metric format."""
    midx = pd.MultiIndex.from_tuples(
        [("a", "b", "c", "d")], names=DataModel.YEAR_IDX_NAMES
    )
    df = pd.DataFrame({f"{name} {unit}": 1}, index=midx)
    df.columns.name = DataModel.METRIC
    if name:
        df.attrs["name"] = name
    if unit:
        df.attrs["unit"] = unit
    return df


@pytest.mark.unit
@pytest.mark.parametrize(
    ("axis", "pos", "expected_index"),
    [
        (
            0,
            0,
            pd.MultiIndex.from_tuples(
                [("VAL", "1"), ("VAL", "2")], names=["NAME", "index"]
            ),
        ),
        (
            1,
            0,
            pd.MultiIndex.from_tuples(
                [("VAL", "A"), ("VAL", "B")], names=["NAME", "columns"]
            ),
        ),
        (
            0,
            1,
            pd.MultiIndex.from_tuples(
                [("1", "VAL"), ("2", "VAL")], names=["index", "NAME"]
            ),
        ),
        (
            1,
            1,
            pd.MultiIndex.from_tuples(
                [("A", "VAL"), ("B", "VAL")], names=["columns", "NAME"]
            ),
        ),
    ],
)
def test_insert_index_level(axis, pos, expected_index, simple_data_frame):
    """
    Test the insert_index_level function.

    This function tests the behavior of the insert_index_level function
    to ensure it correctly modifies the index or columns of a DataFrame
    or Series based on the provided parameters.
    """
    result = insert_index_level(simple_data_frame, "VAL", "NAME", axis, pos)
    idx = result.index if axis == 0 else result.columns
    assert idx.equals(expected_index)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("n", "r", "expected"),
    [
        (10, 0.07, 0.14237750272736466),
        (15, 0.03, 0.08376658046228799),
        (1, 0.05, 1.049999999999999),  # edge case: one year (seems odd to me!)
        (20, 0, 0.05),  # edge case: zero rate
    ],
)
def test_calculate_cost_annuity(n, r, expected):
    """
    Test the calculate_cost_annuity function.

    This test verifies that the calculate_cost_annuity function
    correctly calculates the annuity factor for various input values
    of asset lifetime and discount rate. The expected results are
    compared against the actual output from the function.
    """
    annuity = calculate_cost_annuity(n, r)
    assert annuity == expected


@pytest.mark.unit
@pytest.mark.parametrize(
    ("n", "r", "expected"),
    [
        (
            20,
            pd.Series([0.01, 0.02, 0.03, 0.04, 0.06]),
            pd.Series(
                [
                    0.05541531489055132,
                    0.061156718125290346,
                    0.06721570759685909,
                    0.07358175032862885,
                    0.0871845569768514,
                ]
            ),
        )
    ],
)
def test_calculate_cost_annuity_series(n, r, expected):
    """
    Test the calculate_cost_annuity function for error cases.

    This test verifies that the calculate_cost_annuity function raises
    the appropriate exceptions when invalid input values are provided.
    """
    annuity = calculate_cost_annuity(n, r)
    pd.testing.assert_series_equal(annuity, expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("n", "r", "expected"),
    [(0, 0.05, ZeroDivisionError)],  # error: zero year
    ids=["zero_year_raises_ZeroDivisionError"],
)
def test_calculate_cost_annuity_fails(n, r, expected):
    """
    Test the calculate_cost_annuity function for error cases.

    This test verifies that the calculate_cost_annuity function raises
    the appropriate exceptions when invalid input values are provided.
    In this specific case, it checks that providing a lifetime of zero
    years raises a ZeroDivisionError.
    """
    with pytest.raises(expected):
        calculate_cost_annuity(n, r)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("input_str", "expected_output"),
    [
        # Happy path tests
        pytest.param("prefix (C)", "C", id="single_unit"),
        pytest.param("prefix (°C)", "°C", id="unit_with_unicode"),
        pytest.param("prefix (kW/€)", "kW/€", id="unit_with_slash"),
        pytest.param("prefix (MW) suffix", "MW", id="unit_with_prefix_and_suffix"),
        pytest.param("(TWh) suffix", "TWh", id="unit_with_suffix"),
        pytest.param("No unit here", "", id="no_parentheses"),
        pytest.param("Multiple ($) in (€) string", "€", id="multiple_parentheses"),
        pytest.param("Empty parentheses ()", "", id="empty_parentheses"),
        pytest.param("Nested (parentheses (inner))", "inner", id="nested_parentheses"),
        pytest.param("Trailing (unit) text", "unit", id="trailing_text"),
        pytest.param("Only opening (parenthesis", "", id="only_opening_parenthesis"),
        pytest.param("Only closing parenthesis)", "", id="only_closing_parenthesis"),
        pytest.param("Mismatched )parentheses(", "", id="mismatched_parentheses"),
    ],
)
def test_get_unit(input_str, expected_output):
    """
    Test the get_unit function.

    This test verifies that the get_unit function correctly extracts the
    unit from a given input string. The unit must be enclosed in
    parentheses, and if multiple sets of parentheses are present, the
    last one will be returned. The test covers various scenarios,
    including happy paths, edge cases, and error cases.
    """
    result = get_unit(input_str)
    assert result == expected_output


@pytest.mark.unit
@pytest.mark.parametrize(
    ("bus_a", "bus_b", "expected"),
    [
        # Happy path tests
        pytest.param("AT0 1 H2", "AT0 1 H2", "local", id="same_region"),
        pytest.param("AT0 1 H2", "AT0 2 H2", "domestic", id="same_country"),
        pytest.param("AT0 1 H2", "DE0 1 H2", "foreign", id="different_country"),
        # Edge cases
        pytest.param("AT0 1 H2", "AT0 1 AC", "local", id="different_carrier"),
        pytest.param("", "AT0 1 H2", "", id="empty_bus_a"),
        pytest.param("AT0 1 H2", "", "", id="empty_bus_b"),
        pytest.param("", "", "", id="empty_both_buses"),
        # Error cases
        pytest.param("InvalidBus1", "AT0 1 H2", "", id="invalid_bus_a"),
        pytest.param("AT0 1 H2", "InvalidBus2", "", id="invalid_bus_b"),
        pytest.param("InvalidBus1", "InvalidBus2", "", id="invalid_both_buses"),
    ],
)
def test_get_trade_type(bus_a, bus_b, expected):
    """
    Test the get_trade_type function.

    This test verifies that the get_trade_type function correctly
    determines the trade type between two buses based on their region
    substrings. The function should return one of the following trade
    types: 'local', 'domestic', 'foreign', or an empty string.
    """
    result = get_trade_type(str(bus_a), str(bus_b))
    assert result == expected


@pytest.mark.unit
@pytest.mark.parametrize(
    ("trade_type", "buses", "expected"),
    [
        # Happy path tests
        pytest.param(
            "local", ("bus0", "bus1"), pd.Series([True, False, False]), id="local_trade"
        ),
        pytest.param(
            "domestic",
            ("bus0", "bus1"),
            pd.Series([False, True, False]),
            id="domestic_trade",
        ),
        pytest.param(
            "foreign",
            ("bus0", "bus1"),
            pd.Series([False, False, True]),
            id="foreign_trade",
        ),
        # Error cases
        pytest.param("invalid", ("bus0", "bus1"), ValueError, id="invalid_trade_type"),
        pytest.param("local", ("bus0", "bus1"), KeyError, id="missing_bus1_column"),
    ],
)
def test_trade_mask(trade_type, buses, expected, df_buses):
    """
    Test the trade_mask function.

    This test verifies that the trade_mask function correctly generates
    a mask for different trade types based on the provided component
    data frame and bus identifiers. The expected results are compared
    against the actual output from the function.
    """
    # sourcery skip: no-conditionals-in-tests
    if isinstance(expected, pd.Series):
        result = trade_mask(df_buses, str(trade_type), buses)
        pd.testing.assert_series_equal(result, expected)
    else:
        with pytest.raises(expected):
            trade_mask(df_buses.drop("bus1", axis=1), str(trade_type), buses)


# @pytest.mark.unit
# @pytest.mark.parametrize(
#     ("value", "level", "expected"),
#     [
#         # Happy path tests
#         pytest.param(
#             "C",
#             "idx1",
#             pd.DataFrame(
#                 {"col": [3, 5, 7]},
#                 index=pd.MultiIndex.from_product(
#                     [["C"], ["1", "2", "3"]], names=["idx1", "idx2"]
#                 ),
#             ),
#             id="replace_level_1",
#         ),
#         pytest.param(
#             "4",
#             "idx2",
#             pd.DataFrame(
#                 {"col": [3, 12]},
#                 index=pd.MultiIndex.from_product(
#                     [["A", "B"], ["4"]], names=["idx1", "idx2"]
#                 ),
#             ),
#             id="replace_level_2",
#         ),
#         # Error cases
#         pytest.param("C", "invalid", KeyError, id="invalid_level"),
#     ],
# )
# def test_replace_index_level_values(value, level, expected, df_multi_index):
#     """Test the replace_index_level_values function.
#
#     This test verifies that the replace_index_level_values function
#     correctly replaces index level values in a DataFrame with a
#     specified value. The expected results are compared against the
#     actual output from the function for both happy path and error
#     cases.
#     """
#     # sourcery skip: no-conditionals-in-tests
#     if isinstance(expected, pd.DataFrame):
#         result = replace_index_level_values(df_multi_index, str(value), str(level))
#         pd.testing.assert_frame_equal(result, expected)
#     else:
#         with pytest.raises(expected):
#             replace_index_level_values(df_multi_index, str(value), str(level))


@pytest.mark.unit
@pytest.mark.parametrize(
    ("exclude", "kwargs", "expected"),
    [
        # Happy path tests
        pytest.param(
            False,
            {"idx1": "A", "idx2": "1"},
            {("A", "1"): {"col": 0}},
            id="idx1_is_A_and_idx2_is_1",
        ),
        pytest.param(
            False,
            {"idx2": ["1", "3"]},
            {
                ("A", "1"): {"col": 0},
                ("A", "3"): {"col": 2},
                ("B", "1"): {"col": 3},
                ("B", "3"): {"col": 5},
            },
            id="idx2_is_1_or_3",
        ),
        pytest.param(
            True,
            {"idx1": "A", "idx2": ["1", "2"]},
            {
                ("A", "3"): {"col": 2},
                ("B", "1"): {"col": 3},
                ("B", "2"): {"col": 4},
                ("B", "3"): {"col": 5},
            },
            id="idx1_is_not_A_and_idx2_is_not_1_or_2",
        ),
        # edge cases
        pytest.param(
            True,
            {"idx1": ["A", "B"], "idx2": ["1", "2", "3"]},
            pd.DataFrame(
                columns=["col"],
                index=pd.MultiIndex.from_tuples([], names=["idx1", "idx2"]),
                dtype=np.int64,
            ),
            id="exclude_all",
        ),
        pytest.param(
            False,
            {"idx1": ["A", "B"], "idx2": ["1", "2", "3"]},
            {
                ("A", "1"): {"col": 0},
                ("A", "2"): {"col": 1},
                ("A", "3"): {"col": 2},
                ("B", "1"): {"col": 3},
                ("B", "2"): {"col": 4},
                ("B", "3"): {"col": 5},
            },
            id="include_all",
        ),
        # Error cases
        pytest.param(
            False,
            {"invalid": None},
            pd.errors.UndefinedVariableError,
            id="df_invalid_key",
        ),
        pytest.param(False, {}, ValueError, id="empty_query"),
    ],
)
def test_filter_by_data_frame(exclude, kwargs, expected, df_multi_index):
    """
    Test the filter_by function for DataFrame filtering.

    This test verifies that the filter_by function correctly filters
    a DataFrame based on specified index values and conditions.
    It checks both the expected output when filtering is successful
    and verifies that the appropriate exceptions are raised for
    invalid inputs.
    """
    # sourcery skip: no-conditionals-in-tests
    if not isinstance(expected, pd.DataFrame) and expected in (
        pd.errors.UndefinedVariableError,
        ValueError,
    ):
        with pytest.raises(expected):
            filter_by(df_multi_index, bool(exclude), **dict(kwargs))
    else:
        result = filter_by(df_multi_index, bool(exclude), **dict(kwargs))
        if isinstance(expected, dict):
            expected = pd.DataFrame.from_dict(expected, orient="index")
            expected.index.names = ["idx1", "idx2"]
        # skipping index type checking, because I couldn't construct the
        # expected index type.
        pd.testing.assert_frame_equal(result, expected, check_index_type=False)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("exclude", "kwargs", "expected"),
    [
        # Happy path tests
        pytest.param(
            False,
            {"idx1": "A", "idx2": "1"},
            {("A", "1"): 0},
            id="idx1_is_A_and_idx2_is_1",
        ),
        pytest.param(
            False,
            {"idx2": ["1", "3"]},
            {
                ("A", "1"): 0,
                ("A", "3"): 2,
                ("B", "1"): 3,
                ("B", "3"): 5,
            },
            id="idx2_is_1_or_3",
        ),
        pytest.param(
            True,
            {"idx1": "A", "idx2": ["1", "2"]},
            {("A", "3"): 2, ("B", "1"): 3, ("B", "2"): 4, ("B", "3"): 5},
            id="idx1_is_not_A_and_idx2_is_not_1_or_2",
        ),
        # edge cases
        pytest.param(
            True,
            {"idx1": ["A", "B"], "idx2": ["1", "2", "3"]},
            pd.Series(
                data=[],
                index=pd.MultiIndex.from_tuples([], names=["idx1", "idx2"]),
                dtype=np.int64,
                name="col",
            ),
            id="exclude_all",
        ),
        pytest.param(
            False,
            {"idx1": ["A", "B"], "idx2": ["1", "2", "3"]},
            {
                ("A", "1"): 0,
                ("A", "2"): 1,
                ("A", "3"): 2,
                ("B", "1"): 3,
                ("B", "2"): 4,
                ("B", "3"): 5,
            },
            id="include_all",
        ),
        # Error cases
        pytest.param(
            False,
            {"invalid": None},
            pd.errors.UndefinedVariableError,
            id="df_invalid_key",
        ),
    ],
)
def test_filter_by_series(exclude, kwargs, expected, ser_multi_index):
    """
    Test the filter_by function for Series filtering.

    This test verifies that the filter_by function correctly filters
    a Series based on specified index values and conditions.
    It checks both the expected output when filtering is successful
    and verifies that the appropriate exceptions are raised for
    invalid inputs.
    """
    # sourcery skip: no-conditionals-in-tests
    if (
        not isinstance(expected, pd.Series)
        and expected == pd.errors.UndefinedVariableError
    ):
        with pytest.raises(expected):
            filter_by(ser_multi_index, bool(exclude), **dict(kwargs))
    else:
        result = filter_by(ser_multi_index, bool(exclude), **dict(kwargs))
        if isinstance(expected, dict):
            expected = pd.Series(expected, name="col")
            expected.index.names = ["idx1", "idx2"]
        # couldn't construct the expected index type.
        pd.testing.assert_series_equal(result, expected, check_index_type=False)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("data_input", "snapshots", "data_expected"),
    [
        # Happy path tests
        pytest.param(
            {"metric": [8760, 17520]},
            pd.Index(["2021-01-01", "2022-01-01"]),
            [[1.0, 1.0], [2.0, 2.0]],
            id="single_column_df",
        ),
        pytest.param(
            [8760, 17520],
            pd.Index(["2021-01-01", "2022-01-01"]),
            [[1.0, 1.0], [2.0, 2.0]],
            id="series",
        ),
        # Edge cases
        pytest.param(
            {"metric": [0, 0]},
            pd.Index(["2021-01-01", "2022-01-01"]),
            [[0.0, 0.0], [0.0, 0.0]],
            id="zero_values",
        ),
        pytest.param(
            {"metric": [8760]},
            pd.Index(["2021-01-01"]),
            [[1.0]],
            id="single_value",
        ),
        pytest.param(
            {"metric": [8760, 17520]},
            pd.Index([]),
            [[], []],
            id="empty_snapshots",
        ),
        pytest.param({"metric": []}, pd.Index(["2021-01-01"]), [], id="empty_values"),
        pytest.param(
            {"A": [8760, 17520], "B": [8760, 17520]},
            pd.Index(["2021-01-01", "2022-01-01"]),
            NotImplementedError,
            id="multi_column_df",
        ),
    ],
)
def test_expand_to_time_series(data_input, snapshots, data_expected):
    """
    Test the expand_to_time_series function.

    This test verifies that the expand_to_time_series function correctly
    converts aggregated values into a time series format based on the
    provided snapshots. It checks both the expected output for valid
    inputs and the handling of edge cases.
    """
    # sourcery skip: no-conditionals-in-tests
    if data_expected is NotImplementedError:
        with pytest.raises(NotImplementedError):
            expand_to_time_series(pd.DataFrame(data_input), snapshots)
    else:
        if isinstance(data_input, dict):
            df_or_ser = pd.DataFrame(data_input)
        else:
            df_or_ser = pd.Series(data_input, name="metric")
        expected = pd.DataFrame(data=data_expected, columns=snapshots, dtype=float)
        result = expand_to_time_series(df_or_ser, snapshots)
        pd.testing.assert_frame_equal(result, expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("values", "ascending", "expected"),
    [
        # Happy path tests
        pytest.param(
            ("d", "c", "b"),
            False,
            pd.DataFrame({"A": ["a", "b", "c", "d"]}, index=[0, 1, 2, 3]),
            id="descending",
        ),
        pytest.param(
            ("d", "c", "b"),
            True,
            pd.DataFrame({"A": ["d", "c", "b", "a"]}, index=[3, 2, 1, 0]),
            id="ascending",
        ),
        pytest.param(
            ("c", "b"),
            False,
            pd.DataFrame({"A": ["a", "d", "b", "c"]}, index=[0, 3, 1, 2]),
            id="mixed_descending",
        ),
        pytest.param(
            ("c", "b"),
            True,
            pd.DataFrame({"A": ["c", "b", "a", "d"]}, index=[2, 1, 0, 3]),
            id="mixed_ascending",
        ),
        # Edge cases
        pytest.param(
            (),
            False,
            pd.DataFrame({"A": ["a", "b", "c", "d"]}, index=[0, 1, 2, 3]),
            id="empty_values",
        ),
        pytest.param(
            ("c",),
            False,
            pd.DataFrame({"A": ["a", "b", "d", "c"]}, index=[0, 1, 3, 2]),
            id="single_value_descending",
        ),
        pytest.param(
            ("c",),
            True,
            pd.DataFrame({"A": ["c", "a", "b", "d"]}, index=[2, 0, 1, 3]),
            id="single_value_ascending",
        ),
    ],
)
def test_custom_sort_happy_and_edge_cases(values, ascending, expected, df_sort):
    """
    Test the custom_sort function with various inputs.

    This test verifies that the custom_sort function correctly sorts
    a DataFrame based on specified values and order. It includes
    happy path tests, edge cases, and ensures that the function
    behaves as expected under different scenarios.
    """
    result = ESMChart.custom_sort(df_sort, "A", values, bool(ascending))
    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("values", "by", "expected"),
    [
        pytest.param(("a", "c", "b"), "invalid", KeyError, id="invalid_column"),
    ],
)
def test_custom_sort_error_cases(values, by, expected, df_sort):
    """
    Test the behavior of custom_sort with invalid input.

    This test checks that custom_sort raises the appropriate exception
    when an invalid column name is specified. It ensures that the
    function behaves correctly in error scenarios.
    """
    with pytest.raises(expected):
        ESMChart.custom_sort(df_sort, by, values, True)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("df", "year", "expected"),
    [
        # Happy path tests
        pytest.param(
            pd.DataFrame(columns=pd.to_datetime(["2015", "2016", "2017"])),
            2024,
            pd.DataFrame(columns=pd.to_datetime(["2024", "2024", "2024"])),
            id="year",
        ),
        # Edge cases
        pytest.param(
            pd.DataFrame(),
            2024,
            pd.DataFrame(),
            id="empty_dataframe",
        ),
        pytest.param(
            pd.DataFrame(columns=pd.to_datetime(["2015"])),
            2024,
            pd.DataFrame(columns=pd.to_datetime(["2024"])),
            id="single_column",
        ),
        pytest.param(
            pd.DataFrame(columns=["A", "B", "C"]),
            2024,
            pd.DataFrame(columns=["A", "B", "C"]),
            id="non_datetime_columns",
        ),
    ],
)
def test_fix_snapshots(df, year, expected):
    """
    Test the fix_snapshots function with various inputs.

    This test verifies that the fix_snapshots function correctly adjusts
    the timestamps in the DataFrame column labels to the specified year.
    It includes happy path tests, edge cases, and ensures that the
    function handles different scenarios appropriately.
    """
    result = ESMTimeSeriesChart.fix_snapshots(df, year)
    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("index", "names", "expected"),
    [
        # Happy path tests
        pytest.param(
            pd.MultiIndex.from_tuples(
                [("A", "B", "AT0 0 H2")], names=["lvl1", "lvl2", "lvl3"]
            ),
            ["lvl1", "lvl2", "loc", "carr"],
            pd.MultiIndex.from_tuples(
                [("A", "B", "AT0 0", "H2")], names=["lvl1", "lvl2", "loc", "carr"]
            ),
            id="single_entry",
        ),
        pytest.param(
            pd.MultiIndex.from_tuples(
                [("A", "B", "DE1 1 AC"), ("C", "D", "DE1 2 DC")],
                names=["lvl1", "lvl2", "lvl3"],
            ),
            ["lvl1", "lvl2", "loc", "carr"],
            pd.MultiIndex.from_tuples(
                [("A", "B", "DE1 1", "AC"), ("C", "D", "DE1 2", "DC")],
                names=["lvl1", "lvl2", "loc", "carr"],
            ),
            id="multiple_entries",
        ),
        # Edge cases
        pytest.param(
            pd.MultiIndex.from_tuples(
                [("A", "B", "CH4")], names=["lvl1", "lvl2", "lvl3"]
            ),
            ["lvl1", "lvl2", "loc", "carr"],
            pd.MultiIndex.from_tuples(
                [("A", "B", "", "CH4")],
                names=["lvl1", "lvl2", "loc", "carr"],
            ),
            id="no_location",
        ),
        pytest.param(
            pd.MultiIndex.from_tuples(
                [("A", "B", "AT0 0")], names=["lvl1", "lvl2", "lvl3"]
            ),
            ["lvl1", "lvl2", "loc", "carr"],
            pd.MultiIndex.from_tuples(
                [("A", "B", "AT0 0", "")], names=["lvl1", "lvl2", "loc", "carr"]
            ),
            id="no_carrier",
        ),
        pytest.param(
            pd.MultiIndex.from_tuples([("A", "B", "")], names=["lvl1", "lvl2", "lvl3"]),
            ["lvl1", "lvl2", "loc", "carr"],
            pd.MultiIndex.from_tuples(
                [("A", "B", "", "")], names=["lvl1", "lvl2", "loc", "carr"]
            ),
            id="emtpy_string",
        ),
        pytest.param(
            pd.MultiIndex.from_tuples([("AT0 0 AC",)], names=["lvl1"]),
            ["loc", "carr"],
            pd.MultiIndex.from_tuples([("AT0 0", "AC")], names=["loc", "carr"]),
            id="no_prefix",
        ),
        pytest.param(
            pd.MultiIndex.from_tuples(
                [("A", "AT0 0 H2", "B")], names=["lvl1", "lvl2", "lvl3"]
            ),
            ["lvl1", "loc", "carr", "lvl2"],
            pd.MultiIndex.from_tuples(
                [("A", "AT0 0 H2", "", "B")],
                names=[
                    "lvl1",
                    "loc",
                    "carr",
                    "lvl2",
                ],
            ),
            id="with_suffix",
        ),
    ],
)
def test_split_location_carrier(index, names, expected):
    """
    Test the split_location_carrier function with various inputs.

    This test verifies that the split_location_carrier function
    correctly splits the location and carrier from the innermost
    level of a MultiIndex. It includes happy path tests, edge
    cases, and ensures that the function handles different
    scenarios appropriately.
    """
    result = split_location_carrier(index, names)
    pd.testing.assert_index_equal(result, expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("df", "to_unit", "expected"),
    [
        # Happy path tests
        pytest.param(
            pd.DataFrame(
                {
                    "base (Wh)": [-1, 0, 1],
                    "kilo (kWh)": [-1, 0, 1],
                    "mega (MWh)": [-1, 0, 1],
                    "giga (GWh)": [-1, 0, 1],
                    "terra (TWh)": [-1, 0, 1],
                },
                dtype=float,
            ),
            "TWh",
            pd.DataFrame(
                {
                    "base (TWh)": [-1e-12, 0e-12, 1e-12],
                    "kilo (TWh)": [-1e-9, 0e-9, 1e-9],
                    "mega (TWh)": [-1e-6, 0e-6, 1e-6],
                    "giga (TWh)": [-1e-3, 0e-3, 1e-3],
                    "terra (TWh)": [-1e0, 0e0, 1e0],
                },
                dtype=float,
            ),
            id="energy_up",
        ),
        pytest.param(
            pd.DataFrame(
                {
                    "base (W)": [-1, 0, 1],
                    "kilo (kW)": [-1, 0, 1],
                    "mega (MW)": [-1, 0, 1],
                    "giga (GW)": [-1, 0, 1],
                    "terra (TW)": [-1, 0, 1],
                },
                dtype=float,
            ),
            "W",
            pd.DataFrame(
                {
                    "base (W)": [-1, 0, 1],
                    "kilo (W)": [-1e3, 0e3, 1e3],
                    "mega (W)": [-1e6, 0e6, 1e6],
                    "giga (W)": [-1e9, 0e9, 1e9],
                    "terra (W)": [-1e12, 0e12, 1e12],
                },
                dtype=float,
            ),
            id="power_down",
        ),
        # Edge cases
        pytest.param(
            pd.DataFrame({"base (W)": [None, np.nan]}),
            "MW",
            pd.DataFrame({"base (MW)": [None, np.nan]}),
            id="nan_values",
        ),
        # Error cases
        pytest.param(
            pd.DataFrame({"base (W)": [1]}), "MWh", ValueError, id="power_to_energy"
        ),
        pytest.param(
            pd.DataFrame({"base (kWh)": [1]}), "TW", ValueError, id="energy_to_power"
        ),
        pytest.param(
            pd.DataFrame({"base (kWh)": [1], "base (kW)": [1]}),
            "TWh",
            ValueError,
            id="mixed_power_and_energy",
        ),
        pytest.param(
            pd.DataFrame({"base (kW)": [1]}),
            "",
            KeyError,
            id="empty_to_unit",
        ),
        pytest.param(
            pd.DataFrame({"base (Wh)": [1]}), "invalid", KeyError, id="invalid_unit"
        ),
    ],
)
def test_scale(df, to_unit, expected):
    """
    Test the behavior of the scale function with various inputs.

    This test verifies that the scale function correctly scales the
    metric columns of a DataFrame to the specified unit. It includes
    happy path tests, edge cases, and ensures that the function
    handles different scenarios appropriately.
    """
    if not isinstance(expected, pd.DataFrame):
        with pytest.raises(expected):
            scale(df, to_unit)
    else:
        result = scale(df, to_unit)
        pd.testing.assert_frame_equal(result, expected)
        assert result.attrs.get("unit", "") == to_unit


@pytest.mark.unit
@pytest.mark.parametrize(
    ("mapper", "level", "agg", "expected"),
    [
        # Happy path tests
        pytest.param(
            {"A": "C", "B": "C"},
            "idx1",
            "sum",
            pd.DataFrame({"col": {("C", "1"): 3.0, ("C", "2"): 5.0, ("C", "3"): 7.0}}),
            id="rename_agg_sum",
        ),
        # Edge cases
        pytest.param({}, "idx1", "mean", None, id="empty_mapper"),
        pytest.param({"stinky": "fish"}, "idx1", "mean", None, id="fishy_mapper"),
        pytest.param({"1": "1", "2": "2"}, "idx2", "max", None, id="identity_mapper"),
        pytest.param(
            {"A": "C", "B": "C"},
            "idx1",
            ["min", "max"],
            pd.DataFrame(
                {
                    ("col", "min"): {("C", "1"): 0, ("C", "2"): 1, ("C", "3"): 2},
                    ("col", "max"): {("C", "1"): 3, ("C", "2"): 4, ("C", "3"): 5},
                }
            ),
            id="multiple_aggregations",
        ),
        # Error cases
        pytest.param({"A": "C"}, "invalid", "sum", KeyError, id="invalid_level"),
        pytest.param({"A": "C"}, "idx1", "invalid", AttributeError, id="invalid_agg"),
    ],
)
def test_apply_mapping(mapper, level, agg, expected, df_multi_index):
    """Test the apply_mapping functions."""
    # sourcery skip: no-conditionals-in-tests
    expected_errors = (KeyError, AttributeError)
    if not isinstance(expected, pd.DataFrame) and expected in expected_errors:
        with pytest.raises(expected):
            rename_aggregate(df_multi_index, mapper, level, agg)
    else:
        # simplify parameters
        expected = df_multi_index if expected is None else expected
        expected.index.names = ["idx1", "idx2"]
        result = rename_aggregate(df_multi_index, mapper, level, agg)
        # the apply_mapping function does not preserve data types!
        pd.testing.assert_frame_equal(result, expected, check_dtype=False)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("df", "limit", "drop", "expected"),
    [
        # Happy path tests
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-0.1, -0.5, -1.0]}),
            0.5,
            True,
            pd.DataFrame({"A": [0.5, 1.0], "B": [-0.5, -1.0]}, index=[1, 2]),
            id="replace_drop",
        ),
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-0.1, -0.5, -1.0]}),
            0.5,
            False,
            pd.DataFrame({"A": [np.nan, 0.5, 1.0], "B": [np.nan, -0.5, -1.0]}),
            id="replace_keep",
        ),
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            0.5,
            False,
            pd.DataFrame({"A": [np.nan, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            id="replace_drop_no_rows",
        ),
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            -0.6,
            True,
            pd.DataFrame({"A": [np.nan, 1.0], "B": [-1.0, -1.0]}, index=[0, 2]),
            id="negative_limit",
        ),
        # Edge cases
        pytest.param(
            pd.DataFrame({"A": [0.0, 0.0, 0.0], "B": [0.0, 0.0, 0.0]}),
            0.1,
            True,
            pd.DataFrame({"A": {}, "B": {}}),
            id="replace_drop_all",
        ),
        pytest.param(
            pd.DataFrame({"A": [0.0, 0.0, 0.0], "B": [0.0, 0.0, 0.0]}),
            0.1,
            False,
            pd.DataFrame(
                {"A": [np.nan, np.nan, np.nan], "B": [np.nan, np.nan, np.nan]}
            ),
            id="replace_keep_all",
        ),
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            0.1,
            True,
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            id="replace_none",
        ),
        pytest.param(pd.DataFrame(), 0.0, True, pd.DataFrame(), id="empty_input"),
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            np.inf,
            True,
            pd.DataFrame({"A": {}, "B": {}}),
            id="infinity_limit",
        ),
        pytest.param(
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            np.nan,
            True,
            pd.DataFrame({"A": [0.1, 0.5, 1.0], "B": [-1.0, -0.5, -1.0]}),
            id="nan_limit",
        ),
    ],
)
def test_apply_cutoff(df, limit, drop, expected):
    """Tests the apply_cutoff function."""
    result = apply_cutoff(df, limit, drop)
    pd.testing.assert_frame_equal(result, expected, check_index_type=False)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("df", "keep_regions", "nice_names", "expected"),
    [
        # Happy path tests
        pytest.param(
            df_locations(),
            ("AT",),
            True,
            pd.DataFrame(
                {
                    "a": {
                        ("", "A"): 6,
                        ("AT - Lower Austria", "A"): 1,
                        ("AT - Vienna", "A"): 2,
                        ("AT - Vienna", "B"): 3,
                        ("Austria", "A"): 3,
                        ("Austria", "B"): 3,
                        ("Europe", "A"): 12,
                        ("Europe", "B"): 3,
                        ("France", "A"): 0,
                        ("Germany", "A"): 9,
                    }
                },
            ),
            id="AT_nice_names",
        ),
        pytest.param(
            df_locations(),
            ("AT",),
            False,
            pd.DataFrame(
                {
                    "a": {
                        ("", "A"): 6,
                        ("AT", "A"): 3,
                        ("AT", "B"): 3,
                        ("AT0 1", "A"): 1,
                        ("AT0 2", "A"): 2,
                        ("AT0 2", "B"): 3,
                        ("DE", "A"): 9,
                        ("EU", "A"): 12,
                        ("EU", "B"): 3,
                        ("FR", "A"): 0,
                    }
                }
            ),
            id="AT_country_codes",
        ),
        pytest.param(
            df_locations(),
            ("DE",),
            True,
            pd.DataFrame(
                {
                    "a": {
                        ("", "A"): 6,
                        ("Austria", "A"): 3,
                        ("Austria", "B"): 3,
                        ("DE - Baden Wurttemberg", "A"): 4,
                        ("DE - Midwest", "A"): 5,
                        ("Europe", "A"): 12,
                        ("Europe", "B"): 3,
                        ("France", "A"): 0,
                        ("Germany", "A"): 9,
                    }
                }
            ),
            id="DE_nice_names",
        ),
        # Edge cases
        pytest.param(
            df_locations().iloc[[0], :],
            ("AT",),
            True,
            pd.DataFrame(
                {
                    "a": {
                        ("Europe", "A"): 0,
                        ("France", "A"): 0,
                    }
                }
            ),
            id="single_entry",
        ),
        pytest.param(
            df_locations(),
            (),
            True,
            pd.DataFrame(
                {
                    "a": {
                        ("", "A"): 6,
                        ("Austria", "A"): 3,
                        ("Austria", "B"): 3,
                        ("Europe", "A"): 12,
                        ("Europe", "B"): 3,
                        ("France", "A"): 0,
                        ("Germany", "A"): 9,
                    }
                }
            ),
            id="no_keep_regions",
        ),
        # Error cases
        pytest.param(pd.DataFrame(), ("AT",), True, KeyError, id="empty_data"),
    ],
)
def test_aggregate_locations(df, keep_regions, nice_names, expected):
    """Test the aggregation logic for locations aka nodes."""
    # sourcery skip: no-conditionals-in-tests
    if not isinstance(expected, pd.DataFrame):
        with pytest.raises(expected):
            aggregate_locations(df, keep_regions, nice_names)
    else:
        result = aggregate_locations(df, keep_regions, nice_names).sort_index()
        expected.index.names = ["location", "carrier"]
        pd.testing.assert_frame_equal(result, expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("df", "expected"),
    [
        pytest.param(df_metric(), None, id="valid_input"),
        pytest.param(pd.Series(data=[1, 2, 3]), AssertionError, id="input_series"),
        pytest.param(
            df_metric().rename_axis(columns="invalid"),
            AssertionError,
            id="missing_column_label",
        ),
        pytest.param(
            df_metric().rename_axis(index=["invalid"] + DataModel.IDX_NAMES),
            AssertionError,
            id="missing_index_label",
        ),
        pytest.param(df_metric("", "(unit)"), AssertionError, id="missing_attrs_name"),
        pytest.param(df_metric("name", ""), AssertionError, id="missing_attrs_unit"),
        pytest.param(
            df_metric("name", "unit"), AssertionError, id="missing_unit_braces"
        ),
        pytest.param(
            pd.concat([df_metric(), df_metric()], axis=1),
            AssertionError,
            id="multiple_metrics",
        ),
    ],
)
def test_verify_metric_format(df, expected):
    """Test the verify_metric_format function."""
    if expected is AssertionError:
        with pytest.raises(expected):
            verify_metric_format(df)
    else:
        verify_metric_format(df)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("x", "expected"),
    [
        pytest.param(0.0, "0.0"),
        pytest.param(0.5, "0.5"),
        pytest.param(1.0, "1.0"),
        pytest.param(1.5, "1.5"),
        pytest.param(10.0, "10"),
        pytest.param(10.05, "10"),
        pytest.param(10.09, "10"),
        pytest.param(1.4499999, "1.4"),
        pytest.param(1.499999, "1.5"),
        pytest.param(0.15, "0.2"),
        pytest.param(9.499999, "9.5"),
        pytest.param(9.5, "10"),
        pytest.param(10.449999, "10"),
        pytest.param(10.49999, "11"),
        pytest.param(10.5, "11"),
        pytest.param(99.49, "99"),
        pytest.param(99.5, "100"),
    ],
)
def test_prettify_numer(x, expected):
    result = prettify_number(x)
    assert result == expected
