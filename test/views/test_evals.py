"""Tests to compare the Toolbox with esmtools evaluation results."""

import json
from itertools import product
from pathlib import Path

import pandas as pd
import pytest

RESULTS_PATH = Path(r"C:\Storage\pypsa-eur-sec\results")
DIR_JSON_OLD = (
    # RESULTS_PATH / "esm_run" / "evaluation_tests-comparison-base" / "JSONs_LC"
    RESULTS_PATH / "esm_run" / "evaluation" / "JSONs_LC"
)
DIR_JSON_NEW = RESULTS_PATH / "esm_run" / "evaluation" / "JSON"


class Constants:
    """Container class to collect constants used in tests."""

    PRECISION_TWH: int = 9
    TOLERANCE_TWH_TIME_SERIES: int = 10  # TWh
    LOCATIONS: tuple = (
        "AL",
        "AT",
        "AT0 0",
        "AT0 1",
        "AT0 2",
        "AT0 3",
        "AT0 4",
        "AT0 5",
        "AT0 6",
        "AT0 7",
        "AT0 8",
        "AT0 9",
        "BA",
        "BE",
        "BG",
        "CH",
        "CZ",
        "DE",
        "DK",
        "EE",
        "ES",
        "EU",
        "FI",
        "FR",
        "GB",
        "GR",
        "HR",
        "HU",
        "IE",
        "IT",
        "LT",
        "LU",
        "LV",
        "ME",
        "MK",
        "NL",
        "NO",
        "PL",
        "PT",
        "RO",
        "RS",
        "SE",
        "SI",
        "SK",
    )
    YEARS: tuple = (
        "2020",
        "2030",
        "2040",
        "2050",
    )

    FILE_NAMES_TIME_SERIES = [
        "elec_prod_dem_time_{year}_{location}",
        "heat_prod_dem_time_{year}_{location}",
        "hydrogen_prod_dem_time_{year}_{location}",
        "methane_prod_dem_time_{year}_{location}",
        # missing in old evaluations:
        "hydro_phs_storage_operation_{year}_{location}",
        "hydro_res_storage_operation_{year}_{location}",
    ]
    FILE_NAMES_BARCHARTS = [
        "biomass_balance_{location}",
        "co2_emissions_{location}",
        "eKapas_{location}",
        "elec_balance_{location}",
        "elec_demand_{location}",
        "final_demand_{location}",
        "fKapas_{location}",
        "heat_balance_{location}",
        "heatKapas_{location}",
        "hydrogen_balance_{location}",
        "ind_sec_tot_{location}",
        "methane_balance_{location}",
        "trans_sec_tot_{location}",
        "building_sec_tot_{location}",
        "building_sec_heat_{location}",
        "StoreVol_{location}",
        "elec_storage_volumes_{location}",
        # missing in old Toolbox implementation
        "elec_production_{location}",
        "curtailment_{location}",
        "capacity_factor_{location}",
        "heat_mix_{location}",
    ]
    FILE_NAMES_GROUPED_BARCHARTS = [
        "building_sec_heat_{location}",
        "ind_sec_dem_{location}",
        "trans_sec_dem_{location}",
        "final_demand_sec_{location}",
    ]

    # from JSONs_LC directory. However, there are more HTMLs
    all_jsons = [
        # maps
        "HydrogenTransmission_{year}",
        "MethaneTransmission_{year}",
        "PowerTransmission_{year}",
        # sankey
        "Sankey_detail_transformation_{year}",
        "Sankey_{year}",
        # barcharts
        "StoreVol",
        "biomass_balance",
        "building_sec_tot",
        "co2_emissions",
        "eKapas",
        "elec_balance",
        "elec_demand",
        "elec_storage_volumes",
        "fKapas",
        "final_demand",
        "heatKapas",
        "heat_balance",
        "hydrogen_balance",
        "ind_sec_tot",
        "methane_balance",
        # grouped barcharts
        "trans_sec_dem",
        "final_demand_sec",
        "ind_sec_dem",
        "building_sec_heat",
        # scatter plots
        "elec_prod_dem_time_{year}",
        "hydrogen_prod_dem_time_{year}",
        "heat_prod_dem_time_{year}",
        "methane_prod_dem_time_{year}",
    ]


def _barchart_file_names():
    """
    Construct all file names that have only location encoded.

    Some file names are marked as xfail. These tests return green
    if they fail and red if they succeed. This is because some
    evaluations have bugs and the tests should pass nevertheless.
    """
    for location, file_name_template in product(
        Constants.LOCATIONS, Constants.FILE_NAMES_BARCHARTS
    ):
        file_name = file_name_template.format(location=location)
        reason = False
        strict = True
        if file_name in (
            "heatKapas_BA",
            "heatKapas_DK",
            "heatKapas_EU",
            "heatKapas_FI",
            "heatKapas_GB",
            "heatKapas_IE",
            "heatKapas_LT",
            "heatKapas_NO",
            "heatKapas_RS",
        ):
            reason = (
                "Heat capacities in the old Toolbox have a bug. "
                "p_nom_opt values for Link ports > 1 "
                "are not multiplied with their respective port "
                "efficiency."
            )
        elif file_name in ("elec_balance_AL", "elec_demand_AL"):
            reason = (
                "The AC split has a Bug in the old Toolbox "
                "implementation. The carrier 'domestic homes and "
                "trade' is counted, although positive "
                "demand should be clipped to zero."
            )
        elif file_name in (
            "final_demand_EU",
            "building_sec_tot_EU",
            "building_sec_heat_EU",
        ):
            reason = (
                "Final Demand in old Toolbox shows unreasonably "
                "low values for EU region."
            )
        elif file_name in ("biomass_balance_BA", "biomass_balance_EU"):
            reason = (
                "Electricity production from biomass uses a wrong "
                "efficiency for AC parts. See the bugfix .ipynb "
                "for more information. This bug exists for all "
                "regions, but only BA and EU have amounts bigger "
                "than the configured cutoff. As a result, the bug "
                "only appears for those regions."
            )
        elif file_name.startswith("heat_balance_"):
            reason = (
                "The Toolbox function 'get_district_heat_balances_data' "
                "returns imprecise results that make the tests fail. "
                "See the bugfix_heat-balances.ipynb for additional info."
            )
            strict = False
        marks = pytest.mark.xfail(reason=reason, strict=strict) if reason else ()
        yield pytest.param(file_name, marks=marks)


def _grouped_barchart_file_names():
    """Construct all time series file names."""
    for location, file_name_template in product(
        Constants.LOCATIONS, Constants.FILE_NAMES_GROUPED_BARCHARTS
    ):
        file_name = file_name_template.format(location=location)
        reason = False
        if file_name == "final_demand_sec_EU":
            reason = (
                "Final Demand for Europe has unreasonably low values in old"
                "Toolbox implementation."
            )
        elif file_name == "final_demand_sec_AL":
            reason = (
                "Final Demand for Albania is wrong, because small "
                "amounts of electricity demand from carrier 'domestic_homes_and_trade'"
                "are counted as production in the Toolbox implementation."
                "See bugfix_missing-clipping-in-ac-split.ipynb for more info."
            )
        marks = pytest.mark.xfail(reason=reason, strict=True) if reason else ()
        yield pytest.param(file_name, marks=marks)


def _time_series_file_names():
    """Construct all time series file names."""
    for location, year, file_name in product(
        Constants.LOCATIONS, Constants.YEARS, Constants.FILE_NAMES_TIME_SERIES
    ):
        file_name = file_name.format(location=location, year=year)
        reason = False
        strict = True

        if file_name == "elec_prod_dem_time_2040_LV":
            reason = (
                "The Toolbox wrongly counts the 'MISS.*' Link at port 1 as Base Load."
            )
        elif file_name.endswith(
            ("AT", "DE", "DK", "ES", "GB", "IT")
        ) and file_name.startswith("elec_prod_dem_time_"):
            if file_name in (
                "elec_prod_dem_time_2040_DK",
                "elec_prod_dem_time_2050_DK",
            ):
                strict = False  # because DK 2040 and 2050 do not show differences
            reason = (
                "The Toolbox aggregates nodes to countries before "
                "clipping, while esmtools evaluations per node and "
                "clips per node before aggregation. This leads to "
                "differences in the base load."
            )
        elif file_name.startswith("methane_prod_dem_time_") and file_name.endswith(
            ("BG", "ES", "GR", "IT", "PT", "EU")
        ):
            reason = (
                "The Toolbox either misses the technologies ["
                "'residential urban decentral gas boiler', "
                "'services urban decentral gas boiler'"
                "] in the time series, or wrongly adds the technologies "
                "in the balances. I decided to include them in both, the "
                "time series and the balances. Hence, the respective "
                "tests are marked as xfail. "
                "See bugfix_methane_time_series.ipynb for additional info."
            )
        elif file_name in (
            "methane_prod_dem_time_2040_AT",
            "methane_prod_dem_time_2040_DE",
        ):
            reason = (
                "Small differences in the Toolbox implementation for "
                "CNG technologies cause differences. See "
                "bugfix_methane-transport-AT-2030.ipynb for additional "
                "info."
            )

        marks = pytest.mark.xfail(reason=reason, strict=strict) if reason else ()
        yield pytest.param(file_name, marks=marks)


def _fetch_json_data(file_name):
    fp_old = Path(DIR_JSON_OLD) / f"{file_name}.json"
    if not fp_old.is_file():
        pytest.skip(f"File not found: {fp_old}")

    fp_new = Path(DIR_JSON_NEW) / f"{file_name}.json"
    assert fp_new.is_file()

    with fp_old.open("r", encoding="utf-8") as fh:
        json_old = json.load(fh)
    with fp_new.open("r", encoding="utf-8") as fh:
        json_new = json.load(fh)

    # old grouped bar charts have no names
    old = sorted(json_old["data"], key=lambda x: x.get("name", ""))
    new = sorted(json_new["data"], key=lambda x: x["name"])

    return old, new


def _as_data_frame(data):
    # cannot compare results by xaxis because xaxis enumeration is different =(
    # idx = [(c[0], d.get("xaxis").replace("x11", "x2")) for d in data for c in d["customdata"]]
    idx = [c[0] for d in data for c in d["customdata"]]
    y = [y for d in data for y in d.get("y")]
    x = [x for d in data for x in d.get("x")]
    df = pd.DataFrame(
        index=idx,
        data=zip(y, x, strict=True),
        columns=["data", "year"],
    )
    df.index.name = "category"
    df = df.set_index("year", append=True).sort_index()
    return df


def _as_data_frame_ts(data):
    index = [d["name"] for d in data]
    # assuming all snapshots are the same
    columns = data[0]["x"] if data else []
    values = [d["y"] for d in data]
    new = pd.DataFrame(index=index, columns=columns, data=values)

    return new


@pytest.mark.migration
@pytest.mark.parametrize("file_name", _barchart_file_names())
def test_compare_barchart_data(file_name: str):
    """Compare old and new JSON files for equal total sums."""
    old_data, new_data = _fetch_json_data(file_name)
    drop_traces = ("Upper Sum", "Lower Sum", "Sum", "No data")

    old_clean = [d for d in old_data if d["name"] not in drop_traces]
    old = _as_data_frame(old_clean)

    new_clean = [d for d in new_data if d["name"] not in drop_traces]
    new = _as_data_frame(new_clean)

    # Remove empty entries or data frames cannot be compared
    old = old[old.abs() > 0].dropna()
    new = new[new.abs() > 0].dropna()

    # The old Toolbox has Electricity (OCGT) instead of the usual and
    # more consistently used Electricity OCGT in the methane balance
    # mappings. I do not want to replicate inconsistencies, so I add
    # the following hotfix for this inconsistency so that the test
    # may pass.
    old = old.rename(
        {"Electricity (OCGT)": "Electricity OCGT"}, level="category"
    ).sort_index()

    # combine duplicated index entries or data frames cannot be compared
    old = old.groupby(old.index.names).sum()
    new = new.groupby(new.index.names).sum()

    # new.index.difference(old.index)
    # old.index.difference(new.index)

    absolute_tolerance = 1e-2  # this is below visibility in the HTML

    # Transpose to improve readability of AssertionError
    pd.testing.assert_frame_equal(old.T, new.T, atol=absolute_tolerance)


@pytest.mark.migration
@pytest.mark.parametrize("file_name", _grouped_barchart_file_names())
def test_compare_grouped_barchart_data(file_name):
    """Compare old and new JSON files."""
    old_data, new_data = _fetch_json_data(file_name)

    old_data_clean = [
        d for d in old_data if d.get("mode", "") not in ("text", "markers+text")
    ]
    old = _as_data_frame(old_data_clean)

    new_data_clean = [d for d in new_data if "Sum" not in d.get("name", "")]
    new = _as_data_frame(new_data_clean)

    # drop emtpy entries or data frames cannot be compared
    old = old[old > 0].dropna()
    new = new[new > 0].dropna()

    # combine duplicated index entries or data frames cannot be compared
    old = old.groupby(old.index.names).sum()
    new = new.groupby(new.index.names).sum()

    # new.index.difference(old.index)
    # old.index.difference(new.index)

    absolute_tolerance = 1e-2
    pd.testing.assert_frame_equal(old.T, new.T, atol=absolute_tolerance)


@pytest.mark.migration
@pytest.mark.parametrize("file_name", _time_series_file_names())
def tests_compare_timeseries_data(file_name):
    """Compare time series charts."""
    old_data, new_data = _fetch_json_data(file_name)

    old_data_clean = [d for d in old_data if d["name"] != "no data"]
    old = _as_data_frame_ts(old_data_clean)
    new = _as_data_frame_ts(new_data)

    # old data has duplicated snapshot values. we only use snapshots
    # from the new data set without duplicated snapshot values to reduce
    # test performance.
    old = old.filter(new.columns)

    # esmtools differentiates between import from outside the EU and
    # net imports. We squash them again for a comparison with the
    # Toolbox.
    new = new.rename({"Import Global": "Net Import"}).groupby(level=0).sum()

    # old data is rounded to second decimal place. new data is not,
    # and I do not like to introduce rounding biases in the evaluation.
    # We round new data series here for the comparison with old results.
    new = new.round(2)

    # We also need to drop all zero time series after rounding,
    # because that is what the old Toolbox does before JSON export.
    # (we also need to do the same for the old values, because the
    # Toolbox does not do this consistently for all time series, e.g.
    # AC and gas use different rounding)
    new = new.loc[new.abs().gt(0.0).any(axis=1)]
    old = old.loc[old.abs().gt(0.0).any(axis=1)]

    # We also need to fill NaNs with zeros (0.0) for a comparison
    new = new.fillna(0.0).sort_index()

    # The Toolbox rounds before clipping, while esmtools only clips
    # values using the cut off. Unfortunately, some small differences
    # cannot be avoided for very small values. Consider 0.05 being
    # rounded to 0.1, but cut off by a threshold of 0.1. There are a
    # edge cases I cannot avoid easily.
    if (
        file_name.startswith("hydrogen_prod_dem_time")
        and "Miscellaneous" in new.index
        and all(new.loc["Miscellaneous"].abs() < 0.05)  # MWh
        and "Miscellaneous" not in old.index
    ):
        new = new.drop("Miscellaneous")

    # new.index.difference(old.index)
    # old.index.difference(new.index)

    absolute_tolerance = 1e-1  # MWh
    pd.testing.assert_frame_equal(old.T, new.T, atol=absolute_tolerance)


@pytest.mark.migration
@pytest.mark.parametrize("file_name", _barchart_file_names())
def test_compare_barchart_layouts(file_name):
    """Compare old and new JSON files for equal total sums."""
    old_data, new_data = _fetch_json_data(file_name)
