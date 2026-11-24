# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the diagonal LT fallback of the stochastic grid analysis
(`stochastic_grid_scenarios.diagonal_lt_fallback`): a diagonal ST cell
(portfolio g on grid g) whose solve failed with the current `solve_fallback`
ladder reads its LT network instead. That is valid only if LT and ST relax the
dispatch alike (load shedding, CO2-cap slack).
"""

import sys

import pytest

sys.path.insert(0, ".")

from scripts._helpers import (  # noqa: E402
    check_lt_st_relaxations,
    lt_replaces_st,
    write_st_failure_marker,
)

SHEDDING = {"enable": True, "carriers": {"low voltage": 4000}}
SLACK_ON = {"enable": True, "penalty": 1e4}
SLACK_OFF = {"enable": False, "penalty": 1e4}
LADDER = [{"BarHomogeneous": 1}]


def _options(shedding: dict, slack: dict) -> dict:
    return {"load_shedding": shedding, "co2_cap_slack": slack}


@pytest.mark.parametrize(
    "lt, st, allow_shedding, allow_slack, error",
    [
        (_options(SHEDDING, SLACK_ON), _options(SHEDDING, SLACK_ON), False, False, None),
        (_options(SHEDDING, SLACK_OFF), _options(SHEDDING, {"enable": False}), False, False, None),
        (_options({"enable": False}, SLACK_OFF), _options(SHEDDING, SLACK_OFF), False, False, "load shedding"),
        (_options({"enable": False}, SLACK_OFF), _options(SHEDDING, SLACK_OFF), True, False, None),
        (_options(SHEDDING, SLACK_OFF), _options(SHEDDING, SLACK_ON), False, False, "CO2-cap slack"),
        (_options(SHEDDING, SLACK_OFF), _options(SHEDDING, SLACK_ON), False, True, None),
        (_options(SHEDDING, SLACK_ON), _options(SHEDDING, {**SLACK_ON, "penalty": 1e3}), False, False, "CO2-cap slack"),
    ],
)
def test_check_lt_st_relaxations(lt, st, allow_shedding, allow_slack, error):
    if error is None:
        check_lt_st_relaxations(lt, st, allow_shedding, allow_slack)
    else:
        with pytest.raises(ValueError, match=error):
            check_lt_st_relaxations(lt, st, allow_shedding, allow_slack)


@pytest.mark.parametrize(
    "st_exists, marker_ladder, replaced",
    [
        (False, LADDER, True),  # failed with every configured tier
        (False, None, False),  # never failed (or killed before writing a marker)
        (False, [], False),  # failed with another ladder: solve again
        (True, LADDER, False),  # a later solve succeeded
    ],
)
def test_lt_replaces_st(tmp_path, st_exists, marker_ladder, replaced):
    st = tmp_path / "portfolio-a_on-a_st.nc"
    if st_exists:
        st.touch()
    if marker_ladder is not None:
        write_st_failure_marker(
            str(st), status="ok", condition="suboptimal", relative_gap=float("inf"),
            solve_fallback=marker_ladder,
        )
    assert lt_replaces_st(str(st), LADDER) == replaced
