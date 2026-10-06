# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/pypsa-de/build_grid_topology.py, in
particular the scenario-specific scaling of line num_parallel/reactance/
resistance and the dropping of zeroed corridors in the deterministic/eev
variants, versus the canvas-impedance behaviour of the native stochastic
variant.

Lines carry a standard ``type``, so ``n.calculate_dependent_values()``
recomputes ``x``/``r`` from ``num_parallel`` - exactly as the solver does -
which is why the tests assert ``x``/``r`` *after* that recompute.
"""

import importlib.util
import pathlib
import sys

import numpy as np
import pandas as pd
import pypsa
import pytest

sys.path.insert(0, ".")
sys.path.append("./scripts")

_spec = importlib.util.spec_from_file_location(
    "build_grid_topology", "scripts/pypsa-de/build_grid_topology.py"
)
bgt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bgt)


def _canvas():
    """
    Minimal fully-built canvas: a 3-bus AC triangle (typed lines) plus one DC
    link. x/r are recomputed from the type and num_parallel.
    """
    n = pypsa.Network()
    n.add(
        "LineType", "t", r_per_length=0.01, x_per_length=0.1, c_per_length=10.0,
        f_nom=50.0, i_nom=1.0,
    )
    for b in ["A", "B", "C"]:
        n.add("Bus", b, v_nom=380.0, carrier="AC")
    kw = dict(type="t", carrier="AC")
    n.add("Line", "L1", bus0="A", bus1="B", length=100.0, num_parallel=2.0, s_nom=200.0, **kw)
    n.add("Line", "L2", bus0="B", bus1="C", length=100.0, num_parallel=4.0, s_nom=400.0, **kw)
    n.add("Line", "L3", bus0="A", bus1="C", length=100.0, num_parallel=8.0, s_nom=800.0, **kw)
    n.add("Link", "K1", bus0="A", bus1="C", p_nom=500.0, carrier="DC")
    n.calculate_dependent_values()
    return n


def _write_csvs(tmp_path, lines, links):
    tmp_path = pathlib.Path(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    lp = tmp_path / "lines.csv"
    kp = tmp_path / "links.csv"
    pd.Series(lines, name="s_nom").rename_axis("name").to_csv(lp)
    pd.Series(links, name="p_nom").rename_axis("name").to_csv(kp)
    return {"lines": str(lp), "links": str(kp)}


@pytest.mark.parametrize(
    "s_nom_new, frac, dropped",
    [
        ({"L1": 100.0, "L2": 400.0}, {"L1": 0.5}, []),  # L1 halved, L2 unchanged
        ({"L1": 50.0}, {"L1": 0.25}, []),  # quartered
        ({"L1": 200.0}, {"L1": 1.0}, []),  # unchanged
        ({"L1": 0.0, "L2": 200.0}, {"L2": 0.5}, ["L1"]),  # zeroed -> dropped
    ],
)
def test_set_line_capacity(s_nom_new, frac, dropped):
    n = _canvas()
    x0 = n.lines["x"].copy()
    r0 = n.lines["r"].copy()
    np0 = n.lines["num_parallel"].copy()

    bgt._set_line_capacity(n, pd.Series(s_nom_new))
    n.calculate_dependent_values()  # as the solver would

    for name in dropped:
        assert name not in n.lines.index
    for name, f in frac.items():
        assert n.lines.at[name, "num_parallel"] == pytest.approx(np0[name] * f)
        assert n.lines.at[name, "x"] == pytest.approx(x0[name] / f)  # impedance up as capacity down
        assert n.lines.at[name, "r"] == pytest.approx(r0[name] / f)
        assert n.lines.at[name, "s_nom"] == pytest.approx(s_nom_new[name])
        assert not n.lines.at[name, "s_nom_extendable"]


def test_deterministic_scales_lines_keeps_links(tmp_path):
    n = _canvas()
    x0 = n.lines["x"].copy()
    cfg = _write_csvs(tmp_path, {"L1": 100.0, "L3": 0.0}, {"K1": 250.0})
    bgt.build_deterministic_topology(n, cfg, clusters="x")
    n.calculate_dependent_values()

    assert n.lines.at["L1", "num_parallel"] == pytest.approx(1.0)  # 2 * 0.5
    assert n.lines.at["L1", "x"] == pytest.approx(x0["L1"] * 2.0)  # halved -> doubled
    assert "L3" not in n.lines.index  # zeroed -> dropped
    assert "L2" in n.lines.index and n.lines.at["L2", "x"] == pytest.approx(x0["L2"])
    # DC link: p_nom overridden, no impedance scaling (not a passive branch)
    assert n.links.at["K1", "p_nom"] == pytest.approx(250.0)
    assert not n.links.at["K1", "p_nom_extendable"]


def test_eev_blends_and_scales(tmp_path):
    n = _canvas()
    x0 = n.lines["x"].copy()
    cfg_a = _write_csvs(tmp_path / "a", {"L1": 200.0, "L3": 0.0}, {"K1": 500.0})
    cfg_b = _write_csvs(tmp_path / "b", {"L1": 0.0, "L3": 0.0}, {"K1": 0.0})
    scenarios = {
        "a": {"probability": 0.5, **cfg_a},
        "b": {"probability": 0.5, **cfg_b},
    }
    bgt.build_eev_topology(n, scenarios, clusters="x")
    n.calculate_dependent_values()

    # blended s_nom(L1) = 0.5*200 + 0.5*0 = 100 -> frac 0.5 -> x doubled
    assert n.lines.at["L1", "s_nom"] == pytest.approx(100.0)
    assert n.lines.at["L1", "x"] == pytest.approx(x0["L1"] * 2.0)
    assert "L3" not in n.lines.index  # zeroed in every scenario -> dropped
    assert n.links.at["K1", "p_nom"] == pytest.approx(250.0)


def test_stochastic_keeps_canvas_impedance(tmp_path):
    for sub, s1 in {"a": 100.0, "b": 0.0}.items():
        _write_csvs(tmp_path / sub, {"L1": s1}, {"K1": 500.0})
    n = _canvas()
    x_canvas = n.lines.at["L1", "x"]
    scenarios = {
        "a": {"probability": 0.5, "lines": str(tmp_path / "a" / "lines.csv"),
              "links": str(tmp_path / "a" / "links.csv")},
        "b": {"probability": 0.5, "lines": str(tmp_path / "b" / "lines.csv"),
              "links": str(tmp_path / "b" / "links.csv")},
    }
    bgt.build_stochastic_topology(n, scenarios, clusters="x")
    n.calculate_dependent_values()

    # canvas num_parallel (-> canvas x) preserved for every scenario slice
    np_by_scen = n.lines.xs("L1", level="name")["num_parallel"]
    x_by_scen = n.lines.xs("L1", level="name")["x"]
    assert np.allclose(np_by_scen.to_numpy(), 2.0)
    assert np.allclose(x_by_scen.to_numpy(), x_canvas)
    assert x_by_scen.std() == 0.0
    # but s_nom still overridden per scenario, and nothing dropped
    assert n.lines.xs("L1", level="name")["s_nom"].to_dict() == pytest.approx(
        {"a": 100.0, "b": 0.0}
    )
