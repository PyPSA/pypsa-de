# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/pypsa-de/build_grid_topology.py, in
particular the scenario-specific scaling of line num_parallel/reactance/
resistance. The deterministic/eev variants also drop zeroed corridors; the
stochastic variant scales impedance per scenario slice and keeps zeroed
corridors (they are dropped per scenario at solve time by add_scenario_kvl).

Lines carry a standard ``type``, so ``n.calculate_dependent_values()``
recomputes ``x``/``r`` from ``num_parallel`` - exactly as the solver does -
which is why the tests assert ``x``/``r`` *after* that recompute.
"""

import importlib.util
import pathlib
import sys

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


def _canvas_with_countries():
    """
    DE + foreign (FR) buses, with DE-internal, foreign AC, and foreign DC
    branches carrying build_years, for the outside_de_grid feature.
    """
    n = pypsa.Network()
    n.add(
        "LineType", "t", r_per_length=0.01, x_per_length=0.1, c_per_length=10.0,
        f_nom=50.0, i_nom=1.0,
    )
    for b, c in [("DE1", "DE"), ("DE2", "DE"), ("FR1", "FR"), ("FR2", "FR")]:
        n.add("Bus", b, v_nom=380.0, carrier="AC", country=c)
    kw = dict(type="t", length=100.0, carrier="AC", s_nom_extendable=True)
    n.add("Line", "de", bus0="DE1", bus1="DE2", num_parallel=2.0, s_nom=200.0, build_year=2020, **kw)
    n.add("Line", "fr_old", bus0="FR1", bus1="FR2", num_parallel=2.0, s_nom=200.0, build_year=2025, **kw)
    n.add("Line", "fr_new", bus0="FR1", bus1="FR2", num_parallel=4.0, s_nom=400.0, build_year=2035, **kw)
    n.add("Link", "fr_dc_old", bus0="FR1", bus1="FR2", p_nom=500.0, carrier="DC", build_year=2025, p_nom_extendable=True)
    n.add("Link", "fr_dc_new", bus0="FR1", bus1="FR2", p_nom=700.0, carrier="DC", build_year=2040, p_nom_extendable=True)
    n.calculate_dependent_values()
    return n


def test_outside_de_grid_endogenous_is_noop():
    n = _canvas_with_countries()
    before = n.lines[["s_nom", "s_nom_extendable"]].copy()
    bgt.apply_outside_de_grid(n, "endogenous")
    assert n.lines[["s_nom", "s_nom_extendable"]].equals(before)
    assert n.links["p_nom_extendable"].all()


def test_outside_de_grid_exogen_fixes_foreign_grid():
    n = _canvas_with_countries()
    bgt.apply_outside_de_grid(n, "2030_exogen")

    # foreign committed-by-2030 branches: fixed, non-extendable
    assert n.lines.at["fr_old", "s_nom"] == pytest.approx(200.0)
    assert not n.lines.at["fr_old", "s_nom_extendable"]
    assert n.links.at["fr_dc_old", "p_nom"] == pytest.approx(500.0)
    assert not n.links.at["fr_dc_old", "p_nom_extendable"]
    # foreign not-yet-committed: AC line dropped, DC link floored to 0
    assert "fr_new" not in n.lines.index
    assert n.links.at["fr_dc_new", "p_nom"] == pytest.approx(0.0)
    assert not n.links.at["fr_dc_new", "p_nom_extendable"]
    # DE-internal branch untouched (still extendable)
    assert n.lines.at["de", "s_nom_extendable"]


def test_outside_de_grid_rejects_bad_spec():
    n = _canvas_with_countries()
    with pytest.raises(ValueError):
        bgt.apply_outside_de_grid(n, "2030_optimal")


def test_stochastic_scales_impedance_per_scenario(tmp_path):
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

    np_by_scen = n.lines.xs("L1", level="name")["num_parallel"]
    x_by_scen = n.lines.xs("L1", level="name")["x"]
    # scenario "a": s_nom 200 -> 100 (frac 0.5) scales num_parallel 2 -> 1, x doubled
    assert np_by_scen["a"] == pytest.approx(1.0)
    assert x_by_scen["a"] == pytest.approx(x_canvas * 2.0)
    # scenario "b": s_nom -> 0 keeps canvas impedance (dropped per scenario at solve)
    assert np_by_scen["b"] == pytest.approx(2.0)
    assert x_by_scen["b"] == pytest.approx(x_canvas)
    # impedance is genuinely scenario-specific, and nothing is dropped here
    assert x_by_scen.std() > 0.0
    assert n.lines.xs("L1", level="name")["s_nom"].to_dict() == pytest.approx(
        {"a": 100.0, "b": 0.0}
    )
