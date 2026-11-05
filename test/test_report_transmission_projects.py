# SPDX-FileCopyrightText: Contributors to PyPSA-DE
# SPDX-License-Identifier: MIT

from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

from scripts.report_transmission_projects import (
    active_capacity,
    add_region_mapping,
    build_ledger,
    endpoint_order,
    identify_source,
    load_regions,
    locate_buses,
    project_metadata,
    visible_projects,
)


@pytest.fixture
def project_inputs(tmp_path):
    buses = pd.DataFrame(
        {"country": ["DE", "DE", "FR"], "x": [8.0, 10.0, 6.0], "y": [50.0, 52.0, 49.0]},
        index=["a", "b", "c"],
    )
    lines = pd.DataFrame(
        {
            "bus0": ["a"],
            "bus1": ["b"],
            "s_nom": [100.0],
            "build_year": [0],
            "lifetime": [np.inf],
        },
        index=["old"],
    )
    links = pd.DataFrame(
        columns=["carrier", "p_nom", "bus0", "bus1", "build_year", "lifetime"]
    )
    for name in ("new_buses", "new_links", "adjust_links"):
        pd.DataFrame().to_csv(tmp_path / f"{name}.csv")
    pd.DataFrame({"build_year": [1990], "lifetime": [38]}, index=["old"]).to_csv(
        tmp_path / "adjust_lines.csv"
    )
    pd.DataFrame(
        {
            "bus0": ["a"],
            "bus1": ["b"],
            "s_nom": [250.0],
            "build_year": [2028],
            "lifetime": [np.inf],
            "project_id": ["NEP42"],
            "source_file": ["nep/upgraded_lines.csv"],
            "replaced_asset": ["old"],
        },
        index=["old_upgraded"],
    ).to_csv(tmp_path / "new_lines.csv")
    return SimpleNamespace(buses=buses, lines=lines, links=links), tmp_path


def test_upgrade_net_and_build_year_boundary(project_inputs):
    n, path = project_inputs
    ledger, _, _ = build_ledger(n, path, pd.DataFrame(), 2025, 2028)
    row = ledger.iloc[0]
    assert row.gross_added == 250
    assert row.retired == 100
    assert row.net_added == 150
    assert row.project_id == "NEP42"
    assert row.provenance == "recorded"
    later, _, _ = build_ledger(n, path, pd.DataFrame(), 2028, 2030)
    assert later.net_added.sum() == 0
    assert not later.changed.any()


def test_zero_nominal_link_and_interconnector(project_inputs):
    n, path = project_inputs
    pd.DataFrame(
        {
            "bus0": ["a", "a"],
            "bus1": ["c", "c"],
            "p_nom": [0.0, 700.0],
            "build_year": [2027, 2027],
        },
        index=["zero", "positive"],
    ).to_csv(path / "new_links.csv")
    ledger, _, _ = build_ledger(n, path, pd.DataFrame(), 2025, 2030)
    dc = ledger.set_index("asset").loc[["zero", "positive"]]
    assert dc.category.eq("DC_interconnectors").all()
    assert dc.net_added.sum() == 700
    assert not dc.loc["zero", "changed"]


def test_missing_retirement_rejected(project_inputs):
    n, path = project_inputs
    pd.DataFrame().to_csv(path / "adjust_lines.csv")
    with pytest.raises(ValueError, match="missing retirement"):
        build_ledger(n, path, pd.DataFrame(), 2025, 2030)


def test_lifetime_and_reverse_interval(project_inputs):
    assert (
        active_capacity(
            pd.Series(dict(p_nom=500, build_year=2020, lifetime=10)), "p_nom", 2030
        )
        == 0
    )
    n, path = project_inputs
    with pytest.raises(ValueError, match="greater"):
        build_ledger(n, path, pd.DataFrame(), 2030, 2025)


def test_legacy_upgrade_provenance_requires_unique_match():
    row = pd.Series(dict(tags='"proj_id"=>"10"', build_year=2028))
    catalog = pd.DataFrame(
        [
            dict(
                project_id="42",
                source_file="x/upgraded_links.csv",
                project_action="upgraded_links",
                **row,
            )
        ]
    )
    assert identify_source("old_upgraded", row, "upgraded_links", catalog) == (
        "x/upgraded_links.csv",
        "42",
        "unique legacy match",
    )
    duplicate = pd.concat([catalog, catalog])
    assert (
        identify_source("old_upgraded", row, "upgraded_links", duplicate)[2]
        == "ambiguous"
    )


def test_within_nuts2_project_retained_and_model_collapse(project_inputs):
    n, path = project_inputs
    ledger, _, _ = build_ledger(n, path, pd.DataFrame(), 2025, 2030)
    locations = pd.DataFrame({"region": ["DE11", "DE11", "FR10"]}, index=n.buses.index)
    ledger = add_region_mapping(
        ledger, locations, pd.Series({"a": "DE11", "b": "DE11"})
    )
    assert ledger.within_nuts2.all()
    assert ledger.model_mapping.eq("collapsed within model bus").all()
    assert ledger.net_added.sum() == 150


def test_upgrade_builder_preserves_original_identity_and_supplied_capacity():
    from scripts.build_transmission_projects import get_upgraded_lines

    original = pd.DataFrame(
        {"bus0": ["a"], "bus1": ["b"], "p_nom": [1000.0]}, index=["old"]
    )
    n = SimpleNamespace(df=lambda component: original)
    project = pd.DataFrame(
        {
            "p_nom": [2000.0],
            "project_id": ["TYNDP42"],
            "source_file": ["tyndp/upgraded_links.csv"],
            "replaced_asset": ["old"],
        },
        index=["TYNDP42"],
    )
    result = get_upgraded_lines("Link", n, project, pd.Series({"TYNDP42": "old"}))
    assert result.loc["old_upgraded", "project_id"] == "TYNDP42"
    assert result.loc["old_upgraded", "replaced_asset"] == "old"
    assert result.loc["old_upgraded", "p_nom"] == 2000


@pytest.mark.parametrize("resolution", ["nuts2", "nuts3", "grid"])
def test_display_resolution_keeps_capacity_and_grid_coordinates(
    project_inputs, resolution
):
    n, path = project_inputs
    # Both physical endpoints are inside one administrative region.
    regions = gpd.GeoDataFrame(
        geometry=[box(0, 40, 20, 60)],
        index=["DE11" if resolution == "nuts2" else "DE111"],
        crs=4326,
    ).to_crs(3035)
    locations = locate_buses(n.buses, regions, resolution)
    ledger, _, _ = build_ledger(n, path, pd.DataFrame(), 2025, 2030)
    ledger = add_region_mapping(ledger, locations, resolution=resolution)
    assert ledger.net_added.sum() == 150
    if resolution == "grid":
        physical = gpd.GeoSeries(
            gpd.points_from_xy(n.buses.x, n.buses.y), crs=4326
        ).to_crs(3035)
        np.testing.assert_allclose(locations.x, physical.x)
        np.testing.assert_allclose(locations.y, physical.y)
        assert ledger.region0.iloc[0] == "a"
        assert not ledger.within_region.any()
    else:
        assert ledger.within_region.all()
        assert locations.loc["a", "x"] == locations.loc["b", "x"]


def test_wrong_nuts_level_rejected(tmp_path):
    path = tmp_path / "nuts2.geojson"
    gpd.GeoDataFrame(
        {"NUTS_ID": ["DE11"], "LEVL_CODE": [2]}, geometry=[box(8, 50, 10, 52)], crs=4326
    ).to_file(path)
    with pytest.raises(ValueError, match="NUTS3"):
        load_regions(path, level=3)


def test_cli_default_intervals_have_separate_outputs(monkeypatch, tmp_path):
    from scripts import report_transmission_projects as report

    calls = []
    monkeypatch.setattr(report, "generate_report", lambda *args: calls.append(args))
    monkeypatch.setattr(
        "sys.argv",
        [
            "report",
            "--base-network",
            "base.nc",
            "--projects",
            "projects",
            "--output-dir",
            str(tmp_path),
        ],
    )
    report.main()
    assert [(args[3], args[4]) for args in calls] == [(2025, 2030), (2030, 2035)]
    assert [args[5] for args in calls] == [
        tmp_path / "2025_until_2030",
        tmp_path / "2030_until_2035",
    ]


def test_cli_requires_both_custom_years(monkeypatch):
    from scripts import report_transmission_projects as report

    monkeypatch.setattr(
        "sys.argv",
        [
            "report",
            "--base-network",
            "base.nc",
            "--projects",
            "projects",
            "--start-year",
            "2025",
        ],
    )
    with pytest.raises(SystemExit) as error:
        report.main()
    assert error.value.code == 2


def test_nep_metadata_uses_reference_and_source_endpoints():
    row = SimpleNamespace(
        project_id="TTG-P67M102",
        source_file="nep/new_lines.csv",
        tags='{"url":"https://example.org/nep.pdf", "origin":"NEP", "NEP-ID":"TTG-P67", "location1":"Altheim"}',
    )
    catalog = pd.DataFrame(
        [
            dict(
                project_id=row.project_id,
                source_file=row.source_file,
                bus0="Simbach",
                bus1="Altheim",
            )
        ]
    )
    metadata = project_metadata(row, catalog)
    assert metadata == dict(
        source_plan="NEP",
        reference_id="TTG-P67",
        project_name="Simbach → Altheim (endpoints)",
        source_url="https://example.org/nep.pdf",
        source_page=None,
        source_endpoint0_name="Simbach",
        source_endpoint1_name="Altheim",
        endpoint0_name="Simbach",
        endpoint1_name="Altheim",
    )


def test_tyndp_hstore_and_manual_metadata():
    row = SimpleNamespace(
        project_id="TYNDP2020_38",
        source_file="tyndp2020/new_links.csv",
        tags='"url"=>"https://example.org/176", "tyndp2020_proj_id"=>"176"',
    )
    metadata = project_metadata(row, pd.DataFrame())
    assert metadata["source_plan"] == "TYNDP 2020"
    assert metadata["reference_id"] == "176"
    assert metadata["project_name"] == "Name not supplied"
    row = SimpleNamespace(
        project_id="TYNDP2024_153",
        source_file="manual/new_links.csv",
        tags="{name:France-Alderney-Britain,url:https://example.org/portfolio.xlsx}",
    )
    metadata = project_metadata(row, pd.DataFrame())
    assert metadata["source_plan"] == "TYNDP 2024 (manual CSV)"
    assert metadata["project_name"] == "France-Alderney-Britain"
    assert metadata["source_url"] == "https://example.org/portfolio.xlsx"


def test_visibility_excludes_display_and_model_collapsed_branches():
    ledger = pd.DataFrame(
        {
            "in_scope": [True, True, True, True, False, True],
            "changed": [True, True, True, True, True, False],
            "within_region": [True, False, False, False, False, False],
            "model_mapping": [
                "between model buses",
                "collapsed within model bus",
                "between model buses",
                "unmapped",
                "between model buses",
                "between model buses",
            ],
        }
    )
    assert visible_projects(ledger).tolist() == [False, False, True, True, False, False]
    assert visible_projects(ledger.drop(columns="model_mapping")).tolist() == [
        False,
        True,
        True,
        True,
        False,
        False,
    ]


def test_undated_project_retained_without_assuming_commissioning(project_inputs):
    n, path = project_inputs
    pd.DataFrame(
        {"bus0": ["a"], "bus1": ["c"], "p_nom": [1000.0], "build_year": [np.nan]},
        index=["undated"],
    ).to_csv(path / "new_links.csv")
    ledger, _, _ = build_ledger(n, path, pd.DataFrame(), 2025, 2030)
    row = ledger.set_index("asset").loc["undated"]
    assert not row.commissioning_known
    assert pd.isna(row.capacity_start) and pd.isna(row.net_added)
    assert not row.changed


def test_new_source_url_and_tyndp_substation_labels():
    row = SimpleNamespace(
        project_id="TYNDP2024_1014",
        source_file="tyndp2024draft/new_links.csv",
        tags='{"legacy_source_url":"https://example.org/old", "source_url":"https://example.org/current", "substation_from":"Verderio", "substation_to":"Bonaduz"}',
    )
    catalog = pd.DataFrame(
        [
            {
                "project_id": row.project_id,
                "source_file": row.source_file,
                "bus0": np.nan,
                "bus1": np.nan,
            }
        ]
    )
    metadata = project_metadata(row, catalog)
    assert metadata["source_url"] == "https://example.org/current"
    assert metadata["source_page"] is None
    assert metadata["project_name"] == "Verderio → Bonaduz (endpoints)"


def test_legend_wrap_preserves_long_urls_within_available_width():
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath

    from scripts.report_transmission_projects import wrap_legend

    url = "https://example.org/" + "long-path-" * 60
    wrapped = wrap_legend(url, width=500)
    assert wrapped.replace("\n", "") == url
    font = FontProperties(family="DejaVu Sans", size=7)
    assert all(
        TextToPath().get_text_width_height_descent(line, font, False)[0] <= 500
        for line in wrapped.splitlines()
    )


def test_source_page_from_tags_and_lookup_override(tmp_path):
    from scripts.report_transmission_projects import apply_source_links

    tags = '{"source_url": "https://example.org/nep.pdf", "source_pdf_page": 186}'
    rows = [
        SimpleNamespace(project_id="NEP_A", source_file="nep2025/x.csv", tags=tags),
        SimpleNamespace(project_id="P43", source_file="nep/new_lines.csv", tags="{}"),
    ]
    metadata = pd.DataFrame([project_metadata(row, pd.DataFrame()) for row in rows])
    assert metadata.source_page.tolist()[0] == 186
    pd.DataFrame(
        {
            "dataset": ["nep"],
            "project_id": ["P43"],
            "source_url": ["https://example.org/confirmation.pdf"],
            "source_pdf_page": [117],
        }
    ).to_csv(tmp_path / "links.csv", index=False)
    ledger = pd.DataFrame(
        {"project_id": ["NEP_A", "P43"], "source_file": [r.source_file for r in rows]}
    )
    linked = apply_source_links(metadata, ledger, tmp_path / "links.csv")
    assert linked.source_url.tolist() == [
        "https://example.org/nep.pdf",
        "https://example.org/confirmation.pdf",
    ]
    assert linked.source_page.tolist() == [186, 117]


def test_continuation_pages_have_only_legend_and_preserve_all_projects(
    monkeypatch, project_inputs
):
    import re

    from scripts.report_transmission_projects import draw_report

    n, path = project_inputs
    ledger, _, _ = build_ledger(n, path, pd.DataFrame(), 2025, 2030)
    locations = pd.DataFrame(
        {
            "region": ["DE11", "DE12", "FR10"],
            "x": [0.0, 1000.0, 2000.0],
            "y": [0.0, 1000.0, 2000.0],
        },
        index=n.buses.index,
    )
    ledger = add_region_mapping(ledger, locations)
    ledger = pd.concat(
        [
            ledger.assign(asset=f"asset_{i}", project_id=f"project_{i}")
            for i in range(30)
        ],
        ignore_index=True,
    )
    ledger = ledger.assign(
        source_plan="NEP",
        reference_id="P1",
        project_name="Example endpoints",
        source_url="https://example.org/project.pdf",
        source_page=135,
        endpoint0_name="Start",
        endpoint1_name="End",
    )
    omitted = pd.concat(
        [
            ledger.iloc[[0]].assign(
                asset=f"omitted_{i}", project_id="omitted-project", within_region=True
            )
            for i in range(2)
        ],
        ignore_index=True,
    )
    ledger = pd.concat([ledger, omitted], ignore_index=True)
    regions = gpd.GeoDataFrame(
        {"NAME_LATN": ["Example region"]},
        geometry=[box(-100, -100, 1100, 1100)],
        index=["DE11"],
        crs=3035,
    )
    pages = []

    class CapturePDF:
        def __init__(self, path):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def savefig(self, fig):
            assert "DE11\nExample region" in [
                text.get_text() for text in fig.axes[0].texts
            ]
            map_labels = [
                text for text in fig.axes[0].texts if text.get_text().isdigit()
            ]
            assert len(map_labels) == 30
            assert all(text.get_text().isdigit() for text in map_labels)
            map_colors = {text.get_text(): text.get_color() for text in map_labels}
            assert len(set(map_colors.values())) == 30
            for patch, label in zip(fig.axes[0].patches, map_labels, strict=True):
                np.testing.assert_allclose(patch.get_edgecolor()[:3], label.get_color())
            for text in fig.texts:
                if text.get_url():
                    assert text.get_text() == "p. 135 \u2197"
                    assert text.get_url() == "https://example.org/project.pdf#page=135"
                if text.get_text().startswith("["):
                    key = text.get_text().split("]", 1)[0][1:]
                    assert text.get_color() == map_colors[key]
            pages.append(
                (
                    fig.axes[0].get_visible(),
                    "\n".join(text.get_text() for text in fig.texts),
                )
            )

    monkeypatch.setattr("matplotlib.backends.backend_pdf.PdfPages", CapturePDF)
    draw_report(
        ledger,
        pd.DataFrame(),
        locations,
        regions,
        "AC_internal",
        2025,
        2030,
        "DE",
        path / "report.pdf",
    )
    assert len(pages) > 1
    assert [visible for visible, _ in pages] == [True] + [False] * (len(pages) - 1)
    ids = re.findall(r"\bproject_\d+\b", "\n".join(text for _, text in pages))
    assert len(ids) == 30 and set(ids) == {f"project_{i}" for i in range(30)}
    legend_text = " ".join(" ".join(text for _, text in pages).split())
    assert "Folder: nep" in legend_text
    assert "https://" not in legend_text
    assert legend_text.count("p. 135") == 30
    assert legend_text.count("omitted-project (2)") == 1


@pytest.mark.parametrize("reverse", [False, True])
def test_endpoint_names_follow_physical_bus_order(reverse):
    source_geometry = "LINESTRING (11.122879 50.998664, 9.746033 50.917545)"
    buses = pd.DataFrame(
        {"x": [11.122879, 9.746033], "y": [50.998664, 50.917545]}, index=["V", "M"]
    )
    if reverse:
        buses = buses.iloc[::-1]
    order = endpoint_order(source_geometry, buses)
    assert order == ("reversed" if reverse else "same")
    row = SimpleNamespace(
        project_id="50HzT-P37M25a",
        source_file="nep2023confirmed/upgraded_lines.csv",
        tags='{"origin":"NEP", "location0":"Vieselbach", "location1":"Mecklar"}',
        endpoint_order=order,
    )
    metadata = project_metadata(row, pd.DataFrame())
    expected = ["Mecklar", "Vieselbach"] if reverse else ["Vieselbach", "Mecklar"]
    assert [metadata["endpoint0_name"], metadata["endpoint1_name"]] == expected
    assert metadata["project_name"] == " → ".join(expected) + " (endpoints)"
    assert metadata["source_endpoint0_name"] == "Vieselbach"
    assert metadata["source_endpoint1_name"] == "Mecklar"


def test_ambiguous_endpoint_order_does_not_assign_source_names_to_buses():
    buses = pd.DataFrame({"x": [10.0, 10.0], "y": [51.0, 51.0]})
    assert endpoint_order("LINESTRING (9 51, 11 51)", buses) == "unknown"
    row = SimpleNamespace(
        project_id="P",
        source_file="nep/new_lines.csv",
        tags='{"location0":"A", "location1":"B"}',
        endpoint_order="unknown",
        bus0="bus-a",
        bus1="bus-b",
    )
    metadata = project_metadata(row, pd.DataFrame())
    assert metadata["endpoint0_name"] == "bus-a"
    assert metadata["endpoint1_name"] == "bus-b"
    assert "order unresolved" in metadata["project_name"]


@pytest.mark.parametrize("shared_endpoint", [False, True])
def test_offshore_endpoints_create_unique_buses_and_reuse_existing(shared_endpoint):
    import pypsa
    from shapely.geometry import LineString

    from scripts.build_transmission_projects import connect_new_lines

    n = pypsa.Network()
    n.add("Bus", "land", x=0.0, y=0.0, carrier="AC", country="DE")
    end = (2.0, 2.0) if shared_endpoint else (5.0, 5.0)
    lines = pd.DataFrame(
        {"geometry": [LineString([(2.0, 2.0), end])]}, index=["project"]
    )
    offshore = gpd.GeoDataFrame(
        {"country": ["DE"]}, geometry=[box(1, 1, 6, 6)], crs=4326
    )
    connected, new_buses = connect_new_lines(
        lines, n, pd.DataFrame(), offshore, distance_upper_bound=0.4
    )
    assert new_buses.index.is_unique
    assert len(new_buses) == (1 if shared_endpoint else 2)
    assert (
        connected.loc["project", "bus0"] == connected.loc["project", "bus1"]
    ) == shared_endpoint
    assert tuple(n.buses.loc[connected.loc["project", "bus1"], ["x", "y"]]) == end
    # Reusing the same offshore sites in another project must not export duplicates.
    again = lines.rename(index={"project": "other"})
    connected_again, exported = connect_new_lines(
        again, n, new_buses, offshore, distance_upper_bound=0.4
    )
    assert len(exported) == len(new_buses)
    assert exported.index.is_unique
    assert connected_again.loc["other", "bus0"] == connected.loc["project", "bus0"]


def test_upgrades_match_original_network_not_new_project(monkeypatch, tmp_path):
    import pypsa

    from scripts import build_transmission_projects as build

    n = pypsa.Network()
    n.add(
        "Bus",
        ["a", "b", "c", "d"],
        x=[8.0, 9.0, 8.1, 9.1],
        y=[50.0, 51.0, 50.1, 51.1],
        carrier="AC",
        v_nom=380,
    )
    n.add("Line", "existing", bus0="a", bus1="b", s_nom=1000, v_nom=380)
    candidates = {"Line": n.lines.copy(), "Link": n.links.copy()}
    n.add("Line", "proposed", bus0="c", bus1="d", s_nom=2000)
    monkeypatch.setattr(build, "n", n, raising=False)
    source = tmp_path / "source"
    source.mkdir()
    pd.DataFrame(
        {
            "x0": [8.1],
            "y0": [50.1],
            "x1": [9.1],
            "y1": [51.1],
            "build_year": [2030],
            "project_status": ["confirmed"],
            "v_nom": [380.0],
            "num_parallel": [float("nan")],
        },
        index=["upgrade"],
    ).to_csv(source / "upgraded_lines.csv")
    result = build.add_projects(
        n,
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        box(0, 40, 20, 60),
        gpd.GeoDataFrame(geometry=[], crs=4326),
        source,
        "source",
        status=["confirmed"],
        upgrade_candidates=candidates,
    )
    assert result[0].index.tolist() == ["existing_upgraded"]
    assert result[2].index.tolist() == ["existing"]


def test_split_corridor_upgrade_is_chained_and_unmatched_is_logged(
    monkeypatch, tmp_path, caplog
):
    import pypsa

    from scripts import build_transmission_projects as build

    n = pypsa.Network()
    n.add(
        "Bus",
        ["a", "m", "b", "c", "d", "e", "f"],
        x=[8.0, 8.5, 9.0, 10.0, 10.1, 11.0, 11.1],
        y=[50.0, 50.5, 51.0, 50.0, 50.05, 50.0, 50.05],
        carrier="AC",
        v_nom=220,
    )
    n.add(
        "Line",
        ["am", "mb", "cd", "ef"],
        bus0=["a", "m", "c", "e"],
        bus1=["m", "b", "d", "f"],
        length=[70.0, 70.0, 8.0, 8.0],
        v_nom=220,
        s_nom=500,
    )
    monkeypatch.setattr(build, "n", n, raising=False)
    source = tmp_path / "source"
    source.mkdir()
    pd.DataFrame(
        {
            "x0": [8.0, 10.0, 11.0, 15.0],
            "y0": [50.0, 50.0, 50.0, 45.0],
            "x1": [9.0, 10.1, 11.1, 16.0],
            "y1": [51.0, 50.05, 50.05, 46.0],
            "length": [150.0, 8.0, 8.0, 100.0],
            "build_year": [2030] * 4,
            "project_status": ["confirmed"] * 4,
            "v_nom": [220.0] * 4,
            "num_parallel": [float("nan")] * 4,
        },
        index=["chain", "short1", "short2", "lost"],
    ).to_csv(source / "upgraded_lines.csv")
    empty = [pd.DataFrame()] * 5
    result = build.add_projects(
        n,
        *empty,
        box(0, 40, 20, 60),
        gpd.GeoDataFrame(geometry=[], crs=4326),
        source,
        "source",
        status=["confirmed"],
    )
    upgraded = result[0].set_index("project_id")
    assert sorted(upgraded.loc["chain", "replaced_asset"]) == ["am", "mb"]
    # chained branches keep their own length, not the project's total length
    assert upgraded.loc["chain", "length"].tolist() == [70.0, 70.0]
    assert "lost" not in upgraded.index
    assert sorted(result[2].index) == ["am", "cd", "ef", "mb"]
    assert "['lost']" in caplog.text


def test_upgrades_of_one_branch_are_staged_and_counted_once(tmp_path):
    from scripts.build_transmission_projects import stage_upgrades

    upgrades = pd.DataFrame(
        {
            "bus0": ["a"] * 3,
            "bus1": ["b"] * 3,
            "s_nom": [300.0, 200.0, 250.0],
            "build_year": [2035, 2028, 2028],
            "project_id": ["NEP2", "NEP1", "TYNDP1"],
            "replaced_asset": ["old"] * 3,
        },
        index=["old_upgraded"] * 3,
    )
    adjust = pd.DataFrame(
        {"build_year": [1990] * 3, "lifetime": [45, 38, 38]}, index=["old"] * 3
    )
    staged, adjust = stage_upgrades(upgrades, adjust)
    assert staged.index.tolist() == ["old_upgraded", "old_upgraded2"]
    assert staged.project_id.tolist() == ["NEP1", "NEP2"]
    assert staged.replaced_asset.tolist() == ["old", "old_upgraded"]
    assert staged.lifetime.tolist() == [7, np.inf]
    assert adjust.lifetime.tolist() == [38]

    buses = pd.DataFrame(
        {"country": ["DE", "DE"], "x": [8.0, 10.0], "y": [50.0, 52.0]},
        index=["a", "b"],
    )
    lines = pd.DataFrame(
        {"bus0": ["a"], "bus1": ["b"], "s_nom": [100.0], "build_year": [0]},
        index=["old"],
    ).assign(lifetime=np.inf)
    links = pd.DataFrame(
        columns=["carrier", "p_nom", "bus0", "bus1", "build_year", "lifetime"]
    )
    n = SimpleNamespace(buses=buses, lines=lines, links=links)
    for name in ("new_buses", "new_links", "adjust_links"):
        pd.DataFrame().to_csv(tmp_path / f"{name}.csv")
    staged.to_csv(tmp_path / "new_lines.csv")
    adjust.to_csv(tmp_path / "adjust_lines.csv")
    totals = [
        build_ledger(n, tmp_path, pd.DataFrame(), *interval)[0].net_added.sum()
        for interval in ((2025, 2030), (2030, 2040))
    ]
    assert totals == [100, 100]


def test_project_twin_is_matched_by_source_geometry_not_snapped_buses(monkeypatch):
    from shapely.geometry import LineString

    from scripts import build_transmission_projects as build

    buses = pd.DataFrame({"x": [8.0, 9.0], "y": [50.0, 51.0]}, index=["a", "b"])
    monkeypatch.setattr(build, "n", SimpleNamespace(buses=buses), raising=False)
    source = LineString([(8.2, 50.2), (9.2, 51.2)])
    added = pd.DataFrame(
        {"bus0": ["a"], "bus1": ["b"], "p_nom": [2000.0], "geometry": [source]},
        index=["NEP_link"],
    )
    twin = pd.DataFrame(
        {"p_nom": [2000.0], "geometry": [LineString(source.coords[::-1])]},
        index=["TYNDP_link"],
    )
    match = build.find_closest_lines(added, twin, distance_upper_bound=0.1)
    assert match.to_dict() == {"TYNDP_link": "NEP_link"}


@pytest.mark.parametrize(
    "v_nom, expected", [(380.0, ""), (220.0, "Al/St 240/40 2-bundle 220.0")]
)
def test_upgrade_keeps_base_line_type_only_at_unchanged_voltage(v_nom, expected):
    from scripts.build_transmission_projects import get_upgraded_lines

    original = pd.DataFrame(
        {
            "bus0": ["a"],
            "bus1": ["b"],
            "v_nom": [220.0],
            "type": ["Al/St 240/40 2-bundle 220.0"],
            "s_nom": [500.0],
        },
        index=["old"],
    )
    n = SimpleNamespace(df=lambda component: original)
    project = pd.DataFrame({"v_nom": [v_nom], "type": [np.nan]}, index=["NEP1"])
    result = get_upgraded_lines("Line", n, project, pd.Series({"NEP1": "old"}))
    assert result["type"].fillna("").iloc[0] == expected
