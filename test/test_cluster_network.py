# SPDX-FileCopyrightText: Contributors to PyPSA-DE
# SPDX-License-Identifier: MIT

import geopandas as gpd
import pandas as pd
import pypsa
import pytest
from shapely.geometry import MultiPolygon, box

from scripts.cluster_network import clustering_for_n_clusters, update_bus_coordinates


def active(lines, year):
    return lines.s_nom[(lines.build_year <= year) & (year < lines.build_year + lines.lifetime)].sum()


@pytest.mark.parametrize("year", [2025, 2027, 2028, 2032, 2033, 2050])
def test_clustering_keeps_active_line_capacity(year):
    """Parallel lines retired by upgrades in different years keep their dates."""
    n = pypsa.Network()
    n.add("Bus", ["a1", "a2", "b1", "b2"], v_nom=380.0, x=[0, 0, 1, 1], y=[0, 0.1, 0, 0.1])
    n.add(
        "Line",
        ["old", "early", "late"],
        bus0=["a1", "a1", "a2"],
        bus1=["b1", "b1", "b2"],
        x=0.1,
        length=100,
        s_nom=[500.0, 1000.0, 3000.0],
        build_year=[0, 1990, 1990],
        lifetime=[float("inf"), 38.0, 43.0],
    )
    busmap = pd.Series({"a1": "A", "a2": "A", "b1": "B", "b2": "B"})
    clustered = clustering_for_n_clusters(n, busmap).n
    assert active(clustered.lines, year) == pytest.approx(active(n.lines, year))


def test_bus_coordinates_follow_busmap(tmp_path):
    """A bus in a minor part of its region keeps that region's coordinates."""
    shapes = gpd.GeoDataFrame(
        {"admin": ["A", "B"], "country": "DE"},
        geometry=[MultiPolygon([box(0, 0, 1, 1), box(1.2, 0, 1.3, 0.1)]), box(1.31, 0, 2, 1)],
        crs="EPSG:4326",
    )
    shapes.to_file(tmp_path / "admin.geojson")
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], x=[1.29, 1.6], y=[0.05, 0.5])
    update_bus_coordinates(n, pd.Series({"a": "A", "b": "B"}), tmp_path / "admin.geojson")
    assert n.buses.at["a", "x"] < 1 and n.buses.at["b", "x"] > 1.31
