from pathlib import Path

import pypsa
import pytest


@pytest.fixture(scope="module")
def postnetwork():
    """Load the test data network."""
    csv_folder = (
        Path().resolve()
        / "data"
        / "pypsa"
        / "examples"
        / "ac-dc-meshed"
        / "results-lopf"
    )  # or "results-lpf"
    n = pypsa.Network(csv_folder)
    country_loc = [
        ("UK", "UK0 0"),
        ("UK", "UK0 1"),
        ("UK", "UK0 1"),
        ("UK", "UK0 2"),
        ("DE", "DE0 0"),
        ("DE", "DE0 0"),
        ("DE", "DE0 1"),
        ("NO", "NO0 0"),
        ("NO", "NO0 0"),
    ]
    n.buses["country"] = [c[0] for c in country_loc]
    n.buses["location"] = [c[1] for c in country_loc]
    n.links_t.p_set = n.links_t.p_set.drop(columns=n.links_t.p_set.columns)
    return n


def test_postnetwork(postnetwork):
    """"""
    print(postnetwork)
    assert True
