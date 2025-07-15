# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import getpass
import pathlib
import subprocess
import zipfile
from functools import reduce
from shutil import unpack_archive
from urllib.request import urlretrieve
from uuid import uuid4

import geopandas as gpd
import pandas as pd
import pypsa
import pytest
import requests
import yaml


@pytest.fixture(scope="function")
def scigrid_network():
    return pypsa.examples.scigrid_de(from_master=True)


@pytest.fixture(scope="function")
def ac_dc_network():
    return pypsa.examples.ac_dc_meshed(from_master=True)


@pytest.fixture(scope="session")
def config():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config", "config.default.yaml")
    with open(path_config) as file:
        config_dict = yaml.safe_load(file)
    return config_dict


@pytest.fixture(scope="function")
def buses_dataframe():
    return pd.DataFrame(
        {
            "bus_id": [5231, 5232],
            "voltage": [380.0, 380.0],
            "dc": ["f", "f"],
            "symbol": ["Substation", "Substation"],
            "under_construction": ["f", "f"],
            "x": [6.8884, 6.8894],
            "y": [45.6783, 45.6793],
            "country": ["IT", "IT"],
            "geometry": ["POINT (6.8884 45.6783)", "POINT (6.8894 45.6793)"],
        },
        index=[5090, 5091],
    )


@pytest.fixture(scope="function")
def converters_dataframe():
    return pd.DataFrame(
        {
            "converter_id": "convert_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def lines_dataframe():
    return pd.DataFrame(
        {
            "line_id": "line_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "circuits": 1.0,
            "length": 1000.0,
            "underground": "t",
            "under_construction": "f",
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def links_dataframe():
    return pd.DataFrame(
        {
            "link_id": "link_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "p_nom": 600.0,
            "length": 1000.0,
            "underground": "t",
            "under_construction": "f",
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def transformers_dataframe():
    return pd.DataFrame(
        {
            "transformer_id": "transf_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def download_natural_earth(tmpdir):
    url = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries_deu.zip"
    directory_to_extract_to = pathlib.Path(tmpdir, "folder")
    zipped_filename = "ne_10m_admin_0_countries_deu.zip"
    path_to_zip_file, headers = urlretrieve(url, zipped_filename)
    with zipfile.ZipFile(path_to_zip_file, "r") as zip_ref:
        zip_ref.extractall(directory_to_extract_to)
    natural_earth_shape_file_path = pathlib.Path(
        directory_to_extract_to, "ne_10m_admin_0_countries_deu.shp"
    )
    yield natural_earth_shape_file_path
    pathlib.Path(zipped_filename).unlink(missing_ok=True)
    pathlib.Path(natural_earth_shape_file_path).unlink(missing_ok=True)


@pytest.fixture(scope="function")
def download_eez(tmpdir):
    name = str(uuid4())[:8]
    org = str(uuid4())[:8]
    zipped_filename = "World_EEZ_v12_20231025_LR.zip"
    response = requests.post(
        "https://www.marineregions.org/download_file.php",
        params={"name": zipped_filename},
        data={
            "name": name,
            "organisation": org,
            "email": f"{name}@{org}.org",
            "country": "Germany",
            "user_category": "academia",
            "purpose_category": "Research",
            "agree": "1",
        },
    )
    zipped_filename_path = pathlib.Path(tmpdir, zipped_filename)
    with open(zipped_filename_path, "wb") as f:
        f.write(response.content)
    unpack_archive(zipped_filename_path, tmpdir)
    output_path = pathlib.Path(
        tmpdir, "World_EEZ_v12_20231025_LR", "eez_v12_lowres.gpkg"
    )
    yield output_path
    pathlib.Path(output_path).unlink(missing_ok=True)


@pytest.fixture(scope="function")
def italy_shape(download_natural_earth, tmpdir):
    shape_file = gpd.read_file(download_natural_earth)
    fieldnames = (
        shape_file[x].where(lambda s: s != "-99")
        for x in ("ISO_A2", "WB_A2", "ADM0_A3")
    )
    shape_file["name"] = reduce(
        lambda x, y: x.fillna(y), fieldnames, next(fieldnames)
    ).str[:2]
    italy_shape_file = shape_file.loc[shape_file.name.isin(["IT"])]
    italy_shape_file_path = pathlib.Path(tmpdir, "italy_shape.geojson")
    italy_shape_file.to_file(italy_shape_file_path, driver="GeoJSON")
    yield italy_shape_file_path


def pytest_addoption(parser) -> None:
    """Register command line arguments."""
    parser.addoption(
        "--result-path", action="store", help="Path to the ESM results folder."
    )


def pytest_configure(config) -> None:
    """Add environment info to HTML report."""
    # The Environment section seems to be broken
    # with pytest-metadata >= 3.x. We must use version < 3
    # https://github.com/pytest-dev/pytest-html/issues/683
    config._metadata["Username"] = getpass.getuser()
    markers = config.option.markexpr or "No Markers"
    config._metadata["Pytest markers"] = markers
    config._metadata["Src directory"] = config.rootdir
    result_path = _parse_result_path(config)
    config._metadata["results path"] = result_path
    config._metadata["Run name"] = result_path.parent
    config._metadata["Git branch"] = _parse_git_branch()
    config._metadata["Git HEAD"] = _parse_git_head_hash()


def _parse_result_path(config) -> pathlib.Path:
    """Parse the esm result path argument from CLI."""
    result_path_arg = [
        s.split("=", maxsplit=1)[-1]
        for s in config.invocation_params.args
        if s.startswith("--result-path=")
    ]
    result_path_default = pathlib.Path(config.rootdir) / "results"
    return pathlib.Path(result_path_arg[0]) if result_path_arg else result_path_default


def _parse_git_branch() -> str:
    """Return the active Git branch."""
    sp = subprocess.run("git branch", capture_output=True)
    lines = sp.stdout.decode("utf-8").split("\n")
    branch = [line for line in lines if line.startswith("*")]
    return branch[0] if branch else "No branch detected."


def _parse_git_head_hash() -> str:
    """Return the has of the current HEAD."""
    sp = subprocess.run("git rev-parse HEAD", capture_output=True)
    return sp.stdout.decode("utf-8").strip() or "No HEAD detected."
