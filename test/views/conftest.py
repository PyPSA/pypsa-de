"""Fixtures in this file are available to test functions as arguments.

Non-fixture functions are functions (hooks), that are
called during different stages of a test runs, e.g.:
  - during test collection
  - before test execution
  - after test execution

Here is the example section from the pytest docs:
https://docs.pytest.org/en/7.1.x/example/simple.html
"""

import getpass
import subprocess
from pathlib import Path

import pytest

from esmtools.fileio import read_networks


@pytest.fixture(scope="session")
def result_path(pytestconfig) -> Path:
    """Retrieve the results path from CLI.

    The default path assumes, that we run from inside a cloned
    pypsa-eur-sec repository.

    Note, that we cannot directly access
    the run_path (project root), because we want to run the tests on
    copied results folder as well. The run_path does not exist anymore
    after copying the results.
    """
    default_path = pytestconfig.rootpath / "tests" / "data"
    result_path = pytestconfig.getoption("result_path")
    return Path(result_path) if result_path else default_path


@pytest.fixture(scope="session")
def eval_path(result_path: Path) -> Path:
    """Retrieve the evaluation path from CLI argument."""
    return Path(result_path) / "esm_run" / "evaluation"


@pytest.fixture(scope="session")
def json_path(eval_path: Path) -> Path:
    """Build the JSON result path."""
    return eval_path / "JSONs_LC"


@pytest.fixture(scope="session")
def networks(result_path: Path) -> dict:
    """Load the postnetwork."""
    return read_networks(result_path)


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
    config._metadata["ESM result path"] = result_path
    config._metadata["ESM run comments"] = _read_comments_file(result_path)
    config._metadata["Git branch"] = _parse_git_branch()
    config._metadata["Git HEAD"] = _parse_git_head_hash()


def _parse_result_path(config) -> Path:
    """Parse the esm result path argument from CLI."""
    result_path_arg = [
        s.split("=", maxsplit=1)[-1]
        for s in config.invocation_params.args
        if s.startswith("--result-path=")
    ]
    result_path_default = Path(config.rootdir) / "results"
    return Path(result_path_arg[0]) if result_path_arg else result_path_default


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


def _read_comments_file(file_path: Path) -> str:
    """Read comments from run comment.txt file."""
    file_path_comment = file_path / "comment.txt"
    if not file_path_comment.is_file():
        return f"No comment file at {file_path_comment}."
    with Path.open(file_path_comment, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
        for nrow, line in enumerate(lines):
            if line.strip() == "### Comment:":
                return lines[nrow + 1]
    return f"No comment section found in {file_path_comment}."
