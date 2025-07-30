"""
Fixtures in this file are available to test functions as arguments.

Non-fixture functions are functions (hooks), that are
called during different stages of a test runs, e.g.:
  - during test collection
  - before test execution
  - after test execution

Here is the example section from the pytest docs:
https://docs.pytest.org/en/7.1.x/example/simple.html
"""

from pathlib import Path

import pytest

from evals.fileio import read_networks


@pytest.fixture(scope="session")
def result_path(pytestconfig) -> Path:
    """
    Retrieve the results path from CLI.

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
    return eval_path / "JSON"


@pytest.fixture(scope="session")
def networks(result_path: Path) -> dict:
    """Load the network."""
    return read_networks(result_path)


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
