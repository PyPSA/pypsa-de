# -*- coding: utf-8 -*-
"""Command Line Interface to run evaluations.

Commands may be run from anywhere as long as the virtual environment is
activated and the esmtools project is installed.

Examples
--------
``` shell
# run a single evaluation by name
run_eval "/opt/data/esm/results" -n "eval_capacity_factor"
```

``` shell
# run multiple evaluations by name
run_eval "/opt/data/esm/results" -n "eval_capacity_factor" -n "eval_transmission_grid"
```

``` shell
# run all evaluations
run_eval "/opt/data/esm/results"
```
"""

import logging
import sys
from time import time

import click

logging.basicConfig(
    level=logging.INFO,
    format="{levelname} - {name} - {message}",
    datefmt="%Y-%m-%d %H:%M",
    style="{",
)
logger = logging.getLogger(__name__)


@click.command()
@click.argument("result_path", type=click.Path(exists=True), required=True)
@click.option(
    "--sub_directory",
    "-s",
    type=str,
    required=False,
    default="networks",
)
@click.option("--names", "-n", multiple=True, required=False, default=[])
@click.option(
    "--config_override", "-c", type=str, multiple=False, required=False, default=""
)
def run_eval(
    result_path: click.Path, sub_directory: str, names: list, config_override: str
) -> None:
    r"""Execute evaluation functions from the evals module.

    Find evaluation functions must be registered under
    evals.\__init__.\__all__ to be exposed and ultimately be found
    by this function. Keep that in mind when adding new evaluations.

    All evaluation function are expected to expose the same interface.
    The evaluation function arguments are listed in the evals module
    [reference section](evals/index.md).

    Parameters
    ----------
    result_path
        The path to the result folder, usually ./pypsa-eur-sec/results.
        Note, that running on copied result folders might fail
        due to missing resource files.
    sub_directory
        The subdirectory with the postnetwork files.
    names
        A list of evaluation name, e.g. "eval_electricity_amounts",
        optional. Defaults to running all evaluations from
        evals.__all__.
    config_override
        A path to a config.toml file with the same section as
        the config.defaults.toml used to override configurations
        used by view functions.

    Returns
    -------
    :
        Exits the program with the number of failed evaluations as exit
        code.
    """
    import evals.views as views
    from evals.fileio import read_networks, read_views_config

    eval_functions = [
        getattr(views, fn) for fn in views.__all__ if (not names or fn in names)
    ]
    n_evals = len(eval_functions)

    if n_evals == 0:
        sys.exit(f"Found no evaluation functions named: {names}")
    logger.info(f"Selected {n_evals} evaluation functions.")

    networks = read_networks(result_path, sub_directory=sub_directory)

    fails = []
    run_start = time()
    for i, func in enumerate(eval_functions, start=1):
        logger.info(f"({i}/{n_evals}) Start {func.__name__}...")
        eval_start = time()
        try:
            config = read_views_config(func, config_override)
            func(result_path=result_path, networks=networks, config=config)
        except Exception:  # noqa
            logger.exception(f"Exception during {func.__name__}.", exc_info=True)
            fails.append(func.__name__)
        else:
            logger.info(
                f"Executing {func.__name__} took {time() - eval_start:.2f} seconds."
            )
        finally:
            logger.info(f"Finish {func.__name__}.")

    logger.info(
        f"Full run took {time() - run_start:.2f} seconds."
        f"\nNumber of Errors: {len(fails)} {fails or ''}"
    )
    sys.exit(len(fails))


@click.command()
def run_tests() -> None:
    """Run test suite in a dev environment."""
    # delayed import to skip dependency in production environments
    import pytest

    rc = pytest.main()
    sys.exit(rc)


if __name__ == "__main__":
    # useful for debugging. Use run_eval command for production.
    args = ["../results/evals-dev"]  # this is the copied folder in local repo results
    # cp -r /mnt/storage/pypsa-at-AT10-365H results/
    # args = ["/mnt/storage/pypsa-at-AT10-365H"]  # This is the folder on the file share
    # args.extend(["-n", "view_heat_capacity"])
    args.extend(["-n", "view_hydrogen_balance"])

    # args.extend(["-c", "config.override.toml"])
    run_eval(args)
