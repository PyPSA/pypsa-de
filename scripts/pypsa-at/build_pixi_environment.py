"""A script to rebuild the pixi environment after upstream pinned version upgrades."""

import shutil
from pathlib import Path
from shutil import move
from subprocess import run


def main():
    # backup existing pixi files
    filenames = ("pixi.toml", "pixi.lock")
    for filename in filenames:
        fp = Path(filename)
        if fp.is_file():
            move(fp, f"{fp.stem}.bak{fp.suffix}")

    # purge existing installation
    if Path(".pixi").is_dir():
        shutil.rmtree(".pixi")

    # re-create pixi files from environment.yaml
    run(
        ["pixi", "init", "--import", "envs/environment.yaml", "--platform", "linux-64"],
        check=True,
    )

    # correct the environment name
    pixi_toml = Path("pixi.toml")
    pixi_toml.write_text(pixi_toml.read_text().replace("pypsa-de", "pypsa-at"))

    run(["pixi", "install"], check=True)

    # add pypsa-at packages
    packages = [
        "Pygments",
        "click",
        "folium",
        "frozendict",
        "markdown",
        "mkdocs",
        "mkdocs-marimo",
        "mkdocs-material",
        "mkdocs-material-extensions",
        "mkdocstrings-python",
        "mknotebooks",
        "pixi-pycharm",
        "ruff-lsp",
        "plotly",
        "pymdown-extensions",
        "pytest",
        "pytest-html",
        "pytest-cov",
        "pytest-xdist",
        "pytest-metadata",
        "xlsxwriter",
        "git-delta",
        "gitpython",
        "pandas-stubs",
    ]
    run(["pixi", "add"] + packages, check=True)

    pypi_packages = ["highspy", "xpress", "tsam", "mkdocs-badges"]
    run(["pixi", "add", "--pypi"] + pypi_packages, check=True)

    run(["pixi", "shell"])

    # format prompt to exclude conda env name
    run(["pixi config set shell.change-ps1", "false"])


if __name__ == "__main__":
    main()
