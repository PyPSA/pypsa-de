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

    # purge existing packages
    shutil.rmtree(".pixi")

    # re-create pixi files from environment.yaml
    run(["pixi", "init", "--import", "envs/environment.yaml"], check=True)
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
        "plotly",
        "pymdown-extensions",
        "pytest",
        "xlsxwriter",
    ]
    run(["pixi", "add"] + packages, check=True)
    # print(output.stdout)

    pypi_packages = ["highspy", "xpress", "tsam", "mkdocs-badges"]
    run(["pixi", "add", "--pypi"] + pypi_packages, check=True)
    # print(output.stdout)


if __name__ == "__main__":
    main()
