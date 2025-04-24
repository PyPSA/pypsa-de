"""A script to rebuild the pixi environment after upstream pinned version upgrades."""

from subprocess import run
from shutil import move


def main():
    # backup existing pixi files
    move("pixi.toml", "pixi.toml.bak")
    move("pixi.lock", "pixi.lock.bak")

    # re-create pixi files from environment.yaml
    run(["pixi", "init", "--import", "envs/linux-pinned.yaml"], check=True)
    run(["pixi", "install"], check=True)

    # add pypsa-at packages
    packages = ['Pygments',
                'click',
                'folium',
                'frozendict',
                'markdown',
                'mkdocs',
                'mkdocs-marimo',
                'mkdocs-material',
                'mkdocs-material-extensions',
                'mkdocstrings-python',
                'mknotebooks',
                'pixi-pycharm',
                'plotly',
                'pymdown-extensions',
                'pytest',
                'xlsxwriter']
    output = run(["pixi", "add"] + packages, check=True, capture_output=True, text=True)
    print(output.stdout)

if __name__ == '__main__':
    main()
