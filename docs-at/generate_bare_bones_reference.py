# -*- coding: utf-8 -*-
import pathlib
from importlib import resources


def py_to_md_file(subdir):
    """"""
    directory = SRC_DIR / subdir
    target = pathlib.Path(r"C:\Storage\esmtools\docs\reference") / subdir

    for p in directory.iterdir():
        if p.suffix != ".py":
            continue
        if p.stem == "__init__" and not subdir:  # special index file
            continue

        if p.stem == "__init__":
            new = target / "index.md"
            module = f"esmtools.{subdir}" if subdir else "esmtools"
        else:
            new = target / (p.stem + ".md")
            module = f"esmtools.{subdir}.{p.stem}" if subdir else f"esmtools.{p.stem}"
        # print(str(new), f"::: {module}")
        target.mkdir(exist_ok=True)
        new.touch(exist_ok=True)
        new.write_text(f"::: {module}")


if __name__ == "__main__":
    SRC_DIR = resources.files("pypsa-at") / "evals"

    py_to_md_file("")
    py_to_md_file("plots")
    py_to_md_file("views")
    py_to_md_file("data")
