import pathlib


def py_to_md_file(subdir):
    """Generate all required markdown files in a subdirectory."""
    directory = SRC_DIR / subdir
    target = ROOTDIR / "docs-at" / "reference" / "evals" / subdir

    for p in directory.iterdir():
        if p.suffix != ".py":
            continue
        if p.stem == "__init__" and not subdir:  # special index file
            continue

        if p.stem == "__init__":
            new = target / "index.md"
            module = f"evals.{subdir}" if subdir else "evals"
        else:
            new = target / (p.stem + ".md")
            module = f"evals.{subdir}.{p.stem}" if subdir else f"evals.{p.stem}"
        # print(str(new), f"::: {module}")
        target.mkdir(exist_ok=True)
        new.touch(exist_ok=True)
        new.write_text(f"::: {module}")


if __name__ == "__main__":
    ROOTDIR = pathlib.Path(".").resolve()  # assuming CDW pypsa-at
    assert ROOTDIR.stem == "pypsa-at", "Must run from repo root 'pypsa-at'"
    SRC_DIR = ROOTDIR / "evals"
    py_to_md_file("")
    py_to_md_file("plots")
    py_to_md_file("views")
    py_to_md_file("data")
