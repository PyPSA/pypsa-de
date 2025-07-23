import pathlib


def py_to_md_file(subdir):
    """Generate all required markdown files in a subdirectory."""
    source = SRC_DIR / subdir
    package = SRC_DIR.stem
    target = ROOTDIR / "docs-at" / "reference" / package / subdir

    for p in source.iterdir():
        if p.suffix != ".py":
            continue
        if p.stem == "__init__" and not subdir:  # special index file
            continue

        if p.stem == "__init__":
            new = target / "index.md"
            module = f"{package}.{subdir}" if subdir else package
        else:
            new = target / (p.stem + ".md")
            module = f"{package}.{subdir}.{p.stem}" if subdir else f"{package}.{p.stem}"

        target.mkdir(exist_ok=True)
        new.touch(exist_ok=True)
        new.write_text(f"::: {module}")  # overwrites all contents


if __name__ == "__main__":
    ROOTDIR = pathlib.Path(".").resolve()  # assuming CDW pypsa-at
    assert ROOTDIR.stem == "pypsa-at", "Must run from repo root 'pypsa-at'"
    SRC_DIR = ROOTDIR / "evals"
    py_to_md_file("")
    py_to_md_file("plots")
    py_to_md_file("views")
    py_to_md_file("data")
