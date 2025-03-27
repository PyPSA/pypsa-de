---
hide:
  - toc
---

|python|3.12|
R|codestyle|black|.dark-badge|
R|linter|ruff|.tools-badge|
R|package manager|pixi|.tools-badge|

[black]: https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html
[ruff]: https://docs.astral.sh/ruff/linter/
[pixi]: https://pixi.sh/latest/

# pypsa-at Documentation

This site hosts the documentation for `esmtools`, a Python module
designed to evaluate energy system model results based on
[PyPSA](https://pypsa.readthedocs.io/en/latest/).
The goal of this project is to provide a framework for
analyzing and visualizing energy model outputs,
using Python. The `esmtools` pacakge heavily relies on
[pypsa.statistics](https://pypsa.readthedocs.io/en/latest/api/statistics.html)
to extract views from solved network files.
The plotting library used for visualizations is [Plotly](https://plotly.com/python/)
and results may be exported in other formats such as Excel, CSV or JSON files.


## Table Of Contents

The documentation follows the best practice for
project documentation as described
in the [Diátaxis documentation framework](https://diataxis.fr/)
and consists of four separate parts:

1. [How-To Guides](how-to-guides/index.md)
2. [Tutorials](tutorials/index.md)
3. [Reference](reference/index.md)
4. [Explanation](explanations/index.md)

Quickly find what you're looking for depending on
your use case by looking at the different pages.

This documentation is built using [MkDocs](https://www.mkdocs.org/),
[mkdocstrings](https://mkdocstrings.github.io/python/), and the
[Material for MkDocs theme](https://squidfunk.github.io/mkdocs-material/getting-started/).
