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

This site hosts the documentation for `pypsa-at`


first, build the scenarios using the public DB access
``` sh
snakemake build_scenarios --configfile=config/config.public.yaml -f
```

Start the analysis:
``` sh
snakemake ariadne_all --configfile=config/config.public.yaml -f
```

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
