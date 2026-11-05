# NEP panel data for the grid-delay scenario weights

Copy of the `nep_panel` outputs (v0.13.0, 2026-11-04; repository `general_agent_jobs/grid delay/nep_panel`).
Definitions: `nep_panel/docs/data.md`. Analysis: `nep_panel/reports/2026-11-04_probability_inputs/report.md`.

| File | Written by | Content |
|---|---|---|
| `panel.csv` | nep_panel `build.py` | one row per line and report (NEP 2012-2025, implementation reports 2018-2022): status, target year |
| `lines.csv` | nep_panel `build.py` | one row per line: outcome, commissioning date and precision |
| `delay_table.csv`, `cancellation.csv` | nep_panel `analysis.py` | Tables 1 and 2, used by the explorer's self-check |
| `delay_pmf.csv` | `notebooks/nep_panel_delay_explorer.ipynb` | delay distribution per status: status, delay (years, `inf` = beyond all observed), probability |
| `cancel_rate.csv` | same | probability of cancellation before commissioning per status, with counts |
| `delay_correlation.csv` | same | delay correlation: level (period, project), rho, 90 % interval (lo, hi), groups |

`notebooks/build_probabilities_from_nep_panel.ipynb` reads the last three. `data/pypsa-de/` is allowlisted in `.gitignore`,
so these files can be versioned.
