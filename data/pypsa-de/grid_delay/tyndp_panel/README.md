# TYNDP panel data for the grid-delay correlation comparison

Copy of `Results/panel_long.csv` from the `transmission-timeline` repository
(github.com/damianoalessi/transmission-timeline, written by `scripts/export_panel.py`, local copy of 2026-11-03;
file not yet committed upstream). MIT License, subject to the terms of the original ENTSO-E publications.

| File | Content |
|---|---|
| `panel_long.csv` | one row per TYNDP investment and report (2010-2026): status (1 consideration, 2 planned not permitting, 3 permitting, 4 construction, 5 commissioned, 6 cancelled), expected commissioning year `Comm_Year`, `Present` (0 = imputed), project, countries, `DE_scope`, `Wind_Connection` |

Read by `notebooks/nep_panel_delay_explorer.ipynb` (Table 8) to compare delay correlations with the NEP panel.
Not used by the probability builder.
