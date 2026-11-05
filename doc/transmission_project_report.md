# Exogenous transmission project audit

`scripts/report_transmission_projects.py` creates four DIN A4 PDF reports for
any forward year difference, independently of Snakemake and optimization:

- German internal AC lines (MVA).
- AC interconnectors touching Germany (MVA).
- German internal DC links (MW).
- DC interconnectors touching Germany (MW).

By default, maps show NUTS2 boundaries and connect representative points of endpoint
regions. Every displayed branch has its own colored number, linked to a
project legend with net capacity change, original CSV ID, source plan,
commissioning year, status, replacement details, and available TYNDP project /
investment IDs or NEP locations. Parallel projects have separate curves.
Projects within one displayed region are omitted from maps, legends,
corridor tables, and displayed totals. With a supplied busmap, branches known
to collapse within a model bus are also omitted at every display resolution.
All accepted projects remain in the physical audit ledger. Offshore endpoints outside NUTS2 retain their coordinates.

Earlier reports used circles for endpoints within a region. These projects
are now hidden. For example, `TTG-P67M102M103aM103b` connects Simbach and
Altheim inside NUTS2 region DE22 and is omitted from the NUTS2 display.
A shared display region implies model removal only when that region matches
the actual model clustering; an optional busmap identifies model collapse
independently of the chosen map resolution.
Its tags say `Ersatzneubau`, but the source lists it in `new_lines.csv`, so
the workflow/report treats it as a new line with zero associated retirement.
Construction-type tags alone do not establish a matched replacement asset.

Each legend entry now shows the source plan (NEP, TYNDP year, or other
dataset), the source project ID, and the supplied project name. Where no
official name is supplied, available endpoint names are explicitly marked
`(endpoints)`; absent names are marked `Name not supplied`. NEP measure IDs,
construction type, and TYNDP investment IDs are shown when supplied. Existing
Project labels show only numbers; each project connection and its legend
entry share a distinct color. Legends include the transmission-project folder
name (for example, `nep2023confirmed`), and fields are separated with `|` and
wrapped only when the available row width is exhausted. URLs and hyperlink
annotations are omitted from PDFs. Source URLs remain available in the ledger,
alongside `source_plan`, `reference_id`, and `project_name`.

Endpoint names in the legend follow physical `bus0` / `bus1` and their region
order. Source geometry is compared with physical bus coordinates before
region snapping, so upgrades matched in reverse are labelled consistently.
The ledger records `endpoint_order`, aligned `endpoint0_name` /
`endpoint1_name`, and original `source_endpoint0_name` / `source_endpoint1_name`;
raw source tags are preserved. Original `location0/1` tags are omitted from
the legend to avoid mixing the two orders. If orientation cannot be resolved,
the report uses bus IDs and explicitly marks source endpoint order unresolved.

There are exactly four PDF **files per interval**. If a legend needs more room, its PDF has
additional A4 sheets containing only the remaining legend entries. The map
appears on the first page only. Compact entries wrap to the full available
page width, and continuation pages use the full height below the header.
This avoids shrinking or omitting project information. Blank categories still
produce an A4 report. Page dimensions remain 210 × 297 mm.
At the end of each legend, a compact omitted-project list shows
`project_id (count)` for changed projects hidden within a displayed region
or known model bus. Counts refer to omitted physical branch rows for that
project in the report interval and category. Empty lists are omitted.

## Standalone usage

Run from the repository root in the existing project environment. Point all
inputs at the **same run**; `base.nc` must be the physical network before
projects are attached, rather than `base_extended.nc` or a clustered network.

```bash
run_resources=resources/20261027-NUTS2-365h-2025EXO-2030EXO-2035EXO-10-50-40-new-plots/KN2045_Bal_v5_abroad_capa_freeze

pixi run python scripts/report_transmission_projects.py \
  --base-network "$run_resources/networks/base.nc" \
  --projects "$run_resources/transmission_projects" \
  --nuts2 data/nuts/NUTS_RG_01M_2021_4326_LEVL_2.geojson \
  --source-dir data/transmission_projects/nep \
  --source-dir data/transmission_projects/manual \
  --source-dir data/transmission_projects/tyndp2020 \
  --output-dir results/transmission_project_audit/nuts2
```

Use `--start-year START --end-year END` for a custom interval; the years need not be planning horizons.
Omitting both year flags exports **2025 until 2030** and **2030 until 2035**,
with four separate PDFs and their audit tables per interval. Default directories
are `results/transmission_project_audit/RESOLUTION/2025_until_2030/` and
`results/transmission_project_audit/RESOLUTION/2030_until_2035/`. PDF names follow
`AC_internal_2025_until_2030.pdf`. With default intervals, `--output-dir` selects
the parent folder; each interval gets its own subfolder. Supply both year flags
for a custom single interval; its explicit output directory is used directly.
All generated outputs belong under the git-ignored `results/` directory.
For example, `--output-dir results/projects_including_planned_not_yet_permitted_under_consideration`
keeps a standalone analysis outside the scenario results tree. Include the
actual CSV statuses `planned_not_yet_permitting` and `under_consideration` in
the configuration before rebuilding the processed inputs; the output name
does not control filtering.
Use `--country` to change the default `DE` scope. Foreign-to-foreign projects
remain in the ledger but are excluded from the four report totals.

This reads the **processed** five project CSVs (`new_lines`, `new_links`,
`adjust_lines`, `adjust_links`, `new_buses`). They contain the actual matched,
filtered projects and nominal ratings the attachment step consumes. It does
not run project preprocessing. After changing a raw project CSV or project
configuration, rebuild these inputs first; reporting old processed files will
describe the old run. Raw source directories supply provenance and a source
inventory, and do not override processed capacities.

Each legend entry has a clickable source link in a right-hand column. With a
`source_pdf_page` tag it reads `p. 186 ↗` and opens that PDF page
(`#page=186`, the PDF page index, not the printed page number); otherwise it
reads `source ↗`. Whether a viewer honours the page jump depends on the viewer
(browser viewers and Adobe do). To add or correct links without editing the
project CSVs, pass `--source-links` with a CSV of `dataset` (source folder
name), `project_id`, `source_url` and `source_pdf_page`.

## Plot resolution

Choose `--resolution nuts2` (default), `--resolution nuts3`, or
`--resolution grid`. NUTS modes connect region representative points; grid
mode connects the original matched bus coordinates and shows the physical
base-network branches. Grid mode draws endpoint connections, not surveyed
line routes. Physical accounting stays unchanged in the ledger. Displayed totals can
differ by resolution because within-region branches are excluded.

The script selects the local NUTS2021 shapes automatically. Override them
with `--regions PATH` (the old `--nuts2` flag remains an alias). Supply NUTS3
polygons for NUTS3 mode, and NUTS2 polygons for NUTS2 or the grid backdrop.
Administrative regions retain their NUTS codes and names as muted labels.
NUTS2 and grid maps label the home-country regions; NUTS3 maps label project
endpoint regions to limit crowding. Foreign endpoint regions are labelled
when available. Project connections keep their separate colored number labels.

For example, use the already built upgrades-enabled project inputs:

```bash
pixi run python scripts/report_transmission_projects.py \
  --base-network "$run_resources/networks/base.nc" \
  --projects results/transmission_project_audit/with_upgrades/processed_projects \
  --resolution nuts3
```

Replace `nuts3` with `grid` for physical bus coordinates. Default outputs are
separated by resolution under `results/transmission_project_audit/RESOLUTION/START_until_END/`.
The ledger includes `region0`, `region1`, `within_region`, and `shown`, plus
`nuts2_0/1` or `nuts3_0/1` for the chosen administrative mode. Corridor files
are named `nuts2_corridors.csv`, `nuts3_corridors.csv`, or `grid_corridors.csv`.
Grid corridor endpoints are physical bus IDs. `summary.json` records the
resolution and shape file. Existing explicit output directories still work;
use separate directories when comparing resolutions.

## Optional workflow target

Run `pixi run snakemake report_transmission_projects_default --cores 1`
to export both default intervals for all three resolutions as separate files.

`report_transmission_projects` in `rules/pypsa-de/stochastic_grid.smk` uses the
current run's base network, processed CSVs, enabled source directories, and
NUTS2021 shapes. Request a concrete output, for example:

```bash
pixi run snakemake \
  results/<run-prefix>/<run-name>/transmission_project_audit/nuts2/2025_until_2030/summary.json \
  --cores 1 --configfile config/config.de.yaml
```

Replace the run path with the path produced by your configuration, and use
`nuts3` or `grid` in place of `nuts2` to select another resolution. The target
builds only the required preprocessing inputs and report. It does not require
solved networks or enabling stochastic grid scenarios, and is not added to
the default analysis target. The audit describes project attachment inputs;
if `transmission_projects.enable` is false, those inputs are not attached to
the model even if explicitly built for this report.

## Capacity accounting and exports

The raw `tyndp2024draft/new_links.csv` row `TYNDP2024_660` is corrected to
match the existing NEP DC2 modelling coordinates for Osterath–Philippsburg.
Normal duplicate detection then retains NEP DC2 once. These reuse existing
modelling coordinates, not independently surveyed positions. The original row
and correction rationale are preserved under
`results/transmission_project_audit/with_upgrades/data_corrections/`.
There is no preprocessing coordinate override.

Offshore bus creation uses distinct endpoint IDs and reuses previously
created nearby buses, preventing duplicate bus exports and accidental
connections of a project's two different offshore ends to the same bus.
Upgrade matching uses the original base assets, rather than proposed lines
temporarily added for duplicate detection. Blank line types use the same
default as missing line types.

Capacity at year `Y` is the processed nominal rating when
`build_year <= Y < build_year + lifetime`, and zero otherwise. Absent project
lifetimes are treated as unlimited. Start-year commissioning is already in
the baseline; end-year commissioning is included in the change.

For upgrades, the original branch comes from `base.nc` with its retirement
date updated using `adjust_lines.csv` / `adjust_links.csv`. The report records
replacement rating and retired rating separately. For example, a 1,000 MW
link replaced by a 2,000 MW link contributes +2,000 MW gross, 1,000 MW retired,
and +1,000 MW net. Missing retirement references and duplicate IDs fail rather
than silently overstate capacity. A processed link with `p_nom: 0` contributes
zero added capacity, including when `new_link_capacity: zero` is configured.

Several upgrades of one existing branch are staged by
`build_transmission_projects.py`: `<branch>_upgraded`, `<branch>_upgraded2`, ...
in build-year order. Each stage replaces the previous one (its `replaced_asset`)
and ends with a finite `lifetime` in the next stage's build year. Only the
original branch has a retirement adjustment. A stage's own row records its
end of life, so the report counts no extra retirement for a staged
replacement. Of two upgrades in the same build year, the first in include
order is kept and the other is logged.

The reported difference is the endpoint capacity change, not the sum of all
commissioning events that might occur and retire inside an interval.

Projects without a usable commissioning year remain in the ledger with
`commissioning_known: false` and unknown (blank/NaN) interval capacities.
They are excluded from dated totals rather than assigned an assumed year.
The build/report logs identify these projects, and `summary.json` counts
undated projects touching the selected country in each category.

- `project_ledger.csv`: every accepted project, including unchanged and
  foreign projects; physical network ID; original source ID/file; provenance
  confidence; tags; endpoints and countries; NUTS2 mapping; nominal, start,
  end, gross, retired, and net capacities; map legend key.
- `<resolution>_corridors.csv`: summed changes per unordered endpoint pair and carrier /
  scope, with contributing project IDs. Includes only displayed projects.
- `source_audit.csv`: inventory of the supplied raw source directories and
  whether each source project could be identified in processed output.
  Absence alone cannot distinguish status filtering, configured skipping,
  geographic exclusion, duplicate detection, or failed matching: consult the
  corresponding `build_transmission_projects.log` and run configuration.
- `summary.json`: displayed category totals, excluded project counts and net
  ratings, years, units, and input paths. AC and DC
  totals remain separate because MVA and MW are different quantities.

The builder now preserves `project_id`, `source_file`, `project_action`, and
`replaced_asset` through upgrade renaming. For older artifacts, new projects
are identified by source CSV row ID; upgrade identity is recovered only from
a unique action / tags / build-year match. Unresolved or ambiguous provenance
is explicitly marked. The complete tags remain in the ledger, including URLs.

AC values are nominal apparent power, without `s_max_pu` derating or DLR.
DC values are nominal active power, without conversion losses. Neither is a
power flow or import limit, and sums across consecutive branches are not
end-to-end transfer capability. Lengths are exported but totals are not MW-km.

## Trace into a clustered model

NUTS2 display mapping is geographic and does not establish final model
membership. To trace bus aggregation, optionally supply `--busmap` with a
**composed physical-to-final-model busmap** (index column + one value column).
The ledger then includes `model_bus0`, `model_bus1`, and `model_mapping`, which
distinguishes between-model-bus, collapsed-within-bus, and unmapped projects.
For this repository's two-stage mapping, compose the simplification and final
clustering maps before passing them:

```python
import pandas as pd

simplification = pd.read_csv("<resources>/busmap_base_s.csv", index_col=0, dtype=str).squeeze("columns")
clustering = pd.read_csv("<resources>/busmap_base_s_adm.csv", index_col=0, dtype=str).squeeze("columns")
simplification.map(clustering).to_csv("/tmp/physical_to_model.csv")
```

The physical ledger includes within-model-bus projects for auditing, but
plots and displayed totals exclude known collapsed projects. Unmapped
projects remain visible and explicitly marked as unmapped in the ledger. Later voltage
scaling, clustering, line removal, endogenous investment, and stochastic
scenario overrides can change the final topology. In particular, the
`*_exogen` scenario generator currently gates nominal capacities by build year
and does not apply lifetime retirement in `branch_table`; this report's net
replacement accounting can therefore differ from a scenario override if the
pristine network still contains both the old and replacement branches. This
report exposes the attachment-stage inputs and retirement assumptions; it
does not claim to reconcile a solved or overridden model network.

The default configuration currently skips `upgraded_lines` and
`upgraded_links`. The source inventory helps make that absence visible; to
include upgrades, change `transmission_projects.skip` and rebuild inputs.

The upgrades-enabled example is saved under
`results/transmission_project_audit/with_upgrades/`: `effective_config.yaml`
records `transmission_projects.skip: []`, `processed_projects/` holds freshly
built inputs, the build log records matching/exclusions, and `2030-2025/`
contains the four PDFs and ledgers. It uses the current raw CSVs and leaves
the main configuration unchanged. Repeat it from the repository root with:

```bash
PYTHONPATH=. pixi run python results/transmission_project_audit/with_upgrades/run.py
```

The example runner now respects the current `transmission_projects` settings
instead of forcing upgrades on. Inspect `effective_config.yaml`, the source
audit, and build log for the exact datasets, accepted projects and exclusions.
`input_checksums.json` records the input files used for the latest refresh.

The latest exports for both default intervals and all three resolutions are
under `results/transmission_project_audit/with_upgrades/<resolution>/`, in
`2025_until_2030/` and `2030_until_2035/`. Earlier report folders are retained.
