# Report: scenario ordering, per-scenario `s_nom`, and the phantom cycle

Status 2026-09-30. Branch `stoch-opt-grid`. PyPSA 1.3.0 (pixi env). Builds on
`ac_lines_kvl_check_HANDOFF.md` (findings F1-F6) and
`ac_lines_kvl_check.ipynb`. Companion notebook:
`notebooks/stochastic_grid_reactance_effects.ipynb` (runs clean end-to-end;
outputs not committed - no `nbconvert`/`nbclient` in env, same as A3).

## BLUF

In the **real pypsa-de stochastic-grid workflow the per-scenario reactance is
shared** (the CSV overrides touch only `s_nom`, never `x`). A stochastic program
is a sum over scenarios, so it is inherently order-invariant; the one thing that
makes PyPSA order-dependent is that it reads the KVL from `scenarios[0]` (F3) -
and with `x` identical across scenarios that dependence vanishes too. Therefore
**scenario ordering is a no-op and "strongest grid first" cannot help** (verified
on the real 49- and 27-cluster networks). The actual error is the **phantom
cycle**: a line at `s_nom=0` keeps its finite built-line reactance, which forces
its two endpoints to equal voltage angle - a spurious constraint that distorts
flows on every other line in its cycle. It is always pessimistic, worst for
double-circuit and meshed corridors, zero for radial ones, and not removable by
ordering. The clean fix is to model uncertain corridors as links (or drop them
from the KVL via `x=inf`), bounded by per-topology deterministic solves.

## Findings

- **F7 - Reactance is shared across scenarios in the real workflow.**
  `build_grid_scenario_csvs.py` writes only the capacity column; `_apply_overrides`
  sets `s_nom` + `s_nom_extendable=False` and never touches `x`. Verified on
  `base_s_49__none_2035_topology-stochastic.nc`: the 7 lines whose `s_nom` toggles
  0 <-> 5265 MW have `x` std = **0.0** across scenarios. (`n.lines["x"]` is the
  series reactance in Ohm - coordinates live on buses, `n.buses["x"]` = longitude.)
- **F8 - The stochastic solve is exactly order-invariant.** Structural argument
  (sum over scenarios; only `scenarios[0]` KVL weights break symmetry; those are
  identical here) + empirical: swapping scenario order gives bit-identical
  objective and flows in the toy (shared `x`) and identical expected backup on the
  real 27-cluster grid (3 orderings -> 14.87 TWh each). **"Strongest grid first"
  is a no-op.** (Rejects the mitigation hypothesis, extends handoff F6.)
- **F9 - The real error is the phantom cycle, not ordering.** `s_nom=0` with
  finite `x` implies flow 0 *and* `theta_i = theta_j`. The angle-lock is a
  spurious coupling the truly-absent grid (`x=inf`) would not have. The
  fully-built scenario is correct (its `x` is the real built-line value); the
  *zeroed* scenarios carry the error.
- **F10 - Three regimes (toy, deterministic, VOLL = 1000, load = 100 MW).**
  phantom vs correct cost: **bridge/radial = 100000 vs 100000 (error 0, safe)**;
  **mesh = 100000 vs 50050 (error ~50k)**; **parallel-reinforce = 100000 vs 100
  (error ~99.9k, catastrophic)**. The phantom bites only when the toggled line is
  in a cycle; zeroing one circuit of a double-circuit corridor angle-locks the
  parallel circuit to zero too.
- **F11 - Three-scenario toy (zero/mod/high), expected system cost.** Ground
  truth (per-topology, correct `x`) = **1750**. As-is stochastic (shared built
  `x`) = **3730**, order-invariant; the `zero` scenario is a phantom catastrophe
  (10000 vs correct 5050) and `mod` carries a reactance mismatch (1090 vs 100).
  "Strongest-first" scaled-x = 3730 (= as-is). Only **weakest-first (`x=inf`
  first) = 1750** recovers the correct cost, by dropping the line from the shared
  KVL (link-like). A large-but-finite `x` does **not** remove the phantom.
- **F12 - Real 27-cluster grid (flow-disruption proxy).** All 18 DE-internal AC
  corridors are in cycles; **every one shows phantom > correct** (always
  pessimistic). Zeroing one representative corridor (line 27) distorts
  energy-weighted mean flows by >50 MW on **12 of 18** other corridors (top:
  line 31 by ~2.3 GW, line 25 by ~1.8 GW). Six corridor *pairs* share identical
  phantom cost - these are double circuits (the catastrophic parallel regime).
- **F13 - Weakest-first works on the real grid.** 3-scenario stochastic on line
  27: as-is expected backup 14.87 TWh (any ordering); weakest-first (`x=inf`)
  2.32 TWh, with the `zero` scenario matching the ground-truth 6.68 TWh exactly.

## Risks / caveats

- **R1 - Magnitude metric is a proxy, not literal cost.** The saved
  sector-coupled networks (PyPSA 1.1.2) do **not** re-solve under 1.3.0 (verified:
  Gurobi and HiGHS both infeasible, even unmodified). So the real-grid numbers
  come from an AC-only redispatch testbed: real topology/reactances/capacities,
  nodal injections reconstructed from the stored flows, backup priced at
  1 EUR/MWh (free curtailment). This **forbids cheap generation redispatch and so
  overstates absolute cost** (per-corridor backup up to ~82 TWh vs 114.7 TWh DE
  throughput - clearly an upper bound). Faithful: the direction (always
  pessimistic), the flow-distortion pattern, order-invariance, and the mitigation
  ranking. The clean *cost* bound is the toy: up to ~100% of corridor value.
- **R2 - 1 h resolution does not fix it, likely worsens the worst case.** The
  phantom is a per-snapshot structural constraint, present every hour. Hourly
  resolution resolves the most congested hours (peak load + low wind) that
  representative snapshots smooth away, where redispatch headroom is tightest - so
  the relative error should grow, not shrink. With endogenous investment the
  pessimism biases the shared decision toward over-building, compounding in a
  myopic pathway. The 126-snapshot figure is a conservative estimate.

## Mitigations (ranked)

Key constraint: **PyPSA's native stochastic KVL cannot hold a per-scenario
reactance.** Writing the correct `x` per scenario is ignored except for
`scenarios[0]` (F3), so "scale `x` with `s_nom`" is **not deliverable** in one
stochastic solve.

- **O4 - Model uncertain corridors as `Link`s** (or put the absent topology first
  with `x=inf` so the corridor leaves the shared KVL). Removes the phantom
  entirely; trades a large *pessimistic* error for a small *optimistic* one (free
  routing). Best general-purpose fix when uncertain corridors are few. [recommended]
- **O2 - Per-topology deterministic solves** with correct `x`. Exact physics,
  loses stochastic coupling; use to bound/validate (the ground truth here).
- **O5 - Make only radial/bridge corridors uncertain.** Phantom error is exactly
  zero off-cycle (F10). Check each toggling line against `n.cycle_matrix()`.
- **O1 - Deliberate `scenarios[0]`** - meaningful only combined with O4 (choosing
  the `x=inf`/link topology first); with shared `x` it does nothing.
- **O6 - Accept and document** where corridors are few and meshed-but-not-
  parallel: bounded, always pessimistic (conservative on the weak-grid scenario).

**Recommended:** O4 + O2 - links inside the stochastic solve for the handful of
genuinely uncertain corridors, with per-topology deterministic solves alongside to
bound the idealisation error. Do not rely on ordering or on per-scenario `x`.

## Next actions

- **A5 - Decide the corridor treatment** in `build_grid_topology.py`: convert the
  toggling AC lines to links for the stochastic branch (O4), or restrict toggling
  to off-cycle corridors (O5). Document the choice and the residual bias.
- **A6 - Add per-topology deterministic bounds** to the comparison exports (O2) so
  every stochastic run ships with its correct-physics envelope.
- **A7 - Write the methods caveat** (phantom cycle + shared reactance + order
  no-op) into the paper/docs. Supersedes handoff A4 with concrete numbers.
- **A8 (optional) - Re-run the magnitude on a freshly-built 1.3.0 prenetwork** to
  replace the AC-only proxy (R1) with a full co-optimised cost, once a current
  prenetwork is available.

## Related files

- Notebook: `notebooks/stochastic_grid_reactance_effects.ipynb` (Sections 0-5 +
  Verdict).
- Prior work: `notebooks/ac_lines_kvl_check.ipynb`,
  `notebooks/ac_lines_kvl_check_HANDOFF.md`.
- Overrides: `scripts/pypsa-de/build_grid_scenario_csvs.py` (writes only `s_nom`),
  `scripts/pypsa-de/build_grid_topology.py` (`_apply_overrides`).
