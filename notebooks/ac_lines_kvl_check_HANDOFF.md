# Handoff: AC lines, s_nom=0 and the stochastic KVL

Status as of 2026-09-24. Branch `stoch-opt-grid`. Verified against PyPSA 1.3.0
(pixi env). Deliverable: `notebooks/ac_lines_kvl_check.ipynb` (runs clean
end-to-end; outputs not saved - see A3).

## Context

Stochastic grid-scenario workflow feeds different topologies to one stochastic
network by overriding `s_nom` per scenario (see `scripts/pypsa-de/
build_grid_topology.py` `_apply_overrides`, which sets `s_nom` and
`s_nom_extendable=False`). Question that started this: is it risky to set an AC
line to `s_nom=0` in one run and raise it in another?

## Findings (all verified in code and/or the notebook)

- **F1 - `s_nom=0` does NOT drop a line from the KVL.** Cycle basis is built
  from topology (bus0/bus1), dropped only for infinite impedance, not
  capacity. `pypsa/network/power_flow.py:709-721` ("Cycles with infinite
  impedance are skipped"); active mask ignores `s_nom`
  (`pypsa/components/descriptors.py:151-169`). Notebook Section 1-2.
- **F2 - Stochastic model is structurally identical across scenarios.** Cycle
  basis taken from `scenarios[0]` only and shared
  (`pypsa/network/power_flow.py:711-713`). Notebook Section 3.
- **F3 - KVL *weights* (reactances) also come from `scenarios[0]`**, not just
  the structure (`pypsa/networks.py:1353-1363`, `cycle_matrix`). Proved by
  giving one line different `x` per scenario and swapping order - coefficient
  always follows the first scenario. Notebook Section 4.
- **F4 - Phantom cycle.** A zeroed line keeps its (finite-x) cycle, so its KVL
  constraint still couples the remaining lines - a spurious constraint the
  truly-absent grid would not have. `s_max_pu=0` behaves identically to
  `s_nom=0` here (NOT safer, contrary to first hypothesis). Only `x=inf`
  removes it, which is unavailable per-scenario (F2/F3). Notebook Section 1.
- **F5 - Reactance-not-scaling can be large.** Raising `s_nom` without cutting
  `x` makes flow refuse the reinforced corridor. Toy 3-bus loop: naive
  reinforcement captures 0% of the value; error ~2475 EUR (~96% of base cost).
  Always pessimistic in direction. Deterministic isolation. Notebook Section 5.
- **F6 - Scenario ordering is NOT a fix.** Putting the most-expanded grid first
  corrects that scenario but distorts the others (small corridors inherit
  too-low reactance -> spurious bottleneck); expected cost can worsen. Error is
  relocated, not removed. Notebook Section 6. This answers the "expanded-first"
  mitigation idea: rejected.

## Options for correct-ish physics (none is free)

- **O1 - Deliberate `scenarios[0]`.** Documented lesser-evil: pick the topology
  whose reactances matter most; accept bounded error elsewhere.
- **O2 - Separate per-topology deterministic solves.** Correct reactances per
  grid; loses the joint stochastic coupling. Use to validate/bound the
  stochastic run.
- **O3 - Links for switchable corridors.** No reactance/KVL at all; correct
  ONLY where the device is genuinely controllable (HVDC/PST). For passive AC
  lines it over-idealises (free routing) and corrupts the built-scenario
  physics - not a general phantom-cycle fix.

## Open questions / next actions

- **Q1 - Which corridors in the real workflow actually toggle across
  scenarios, and are any radial/bridge-like?** That is where F4/F5 bite
  hardest. A2 would quantify.
- **A1 - Decide `scenarios[0]` policy** in `build_grid_topology.py` (stochastic
  branch) and document it. Currently order is whatever the scenarios dict
  yields - make it explicit.
- **A2 - Port the toy sensitivity (Sections 5-6) to the real 49-cluster
  network:** compare stochastic flows/cost vs per-topology deterministic solves
  for the toggled corridors; report the magnitude for the actual grid.
- **A3 - `pixi add nbconvert`** so the notebook can be committed with tables +
  figure pre-rendered (nbconvert not in env; I validated by exec, outputs are
  empty on disk).
- **A4 - Write the methods caveat** (phantom cycle + scenario-0 reactances)
  into the paper/docs once A2 gives numbers.

## Related files

- Notebook: `notebooks/ac_lines_kvl_check.ipynb` (Sections 1-6 + Verdict).
- Scratchpad notes (session-local, copy out if wanted): `kvl_stochastic_answer.md`,
  `phantom_cycle_explained.md`, `reactance_not_scaling.md`.
- Overrides: `scripts/pypsa-de/build_grid_topology.py`,
  `scripts/pypsa-de/build_grid_scenario_csvs.py`.
- Iteration skip (non-extendable lines): `scripts/solve_network.py:1592`.
