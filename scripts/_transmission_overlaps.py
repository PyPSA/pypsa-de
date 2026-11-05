# SPDX-FileCopyrightText: Contributors to PyPSA-DE
#
# SPDX-License-Identifier: MIT

"""
Reconcile competing plan descriptions on matched model branches before staging.

A collision is a modelling conflict, never proof of physical project identity.
Reviewed shared scope uses configured include order only on the common branch.
Unreviewed cross-plan collisions use the same conservative provisional choice,
but are always queued for scope review. Explicit successive work is exempt.
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def project_key(row):
    ident = str(row.get("project_id", ""))
    if ident.startswith("TYNDP"):
        return "TYNDP:" + ident.rsplit("_", 1)[-1]
    try:
        tags = json.loads(row.get("tags") or "{}")
    except (ValueError, TypeError):
        tags = {}
    project = re.search(r"P\d+[a-z]?", str(tags.get("project_id", ident)))
    measure = tags.get("measure_id") or (
        re.search(r"_?(M\d+[\w]*)$", ident).group(1)
        if re.search(r"_?(M\d+[\w]*)$", ident)
        else ""
    )
    return f"NEP:{project[0]}:{measure}" if project and measure else ident


def reconcile(branches, adjustments, decisions=None):
    """
    Return reconciled branches, refreshed retirements and a complete audit.

    Include order is the original concatenation order, not commissioning order.
    Unshared sections of a chained project remain. No source CSV is changed.
    """
    columns = [
        "base_asset",
        "kept_project",
        "other_project",
        "kept_plan",
        "other_plan",
        "kept_key",
        "other_key",
        "kept_year",
        "other_year",
        "action",
        "review_status",
        "decision_id",
        "basis",
        "source_file",
    ]
    if (
        branches.empty
        or "replaced_asset" not in branches
        or "source_file" not in branches
    ):
        return branches, adjustments, pd.DataFrame(columns=columns)
    decisions = pd.DataFrame() if decisions is None else decisions
    lookup = {}
    for row in decisions.to_dict("records"):
        pair = frozenset((row["key_a"], row["key_b"]))
        if pair in lookup:
            raise ValueError("Duplicate overlap decision: " + str(pair))
        if row["policy"] not in ("shared_branch_precedence", "separate_successive"):
            raise ValueError("Unknown overlap policy: " + row["policy"])
        if not row.get("evidence") or not row.get("basis"):
            raise ValueError("Overlap decision needs local evidence and basis")
        lookup[pair] = row
    keep = np.ones(len(branches), dtype=bool)
    events = []
    ups = branches[branches.replaced_asset.notna()].copy()
    ups["_position"] = np.flatnonzero(branches.replaced_asset.notna())
    ups["_plan"] = ups.source_file.map(lambda p: Path(str(p)).parent.name)
    for asset, group in ups.groupby("replaced_asset", sort=False):
        # All stages from the authoritative plan remain; within-plan scope is
        # not adjudicated by a cross-dataset reconciliation rule.
        first = group.iloc[0]
        authority = group[group["_plan"].eq(first["_plan"])]
        for _, other in group.iterrows():
            if first["_plan"] == other["_plan"]:
                continue
            b = project_key(other)
            # Prefer the reviewed project relation. Otherwise attribute the
            # conflict to the closest dated stage of the authoritative plan,
            # rather than always naming its first (possibly unrelated) stage.
            reviewed = authority[
                authority.apply(
                    lambda candidate: frozenset((project_key(candidate), b)) in lookup,
                    axis=1,
                )
            ]
            candidates = reviewed if not reviewed.empty else authority
            order = (
                (candidates.build_year - float(other.build_year)).abs().fillna(np.inf)
            )
            primary = candidates.iloc[int(np.argmin(order.to_numpy()))]
            a = project_key(primary)
            decision = lookup.get(frozenset((a, b)), {})
            separate = decision.get("policy") == "separate_successive"
            if not separate:
                keep[int(other._position)] = False
            events.append(
                dict(
                    base_asset=asset,
                    kept_project=primary.project_id,
                    other_project=other.project_id,
                    kept_plan=primary._plan,
                    other_plan=other._plan,
                    kept_key=a,
                    other_key=b,
                    kept_year=primary.build_year,
                    other_year=other.build_year,
                    action="separate_successive"
                    if separate
                    else "shared_branch_precedence",
                    review_status="reviewed_model_approximation"
                    if decision and not separate
                    else "reviewed_distinct_work"
                    if separate
                    else "unresolved_scope",
                    decision_id=decision.get("decision_id", ""),
                    basis=decision.get(
                        "basis",
                        "No physical identity inferred: first included plan retained provisionally on shared model branch",
                    ),
                    source_file=str(other.source_file),
                )
            )
    result = branches.iloc[keep].copy()
    # Retirement was originally computed before plan reconciliation. Rebuild
    # it from retained rows so a suppressed earlier date cannot retire the base.
    active = result[result.replaced_asset.notna()]
    if not active.empty:
        years = active.groupby("replaced_asset").build_year.min()
        adjustments = adjustments.groupby(level=0).first().loc[years.index].copy()
        adjustments["build_year"] = 1990
        adjustments["lifetime"] = years - 1990
    else:
        adjustments = adjustments.iloc[:0].copy()
    return result, adjustments, pd.DataFrame(events, columns=columns)
