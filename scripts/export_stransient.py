"""
Utility script to extract network components from PyPSA-DE exports and format
them into STRANSIENT-compatible CSV files.

This can be run standalone via the CLI, or automatically as part of a Snakemake workflow.
"""

import argparse
from pathlib import Path

import pandas as pd


def export_stransient(base: Path, out: Path):
    """
    Parses PyPSA network CSV exports and translates them into the STRANSIENT format.

    Parameters
    ----------
    base : Path
        Input directory containing standard PyPSA exports (buses.csv, lines.csv, generators.csv, loads.csv).
    out : Path
        Output directory to save the formatted STRANSIENT CSV files.
    """
    out.mkdir(parents=True, exist_ok=True)
    bus_df = pd.read_csv(base / "buses.csv").rename(
        columns={"name": "bus_id", "v_nom": "vn_kv"}
    )
    bus_df["vm_pu"] = bus_df["v_mag_pu_set"]
    strans_bus = bus_df[["bus_id", "vn_kv", "type", "vm_pu"]].copy()
    strans_bus["area"] = "DE"
    strans_bus.to_csv(out / "stransient_bus.csv", index=False)
    lines = pd.read_csv(base / "lines.csv").rename(columns={"name": "branch_id"})
    lines["type"] = "AC_line"
    lines[
        ["branch_id", "bus0", "bus1", "r_pu", "x_pu", "length", "i_nom", "type"]
    ].to_csv(out / "stransient_branch.csv", index=False)
    gens = pd.read_csv(base / "generators.csv")
    q_source = (
        "q_nom"
        if "q_nom" in gens.columns
        else "q_set"
        if "q_set" in gens.columns
        else None
    )
    rename_map = {"name": "gen_id", "p_nom": "p_max_mw", "carrier": "type"}
    if q_source is not None:
        rename_map[q_source] = "q_max_mvar"
    else:
        gens["q_max_mvar"] = 0.0
    gens = gens.rename(columns=rename_map)
    required_cols = ["gen_id", "bus", "p_max_mw", "q_max_mvar", "type"]
    gens[required_cols].to_csv(out / "stransient_gen.csv", index=False)
    loads = pd.read_csv(base / "loads.csv")
    p_source = (
        "p_mw"
        if "p_mw" in loads.columns
        else "p_set"
        if "p_set" in loads.columns
        else None
    )
    q_source = (
        "q_mvar"
        if "q_mvar" in loads.columns
        else "q_set"
        if "q_set" in loads.columns
        else None
    )
    load_rename = {"name": "load_id"}
    if p_source:
        load_rename[p_source] = "p_mw"
    else:
        loads["p_mw"] = 0.0
    if q_source:
        load_rename[q_source] = "q_mvar"
    else:
        loads["q_mvar"] = 0.0
    loads = loads.rename(columns=load_rename)
    loads[["load_id", "bus", "p_mw", "q_mvar"]].to_csv(
        out / "stransient_load.csv", index=False
    )
    print("wired exports done")


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake = globals()["snakemake"]
        base_dir = Path(snakemake.input.exports_dir)
        out_dir = Path(snakemake.output.stransient_dir)
        export_stransient(base_dir, out_dir)
    else:
        parser = argparse.ArgumentParser(
            description="Export STRANSIENT grids from PyPSA"
        )
        parser.add_argument(
            "--exports-dir",
            type=Path,
            default=Path(
                "results/20260114_limit_cross_border_flows/KN2045_Mix/exports"
            ),
            help="Directory containing PyPSA export CSVs",
        )
        parser.add_argument(
            "--out-dir",
            type=Path,
            default=None,
            help="Output directory. Defaults to <exports_dir>/../stransient",
        )
        args = parser.parse_args()

        base_dir = args.exports_dir
        out_dir = args.out_dir if args.out_dir else base_dir.parent / "stransient"
        export_stransient(base_dir, out_dir)
