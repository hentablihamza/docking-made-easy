"""Check the committed EGFR example outputs against each other.

Run in CI, where neither AutoDock Vina nor the network is available. It cannot
rerun the docking, but it can prove that the numbers quoted in the README and
in report.json are the numbers in the CSV files, and that the committed box is
exactly the box the code computes from the real 1M17 ligand coordinates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "examples" / "egfr_1m17"
FIXTURE = ROOT / "tests" / "fixtures" / "mini_complex.pdb"

sys.path.insert(0, str(ROOT / "src"))

from docking_made_easy import Structure, box_from_ligand  # noqa: E402

failures: list[str] = []


def check(condition: bool, message: str) -> None:
    print(("  ok   " if condition else "  FAIL ") + message)
    if not condition:
        failures.append(message)


def main() -> int:
    report = json.loads((DATA / "report.json").read_text())
    redock = pd.read_csv(DATA / "redock.csv")
    poses = pd.read_csv(DATA / "poses.csv")
    summary = pd.read_csv(DATA / "summary.csv")

    print("Box")
    ligand = Structure.from_file(FIXTURE).select_ligand("AQ4")
    box = box_from_ligand(ligand, padding=4.0)
    for key, value in box.to_dict().items():
        if key == "source":
            continue
        check(
            abs(value - report["box"][key]) < 1e-6,
            f"{key} recomputed from the crystal ligand = {value:.3f}",
        )

    print("\nRedocking")
    check(len(redock) == 9, f"{len(redock)} poses returned")
    check(
        abs(redock.iloc[0]["rmsd_to_crystal"] - report["redock"]["top_pose_rmsd"]) < 1e-3,
        f"top pose RMSD {redock.iloc[0]['rmsd_to_crystal']:.2f} A matches report.json",
    )
    check(
        redock.iloc[0]["rmsd_to_crystal"] < 2.0,
        "the top-ranked pose reproduces the crystal pose (< 2 A)",
    )
    check(
        abs(redock.iloc[0]["affinity"] - report["redock"]["top_affinity"]) < 1e-3,
        f"top affinity {redock.iloc[0]['affinity']:.2f} kcal/mol matches report.json",
    )
    check(
        list(redock["affinity"]) == sorted(redock["affinity"]),
        "poses are ranked best affinity first",
    )

    print("\nScreen")
    check(summary["ligand"].nunique() == len(summary), "one row per ligand in summary.csv")
    for _, row in summary.iterrows():
        best = poses[poses["ligand"] == row["ligand"]]["affinity"].min()
        check(
            abs(best - row["best_affinity"]) < 1e-6,
            f"{row['ligand']}: best affinity {best:.2f} kcal/mol agrees with poses.csv",
        )
    inhibitors = summary[~summary["ligand"].isin(["aspirin", "glucose"])]["best_affinity"].max()
    non_binders = summary[summary["ligand"].isin(["aspirin", "glucose"])]["best_affinity"].min()
    check(
        inhibitors < non_binders,
        f"every kinase inhibitor scores better than every non-binder "
        f"({inhibitors:.2f} vs {non_binders:.2f} kcal/mol)",
    )

    print()
    if failures:
        print(f"{len(failures)} check(s) failed.")
        return 1
    print("All example outputs are consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
