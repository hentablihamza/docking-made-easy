"""Draw the figure in docs/egfr_1m17_docking.png from the committed example outputs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent
DATA = HERE / "egfr_1m17"
OUT = HERE.parent / "docs" / "egfr_1m17_docking.png"

BLUE = "#2b6cb0"
GREY = "#a0aec0"
RED = "#c53030"
GREEN = "#2f855a"


def main() -> None:
    report = json.loads((DATA / "report.json").read_text())
    redock = pd.read_csv(DATA / "redock.csv")
    summary = pd.read_csv(DATA / "summary.csv")

    fig, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # -- redocking -----------------------------------------------------------
    colours = [GREEN if v <= 2.0 else GREY for v in redock["rmsd_to_crystal"]]
    left.scatter(
        redock["rmsd_to_crystal"], redock["affinity"],
        c=colours, s=90, edgecolor="white", zorder=3,
    )
    left.axvline(2.0, color=RED, ls="--", lw=1.2, zorder=2)
    left.text(2.1, redock["affinity"].max(), "2 Å success\nthreshold",
              color=RED, fontsize=9, va="top")
    top = redock.iloc[0]
    left.annotate(
        f"top-ranked pose\n{top['rmsd_to_crystal']:.2f} Å, {top['affinity']:.2f} kcal/mol",
        xy=(top["rmsd_to_crystal"], top["affinity"]),
        xytext=(top["rmsd_to_crystal"] + 1.6, top["affinity"] + 0.03),
        fontsize=9, va="center",
        arrowprops=dict(arrowstyle="->", color="#4a5568", lw=1),
    )
    left.set_xlabel("RMSD to the crystal pose (Å, symmetry-corrected)")
    left.set_ylabel("Vina affinity (kcal/mol)")
    left.set_title("Redocking erlotinib into 1M17", fontsize=12, loc="left")
    left.grid(alpha=0.25, zorder=0)

    # -- the small screen ----------------------------------------------------
    summary = summary.sort_values("best_affinity", ascending=False)
    colours = [GREY if name in ("aspirin", "glucose") else BLUE for name in summary["ligand"]]
    bars = right.barh(summary["ligand"], summary["best_affinity"], color=colours, zorder=3)
    for bar, value in zip(bars, summary["best_affinity"]):
        right.text(value + 0.12, bar.get_y() + bar.get_height() / 2,
                   f"{value:.2f}", va="center", ha="left", color="white", fontsize=9)
    right.set_xlabel("best Vina affinity (kcal/mol)")
    right.set_title("Six kinase inhibitors and two non-binders (grey)",
                    fontsize=12, loc="left")
    right.set_xlim(summary["best_affinity"].min() - 1.0, 0)
    right.grid(axis="x", alpha=0.25, zorder=0)

    fig.suptitle(
        f"EGFR kinase (PDB {report['pdb_id']}), box from the co-crystal ligand, "
        f"exhaustiveness {report['exhaustiveness']}",
        fontsize=11, y=0.99, x=0.01, ha="left", color="#4a5568",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
