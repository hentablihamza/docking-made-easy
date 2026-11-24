"""Worked example: EGFR kinase (PDB 1M17) with erlotinib.

Runs the whole workflow and writes everything under ``examples/egfr_1m17/``:

1. download 1M17 and look at what is inside,
2. build the docking box from the co-crystal ligand AQ4 (erlotinib),
3. prepare the receptor (ligand/water removal, missing atoms, protonation, typing),
4. redock erlotinib and measure the RMSD to the crystal pose,
5. dock a small set of EGFR inhibitors and non-binders, and rank them.

    python examples/egfr_1m17_example.py

Needs AutoDock Vina, Meeko and PDBFixer. `dme doctor` says what is missing.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from docking_made_easy import (
    Structure,
    box_from_ligand,
    check_box,
    dock_many,
    pocket_residues,
    prepare_ligands,
    prepare_receptor,
    poses_to_sdf,
    redock,
    results_table,
    save_table,
    summary_table,
    write_complex,
)

HERE = Path(__file__).parent
OUT = HERE / "egfr_1m17"

PDB_ID = "1M17"
LIGAND_RES = "AQ4"              # erlotinib, as deposited in 1M17
ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"

#: Five quinazoline EGFR inhibitors, one unrelated kinase inhibitor, and two
#: molecules that have no business binding a kinase hinge.
LIBRARY = [
    ("erlotinib", ERLOTINIB),
    ("gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"),
    ("lapatinib", "CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1"),
    ("afatinib", "CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1"),
    ("vandetanib", "COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCC1CCN(C)CC1"),
    ("imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
]

EXHAUSTIVENESS = 16
SEED = 42


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {"pdb_id": PDB_ID, "exhaustiveness": EXHAUSTIVENESS, "seed": SEED}
    started = time.perf_counter()

    # 1. What is in the structure? ------------------------------------------
    structure = Structure.fetch(PDB_ID)
    summary = structure.summary()
    print(json.dumps(summary, indent=2))
    report["structure"] = summary

    ligand = structure.select_ligand(LIGAND_RES)
    print(f"\nCo-crystal ligand: {ligand.key} ({ligand.n_heavy_atoms} heavy atoms)")

    # 2. The box, straight from the ligand ----------------------------------
    box = box_from_ligand(ligand, padding=4.0)
    receptor_atoms = structure.receptor_atoms()
    box_report = check_box(box, ligand, receptor_atoms)
    pocket = pocket_residues(receptor_atoms, ligand, cutoff=5.0)
    print(f"Box: {box}")
    print(f"Pocket residues within 5 A ({len(pocket)}): "
          + ", ".join(name for name, _ in pocket))
    box.write_config(OUT / "box.txt")
    box.write_json(OUT / "box.json")
    (OUT / "box.pml").write_text(box.to_pymol())
    report["box"] = box.to_dict()
    report["box_check"] = box_report
    report["pocket_residues"] = [{"residue": n, "distance": round(d, 2)} for n, d in pocket]

    # 3. Receptor preparation ------------------------------------------------
    prep = prepare_receptor(PDB_ID, out_dir=OUT / "receptor", name="1m17", structure=structure)
    print("\n" + prep.summary())
    report["receptor"] = {
        "method": prep.method,
        "ph": prep.ph,
        "atoms": prep.n_atoms,
        "hydrogens": prep.n_hydrogens,
        "notes": prep.notes,
    }

    # 4. Redocking: does the setup reproduce the crystal pose? ---------------
    print("\nRedocking the co-crystal ligand...")
    check = redock(
        prep.pdbqt, ligand, ERLOTINIB, box,
        out_dir=OUT / "redock", exhaustiveness=EXHAUSTIVENESS, seed=SEED,
    )
    print(check.to_dataframe().to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(check.summary())
    save_table(check.to_dataframe(), OUT / "redock.csv")
    report["redock"] = {
        "top_pose_rmsd": round(check.top_pose_rmsd, 3),
        "best_rmsd": round(check.best_rmsd, 3),
        "top_affinity": round(check.top_affinity, 3),
        "success": check.success,
        "rmsd_per_mode": [round(v, 3) for v in check.rmsd_per_mode],
        "seconds": round(check.docking.seconds, 1),
    }
    poses_to_sdf(
        check.docking.out_pdbqt, OUT / "erlotinib_redocked.sdf",
        affinities=[p.affinity for p in check.docking.poses], name="erlotinib_redocked",
    )

    # 5. A small screen ------------------------------------------------------
    smi = OUT / "library.smi"
    smi.write_text("".join(f"{smiles} {name}\n" for name, smiles in LIBRARY))
    prepared, failures = prepare_ligands(smi, out_dir=OUT / "ligands", seed=SEED)
    for name, error in failures:
        print(f"  skipped {name}: {error}")

    print(f"\nDocking {len(prepared)} ligands (exhaustiveness {EXHAUSTIVENESS})...")
    results, dock_failures = dock_many(
        prep.pdbqt,
        [(item.name, item.pdbqt) for item in prepared],
        box,
        out_dir=OUT / "poses",
        n_workers=2,
        exhaustiveness=EXHAUSTIVENESS,
        seed=SEED,
        progress=lambda name, done, total: print(f"  [{done}/{total}] {name}", flush=True),
    )
    for name, error in dock_failures:
        print(f"  FAILED {name}: {error}")

    table = results_table(results)
    summary_frame = summary_table(results)
    save_table(table, OUT / "poses.csv")
    save_table(summary_frame, OUT / "summary.csv")
    print()
    print(summary_frame.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    report["screen"] = json.loads(summary_frame.to_json(orient="records"))

    best = min(results, key=lambda r: r.best_affinity)
    poses_to_sdf(
        best.out_pdbqt, OUT / f"{best.ligand}_best_pose.sdf", only_best=True,
        affinities=[p.affinity for p in best.poses], name=best.ligand,
    )

    report["seconds_total"] = round(time.perf_counter() - started, 1)
    (OUT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nEverything written to {OUT}")


if __name__ == "__main__":
    main()
