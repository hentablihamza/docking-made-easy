"""The ``dme`` command line."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .version import __version__
from .errors import DockingError
from .box import (
    Box,
    box_from_center,
    box_from_ligand,
    box_from_residues,
    check_box,
    pocket_residues,
)
from .structure import Structure


def _structure(source: str) -> Structure:
    if len(source) == 4 and not Path(source).exists():
        return Structure.fetch(source)
    return Structure.from_file(source)


def _resolve_box(args, structure: Structure | None = None, ligand=None) -> Box:
    if getattr(args, "config", None):
        return Box.from_config(args.config)
    if getattr(args, "center", None):
        size = args.size if args.size else [22.5, 22.5, 22.5]
        if len(size) == 1:
            size = size * 3
        return box_from_center(args.center, size)
    if ligand is None and structure is not None:
        ligand = structure.select_ligand(getattr(args, "ligand_res", None))
    if ligand is None:
        raise DockingError("Give a box: --ligand-res, --center or --config.")
    box = box_from_ligand(
        ligand,
        padding=args.padding,
        mode=getattr(args, "box_mode", "extent"),
        cubic=getattr(args, "cubic", False),
    )
    if args.size:
        size = args.size * 3 if len(args.size) == 1 else args.size
        box = box_from_center(box.center, size)
    return box


# -- commands ----------------------------------------------------------------


def cmd_info(args) -> int:
    structure = _structure(args.source)
    summary = structure.summary()
    if args.json:
        print(json.dumps(summary, indent=2))
        return 0
    print(f"Structure: {args.source}")
    print(f"  atoms            : {summary['n_atoms']}")
    print(f"  chains           : {', '.join(summary['chains'])}")
    print(f"  polymer residues : {summary['n_polymer_residues']}")
    print(f"  waters           : {summary['n_waters']}")
    if summary["ions"]:
        print(f"  ions             : {', '.join(summary['ions'])}")
    if summary["ligands"]:
        print("  ligand candidates:")
        for lig in summary["ligands"]:
            print(f"    {lig['key']:<16} {lig['n_heavy_atoms']} heavy atoms")
    else:
        print("  ligand candidates: none (set the box with --center)")
    return 0


def cmd_box(args) -> int:
    structure = _structure(args.source)
    ligand = None
    if args.residues:
        selected = [structure.select_ligand(r, min_heavy_atoms=1) for r in args.residues]
        box = box_from_residues(selected, padding=args.padding, cubic=args.cubic)
    else:
        if not args.center and not args.config:
            ligand = structure.select_ligand(args.ligand_res)
        box = _resolve_box(args, structure, ligand)

    print(f"Box source : {box.source}")
    print(f"center_x = {box.center_x:.3f}")
    print(f"center_y = {box.center_y:.3f}")
    print(f"center_z = {box.center_z:.3f}")
    print(f"size_x = {box.size_x:.3f}")
    print(f"size_y = {box.size_y:.3f}")
    print(f"size_z = {box.size_z:.3f}")
    print(f"volume   : {box.volume:,.0f} A^3")

    if ligand is not None:
        report = check_box(box, ligand, structure.receptor_atoms())
        print(f"ligand atoms inside: {report['ligand_atoms_inside'] * 100:.0f}%")
        print(f"receptor heavy atoms inside: {report.get('receptor_heavy_atoms_in_box', 0)}")
        if args.pocket:
            residues = pocket_residues(structure.receptor_atoms(), ligand, cutoff=args.cutoff)
            print(f"pocket residues within {args.cutoff:.1f} A ({len(residues)}):")
            print("  " + ", ".join(f"{name}" for name, _ in residues))
        for warning in report["warnings"]:
            print(f"warning: {warning}")
    else:
        for warning in box.warnings():
            print(f"warning: {warning}")

    if args.out:
        box.write_config(args.out)
        print(f"Wrote {args.out}")
    if args.json_out:
        box.write_json(args.json_out)
        print(f"Wrote {args.json_out}")
    if args.pymol:
        Path(args.pymol).write_text(box.to_pymol())
        print(f"Wrote {args.pymol} (run `pymol {args.source} {args.pymol}`)")
    return 0


def cmd_prep_receptor(args) -> int:
    from .receptor import prepare_receptor

    result = prepare_receptor(
        args.source,
        out_dir=args.out_dir,
        name=args.name,
        keep_waters=args.keep_waters,
        keep_ions=args.keep_ions,
        keep_residues=args.keep_residues or (),
        chains=args.chains or (),
        ph=args.ph,
        method=args.method,
        add_missing_residues=args.add_missing_residues,
        keep_oxt=args.keep_oxt,
    )
    print(result.summary())
    return 0


def cmd_prep_ligand(args) -> int:
    from .ligand import prepare_ligand, prepare_ligands

    if args.file:
        prepared, failures = prepare_ligands(
            args.file, out_dir=args.out_dir, seed=args.seed, ph=None if args.no_ph else args.ph
        )
        print(f"Prepared {len(prepared)} ligand(s) in {args.out_dir}")
        for item in prepared:
            print(f"  {item.name:<24} {item.n_heavy_atoms:>3} heavy atoms, "
                  f"{item.n_rotatable_bonds:>2} rotatable bonds")
        for name, error in failures:
            print(f"  FAILED {name}: {error}")
        return 0
    item = prepare_ligand(
        smiles=args.smiles,
        name=args.name,
        out_dir=args.out_dir,
        seed=args.seed,
        ph=None if args.no_ph else args.ph,
    )
    print(f"{item.pdbqt}")
    print(f"  SMILES          : {item.smiles}")
    print(f"  heavy atoms     : {item.n_heavy_atoms}")
    print(f"  rotatable bonds : {item.n_rotatable_bonds}")
    return 0


def cmd_dock(args) -> int:
    from .dock import dock, dock_many
    from .ligand import prepare_ligand, prepare_ligands
    from .results import results_table, save_table, summary_table, poses_to_sdf

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.config or args.center:
        box = _resolve_box(args)
    else:
        if not args.structure:
            raise DockingError(
                "Give --structure (a PDB ID or file with the co-crystal ligand), "
                "or set the box with --center/--config."
            )
        box = _resolve_box(args, _structure(args.structure))
    print(f"Box: {box}")

    ligands: list[tuple[str, Path]] = []
    if args.ligand_pdbqt:
        ligands = [(Path(p).stem, Path(p)) for p in args.ligand_pdbqt]
    elif args.smiles:
        item = prepare_ligand(
            smiles=args.smiles, name=args.name, out_dir=out_dir / "ligands", seed=args.seed
        )
        ligands = [(item.name, item.pdbqt)]
    elif args.ligands:
        prepared, failures = prepare_ligands(args.ligands, out_dir=out_dir / "ligands", seed=args.seed)
        for name, error in failures:
            print(f"  skipped {name}: {error}")
        ligands = [(item.name, item.pdbqt) for item in prepared]
    else:
        raise DockingError("Give --smiles, --ligands or --ligand-pdbqt.")

    def progress(name: str, done: int, total: int) -> None:
        print(f"  [{done}/{total}] {name}", flush=True)

    print(f"Docking {len(ligands)} ligand(s) with exhaustiveness {args.exhaustiveness}...")
    results, failures = dock_many(
        args.receptor,
        ligands,
        box,
        out_dir=out_dir / "poses",
        n_workers=args.workers,
        progress=progress,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        seed=args.seed,
        cpu=args.cpu,
        overwrite=not args.reuse,
    )
    for name, error in failures:
        print(f"  FAILED {name}: {error}")
    if not results:
        return 1

    table = results_table(results)
    summary = summary_table(results)
    print()
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    poses_csv = args.out or out_dir / "poses.csv"
    save_table(table, poses_csv)
    save_table(summary, out_dir / "summary.csv")
    print(f"\nWrote {poses_csv} and {out_dir / 'summary.csv'}")

    if args.sdf:
        sdf_dir = out_dir / "best_poses"
        sdf_dir.mkdir(parents=True, exist_ok=True)
        for result in results:
            poses_to_sdf(
                result.out_pdbqt,
                sdf_dir / f"{result.ligand}.sdf",
                only_best=not args.all_poses,
                affinities=[p.affinity for p in result.poses],
                name=result.ligand,
            )
        print(f"Wrote poses to {sdf_dir}")
    return 0


def cmd_parse(args) -> int:
    from .results import parse_vina_pdbqt, poses_to_dataframe, pose_to_pdb, poses_to_sdf, save_table

    poses = parse_vina_pdbqt(args.out_pdbqt, ligand_name=args.name)
    frame = poses_to_dataframe(poses)
    print(frame.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if args.csv:
        save_table(frame, args.csv)
        print(f"Wrote {args.csv}")
    if args.sdf:
        poses_to_sdf(
            args.out_pdbqt, args.sdf, only_best=args.best_only,
            affinities=[p.affinity for p in poses], name=args.name,
        )
        print(f"Wrote {args.sdf}")
    if args.pdb:
        pose_to_pdb(args.out_pdbqt, args.pdb, mode=args.mode)
        print(f"Wrote {args.pdb}")
    return 0


def cmd_redock(args) -> int:
    from .receptor import prepare_receptor
    from .results import save_table
    from .validate import redock

    out_dir = Path(args.out_dir)
    structure = _structure(args.source)
    ligand = structure.select_ligand(args.ligand_res)
    print(f"Co-crystal ligand: {ligand.key} ({ligand.n_heavy_atoms} heavy atoms)")

    if args.receptor:
        receptor_pdbqt = Path(args.receptor)
    else:
        print("Preparing the receptor...")
        prep = prepare_receptor(
            args.source, out_dir=out_dir / "receptor", name="receptor",
            ph=args.ph, structure=structure, method=args.method,
        )
        print(prep.summary())
        receptor_pdbqt = prep.pdbqt

    box = _resolve_box(args, structure, ligand)
    print(f"Box: {box}")

    result = redock(
        receptor_pdbqt,
        ligand,
        args.smiles,
        box,
        out_dir=out_dir,
        threshold=args.threshold,
        seed=args.seed,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
    )
    print()
    print(result.to_dataframe().to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print()
    print(result.summary())
    if args.csv:
        save_table(result.to_dataframe(), args.csv)
        print(f"Wrote {args.csv}")
    return 0 if result.success else 2


def cmd_doctor(args) -> int:
    from .tools import check_tools

    statuses = check_tools()
    missing = 0
    for status in statuses:
        mark = "ok  " if status.available else "MISS"
        extra = status.version or status.path or ""
        print(f"[{mark}] {status.name}  {extra}")
        if not status.available:
            missing += 1
            print(f"       {status.hint}")
    if missing:
        print(f"\n{missing} tool(s) missing. Box calculation, parsing and ligand prep "
              "still work without Vina; docking does not.")
    else:
        print("\nEverything docking-made-easy needs is installed.")
    return 0


# -- parser ------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dme",
        description="AutoDock Vina docking with protein prep, a box from the "
                    "co-crystal ligand, and result parsing.",
    )
    parser.add_argument("--version", action="version", version=f"docking-made-easy {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_box_options(p: argparse.ArgumentParser) -> None:
        p.add_argument("--ligand-res", help="co-crystal ligand, e.g. AQ4 or AQ4:A:999")
        p.add_argument("--padding", type=float, default=4.0,
                       help="angstrom added around the ligand (default: 4)")
        p.add_argument("--box-mode", choices=["extent", "centroid"], default="extent")
        p.add_argument("--cubic", action="store_true", help="make the box a cube")
        p.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"),
                       help="set the centre manually")
        p.add_argument("--size", type=float, nargs="+", metavar="S",
                       help="set the size manually (one value or three)")
        p.add_argument("--config", help="read the box from a Vina config file")

    p_info = sub.add_parser("info", help="what is in this structure?")
    p_info.add_argument("source", help="PDB ID or .pdb file")
    p_info.add_argument("--json", action="store_true")
    p_info.set_defaults(func=cmd_info)

    p_box = sub.add_parser("box", help="compute the docking box from the co-crystal ligand")
    p_box.add_argument("source", help="PDB ID or .pdb file")
    add_box_options(p_box)
    p_box.add_argument("--residues", nargs="+", help="build the box from these residues instead")
    p_box.add_argument("--pocket", action="store_true", help="list the pocket residues")
    p_box.add_argument("--cutoff", type=float, default=5.0)
    p_box.add_argument("-o", "--out", help="write a Vina config file")
    p_box.add_argument("--json-out", help="write the box as JSON")
    p_box.add_argument("--pymol", help="write a PyMOL script that draws the box")
    p_box.set_defaults(func=cmd_box)

    p_rec = sub.add_parser("prep-receptor", help="PDB to receptor PDBQT")
    p_rec.add_argument("source", help="PDB ID or .pdb file")
    p_rec.add_argument("-o", "--out-dir", default="receptor")
    p_rec.add_argument("--name", default="receptor")
    p_rec.add_argument("--ph", type=float, default=7.4)
    p_rec.add_argument("--method", choices=["auto", "meeko", "openbabel"], default="auto")
    p_rec.add_argument("--keep-waters", action="store_true")
    p_rec.add_argument("--keep-ions", action="store_true")
    p_rec.add_argument("--keep-residues", nargs="+", help="keep these HETATM residues (cofactors)")
    p_rec.add_argument("--chains", nargs="+", help="keep only these chains")
    p_rec.add_argument("--add-missing-residues", action="store_true",
                       help="let PDBFixer model missing loops (a guess, off by default)")
    p_rec.add_argument("--keep-oxt", action="store_true")
    p_rec.set_defaults(func=cmd_prep_receptor)

    p_lig = sub.add_parser("prep-ligand", help="SMILES or SDF to ligand PDBQT")
    p_lig.add_argument("smiles", nargs="?", help="a SMILES string")
    p_lig.add_argument("-f", "--file", help="a .smi/.csv/.sdf file of ligands")
    p_lig.add_argument("--name", default="ligand")
    p_lig.add_argument("-o", "--out-dir", default="ligands")
    p_lig.add_argument("--seed", type=int, default=42)
    p_lig.add_argument("--ph", type=float, default=7.4)
    p_lig.add_argument("--no-ph", action="store_true", help="keep the SMILES protonation as given")
    p_lig.set_defaults(func=cmd_prep_ligand)

    p_dock = sub.add_parser("dock", help="run Vina on one ligand or a library")
    p_dock.add_argument("-r", "--receptor", required=True, help="receptor PDBQT")
    p_dock.add_argument("--smiles", help="dock a single SMILES")
    p_dock.add_argument("--ligands", help="a .smi/.csv/.sdf file of ligands")
    p_dock.add_argument("--ligand-pdbqt", nargs="+", help="already-prepared ligand PDBQT files")
    p_dock.add_argument("--structure", help="PDB ID or file holding the co-crystal ligand")
    add_box_options(p_dock)
    p_dock.add_argument("--name", default="ligand")
    p_dock.add_argument("--exhaustiveness", type=int, default=8)
    p_dock.add_argument("--num-modes", type=int, default=9)
    p_dock.add_argument("--seed", type=int, default=42)
    p_dock.add_argument("--cpu", type=int, help="cores per Vina run")
    p_dock.add_argument("--workers", type=int, default=1, help="ligands docked at the same time")
    p_dock.add_argument("--reuse", action="store_true", help="keep existing pose files")
    p_dock.add_argument("--sdf", action="store_true", help="also write poses as SDF")
    p_dock.add_argument("--all-poses", action="store_true", help="write every mode, not just the best")
    p_dock.add_argument("-o", "--out", help="poses CSV (default: <out-dir>/poses.csv)")
    p_dock.add_argument("--out-dir", default="docking")
    p_dock.set_defaults(func=cmd_dock)

    p_parse = sub.add_parser("parse", help="Vina output PDBQT to a table, SDF or PDB")
    p_parse.add_argument("out_pdbqt")
    p_parse.add_argument("--name")
    p_parse.add_argument("--csv")
    p_parse.add_argument("--sdf")
    p_parse.add_argument("--pdb")
    p_parse.add_argument("--mode", type=int, default=1, help="which pose to write as PDB")
    p_parse.add_argument("--best-only", action="store_true")
    p_parse.set_defaults(func=cmd_parse)

    p_re = sub.add_parser("redock", help="redock the co-crystal ligand and measure the RMSD")
    p_re.add_argument("source", help="PDB ID or .pdb file")
    p_re.add_argument("--smiles", required=True, help="SMILES of the co-crystal ligand")
    add_box_options(p_re)
    p_re.add_argument("-r", "--receptor", help="an already-prepared receptor PDBQT")
    p_re.add_argument("--method", choices=["auto", "meeko", "openbabel"], default="auto")
    p_re.add_argument("--ph", type=float, default=7.4)
    p_re.add_argument("--threshold", type=float, default=2.0)
    p_re.add_argument("--exhaustiveness", type=int, default=8)
    p_re.add_argument("--num-modes", type=int, default=9)
    p_re.add_argument("--seed", type=int, default=42)
    p_re.add_argument("--csv")
    p_re.add_argument("-o", "--out-dir", default="redock")
    p_re.set_defaults(func=cmd_redock)

    p_doc = sub.add_parser("doctor", help="check which external tools are installed")
    p_doc.set_defaults(func=cmd_doctor)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except DockingError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:  # pragma: no cover
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
