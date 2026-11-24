# When docking goes wrong

A catalogue of the errors this workflow actually produces, and what each one means.

## Setup and tools

**`vina is not available`**
The binary is not on your `PATH`. `conda install -c conda-forge vina`, or download the release binary and point `DME_VINA` at it:

```bash
export DME_VINA=/opt/autodock/vina_1.2.7_linux_x86_64
dme doctor
```

**`pdbfixer is not available`**
Receptor preparation needs it: `conda install -c conda-forge pdbfixer openmm`. If you cannot install it, protonate the receptor in another tool and pass the result:

```python
prepare_receptor("1M17", protonated="my_receptor_with_hydrogens.pdb")
```

**MGLTools / `prepare_receptor4.py` / ADFR**
You do not need them. PDBFixer plus Meeko replaces them, and Open Babel (`--method openbabel`) is the fallback. MGLTools is Python 2 and should be left alone.

## Structure problems

**`No ATOM/HETATM records found`**
You downloaded mmCIF. Get the `.pdb` format, or convert it first. `dme info 1M17` downloads the right thing for you.

**`No co-crystal ligand found`**
The entry is apo, or the ligand is smaller than the 6 heavy atoms required. Build the box from pocket residues (`dme box apo.pdb --residues HIS41 CYS145`) or give `--center`/`--size`.

**`'LIG' matches 3 residues (LIG:A:1, LIG:B:1, LIG:C:1)`**
Several copies, one per chain. Pick one: `--ligand-res LIG:A:1`. Usually you also want `--chains A`.

**`Several ligands are equally large`**
Two different HETATM groups with the same heavy-atom count. Name the one you mean.

**`unable to build rdkit mol for residue X corresponding to key A:nnn` (Meeko)**
Meeko's residue templates could not match that residue. Common causes, in order:

1. a terminal OXT on a ragged chain end — handled automatically (the OXT is removed; `--keep-oxt` keeps it),
2. a modified or non-standard residue — PDBFixer replaces the ones it knows; otherwise delete it or supply a template,
3. a covalently bound ligand or a crystallographic artefact clashing with the residue.

With `--method auto` the run falls back to Open Babel and tells you. `--method meeko` makes it an error instead.

**Missing loops.** Gaps in the chain are left unmodelled by default and reported as a note. `--add-missing-residues` lets PDBFixer build them, but a loop modelled from nothing is a guess; never let one line the pocket.

## Ligand problems

**`RDKit could not parse the SMILES`**
Check the ring closures, the brackets and the aromatic lowercase atoms. Paste the SMILES into RDKit on its own before blaming the pipeline.

**`RDKit could not generate a 3D conformer`**
Strained macrocycles and cages. Try `--no-ph`, or build the 3D structure elsewhere and pass an SDF.

**Protonation.** The built-in model is deliberately simple: carboxylic, sulfonic and phosphonic acids are deprotonated, aliphatic amines are protonated, everything else is left alone. For tautomers, zwitterions or a pH far from 7.4, prepare the ligand properly and dock the SDF.

## Docking problems

**`Vina reported success but wrote no output`**
The ligand is entirely outside the box. Check the centre: `dme box <structure> --ligand-res <LIG>`.

**Every score is around -3 to -5 kcal/mol, for everything**
The box is in the solvent, or the receptor PDBQT is empty or untyped. `check_box()` warns when fewer than 50 receptor heavy atoms sit inside the box.

**Redocking fails (RMSD > 2 Å)**

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| the crystal pose is among the modes but not ranked first | scoring, not setup | accept it; report the rank, and consider rescoring |
| no pose anywhere near the crystal one | box too small or off-centre | raise `--padding`, check `--pocket` output |
| the pose is flipped | a bridging water was removed, or a histidine is protonated the wrong way | `--keep-waters`, or fix the protonation |
| poses move between seeds | the search has not converged | `--exhaustiveness 32` or higher |
| the ligand is large and flexible | > 10 rotatable bonds | more exhaustiveness, and expect a worse RMSD |

**A screen is slow**
Vina already multithreads one ligand. `--workers N` splits the cores between N ligands, which is faster for a library of small ligands. Shrink the box before raising exhaustiveness: cost grows with the volume.

**Results differ between machines**
Vina is deterministic for a given seed, binary version and thread count. Pin `--seed` *and* `--cpu`, and record the Vina version.

## Reading the numbers

- `affinity` is kcal/mol, more negative being better. It is an estimate from an empirical scoring function and it correlates only loosely with measured affinity.
- `rmsd_lb`/`rmsd_ub` are distances from each mode to the **best mode**, not to anything experimental.
- The only experimental comparison in this package is `rmsd_to_crystal`, from `dme redock`.
- Differences under ~0.5 kcal/mol between ligands are noise. Do not rank a library on the third decimal.
