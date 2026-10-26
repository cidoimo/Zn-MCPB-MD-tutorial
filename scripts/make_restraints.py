#!/usr/bin/env python3
"""Write Amber DISANG restraint files for the Zn and Mg sites.

Zn, Mg, nucleotide and the Zn-bound water are found automatically (by name / distance),
so the same script works on the prep, MCPB and solvated topologies.

Modes
  prep : Zn ligands (2.30 A) + Zn-water O (2.05 A) attractive,
         repulsive (> --rmin) on EVERY other N/O/S within --shell A of Zn,
         Mg ligands (2.10 A)                       -> GB pre-organization
  mg   : Mg ligands only (k = --rk)                -> heating / equilibration
  znow : Zn-water O flat bottom 1.95-2.20 A        -> whole MD (bonded model)

Usage
  python make_restraints.py prep.prmtop prep.inpcrd --mode prep \
      --zn-ligands "71@SG 73@SG 74@SG" --mg-ligands "gtp@O1G gtp@O2B 19@OG" > site.rst
"""
import argparse
import sys

import numpy as np
import parmed as pmd

ap = argparse.ArgumentParser()
ap.add_argument("prmtop")
ap.add_argument("coords")
ap.add_argument("--mode", choices=["prep", "mg", "znow"], required=True)
ap.add_argument("--zn-ligands", default="")
ap.add_argument("--mg-ligands", default="")
ap.add_argument("--rk", type=float, default=100.0)
ap.add_argument("--shell", type=float, default=8.0)
ap.add_argument("--rmin", type=float, default=3.4)
ap.add_argument("--znow", default="1.95,2.20", help="flat-bottom range for Zn-water O")
a = ap.parse_args()

p = pmd.load_file(a.prmtop, xyz=a.coords)
xyz = np.array(p.coordinates)


def mask(m):
    """'71@SG' -> ':71@SG' ; 'gtp@O1G' -> ':gtp@O1G' ; '@123' kept."""
    return m if m.startswith(("@", ":")) else ":" + m


def one(m):
    sel = list(pmd.amber.AmberMask(p, mask(m)).Selected())
    if len(sel) != 1:
        sys.exit(f"ERROR: mask {mask(m)!r} selects {len(sel)} atoms")
    return sel[0]


def name(i):
    at = p.atoms[i]
    return f":{at.residue.idx + 1}{at.residue.name}@{at.name}"


def rst(i, j, r1, r2, r3, r4, k2, k3, label):
    print(f"# {label}: {name(i)} -- {name(j)}")
    print(f" &rst iat={i + 1},{j + 1}, r1={r1:.2f}, r2={r2:.2f}, r3={r3:.2f}, r4={r4:.2f}, "
          f"rk2={k2:.1f}, rk3={k3:.1f}, /")


zn = [at.idx for at in p.atoms if at.name.upper() == "ZN"]
if len(zn) != 1:
    sys.exit(f"ERROR: found {len(zn)} Zn atoms")
zn = zn[0]
wat_o = [at.idx for at in p.atoms if at.element_name == "O" and len(at.residue.atoms) == 3
         and {x.element_name for x in at.residue.atoms} == {"O", "H"}]
wat = min(wat_o, key=lambda j: np.linalg.norm(xyz[j] - xyz[zn])) if wat_o else None
if wat is not None and np.linalg.norm(xyz[wat] - xyz[zn]) > 2.8:
    wat = None
sys.stderr.write(f"Zn {name(zn)}; coordinated water: "
                 f"{name(wat) + f' ({np.linalg.norm(xyz[wat] - xyz[zn]):.2f} A)' if wat is not None else 'none'}\n")

if a.mode == "prep":
    lig = [one(m) for m in a.zn_ligands.split()]
    for j in lig:
        rst(zn, j, 0.0, 2.25, 2.35, 3.30, a.rk, a.rk, "Zn ligand")
    if wat is not None:
        rst(zn, wat, 0.0, 2.00, 2.10, 3.05, a.rk, a.rk, "Zn water")
    n = 0
    for j, at in enumerate(p.atoms):
        if j in lig or j in (zn, wat) or at.element_name not in ("N", "O", "S"):
            continue
        if np.linalg.norm(xyz[j] - xyz[zn]) <= a.shell:
            rst(zn, j, a.rmin - 1.0, a.rmin, 99.0, 100.0, a.rk, 0.0, f"non-ligand > {a.rmin} A")
            n += 1
    sys.stderr.write(f"Zn: {len(lig)} ligands + {'1 water' if wat is not None else 'no water'} "
                     f"+ {n} repulsive restraints\n")

if a.mode in ("prep", "mg") and a.mg_ligands:
    mg = [at.idx for at in p.atoms if at.name.upper() == "MG"]
    if len(mg) != 1:
        sys.exit(f"ERROR: found {len(mg)} Mg atoms")
    for m in a.mg_ligands.split():
        rst(mg[0], one(m), 0.0, 2.05, 2.15, 3.10, a.rk, a.rk, "Mg ligand")

if a.mode == "znow":
    if wat is None:
        sys.exit("ERROR: no water within 2.8 A of Zn")
    lo, hi = map(float, a.znow.split(","))
    rst(zn, wat, lo - 0.45, lo, hi, hi + 0.5, 50.0, 50.0, "Zn-water flat bottom")
