#!/usr/bin/env python3
"""Fail unless the Zn first shell (< --cut A, heavy atoms) is exactly the expected ligands
(+ the coordinated water if present).
Usage: python check_zn_shell.py X.prmtop X.rst7 --ligands "71@SG 73@SG 74@SG" [--water]
"""
import argparse
import sys

import numpy as np
import parmed as pmd

ap = argparse.ArgumentParser()
ap.add_argument("prmtop")
ap.add_argument("coords")
ap.add_argument("--ligands", required=True)
ap.add_argument("--water", action="store_true")
ap.add_argument("--cut", type=float, default=2.8)
a = ap.parse_args()

p = pmd.load_file(a.prmtop, xyz=a.coords)
xyz = np.array(p.coordinates)
sel = lambda m: list(pmd.amber.AmberMask(p, m if m.startswith(":") else ":" + m).Selected())
zn = [at.idx for at in p.atoms if at.name.upper() == "ZN"][0]
expected = {i for m in a.ligands.split() for i in sel(m)}
near = []
for j, at in enumerate(p.atoms):
    if j == zn or at.element_name == "H":
        continue
    d = np.linalg.norm(xyz[j] - xyz[zn])
    if d <= a.cut + 0.6:
        near.append((d, j))
waters = [j for d, j in near if d <= a.cut and len(p.atoms[j].residue.atoms) == 3
          and p.atoms[j].element_name == "O"]
if a.water and waters:
    expected.add(min(waters, key=lambda j: np.linalg.norm(xyz[j] - xyz[zn])))
found = {j for d, j in near if d <= a.cut}
for d, j in sorted(near):
    at = p.atoms[j]
    tag = "LIGAND" if j in expected else ("INTRUDER" if d <= a.cut else "")
    print(f"{at.residue.idx + 1:5d} {at.residue.name:4s} {at.name:4s} {d:5.2f} A  {tag}")
if found == expected:
    print(f"OK: Zn first shell = expected ligands ({len(expected)})")
else:
    print(f"FAIL: missing {len(expected - found)}, intruders {len(found - expected)}")
    sys.exit(1)
