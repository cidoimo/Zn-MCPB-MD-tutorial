#!/usr/bin/env python3
"""Post-MCPB fixes of the solvated topology (bonded Zn + optional 12-6-4 Mg).

1. Zn-bound water H atoms get a small LJ (r 0.6 A, eps 0.0157): TIP3P H has none and
   collapses onto the thiolate S.
2. L-Zn-O-H dihedrals deleted (zero barrier anyway): removes the 1-4 S...H attraction.
3. If the Seminario Zn-O(water) force constant is < --kmin (weak water, as in PA5535):
     Zn-O bond     -> k = --k-znow, r0 = --r0-znow
     L-Zn-O angles -> k = 35, theta0 = current geometry   (MCPB empirical scheme)
     Zn-O-H angles -> k = 50, theta0 = current geometry
4. Mg (if present) -> 12-6-4 (Li & Merz), polarizabilities of the MCPB atom types added.
5. Writes NAME_final.prmtop and site_final.sh (atom/residue ids for restraints/analysis).

Usage: python fix_topology.py NAME_solv.prmtop NAME_solv.inpcrd NAME_final.prmtop
"""
import argparse
import math
import os

import numpy as np
import parmed as pmd
from parmed.tools import actions as A

ap = argparse.ArgumentParser()
ap.add_argument("prmtop")
ap.add_argument("inpcrd")
ap.add_argument("out")
ap.add_argument("--kmin", type=float, default=20.0)
ap.add_argument("--k-znow", type=float, default=150.0)
ap.add_argument("--r0-znow", type=float, default=2.05)
a = ap.parse_args()

p = pmd.load_file(a.prmtop, xyz=a.inpcrd)
xyz = np.array(p.coordinates)
at = p.atoms
m = lambda i: f"@{i + 1}"


def run(action, *args):
    act = action(p, *args)
    print(f"  {action.__name__} {' '.join(map(str, args))}")
    act.execute()


def angle(i, j, k):
    v1, v2 = xyz[i] - xyz[j], xyz[k] - xyz[j]
    return math.degrees(math.acos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))


zn = [x.idx for x in at if x.name.upper() == "ZN"]
assert len(zn) == 1, f"found {len(zn)} Zn"
zn = zn[0]
ligs = sorted({b.atom1.idx if b.atom2.idx == zn else b.atom2.idx
               for b in p.bonds if zn in (b.atom1.idx, b.atom2.idx)})
water = [j for j in ligs if at[j].element_name == "O" and len(at[j].residue.atoms) == 3]
others = [j for j in ligs if j not in water]
print(f"Zn {at[zn].residue.name}{at[zn].residue.idx + 1}, ligands: "
      + ", ".join(f"{at[j].residue.name}{at[j].residue.idx + 1}@{at[j].name}" for j in ligs))

if water:
    o = water[0]
    hs = [x.idx for x in at[o].residue.atoms if x.element_name == "H"]
    print("1. LJ on the Zn-bound water H")
    run(A.addLJType, ",".join(m(h) for h in hs), "radius", 0.6, "epsilon", 0.0157)

    print("2. delete L-Zn-O-H dihedrals")
    for l in others:
        for h in hs:
            if any({d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx} == {l, zn, o, h}
                   for d in p.dihedrals):
                run(A.deleteDihedral, m(l), m(zn), m(o), m(h))

    b = next(b for b in p.bonds if {b.atom1.idx, b.atom2.idx} == {zn, o})
    print(f"3. Zn-O(water): k = {b.type.k:.1f}, r0 = {b.type.req:.3f}")
    if b.type.k < a.kmin:
        print(f"   k < {a.kmin}: weak water -> hybrid parameters")
        run(A.setBond, m(zn), m(o), a.k_znow, a.r0_znow)
        for l in others:
            run(A.setAngle, m(l), m(zn), m(o), 35.0, round(angle(l, zn, o), 2))
        for h in hs:
            run(A.setAngle, m(zn), m(o), m(h), 50.0, round(angle(zn, o, h), 2))
    else:
        print("   Seminario value kept")

mg = [x.idx for x in at if x.name.upper() == "MG"]
if mg:
    print("4. Mg -> 12-6-4")
    base = os.path.join(os.environ["AMBERHOME"], "dat/leap/parm/lj_1264_pol.dat")
    pol = {}
    for l in open(base):
        t = l.split()
        if len(t) >= 2:
            try:
                pol[t[0]] = float(t[1])
            except ValueError:
                pass
    ref = {"S": pol.get("SH", pol.get("S")), "O": pol.get("OW", pol.get("O")),
           "N": pol.get("N", pol.get("NB")), "H": pol.get("HW", 0.0)}
    # a new MCPB type gets the polarizability of a standard type with the SAME LJ
    # parameters (12-6-4 requires it), else a per-element default
    lj = {}
    for x in at:
        lj.setdefault(x.type, (round(x.rmin, 4), round(x.epsilon, 4)))
    extra = {}
    for t, (r, e) in lj.items():
        if t in pol:
            continue
        twin = [u for u, v in lj.items() if v == (r, e) and u in pol]
        el = next(x.element_name for x in at if x.type == t)
        extra[t] = pol[twin[0]] if twin else (0.0 if el in ("Zn", "Mg") else ref.get(el, 0.0))
    with open("lj_1264_pol_mcpb.dat", "w") as f:
        f.write(open(base).read().rstrip("\n") + "\n")
        for t, v in extra.items():
            f.write(f"{t:<8s}{v:8.3f}\n")
    print(f"   added polarizabilities for: {', '.join(extra) or 'none'}")
    run(A.add12_6_4, f"@{mg[0] + 1}", "watermodel", "TIP3P", "polfile", "lj_1264_pol_mcpb.dat")

print("Zn site after fixes:")
for bd in p.bonds:
    if zn in (bd.atom1.idx, bd.atom2.idx):
        j = bd.atom2.idx if bd.atom1.idx == zn else bd.atom1.idx
        print(f"  Zn-{at[j].residue.name}{at[j].residue.idx + 1}@{at[j].name}: "
              f"k {bd.type.k:7.1f}  r0 {bd.type.req:.3f}  now {np.linalg.norm(xyz[j] - xyz[zn]):.2f}")
A.outparm(p, a.out).execute()

with open("site_final.sh", "w") as f:
    f.write(f"ZN_ATOM={zn + 1}\nZN_RESNAME={at[zn].residue.name}\nZN_RES={at[zn].residue.idx + 1}\n")
    if water:
        f.write(f"WAT_O_ATOM={water[0] + 1}\nWAT_RESNAME={at[water[0]].residue.name}\n"
                f"WAT_RES={at[water[0]].residue.idx + 1}\n")
    if mg:
        f.write(f"MG_ATOM={mg[0] + 1}\nMG_RES={at[mg[0]].residue.idx + 1}\n")
print(f"{a.out} and site_final.sh written")
