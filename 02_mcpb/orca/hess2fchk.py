#!/usr/bin/env python3
"""ORCA .hess (from Opt+Freq) -> minimal Gaussian fchk readable by MCPB.py -s 2 (Seminario).

MCPB reads only 'Current cartesian coordinates' (bohr) and 'Cartesian Force Constants'
(lower triangle, Hartree/bohr^2). ORCA .hess uses the same units.
Usage: python orca/hess2fchk.py PA5535_small_optfreq.hess PA5535_small_opt.fchk [--com PA5535_small_opt.com]
"""
import argparse
import sys

import numpy as np

Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15, "ZN": 30, "Zn": 30, "MG": 12, "Mg": 12,
     "FE": 26, "Fe": 26, "CU": 29, "Cu": 29, "NA": 11, "Na": 11, "CL": 17, "Cl": 17}


def blocks(lines):
    """Yield (keyword, start index) for every $block."""
    for i, l in enumerate(lines):
        if l.startswith("$"):
            yield l.strip(), i


def read_hess(fname):
    lines = open(fname).read().splitlines()
    pos = dict(blocks(lines))
    if "$hessian" not in pos or "$atoms" not in pos:
        sys.exit("ERROR: $hessian or $atoms block not found - did the Freq job finish?")

    # $hessian: dim, then column blocks (header with column indices, then dim rows)
    i = pos["$hessian"] + 1
    n = int(lines[i].split()[0])
    H = np.zeros((n, n))
    i += 1
    while i < len(lines):
        t = lines[i].split()
        if not t:
            i += 1
            continue
        if t[0].startswith("$"):
            break
        cols = [int(c) for c in t]                     # header line
        for r in range(n):
            row = lines[i + 1 + r].split()
            H[int(row[0]), cols] = [float(v) for v in row[1:1 + len(cols)]]
        i += 1 + n
        if cols[-1] == n - 1:
            break

    # $atoms: natoms, then "El mass x y z" (bohr)
    i = pos["$atoms"] + 1
    nat = int(lines[i].split()[0])
    els, xyz = [], []
    for l in lines[i + 1:i + 1 + nat]:
        t = l.split()
        els.append(t[0])
        xyz.append([float(v) for v in t[2:5]])
    if 3 * nat != n:
        sys.exit(f"ERROR: {nat} atoms but Hessian dimension {n}")
    return els, np.array(xyz), H


def fchk_array(f, title, vals, kind="R"):
    f.write(f"{title:<43s}{kind}   N={len(vals):12d}\n")
    for k in range(0, len(vals), 5):
        f.write("".join(f"{v:16.8E}" for v in vals[k:k + 5]) + "\n")


ap = argparse.ArgumentParser()
ap.add_argument("hess")
ap.add_argument("fchk")
ap.add_argument("--com", help="MCPB small_opt.com: check element order")
a = ap.parse_args()

els, xyz, H = read_hess(a.hess)
asym = np.abs(H - H.T).max()
H = 0.5 * (H + H.T)

if a.com:
    com_els = []
    started = False
    for l in open(a.com):
        t = l.split()
        if not started:
            started = len(t) == 2 and t[0].lstrip("-").isdigit() and t[1].isdigit()
            continue
        if not t:
            break
        com_els.append(t[0].upper())
    if [e.upper() for e in els] != com_els:
        sys.exit("ERROR: element order in .hess differs from the .com file")

with open(a.fchk, "w") as f:
    f.write("MCPB small model - converted from ORCA\n")
    f.write("Freq      RB3LYP                                                      6-31G(d)\n")
    f.write(f"{'Number of atoms':<43s}I     {len(els):12d}\n")
    f.write(f"{'Atomic numbers':<43s}I   N={len(els):12d}\n")
    zs = [Z.get(e, Z.get(e.capitalize(), 0)) for e in els]
    for k in range(0, len(zs), 6):
        f.write("".join(f"{z:12d}" for z in zs[k:k + 6]) + "\n")
    fchk_array(f, "Current cartesian coordinates", xyz.ravel())
    fchk_array(f, "Cartesian Force Constants", H[np.tril_indices(len(H))])

ev = np.linalg.eigvalsh(H)
print(f"{a.fchk}: {len(els)} atoms, Hessian {len(H)}x{len(H)}, max asymmetry {asym:.2e}")
print(f"lowest 8 Hessian eigenvalues (6 ~0 expected at a minimum): {np.round(ev[:8], 5)}")
