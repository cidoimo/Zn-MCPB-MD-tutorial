#!/usr/bin/env python3
"""Merz-Kollman ESP of the large model with PySCF -> Gaussian-like log readable by MCPB.py -s 3.

Geometry: ORCA H-only optimized large model (PA5535_large_hopt.xyz), atom order = large_mk.com.
Level: B3LYP (Gaussian VWN) / 6-31G*, same as MCPB default.
Grid: Singh-Kollman, 4 shells (1.4, 1.6, 1.8, 2.0 x MK radius), 6 points/A^2 (= IOp(6/33=2,6/42=6)).
Metal radius taken from the ReadRadii line at the end of large_mk.com.
Usage: python orca/esp_mk.py PA5535 [--xyz PA5535_large_hopt.xyz] [--np 8]
"""
import argparse
import os
import re
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("gname", nargs="?", default="PA5535")
ap.add_argument("--xyz", help="optimized large model (default <gname>_large_hopt.xyz)")
ap.add_argument("--np", type=int, default=8)
ap.add_argument("--basis", default="6-31g*")
ap.add_argument("--xc", default="b3lypg")
ap.add_argument("--density", type=float, default=6.0, help="points per A^2")
a = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", str(a.np))

from pyscf import dft, gto, lib  # noqa: E402
from pyscf.df import incore  # noqa: E402

lib.num_threads(a.np)
BOHR = 0.52917721092
MK = {"H": 1.20, "C": 1.50, "N": 1.50, "O": 1.40, "S": 1.75, "P": 1.80,
      "F": 1.35, "CL": 1.70, "BR": 1.85}   # Gaussian MK defaults (A)

# --- charge, multiplicity, element order and metal radii from large_mk.com
com = open(f"{a.gname}_large_mk.com").read().splitlines()
i = next(k for k, l in enumerate(com) if re.match(r"^\s*-?\d+\s+\d+\s*$", l))
chg, mult = map(int, com[i].split())
els = []
i += 1
while com[i].strip():
    els.append(com[i].split()[0])
    i += 1
for l in com[i:]:
    t = l.split()
    if len(t) >= 2 and t[0].upper() not in MK:
        try:
            MK[t[0].upper()] = float(t[1])
            print(f"radius {t[0]} = {t[1]} A (ReadRadii)")
        except ValueError:
            pass

# --- geometry
xyzf = a.xyz or f"{a.gname}_large_hopt.xyz"
xl = open(xyzf).read().splitlines()
nat = int(xl[0])
geo = [(l.split()[0], np.array(list(map(float, l.split()[1:4])))) for l in xl[2:2 + nat]]
if [g[0].upper() for g in geo] != [e.upper() for e in els]:
    sys.exit("ERROR: element order in the xyz differs from large_mk.com")
xyz = np.array([g[1] for g in geo])
for e in set(els):
    if e.upper() not in MK:
        sys.exit(f"ERROR: no MK radius for {e}")
rad = np.array([MK[e.upper()] for e in els])


# --- Singh-Kollman points
def sphere(n):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5 ** 0.5) * k
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


pts = []
for f in (1.4, 1.6, 1.8, 2.0):
    for j in range(nat):
        r = f * rad[j]
        p = xyz[j] + r * sphere(max(int(4 * np.pi * r * r * a.density), 1))
        d = np.linalg.norm(p[:, None, :] - xyz[None, :, :], axis=2)
        keep = np.all(d >= f * rad[None, :] - 1e-8, axis=1)
        pts.append(p[keep])
pts = np.vstack(pts)
print(f"{nat} atoms, {len(pts)} ESP points")

# --- SCF
mol = gto.M(atom=[(e, x) for e, x in zip(els, xyz)], basis=a.basis, charge=chg,
            spin=mult - 1, unit="Angstrom", verbose=4, cart=True)   # 6-31G* = 6D in Gaussian
mf = dft.RKS(mol) if mult == 1 else dft.UKS(mol)
mf.xc = a.xc
mf.grids.level = 4
mf.conv_tol = 1e-9
mf.max_cycle = 200
mf.chkfile = f"{a.gname}_large_mk.chk"
if os.path.exists(mf.chkfile):          # restart from a previous run
    mf.init_guess = "chkfile"
e = mf.kernel()
if not mf.converged:
    mf = mf.newton()
    e = mf.kernel()
    if not mf.converged:
        sys.exit("ERROR: SCF not converged")
dm = mf.make_rdm1()
if dm.ndim == 3:
    dm = dm[0] + dm[1]

# --- ESP at points (a.u.)
pb = pts / BOHR
vnuc = np.zeros(len(pb))
for z, c in zip(mol.atom_charges(), mol.atom_coords()):
    vnuc += z / np.linalg.norm(pb - c, axis=1)
vele = np.zeros(len(pb))
for s in range(0, len(pb), 200):
    fake = gto.fakemol_for_charges(pb[s:s + 200])
    fake.cart = mol.cart          # point charges: same orbital convention as mol
    ints = incore.aux_e2(mol, fake, intor="int3c2e")
    vele[s:s + 200] = np.einsum("ijp,ij->p", ints, dm)
esp = vnuc - vele

# --- Gaussian-like log
out = f"{a.gname}_large_mk.log"
with open(out, "w") as f:
    f.write(" Merz-Kollman ESP computed with PySCF (converted for MCPB.py)\n")
    f.write(f" SCF Done:  E(R{a.xc.upper()}) = {e:.9f}\n")
    f.write(" Electrostatic Properties Using The SCF Density:\n")
    for k, x in enumerate(xyz, 1):
        f.write(f"{'      Atomic Center %4d is at' % k:<32s}{x[0]:10.6f}{x[1]:10.6f}{x[2]:10.6f}\n")
    for k, x in enumerate(pts, nat + 1):
        f.write(f"{'     ESP Fit Center %4d is at' % k:<32s}{x[0]:10.6f}{x[1]:10.6f}{x[2]:10.6f}\n")
    f.write(" **********************************************************************\n")
    f.write("\n")
    f.write("            Electrostatic Properties (Atomic Units)\n")
    f.write("\n")
    f.write(" -----------------------------------------------------------------\n")
    f.write("    Center     Electric         -------- Electric Field --------\n")
    f.write("               Potential          X             Y             Z\n")
    f.write(" -----------------------------------------------------------------\n")
    for k in range(1, nat + 1):
        f.write(f" {k:6d} Atom    {0.0:14.6f}\n")
    for k, v in enumerate(esp, nat + 1):
        f.write(f" {k:6d} Fit     {v:14.6f}\n")
    f.write(" -----------------------------------------------------------------\n")
    f.write(" Normal termination of PySCF ESP (MCPB converter).\n")
print(f"{out} written: E = {e:.6f} Eh, ESP range {esp.min():.4f} .. {esp.max():.4f} a.u.")
