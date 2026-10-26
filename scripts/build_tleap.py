#!/usr/bin/env python3
"""Build the final tleap input from the one written by MCPB.py -s 4.

- (optional) appends nucleotide + Mg (from the minimized prep structure) to the MCPB pdb
- loads GTP.prep / frcmod.phos (Meagher) with the 'gtp' unit alias and the O3 atom type
- loads frcmod.water (H-O-H angle of the Zn-bound water, missing in the MCPB frcmod)
- keeps all MCPB bond commands; replaces solvation / ions / save commands
Usage:
  python build_tleap.py NAME [--ligands GTP_MG_min.pdb] [--water-frcmod frcmod.water]
                        [--box 12.0] [--ion Na+|Cl-|none]
  --ion none  -> dry topology only (used to compute the net charge)
"""
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--ligands")
ap.add_argument("--water-frcmod")
ap.add_argument("--box", type=float, default=12.0)
ap.add_argument("--ion", default="none")
a = ap.parse_args()
N = a.name

# ---- pdb: MCPB site system (+ nucleotide/Mg)
site = [l.rstrip("\n") for l in open(f"{N}_mcpbpy.pdb") if not l.startswith(("END", "CONECT"))]
pdb_out = f"{N}_full.pdb"
if a.ligands:
    atoms = [l for l in site if l.startswith(("ATOM", "HETATM"))]
    resi = max(int(l[22:26]) for l in atoms)
    serial = max(int(l[6:11]) for l in atoms)
    if not site[-1].startswith("TER"):
        site.append("TER")
    prev = None
    for l in open(a.ligands):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[17:20], l[22:26])
        if key != prev:
            if prev is not None:
                site.append("TER")
            resi += 1
            prev = key
        serial += 1
        site.append(f"HETATM{serial:5d}{l[11:22]}{resi:4d}{l[26:].rstrip()}")
    site.append("TER")
open(pdb_out, "w").write("\n".join(site + ["END"]) + "\n")

# ---- tleap input
lines = open(f"{N}_tleap.in").read().splitlines()
keep = [l for l in lines if not re.match(r"^\s*(savepdb|saveamberparm|solvate\w*|addions\w*|quit)\b", l, re.I)]
i = next(k for k, l in enumerate(keep) if re.search(r"=\s*loadpdb", l, re.I))
extra = []
if a.water_frcmod:
    extra.append(f"loadamberparams {a.water_frcmod}")
if a.ligands:
    extra += ['addAtomTypes { { "O3" "O" "sp2" } }', "loadamberprep GTP.prep",
              "loadamberparams frcmod.phos", "GTP = gtp"]
keep[i] = re.sub(r"loadpdb\s+\S+", f"loadpdb {pdb_out}", keep[i])
mol = keep[i].split("=")[0].strip()
keep = keep[:i] + extra + keep[i:]
keep += [f"charge {mol}", f"saveamberparm {mol} {N}_dry.prmtop {N}_dry.inpcrd"]
if a.ion != "none":
    keep += [f"solvatebox {mol} TIP3PBOX {a.box}", f"addions {mol} {a.ion} 0",
             f"savepdb {mol} {N}_solv.pdb",
             f"saveamberparm {mol} {N}_solv.prmtop {N}_solv.inpcrd"]
keep.append("quit")
out = f"{N}_build.in"
open(out, "w").write("\n".join(keep) + "\n")
print(f"{pdb_out} and {out} written (ion: {a.ion})")
