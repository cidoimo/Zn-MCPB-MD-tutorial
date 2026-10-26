#!/usr/bin/env python3
"""Clean an AF3 model for Amber and write the residue numbering of the hetero groups.

- Zn-bound Cys -> CYM (HG removed)
- keeps ZN, MG, the nucleotide (GTP/GDP) and ONLY the water closest to Zn (oxygen only)
- drops CONECT/TER, renumbers atoms; hetero residues after the last protein residue:
  nucleotide, ZN, MG, WAT (only those present)
- writes site_ids.sh (residue numbers used by the following steps)

Usage: python clean_model.py model.pdb out.pdb --cys "71 73 74" [--ids site_ids.sh]
"""
import argparse
import numpy as np

NUC = {"GTP", "GDP", "GNP", "GSP"}
WAT = {"WAT", "HOH", "H2O", "TIP"}

ap = argparse.ArgumentParser()
ap.add_argument("model")
ap.add_argument("out")
ap.add_argument("--cys", default="", help="Zn-bound Cys residue numbers")
ap.add_argument("--ids", default="site_ids.sh")
a = ap.parse_args()
zn_cys = {int(x) for x in a.cys.split()}

prot, het = [], {}
for l in open(a.model):
    if not l.startswith(("ATOM", "HETATM")):
        continue
    resn, resi, name = l[17:20].strip(), int(l[22:26]), l[12:16].strip()
    xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if l.startswith("ATOM") and resn not in WAT:
        if resn == "CYS" and resi in zn_cys:
            if name == "HG":
                continue
            resn = "CYM"
        prot.append((l[12:16], resn, resi, l[26:].rstrip("\n")))
    else:
        key = "NUC" if resn in NUC else ("WAT" if resn in WAT else resn)
        het.setdefault(key, {}).setdefault(resi, []).append((l, resn, name, xyz))

if "ZN" not in het:
    raise SystemExit("ERROR: no ZN in the model")
zn_xyz = list(het["ZN"].values())[0][0][3]

# coordinating water = water oxygen closest to Zn
if "WAT" in het:
    best = min((np.linalg.norm(at[3] - zn_xyz), resi)
               for resi, ats in het["WAT"].items() for at in ats if at[2].startswith("O"))
    print(f"coordinated water: residue {best[1]}, Zn-O {best[0]:.2f} A")
    het["WAT"] = {best[1]: [at for at in het["WAT"][best[1]] if at[2].startswith("O")]}

last = max(r[2] for r in prot)
out, serial, ids = [], 1, {"LAST_PROT": last}


def emit(rec, aname, resn, resi, tail):
    global serial
    out.append(f"{rec}{serial:5d} {aname} {resn:>3} A{resi:4d}{tail}")
    serial += 1


for r in prot:
    emit("ATOM  ", *r)
out.append("TER")
resi = last
for key in ("NUC", "ZN", "MG", "WAT"):
    if key not in het:
        continue
    resi += 1
    ids[f"{key}_RES"] = resi
    for ats in het[key].values():
        for l, resn, name, xyz in ats:
            resn = "WAT" if key == "WAT" else resn
            emit("HETATM", l[12:16], resn, resi, l[26:].rstrip("\n"))
        break                      # one residue per group
    out.append("TER")
out.append("END")
open(a.out, "w").write("\n".join(out) + "\n")
with open(a.ids, "w") as f:
    for k, v in ids.items():
        f.write(f"{k}={v}\n")
print(f"{a.out}: {serial - 1} atoms; " + ", ".join(f"{k}={v}" for k, v in ids.items()))
