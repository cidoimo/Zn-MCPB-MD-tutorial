#!/usr/bin/env python3
"""Turn the Gaussian inputs written by MCPB.py -s 1 into ORCA inputs.

  PA5535_small_opt.com -> PA5535_small_optfreq.inp   (B3LYP/G 6-31G*, Opt + analytic Freq)
  PA5535_large_mk.com  -> PA5535_large_hopt.inp      (H-only optimization, heavy atoms fixed)

Atom order is kept exactly as in the .com files (it must match the MCPB fingerprints).
Usage: python orca/com2orca.py PA5535 [--np 8] [--maxcore 3000]
"""
import argparse
import re


def read_com(fname):
    """Return charge, multiplicity, atoms [(el, x, y, z, flag)], trailing lines."""
    lines = open(fname).read().splitlines()
    i = 0
    # skip Link0, route, title until the charge/multiplicity line
    while not re.match(r"^\s*-?\d+\s+\d+\s*$", lines[i]):
        i += 1
    chg, mult = map(int, lines[i].split())
    atoms = []
    i += 1
    while i < len(lines) and lines[i].strip():
        t = lines[i].split()
        if len(t) == 5:      # element flag x y z  (large model: 0 = free, -1 = frozen)
            atoms.append((t[0], float(t[2]), float(t[3]), float(t[4]), int(t[1])))
        else:                # element x y z
            atoms.append((t[0], float(t[1]), float(t[2]), float(t[3]), 0))
        i += 1
    return chg, mult, atoms, lines[i:]


def xyz_block(chg, mult, atoms):
    out = [f"* xyz {chg} {mult}"]
    out += [f"  {el:<3s} {x:12.6f} {y:12.6f} {z:12.6f}" for el, x, y, z, _ in atoms]
    out.append("*")
    return "\n".join(out)


ap = argparse.ArgumentParser()
ap.add_argument("gname", nargs="?", default="PA5535")
ap.add_argument("--np", type=int, default=8)
ap.add_argument("--maxcore", type=int, default=3000, help="MB per core")
a = ap.parse_args()
common = f"%pal nprocs {a.np} end\n%maxcore {a.maxcore}\n"

# small model: optimization + analytic Hessian (same level as MCPB default)
chg, mult, atoms, _ = read_com(f"{a.gname}_small_opt.com")
with open(f"{a.gname}_small_optfreq.inp", "w") as f:
    f.write("! B3LYP/G 6-31G* TightOpt Freq DefGrid3 TightSCF\n")
    f.write(common)
    f.write(xyz_block(chg, mult, atoms) + "\n")
print(f"{a.gname}_small_optfreq.inp: {len(atoms)} atoms, charge {chg}, mult {mult}")

# large model: optimize H only (MCPB large_opt 1), heavy atoms fixed
chg, mult, atoms, _ = read_com(f"{a.gname}_large_mk.com")
heavy = [i for i, at in enumerate(atoms) if at[0].upper() != "H"]
with open(f"{a.gname}_large_hopt.inp", "w") as f:
    f.write("! B3LYP/G 6-31G* Opt DefGrid3 TightSCF\n")
    f.write(common)
    f.write("%geom\n  Constraints\n")
    for i in heavy:
        f.write(f"    {{ C {i} C }}\n")
    f.write("  end\nend\n")
    f.write(xyz_block(chg, mult, atoms) + "\n")
print(f"{a.gname}_large_hopt.inp: {len(atoms)} atoms ({len(heavy)} fixed), charge {chg}, mult {mult}")
