#!/usr/bin/env python3
"""Mean / sd / min / max of every column of the .dat files written by monitor.cpptraj.
Usage: python summarize_site.py [file.dat ...]   (default: all *.dat in the folder)
"""
import glob
import sys

import numpy as np

files = sys.argv[1:] or sorted(glob.glob("*.dat"))
for f in files:
    try:
        head = open(f).readline().split()[1:]
        d = np.atleast_2d(np.loadtxt(f, comments="#"))[:, 1:]
    except (OSError, ValueError, IndexError):
        continue
    print(f"\n{f}  ({len(d)} frames)")
    print(f"{'':22s}{'mean':>8s}{'sd':>8s}{'min':>8s}{'max':>8s}")
    for n, c in zip(head, d.T):
        print(f"{n:22s}{c.mean():8.2f}{c.std():8.2f}{c.min():8.2f}{c.max():8.2f}")
