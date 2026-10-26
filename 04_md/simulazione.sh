#!/bin/bash
# Usage (from a replica folder, e.g. 04_md/rep1):
#   ../simulazione.sh ../PA5535_final.prmtop ../PA5535_final.inpcrd
# Inputs in ../input ; restraint files ../input/site_eq.rst, ../input/site_prod.rst
# Total: 4 min + 4 NVT (2 ns) + NPT (2 ns) + GaMD 150 ns (30 cMD + 50 GaMD eq + 70 GaMD prod)

set -e
export LC_ALL=C
NP=${NP:-8}
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}

#minimizzazione (4)
mpirun -np $NP pmemd.MPI -O -i ../input/em.in -p $1 -c $2 -ref $2 -suffix em1

for i in `seq 2 4`
	do
		x=$(($i-1))
		mpirun -np $NP pmemd.MPI -O -i ../input/em$i.in -p $1 -c restrt.em$x -ref restrt.em$x -suffix em$i
done

#equilibratura --> Ensemble NVT (4)
pmemd.cuda -O -i ../input/nvt.in -p $1 -c restrt.em4 -ref restrt.em4 -suffix nvt1

for i in `seq 2 4`
	do
		x=$(($i-1))
		pmemd.cuda -O -i ../input/nvt$i.in -p $1 -c restrt.nvt$x -ref restrt.nvt$x -suffix nvt$i
done

#equilibratura --> Ensemble NPT (1)
pmemd.cuda -O -i ../input/npt.in -p $1 -c restrt.nvt4 -suffix npt

#GaMD run 1: 30 ns cMD + 50 ns GaMD equilibration (writes gamd-restart.dat)
pmemd.cuda -O -i ../input/firstGaMDrun.in -p $1 -c restrt.npt -inf mdinfo.gamd01 -o mdout.gamd01 -r restrt.gamd01 -x mdcrd.gamd01 -gamd gamd01.log

#GaMD run 2: 70 ns GaMD production (reads gamd-restart.dat of run 1)
pmemd.cuda -O -i ../input/secondGaMDrun.in -p $1 -c restrt.gamd01 -inf mdinfo.gamd02 -o mdout.gamd02 -r restrt.gamd02 -x mdcrd.gamd02 -gamd gamd02.log
