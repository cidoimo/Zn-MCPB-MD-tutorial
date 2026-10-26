#!/usr/bin/env bash
# STEP 5 - metal-site monitoring of one replica
# Usage: ./analyze_site.sh ../04_md/rep1
source "$(dirname "$0")/../scripts/common.sh"
REP=$(cd "${1:?replica folder}" && pwd)
source "$ROOT/04_md/site_final.sh"
TOP=$ROOT/04_md/${NAME}_final.prmtop
OUT=$ROOT/05_analysis/$(basename "$REP"); mkdir -p "$OUT"; cd "$OUT"

{
  echo "parm $TOP"
  for t in mdcrd.nvt1 mdcrd.nvt2 mdcrd.nvt3 mdcrd.nvt4 mdcrd.npt mdcrd.gamd01 mdcrd.gamd02; do
    [[ -s $REP/$t ]] && echo "trajin $REP/$t"
  done
  echo "autoimage"
  for l in $ZN_LIGAND_ATOMS; do
    echo "distance zn_${l//@/_} @$ZN_ATOM :$l out zn_site.dat"
  done
  [[ -n ${WAT_O_ATOM:-} ]] && echo "distance zn_water @$ZN_ATOM @$WAT_O_ATOM out zn_site.dat"
  echo "watershell @$ZN_ATOM out zn_solvent.dat lower 2.6 upper 3.5"
  if [[ -n ${MG_ATOM:-} ]]; then
    for l in $MG_LIGAND_ATOMS; do echo "distance mg_${l//@/_} @$MG_ATOM :$l out mg_site.dat"; done
    echo "watershell @$MG_ATOM out mg_solvent.dat lower 2.6 upper 3.5"
    echo "nativecontacts @$MG_ATOM (:ASP@OD1,OD2|:GLU@OE1,OE2) distance 2.6 out mg_carboxylates.dat"
  fi
  echo "rms bb_rmsd first @CA,C,N out rmsd.dat"
  echo "atomicfluct :1-$((ZN_RES-1))@CA byres out rmsf.agr"
  echo "run"
} > monitor.cpptraj
cpptraj -i monitor.cpptraj > monitor.log 2>&1 || { tail -20 monitor.log; die "cpptraj failed"; }
$PY "$SCRIPTS/summarize_site.py" zn_site.dat zn_solvent.dat mg_site.dat mg_solvent.dat mg_carboxylates.dat rmsd.dat
echo "Data in $OUT"
