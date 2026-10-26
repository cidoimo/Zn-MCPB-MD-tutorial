#!/usr/bin/env bash
# STEP 2 - MCPB.py with ORCA (QM) + PySCF (Merz-Kollman ESP)
# Usage: ./run_mcpb.sh mol2 | 1 | qm | 234 | all
source "$(dirname "$0")/../scripts/common.sh"
cd "$ROOT/02_mcpb"
G=$NAME
PDB=${G}_fixed_H.pdb

mol2() {
  step "2.1 mol2 files + MCPB input"
  [[ -f $PDB ]] || die "$PDB missing: run 01_prep/run_prep.sh"
  resn() { awk -v r="$1" '/^(ATOM|HETATM)/{n=substr($0,18,3); gsub(/ /,"",n); if(n==r) print}' "$PDB"; }
  resn ZN > ZN.pdb
  [[ $(wc -l < ZN.pdb) -eq 1 ]] || die "expected 1 ZN atom in $PDB"
  metalpdb2mol2.py -i ZN.pdb -o ZN.mol2 -c 2
  NAA=""
  if [[ $HAS_ZN_WATER == 1 ]]; then
    resn WAT > WAT_H.pdb
    [[ $(wc -l < WAT_H.pdb) -eq 3 ]] || die "expected 3 WAT atoms in $PDB"
    antechamber -fi pdb -fo mol2 -i WAT_H.pdb -o WAT_pre.mol2 -at amber -c bcc -nc 0 -pf y
    awk 'BEGIN{a=0} /^@<TRIPOS>ATOM/{a=1;print;next} /^@<TRIPOS>/{a=0}
         a && NF>=6 {if($2=="O")$6="OW"; else if($2~/^H/)$6="HW";
         printf "%7d %-8s %10.4f %10.4f %10.4f %-6s %4d %-6s %10.6f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9; next} {print}' \
         WAT_pre.mol2 > WAT.mol2
    NAA="naa_mol2files    WAT.mol2"
  fi
  cp "$AMBERHOME/dat/leap/parm/frcmod.tip3p" .
  ZNID=$(awk '{print substr($0,7,5)+0}' ZN.pdb)
  cat > $G.in <<EOM
original_pdb     $PDB
group_name       $G
cut_off          2.8
ion_ids          $ZNID
ion_mol2files    ZN.mol2
$NAA
frcmod_files     frcmod.tip3p
force_field      ff14SB
water_model      TIP3P
software_version g16
large_opt        1
EOM
  cat $G.in
}

s1() {
  step "2.2 MCPB.py step 1 (models + QM inputs)"
  MCPB.py -i $G.in -s 1 > mcpb_s1.log 2>&1 || { tail -20 mcpb_s1.log; die "MCPB step 1 failed"; }
  echo "LINK lines (one per Zn ligand):"; grep LINK ${G}_standard.fingerprint
  echo "QM charge / multiplicity:"
  grep -hE "^-?[0-9]+ +[0-9]+ *$" ${G}_small_opt.com ${G}_large_mk.com
}

qm() {
  step "2.3 QM: ORCA (small Opt+Freq, large H-opt) + PySCF (MK ESP)"
  export PATH="$(dirname "$QM_PYTHON"):$PATH"    # mpirun 4.1.8 of the md-qm env, for ORCA
  [[ -x $ORCA_DIR/orca ]] || die "ORCA not found in ORCA_DIR=$ORCA_DIR"
  command -v mpirun >/dev/null || die "mpirun not in PATH (needed by ORCA %pal)"
  $PY orca/com2orca.py $G --np $NP --maxcore $MAXCORE
  "$ORCA_DIR/orca" ${G}_small_optfreq.inp > ${G}_small_optfreq.out
  grep -q "ORCA TERMINATED NORMALLY" ${G}_small_optfreq.out || die "small_optfreq failed"
  echo "imaginary modes: $(grep -c 'imaginary mode' ${G}_small_optfreq.out || true)"
  $PY - ${G}_small_optfreq.xyz <<'PYEOF'
import sys, numpy as np
l = open(sys.argv[1]).read().splitlines(); n = int(l[0])
at = [(x.split()[0], np.array(list(map(float, x.split()[1:4])))) for x in l[2:2 + n]]
zn = [c for e, c in at if e.upper() == "ZN"][0]
for i, (e, c) in enumerate(at, 1):
    d = np.linalg.norm(c - zn)
    if e.upper() in ("S", "O", "N") and d < 3.2:
        print(f"  {i:3d} {e:2s} Zn-{e} {d:5.2f} A {'BOUND' if d < 2.6 else 'NOT BOUND'}")
PYEOF
  $PY orca/hess2fchk.py ${G}_small_optfreq.hess ${G}_small_opt.fchk --com ${G}_small_opt.com
  "$ORCA_DIR/orca" ${G}_large_hopt.inp > ${G}_large_hopt.out
  grep -q "ORCA TERMINATED NORMALLY" ${G}_large_hopt.out || die "large_hopt failed"
  $PY orca/esp_mk.py $G --xyz ${G}_large_hopt.xyz --np $NP > esp_mk.out 2>&1 || { tail esp_mk.out; die "ESP failed"; }
  tail -1 esp_mk.out
}

s234() {
  step "2.4 MCPB.py steps 2 (Seminario), 3 (RESP), 4 (tleap input)"
  for s in 2 3 4; do
    MCPB.py -i $G.in -s $s > mcpb_s$s.log 2>&1 || { tail -20 mcpb_s$s.log; die "MCPB step $s failed"; }
  done
  echo "Metal bonds / angles (k must be > 0):"
  awk '/^BOND/,/^ *$/' ${G}_mcpbpy.frcmod | grep -E "M[0-9]"
  awk '/^ANGL/,/^ *$/' ${G}_mcpbpy.frcmod | grep -E "M[0-9]"
  echo "RESP charges of the site residues (sum must be the site charge):"
  for f in $(grep -l "@<TRIPOS>MOLECULE" *.mol2 | grep -vE "^(ZN|WAT|WAT_pre)\.mol2$"); do
    awk '/@<TRIPOS>ATOM/{a=1;next}/@<TRIPOS>/{a=0}a{s+=$9}END{printf "  %-10s %8.3f\n", FILENAME, s}' $f
  done | tee charges.txt
  awk '{s+=$2}END{printf "  TOTAL      %8.3f\n", s}' charges.txt
  cp ${G}_mcpbpy.pdb ${G}_tleap.in ${G}_mcpbpy.frcmod frcmod.tip3p *.mol2 "$ROOT/03_system/"
  echo "Done. Next: cd ../03_system && ./build_system.sh"
}

case "${1:-}" in
  mol2) mol2 ;;  1) s1 ;;  qm) qm ;;  234) s234 ;;
  all)  mol2; s1; qm; s234 ;;
  *) echo "Usage: $0 mol2|1|qm|234|all"; exit 1 ;;
esac
