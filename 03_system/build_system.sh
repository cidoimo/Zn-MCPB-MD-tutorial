#!/usr/bin/env bash
# STEP 3 - solvated system: MCPB site + nucleotide/Mg, post-MCPB fixes, MD restraint files
source "$(dirname "$0")/../scripts/common.sh"
cd "$ROOT/03_system"
G=$NAME
for f in ${G}_mcpbpy.pdb ${G}_tleap.in ${G}_mcpbpy.frcmod; do
  [[ -f $f ]] || die "$f missing: run 02_mcpb/run_mcpb.sh 234"
done

step "3.1 H-O-H angle of the Zn-bound water (missing in the MCPB frcmod)"
WFRC=""
if [[ $HAS_ZN_WATER == 1 ]]; then
  OT=""
  for f in *.mol2; do
    case $f in ZN.mol2|WAT.mol2|WAT_pre.mol2) continue;; esac
    n=$(awk '/@<TRIPOS>ATOM/{a=1;next}/@<TRIPOS>/{a=0}a&&NF>5{c++}END{print c+0}' $f)
    if [[ $n == 3 ]]; then OT=$(awk '/@<TRIPOS>ATOM/{a=1;next}/@<TRIPOS>/{a=0}a&&$2=="O"{print $6}' $f); WRES=${f%.mol2}; fi
  done
  [[ -n $OT ]] || die "Zn-bound water mol2 not found"
  printf "Zn-bound water: H-O-H angle (TIP3P geometry)\nANGLE\nHW-%s-HW   100.000   104.520    TIP3P angle\n\n" "$OT" > frcmod.water
  echo "water residue $WRES, O type $OT"; WFRC="--water-frcmod frcmod.water"
fi

LIG=""; [[ $HAS_GTP == 1 ]] && LIG="--ligands GTP_MG_min.pdb"

step "3.2 tleap pass 1 (dry, net charge)"
$PY "$SCRIPTS/build_tleap.py" $G $LIG $WFRC --ion none
tleap -s -f ${G}_build.in > tleap_dry.out
grep -q "Errors = 0" tleap_dry.out || { grep -B2 -A4 -i "error" tleap_dry.out | head -40; die "tleap errors"; }
Q=$($PY -c "import parmed as p; print(round(sum(a.charge for a in p.load_file('${G}_dry.prmtop').atoms)))")
ION=$([[ $Q -lt 0 ]] && echo "Na+" || echo "Cl-")
echo "net charge $Q -> neutralize with $ION"

step "3.3 tleap pass 2 (solvated, box $BOX A)"
$PY "$SCRIPTS/build_tleap.py" $G $LIG $WFRC --box $BOX --ion "$ION"
tleap -s -f ${G}_build.in > tleap_solv.out
grep -E "Errors =" tleap_solv.out
grep -q "Errors = 0" tleap_solv.out || { grep -B2 -A4 -i "error" tleap_solv.out | head -40; die "tleap errors"; }

step "3.4 post-MCPB fixes (water H LJ, L-Zn-O-H dihedrals, weak Zn-O, Mg 12-6-4)"
$PY "$SCRIPTS/fix_topology.py" ${G}_solv.prmtop ${G}_solv.inpcrd ${G}_final.prmtop | tee fix_topology.log
cp ${G}_solv.inpcrd ${G}_final.inpcrd

step "3.5 restraint files for MD (04_md/input)"
mkdir -p "$ROOT/04_md/input"
: > "$ROOT/04_md/input/site_eq.rst"; : > "$ROOT/04_md/input/site_prod.rst"
if [[ $HAS_ZN_WATER == 1 ]]; then
  $PY "$SCRIPTS/make_restraints.py" ${G}_final.prmtop ${G}_final.inpcrd --mode znow > zn_ow.rst
  cat zn_ow.rst >> "$ROOT/04_md/input/site_eq.rst"
  cat zn_ow.rst >> "$ROOT/04_md/input/site_prod.rst"
fi
if [[ $HAS_GTP == 1 ]]; then
  $PY "$SCRIPTS/make_restraints.py" ${G}_final.prmtop ${G}_final.inpcrd --mode mg \
      --mg-ligands "$MG_LIGAND_ATOMS" --rk 10 >> "$ROOT/04_md/input/site_eq.rst"
fi
cat "$ROOT/04_md/input/site_eq.rst"
cp ${G}_final.prmtop ${G}_final.inpcrd site_final.sh "$ROOT/04_md/"
echo "Done. Next: cd ../04_md && see README (simulazione.sh)"
