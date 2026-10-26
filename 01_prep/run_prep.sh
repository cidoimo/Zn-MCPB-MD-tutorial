#!/usr/bin/env bash
# STEP 1 - clean the AF3 model and pre-organize the metal sites
#          (nonbonded Zn/Mg, GB minimization with flat-bottom + repulsive restraints)
source "$(dirname "$0")/../scripts/common.sh"
cd "$ROOT/01_prep"

if [[ $HAS_GTP == 1 ]]; then
  for f in GTP.prep frcmod.phos; do [[ -f $f ]] || die "$f missing in 01_prep (see README, Requirements)"; done
fi

step "1.1 clean model"
$PY "$SCRIPTS/clean_model.py" "$ROOT/$MODEL" ${NAME}_clean.pdb --cys "$ZN_CYS" --ids site_ids.sh
source site_ids.sh

step "1.2 nonbonded topology"
{
  echo "source leaprc.protein.ff14SB"
  echo "source leaprc.water.tip3p"
  if [[ $HAS_GTP == 1 ]]; then
    echo 'addAtomTypes { { "O3" "O" "sp2" } }'
    echo "loadamberprep GTP.prep"
    echo "loadamberparams frcmod.phos"
    echo "GTP = gtp"
  fi
  echo "set default PBradii mbondi3"
  echo "mol = loadpdb ${NAME}_clean.pdb"
  echo "check mol"
  echo "charge mol"
  echo "saveamberparm mol prep.prmtop prep.inpcrd"
  echo "quit"
} > prep_tleap.in
tleap -s -f prep_tleap.in > prep_tleap.out
grep -E "Errors =" prep_tleap.out
grep -q "Errors = 0" prep_tleap.out || die "tleap errors, see 01_prep/prep_tleap.out"
cpptraj -p prep.prmtop <<< "resinfo :$((LAST_PROT-1))-$((LAST_PROT+4))"

step "1.3 restraints"
MGL=""; [[ $HAS_GTP == 1 ]] && MGL="$MG_LIGAND_ATOMS"
$PY "$SCRIPTS/make_restraints.py" prep.prmtop prep.inpcrd --mode prep \
    --zn-ligands "$ZN_LIGAND_ATOMS" --mg-ligands "$MGL" > site.rst

step "1.4 GB minimization"
for i in 1 2; do
  wt=$([[ $i == 1 ]] && echo 10.0 || echo 1.0)
  cat > min$i.in <<EOM
GB restrained minimization $i (backbone $wt kcal/mol/A^2 outside free regions)
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=0, igb=8, cut=999.0, saltcon=0.15,
  ntr=1, restraint_wt=$wt, restraintmask='@CA,C,N,O & !(:$FREE_REGIONS)',
  nmropt=1, ntpr=200,
 /
 &wt type='END' /
DISANG=site.rst
EOM
done
MIN="mpirun -np $NP sander.MPI"
$MIN -O -i min1.in -p prep.prmtop -c prep.inpcrd -ref prep.inpcrd -r min1.rst7 -o min1.out
$MIN -O -i min2.in -p prep.prmtop -c min1.rst7   -ref min1.rst7   -r min2.rst7 -o min2.out
grep -ciE "nan|\*\*\*\*" min2.out | grep -q "^0$" || die "NaN in min2.out"

step "1.5 Zn first shell"
WFLAG=""; [[ $HAS_ZN_WATER == 1 ]] && WFLAG="--water"
$PY "$SCRIPTS/check_zn_shell.py" prep.prmtop min2.rst7 --ligands "$ZN_LIGAND_ATOMS" $WFLAG \
  || die "Zn shell wrong: do not run MCPB.py (widen FREE_REGIONS or check the model)"

step "1.6 split: Zn system for MCPB / nucleotide+Mg for later"
STRIP=""
[[ $HAS_GTP == 1 ]] && STRIP=":${NUC_RES},${MG_RES}"
if [[ -n $STRIP ]]; then
  printf "parm prep.prmtop\ntrajin min2.rst7\nstrip %s\ntrajout %s_site_min.pdb pdb nobox\nrun\n" "$STRIP" "$NAME" > split1.cpptraj
  printf "parm prep.prmtop\ntrajin min2.rst7\nstrip !(%s)\ntrajout GTP_MG_min.pdb pdb nobox\nrun\n" "$STRIP" > split2.cpptraj
  cpptraj -i split1.cpptraj > split1.log
  cpptraj -i split2.cpptraj > split2.log
  [[ $(grep -cE " gtp | MG " GTP_MG_min.pdb) -gt 0 ]] || die "GTP_MG_min.pdb empty"
  cp GTP_MG_min.pdb GTP.prep frcmod.phos "$ROOT/03_system/"
else
  printf "parm prep.prmtop\ntrajin min2.rst7\ntrajout %s_site_min.pdb pdb nobox\nrun\n" "$NAME" > split1.cpptraj
  cpptraj -i split1.cpptraj > split1.log
fi
pdb4amber -i ${NAME}_site_min.pdb -o ${NAME}_fixed_H.pdb > pdb4amber.log 2>&1
grep -qE " gtp | MG " ${NAME}_fixed_H.pdb && die "nucleotide/Mg still in ${NAME}_fixed_H.pdb"
cp ${NAME}_fixed_H.pdb "$ROOT/02_mcpb/"
echo "Done. Next: cd ../02_mcpb && ./run_mcpb.sh mol2"
