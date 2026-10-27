# Zn-MCPB-MD: bonded Zn sites with MCPB.py, ORCA and GaMD
© cidoimo - Davide Pietrafesa 

Pipeline and tutorial for the MD of Zn proteins from *Pseudomonas aeruginosa* (Zur regulon), worked out on **PA5535** (ZigA-like COG0523 GTPase, closed state with Mg·GTP), on Plutone server (University of Rome Tor Vergata):

- **Zn²⁺** at the CxCC motif (Cys71, Cys73, Cys74) + one water → **bonded model** (MCPB.py)
- QM with **ORCA** (geometry, Hessian) and **PySCF** (Merz–Kollman ESP) → no Gaussian needed
- **Mg²⁺·GTP**: Mg nonbonded (Li–Merz **12-6-4**), GTP from Meagher et al.
- MD on GPU with **pmemd.cuda**: 4 minimizations, 4 NVT, 1 NPT, **150 ns GaMD**

Everything protein-specific is in **`config.sh`**: the same scripts are meant to be reused for **PA5532** (see [section 9](#9-a-new-protein-pa5532)).

```
00_model ─► 01_prep ──────────────► 02_mcpb ─────────────────────► 03_system ──────────► 04_md ──────► 05_analysis
 AF3 model   clean, nonbonded       mol2 → MCPB 1 → ORCA/PySCF     tleap, post-MCPB      min/NVT/NPT   Zn/Mg site
             GB min + restraints    → MCPB 2,3,4                   fixes, restraints     GaMD 150 ns   monitoring
```

---

## Contents

1. [Requirements and setup on plutone](#1-requirements-and-setup-on-plutone)
2. [Repository layout](#2-repository-layout)
3. [config.sh](#3-configsh)
4. [Step 1 – Model cleaning and site pre-organization](#4-step-1--model-cleaning-and-site-pre-organization)
5. [Step 2 – MCPB.py with ORCA](#5-step-2--mcpbpy-with-orca)
6. [Step 3 – Solvated system and post-MCPB fixes](#6-step-3--solvated-system-and-post-mcpb-fixes)
7. [Step 4 – Equilibration and GaMD on GPU](#7-step-4--equilibration-and-gamd-on-gpu)
8. [Step 5 – Site monitoring](#8-step-5--site-monitoring)
9. [A new protein (PA5532)](#9-a-new-protein-pa5532)
10. [Troubleshooting](#10-troubleshooting)
11. [References](#12-references)

---

## 1. Requirements and setup on plutone

| Software | Where | Notes |
|---|---|---|
| AmberTools / Amber 22 | `/home/dave/programs/amber22` | `source .../amber.sh` (done by the scripts via `AMBER_SH` in `config.sh`); `pmemd.MPI`, `pmemd.cuda` |
| conda env **pyscf** (shared) | `/home/dave/programs/miniconda3/envs/pyscf` | python 3.11, pyscf, parmed, openmpi 4.1.8 |
| ORCA 6.1.1 (shared) | `/home/dave/programs/orca6` | |
| GTP.prep, frcmod.phos | AMBER Parameter Database (Manchester), "GTP – revised phosphate parameters" | put both in `01_prep/` (not in the repo, third-party files) |
| GPU | `CUDA_DEVICE` in `config.sh` | for `pmemd.cuda` |

### 1.1 Shared software on plutone

The conda env (`pyscf`: python 3.11, pyscf, parmed, openmpi 4.1.8) and ORCA 6.1.1 are mine's installs, shared read-only. `config.sh` already points to them:
```bash
QM_PYTHON=/home/dave/programs/miniconda3/envs/pyscf/bin/python
ORCA_DIR=/home/dave/programs/orca6
```
Nothing to install or activate. Check once:
```bash
/home/dave/programs/miniconda3/envs/pyscf/bin/python -c "import pyscf, parmed; print('ok')"
ls /home/dave/programs/orca6/orca
```
If one of them fails, ask me (permissions). Never `conda install` into that env from your account.

The scripts call that python directly with `PYTHONPATH` unset (so it never mixes with Amber's python 3.12 packages) and put its `bin` in `PATH` only for ORCA (`mpirun` 4.1.8). `envs/md-qm.yml` describes the same env, as a backup for another machine.

## 2. Repository layout

```
config.sh                 everything protein-specific
scripts/                  shared code (do not edit)
  common.sh                 env, LC_ALL=C, paths
  clean_model.py            CYS->CYM, keep coordinating water, residue numbering
  make_restraints.py        DISANG files (prep / mg / znow), auto-detects Zn, Mg, water
  check_zn_shell.py         stops the pipeline if the Zn first shell is wrong
  build_tleap.py            final tleap input from the MCPB one (+ GTP/Mg, water angle, ions)
  fix_topology.py           post-MCPB fixes (water H LJ, dihedrals, weak Zn-O, Mg 12-6-4)
  summarize_site.py         mean/sd/min/max of cpptraj .dat files
00_model/                 AF3 model goes here (not in git: see 00_model/README.md)
01_prep/run_prep.sh       step 1
02_mcpb/run_mcpb.sh       step 2   (orca/: com2orca.py, hess2fchk.py, esp_mk.py)
03_system/build_system.sh step 3
04_md/input/*.in          minimization, NVT, NPT, GaMD inputs
04_md/simulazione.sh      step 4 (one replica)
05_analysis/analyze_site.sh  step 5
envs/md-qm.yml
```

## 3. config.sh

```bash
NAME=PA5535
MODEL=00_model/PA5535_AF3_closed.pdb   # put the model here first (not in git)
ZN_CYS="71 73 74"                    # Cys -> CYM
ZN_LIGAND_ATOMS="71@SG 73@SG 74@SG"  # protein atoms bound to Zn
HAS_ZN_WATER=1                       # one water bound to Zn
HAS_GTP=1                            # GTP + Mg present
MG_LIGAND_ATOMS="gtp@O1G gtp@O2B 19@OG"
FREE_REGIONS="15-21,69-76"           # unrestrained during GB minimization
NP=8; MAXCORE=3000; CUDA_DEVICE=0
BOX=12.0
```
Residue numbers are those of the model (the protein must start at residue 1).

## 4. Step 1 – Model cleaning and site pre-organization

```bash
cd 01_prep            # GTP.prep and frcmod.phos must be here
./run_prep.sh
```

**Why.** The AF3 site cannot go straight into MCPB.py. In PA5535: Zn–S 1.6/1.9/4.1 Å (Cys74 not bound), water H pointing at Zn, deformed phosphates, pLDDT ~45–50 on the CxCC loop. With `cut_off 2.8` MCPB would miss Cys74 and fit RESP charges on a distorted model.

**What it does.**

| Sub-step | Action |
|---|---|
| 1.1 | `clean_model.py`: Cys71/73/74 → CYM (HG removed); keeps ZN, MG, GTP and only the water closest to Zn (O only, tleap rebuilds H); hetero residues after the protein: gtp 401, ZN 402, MG 403, WAT 404 → `site_ids.sh` |
| 1.2 | nonbonded topology: ff14SB, TIP3P ions (Zn²⁺, Mg²⁺), GTP Meagher (`GTP = gtp` alias, `O3` atom type), mbondi3 |
| 1.3 | `make_restraints.py --mode prep`: Zn–SG ×3 (2.30 Å), Zn–O(water) (2.05 Å), **repulsive restraints (> 3.4 Å) on every other N/O/S within 8 Å of Zn**, Mg–O1G/O2B/Ser19 OG (2.10 Å) |
| 1.4 | 2 × 5000-step GB minimization (igb 8), backbone restrained (10 → 1) outside `FREE_REGIONS` |
| 1.5 | `check_zn_shell.py`: the first shell (< 2.8 Å) must be exactly 3 SG + water O, otherwise **the pipeline stops** |
| 1.6 | split: protein + Zn + water → `02_mcpb/PA5535_fixed_H.pdb`; gtp + Mg → `03_system/GTP_MG_min.pdb` |

Without the repulsive restraints the nonbonded Zn pulls in Glu242 (6 Å → 1.8 Å) and the backbone carbonyls of Cys71/73, and MCPB then builds a 5–6-coordinate site with the wrong charge (−2 instead of −1).

## 5. Step 2 – MCPB.py with ORCA

```bash
cd ../02_mcpb
./run_mcpb.sh mol2     # ZN.mol2, WAT.mol2 (OW/HW), PA5535.in with ion_ids
./run_mcpb.sh 1        # small/large models + Gaussian-style inputs
./run_mcpb.sh qm       # ORCA + PySCF  (hours: run with nohup)
./run_mcpb.sh 234      # Seminario, RESP, tleap input
```
or `nohup ./run_mcpb.sh all > mcpb.log 2>&1 &`.

**Checks after `1`:**
- `LINK` lines = **4** (3 × Zn–SG, 1 × Zn–O)
- QM charge/multiplicity = **`-1 1`** for small and large model

**QM (`qm`)** – MCPB.py reads only Gaussian or GAMESS files, so:

| File MCPB needs | Produced by | Content |
|---|---|---|
| `PA5535_small_opt.fchk` | ORCA `! B3LYP/G 6-31G* TightOpt Freq` → `hess2fchk.py` | optimized coordinates + Cartesian Hessian |
| `PA5535_large_mk.log` | ORCA H-only optimization → PySCF `esp_mk.py` | Merz–Kollman ESP (4 shells 1.4–2.0 × MK radius, 6 pts/Å², Zn radius from the `.com`) |

Atom order is kept identical to the MCPB `.com` files. Same level of theory as the MCPB default (B3LYP with Gaussian VWN, 6-31G\*).

PA5535 timing on 8 cores: small Opt+Freq (28 atoms) ~3–5 h, large H-opt (53 atoms) ~1 h, ESP ~20 min.

**Checks during `qm`:** 0 imaginary frequencies; all ligands `BOUND` in the optimized small model.

**Checks after `234`:** no negative k in the metal BOND/ANGL lines; site charges sum to the site charge (PA5535: CM1 −0.817, CM2 −0.789, CM3 −0.661, WT1 +0.017, ZN1 +1.250 → **−1.000**).

## 6. Step 3 – Solvated system and post-MCPB fixes

```bash
cd ../03_system
./build_system.sh
```

| Sub-step | Action |
|---|---|
| 3.1 | `frcmod.water`: **HW–Y4–HW angle** (TIP3P 104.52°) – the MCPB frcmod lacks it and tleap stops (`Could not find angle parameter HW - Y4 - HW`) |
| 3.2 | tleap pass 1 (dry) → net charge (PA5535: −24) → choice of counter-ion |
| 3.3 | tleap pass 2: TIP3P box `BOX` Å, neutralized |
| 3.4 | `fix_topology.py` (below) → `PA5535_final.prmtop` |
| 3.5 | restraint files `04_md/input/site_eq.rst` (Mg + Zn–O) and `site_prod.rst` (Zn–O) |

**Post-MCPB fixes (all found the hard way on PA5535):**

| Problem | Symptom | Fix |
|---|---|---|
| TIP3P H has no LJ; S(−0.8)···H(+0.4) of the bound water | minimization explodes, H on top of Cys73 SG (0.00 Å), E ≈ −10⁷ | LJ on the bound-water H: r = 0.6 Å, ε = 0.0157 kcal/mol |
| 1-4 electrostatics S–Zn–O–H through zero-barrier dihedrals | water pushed onto Zn | the 6 L–Zn–O–H dihedrals deleted |
| Seminario Zn–O(water) **k = 0.4** kcal mol⁻¹ Å⁻² (water barely bound at QM level, Zn–O 2.37 Å) | Zn–O collapses to 1.25 Å | if k < 20: Zn–O **k = 150, r₀ = 2.05 Å**; S–Zn–O angles k = 35, Zn–O–H k = 50, θ₀ from geometry (MCPB empirical scheme) |
| Mg 12-6 binds carboxylates too strongly | — | Mg → **12-6-4** (polarizabilities of the MCPB atom types added automatically) |

Zn–S parameters stay **Seminario** (k 76.7 / 56.8 / 57.9, r₀ 2.28 / 2.35 / 2.34 Å); the MCPB empirical formula at the same distances gives 77.6 / 52.3 / 54.3, i.e. the Seminario Zn–S values are consistent.

The `04_md/input` restraint files: Zn–O flat bottom **1.95–2.20 Å** (k 50) in all stages; Mg–O1G/O2B/Ser19 OG (k 10) only in minimization/NVT/NPT.

## 7. Step 4 – Equilibration and GaMD on GPU

```bash
cd ../04_md
mkdir rep1 && cd rep1
nohup ../simulazione.sh ../PA5535_final.prmtop ../PA5535_final.inpcrd > sim.log 2>&1 &
```
Replicas: `rep2`, `rep3` with the same command (`ig=-1` → different seeds).

| Input | Engine | Ensemble | Length | Positional restraints | Site restraints |
|---|---|---|---|---|---|
| `em.in` | pmemd.MPI | min | 5000 | solute heavy 10 | eq |
| `em2.in` | pmemd.MPI | min | 5000 | solute heavy 5 | eq |
| `em3.in` | pmemd.MPI | min | 5000 | backbone 1 | eq |
| `em4.in` | pmemd.MPI | min | 5000 | – | eq |
| `nvt.in` | pmemd.cuda | NVT, 0→300 K | 500 ps | solute heavy 10 | eq |
| `nvt2-4.in` | pmemd.cuda | NVT 300 K | 3 × 500 ps | backbone 5 / 2 / 1 | eq |
| `npt.in` | pmemd.cuda | NPT 1 bar (MC barostat) | 2 ns | – | eq |
| `firstGaMDrun.in` | pmemd.cuda | NVT | 80 ns = **30 ns cMD + 50 ns GaMD eq.** | – | prod (Zn–O) |
| `secondGaMDrun.in` | pmemd.cuda | NVT | **70 ns GaMD production** | – | prod (Zn–O) |

All inputs: dt 2 fs, SHAKE on H, cutoff 10 Å, Langevin thermostat, `lj1264=1` (Mg 12-6-4: **remove it if the system has no Mg**).

**GaMD** (dual boost, `igamd=3`, threshold `iE=1`, σ₀P = σ₀D = 6.0 kcal/mol):

| Run | `ntcmd` | `nteb` | `irest_gamd` | Notes |
|---|---|---|---|---|
| 1 | 15 000 000 (30 ns) | 25 000 000 (50 ns) | 0 | `ntcmdprep`/`ntebprep` 2 ns, `ntave` 50 000; writes `gamd-restart.dat` |
| 2 | 0 | 0 | 1 | reads `gamd-restart.dat` from the same folder: run it **in the same folder** as run 1 |

Analyse only run 2 (and, if needed, the last part of run 1) for reweighting: `gamd02.log` with PyReweighting.

## 8. Step 5 – Site monitoring

```bash
cd ../../05_analysis
./analyze_site.sh ../04_md/rep1
```
Runs on everything already written (also while the MD is running).

**PA5535 reference values** (test MD, 1 + 7 ns, 12-6-4 Mg):

| Quantity | Mean ± sd |
|---|---|
| Zn–SG 71 / 73 / 74 | 2.20 ± 0.07 / 2.24 ± 0.08 / 2.26 ± 0.08 Å |
| Zn–O(water) | 1.93 ± 0.04 Å |
| Glu242–Zn | 6.9 ± 0.4 Å |
| Solvent waters < 2.6 Å from Zn | 0.2 (transient) |
| Mg–O1G / O2B / Ser19 OG | 1.88 / 1.91 / 2.16 Å |
| Mg–carboxylate contacts | 3 (Asp47 bidentate + Glu99), no water in the Mg shell |
| Backbone RMSD | 1.5 ± 0.2 Å |

The Mg coordination by Asp47 and Glu99 comes from the AF3 closed model and is stable with both 12-6 and 12-6-4 Mg: it is a known **limitation** (in P-loop GTPases the Walker B Asp usually binds Mg through a water). It does not affect the Zn site.

## 9. A new protein (PA5532)

The procedure is the same if the protein has a Zn site of the same kind. Change only `config.sh`:

1. `NAME`, `MODEL` (AF3 model with Zn; place the coordinating water by hand if the 4th ligand is a water)
2. `ZN_CYS`, `ZN_LIGAND_ATOMS`: Zn ligands. For **His** ligands put the atom (`120@NE2` or `120@ND1`) in `ZN_LIGAND_ATOMS` and set the tautomer by hand in the model before step 1: the N that is **not** bound to Zn carries the H, so Zn on NE2 → `HID`, Zn on ND1 → `HIE`
3. `HAS_ZN_WATER` = 0 if Zn has no water
4. `HAS_GTP` = 0 if there is no nucleotide/Mg (and remove `lj1264=1` from `04_md/input/*.in`)
5. `FREE_REGIONS`: the loops around the metal sites

Then run steps 1–5 as above. **Do not reuse PA5535 parameters**: every site needs its own MCPB/QM run. The post-MCPB fixes are applied automatically only when needed (e.g. the Zn–O override only if the Seminario k < 20).

Things to check by eye for every new protein: LINK lines, QM charge, `BOUND` ligands after the ORCA optimization, metal force constants in the frcmod, total site charge, first minimization of step 4.

## 10. Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `Atom .R<GTP 401>.A<PA 1> does not have a type` | unit in `GTP.prep` is `gtp` | `GTP = gtp` (done by the scripts) |
| `mask ':GTP@O1G' selects 0 atoms` | residue is `gtp` (lowercase) | use `gtp@...` or residue numbers |
| `could not convert string to float: '-9,0000'` | Italian locale in awk/printf | `export LC_ALL=C` (done in `common.sh`) |
| `Could not perform modeling without setting $AMBERHOME` | amber.sh not sourced | check `AMBER_SH` in `config.sh` |
| `No module named parmed` / `externally-managed-environment` | wrong python | set `QM_PYTHON` (section 1.1) |
| MCPB charge −2, 5–6 LINK lines | Glu242 / backbone O pulled onto Zn | repulsive restraints of step 1 (automatic) |
| `Cartesian orbitals for mol and spherical orbitals for auxmol` | old `esp_mk.py` | use the one in this repo |
| `Could not find angle parameter HW - Y4 - HW` | MCPB frcmod | `frcmod.water` (automatic) |
| `LJ12_6_4Error: Could not find parameters for ATOM_TYPE Y1` | MCPB types not in the 12-6-4 table | automatic in `fix_topology.py` |
| minimization explodes on Cys SG | water H without LJ | automatic in `fix_topology.py` |
| `pmemd.MPI` segfault in minimization | usually exploding forces (see above); otherwise MPI | `mpirun -np 8 sander.MPI` instead of `pmemd.MPI` in `simulazione.sh` |
| `recv ... Connection reset by peer` at the end of a step | OpenMPI teardown message | harmless if `mdout` ends with timings |
| `Unit 5 Error on OPEN: xxx.in` | input deleted / wrong path | inputs are in `04_md/input`, run from `04_md/repN` |

## 11. References

Abramson, J.; Adler, J.; Dunger, J.; et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* **2024**, *630*, 493–500. DOI: 10.1038/s41586-024-07487-w

Jorgensen, W. L.; Chandrasekhar, J.; Madura, J. D.; Impey, R. W.; Klein, M. L. Comparison of simple potential functions for simulating liquid water. *J. Chem. Phys.* **1983**, *79* (2), 926–935. DOI: 10.1063/1.445869

Li, P.; Merz, K. M., Jr. Taking into Account the Ion-Induced Dipole Interaction in the Nonbonded Model of Ions. *J. Chem. Theory Comput.* **2014**, *10* (1), 289–297. DOI: 10.1021/ct400751u

Li, P.; Merz, K. M., Jr. MCPB.py: A Python Based Metal Center Parameter Builder. *J. Chem. Inf. Model.* **2016**, *56* (4), 599–604. DOI: 10.1021/acs.jcim.5b00674

Li, P.; Merz, K. M., Jr. Metal Ion Modeling Using Classical Mechanics. *Chem. Rev.* **2017**, *117* (3), 1564–1686. DOI: 10.1021/acs.chemrev.6b00440

Maier, J. A.; Martinez, C.; Kasavajhala, K.; Wickstrom, L.; Hauser, K. E.; Simmerling, C. ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB. *J. Chem. Theory Comput.* **2015**, *11* (8), 3696–3713. DOI: 10.1021/acs.jctc.5b00255

Meagher, K. L.; Redman, L. T.; Carlson, H. A. Development of polyphosphate parameters for use with the AMBER force field. *J. Comput. Chem.* **2003**, *24* (9), 1016–1025. DOI: 10.1002/jcc.10262

Miao, Y.; Feher, V. A.; McCammon, J. A. Gaussian Accelerated Molecular Dynamics: Unconstrained Enhanced Sampling and Free Energy Calculation. *J. Chem. Theory Comput.* **2015**, *11* (8), 3584–3595. DOI: 10.1021/acs.jctc.5b00436

Neese, F. Software update: The ORCA program system—Version 5.0. *WIREs Comput. Mol. Sci.* **2022**, *12* (5), e1606. DOI: 10.1002/wcms.1606

Nguyen, H.; Roe, D. R.; Simmerling, C. Improved Generalized Born Solvent Model Parameters for Protein Simulations. *J. Chem. Theory Comput.* **2013**, *9* (4), 2020–2034. DOI: 10.1021/ct3010485

Seminario, J. M. Calculation of Intramolecular Force Fields from Second-Derivative Tensors. *Int. J. Quantum Chem.* **1996**, *60* (7), 1271–1277. DOI: 10.1002/(SICI)1097-461X(1996)60:7<1271::AID-QUA8>3.0.CO;2-W

Singh, U. C.; Kollman, P. A. An approach to computing electrostatic charges for molecules. *J. Comput. Chem.* **1984**, *5* (2), 129–145. DOI: 10.1002/jcc.540050204

Sun, Q.; Zhang, X.; Banerjee, S.; et al. Recent developments in the PySCF program package. *J. Chem. Phys.* **2020**, *153* (2), 024109. DOI: 10.1063/5.0006074
