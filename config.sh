# ============================================================================
#  Project configuration - the ONLY file to edit for a new protein
#  Every script does: source "$ROOT/config.sh"
# ============================================================================

NAME=PA5535                                  # prefix for all output files
MODEL=00_model/PA5535_AF3_closed.pdb         # input model (protein + ZN [+ GTP + MG] + coordinating water)

# --- Zn site ---------------------------------------------------------------
ZN_CYS="71 73 74"                    # Cys bound to Zn -> become CYM (thiolate)
ZN_LIGAND_ATOMS="71@SG 73@SG 74@SG"  # protein atoms bound to Zn (His: e.g. 120@NE2)
HAS_ZN_WATER=1                       # 1 = one water coordinated to Zn (in the model as WAT/HOH)

# --- Nucleotide / Mg (set HAS_GTP=0 if absent) -----------------------------
HAS_GTP=1
MG_LIGAND_ATOMS="gtp@O1G gtp@O2B 19@OG"   # Mg ligands for equilibration restraints

# --- Pre-organization (GB minimization) -------------------------------------
FREE_REGIONS="15-21,69-76"           # residues NOT restrained (P-loop, Zn loop)

# --- Resources --------------------------------------------------------------
NP=8                                 # CPU cores (GB minimization, ORCA, PySCF)
MAXCORE=3000                         # ORCA memory per core (MB)
CUDA_DEVICE=0                        # GPU used by pmemd.cuda

# --- System and MD ----------------------------------------------------------
BOX=12.0                             # TIP3P buffer (A)
PROD_NS=200                          # production length per replica (ns)
SEG_NS=50                            # production segment length (ns), restartable
REPLICAS=3                           # independent replicas (different seeds)

# --- Software (shared installs on plutone) ----------------------------------
AMBER_SH=/home/dave/programs/amber22/amber.sh
ORCA_DIR=/home/dave/programs/orca6         # shared ORCA 6.1.1 install
QM_ENV=pyscf                         # Dave's conda env (pyscf + parmed + openmpi), shared read-only
QM_PYTHON=${QM_PYTHON:-/home/dave/programs/miniconda3/envs/$QM_ENV/bin/python}   # python of that env
