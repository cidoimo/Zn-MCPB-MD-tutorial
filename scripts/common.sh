# sourced by every step script
set -euo pipefail
export LC_ALL=C                      # decimal point, never comma (Italian locale breaks mol2/pdb)
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$ROOT/config.sh"
[[ -f $AMBER_SH ]] && source "$AMBER_SH"
[[ -n ${AMBERHOME:-} ]] || { echo "AMBERHOME not set: check AMBER_SH in config.sh"; exit 1; }
SCRIPTS=$ROOT/scripts
[[ -x $QM_PYTHON ]] || { echo "QM_PYTHON not found ($QM_PYTHON): see envs/README"; exit 1; }
# env python without Amber's PYTHONPATH (Amber ships its own python 3.12 packages)
PY="env -u PYTHONPATH $QM_PYTHON"
step() { echo; echo "=================== $* ==================="; }
die()  { echo "ERROR: $*"; exit 1; }
