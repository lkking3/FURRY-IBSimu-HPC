#!/bin/bash
#SBATCH -J ibsimu_pipeline
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH -t 00:20:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Lightweight ORCHESTRATOR: prepares and submits sweep -> reduce -> (optional) refine.

set -euo pipefail
# ---- debug / trace (set DEBUG=1 to enable) ----
DEBUG="${DEBUG:-0}"
if [[ "$DEBUG" == "1" ]]; then set -x; fi
trap 'echo "[orchestrator] ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR


umask 002

# ---- LMOD init (Omega) ----
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
[ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash || true

# =====================================================================================
# USER CONFIG
# =====================================================================================

# Resources for Stage-1 compute tasks (per array element)
SWEEP_CPUS=2
SWEEP_MEM=2GB
SWEEP_TIME=24:00:00
SWEEP_PART=medium       # comment out to skip setting partition

# Array concurrency (how many Stage-1 tasks at once)
ARRAY_CONC=5

# Results root & run name prefix
RESULTS_DIR="${RESULTS_DIR:-$PWD/results}"
RUN_PREFIX="two_grid_Convergence_Study_Idealized_Grid_offset"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"

# ---- Stage-1 verticality sampling knobs (passed to compute jobs) ----
TRACE_ENABLE=0            # 0=metrics only; 1=also write CSV tracers
VERT_PAD_M=5e-4
VERT_XLEN_M=5e-3
VERT_YPAD_M=5e-4
VERT_STEP_M=1e-4
VERT_THRESH_DEG=2.0

# ---- Base geometry (non-swept defaults) ----
H=8e-5
X_LEFT_M=2e-3
X_RIGHT_M="${X_RIGHT_M:-0.02}"
# Physical tube/sample length used for envelope projection (can exceed X_RIGHT_M)
X_RIGHT_PHYS_M="${X_RIGHT_PHYS_M:-0.55}"
YBOX_M=0.086
ACCEL_OFF_Y_M=0.0784
SCREEN_OFF_Y_M=0.0784

# Offset sweep controls (defaults = single value at the base offsets)
ACCEL_OFF_START_M="${ACCEL_OFF_START_M:-$ACCEL_OFF_Y_M}"
ACCEL_OFF_STEP_M="${ACCEL_OFF_STEP_M:-0.00005}"
ACCEL_OFF_COUNT="${ACCEL_OFF_COUNT:-10}"
# Potentials & downstream slices
VS_V=0
VA_V=-10000
SAMPLE_V=$VA_V
TUBE_WALL_T_M=0.0002
SAMPLE_PLATE_T_M=0.0002
SAMPLE_PROFILE_X_M="${SAMPLE_PROFILE_X_M:-0.03}"

# Auto-plot sample profile JSONs after each run
PLOT_PROFILE="${PLOT_PROFILE:-0}"
PLOT_PROFILE_Y="${PLOT_PROFILE_Y:-I_A}"
PLOT_PROFILE_X_UNITS="${PLOT_PROFILE_X_UNITS:-mm}"
PLOT_PROFILE_EXT="${PLOT_PROFILE_EXT:-png}"

RUN_SOLVE=1
WRITE_PNG=1
PNG_NAME="two_grid_geom.png"
CLOSE_PAD_M=5e-4
CLOSE_BASE_PX=800
PLOT_FONT_PX=${PLOT_FONT_PX:-22}

# Ion / plasma toggles
# Stage 1: field-only (fast) by default
# Stage 2: enable ions manually or by setting this to 1
ENABLE_IONS_STAGE1=${ENABLE_IONS_STAGE1:-1}
ENABLE_IONS_STAGE2=${ENABLE_IONS_STAGE2:-1}   # set to 1 when you want Stage-2 with ions

# Plasma defaults (match C++ defaults; you can override here or via sbatch --export)
PLASMA_NI_M3=${PLASMA_NI_M3:-2e18}
PLASMA_TE_EV=${PLASMA_TE_EV:-3.7}

# Beam / SC tuning knobs
ION_J_SCALE=${ION_J_SCALE:-1.0} #1.0 = nominal Bohm current
SC_FACTOR=${SC_FACTOR:-0.0005} # 1.0 = full space charge, <1.0 = partially neutralized

# =====================================================================================
# SWEEP LISTS
# =====================================================================================
# Geometry
AP_RAD_LIST=(0.0015)               # 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005)
GRID_T_LIST=(0.005)                #0.004 0.005 0.006 0.007 0.008 0.009 0.01)
GAP_LIST=(0.001)

# Chamfers (8 independent controls)
SCR_UP_DEPTH_LIST=(0.0)
SCR_UP_ANGLE_LIST=(0.0)
SCR_DN_DEPTH_LIST=(0.0) #0.002 0.003 0.004)
SCR_DN_ANGLE_LIST=(0) #20 45)
ACC_UP_DEPTH_LIST=(0.0) #0.002 0.003 0.004)
ACC_UP_ANGLE_LIST=(0) #20 45)
ACC_DN_DEPTH_LIST=(0.0)
ACC_DN_ANGLE_LIST=(0.0)

# Offset sweeps
make_seq() {
  python3 - "$1" "$2" "$3" <<'PY'
import sys, math
start = float(sys.argv[1])
step = float(sys.argv[2])
try:
    count = int(float(sys.argv[3]))
except Exception:
    count = 1
if count < 1:
    count = 1
vals = [start + i * step for i in range(count)]
print(" ".join(f"{v:.12g}" for v in vals))
PY
}

ACCEL_OFF_LIST=( $(make_seq "$ACCEL_OFF_START_M" "$ACCEL_OFF_STEP_M" "$ACCEL_OFF_COUNT") )

# =====================================================================================
# Stage-2 control
# =====================================================================================
RUN_STAGE2="${RUN_STAGE2:-1}"   # 0=no refine; 1=submit refine after reducer

# FIX: set an explicit knob list. If you want auto-pick later, set OPT_KNOBS="".
# (current default: chamfer-only tweaks)
OPT_KNOBS="${OPT_KNOBS:-ACCEL_OFF_Y_M}"

# These are ignored when OPT_KNOBS is non-empty, but keep sane defaults:
MAX_KNOBS="${MAX_KNOBS:-3}"
VERT_METRIC="${VERT_METRIC:-mean_angle_deg}"   # or e.g. mean_angle_deg
TOPK="${TOPK:-10}"

# Default delta floors used when reducer suggests step sizes
DELTA_FLOORS_JSON="${DELTA_FLOORS_JSON:-$(
  cat <<'JSON'
{
  "AP_RAD_M":         2.0e-4,
  "GAP_M":            2.0e-4,
  "GRID_T_M":         1.0e-3,
  "ACCEL_OFF_Y_M":    0.00001,
  "SCREEN_OFF_Y_M":   5.0e-4,
  "SCR_UP_DEPTH_M":   5.0e-4,
  "SCR_DN_DEPTH_M":   5.0e-4,
  "ACC_UP_DEPTH_M":   5.0e-4,
  "ACC_DN_DEPTH_M":   5.0e-4,
  "SCR_UP_ANGLE_DEG": 2.0,
  "SCR_DN_ANGLE_DEG": 2.0,
  "ACC_UP_ANGLE_DEG": 2.0,
  "ACC_DN_ANGLE_DEG": 2.0
}
JSON
)}"

# Encode JSON safely for sed injection later
DF_B64="$(printf '%s' "$DELTA_FLOORS_JSON" | base64 -w0 2>/dev/null || printf '%s' "$DELTA_FLOORS_JSON" | base64)"

# -----------------------------------------------------------------------------
# Stage-2 optimizer hyperparameters (used by stage2_opt.py)
# -----------------------------------------------------------------------------
# Choose optimizer mode:
#   smart = feedback-driven hill-climber (stage2_opt.py)
#   grid  = legacy Cartesian 3^N search (inline Python)
STAGE2_MODE="${STAGE2_MODE:-smart}"

# Outer loop + stopping criteria
ST2_MAX_ITERS="${ST2_MAX_ITERS:-40}"      # max outer iterations
ST2_NO_IMPROVE="${ST2_NO_IMPROVE:-5}"    # stop after N iters with no better score
ST2_BATCH_SIZE="${ST2_BATCH_SIZE:-5}"    # concurrent srun evaluations
ST2_SEED="${ST2_SEED:-42069}"                # RNG seed (0 = deterministic per campaign)

# Stage-2 scoring configuration
# mode: lex (tuple ordering) or weighted (single scalar)
ST2_SCORE_MODE="${ST2_SCORE_MODE:-weighted}"
ST2_I_FIELD="${ST2_I_FIELD:-currents.I_ag_out_A|I_ag_out_A}"
ST2_DIV_FIELD="${ST2_DIV_FIELD:-DIVERGENCE_ANGLE_DEG|divergence_angle_deg}"
ST2_Y_FIELD="${ST2_Y_FIELD:-deflection.y_mean_pred_500mm_m|y_mean_pred_500mm_m|deflection.y_mean_right_m}"
ST2_W_I="${ST2_W_I:-0.5}"
ST2_W_DIV="${ST2_W_DIV:-0.1}"
ST2_W_Y="${ST2_W_Y:-1.0}"
ST2_I_REF_A="${ST2_I_REF_A:-0.004}"   # reference current for scoring (e.g. 4 mA)
ST2_DIV_REF_DEG="${ST2_DIV_REF_DEG:-1.0}"
ST2_Y_REF_M="${ST2_Y_REF_M:-1.0e-3}"
ST2_Y_ABS="${ST2_Y_ABS:-1}"

# Legacy weights below are unused but kept for compatibility.
ST2_W_VERT="${ST2_W_VERT:-1.0}"
ST2_W_COLL="${ST2_W_COLL:-1.0}"
ST2_W_CURR="${ST2_W_CURR:-1.0}"   # weight for "distance from 4A" term

# Current target: we want a grid that *could* reach about 4 A total
ST2_I_TARGET_A="${ST2_I_TARGET_A:-4.0}"    # target total current for full grid [A]

# Minimum acceptable total current as a fraction of target.
# e.g. 0.25 -> require at least 1 A potential total current to not discard.
ST2_I_MIN_FRAC="${ST2_I_MIN_FRAC:-0.66}"

# Optional explicit floor in amps; if unset, we derive from fraction * target.
ST2_I_MIN_TOTAL_A="${ST2_I_MIN_TOTAL_A:-}"

# (keep existing export line, just extend it)
export RUN_STAGE2 OPT_KNOBS MAX_KNOBS VERT_METRIC TOPK \
       STAGE2_MODE ST2_MAX_ITERS ST2_NO_IMPROVE ST2_BATCH_SIZE ST2_SEED \
       ST2_SCORE_MODE ST2_I_FIELD ST2_DIV_FIELD ST2_Y_FIELD \
       ST2_W_I ST2_W_DIV ST2_W_Y ST2_I_REF_A ST2_DIV_REF_DEG ST2_Y_REF_M ST2_Y_ABS \
       ST2_W_VERT ST2_W_COLL ST2_W_CURR \
       ST2_I_TARGET_A ST2_I_MIN_FRAC ST2_I_MIN_TOTAL_A \
       ST2_THETA_REF_DEG ST2_R_REF_M ST2_I_REF_A ST2_I_MIN_A

# =====================================================================================
# INTERNALS (no edits)
# =====================================================================================

PIPE_ROOT="$PWD/.pipeline_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$PIPE_ROOT"
echo "[orchestrator] workspace: $PIPE_ROOT"
cp -f plot_runlog_compact.py "$PIPE_ROOT/plot_runlog_compact.py"

# ---------- helpers ----------
arrlit(){ local -n A="$1"; printf "("; for v in "${A[@]}"; do printf " %s" "$v"; done; printf " )"; }

# Escape strings for safe sed replacement (we use '#' delimiter)
sed_escape() {
  printf '%s' "$1" | sed -e 's/[\/&]/\\&/g' -e 's/#/\\#/g'
}
arrlit_escaped(){ local -n A="$1"; sed_escape "$(arrlit "$1")"; }

# =====================================================================================
# Write & fill sweep_lists.env
# =====================================================================================
cat > "$PIPE_ROOT/sweep_lists.env" <<'EOF_ENV'
# shellcheck shell=bash
# Populated by orchestrator.

AP_RAD_LIST=__AP_RAD_LIST__
GRID_T_LIST=__GRID_T_LIST__
GAP_LIST=__GAP_LIST__
SCR_UP_DEPTH_LIST=__SCR_UP_DEPTH_LIST__
SCR_UP_ANGLE_LIST=__SCR_UP_ANGLE_LIST__
SCR_DN_DEPTH_LIST=__SCR_DN_DEPTH_LIST__
SCR_DN_ANGLE_LIST=__SCR_DN_ANGLE_LIST__
ACC_UP_DEPTH_LIST=__ACC_UP_DEPTH_LIST__
ACC_UP_ANGLE_LIST=__ACC_UP_ANGLE_LIST__
ACC_DN_DEPTH_LIST=__ACC_DN_DEPTH_LIST__
ACC_DN_ANGLE_LIST=__ACC_DN_ANGLE_LIST__
ACCEL_OFF_LIST=__ACCEL_OFF_LIST__

export H=__H__
export X_LEFT_M=__X_LEFT_M__
export X_RIGHT_M=__X_RIGHT_M__
export X_RIGHT_PHYS_M=__X_RIGHT_PHYS_M__
export YBOX_M=__YBOX_M__
export ACCEL_OFF_Y_M=__ACCEL_OFF_Y_M__
export SCREEN_OFF_Y_M=__SCREEN_OFF_Y_M__

export VS_V=__VS_V__
export VA_V=__VA_V__
export SAMPLE_V=__SAMPLE_V__
export TUBE_WALL_T_M=__TUBE_WALL_T_M__
export SAMPLE_PLATE_T_M=__SAMPLE_PLATE_T_M__
export SAMPLE_PROFILE_X_M=__SAMPLE_PROFILE_X_M__

export RUN_SOLVE=__RUN_SOLVE__
export WRITE_PNG=__WRITE_PNG__
export PNG_NAME="__PNG_NAME__"
export CLOSE_PAD_M=__CLOSE_PAD_M__
export CLOSE_BASE_PX=__CLOSE_BASE_PX__
export PLOT_FONT_PX=__PLOT_FONT_PX__

export ENABLE_IONS=__ENABLE_IONS__       # force field-only / ions for Stage-1
export PLASMA_NI_M3=__PLASMA_NI_M3__
export PLASMA_TE_EV=__PLASMA_TE_EV__
export ION_J_SCALE=__ION_J_SCALE__
export SC_FACTOR=__SC_FACTOR__
export RESULTS_DIR="__RESULTS_DIR__"
export RUN_PREFIX="__RUN_PREFIX__"
export RUN_STAMP="__RUN_STAMP__"
export PLOT_PROFILE=__PLOT_PROFILE__
export PLOT_PROFILE_Y="__PLOT_PROFILE_Y__"
export PLOT_PROFILE_X_UNITS="__PLOT_PROFILE_X_UNITS__"
export PLOT_PROFILE_EXT="__PLOT_PROFILE_EXT__"

export TRACE_ENABLE=__TRACE_ENABLE__
export VERT_PAD_M=__VERT_PAD_M__
export VERT_XLEN_M=__VERT_XLEN_M__
export VERT_YPAD_M=__VERT_YPAD_M__
export VERT_STEP_M=__VERT_STEP_M__
export VERT_THRESH_DEG=__VERT_THRESH_DEG__
EOF_ENV

# Pre-escape lists
AP_R="$(arrlit_escaped AP_RAD_LIST)"
GR_T="$(arrlit_escaped GRID_T_LIST)"
GAP_L="$(arrlit_escaped GAP_LIST)"
SU_D="$(arrlit_escaped SCR_UP_DEPTH_LIST)"
SU_A="$(arrlit_escaped SCR_UP_ANGLE_LIST)"
SD_D="$(arrlit_escaped SCR_DN_DEPTH_LIST)"
SD_A="$(arrlit_escaped SCR_DN_ANGLE_LIST)"
AU_D="$(arrlit_escaped ACC_UP_DEPTH_LIST)"
AU_A="$(arrlit_escaped ACC_UP_ANGLE_LIST)"
AD_D="$(arrlit_escaped ACC_DN_DEPTH_LIST)"
AD_A="$(arrlit_escaped ACC_DN_ANGLE_LIST)"
AO_L="$(arrlit_escaped ACCEL_OFF_LIST)"

# Pre-escape scalars
E_H="$(sed_escape "$H")"
E_XL="$(sed_escape "$X_LEFT_M")"
E_XR="$(sed_escape "$X_RIGHT_M")"
E_XRP="$(sed_escape "$X_RIGHT_PHYS_M")"
E_YB="$(sed_escape "$YBOX_M")"
E_OFF="$(sed_escape "$ACCEL_OFF_Y_M")"
E_SOFF="$(sed_escape "$SCREEN_OFF_Y_M")"
E_VS="$(sed_escape "$VS_V")"
E_VA="$(sed_escape "$VA_V")"
E_SAMP="$(sed_escape "$SAMPLE_V")"
E_TW="$(sed_escape "$TUBE_WALL_T_M")"
E_EP="$(sed_escape "$SAMPLE_PLATE_T_M")"
E_SPX="$(sed_escape "${SAMPLE_PROFILE_X_M:-nan}")"
E_RUNSOLVE="$(sed_escape "$RUN_SOLVE")"
E_WRITEPNG="$(sed_escape "$WRITE_PNG")"
E_PNGNAME="$(sed_escape "$PNG_NAME")"
E_CP="$(sed_escape "$CLOSE_PAD_M")"
E_CB="$(sed_escape "$CLOSE_BASE_PX")"
E_RDIR="$(sed_escape "$RESULTS_DIR")"
E_RPFX="$(sed_escape "$RUN_PREFIX")"
E_RSTAMP="$(sed_escape "$RUN_STAMP")"
E_PPROF="$(sed_escape "$PLOT_PROFILE")"
E_PPY="$(sed_escape "$PLOT_PROFILE_Y")"
E_PPX="$(sed_escape "$PLOT_PROFILE_X_UNITS")"
E_PPE="$(sed_escape "$PLOT_PROFILE_EXT")"
E_TRACE="$(sed_escape "$TRACE_ENABLE")"
E_VPAD="$(sed_escape "$VERT_PAD_M")"
E_PFPX="$(sed_escape "$PLOT_FONT_PX")"
E_VXLEN="$(sed_escape "$VERT_XLEN_M")"
E_VYPAD="$(sed_escape "$VERT_YPAD_M")"
E_VSTEP="$(sed_escape "$VERT_STEP_M")"
E_VTH="$(sed_escape "$VERT_THRESH_DEG")"
E_ENIONS1="$(sed_escape "$ENABLE_IONS_STAGE1")"
E_PNI="$(sed_escape "$PLASMA_NI_M3")"
E_PTE="$(sed_escape "$PLASMA_TE_EV")"
E_IONJ="$(sed_escape "$ION_J_SCALE")"
E_SCF="$(sed_escape "$SC_FACTOR")"

sed -i \
  -e "s#__AP_RAD_LIST__#${AP_R}#" \
  -e "s#__GRID_T_LIST__#${GR_T}#" \
  -e "s#__GAP_LIST__#${GAP_L}#" \
  -e "s#__SCR_UP_DEPTH_LIST__#${SU_D}#" \
  -e "s#__SCR_UP_ANGLE_LIST__#${SU_A}#" \
  -e "s#__SCR_DN_DEPTH_LIST__#${SD_D}#" \
  -e "s#__SCR_DN_ANGLE_LIST__#${SD_A}#" \
  -e "s#__ACC_UP_DEPTH_LIST__#${AU_D}#" \
  -e "s#__ACC_UP_ANGLE_LIST__#${AU_A}#" \
  -e "s#__ACC_DN_DEPTH_LIST__#${AD_D}#" \
  -e "s#__ACC_DN_ANGLE_LIST__#${AD_A}#" \
  -e "s#__ACCEL_OFF_LIST__#${AO_L}#" \
  -e "s#__H__#${E_H}#" \
  -e "s#__X_LEFT_M__#${E_XL}#" \
  -e "s#__X_RIGHT_M__#${E_XR}#" \
  -e "s#__X_RIGHT_PHYS_M__#${E_XRP}#" \
  -e "s#__YBOX_M__#${E_YB}#" \
  -e "s#__ACCEL_OFF_Y_M__#${E_OFF}#" \
  -e "s#__SCREEN_OFF_Y_M__#${E_SOFF}#" \
  -e "s#__VS_V__#${E_VS}#" \
  -e "s#__VA_V__#${E_VA}#" \
  -e "s#__SAMPLE_V__#${E_SAMP}#" \
  -e "s#__TUBE_WALL_T_M__#${E_TW}#" \
  -e "s#__SAMPLE_PLATE_T_M__#${E_EP}#" \
  -e "s#__SAMPLE_PROFILE_X_M__#${E_SPX}#" \
  -e "s#__RUN_SOLVE__#${E_RUNSOLVE}#" \
  -e "s#__WRITE_PNG__#${E_WRITEPNG}#" \
  -e "s#__PNG_NAME__#${E_PNGNAME}#" \
  -e "s#__CLOSE_PAD_M__#${E_CP}#" \
  -e "s#__CLOSE_BASE_PX__#${E_CB}#" \
  -e "s#__PLOT_FONT_PX__#${E_PFPX}#" \
  -e "s#__RESULTS_DIR__#${E_RDIR}#" \
  -e "s#__RUN_PREFIX__#${E_RPFX}#" \
  -e "s#__RUN_STAMP__#${E_RSTAMP}#" \
  -e "s#__PLOT_PROFILE__#${E_PPROF}#" \
  -e "s#__PLOT_PROFILE_Y__#${E_PPY}#" \
  -e "s#__PLOT_PROFILE_X_UNITS__#${E_PPX}#" \
  -e "s#__PLOT_PROFILE_EXT__#${E_PPE}#" \
  -e "s#__TRACE_ENABLE__#${E_TRACE}#" \
  -e "s#__VERT_PAD_M__#${E_VPAD}#" \
  -e "s#__VERT_XLEN_M__#${E_VXLEN}#" \
  -e "s#__VERT_YPAD_M__#${E_VYPAD}#" \
  -e "s#__VERT_STEP_M__#${E_VSTEP}#" \
  -e "s#__VERT_THRESH_DEG__#${E_VTH}#" \
  -e "s#__ENABLE_IONS__#${E_ENIONS1}#" \
  -e "s#__PLASMA_NI_M3__#${E_PNI}#" \
  -e "s#__PLASMA_TE_EV__#${E_PTE}#" \
  -e "s#__ION_J_SCALE__#${E_IONJ}#" \
  -e "s#__SC_FACTOR__#${E_SCF}#" \
  "$PIPE_ROOT/sweep_lists.env"

# Source to compute cardinalities
# shellcheck source=/dev/null
source "$PIPE_ROOT/sweep_lists.env"

len(){ local -n A="$1"; echo "${#A[@]}"; }
nR=$(len AP_RAD_LIST); nT=$(len GRID_T_LIST); nG=$(len GAP_LIST)
nSUd=$(len SCR_UP_DEPTH_LIST); nSUa=$(len SCR_UP_ANGLE_LIST)
nSDd=$(len SCR_DN_DEPTH_LIST); nSDa=$(len SCR_DN_ANGLE_LIST)
nAUd=$(len ACC_UP_DEPTH_LIST); nAUa=$(len ACC_UP_ANGLE_LIST)
nADd=$(len ACC_DN_DEPTH_LIST); nADa=$(len ACC_DN_ANGLE_LIST)
nAO=$(len ACCEL_OFF_LIST)
TOTAL=$(( nR*nT*nG*nSUd*nSUa*nSDd*nSDa*nAUd*nAUa*nADd*nADa*nAO ))
echo "[orchestrator] Stage-1 total combinations = $TOTAL"

# =====================================================================================
# Stage-1 array script
# =====================================================================================
cat > "$PIPE_ROOT/sweep_two_grid.slurm" <<'EOF_SWEEP'
#!/bin/bash
#SBATCH -J two_grid_sweep
#SBATCH -n 1
#SBATCH --cpus-per-task=__SWEEP_CPUS__
#SBATCH --mem=__SWEEP_MEM__
#SBATCH -t __SWEEP_TIME__
#SBATCH -o %x.%A_%a.out
#SBATCH -e %x.%A_%a.err
__SWEEP_PART_LINE__

set -euo pipefail
umask 002
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
[ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash || true

# shellcheck source=/dev/null
source "__PIPE_ROOT__/sweep_lists.env"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

pick(){ local idx="$1"; shift; local arr=( "$@" ); echo "${arr[$idx]}"; }

idx=${SLURM_ARRAY_TASK_ID:?}
idx=$(( idx + ${ARRAY_BASE:-0} ))
i=$idx
iR=$(( i % ${#AP_RAD_LIST[@]} ));  i=$(( i / ${#AP_RAD_LIST[@]} ))
iT=$(( i % ${#GRID_T_LIST[@]} ));  i=$(( i / ${#GRID_T_LIST[@]} ))
iG=$(( i % ${#GAP_LIST[@]} ));     i=$(( i / ${#GAP_LIST[@]} ))
iSUd=$(( i % ${#SCR_UP_DEPTH_LIST[@]} )); i=$(( i / ${#SCR_UP_DEPTH_LIST[@]} ))
iSUa=$(( i % ${#SCR_UP_ANGLE_LIST[@]} )); i=$(( i / ${#SCR_UP_ANGLE_LIST[@]} ))
iSDd=$(( i % ${#SCR_DN_DEPTH_LIST[@]} )); i=$(( i / ${#SCR_DN_DEPTH_LIST[@]} ))
iSDa=$(( i % ${#SCR_DN_ANGLE_LIST[@]} )); i=$(( i / ${#SCR_DN_ANGLE_LIST[@]} ))
iAUd=$(( i % ${#ACC_UP_DEPTH_LIST[@]} )); i=$(( i / ${#ACC_UP_DEPTH_LIST[@]} ))
iAUa=$(( i % ${#ACC_UP_ANGLE_LIST[@]} )); i=$(( i / ${#ACC_UP_ANGLE_LIST[@]} ))
iADd=$(( i % ${#ACC_DN_DEPTH_LIST[@]} )); i=$(( i / ${#ACC_DN_DEPTH_LIST[@]} ))
iADa=$(( i % ${#ACC_DN_ANGLE_LIST[@]} )); i=$(( i / ${#ACC_DN_ANGLE_LIST[@]} ))
iAO=$(( i % ${#ACCEL_OFF_LIST[@]} ))

export AP_RAD_M=$(pick "$iR"  "${AP_RAD_LIST[@]}")
export GRID_T_M=$(pick "$iT"  "${GRID_T_LIST[@]}")
export GAP_M=$(pick "$iG"     "${GAP_LIST[@]}")
export SCR_UP_DEPTH_M=$(pick "$iSUd" "${SCR_UP_DEPTH_LIST[@]}")
export SCR_UP_ANGLE_DEG=$(pick "$iSUa" "${SCR_UP_ANGLE_LIST[@]}")
export SCR_DN_DEPTH_M=$(pick "$iSDd" "${SCR_DN_DEPTH_LIST[@]}")
export SCR_DN_ANGLE_DEG=$(pick "$iSDa" "${SCR_DN_ANGLE_LIST[@]}")
export ACC_UP_DEPTH_M=$(pick "$iAUd" "${ACC_UP_DEPTH_LIST[@]}")
export ACC_UP_ANGLE_DEG=$(pick "$iAUa" "${ACC_UP_ANGLE_LIST[@]}")
export ACC_DN_DEPTH_M=$(pick "$iADd" "${ACC_DN_DEPTH_LIST[@]}")
export ACC_DN_ANGLE_DEG=$(pick "$iADa" "${ACC_DN_ANGLE_LIST[@]}")
export ACCEL_OFF_Y_M=$(pick "$iAO" "${ACCEL_OFF_LIST[@]}")

sanitize(){ echo "$1" | sed -e 's/\./p/g' -e 's/-/m/g'; }
rt=$(sanitize "$AP_RAD_M"); gt=$(sanitize "$GRID_T_M"); gp=$(sanitize "$GAP_M")
suD=$(sanitize "$SCR_UP_DEPTH_M"); suA=$(sanitize "$SCR_UP_ANGLE_DEG")
sdD=$(sanitize "$SCR_DN_DEPTH_M"); sdA=$(sanitize "$SCR_DN_ANGLE_DEG")
auD=$(sanitize "$ACC_UP_DEPTH_M");  auA=$(sanitize "$ACC_UP_ANGLE_DEG")
adD=$(sanitize "$ACC_DN_DEPTH_M");  adA=$(sanitize "$ACC_DN_ANGLE_DEG")
ao=$(sanitize "$ACCEL_OFF_Y_M");    so=$(sanitize "$SCREEN_OFF_Y_M")

export RUN_TAG="i${idx}_rad${rt}_t${gt}_gap${gp}_AO${ao}_SO${so}_SU(${suD},${suA})_SD(${sdD},${sdA})_AU(${auD},${auA})_AD(${adD},${adA})"

echo "[$(date)] running..."
srun --ntasks=1 "__PIPE_ROOT__/bin/two_grid_2d"

# Optional: plot sample profile automatically
if [[ "${PLOT_PROFILE:-1}" != "0" ]]; then
  OUTDIR="${RESULTS_DIR}/${RUN_PREFIX}_${RUN_STAMP}"
  if [[ -n "${RUN_TAG:-}" ]]; then OUTDIR="${OUTDIR}_${RUN_TAG}"; fi
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then OUTDIR="${OUTDIR}_j${SLURM_JOB_ID}"; fi

  PROF_JSON="${OUTDIR}/sample_diameter_profile.json"
  if [[ ! -f "${PROF_JSON}" ]]; then
    PROF_JSON="${OUTDIR}/sample_radial_profile.json"
  fi

  if [[ -f "${PROF_JSON}" ]]; then
    ext="${PLOT_PROFILE_EXT:-png}"
    python3 "__PIPE_ROOT__/plot_runlog_compact.py" sample-profile \
      --json "${PROF_JSON}" \
      --out "${OUTDIR}/sample_profile.${ext}" \
      --y "${PLOT_PROFILE_Y:-I_Apm}" \
      --x-units "${PLOT_PROFILE_X_UNITS:-mm}" \
      || echo "[plot] WARNING: failed to generate sample_profile.${ext}"
  fi
fi
EOF_SWEEP

# Robust replace for sweep resources/paths
E_CPUS="$(sed_escape "$SWEEP_CPUS")"
E_MEM="$(sed_escape "$SWEEP_MEM")"
E_TIME="$(sed_escape "$SWEEP_TIME")"
E_ROOT="$(sed_escape "$PIPE_ROOT")"
sed -i \
  -e "s#__SWEEP_CPUS__#${E_CPUS}#" \
  -e "s#__SWEEP_MEM__#${E_MEM}#" \
  -e "s#__SWEEP_TIME__#${E_TIME}#" \
  -e "s#__PIPE_ROOT__#${E_ROOT}#" \
  "$PIPE_ROOT/sweep_two_grid.slurm"

if [[ -n "${SWEEP_PART:-}" ]]; then
  sed -i -e "s#__SWEEP_PART_LINE__#\#SBATCH --partition=$(sed_escape "$SWEEP_PART")#" "$PIPE_ROOT/sweep_two_grid.slurm"
else
  sed -i -e "s#__SWEEP_PART_LINE__##" "$PIPE_ROOT/sweep_two_grid.slurm"
fi

# =====================================================================================
# Reducer (Python + wrapper Slurm)
# =====================================================================================
cat > "$PIPE_ROOT/reduce_best.py" <<'PY'
#!/usr/bin/env python3
import sys, json, os, time, shutil, statistics
from pathlib import Path

RESULTS = Path(os.environ.get("RESULTS_DIR", "./results"))
METRIC  = os.environ.get("VERT_METRIC","mean_tan")           # legacy verticality.* key
COL_METRIC = os.environ.get("VERT_METRIC_COL","col_angle_max_deg")  # column_verticality.* key
ASCEND  = True
TOPK    = int(os.environ.get("TOPK","10"))
EXPLICIT  = os.environ.get("OPT_KNOBS","").strip()
MAX_KNOBS = int(os.environ.get("MAX_KNOBS","4"))
DELTA_FLOORS = json.loads(os.environ.get("DELTA_FLOORS_JSON","{}"))

def loadj(p):
    try: return json.loads(p.read_text())
    except Exception: return None

def score_run(run: Path):
    """
    Score a run using beam diagnostics:
      1) |Y_MEAN_PRED_500MM_M| closest to 0
      2) DIVERGENCE_ANGLE_DEG smallest
      3) I_AG_OUT_A largest

    Disqualify if LOST_TO_SIDEWALLS is true.
    """
    bm = loadj(run/"beam_metrics.json") or {}
    coll = bm.get("collimation") or {}
    if coll.get("lost_to_sidewalls") in (True, 1, "1", "true"):
        return None

    defl = bm.get("deflection") or {}
    cur  = bm.get("currents") or {}

    y_pred = defl.get("y_mean_pred_500mm_m")
    div_deg = bm.get("DIVERGENCE_ANGLE_DEG")
    if div_deg is None:
        div_deg = bm.get("divergence_angle_deg")
    I_ag = cur.get("I_ag_out_A")

    try:
        y_pred = float(y_pred)
        div_deg = float(div_deg)
        I_ag = float(I_ag)
    except Exception:
        return None

    abs_y = abs(y_pred)
    score = (abs_y, div_deg, -I_ag)
    info = {
        "y_mean_pred_500mm_m": y_pred,
        "abs_y_mean_pred_500mm_m": abs_y,
        "divergence_angle_deg": div_deg,
        "I_ag_out_A": I_ag,
        "lost_to_sidewalls": bool(coll.get("lost_to_sidewalls", False)),
    }
    return score, info


def pick_pngs(d:Path):
    mains=list(d.glob("*.png"))
    main=next((p for p in mains if p.name=="two_grid_geom.png"), None)
    if not main:
        main=next((p for p in mains if not p.stem.endswith("_closeup")), None)
    close=next((p for p in mains if p.stem.endswith("_closeup")), None)
    return main,close

def extract_geom(meta):
    out={}
    if not meta: return out
    scr=meta.get("screen_chamfer") or {}
    acc=meta.get("accel_chamfer") or {}
    out.update({k:meta.get(k) for k in [
        "AP_RAD_M","GRID_T_M","GAP_M","ACCEL_OFF_Y_M","SCREEN_OFF_Y_M","H","X_LEFT_M","X_RIGHT_M","YBOX_M",
        "VS_V","VA_V","SAMPLE_V"
    ] if k in meta})
    out.update({
        "SCR_UP_DEPTH_M":scr.get("up_depth"), "SCR_UP_ANGLE_DEG":scr.get("up_angle_deg"),
        "SCR_DN_DEPTH_M":scr.get("dn_depth"), "SCR_DN_ANGLE_DEG":scr.get("dn_angle_deg"),
        "ACC_UP_DEPTH_M":acc.get("up_depth"), "ACC_UP_ANGLE_DEG":acc.get("up_angle_deg"),
        "ACC_DN_DEPTH_M":acc.get("dn_depth"), "ACC_DN_ANGLE_DEG":acc.get("dn_angle_deg"),
    })
    return out

def collect_rows(root: Path):
    rows=[]
    for meta in root.rglob("meta.json"):
        run = meta.parent
        scored = score_run(run)
        if not scored:
            continue
        score, info = scored
        y = info.get("abs_y_mean_pred_500mm_m")
        if y is None:
            continue

        m   = loadj(meta) or {}
        scr = (m.get("screen_chamfer") or {})
        acc = (m.get("accel_chamfer") or {})
        feat = {
          "AP_RAD_M": m.get("AP_RAD_M"), "GRID_T_M": m.get("GRID_T_M"),
          "GAP_M": m.get("GAP_M"), "ACCEL_OFF_Y_M": m.get("ACCEL_OFF_Y_M"), "SCREEN_OFF_Y_M": m.get("SCREEN_OFF_Y_M"),
          "SCR_UP_DEPTH_M": scr.get("up_depth"), "SCR_UP_ANGLE_DEG": scr.get("up_angle_deg"),
          "SCR_DN_DEPTH_M": scr.get("dn_depth"), "SCR_DN_ANGLE_DEG": scr.get("dn_angle_deg"),
          "ACC_UP_DEPTH_M": acc.get("up_depth"), "ACC_UP_ANGLE_DEG": acc.get("up_angle_deg"),
          "ACC_DN_DEPTH_M": acc.get("dn_depth"), "ACC_DN_ANGLE_DEG": acc.get("dn_angle_deg"),
        }
        if any(vv is None for vv in feat.values()):
            continue
        rows.append((feat, float(y), run))
    return rows


def suggest_knobs(rows, best_geom, top):
    if EXPLICIT:
        knob_list=[k.strip() for k in EXPLICIT.split(",") if k.strip()]
    else:
        try:
            import numpy as np
            cols=list(rows[0][0].keys())
            X=np.array([[r[0][c] for c in cols] for r in rows], float)
            y=np.array([r[1] for r in rows], float)
            mu=X.mean(0); sig=X.std(0)+1e-12
            Xz=(X-mu)/sig
            lam=1e-6
            w=np.linalg.solve(Xz.T@Xz + lam*np.eye(Xz.shape[1]), Xz.T@y)
            imp=np.abs(w)
            rank=np.argsort(-imp)
            knob_list=[cols[i] for i in rank[:MAX_KNOBS]]
        except Exception:
            cols=[k for k in top[0].keys() if k.endswith("_M") or k.endswith("_DEG")]
            var={c: (statistics.pstdev([float(r[c]) for r in top if r.get(c) is not None]) if sum(1 for r in top if r.get(c) is not None)>1 else 0.0) for c in cols}
            knob_list=[k for k,_ in sorted(var.items(), key=lambda kv: -kv[1])[:MAX_KNOBS]]

    # physics guardrail: angle only if its depth > 0
    pairs=[("SCR_UP_ANGLE_DEG","SCR_UP_DEPTH_M"),
           ("SCR_DN_ANGLE_DEG","SCR_DN_DEPTH_M"),
           ("ACC_UP_ANGLE_DEG","ACC_UP_DEPTH_M"),
           ("ACC_DN_ANGLE_DEG","ACC_DN_DEPTH_M")]
    for ang,dep in pairs:
        if ang in knob_list and (best_geom.get(dep,0.0)<=0.0):
            knob_list.remove(ang)
            if dep not in knob_list:
                knob_list.append(dep)
    return knob_list

def suggest_deltas(knobs, top, floors):
    deltas={}
    for k in knobs:
        vals=[t.get(k) for t in top if t.get(k) is not None]
        spread=(max(vals)-min(vals))/4.0 if len(vals)>=2 else 0.0
        floor=floors.get(k, 0.0)
        deltas[k]=max(spread, floor)
    return deltas

def main():
    items=[]
    for meta in RESULTS.rglob("meta.json"):
        run=meta.parent
        scored = score_run(run)
        if not scored:
            continue
        s, info = scored
        items.append({"run":run, "meta":loadj(meta), "beam":info, "score":s})
    if not items:
        print("[reduce] no runs with metrics found", file=sys.stderr); sys.exit(2)

    items.sort(key=lambda r: r["score"])
    stamp=time.strftime("%Y%m%d_%H%M%S")
    outdir=RESULTS/"_stage1_best"/stamp
    outdir.mkdir(parents=True, exist_ok=True)

    top=[]
    for r in items[:min(TOPK,len(items))]:
        row={"run_dir":str(r["run"]),
             "score_keys":["abs_y_mean_pred_500mm_m","divergence_angle_deg","-I_ag_out_A"],
             "score_values":[r["beam"]["abs_y_mean_pred_500mm_m"],
                             r["beam"]["divergence_angle_deg"],
                             -r["beam"]["I_ag_out_A"]],
             "beam_metrics":r["beam"]}
        row.update(extract_geom(r["meta"]))
        top.append(row)
    (outdir/"topk.json").write_text(json.dumps({"ranking":"abs(y_mean_pred_500mm_m), divergence_angle_deg, -I_ag_out_A",
                                                "topk":top}, indent=2))

    best=items[0]
    geom=extract_geom(best["meta"])
    (outdir/"best.json").write_text(json.dumps({
        "ranking":"abs(y_mean_pred_500mm_m), divergence_angle_deg, -I_ag_out_A",
        "score_values":[best["beam"]["abs_y_mean_pred_500mm_m"],
                        best["beam"]["divergence_angle_deg"],
                        -best["beam"]["I_ag_out_A"]],
        "run_dir":str(best["run"]), "beam_metrics":best["beam"], "geom":geom
    }, indent=2))

    seed=(outdir/"seed_stage2.env")
    with seed.open("w") as f:
        for k,v in geom.items():
            if v is not None: f.write(f"export {k}={v}\n")

    main_png,close_png=pick_pngs(best["run"])
    if main_png: shutil.copy2(main_png, outdir/main_png.name)
    if close_png: shutil.copy2(close_png, outdir/close_png.name)

    latest=RESULTS/"_stage1_best"/"latest"
    try:
        if latest.exists() or latest.is_symlink(): latest.unlink()
        latest.symlink_to(outdir.name)
    except Exception: pass
    (RESULTS/"_stage1_best"/"BEST_PATH.txt").write_text(str(outdir)+"\n")

    rows=collect_rows(RESULTS)
    top_full=json.loads((outdir/"topk.json").read_text())["topk"]
    knobs=suggest_knobs(rows, geom, top_full if top_full else [])
    deltas=suggest_deltas(knobs, top_full if top_full else [], json.loads(os.environ.get("DELTA_FLOORS_JSON","{}")))

    knobs_env = RESULTS/"_stage1_best"/"latest"/"knobs_stage2.env"
    with knobs_env.open("w") as f:
        f.write("export OPT_KNOBS=" + ",".join(knobs) + "\n")
        for k,d in deltas.items():
            f.write(f"export DELTA_{k}={d}\n")

    print(f"[reduce] best: {best['run']}  score={best['score']}")
    print(f"[reduce] seed: {seed}")
    print(f"[reduce] knobs: {knobs} with deltas {deltas}")
    print(f"[reduce] out : {outdir}")

if __name__=="__main__": main()
PY
chmod +x "$PIPE_ROOT/reduce_best.py"

cat > "$PIPE_ROOT/reduce_best.slurm" <<'EOF_REDUCE'
#!/bin/bash
#SBATCH -J two_grid_reduce
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH -t 00:05:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# SBATCH --partition=medium
set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
[ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash || true
export RESULTS_DIR="__RESULTS_DIR__"
export VERT_METRIC="__VERT_METRIC__"
export TOPK="__TOPK__"
export OPT_KNOBS="__OPT_KNOBS__"
export MAX_KNOBS="__MAX_KNOBS__"
export DELTA_FLOORS_JSON="$(echo '__DELTA_FLOORS_B64__' | base64 -d)"
python3 "__PIPE_ROOT__/reduce_best.py"
EOF_REDUCE

E_RDIR="$(sed_escape "$RESULTS_DIR")"
E_ROOT="$(sed_escape "$PIPE_ROOT")"
E_MET="$(sed_escape "$VERT_METRIC")"
E_TOP="$(sed_escape "$TOPK")"
E_KNB="$(sed_escape "$OPT_KNOBS")"
E_MKN="$(sed_escape "$MAX_KNOBS")"
sed -i \
  -e "s#__RESULTS_DIR__#${E_RDIR}#" \
  -e "s#__PIPE_ROOT__#${E_ROOT}#" \
  -e "s#__VERT_METRIC__#${E_MET}#" \
  -e "s#__TOPK__#${E_TOP}#" \
  -e "s#__OPT_KNOBS__#${E_KNB}#" \
  -e "s#__MAX_KNOBS__#${E_MKN}#" \
  -e "s#__DELTA_FLOORS_B64__#${DF_B64}#" \
  "$PIPE_ROOT/reduce_best.slurm"



# =====================================================================================
# Stage-2 optimizer controller (Python)
# =====================================================================================

if [[ "$RUN_STAGE2" == "1" ]]; then
cat > "$PIPE_ROOT/stage2_opt.py" <<'PY'
#!/usr/bin/env python3
import os, json, time, math, subprocess, random
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

# --- Environment & basic paths ---
RESULTS = Path(os.environ.get("RESULTS_DIR", "./results"))
PIPE_ROOT = Path(os.environ.get("PIPE_ROOT", "."))
RUN_PREFIX = os.environ.get("RUN_PREFIX", "two_grid_refine")
JOBID = os.environ.get("SLURM_JOB_ID", "")

# Reuse RUN_STAMP if provided; otherwise create one for this Stage-2 campaign
stamp = os.environ.get("RUN_STAMP", "").strip()
if not stamp:
    stamp = time.strftime("%Y%m%d_%H%M%S")
    os.environ["RUN_STAMP"] = stamp

# Stage-2 scoring:
#   1) I_ag_out_A (higher is better)
#   2) divergence_angle_deg (lower is better)
#   3) |y_mean_pred_500mm_m| (lower is better)


MAX_ITERS       = int(os.environ.get("ST2_MAX_ITERS", "10"))
NO_IMPROVE_MAX  = int(os.environ.get("ST2_NO_IMPROVE", "3"))
BATCH_SIZE      = int(os.environ.get("ST2_BATCH_SIZE", "3"))   # concurrent evals
RANDOM_SEED     = int(os.environ.get("ST2_SEED", "0") or "0")

random.seed(RANDOM_SEED)

# Aperture packing: approximate max number of beamlets vs aperture radius.
# Given as [radius_m, max_N_aps] translated from your mm table:
# [1mm:2221, 2mm:769, 3mm:379, 4mm:223, 5mm:151, 6mm:109, 7mm:73, 8mm:55]
AP_PACK = [
    (1.0e-3, 2221.0),
    (2.0e-3,  769.0),
    (3.0e-3,  379.0),
    (4.0e-3,  223.0),
    (5.0e-3,  151.0),
    (6.0e-3,  109.0),
    (7.0e-3,   73.0),
    (8.0e-3,   55.0),
]

def estimate_beamlet_count(ap_rad_m: float) -> float:
    """Interpolate max number of apertures as a function of aperture radius."""
    if ap_rad_m is None or ap_rad_m <= 0.0:
        return AP_PACK[0][1]
    # Clamp below/above table
    if ap_rad_m <= AP_PACK[0][0]:
        return AP_PACK[0][1]
    if ap_rad_m >= AP_PACK[-1][0]:
        return AP_PACK[-1][1]
    # Piecewise linear interpolation
    for (r0, n0), (r1, n1) in zip(AP_PACK[:-1], AP_PACK[1:]):
        if r0 <= ap_rad_m <= r1:
            t = (ap_rad_m - r0) / (r1 - r0) if r1 > r0 else 0.0
            return n0 + t * (n1 - n0)
    return AP_PACK[-1][1]

def load_json(p: Path) -> Any:
    try:
        return json.loads(p.read_text())
    except Exception:
        return None

PLOT_PROFILE = os.environ.get("PLOT_PROFILE", "1").strip() != "0"
PLOT_PROFILE_Y = os.environ.get("PLOT_PROFILE_Y", "I_Apm").strip() or "I_Apm"
PLOT_PROFILE_X_UNITS = os.environ.get("PLOT_PROFILE_X_UNITS", "mm").strip() or "mm"
PLOT_PROFILE_EXT = os.environ.get("PLOT_PROFILE_EXT", "png").strip() or "png"

# Stage-2 scoring configuration
SCORE_MODE = os.environ.get("ST2_SCORE_MODE", "lex").strip().lower()
I_FIELD = os.environ.get("ST2_I_FIELD", "I_ag_out_A").strip() or "I_ag_out_A"
DIV_FIELD = os.environ.get("ST2_DIV_FIELD", "divergence_angle_deg").strip() or "divergence_angle_deg"
Y_FIELD = os.environ.get("ST2_Y_FIELD", "y_mean_pred_500mm_m").strip() or "y_mean_pred_500mm_m"
W_I = float(os.environ.get("ST2_W_I", "1.0"))
W_DIV = float(os.environ.get("ST2_W_DIV", "1.0"))
W_Y = float(os.environ.get("ST2_W_Y", "1.0"))
I_REF = float(os.environ.get("ST2_I_REF_A", "1.0"))
DIV_REF = float(os.environ.get("ST2_DIV_REF_DEG", "1.0"))
Y_REF = float(os.environ.get("ST2_Y_REF_M", "1.0e-3"))
PLOT_SCRIPT = PIPE_ROOT / "plot_runlog_compact.py"

def plot_profile(run_dir: Path) -> None:
    if not PLOT_PROFILE:
        return
    if not PLOT_SCRIPT.exists():
        return
    prof = run_dir / "sample_diameter_profile.json"
    if not prof.exists():
        prof = run_dir / "sample_radial_profile.json"
    if not prof.exists():
        return
    out_html = run_dir / f"sample_profile.{PLOT_PROFILE_EXT}"
    cmd = [
        "python3",
        str(PLOT_SCRIPT),
        "sample-profile",
        "--json",
        str(prof),
        "--out",
        str(out_html),
        "--y",
        PLOT_PROFILE_Y,
        "--x-units",
        PLOT_PROFILE_X_UNITS,
    ]
    try:
        subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass
CHAMFER_PAIRS = [
    ("SCR_UP_DEPTH_M", "SCR_UP_ANGLE_DEG"),
    ("SCR_DN_DEPTH_M", "SCR_DN_ANGLE_DEG"),
    ("ACC_UP_DEPTH_M", "ACC_UP_ANGLE_DEG"),
    ("ACC_DN_DEPTH_M", "ACC_DN_ANGLE_DEG"),
]
EPS_ZERO = 1.0e-12

def is_zero(v: float) -> bool:
    return abs(v) <= EPS_ZERO

def is_valid_chamfer(vals: Dict[str, float]) -> bool:
    for dkey, akey in CHAMFER_PAIRS:
        if dkey not in vals or akey not in vals:
            continue
        d = float(vals[dkey])
        a = float(vals[akey])
        if is_zero(d) ^ is_zero(a):
            return False
    return True

def closest_index(vals: List[float], want_zero: bool) -> int:
    idx = None
    best = None
    for i, v in enumerate(vals):
        v = float(v)
        if want_zero and not is_zero(v):
            continue
        if not want_zero and is_zero(v):
            continue
        score = abs(v)
        if best is None or score < best:
            best = score
            idx = i
    return -1 if idx is None else idx

def to_float(v: Any) -> Any:
    try:
        f = float(v)
    except Exception:
        return None
    return f if math.isfinite(f) else None

def fmt_score(score: Any) -> str:
    if isinstance(score, (tuple, list)):
        parts = []
        for x in score:
            if isinstance(x, (int, float)) and math.isfinite(x):
                parts.append(f"{x:.6g}")
            else:
                parts.append(str(x))
        return "(" + ", ".join(parts) + ")"
    if isinstance(score, (int, float)) and math.isfinite(score):
        return f"{score:.6g}"
    return str(score)

def score_key(score: Any) -> Tuple[float, ...]:
    if isinstance(score, (tuple, list)):
        try:
            return tuple(float(x) for x in score)
        except Exception:
            return (math.inf, math.inf, math.inf)
    try:
        return (float(score),)
    except Exception:
        return (math.inf,)

def score_to_list(score: Any) -> List[float]:
    if isinstance(score, (tuple, list)):
        return [float(x) for x in score]
    return [float(score)]


def score_run(run_dir: Path) -> Tuple[Any, Dict[str, Any]]:
    """
    Score a run using beam diagnostics from beam_metrics.json.
    """
    bm = load_json(run_dir / "beam_metrics.json") or {}
    coll = (bm.get("collimation") or {})
    if coll.get("lost_to_sidewalls") in (True, 1, "1", "true"):
        return (math.inf, math.inf, math.inf), {"reason": "lost_to_sidewalls"}

    defl = (bm.get("deflection") or {})
    cur = (bm.get("currents") or {})

    def pick_field(field_spec: str) -> Tuple[Optional[float], Optional[str]]:
        if not field_spec:
            return None, None
        for field in field_spec.split("|"):
            field = field.strip()
            if not field:
                continue
            if "." in field:
                sect, key = field.split(".", 1)
                if sect == "deflection":
                    val = defl.get(key)
                elif sect == "currents":
                    val = cur.get(key)
                else:
                    val = bm.get(field)
            elif field in defl:
                val = defl.get(field)
            elif field in cur:
                val = cur.get(field)
            elif field in bm:
                val = bm.get(field)
            else:
                val = None
            f = to_float(val)
            if f is not None:
                return f, field
        return None, None

    i_val, i_used = pick_field(I_FIELD)
    div_val, div_used = pick_field(DIV_FIELD)
    y_val, y_used = pick_field(Y_FIELD)

    if i_val is None or div_val is None or y_val is None:
        missing = []
        if i_val is None:
            missing.append(f"I:{I_FIELD}")
        if div_val is None:
            missing.append(f"DIV:{DIV_FIELD}")
        if y_val is None:
            missing.append(f"Y:{Y_FIELD}")
        return (math.inf, math.inf, math.inf), {"reason": "missing_metrics", "missing": missing}

    y_term = abs(y_val) if (os.environ.get("ST2_Y_ABS", "1").strip() != "0") else y_val

    info = {
        "I_metric": i_val,
        "DIV_metric": div_val,
        "Y_metric": y_val,
        "I_used": i_used,
        "DIV_used": div_used,
        "Y_used": y_used,
        "Y_term": y_term,
        "I_field": I_FIELD,
        "DIV_field": DIV_FIELD,
        "Y_field": Y_FIELD,
        "lost_to_sidewalls": bool(coll.get("lost_to_sidewalls", False)),
    }

    if SCORE_MODE == "weighted":
        i_ref = I_REF if I_REF > 0.0 else 1.0
        d_ref = DIV_REF if DIV_REF > 0.0 else 1.0
        y_ref = Y_REF if Y_REF > 0.0 else 1.0
        score = (-W_I * (i_val / i_ref)) + (W_DIV * (div_val / d_ref)) + (W_Y * (y_term / y_ref))
        info.update({
            "score_mode": "weighted",
            "weights": {"I": W_I, "DIV": W_DIV, "Y": W_Y},
            "refs": {"I": i_ref, "DIV": d_ref, "Y": y_ref},
        })
        return score, info

    score = (-i_val, div_val, y_term)
    info.update({"score_mode": "lex"})
    return score, info


def parse_refine_dump() -> Dict[str, List[float]]:
    """
    REFINE_DUMP comes from refine_stage2.slurm and encodes lists for each knob.

    Example:
      {"SCR_DN_DEPTH_M":"0.0003 0.0005 0.0007","SCR_DN_ANGLE_DEG":"-4 0 4"}
    """
    raw = os.environ.get("REFINE_DUMP", "{}").strip()
    if not raw or raw == "{}":
        return {}

    s = raw.strip("{}").strip()
    if not s:
        return {}

    d: Dict[str, List[float]] = {}
    # crude parser but sufficient for our simple format
    for pair in s.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if ":" not in pair:
            continue
        k, v = pair.split(":", 1)
        k = k.strip().strip('"')
        v = v.strip().strip('"')
        if not k:
            continue
        vals = []
        for token in v.split():
            try:
                vals.append(float(token))
            except Exception:
                pass
        if vals:
            d[k] = vals
    return d


def make_run_dir(tag: str) -> Path:
    """Predict the results directory name based on two_grid_2d's convention."""
    name = f"{RUN_PREFIX}_{stamp}"
    if tag:
        name += f"_{tag}"
    if JOBID:
        name += f"_j{JOBID}"
    return RESULTS / name


def launch_candidate(tag: str, knob_vals: Dict[str, float]) -> subprocess.Popen:
    """
    Launch a single two_grid_2d evaluation via srun, with env overrides for knobs.

    Returns a Popen handle that we will wait() on.
    """
    env = os.environ.copy()
    env["RUN_TAG"] = tag
    # Ensure we stay in Stage-2 / ions mode
    env["ENABLE_IONS"] = env.get("ENABLE_IONS", "1")
    # Override knobs
    for k, val in knob_vals.items():
        env[k] = f"{val:.12g}"

    rdir = make_run_dir(tag)
    print(f"[stage2] launching {tag} -> {rdir}", flush=True)

    cmd = ["srun", "--ntasks=1", str(PIPE_ROOT / "bin" / "two_grid_2d")]
    return subprocess.Popen(cmd, env=env)


def evaluate_candidate(tag: str, knob_vals: Dict[str, float]) -> Tuple[Tuple[float, float, float], Dict[str, Any]]:
    """
    Launch and wait for a candidate, then read its metrics and score it.
    """
    if not is_valid_chamfer(knob_vals):
        return (math.inf, math.inf, math.inf), {"reason": "invalid_chamfer", "tag": tag, "knobs": knob_vals}
    proc = launch_candidate(tag, knob_vals)
    ret = proc.wait()
    if ret != 0:
        print(f"[stage2] WARNING: run {tag} exited with code {ret}", flush=True)
        return (math.inf, math.inf, math.inf), {"reason": f"exit_{ret}", "tag": tag, "knobs": knob_vals}

    rdir = make_run_dir(tag)
    score, info = score_run(rdir)
    plot_profile(rdir)
    info["run_dir"] = str(rdir)
    info["tag"] = tag
    info["knobs"] = {k: float(v) for k, v in knob_vals.items()}
    print(f"[stage2] result {tag}: score={fmt_score(score)} info={info}", flush=True)
    return score, info


def coordinate_neighbors(
    keys: List[str],
    value_lists: Dict[str, List[float]],
    idx_map: Dict[str, int],
) -> List[Dict[str, int]]:
    """
    Generate neighbors by varying one knob at a time to its other discrete values.
    Returns a list of new idx_map dicts.
    """
    nbrs: List[Dict[str, int]] = []
    handled = set()

    # For chamfer pairs, vary depth+angle together so we don't generate invalid no-op cases.
    for dkey, akey in CHAMFER_PAIRS:
        if dkey in idx_map and akey in idx_map:
            handled.add(dkey)
            handled.add(akey)
            dvals = value_lists[dkey]
            avals = value_lists[akey]
            for di in range(len(dvals)):
                for ai in range(len(avals)):
                    if di == idx_map[dkey] and ai == idx_map[akey]:
                        continue
                    nm = dict(idx_map)
                    nm[dkey] = di
                    nm[akey] = ai
                    nbrs.append(nm)

    # Other knobs vary one at a time.
    for k in keys:
        if k in handled:
            continue
        vals = value_lists[k]
        cur_idx = idx_map[k]
        for j in range(len(vals)):
            if j == cur_idx:
                continue
            nm = dict(idx_map)
            nm[k] = j
            nbrs.append(nm)
    random.shuffle(nbrs)
    return nbrs


def main():
    refine_lists = parse_refine_dump()
    if not refine_lists:
        print("[stage2] REFINE_DUMP empty or invalid; running single seed case.")
        base_tag = "s2_seed"
        score, info = evaluate_candidate(base_tag, {})
        outdir = RESULTS / "_stage2_best"
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "stage2_best.json").write_text(json.dumps({
            "best_score": score_to_list(score),
            "best": info,
            "keys": [],
            "values": {},
            "history": [{"score": score_to_list(score), **info}],
        }, indent=2))
        return

    keys = sorted(refine_lists.keys())
    value_lists: Dict[str, List[float]] = refine_lists

    # Discrete index per knob; start from the middle value (center)
    idx_map: Dict[str, int] = {}
    for k in keys:
        n = len(value_lists[k])
        if n == 0:
            continue
        idx_map[k] = n // 2  # for triplets, this is the center

    # Repair seed indices if they imply a no-op chamfer (depth=0 or angle=0 only).
    for dkey, akey in CHAMFER_PAIRS:
        if dkey not in idx_map or akey not in idx_map:
            continue
        dvals = value_lists[dkey]
        avals = value_lists[akey]
        d = float(dvals[idx_map[dkey]])
        a = float(avals[idx_map[akey]])
        if is_zero(d) ^ is_zero(a):
            if is_zero(d):
                # Prefer no-chamfer if angle has a zero entry; otherwise force depth non-zero.
                ai = closest_index(avals, want_zero=True)
                if ai >= 0:
                    idx_map[akey] = ai
                else:
                    di = closest_index(dvals, want_zero=False)
                    if di >= 0:
                        idx_map[dkey] = di
            else:
                # Angle is zero, depth non-zero.
                di = closest_index(dvals, want_zero=True)
                if di >= 0:
                    idx_map[dkey] = di
                else:
                    ai = closest_index(avals, want_zero=False)
                    if ai >= 0:
                        idx_map[akey] = ai

    def materialize(map_idx: Dict[str, int]) -> Dict[str, float]:
        return {k: value_lists[k][map_idx[k]] for k in keys if k in map_idx}

    # Evaluate the seed once (center point)
    best_idx = dict(idx_map)
    best_vals = materialize(best_idx)
    best_tag = "s2_seed"
    if not is_valid_chamfer(best_vals):
        best_score = (math.inf, math.inf, math.inf)
        best_info = {"reason": "invalid_chamfer_seed", "knobs": best_vals}
    else:
        best_score, best_info = evaluate_candidate(best_tag, best_vals)
    print(f"[stage2] initial best: score={fmt_score(best_score)} tag={best_tag}", flush=True)

    history = []
    history.append({"tag": best_tag, "score": best_score, **best_info})

    no_improve = 0

    for it in range(MAX_ITERS):
        print(f"[stage2] iteration {it}", flush=True)

        neighbors = coordinate_neighbors(keys, value_lists, best_idx)
        neighbors = [nm for nm in neighbors if is_valid_chamfer(materialize(nm))]
        if not neighbors:
            print("[stage2] no neighbors to explore; stopping.")
            break

        improved_this_iter = False
        cand_idx = 0

        while cand_idx < len(neighbors):
            batch = neighbors[cand_idx:cand_idx + BATCH_SIZE]
            cand_idx += len(batch)

            procs = []
            tags = []
            maps = []
            knob_snapshots = []

            # Launch batch
            for j, nm in enumerate(batch):
                tag = f"s2_it{it}_c{cand_idx - len(batch) + j}"
                vals = materialize(nm)
                p = launch_candidate(tag, vals)
                procs.append(p)
                tags.append(tag)
                maps.append(nm)
                knob_snapshots.append(vals)

            # Wait and score batch
            for p, tag, nm, vals in zip(procs, tags, maps, knob_snapshots):
                ret = p.wait()
                if ret != 0:
                    print(f"[stage2] WARNING: run {tag} exited with code {ret}", flush=True)
                    score = (math.inf, math.inf, math.inf)
                    info = {"reason": f"exit_{ret}", "tag": tag, "knobs": vals}
                else:
                    rdir = make_run_dir(tag)
                    score, info = score_run(rdir)
                    plot_profile(rdir)
                    info["run_dir"] = str(rdir)
                    info["tag"] = tag
                    info["knobs"] = vals

                history.append({"tag": tag, "score": score, **info})
                print(f"[stage2] candidate {tag}: score={fmt_score(score)}", flush=True)

                if score < best_score:
                    print(f"[stage2] NEW BEST: {fmt_score(score)} (prev {fmt_score(best_score)}) at {tag}", flush=True)
                    best_score = score
                    best_info = info
                    best_idx = dict(nm)
                    improved_this_iter = True

        if improved_this_iter:
            no_improve = 0
        else:
            no_improve += 1
            print(f"[stage2] no improvement this iteration (streak={no_improve})", flush=True)
            if no_improve >= NO_IMPROVE_MAX:
                print("[stage2] stopping due to no improvement.", flush=True)
                break

    # Write summary
    outdir = RESULTS / "_stage2_best"
    outdir.mkdir(parents=True, exist_ok=True)
    hist_sorted = sorted(history, key=lambda h: score_key(h.get("score", (math.inf, math.inf, math.inf))))
    hist_out = []
    for h in hist_sorted[:50]:
        row = dict(h)
        row["score"] = score_to_list(h.get("score", (math.inf, math.inf, math.inf)))
        hist_out.append(row)

    summary = {
        "best_score": score_to_list(best_score),
        "best": best_info,
        "history": hist_out,
        "keys": keys,
        "values": {k: [float(x) for x in vs] for k, vs in value_lists.items()},
        "score_config": {
            "mode": SCORE_MODE,
            "fields": {"I": I_FIELD, "DIV": DIV_FIELD, "Y": Y_FIELD},
            "weights": {"I": W_I, "DIV": W_DIV, "Y": W_Y},
            "refs": {"I": I_REF, "DIV": DIV_REF, "Y": Y_REF},
            "y_abs": (os.environ.get("ST2_Y_ABS", "1").strip() != "0"),
        },
    }
    (outdir / "stage2_best.json").write_text(json.dumps(summary, indent=2))
    print(f"[stage2] done. Best score={fmt_score(best_score)}, summary at {outdir/'stage2_best.json'}", flush=True)


if __name__ == "__main__":
    main()
PY
chmod +x "$PIPE_ROOT/stage2_opt.py"
fi


# =====================================================================================
# Stage-2 refine (optional)
# =====================================================================================

if [[ "$RUN_STAGE2" == "1" ]]; then
cat > "$PIPE_ROOT/refine_stage2.slurm" <<'EOF_REF2'
#!/bin/bash
#SBATCH -J two_grid_refine
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10GB
#SBATCH -t 24:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# SBATCH --partition=long

set -euo pipefail
umask 002
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
[ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash || true
module load gcc/11.x
module load libibsimu/1.0.6
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Center on Stage-1 best:
# shellcheck source=/dev/null
source "__RESULTS_DIR__/_stage1_best/latest/seed_stage2.env"

# Use explicit knobs if provided, else read auto-selected list
if [[ -z "${OPT_KNOBS:-}" && -f "__RESULTS_DIR__/_stage1_best/latest/knobs_stage2.env" ]]; then
  # shellcheck source=/dev/null
  source "__RESULTS_DIR__/_stage1_best/latest/knobs_stage2.env"
fi
if [[ -z "${OPT_KNOBS:-}" ]]; then
  OPT_KNOBS="SCR_UP_DEPTH_M,SCR_UP_ANGLE_DEG,SCR_DN_DEPTH_M,SCR_DN_ANGLE_DEG,ACC_UP_DEPTH_M,ACC_UP_ANGLE_DEG,ACC_DN_DEPTH_M,ACC_DN_ANGLE_DEG"
fi
echo "[refine] optimizing knobs: ${OPT_KNOBS:-<none>}"

triplet() { python3 - "$1" "$2" <<'PY'
import sys
c=float(sys.argv[1]); d=float(sys.argv[2])
print(f"{c-d} {c} {c+d}")
PY
}

def_delta() { case "$1" in *ANGLE_DEG) echo 2.0 ;; *) echo 1e-4 ;; esac; }

declare -A LISTS
add_list() {
  local name="$1" center="$2" delta_var="DELTA_$1"
  local d="${!delta_var:-$(def_delta "$name")}"
  local arr=( $(triplet "$center" "$d") )
  LISTS["$name"]="${arr[*]}"
}

IFS=',' read -r -a KNOBS <<< "${OPT_KNOBS:-}"
for K in "${KNOBS[@]}"; do
  K=$(echo "$K" | xargs)
  case "$K" in
    SCR_UP_DEPTH_M|SCR_DN_DEPTH_M|ACC_UP_DEPTH_M|ACC_DN_DEPTH_M|SCR_UP_ANGLE_DEG|SCR_DN_ANGLE_DEG|ACC_UP_ANGLE_DEG|ACC_DN_ANGLE_DEG|ACCEL_OFF_Y_M)
      add_list "$K" "${!K}"
      ;;
    *) echo "[refine] WARNING: unknown knob '$K' ignored";;
  esac
done
if [[ ${#LISTS[@]} -eq 0 ]]; then
  echo "[refine] no chamfer knobs selected; falling back to default chamfer list."
  for K in SCR_UP_DEPTH_M SCR_UP_ANGLE_DEG SCR_DN_DEPTH_M SCR_DN_ANGLE_DEG ACC_UP_DEPTH_M ACC_UP_ANGLE_DEG ACC_DN_DEPTH_M ACC_DN_ANGLE_DEG; do
    add_list "$K" "${!K}"
  done
fi

# Base env for two_grid_2d
export H="__H__" X_LEFT_M="__X_LEFT_M__" X_RIGHT_M="__X_RIGHT_M__" YBOX_M="__YBOX_M__"
export ACCEL_OFF_Y_M="${ACCEL_OFF_Y_M:-__ACCEL_OFF_Y_M__}" SCREEN_OFF_Y_M="${SCREEN_OFF_Y_M:-__SCREEN_OFF_Y_M__}"
export X_RIGHT_PHYS_M="__X_RIGHT_PHYS_M__"
export VS_V="__VS_V__" VA_V="__VA_V__" SAMPLE_V="__SAMPLE_V__"
export TUBE_WALL_T_M="__TUBE_WALL_T_M__" SAMPLE_PLATE_T_M="__SAMPLE_PLATE_T_M__"
export SAMPLE_PROFILE_X_M="__SAMPLE_PROFILE_X_M__"
export RUN_SOLVE="__RUN_SOLVE__" WRITE_PNG="__WRITE_PNG__" PNG_NAME="__PNG_NAME__"
export CLOSE_PAD_M="__CLOSE_PAD_M__" CLOSE_BASE_PX="__CLOSE_BASE_PX__"
export RESULTS_DIR="__RESULTS_DIR__" RUN_PREFIX="two_grid_refine"
export RUN_STAMP="__RUN_STAMP__"
export PIPE_ROOT="__PIPE_ROOT__"
export ENABLE_IONS="__ENABLE_IONS__"
export PLASMA_NI_M3="__PLASMA_NI_M3__"
export PLASMA_TE_EV="__PLASMA_TE_EV__"
export TRACE_ENABLE="0" VERT_PAD_M="__VERT_PAD_M__" VERT_XLEN_M="__VERT_XLEN_M__"
export VERT_YPAD_M="__VERT_YPAD_M__" VERT_STEP_M="__VERT_STEP_M__" VERT_THRESH_DEG="__VERT_THRESH_DEG__"
export PLOT_PROFILE="__PLOT_PROFILE__"
export PLOT_PROFILE_Y="__PLOT_PROFILE_Y__"
export PLOT_PROFILE_X_UNITS="__PLOT_PROFILE_X_UNITS__"
export PLOT_PROFILE_EXT="__PLOT_PROFILE_EXT__"


# Pack LISTS into env for python
DUMP="{"
for key in "${!LISTS[@]}"; do DUMP="$DUMP\"$key\":\"${LISTS[$key]}\","; done
DUMP="${DUMP%,}}"
export REFINE_DUMP="$DUMP"

if [[ "${STAGE2_MODE:-smart}" == "grid" ]]; then
  echo "[refine] STAGE2_MODE=grid -> using legacy Cartesian grid search."
  python3 - <<'PY'
import os, itertools, subprocess, time
dump=os.environ.get("REFINE_DUMP","{}").strip("{}")
d={}
if dump and dump!="{}":
    for pair in dump.split(","):
        pair=pair.strip()
        if not pair: continue
        k,v=pair.split(":",1)
        k=k.strip('" ')
        vals=v.strip('" ').split()
        d[k]=vals
keys=sorted(d.keys())
if not keys:
    print("[refine] No knobs selected; running single seed case.")
    combos=[{}]
else:
    combos=[dict(zip(keys,vals)) for vals in itertools.product(*[d[k] for k in keys])]
print(f"[refine] combos: {len(combos)} over knobs {keys}")
def is_zero(v):
    try:
        return abs(float(v)) <= 1.0e-12
    except Exception:
        return False

def valid_chamfer(combo):
    pairs = [
        ("SCR_UP_DEPTH_M","SCR_UP_ANGLE_DEG"),
        ("SCR_DN_DEPTH_M","SCR_DN_ANGLE_DEG"),
        ("ACC_UP_DEPTH_M","ACC_UP_ANGLE_DEG"),
        ("ACC_DN_DEPTH_M","ACC_DN_ANGLE_DEG"),
    ]
    for dkey, akey in pairs:
        if dkey in combo and akey in combo:
            d = combo[dkey]; a = combo[akey]
            if is_zero(d) ^ is_zero(a):
                return False
    return True

combos = [c for c in combos if valid_chamfer(c)]
print(f"[refine] combos after chamfer filter: {len(combos)}")

def predict_outdir(tag):
    base = os.environ.get("RESULTS_DIR", "./results")
    prefix = os.environ.get("RUN_PREFIX", "run")
    stamp = os.environ.get("RUN_STAMP", "").strip()
    if not stamp:
        stamp = time.strftime("%Y%m%d_%H%M%S")
    name = f"{prefix}_{stamp}"
    if tag:
        name += f"_{tag}"
    jobid = os.environ.get("SLURM_JOB_ID", "")
    if jobid:
        name += f"_j{jobid}"
    return os.path.join(base, name)

def plot_profile(tag):
    if os.environ.get("PLOT_PROFILE", "1").strip() == "0":
        return
    pipe_root = os.environ.get("PIPE_ROOT", ".")
    script = os.path.join(pipe_root, "plot_runlog_compact.py")
    if not os.path.exists(script):
        return
    outdir = predict_outdir(tag)
    prof = os.path.join(outdir, "sample_diameter_profile.json")
    if not os.path.exists(prof):
        prof = os.path.join(outdir, "sample_radial_profile.json")
    if not os.path.exists(prof):
        return
    ext = os.environ.get("PLOT_PROFILE_EXT", "png")
    out_html = os.path.join(outdir, f"sample_profile.{ext}")
    y = os.environ.get("PLOT_PROFILE_Y", "I_Apm")
    xunits = os.environ.get("PLOT_PROFILE_X_UNITS", "mm")
    cmd = ["python3", script, "sample-profile", "--json", prof, "--out", out_html,
           "--y", y, "--x-units", xunits]
    try:
        subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass

for i,combo in enumerate(combos):
    env=os.environ.copy()
    tag="ref_i{}".format(i)
    for k,val in combo.items():
        env[k]=val; tag+=f"_{k}={val}"
    env["RUN_TAG"]=tag
    print("[refine] running", tag)
    subprocess.check_call(["srun","--ntasks=1","__PIPE_ROOT__/bin/two_grid_2d"], env=env)
    plot_profile(tag)
PY
else
  echo "[refine] STAGE2_MODE=${STAGE2_MODE:-smart} -> using stage2_opt.py controller."
  python3 "__PIPE_ROOT__/stage2_opt.py"
fi
EOF_REF2


# Fill constants (robust escaping)
fill_ref2() { sed -i -e "s#__$1__#$(sed_escape "$2")#" "$PIPE_ROOT/refine_stage2.slurm"; }
fill_ref2 RESULTS_DIR "$RESULTS_DIR"
fill_ref2 H "$H";               fill_ref2 X_LEFT_M "$X_LEFT_M";  fill_ref2 X_RIGHT_M "$X_RIGHT_M"; fill_ref2 YBOX_M "$YBOX_M"
fill_ref2 ACCEL_OFF_Y_M "$ACCEL_OFF_Y_M"; fill_ref2 SCREEN_OFF_Y_M "$SCREEN_OFF_Y_M"
fill_ref2 X_RIGHT_PHYS_M "$X_RIGHT_PHYS_M"
fill_ref2 VS_V "$VS_V";         fill_ref2 VA_V "$VA_V";          fill_ref2 SAMPLE_V "$SAMPLE_V"
fill_ref2 TUBE_WALL_T_M "$TUBE_WALL_T_M"; fill_ref2 SAMPLE_PLATE_T_M "$SAMPLE_PLATE_T_M"
fill_ref2 SAMPLE_PROFILE_X_M "${SAMPLE_PROFILE_X_M:-nan}"
fill_ref2 RUN_SOLVE "$RUN_SOLVE"; fill_ref2 WRITE_PNG "$WRITE_PNG"; fill_ref2 PNG_NAME "$PNG_NAME"
fill_ref2 CLOSE_PAD_M "$CLOSE_PAD_M"; fill_ref2 CLOSE_BASE_PX "$CLOSE_BASE_PX"
fill_ref2 VERT_PAD_M "$VERT_PAD_M";   fill_ref2 VERT_XLEN_M "$VERT_XLEN_M"
fill_ref2 VERT_YPAD_M "$VERT_YPAD_M"; fill_ref2 VERT_STEP_M "$VERT_STEP_M"; fill_ref2 VERT_THRESH_DEG "$VERT_THRESH_DEG"
fill_ref2 PIPE_ROOT "$PIPE_ROOT"
fill_ref2 ENABLE_IONS "$ENABLE_IONS_STAGE2"
fill_ref2 PLASMA_NI_M3 "$PLASMA_NI_M3"
fill_ref2 PLASMA_TE_EV "$PLASMA_TE_EV"
fill_ref2 RUN_STAMP "$RUN_STAMP"
fill_ref2 PLOT_PROFILE "$PLOT_PROFILE"
fill_ref2 PLOT_PROFILE_Y "$PLOT_PROFILE_Y"
fill_ref2 PLOT_PROFILE_X_UNITS "$PLOT_PROFILE_X_UNITS"
fill_ref2 PLOT_PROFILE_EXT "$PLOT_PROFILE_EXT"
fi


# --- Build once and stage the binary ---
module load gcc/11.x
module load libibsimu/1.0.6

echo "[orchestrator] building two_grid_2d once..."
make -j

mkdir -p "$PIPE_ROOT/bin"
cp -f two_grid_2d "$PIPE_ROOT/bin/two_grid_2d"
chmod +x "$PIPE_ROOT/bin/two_grid_2d"

# =====================================================================================
# SUBMIT: Stage-1 array (chunked) -> reducer -> (optional) Stage-2
# =====================================================================================

# Max tasks per array submission (tune if your site allows bigger arrays)
MAX_ARRAY_CHUNK="${MAX_ARRAY_CHUNK:-1000}"

submit_ids=()
base=0
while (( base < TOTAL )); do
  remain=$(( TOTAL - base ))
  this=$(( remain < MAX_ARRAY_CHUNK ? remain : MAX_ARRAY_CHUNK ))
  # array indices for this chunk are 0..(this-1)
  ARRAY_SPEC="--array=0-$((this-1))%$ARRAY_CONC"
  echo "[orchestrator] sbatch $ARRAY_SPEC (ARRAY_BASE=$base) $PIPE_ROOT/sweep_two_grid.slurm"

  # Pass ARRAY_BASE to the batch via --export
  jid=$(sbatch --parsable $ARRAY_SPEC --export=ALL,ARRAY_BASE="$base" "$PIPE_ROOT/sweep_two_grid.slurm")
  echo "[orchestrator]  chunk @base=$base size=$this  -> jobid: $jid"

  submit_ids+=("$jid")
  base=$(( base + this ))
done

# Build dependency on ALL sweep chunks
dep="afterok"
for j in "${submit_ids[@]}"; do dep="$dep:$j"; done

REDUCE_JOBID=$(sbatch --parsable --dependency="$dep" "$PIPE_ROOT/reduce_best.slurm")
echo "[orchestrator] Reducer jobid      : $REDUCE_JOBID"

if [[ "$RUN_STAGE2" == "1" ]]; then
  REFINE_JOBID=$(sbatch --parsable --dependency=afterok:$REDUCE_JOBID "$PIPE_ROOT/refine_stage2.slurm")
  echo "[orchestrator] Stage-2 refine job : $REFINE_JOBID"
fi

echo "[orchestrator] Submitted. Monitor:"
echo "  squeue -u $USER -o \"%.10i %.5t %12j %8A %8a %20S\""
