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
RUN_PREFIX="two_grid_SC_Ramp_Rework_Test"

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
X_RIGHT_M="${X_RIGHT_M:-0.06}"
# Physical tube/sample length used for envelope projection (can exceed X_RIGHT_M)
X_RIGHT_PHYS_M="${X_RIGHT_PHYS_M:-0.55}"
YBOX_M=0.086
ACCEL_OFF_Y_M=0.0

# Potentials & downstream slices
VS_V=0
VA_V=-10000
SAMPLE_V=$VA_V
TUBE_WALL_T_M=0.0002
SAMPLE_PLATE_T_M=0.0002

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
ENABLE_IONS_STAGE2=${ENABLE_IONS_STAGE2:-0}   # set to 1 when you want Stage-2 with ions

# Plasma defaults (match C++ defaults; you can override here or via sbatch --export)
PLASMA_NI_M3=${PLASMA_NI_M3:-2e18}
PLASMA_TE_EV=${PLASMA_TE_EV:-3.7}

# Beam / SC tuning knobs
ION_J_SCALE=${ION_J_SCALE:-1.0} #1.0 = nominal Bohm current
SC_FACTOR=${SC_FACTOR:-0.1} # 1.0 = full space charge, <1.0 = partially neutralized

# =====================================================================================
# SWEEP LISTS
# =====================================================================================
# Geometry
AP_RAD_LIST=(0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004)                # 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005)
GRID_T_LIST=(0.007)                #0.004 0.005 0.006 0.007 0.008 0.009 0.01)
GAP_LIST=(0.001 0.0012 0.0013 0.0014 0.0016 0.0018 0.002)

# Chamfers (8 independent controls)
SCR_UP_DEPTH_LIST=(0.0)
SCR_UP_ANGLE_LIST=(0.0)
SCR_DN_DEPTH_LIST=(0.0) #0.002 0.003 0.004)
SCR_DN_ANGLE_LIST=(0) #20 45)
ACC_UP_DEPTH_LIST=(0.0) #0.002 0.003 0.004)
ACC_UP_ANGLE_LIST=(0) #20 45)
ACC_DN_DEPTH_LIST=(0.0)
ACC_DN_ANGLE_LIST=(0.0)

# =====================================================================================
# Stage-2 control
# =====================================================================================
RUN_STAGE2="${RUN_STAGE2:-0}"   # 0=no refine; 1=submit refine after reducer

# FIX: set an explicit knob list. If you want auto-pick later, set OPT_KNOBS="".
# (current default: aperture radius + grid thickness + gap)
OPT_KNOBS="${OPT_KNOBS:-AP_RAD_M, GRID_T_M, GAP_M}"

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
  "ACCEL_OFF_Y_M":    5.0e-4,
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
ST2_MAX_ITERS="${ST2_MAX_ITERS:-20}"      # max outer iterations
ST2_NO_IMPROVE="${ST2_NO_IMPROVE:-3}"    # stop after N iters with no better score
ST2_BATCH_SIZE="${ST2_BATCH_SIZE:-5}"    # concurrent srun evaluations
ST2_SEED="${ST2_SEED:-42069}"                # RNG seed (0 = deterministic per campaign)

# Scoring weights (verticality vs collimation vs current)
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
       ST2_W_VERT ST2_W_COLL ST2_W_CURR \
       ST2_I_TARGET_A ST2_I_MIN_FRAC ST2_I_MIN_TOTAL_A \
       ST2_THETA_REF_DEG ST2_R_REF_M ST2_I_REF_A ST2_I_MIN_A

# =====================================================================================
# INTERNALS (no edits)
# =====================================================================================

PIPE_ROOT="$PWD/.pipeline_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$PIPE_ROOT"
echo "[orchestrator] workspace: $PIPE_ROOT"

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

export H=__H__
export X_LEFT_M=__X_LEFT_M__
export X_RIGHT_M=__X_RIGHT_M__
export X_RIGHT_PHYS_M=__X_RIGHT_PHYS_M__
export YBOX_M=__YBOX_M__
export ACCEL_OFF_Y_M=__ACCEL_OFF_Y_M__

export VS_V=__VS_V__
export VA_V=__VA_V__
export SAMPLE_V=__SAMPLE_V__
export TUBE_WALL_T_M=__TUBE_WALL_T_M__
export SAMPLE_PLATE_T_M=__SAMPLE_PLATE_T_M__

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

# Pre-escape scalars
E_H="$(sed_escape "$H")"
E_XL="$(sed_escape "$X_LEFT_M")"
E_XR="$(sed_escape "$X_RIGHT_M")"
E_XRP="$(sed_escape "$X_RIGHT_PHYS_M")"
E_YB="$(sed_escape "$YBOX_M")"
E_OFF="$(sed_escape "$ACCEL_OFF_Y_M")"
E_VS="$(sed_escape "$VS_V")"
E_VA="$(sed_escape "$VA_V")"
E_SAMP="$(sed_escape "$SAMPLE_V")"
E_TW="$(sed_escape "$TUBE_WALL_T_M")"
E_EP="$(sed_escape "$SAMPLE_PLATE_T_M")"
E_RUNSOLVE="$(sed_escape "$RUN_SOLVE")"
E_WRITEPNG="$(sed_escape "$WRITE_PNG")"
E_PNGNAME="$(sed_escape "$PNG_NAME")"
E_CP="$(sed_escape "$CLOSE_PAD_M")"
E_CB="$(sed_escape "$CLOSE_BASE_PX")"
E_RDIR="$(sed_escape "$RESULTS_DIR")"
E_RPFX="$(sed_escape "$RUN_PREFIX")"
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
  -e "s#__H__#${E_H}#" \
  -e "s#__X_LEFT_M__#${E_XL}#" \
  -e "s#__X_RIGHT_M__#${E_XR}#" \
  -e "s#__X_RIGHT_PHYS_M__#${E_XRP}#" \
  -e "s#__YBOX_M__#${E_YB}#" \
  -e "s#__ACCEL_OFF_Y_M__#${E_OFF}#" \
  -e "s#__VS_V__#${E_VS}#" \
  -e "s#__VA_V__#${E_VA}#" \
  -e "s#__SAMPLE_V__#${E_SAMP}#" \
  -e "s#__TUBE_WALL_T_M__#${E_TW}#" \
  -e "s#__SAMPLE_PLATE_T_M__#${E_EP}#" \
  -e "s#__RUN_SOLVE__#${E_RUNSOLVE}#" \
  -e "s#__WRITE_PNG__#${E_WRITEPNG}#" \
  -e "s#__PNG_NAME__#${E_PNGNAME}#" \
  -e "s#__CLOSE_PAD_M__#${E_CP}#" \
  -e "s#__CLOSE_BASE_PX__#${E_CB}#" \
  -e "s#__PLOT_FONT_PX__#${E_PFPX}#" \
  -e "s#__RESULTS_DIR__#${E_RDIR}#" \
  -e "s#__RUN_PREFIX__#${E_RPFX}#" \
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
TOTAL=$(( nR*nT*nG*nSUd*nSUa*nSDd*nSDa*nAUd*nAUa*nADd*nADa ))
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
iADa=$(( i % ${#ACC_DN_ANGLE_LIST[@]} ))

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

sanitize(){ echo "$1" | sed -e 's/\./p/g' -e 's/-/m/g'; }
rt=$(sanitize "$AP_RAD_M"); gt=$(sanitize "$GRID_T_M"); gp=$(sanitize "$GAP_M")
suD=$(sanitize "$SCR_UP_DEPTH_M"); suA=$(sanitize "$SCR_UP_ANGLE_DEG")
sdD=$(sanitize "$SCR_DN_DEPTH_M"); sdA=$(sanitize "$SCR_DN_ANGLE_DEG")
auD=$(sanitize "$ACC_UP_DEPTH_M");  auA=$(sanitize "$ACC_UP_ANGLE_DEG")
adD=$(sanitize "$ACC_DN_DEPTH_M");  adA=$(sanitize "$ACC_DN_ANGLE_DEG")

export RUN_TAG="i${idx}_rad${rt}_t${gt}_gap${gp}_SU(${suD},${suA})_SD(${sdD},${sdA})_AU(${auD},${auA})_AD(${adD},${adA})"

echo "[$(date)] running..."
srun --ntasks=1 "__PIPE_ROOT__/bin/two_grid_2d"
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

def score(m):
    """
    Stage-1 score:

    Prefer the column-based metric if available:
      metrics["column_verticality"][COL_METRIC]
    else fall back to the legacy area metric:
      metrics["verticality"][METRIC]

    Lower score is better.
    """
    m = m or {}
    v_col = m.get("column_verticality") or {}
    v     = m.get("verticality") or {}

    # Column metric first (uses averaged field vectors along columns)
    if v_col:
        try:
            val = float(v_col.get(COL_METRIC, 1e99))
        except Exception:
            val = 1e99
        return val, {"column_verticality": v_col, "verticality": v}

    # Fallback: old near-field area metric
    try:
        val = float(v.get(METRIC, 1e99))
    except Exception:
        val = 1e99
    return val, v


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
        "AP_RAD_M","GRID_T_M","GAP_M","ACCEL_OFF_Y_M","H","X_LEFT_M","X_RIGHT_M","YBOX_M",
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
        met = loadj(run/"metrics.json") or {}

        v_col = met.get("column_verticality") or {}
        v     = met.get("verticality") or {}

        if v_col:
            y = v_col.get(COL_METRIC, None)
        else:
            y = v.get(METRIC, None)

        if y is None:
            continue

        m   = loadj(meta) or {}
        scr = (m.get("screen_chamfer") or {})
        acc = (m.get("accel_chamfer") or {})
        feat = {
          "AP_RAD_M": m.get("AP_RAD_M"), "GRID_T_M": m.get("GRID_T_M"),
          "GAP_M": m.get("GAP_M"), "ACCEL_OFF_Y_M": m.get("ACCEL_OFF_Y_M"),
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
        met = loadj(run/"metrics.json")
        if not met: continue
        s, vfull = score(met)
        if s>=1e98: continue
        items.append({"run":run, "meta":loadj(meta), "metrics":vfull, "score":s})
    if not items:
        print("[reduce] no runs with metrics found", file=sys.stderr); sys.exit(2)

    items.sort(key=lambda r: r["score"], reverse=not ASCEND)
    stamp=time.strftime("%Y%m%d_%H%M%S")
    outdir=RESULTS/"_stage1_best"/stamp
    outdir.mkdir(parents=True, exist_ok=True)

    top=[]
    for r in items[:min(TOPK,len(items))]:
        row={"run_dir":str(r["run"]), "score_metric":METRIC, "score_value":r["score"], "verticality":r["metrics"]}
        row.update(extract_geom(r["meta"]))
        top.append(row)
    (outdir/"topk.json").write_text(json.dumps({"metric":METRIC,"ascending":ASCEND,"topk":top}, indent=2))

    best=items[0]
    geom=extract_geom(best["meta"])
    (outdir/"best.json").write_text(json.dumps({
        "metric":METRIC,"ascending":ASCEND,"score":best["score"],
        "run_dir":str(best["run"]), "verticality":best["metrics"], "geom":geom
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

    print(f"[reduce] best: {best['run']}  score={best['score']:.6e}")
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
from typing import Dict, List, Tuple, Any

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

# Verticality metric name (should match what two_grid_2d wrote into metrics.json)
VERT_METRIC = os.environ.get("VERT_METRIC", "mean_angle_deg")

# Stage-2 scoring parameters (can override via env)
W_VERT = float(os.environ.get("ST2_W_VERT", "1.0"))
W_COLL = float(os.environ.get("ST2_W_COLL", "1.0"))
W_CURR = float(os.environ.get("ST2_W_CURR", "0.3"))

THETA_REF = float(os.environ.get("ST2_THETA_REF_DEG", "2.0"))      # deg
R_REF     = float(os.environ.get("ST2_R_REF_M", "2.0e-3"))         # m
I_TARGET  = float(os.environ.get("ST2_I_TARGET_A", "4.0"))         # A, full-grid target
I_MIN_FRAC = float(os.environ.get("ST2_I_MIN_FRAC", "0.66"))

# If ST2_I_MIN_TOTAL_A is a non-empty, numeric string, use it.
# Otherwise derive from fraction * target.
_raw_min = os.environ.get("ST2_I_MIN_TOTAL_A", "").strip()
if _raw_min:
    try:
        I_MIN_TOTAL = float(_raw_min)
    except ValueError:
        # fall back if user set a non-numeric value
        I_MIN_TOTAL = I_TARGET * I_MIN_FRAC
else:
    I_MIN_TOTAL = I_TARGET * I_MIN_FRAC


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


def score_run(run_dir: Path) -> Tuple[float, Dict[str, Any]]:
    """
    Compute a scalar score from metrics.json + beam_metrics.json,
    using aperture-scaled total beam current.

    Lower score is better.
    """
    met = load_json(run_dir / "metrics.json") or {}
    v   = (met.get("verticality") or {})
    theta = v.get(VERT_METRIC, None)
    if theta is None:
        return math.inf, {"reason": "no_verticality"}

    try:
        theta = float(theta)
    except Exception:
        return math.inf, {"reason": "bad_verticality"}

    # Get aperture radius for this run from meta.json (preferred) or env fallback
    meta = load_json(run_dir / "meta.json") or {}
    ap_rad_m = None
    if isinstance(meta, dict):
        ap_rad_m = meta.get("AP_RAD_M")
    try:
        if ap_rad_m is None:
            ap_rad_m = float(os.environ.get("AP_RAD_M", "1.0e-3"))
        else:
            ap_rad_m = float(ap_rad_m)
    except Exception:
        ap_rad_m = 1.0e-3  # conservative fallback: 1 mm

    bm = load_json(run_dir / "beam_metrics.json") or {}
    coll = (bm.get("collimation") or {})
    sm   = (bm.get("sample") or {})

    has_sample = bool(coll.get("has_sample_beam", False))
    good_single = bool(coll.get("good_single_beam", False))

    # Hard filters: must have a good single beam hitting the sample
    if not has_sample:
        return math.inf, {"theta_deg": theta, "reason": "no_sample_beam"}
    if not good_single:
        return math.inf, {"theta_deg": theta, "reason": "not_single_beam"}

    try:
        I_A   = float(sm.get("I_A", 0.0))          # single-beamlet 3D current [A]
        r_rms = float(sm.get("y_rms_m", 0.0))      # rms radius at sample [m]
    except Exception:
        return math.inf, {"theta_deg": theta, "reason": "bad_beam_numbers"}

    # Estimate how many such apertures could fit on the grid
    n_beamlets = estimate_beamlet_count(ap_rad_m)
    I_total = I_A * n_beamlets

    # Hard filter: require minimum total current capability
    if I_total < I_MIN_TOTAL:
        return math.inf, {
            "theta_deg": theta,
            "I_A": I_A,
            "I_total_A": I_total,
            "n_beamlets": n_beamlets,
            "reason": "total_current_below_min",
        }

    # Normalized components for cost function
    theta_norm = theta / (THETA_REF if THETA_REF > 0 else 1.0)
    r_norm     = r_rms / (R_REF     if R_REF     > 0 else 1.0)

    if I_TARGET > 0.0:
        I_mismatch = abs(I_total - I_TARGET) / I_TARGET
    else:
        I_mismatch = 0.0

    # Composite score:
    #   lower theta (straighter equipotentials)  -> better
    #   smaller beam radius at sample           -> better
    #   total current closer to 4 A             -> better (smaller mismatch)
    score = W_VERT * theta_norm + W_COLL * r_norm + W_CURR * I_mismatch

    info = {
        "theta_deg": theta,
        "theta_norm": theta_norm,
        "r_rms_m": r_rms,
        "r_norm": r_norm,
        "I_A": I_A,
        "I_total_A": I_total,
        "I_mismatch": I_mismatch,
        "n_beamlets": n_beamlets,
        "ap_rad_m": ap_rad_m,
        "has_sample_beam": has_sample,
        "good_single_beam": good_single,
    }
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


def evaluate_candidate(tag: str, knob_vals: Dict[str, float]) -> Tuple[float, Dict[str, Any]]:
    """
    Launch and wait for a candidate, then read its metrics and score it.
    """
    proc = launch_candidate(tag, knob_vals)
    ret = proc.wait()
    if ret != 0:
        print(f"[stage2] WARNING: run {tag} exited with code {ret}", flush=True)
        return math.inf, {"reason": f"exit_{ret}", "tag": tag, "knobs": knob_vals}

    rdir = make_run_dir(tag)
    score, info = score_run(rdir)
    info["run_dir"] = str(rdir)
    info["tag"] = tag
    info["knobs"] = {k: float(v) for k, v in knob_vals.items()}
    print(f"[stage2] result {tag}: score={score:.6g} info={info}", flush=True)
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
    for k in keys:
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
            "best_score": score,
            "best": info,
            "keys": [],
            "values": {},
            "history": [info],
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

    def materialize(map_idx: Dict[str, int]) -> Dict[str, float]:
        return {k: value_lists[k][map_idx[k]] for k in keys if k in map_idx}

    # Evaluate the seed once (center point)
    best_idx = dict(idx_map)
    best_vals = materialize(best_idx)
    best_tag = "s2_seed"
    best_score, best_info = evaluate_candidate(best_tag, best_vals)
    print(f"[stage2] initial best: score={best_score:.6g} tag={best_tag}", flush=True)

    history = []
    history.append({"tag": best_tag, "score": best_score, **best_info})

    no_improve = 0

    for it in range(MAX_ITERS):
        print(f"[stage2] iteration {it}", flush=True)

        neighbors = coordinate_neighbors(keys, value_lists, best_idx)
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
                    score = math.inf
                    info = {"reason": f"exit_{ret}", "tag": tag, "knobs": vals}
                else:
                    rdir = make_run_dir(tag)
                    score, info = score_run(rdir)
                    info["run_dir"] = str(rdir)
                    info["tag"] = tag
                    info["knobs"] = vals

                history.append({"tag": tag, "score": score, **info})
                print(f"[stage2] candidate {tag}: score={score:.6g}", flush=True)

                if score < best_score:
                    print(f"[stage2] NEW BEST: {score:.6g} (prev {best_score:.6g}) at {tag}", flush=True)
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
    summary = {
        "best_score": best_score,
        "best": best_info,
        "history": sorted(history, key=lambda h: h.get("score", math.inf))[:50],
        "keys": keys,
        "values": {k: [float(x) for x in vs] for k, vs in value_lists.items()},
    }
    (outdir / "stage2_best.json").write_text(json.dumps(summary, indent=2))
    print(f"[stage2] done. Best score={best_score:.6g}, summary at {outdir/'stage2_best.json'}", flush=True)


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
    AP_RAD_M|GAP_M|GRID_T_M|ACCEL_OFF_Y_M|SCR_UP_DEPTH_M|SCR_DN_DEPTH_M|ACC_UP_DEPTH_M|ACC_DN_DEPTH_M|SCR_UP_ANGLE_DEG|SCR_DN_ANGLE_DEG|ACC_UP_ANGLE_DEG|ACC_DN_ANGLE_DEG)
      add_list "$K" "${!K}"
      ;;
    *) echo "[refine] WARNING: unknown knob '$K' ignored";;
  esac
done

# Base env for two_grid_2d
export H="__H__" X_LEFT_M="__X_LEFT_M__" X_RIGHT_M="__X_RIGHT_M__" YBOX_M="__YBOX_M__"
export X_RIGHT_PHYS_M="__X_RIGHT_PHYS_M__"
export VS_V="__VS_V__" VA_V="__VA_V__" SAMPLE_V="__SAMPLE_V__"
export TUBE_WALL_T_M="__TUBE_WALL_T_M__" SAMPLE_PLATE_T_M="__SAMPLE_PLATE_T_M__"
export RUN_SOLVE="__RUN_SOLVE__" WRITE_PNG="__WRITE_PNG__" PNG_NAME="__PNG_NAME__"
export CLOSE_PAD_M="__CLOSE_PAD_M__" CLOSE_BASE_PX="__CLOSE_BASE_PX__"
export RESULTS_DIR="__RESULTS_DIR__" RUN_PREFIX="two_grid_refine"
export PIPE_ROOT="__PIPE_ROOT__"
export ENABLE_IONS="__ENABLE_IONS__"
export PLASMA_NI_M3="__PLASMA_NI_M3__"
export PLASMA_TE_EV="__PLASMA_TE_EV__"
export TRACE_ENABLE="0" VERT_PAD_M="__VERT_PAD_M__" VERT_XLEN_M="__VERT_XLEN_M__"
export VERT_YPAD_M="__VERT_YPAD_M__" VERT_STEP_M="__VERT_STEP_M__" VERT_THRESH_DEG="__VERT_THRESH_DEG__"


# Pack LISTS into env for python
DUMP="{"
for key in "${!LISTS[@]}"; do DUMP="$DUMP\"$key\":\"${LISTS[$key]}\","; done
DUMP="${DUMP%,}}"
export REFINE_DUMP="$DUMP"

if [[ "${STAGE2_MODE:-smart}" == "grid" ]]; then
  echo "[refine] STAGE2_MODE=grid -> using legacy Cartesian grid search."
  python3 - <<'PY'
import os, itertools, subprocess
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
for i,combo in enumerate(combos):
    env=os.environ.copy()
    tag="ref_i{}".format(i)
    for k,val in combo.items():
        env[k]=val; tag+=f"_{k}={val}"
    env["RUN_TAG"]=tag
    print("[refine] running", tag)
    subprocess.check_call(["srun","--ntasks=1","__PIPE_ROOT__/bin/two_grid_2d"], env=env)
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
fill_ref2 X_RIGHT_PHYS_M "$X_RIGHT_PHYS_M"
fill_ref2 VS_V "$VS_V";         fill_ref2 VA_V "$VA_V";          fill_ref2 SAMPLE_V "$SAMPLE_V"
fill_ref2 TUBE_WALL_T_M "$TUBE_WALL_T_M"; fill_ref2 SAMPLE_PLATE_T_M "$SAMPLE_PLATE_T_M"
fill_ref2 RUN_SOLVE "$RUN_SOLVE"; fill_ref2 WRITE_PNG "$WRITE_PNG"; fill_ref2 PNG_NAME "$PNG_NAME"
fill_ref2 CLOSE_PAD_M "$CLOSE_PAD_M"; fill_ref2 CLOSE_BASE_PX "$CLOSE_BASE_PX"
fill_ref2 VERT_PAD_M "$VERT_PAD_M";   fill_ref2 VERT_XLEN_M "$VERT_XLEN_M"
fill_ref2 VERT_YPAD_M "$VERT_YPAD_M"; fill_ref2 VERT_STEP_M "$VERT_STEP_M"; fill_ref2 VERT_THRESH_DEG "$VERT_THRESH_DEG"
fill_ref2 PIPE_ROOT "$PIPE_ROOT"
fill_ref2 ENABLE_IONS "$ENABLE_IONS_STAGE2"
fill_ref2 PLASMA_NI_M3 "$PLASMA_NI_M3"
fill_ref2 PLASMA_TE_EV "$PLASMA_TE_EV"
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


