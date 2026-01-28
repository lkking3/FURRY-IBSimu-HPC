#!/usr/bin/env python3
"""
update_runlog_compact.py

Creates/updates a compact runlog CSV where:
- Each row = one run directory under RESULTS_DIR containing meta.json
- Columns are ordered: identity -> assigned variables -> beam outputs

Design goals:
- No system-estimate columns (no N_beamlets scaling, no "sys" perveance)
- Perveance kept simple: CL perveance, AG geometric perveance, and their ratio (normalized perveance)
- Keep only one timing metric: total wall time

Usage examples:
  python update_runlog_compact.py --results-dir results --scan --refresh
  python update_runlog_compact.py --results-dir results --csv results/runlog_compact.csv --scan
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

def _load_json(p: Path) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(p.read_text())
    except Exception:
        return None

def _get(d: Any, path: str, default=None):
    cur = d
    for key in path.split("."):
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur

def _as_float(x) -> Optional[float]:
    if x is None:
        return None
    try:
        v = float(x)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except Exception:
        return None

def _fmt(x):
    if x is None:
        return ""
    if isinstance(x, bool):
        return "1" if x else "0"
    return x

IDENTITY_COLS = [
    "RUN_DIR",
    "runname",
    "RUN_TAG",
    "RUN_STAMP",
]

ASSIGNED_COLS = [
    "VS_V","VA_V","SAMPLE_V",
    "H","AP_RAD_M","GRID_T_M","GAP_M","ACCEL_OFF_Y_M","SCREEN_OFF_Y_M",
    "X_LEFT_M","X_RIGHT_M","X_RIGHT_PHYS_M","TRUNCATED_DRIFT","YBOX_M",
    "TUBE_X_START_M","TUBE_X_END_M","TUBE_WALL_T_M","SAMPLE_PLATE_T_M",
    "SCR_UP_DEPTH_M","SCR_UP_ANGLE_DEG","SCR_DN_DEPTH_M","SCR_DN_ANGLE_DEG",
    "ACC_UP_DEPTH_M","ACC_UP_ANGLE_DEG","ACC_DN_DEPTH_M","ACC_DN_ANGLE_DEG",
    "ENABLE_IONS","ION_ITER_MAX","ION_M_AMU","ION_TP_EV","ION_TT_EV",
    "PLASMA_NI_M3","PLASMA_TE_EV","PLASMA_UP_V",
    "ION_J_SCALE","SC_FACTOR","SC_RAMP_START_M","SC_RAMP_LEN_M",
]

OUTPUT_COLS = [
    "WALL_TOTAL_S",
    "HAS_SAMPLE_BEAM",
    "GOOD_SINGLE_BEAM",
    "LOST_TO_SIDEWALLS",
    "I_PG_IN_A",
    "I_AG_OUT_A",
    "I_SAMPLE_A",
    "Y_RMS_MAX_M",
    "Y_ABSMAX_MAX_M",
    "TUBE_INNER_HALF_M",
    "Y_RMS_SAMPLE_M",
    "Y_ABSMAX_SAMPLE_M",
    "FOOTPRINT_PRED_M",
    "FOOTPRINT_PRED_400MM_M",
    "FOOTPRINT_PRED_500MM_M",
    "FOOTPRINT_PRED_600MM_M",
    "CLEARANCE_PRED_M",
    "HITS_TUBE_PRED",
    "DIVERGENCE_ANGLE_DEG",
    "P_CL_A_PER_V32",
    "P_GEOM_AG_A_PER_V32",
    "P_NORM_AG_OVER_CL",
]

ALL_COLS = IDENTITY_COLS + ASSIGNED_COLS + OUTPUT_COLS

def extract_row(run_dir: Path, results_dir: Path) -> Dict[str, Any]:
    meta = _load_json(run_dir / "meta.json") or {}
    bm   = _load_json(run_dir / "beam_metrics.json") or {}
    tj   = _load_json(run_dir / "timing.json") or {}

    row: Dict[str, Any] = {c: None for c in ALL_COLS}

    try:
        row["RUN_DIR"] = str(run_dir.relative_to(results_dir)).replace("\\", "/")
    except Exception:
        row["RUN_DIR"] = str(run_dir).replace("\\", "/")
    row["runname"]   = meta.get("runname")
    row["RUN_TAG"]   = meta.get("RUN_TAG")
    row["RUN_STAMP"] = meta.get("RUN_STAMP")

    for k in ["VS_V","VA_V","SAMPLE_V","H","AP_RAD_M","GRID_T_M","GAP_M","ACCEL_OFF_Y_M","SCREEN_OFF_Y_M",
              "X_LEFT_M","X_RIGHT_M","X_RIGHT_PHYS_M","TRUNCATED_DRIFT","YBOX_M"]:
        if k in meta:
            row[k] = meta.get(k)

    tube = meta.get("tube") or {}
    row["TUBE_X_START_M"] = tube.get("x_start")
    row["TUBE_X_END_M"]   = tube.get("x_end")
    row["TUBE_WALL_T_M"]  = tube.get("wall_t")
    row["SAMPLE_PLATE_T_M"]= tube.get("endplate_t")
    
    scr = meta.get("screen_chamfer") or {}
    acc = meta.get("accel_chamfer") or {}
    row["SCR_UP_DEPTH_M"]   = scr.get("up_depth")
    row["SCR_UP_ANGLE_DEG"] = scr.get("up_angle_deg")
    row["SCR_DN_DEPTH_M"]   = scr.get("dn_depth")
    row["SCR_DN_ANGLE_DEG"] = scr.get("dn_angle_deg")
    row["ACC_UP_DEPTH_M"]   = acc.get("up_depth")
    row["ACC_UP_ANGLE_DEG"] = acc.get("up_angle_deg")
    row["ACC_DN_DEPTH_M"]   = acc.get("dn_depth")
    row["ACC_DN_ANGLE_DEG"] = acc.get("dn_angle_deg")

    phys = meta.get("physics") or {}
    for k in ["ENABLE_IONS","ION_ITER_MAX","ION_M_AMU","ION_TP_EV","ION_TT_EV",
              "PLASMA_NI_M3","PLASMA_TE_EV","PLASMA_UP_V","ION_J_SCALE","SC_FACTOR",
              "SC_RAMP_START_M","SC_RAMP_LEN_M"]:
        if k in phys:
            row[k] = phys.get(k)

    row["WALL_TOTAL_S"] = tj.get("wall_total_s")

    coll = bm.get("collimation") or {}
    cur  = bm.get("currents") or {}
    sm   = bm.get("sample") or {}

    row["HAS_SAMPLE_BEAM"]   = coll.get("has_sample_beam")
    row["GOOD_SINGLE_BEAM"]  = coll.get("good_single_beam")
    row["LOST_TO_SIDEWALLS"] = coll.get("lost_to_sidewalls")

    row["I_PG_IN_A"]  = cur.get("I_pg_in_A")
    row["I_AG_OUT_A"] = cur.get("I_ag_out_A")
    row["I_SAMPLE_A"] = sm.get("I_A")

    row["Y_RMS_MAX_M"]       = coll.get("y_rms_max_m")
    row["Y_ABSMAX_MAX_M"]    = coll.get("y_absmax_max_m")
    row["TUBE_INNER_HALF_M"] = coll.get("tube_inner_half_m")
    row["Y_RMS_SAMPLE_M"]    = sm.get("y_rms_m")
    row["Y_ABSMAX_SAMPLE_M"] = sm.get("y_absmax_m")

    # optional predicted footprint fields (filled only if your run writes them)
    for out_key, paths in [
        ("FOOTPRINT_PRED_M", ["envelope.footprint_pred_m", "envelope.footprint_m", "FOOTPRINT_PRED_M"]),
        ("FOOTPRINT_PRED_400MM_M", ["envelope.footprint_pred_400mm_m", "FOOTPRINT_PRED_400MM_M"]),
        ("FOOTPRINT_PRED_500MM_M", ["envelope.footprint_pred_500mm_m", "FOOTPRINT_PRED_500MM_M"]),
        ("FOOTPRINT_PRED_600MM_M", ["envelope.footprint_pred_600mm_m", "FOOTPRINT_PRED_600MM_M"]),
        ("CLEARANCE_PRED_M", ["envelope.clearance_pred_m", "CLEARANCE_PRED_M"]),
        ("HITS_TUBE_PRED",   ["envelope.hits_tube_pred", "HITS_TUBE_PRED"]),
        ("DIVERGENCE_ANGLE_DEG", ["envelope.divergence_angle_deg", "DIVERGENCE_ANGLE_DEG"]),
    ]:
        val = None
        for p in paths:
            val = _get(meta, p, None)
            if val is None:
                val = _get(bm, p, None)
            if val is not None:
                break
        row[out_key] = val

    pv = bm.get("perveance") or {}
    beamlet = (pv.get("beamlet") or {})

    P_cl = _as_float(pv.get("P_CL_A_per_V32"))
    P_ag = _as_float(beamlet.get("P_geom_ag_A_per_V32"))

    row["P_CL_A_PER_V32"] = P_cl
    row["P_GEOM_AG_A_PER_V32"] = P_ag
    if P_cl and P_cl > 0 and P_ag is not None:
        row["P_NORM_AG_OVER_CL"] = P_ag / P_cl
    else:
        row["P_NORM_AG_OVER_CL"] = None

    return row

def read_existing(csv_path: Path) -> Dict[str, Dict[str, Any]]:
    if not csv_path.exists():
        return {}
    out: Dict[str, Dict[str, Any]] = {}
    with csv_path.open("r", newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            key = (r.get("RUN_DIR") or "").strip()
            if key:
                out[key] = r
    return out

def write_csv(csv_path: Path, rows: Iterable[Dict[str, Any]]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=ALL_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in ALL_COLS})

def scan_run_dirs(results_dir: Path) -> List[Path]:
    runs = []
    for p in results_dir.rglob("meta.json"):
        d = p.parent
        parts = set(d.parts)
        if "_stage1_best" in parts or "_stage2_best" in parts:
            continue
        runs.append(d)
    runs.sort(key=lambda x: str(x))
    return runs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    csv_path = Path(args.csv).resolve() if args.csv else (results_dir / "runlog_compact.csv")

    if not results_dir.exists():
        raise SystemExit(f"results-dir not found: {results_dir}")

    existing = read_existing(csv_path)
    run_dirs = scan_run_dirs(results_dir) if args.scan else []

    out_rows: List[Dict[str, Any]] = []
    if existing and not args.refresh:
        for key, r in existing.items():
            out_rows.append({c: r.get(c, "") for c in ALL_COLS})

    existing_keys = {r.get("RUN_DIR") for r in out_rows if isinstance(r, dict)}

    for d in run_dirs:
        try:
            key = str(d.relative_to(results_dir)).replace("\\", "/")
        except Exception:
            key = str(d).replace("\\", "/")

        if (not args.refresh) and key in existing_keys:
            continue

        row = extract_row(d, results_dir)

        replaced = False
        for i, rr in enumerate(out_rows):
            if (rr.get("RUN_DIR") or "") == key:
                out_rows[i] = row
                replaced = True
                break
        if not replaced:
            out_rows.append(row)
        existing_keys.add(key)

    def sort_key(r):
        s = r.get("RUN_STAMP")
        return (s is None, str(s), str(r.get("RUN_DIR") or ""))
    out_rows.sort(key=sort_key)

    write_csv(csv_path, out_rows)
    print(f"[runlog] wrote {csv_path}  rows={len(out_rows)}  cols={len(ALL_COLS)}")

if __name__ == "__main__":
    main()
