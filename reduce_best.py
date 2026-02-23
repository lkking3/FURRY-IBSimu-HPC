#!/usr/bin/env python3
"""
Reduce Stage-1 sweep outputs to the best candidate(s) using beam diagnostics.

Ranking (most important first):
  1) |y_mean_pred_500mm_m| closest to 0
  2) divergence_angle_deg smallest
  3) I_ag_out_A largest

Disqualify runs if collimation.lost_to_sidewalls is true.

Env vars:
  RESULTS_DIR       : root results folder (default: ./results)
  TOPK              : how many top rows to keep (default: 10)
  OPT_KNOBS         : comma list to seed Stage-2; if empty, auto-pick
  MAX_KNOBS         : limit for auto-picked knobs (default: 4)
  DELTA_FLOORS_JSON : JSON object of min step sizes (optional)

Outputs under:
  RESULTS_DIR/_stage1_best/<timestamp>/
    - best.json
    - seed_stage2.env
    - topk.json
    - two_grid_geom.png (if present)
    - two_grid_geom_closeup.png (if present)
  RESULTS_DIR/_stage1_best/latest -> symlink to the stamp folder
  RESULTS_DIR/_stage1_best/BEST_PATH.txt
  RESULTS_DIR/_stage1_best/latest/knobs_stage2.env
"""

import json
import math
import os
import shutil
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

RESULTS = Path(os.environ.get("RESULTS_DIR", "./results"))
TOPK = int(os.environ.get("TOPK", "10"))
EXPLICIT = os.environ.get("OPT_KNOBS", "").strip()
MAX_KNOBS = int(os.environ.get("MAX_KNOBS", "4"))
DELTA_FLOORS = json.loads(os.environ.get("DELTA_FLOORS_JSON", "{}"))


def load_json(p: Path) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def to_float(v: Any) -> Optional[float]:
    try:
        f = float(v)
    except Exception:
        return None
    if not math.isfinite(f):
        return None
    return f


def score_run(run: Path) -> Optional[Tuple[Tuple[float, float, float], Dict[str, Any]]]:
    """
    Score a run using beam diagnostics:
      1) |y_mean_pred_500mm_m| closest to 0
      2) divergence_angle_deg smallest
      3) I_ag_out_A largest
    """
    bm = load_json(run / "beam_metrics.json") or {}
    coll = bm.get("collimation") or {}
    if coll.get("lost_to_sidewalls") in (True, 1, "1", "true"):
        return None

    defl = bm.get("deflection") or {}
    cur = bm.get("currents") or {}

    y_pred = to_float(defl.get("y_mean_pred_500mm_m"))
    div_deg = bm.get("DIVERGENCE_ANGLE_DEG")
    if div_deg is None:
        div_deg = bm.get("divergence_angle_deg")
    div_deg = to_float(div_deg)
    i_ag = to_float(cur.get("I_ag_out_A"))

    if y_pred is None or div_deg is None or i_ag is None:
        return None

    abs_y = abs(y_pred)
    score = (abs_y, div_deg, -i_ag)
    info = {
        "y_mean_pred_500mm_m": y_pred,
        "abs_y_mean_pred_500mm_m": abs_y,
        "divergence_angle_deg": div_deg,
        "I_ag_out_A": i_ag,
        "lost_to_sidewalls": bool(coll.get("lost_to_sidewalls", False)),
    }
    return score, info


def pick_pngs(run_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    mains = list(run_dir.glob("*.png"))
    main = next((p for p in mains if p.name == "two_grid_geom.png"), None)
    if not main:
        main = next((p for p in mains if not p.stem.endswith("_closeup")), None)
    close = next((p for p in mains if p.stem.endswith("_closeup")), None)
    return main, close


def extract_geom(meta: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    if not meta:
        return out
    scr = meta.get("screen_chamfer") or {}
    acc = meta.get("accel_chamfer") or {}
    for k in [
        "AP_RAD_M",
        "SCREEN_AP_RAD_M",
        "ACCEL_AP_RAD_M",
        "GRID_T_M",
        "GAP_M",
        "ACCEL_OFF_Y_M",
        "SCREEN_OFF_Y_M",
        "H",
        "X_LEFT_M",
        "X_RIGHT_M",
        "YBOX_M",
        "VS_V",
        "VA_V",
        "SAMPLE_V",
    ]:
        if k in meta:
            out[k] = meta[k]

    out.update(
        {
            "SCR_UP_DEPTH_M": scr.get("up_depth"),
            "SCR_UP_ANGLE_DEG": scr.get("up_angle_deg"),
            "SCR_DN_DEPTH_M": scr.get("dn_depth"),
            "SCR_DN_ANGLE_DEG": scr.get("dn_angle_deg"),
            "ACC_UP_DEPTH_M": acc.get("up_depth"),
            "ACC_UP_ANGLE_DEG": acc.get("up_angle_deg"),
            "ACC_DN_DEPTH_M": acc.get("dn_depth"),
            "ACC_DN_ANGLE_DEG": acc.get("dn_angle_deg"),
        }
    )
    return out


def collect_rows(root: Path):
    rows = []
    for meta in root.rglob("meta.json"):
        run = meta.parent
        scored = score_run(run)
        if not scored:
            continue
        _, info = scored
        y = info.get("abs_y_mean_pred_500mm_m")
        if y is None:
            continue

        m = load_json(meta) or {}
        scr = m.get("screen_chamfer") or {}
        acc = m.get("accel_chamfer") or {}
        feat = {
            "AP_RAD_M": m.get("AP_RAD_M"),
            "SCREEN_AP_RAD_M": m.get("SCREEN_AP_RAD_M"),
            "ACCEL_AP_RAD_M": m.get("ACCEL_AP_RAD_M"),
            "GRID_T_M": m.get("GRID_T_M"),
            "GAP_M": m.get("GAP_M"),
            "ACCEL_OFF_Y_M": m.get("ACCEL_OFF_Y_M"),
            "SCREEN_OFF_Y_M": m.get("SCREEN_OFF_Y_M"),
            "SCR_UP_DEPTH_M": scr.get("up_depth"),
            "SCR_UP_ANGLE_DEG": scr.get("up_angle_deg"),
            "SCR_DN_DEPTH_M": scr.get("dn_depth"),
            "SCR_DN_ANGLE_DEG": scr.get("dn_angle_deg"),
            "ACC_UP_DEPTH_M": acc.get("up_depth"),
            "ACC_UP_ANGLE_DEG": acc.get("up_angle_deg"),
            "ACC_DN_DEPTH_M": acc.get("dn_depth"),
            "ACC_DN_ANGLE_DEG": acc.get("dn_angle_deg"),
        }
        if any(v is None for v in feat.values()):
            continue
        rows.append((feat, float(y), run))
    return rows


def suggest_knobs(rows, best_geom, top):
    if EXPLICIT:
        knob_list = [k.strip() for k in EXPLICIT.split(",") if k.strip()]
    else:
        try:
            import numpy as np

            cols = list(rows[0][0].keys())
            X = np.array([[r[0][c] for c in cols] for r in rows], float)
            y = np.array([r[1] for r in rows], float)
            mu = X.mean(0)
            sig = X.std(0) + 1e-12
            Xz = (X - mu) / sig
            lam = 1e-6
            w = np.linalg.solve(Xz.T @ Xz + lam * np.eye(Xz.shape[1]), Xz.T @ y)
            imp = np.abs(w)
            rank = np.argsort(-imp)
            knob_list = [cols[i] for i in rank[:MAX_KNOBS]]
        except Exception:
            cols = [k for k in top[0].keys() if k.endswith("_M") or k.endswith("_DEG")]
            var = {
                c: (
                    statistics.pstdev([float(r[c]) for r in top if r.get(c) is not None])
                    if sum(1 for r in top if r.get(c) is not None) > 1
                    else 0.0
                )
                for c in cols
            }
            knob_list = [k for k, _ in sorted(var.items(), key=lambda kv: -kv[1])[:MAX_KNOBS]]

    # physics guardrail: angle only if its depth > 0
    pairs = [
        ("SCR_UP_ANGLE_DEG", "SCR_UP_DEPTH_M"),
        ("SCR_DN_ANGLE_DEG", "SCR_DN_DEPTH_M"),
        ("ACC_UP_ANGLE_DEG", "ACC_UP_DEPTH_M"),
        ("ACC_DN_ANGLE_DEG", "ACC_DN_DEPTH_M"),
    ]
    for ang, dep in pairs:
        if ang in knob_list and (best_geom.get(dep, 0.0) <= 0.0):
            knob_list.remove(ang)
            if dep not in knob_list:
                knob_list.append(dep)
    return knob_list


def suggest_deltas(knobs, top, floors):
    deltas = {}
    for k in knobs:
        vals = [t.get(k) for t in top if t.get(k) is not None]
        spread = (max(vals) - min(vals)) / 4.0 if len(vals) >= 2 else 0.0
        floor = floors.get(k, 0.0)
        deltas[k] = max(spread, floor)
    return deltas


def main() -> int:
    if not RESULTS.exists():
        print(f"[reduce] RESULTS_DIR not found: {RESULTS}", file=sys.stderr)
        return 2

    items = []
    for meta in RESULTS.rglob("meta.json"):
        run = meta.parent
        scored = score_run(run)
        if not scored:
            continue
        score, info = scored
        items.append({"run": run, "meta": load_json(meta), "beam": info, "score": score})

    if not items:
        print("[reduce] no runs with usable beam metrics found", file=sys.stderr)
        return 3

    items.sort(key=lambda r: r["score"])
    stamp = time.strftime("%Y%m%d_%H%M%S")
    best_root = RESULTS / "_stage1_best"
    outdir = best_root / stamp
    outdir.mkdir(parents=True, exist_ok=True)

    top = []
    for r in items[: min(TOPK, len(items))]:
        row = {
            "run_dir": str(r["run"]),
            "score_keys": [
                "abs_y_mean_pred_500mm_m",
                "divergence_angle_deg",
                "-I_ag_out_A",
            ],
            "score_values": [
                r["beam"]["abs_y_mean_pred_500mm_m"],
                r["beam"]["divergence_angle_deg"],
                -r["beam"]["I_ag_out_A"],
            ],
            "beam_metrics": r["beam"],
        }
        row.update(extract_geom(r["meta"]))
        top.append(row)

    (outdir / "topk.json").write_text(
        json.dumps(
            {
                "ranking": "abs(y_mean_pred_500mm_m), divergence_angle_deg, -I_ag_out_A",
                "topk": top,
            },
            indent=2,
        )
    )

    best = items[0]
    geom = extract_geom(best["meta"])
    (outdir / "best.json").write_text(
        json.dumps(
            {
                "ranking": "abs(y_mean_pred_500mm_m), divergence_angle_deg, -I_ag_out_A",
                "score_values": [
                    best["beam"]["abs_y_mean_pred_500mm_m"],
                    best["beam"]["divergence_angle_deg"],
                    -best["beam"]["I_ag_out_A"],
                ],
                "run_dir": str(best["run"]),
                "beam_metrics": best["beam"],
                "geom": geom,
            },
            indent=2,
        )
    )

    seed = outdir / "seed_stage2.env"
    with seed.open("w") as f:
        for k, v in geom.items():
            if v is not None:
                f.write(f"export {k}={v}\n")

    main_png, close_png = pick_pngs(best["run"])
    if main_png:
        shutil.copy2(main_png, outdir / main_png.name)
    if close_png:
        shutil.copy2(close_png, outdir / close_png.name)

    # Update latest symlink and path helper.
    try:
        latest = best_root / "latest"
        if latest.exists() or latest.is_symlink():
            latest.unlink()
        os.symlink(stamp, latest)
    except Exception:
        pass
    (best_root / "BEST_PATH.txt").write_text(str(outdir) + "\n")

    rows = collect_rows(RESULTS)
    if rows and top:
        knobs = suggest_knobs(rows, geom, top)
        deltas = suggest_deltas(knobs, top, DELTA_FLOORS)
        knobs_env = best_root / "latest" / "knobs_stage2.env"
        with knobs_env.open("w") as f:
            f.write("export OPT_KNOBS=" + ",".join(knobs) + "\n")
            for k, d in deltas.items():
                f.write(f"export DELTA_{k}={d}\n")

        print(f"[reduce] knobs: {knobs} with deltas {deltas}")

    print(f"[reduce] best: {best['run']}  score={best['score']}")
    print(f"[reduce] seed: {seed}")
    print(f"[reduce] out : {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
