#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path


def load_json(path: Path):
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def to_float(value):
    try:
        if value is None:
            return None
        return float(value)
    except Exception:
        return None


def score_run(run_dir: Path, max_grid_loss_frac: float):
    bm = load_json(run_dir / 'beam_metrics.json') or {}
    meta = load_json(run_dir / 'meta.json') or {}

    currents = bm.get('currents') or {}
    coll = bm.get('collimation') or {}
    defl = bm.get('deflection') or {}
    acc = meta.get('accel_chamfer') or {}

    steer = to_float(defl.get('steer_angle_deg'))
    i_pg = to_float(currents.get('I_pg_in_A'))
    i_ag = to_float(currents.get('I_ag_out_A'))
    grid_loss_frac = to_float(currents.get('grid_loss_frac'))
    grid_tx = to_float(currents.get('grid_transmission_frac'))

    if grid_loss_frac is None and i_pg is not None and i_pg > 0.0 and i_ag is not None:
        grid_tx = max(0.0, min(1.0, i_ag / i_pg))
        grid_loss_frac = max(0.0, 1.0 - grid_tx)

    if steer is None or not math.isfinite(steer):
        return None
    if grid_loss_frac is None:
        return None

    lost_sidewalls = coll.get('lost_to_sidewalls') in (True, 1, '1', 'true')
    safe = (not lost_sidewalls) and (grid_loss_frac <= max_grid_loss_frac)

    return {
        'run_dir': str(run_dir),
        'run_name': run_dir.name,
        'steer_angle_deg': steer,
        'abs_steer_angle_deg': abs(steer),
        'I_pg_in_A': i_pg,
        'I_ag_out_A': i_ag,
        'grid_loss_frac': grid_loss_frac,
        'grid_transmission_frac': grid_tx,
        'lost_to_sidewalls': bool(lost_sidewalls),
        'safe': bool(safe),
        'ACCEL_OFF_Y_M': meta.get('ACCEL_OFF_Y_M'),
        'ACC_DN_DEPTH_M': acc.get('dn_depth'),
        'ACC_DN_ANGLE_DEG': acc.get('dn_angle_deg'),
        'GRID_T_M': meta.get('GRID_T_M'),
        'GAP_M': meta.get('GAP_M'),
    }


def main():
    parser = argparse.ArgumentParser(description='Find the maximum steering angle from two-grid runs without grid loss.')
    parser.add_argument('--results', default='results', help='Results directory root')
    parser.add_argument('--topk', type=int, default=20, help='Number of ranked runs to keep')
    parser.add_argument('--max-grid-loss-frac', type=float, default=1.0e-3,
                        help='Maximum allowed PG-to-AG current loss fraction')
    parser.add_argument('--out', default=None, help='Optional output JSON path')
    args = parser.parse_args()

    results_root = Path(args.results)
    rows = []
    for meta_path in results_root.rglob('meta.json'):
        scored = score_run(meta_path.parent, args.max_grid_loss_frac)
        if scored is not None:
            rows.append(scored)

    safe_rows = [r for r in rows if r['safe']]
    safe_rows.sort(key=lambda r: (-r['abs_steer_angle_deg'], r['grid_loss_frac'], -(r['I_ag_out_A'] or 0.0)))

    best_by_chamfer = {}
    for row in safe_rows:
        key = (row['ACC_DN_DEPTH_M'], row['ACC_DN_ANGLE_DEG'])
        if key not in best_by_chamfer:
            best_by_chamfer[key] = row
    chamfer_rows = list(best_by_chamfer.values())
    chamfer_rows.sort(key=lambda r: (-r['abs_steer_angle_deg'], r['grid_loss_frac'], -(r['I_ag_out_A'] or 0.0)))

    payload = {
        'results_root': str(results_root),
        'max_grid_loss_frac': args.max_grid_loss_frac,
        'num_runs_scanned': len(rows),
        'num_safe_runs': len(safe_rows),
        'best_safe_run': safe_rows[0] if safe_rows else None,
        'top_safe_runs': safe_rows[:args.topk],
        'best_per_accel_chamfer': chamfer_rows[:args.topk],
    }

    out_path = Path(args.out) if args.out else (results_root / 'max_steer_search.json')
    out_path.write_text(json.dumps(payload, indent=2))

    print(f'scanned: {len(rows)} runs')
    print(f'safe:    {len(safe_rows)} runs (grid_loss_frac <= {args.max_grid_loss_frac:g}, no sidewall loss)')
    if not safe_rows:
        print('no safe runs found')
        return 0

    best = safe_rows[0]
    print('best safe run:')
    print(f"  run: {best['run_name']}")
    print(f"  steer_angle_deg: {best['steer_angle_deg']:.6g}")
    print(f"  grid_loss_frac: {best['grid_loss_frac']:.6g}")
    print(f"  accel_off_y_m: {best['ACCEL_OFF_Y_M']}")
    print(f"  acc_dn_depth_m: {best['ACC_DN_DEPTH_M']}")
    print(f"  acc_dn_angle_deg: {best['ACC_DN_ANGLE_DEG']}")
    print(f'json: {out_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
