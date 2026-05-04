#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

try:
    import pandas as pd
    _HAS_PANDAS = True
except Exception:
    pd = None
    _HAS_PANDAS = False

try:
    import plotly.graph_objects as go
    _HAS_PLOTLY = True
except Exception:
    go = None
    _HAS_PLOTLY = False


def _split_floats(text: str):
    out = []
    for tok in (text or '').split(','):
        tok = tok.strip()
        if not tok:
            continue
        out.append(float(tok))
    return out


def _load_csv_lines(args):
    if not _HAS_PANDAS:
        raise SystemExit('pandas is required for --csv mode')

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f'csv not found: {csv_path}')

    df = pd.read_csv(csv_path)
    if args.filter:
        try:
            df = df.query(args.filter, engine='python')
        except Exception as exc:
            raise SystemExit(f'bad --filter expression: {exc}') from exc

    needed = [args.group_col, args.gap_col, args.steer_col]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f'missing required columns: {", ".join(missing)}')

    g = pd.to_numeric(df[args.group_col], errors='coerce')
    gap = pd.to_numeric(df[args.gap_col], errors='coerce')
    steer = pd.to_numeric(df[args.steer_col], errors='coerce')
    work = df.assign(_group=g, _gap=gap, _steer=steer).dropna(subset=['_group', '_gap', '_steer'])
    if work.empty:
        raise SystemExit('no valid numeric rows left after CSV filtering')

    if args.abs_steer:
        work['_score'] = work['_steer'].abs()
    else:
        work['_score'] = work['_steer']

    idx = work.groupby('_group')['_score'].idxmax()
    chosen = work.loc[idx].sort_values('_group', kind='mergesort')
    if chosen.empty:
        raise SystemExit('no rows selected from CSV')

    lines = []
    for _, row in chosen.iterrows():
        radius_mm = row['_group'] * 1000.0
        gap_mm = row['_gap'] * 1000.0
        steer_deg = float(row['_steer'])
        label = f'r={radius_mm:g} mm, gap={gap_mm:g} mm'
        lines.append((steer_deg, label))
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description='Simple geometric steering estimate plot')
    ap.add_argument('--angles-deg', default='1,2,3,4', help='comma-separated steering angles in degrees')
    ap.add_argument('--csv', default=None, help='optional runlog_compact.csv to select one max-steering case per radius')
    ap.add_argument('--filter', default=None, help='optional pandas query applied in --csv mode')
    ap.add_argument('--group-col', default='AP_RAD_M', help='radius column in --csv mode')
    ap.add_argument('--gap-col', default='GAP_M', help='gap column in --csv mode')
    ap.add_argument('--steer-col', default='STEER_ANGLE_DEG', help='steering-angle column in --csv mode')
    ap.add_argument('--abs-steer', action='store_true', help='rank CSV rows by absolute steering angle')
    ap.add_argument('--x0-m', type=float, default=0.0, help='trajectory origin in meters')
    ap.add_argument('--x-max-m', type=float, default=0.55, help='maximum downstream distance in meters')
    ap.add_argument('--npts', type=int, default=300, help='number of x samples')
    ap.add_argument('--x-units', default='mm', choices=['m', 'mm'], help='x-axis units')
    ap.add_argument('--y-units', default='mm', choices=['m', 'mm'], help='y-axis units')
    ap.add_argument('--x-label', default='Range Down Beam Tube (mm)', help='x-axis label')
    ap.add_argument('--y-label', default='Estimated Beam Center Offset (mm)', help='y-axis label')
    ap.add_argument('--vlines-mm', default='', help='comma-separated vertical reference lines in mm')
    ap.add_argument('--title', default='Estimated Beam-Center Trajectory from Steering Angle', help='plot title')
    ap.add_argument('--out', default='steer_estimate.html', help='output html path')
    args = ap.parse_args(argv)

    if not _HAS_PLOTLY:
        raise SystemExit('plotly is required (pip install plotly)')

    if args.csv:
        lines = _load_csv_lines(args)
        angles = [angle for angle, _ in lines]
    else:
        angles = _split_floats(args.angles_deg)
        lines = [(angle, f'{angle:g} deg') for angle in angles]
    if not angles:
        raise SystemExit('no angles provided')
    if args.npts < 2:
        raise SystemExit('--npts must be at least 2')
    if args.x_max_m <= args.x0_m:
        raise SystemExit('--x-max-m must be greater than --x0-m')

    xscale = 1000.0 if args.x_units == 'mm' else 1.0
    yscale = 1000.0 if args.y_units == 'mm' else 1.0

    xs_m = [args.x0_m + i * (args.x_max_m - args.x0_m) / (args.npts - 1) for i in range(args.npts)]

    fig = go.Figure()
    for ang, label in lines:
        slope = math.tan(math.radians(ang))
        ys_m = [(x - args.x0_m) * slope for x in xs_m]
        fig.add_trace(go.Scatter(
            x=[x * xscale for x in xs_m],
            y=[y * yscale for y in ys_m],
            mode='lines',
            name=label,
            hovertemplate=(
                f'{args.x_label}=%{{x:.6g}}<br>'
                f'{args.y_label}=%{{y:.6g}}<br>'
                f'label={label}<br>'
                f'steering={ang:g} deg<extra></extra>'
            ),
        ))

    for xv in _split_floats(args.vlines_mm):
        fig.add_vline(
            x=xv,
            line_width=1,
            line_dash='dot',
            line_color='gray',
            opacity=0.75,
            annotation_text=f'{xv:g} mm',
            annotation_position='top'
        )

    fig.update_layout(
        template='plotly_white',
        title=args.title,
        xaxis_title=args.x_label,
        yaxis_title=args.y_label,
        legend_title='Steering Angle',
        font=dict(size=20),
        title_font=dict(size=30),
        legend=dict(font=dict(size=18), title_font=dict(size=20)),
    )
    fig.update_xaxes(title_font=dict(size=24), tickfont=dict(size=18))
    fig.update_yaxes(title_font=dict(size=24), tickfont=dict(size=18))
    out = Path(args.out)
    fig.write_html(out)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
