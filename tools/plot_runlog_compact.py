#!/usr/bin/env python3
"""
plot_runlog_compact.py

Interactive plots for runlog_compact.csv (scatter + heatmap).
Uses pandas + plotly to generate self-contained HTML files.
"""
import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable, Optional, List, Tuple, Dict, Any

import pandas as pd
import numpy as np

try:
    import plotly.express as px
    import plotly.graph_objects as go
    _HAS_PLOTLY = True
except Exception:
    px = None
    go = None
    _HAS_PLOTLY = False

try:
    from scipy.ndimage import gaussian_filter
except Exception:  # pragma: no cover
    gaussian_filter = None


def _load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"csv not found: {path}")
    df = pd.read_csv(path)
    # Coerce numeric-like columns; leave pure strings intact.
    for col in df.columns:
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.notna().any():
            df[col] = coerced
    return df


def _apply_filter(df: pd.DataFrame, expr: Optional[str]) -> pd.DataFrame:
    if not expr:
        return df
    try:
        return df.query(expr, engine="python")
    except Exception as exc:
        raise SystemExit(f"bad --filter expression: {exc}") from exc


def _safe_name(text: str) -> str:
    out = []
    for ch in text:
        if ch.isalnum() or ch in ("-", "_", "."):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out)


def _load_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _as_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        v = float(value)
    except Exception:
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def _add_suffix_before_ext(path: Path, suffix: str) -> Path:
    name = path.name
    if "." in name:
        stem, ext = name.rsplit(".", 1)
        return path.with_name(f"{stem}{suffix}.{ext}")
    return path.with_name(f"{name}{suffix}")


def _split_list(text: Optional[str]) -> List[str]:
    if not text:
        return []
    return [t.strip() for t in text.split(",") if t.strip()]


def _auto_bins(series: pd.Series, bins: Optional[int]) -> Optional[int]:
    if bins is not None:
        return bins
    if not pd.api.types.is_numeric_dtype(series):
        return None
    nunique = series.nunique(dropna=True)
    if nunique > 50:
        return 50
    return None


def _bin_series(series: pd.Series, bins: Optional[int]) -> pd.Series:
    if bins is None:
        return series
    if not pd.api.types.is_numeric_dtype(series):
        return series
    cut = pd.cut(series, bins=bins)
    return cut.astype(str)


def _gaussian_smooth(z: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0.0:
        return z
    mask = np.isfinite(z).astype(float)
    z0 = np.nan_to_num(z, nan=0.0)
    zf = gaussian_filter(z0, sigma=sigma, mode="nearest")
    wf = gaussian_filter(mask, sigma=sigma, mode="nearest")
    return zf / np.maximum(wf, 1e-12)


_IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".svg", ".pdf"}


def _write_figure(fig: "go.Figure", out: Path) -> Path:
    if not _HAS_PLOTLY:
        raise RuntimeError("plotly is required for HTML/image export")
    if out.suffix.lower() in _IMAGE_EXTS:
        try:
            fig.write_image(out)
            return out
        except Exception as exc:
            fallback = out.with_suffix(".html")
            fig.write_html(fallback)
            print(f"warning: image export failed ({exc}); wrote {fallback} instead")
            return fallback
    fig.write_html(out)
    return out


def _apply_presentation_fonts(fig: "go.Figure") -> None:
    if not _HAS_PLOTLY:
        return
    fig.update_layout(
        font=dict(size=20),
        title_font=dict(size=30),
        legend=dict(font=dict(size=18), title_font=dict(size=20)),
    )
    fig.update_xaxes(title_font=dict(size=24), tickfont=dict(size=18))
    fig.update_yaxes(title_font=dict(size=24), tickfont=dict(size=18))


def _summary(df: pd.DataFrame, top: int) -> None:
    print(f"rows={len(df)} cols={len(df.columns)}")
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    cat_cols = [c for c in df.columns if c not in num_cols]

    print("\nNumeric columns (top by unique count):")
    rows = []
    for col in num_cols:
        s = df[col].dropna()
        if s.empty:
            continue
        rows.append((col, s.nunique(), s.min(), s.max(), s.mean()))
    rows.sort(key=lambda r: r[1], reverse=True)
    for col, nunique, vmin, vmax, mean in rows[:top]:
        print(f"  {col}: unique={nunique} min={vmin:.6g} max={vmax:.6g} mean={mean:.6g}")

    print("\nCategorical columns (top by unique count):")
    rows = []
    for col in cat_cols:
        s = df[col].dropna()
        if s.empty:
            continue
        rows.append((col, s.nunique()))
    rows.sort(key=lambda r: r[1], reverse=True)
    for col, nunique in rows[:top]:
        print(f"  {col}: unique={nunique}")


def _scatter(df: pd.DataFrame, args: argparse.Namespace) -> None:
    if not _HAS_PLOTLY:
        raise SystemExit("plotly is required for scatter plots (pip install plotly)")
    hover = _split_list(args.hover) or None
    cols = [args.x, args.y] + [c for c in (args.color, args.size, args.facet_row, args.facet_col) if c]
    df = df.dropna(subset=[c for c in cols if c in df.columns])

    fig = px.scatter(
        df,
        x=args.x,
        y=args.y,
        color=args.color,
        size=args.size,
        hover_data=hover,
        facet_row=args.facet_row,
        facet_col=args.facet_col,
        title=args.title,
    )
    if args.logx:
        fig.update_xaxes(type="log")
    if args.logy:
        fig.update_yaxes(type="log")
    fig.update_layout(template="plotly_white")
    _apply_presentation_fonts(fig)

    out = Path(args.out) if args.out else Path(f"scatter_{_safe_name(args.x)}_vs_{_safe_name(args.y)}.html")
    _apply_presentation_fonts(fig)
    fig.write_html(out)
    print(f"wrote {out}")
    if args.show:
        fig.show()


def _heatmap(df: pd.DataFrame, args: argparse.Namespace) -> None:
    if not _HAS_PLOTLY:
        raise SystemExit("plotly is required for heatmaps (pip install plotly)")
    xbins = _auto_bins(df[args.x], args.xbins)
    ybins = _auto_bins(df[args.y], args.ybins)

    df = df.dropna(subset=[args.x, args.y, args.z])
    x_key = _bin_series(df[args.x], xbins)
    y_key = _bin_series(df[args.y], ybins)

    pivot = pd.pivot_table(
        df.assign(_x=x_key, _y=y_key),
        index="_y",
        columns="_x",
        values=args.z,
        aggfunc=args.agg,
    )

    z = pivot.values.astype(float)
    if args.smooth:
        if gaussian_filter is None:
            raise SystemExit("scipy is required for --smooth")
        z = _gaussian_smooth(z, args.smooth)

    zsmooth = False if args.zsmooth == "none" else args.zsmooth
    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=pivot.columns.astype(str),
            y=pivot.index.astype(str),
            zsmooth=zsmooth,
            colorbar=dict(title=args.z),
        )
    )
    fig.update_layout(
        title=args.title,
        xaxis_title=x_label,
        yaxis_title=args.y,
        template="plotly_white",
    )
    _apply_presentation_fonts(fig)
    if args.logz:
        fig.update_coloraxes(type="log")

    out = Path(args.out) if args.out else Path(f"heatmap_{_safe_name(args.z)}_by_{_safe_name(args.x)}_{_safe_name(args.y)}.html")
    fig.write_html(out)
    print(f"wrote {out}")
    if args.show:
        fig.show()


def _format_group_value(value: Any, units: str) -> str:
    if value is None:
        return "unknown"
    fv = _as_float(value)
    if fv is None:
        return str(value)
    scale = 1.0e3 if units == "mm" else 1.0
    suffix = " mm" if units == "mm" else " m"
    return f"{fv * scale:.6g}{suffix}"


def _segment_indices(mask: List[bool]) -> List[Tuple[int, int, bool]]:
    if not mask:
        return []
    out: List[Tuple[int, int, bool]] = []
    start = 0
    state = bool(mask[0])
    for i in range(1, len(mask)):
        cur = bool(mask[i])
        if cur != state:
            out.append((start, i, state))
            start = i
            state = cur
    out.append((start, len(mask), state))
    return out


def _steer_curves(df: pd.DataFrame, args: argparse.Namespace) -> None:
    if not _HAS_PLOTLY:
        raise SystemExit("plotly is required for steer-curves plots (pip install plotly)")

    needed = [args.x, args.y, args.group]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"missing required columns: {', '.join(missing)}")

    plot_df = df.copy()
    if args.loss_col not in plot_df.columns:
        plot_df[args.loss_col] = 0.0
    if args.loss_abs_col and args.loss_abs_col not in plot_df.columns:
        plot_df[args.loss_abs_col] = np.nan

    cols = [args.x, args.y, args.group, args.loss_col]
    plot_df = plot_df.dropna(subset=[c for c in cols if c in plot_df.columns])
    if plot_df.empty:
        raise SystemExit("no rows left after filtering and dropping NaNs")

    xvals = pd.to_numeric(plot_df[args.x], errors="coerce")
    yvals = pd.to_numeric(plot_df[args.y], errors="coerce")
    groupvals = plot_df[args.group]
    lossvals = pd.to_numeric(plot_df[args.loss_col], errors="coerce").fillna(0.0)
    plot_df = plot_df.assign(_x=xvals, _y=yvals, _group=groupvals, _loss=lossvals)
    plot_df = plot_df.dropna(subset=["_x", "_y", "_group"])
    if plot_df.empty:
        raise SystemExit("no valid numeric rows left for steer-curves plot")

    if args.abs_y:
        plot_df["_y"] = plot_df["_y"].abs()

    if args.require_no_sidewalls and "LOST_TO_SIDEWALLS" in plot_df.columns:
        keep = ~plot_df["LOST_TO_SIDEWALLS"].astype(str).str.lower().isin(["true", "1"])
        plot_df = plot_df[keep]
    if plot_df.empty:
        raise SystemExit("no rows left after LOST_TO_SIDEWALLS filtering")

    plot_df = plot_df.sort_values(["_group", "_x", "_loss"], kind="mergesort")

    if getattr(args, "debug_summary", False):
        parts = [
            f"rows={len(plot_df)}",
            f"groups={sorted(pd.Series(plot_df['_group']).dropna().unique().tolist())}",
            f"x_min={float(plot_df['_x'].min()):.6g}",
            f"x_max={float(plot_df['_x'].max()):.6g}",
        ]
        if 'GAP_M' in plot_df.columns:
            gaps = sorted(pd.to_numeric(plot_df['GAP_M'], errors='coerce').dropna().unique().tolist())
            parts.append(f"gaps={gaps}")
        print('[steer-curves] ' + ' '.join(parts))

    xscale = 1000.0 if args.x_units == "mm" else 1.0
    x_label = f"{args.x} (mm)" if args.x_units == "mm" else f"{args.x} (m)"
    plot_df["_x_plot"] = plot_df["_x"] * xscale

    colorway = px.colors.qualitative.Plotly if px is not None else ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    fig = go.Figure()

    hover_cols = [c for c in _split_list(args.hover) if c in plot_df.columns]
    groups = list(dict.fromkeys(plot_df["_group"].tolist()))
    for gi, group_value in enumerate(groups):
        gdf = plot_df[plot_df["_group"] == group_value].copy()
        gdf = gdf.sort_values("_x", kind="mergesort")
        if gdf.empty:
            continue
        color = colorway[gi % len(colorway)]
        label = _format_group_value(group_value, args.group_units)
        lossy_mask = (gdf["_loss"] > args.loss_threshold).tolist()
        segments = _segment_indices(lossy_mask)
        custom = np.column_stack([
            gdf[args.group].astype(str).to_numpy(),
            gdf["_loss"].to_numpy(),
        ])
        if args.loss_abs_col and args.loss_abs_col in gdf.columns:
            loss_abs_vals = pd.to_numeric(gdf[args.loss_abs_col], errors="coerce").to_numpy()
        else:
            loss_abs_vals = np.full(len(gdf), np.nan)
        custom = np.column_stack([custom, loss_abs_vals])
        for hover_col in hover_cols:
            custom = np.column_stack([custom, gdf[hover_col].astype(str).to_numpy()])

        hovertemplate = (
            f"{x_label}=%{{x:.6g}}<br>"
            f"{args.y}=%{{y:.6g}}<br>"
            f"{args.group}=%{{customdata[0]}}<br>"
            f"{args.loss_col}=%{{customdata[1]:.3g}}<br>"
            f"{args.loss_abs_col or 'grid_loss_A'}=%{{customdata[2]:.3g}}"
        )
        for hi, hover_col in enumerate(hover_cols, start=3):
            hovertemplate += f"<br>{hover_col}=%{{customdata[{hi}]}}"
        hovertemplate += "<extra></extra>"

        showlegend = True
        for seg_start, seg_end, is_lossy in segments:
            seg = gdf.iloc[seg_start:seg_end]
            fig.add_trace(go.Scatter(
                x=seg["_x_plot"],
                y=seg["_y"],
                mode="lines+markers",
                name=label,
                legendgroup=label,
                showlegend=showlegend,
                line=dict(color=color, width=3, dash=("dash" if is_lossy else "solid")),
                marker=dict(color=color, size=8, symbol=("x" if is_lossy else "circle"), opacity=(0.75 if is_lossy else 0.95)),
                customdata=custom[seg_start:seg_end],
                hovertemplate=hovertemplate,
                opacity=(0.65 if is_lossy else 1.0),
            ))
            showlegend = False

        if args.show_threshold:
            threshold_rows = gdf[gdf["_loss"] > args.loss_threshold]
            if not threshold_rows.empty:
                first_bad = threshold_rows.iloc[0]
                fig.add_vline(
                    x=float(first_bad["_x_plot"]),
                    line_width=1,
                    line_dash="dot",
                    line_color=color,
                    opacity=0.45,
                )

    y_label = args.y_label or (f"abs({args.y})" if args.abs_y else args.y)
    fig.update_layout(
        template="plotly_white",
        title=args.title or f"{args.y} vs {args.x} grouped by {args.group}",
        xaxis_title=(args.x_label or x_label),
        yaxis_title=y_label,
        legend_title=args.group,
    )
    _apply_presentation_fonts(fig)
    if args.logx:
        fig.update_xaxes(type="log")
    if args.logy:
        fig.update_yaxes(type="log")

    out = Path(args.out) if args.out else Path(f"steer_curves_{_safe_name(args.group)}.html")
    fig = fig
    written = _write_figure(fig, out)
    print(f"wrote {written}")
    if args.show:
        fig.show()


def _max_safe_offset_curves(df: pd.DataFrame, args: argparse.Namespace) -> None:
    if not _HAS_PLOTLY:
        raise SystemExit("plotly is required for max-safe-offset-curves plots (pip install plotly)")

    needed = [args.x, args.offset_col, args.group]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"missing required columns: {', '.join(missing)}")

    plot_df = df.copy()
    if args.loss_col not in plot_df.columns:
        plot_df[args.loss_col] = 0.0

    cols = [args.x, args.offset_col, args.group, args.loss_col]
    plot_df = plot_df.dropna(subset=[c for c in cols if c in plot_df.columns])
    if plot_df.empty:
        raise SystemExit("no rows left after filtering and dropping NaNs")

    xvals = pd.to_numeric(plot_df[args.x], errors="coerce")
    offvals = pd.to_numeric(plot_df[args.offset_col], errors="coerce")
    groupvals = plot_df[args.group]
    lossvals = pd.to_numeric(plot_df[args.loss_col], errors="coerce").fillna(np.inf)
    plot_df = plot_df.assign(_x=xvals, _offset=offvals, _group=groupvals, _loss=lossvals)
    plot_df = plot_df.dropna(subset=["_x", "_offset", "_group"])
    if plot_df.empty:
        raise SystemExit("no valid numeric rows left for max-safe-offset-curves plot")

    if args.require_no_sidewalls and "LOST_TO_SIDEWALLS" in plot_df.columns:
        keep = ~plot_df["LOST_TO_SIDEWALLS"].astype(str).str.lower().isin(["true", "1"])
        plot_df = plot_df[keep]
    if plot_df.empty:
        raise SystemExit("no rows left after LOST_TO_SIDEWALLS filtering")

    safe_df = plot_df[plot_df["_loss"] <= args.loss_threshold].copy()
    if safe_df.empty:
        raise SystemExit("no rows satisfy the safe-loss threshold")

    grouped = (
        safe_df.groupby(["_group", "_x"], dropna=False)
        .agg(
            _offset_max=("_offset", "max"),
            _safe_count=("_offset", "size"),
            _loss_max=("_loss", "max"),
        )
        .reset_index()
        .sort_values(["_group", "_x"], kind="mergesort")
    )

    xscale = 1000.0 if args.x_units == "mm" else 1.0
    yscale = 1000.0 if args.y_units == "mm" else 1.0
    x_label = args.x_label or (f"{args.x} (mm)" if args.x_units == "mm" else f"{args.x} (m)")
    y_label = args.y_label or (f"Max Safe Offset (mm)" if args.y_units == "mm" else f"Max Safe Offset (m)")
    grouped["_x_plot"] = grouped["_x"] * xscale
    grouped["_y_plot"] = grouped["_offset_max"] * yscale

    colorway = px.colors.qualitative.Plotly if px is not None else ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    fig = go.Figure()

    groups = list(dict.fromkeys(grouped["_group"].tolist()))
    for gi, group_value in enumerate(groups):
        gdf = grouped[grouped["_group"] == group_value].copy()
        gdf = gdf.sort_values("_x", kind="mergesort")
        if gdf.empty:
            continue
        color = colorway[gi % len(colorway)]
        label = _format_group_value(group_value, args.group_units)
        custom = np.column_stack([
            gdf["_group"].astype(str).to_numpy(),
            gdf["_safe_count"].to_numpy(),
            gdf["_loss_max"].to_numpy(),
            gdf["_offset_max"].to_numpy(),
        ])
        fig.add_trace(go.Scatter(
            x=gdf["_x_plot"],
            y=gdf["_y_plot"],
            mode="lines+markers",
            name=label,
            legendgroup=label,
            line=dict(color=color, width=3),
            marker=dict(color=color, size=8),
            customdata=custom,
            hovertemplate=(
                f"{x_label}=%{{x:.6g}}<br>"
                f"{y_label}=%{{y:.6g}}<br>"
                f"{args.group}=%{{customdata[0]}}<br>"
                f"safe run count=%{{customdata[1]}}<br>"
                f"max safe {args.loss_col}=%{{customdata[2]:.3g}}<br>"
                f"max safe {args.offset_col}=%{{customdata[3]:.6g}}<extra></extra>"
            ),
        ))

    fig.update_layout(
        template="plotly_white",
        title=args.title or f"Maximum safe {args.offset_col} vs {args.x} grouped by {args.group}",
        xaxis_title=x_label,
        yaxis_title=y_label,
        legend_title=args.group,
    )
    _apply_presentation_fonts(fig)
    if args.logx:
        fig.update_xaxes(type="log")
    if args.logy:
        fig.update_yaxes(type="log")

    out = Path(args.out) if args.out else Path(f"max_safe_offset_curves_{_safe_name(args.group)}.html")
    written = _write_figure(fig, out)
    print(f"wrote {written}")
    if args.show:
        fig.show()



def _beam_centerlines(df: pd.DataFrame, args: argparse.Namespace) -> None:
    if not _HAS_PLOTLY:
        raise SystemExit("plotly is required for beam-centerlines plots (pip install plotly)")

    needed = [args.steer_col, args.path_col]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"missing required columns: {', '.join(missing)}")

    plot_df = df.copy()
    plot_df["_steer"] = pd.to_numeric(plot_df[args.steer_col], errors="coerce")
    plot_df = plot_df.dropna(subset=["_steer", args.path_col])
    if plot_df.empty:
        raise SystemExit("no rows left after dropping NaN steering values")

    if args.abs_steer:
        plot_df["_steer"] = plot_df["_steer"].abs()

    if args.require_no_sidewalls and "LOST_TO_SIDEWALLS" in plot_df.columns:
        keep = ~plot_df["LOST_TO_SIDEWALLS"].astype(str).str.lower().isin(["true", "1"])
        plot_df = plot_df[keep]
    if args.max_loss_frac is not None and args.loss_col in plot_df.columns:
        loss = pd.to_numeric(plot_df[args.loss_col], errors="coerce")
        plot_df = plot_df[loss <= args.max_loss_frac]
    if plot_df.empty:
        raise SystemExit("no rows left after filtering")

    targets = []
    for tok in _split_list(args.targets):
        try:
            targets.append(float(tok))
        except Exception:
            raise SystemExit(f"bad steering target: {tok}")
    if not targets:
        raise SystemExit("no steering targets provided")

    hover_cols = [c for c in _split_list(args.hover) if c in plot_df.columns]
    xscale = 1000.0 if args.x_units == "mm" else 1.0
    yscale = 1000.0 if args.y_units == "mm" else 1.0
    x_label = args.x_label or ("Range Down Beam Tube (mm)" if args.x_units == "mm" else "Range Down Beam Tube (m)")
    y_label = args.y_label or ("Beam Radial Center (mm)" if args.y_units == "mm" else "Beam Radial Center (m)")

    colorway = px.colors.qualitative.Plotly if px is not None else ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    fig = go.Figure()
    used_paths = set()
    added = 0

    for ti, target in enumerate(targets):
        cand = plot_df.assign(_delta=(plot_df["_steer"] - target).abs()).sort_values(["_delta", "_steer"], kind="mergesort")
        chosen = None
        for _, row in cand.iterrows():
            run_dir = str(row[args.path_col])
            if run_dir in used_paths:
                continue
            chosen = row
            break
        if chosen is None:
            continue

        run_dir = Path(str(chosen[args.path_col]))
        bm_path = run_dir / args.beam_json
        bm = _load_json(bm_path)
        if not bm:
            continue
        planes = bm.get("planes") or []
        rows = []
        for p in planes:
            x = _as_float((p or {}).get("x_m"))
            y = _as_float((p or {}).get("y_mean_m"))
            if x is None or y is None:
                continue
            rows.append((x, y))
        if not rows:
            continue
        rows.sort(key=lambda r: r[0])
        xs = [r[0] * xscale for r in rows]
        ys = [r[1] * yscale for r in rows]

        steer_actual = float(chosen["_steer"])
        delta = abs(steer_actual - target)
        label = f"target {target:g} deg, actual {steer_actual:.3g} deg"
        color = colorway[ti % len(colorway)]

        custom_rows = []
        for x_raw, y_raw in rows:
            vals = [target, steer_actual, delta, str(run_dir)]
            for col in hover_cols:
                vals.append(str(chosen[col]))
            custom_rows.append(vals)

        hovertemplate = (
            f"{x_label}=%{{x:.6g}}<br>"
            f"{y_label}=%{{y:.6g}}<br>"
            "target steer=%{customdata[0]:.6g} deg<br>"
            "actual steer=%{customdata[1]:.6g} deg<br>"
            "|delta|=%{customdata[2]:.6g} deg<br>"
            "run=%{customdata[3]}"
        )
        for hi, col in enumerate(hover_cols, start=4):
            hovertemplate += f"<br>{col}=%{{customdata[{hi}]}}"
        hovertemplate += "<extra></extra>"

        fig.add_trace(go.Scatter(
            x=xs,
            y=ys,
            mode="lines+markers",
            name=label,
            legendgroup=label,
            line=dict(color=color, width=3),
            marker=dict(color=color, size=7),
            customdata=np.array(custom_rows, dtype=object),
            hovertemplate=hovertemplate,
        ))
        used_paths.add(str(run_dir))
        added += 1

    if added == 0:
        raise SystemExit("no valid beam_metrics plane traces found for the requested steering targets")

    fig.update_layout(
        template="plotly_white",
        title=args.title or "Beam Centerline vs Range Down Beam Tube",
        xaxis_title=x_label,
        yaxis_title=y_label,
        legend_title="Selected steering",
    )
    _apply_presentation_fonts(fig)
    if args.logx:
        fig.update_xaxes(type="log")
    if args.logy:
        fig.update_yaxes(type="log")

    out = Path(args.out) if args.out else Path("beam_centerlines.html")
    written = _write_figure(fig, out)
    print(f"wrote {written}")
    if args.show:
        fig.show()



def _find_latest_profile(search_dir: Path, kind: str) -> Path:
    if not search_dir.exists():
        raise SystemExit(f"search dir not found: {search_dir}")

    patterns = []
    if kind == "diameter":
        patterns = ["sample_diameter_profile.json"]
    else:
        patterns = ["sample_diameter_profile.json"]

    best_path = None
    best_mtime = None
    for pat in patterns:
        for p in search_dir.rglob(pat):
            try:
                mtime = p.stat().st_mtime
            except Exception:
                continue
            if best_mtime is None or mtime > best_mtime:
                best_mtime = mtime
                best_path = p

    if not best_path:
        raise SystemExit(f"no sample profile json found under {search_dir}")
    return best_path


def _load_sample_profile(path: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]], str]:
    if not path.exists():
        raise SystemExit(f"json not found: {path}")
    data = json.loads(path.read_text())
    bins = data.get("bins", [])
    mode = "unknown"
    if any("y_lo_m" in b for b in bins):
        mode = "diameter"
    elif any("r_lo_m" in b for b in bins):
        mode = "radial"

    rows = []
    sample_rad = None
    try:
        sample_rad = float(data.get("sample_rad_m"))
    except Exception:
        sample_rad = None
    for b in bins:
        if mode == "diameter":
            x_lo = float(b.get("y_lo_m", 0.0))
            x_hi = float(b.get("y_hi_m", 0.0))
        else:
            x_lo = float(b.get("r_lo_m", 0.0))
            x_hi = float(b.get("r_hi_m", 0.0))
        x_mid = 0.5 * (x_lo + x_hi)
        in_sample = b.get("in_sample")
        if in_sample is None and mode == "diameter" and sample_rad is not None:
            in_sample = abs(x_mid) <= sample_rad
        if in_sample is None:
            in_sample = True
        rows.append(
            {
                "x_lo_m": x_lo,
                "x_hi_m": x_hi,
                "x_mid_m": x_mid,
                "I_Apm": float(b.get("I_Apm", 0.0)),
                "I_A": float(b.get("I_A", 0.0)),
                "fraction": float(b.get("fraction", 0.0)),
                "count": int(b.get("count", 0)),
                "in_sample": bool(in_sample),
            }
        )
    return data, rows, mode


def _write_profile_csv(rows: List[Dict[str, Any]], path: Path, mode: str) -> None:
    head = "y_lo_m,y_hi_m,y_mid_m" if mode == "diameter" else "r_lo_m,r_hi_m,r_mid_m"
    with path.open("w") as f:
        f.write(head + ",I_Apm,I_A,fraction,count\n")
        for r in rows:
            f.write(
                f"{r['x_lo_m']},{r['x_hi_m']},{r['x_mid_m']},"
                f"{r['I_Apm']},{r['I_A']},{r['fraction']},{r['count']}\n"
            )


def _mirror_profile_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    neg = []
    for r in rows:
        neg.append(
            {
                "x_lo_m": -r["x_hi_m"],
                "x_hi_m": -r["x_lo_m"],
                "x_mid_m": -r["x_mid_m"],
                "I_Apm": r["I_Apm"],
                "I_A": r["I_A"],
                "fraction": r["fraction"],
                "count": r["count"],
            }
        )
    return sorted(neg + rows, key=lambda x: x["x_mid_m"])


def _sample_profile(args: argparse.Namespace) -> None:
    if args.json:
        json_path = Path(args.json)
    else:
        json_path = _find_latest_profile(Path(args.search_dir), args.kind)
        print(f"using {json_path}")
    data, rows, mode = _load_sample_profile(json_path)
    if not rows:
        raise SystemExit("no bins found in sample profile json")

    if args.mirror:
        if mode != "radial":
            print("note: --mirror ignored for diameter profiles")
        else:
            rows = _mirror_profile_rows(rows)

    if args.csv:
        _write_profile_csv(rows, Path(args.csv), mode)

    xs = [r["x_mid_m"] for r in rows]
    ys = [r[args.y] for r in rows]
    xscale = 1000.0 if args.x_units == "mm" else 1.0
    axis_name = "r" if mode == "radial" else "y"
    x_label = f"{axis_name} (mm)" if args.x_units == "mm" else f"{axis_name} (m)"
    xs = [x * xscale for x in xs]
    bin_w = data.get("bin_width_m")
    try:
        bin_w = float(bin_w) if bin_w is not None else None
    except Exception:
        bin_w = None
    bar_w = (bin_w * xscale) if bin_w else None

    if args.out or args.show:
        out = Path(args.out) if args.out else Path(f"sample_profile_{_safe_name(args.y)}.html")
        colors = None
        if mode == "diameter" and not args.line:
            colors = ["crimson" if not r.get("in_sample", True) else "royalblue" for r in rows]
        if _HAS_PLOTLY:
            if args.line:
                fig = go.Figure(data=go.Scatter(x=xs, y=ys, mode="lines+markers"))
            else:
                fig = go.Figure(data=go.Bar(x=xs, y=ys, marker_color=colors))

            title = args.title
            if not title:
                title = f"Sample radial profile ({args.y})"
            fig.update_layout(
                title=title,
                xaxis_title=x_label,
                yaxis_title=args.y,
                template="plotly_white",
                margin=dict(l=70, r=30, t=40, b=60),
                font=dict(size=18),
            )
            fig.update_xaxes(title_font=dict(size=30), tickfont=dict(size=20))
            fig.update_yaxes(title_font=dict(size=30), tickfont=dict(size=20))
            if args.logy:
                fig.update_yaxes(type="log")

            actual = _write_figure(fig, out)
            print(f"wrote {actual}")
            if args.show:
                fig.show()
        else:
            if out.suffix.lower() not in _IMAGE_EXTS:
                out = out.with_suffix(".png")
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
            except Exception as exc:
                raise SystemExit(f"plotly not available and matplotlib missing: {exc}") from exc

            plt.figure(figsize=(8, 4.5))
            if args.line:
                plt.plot(xs, ys, marker="o")
            else:
                if bar_w is not None:
                    plt.bar(xs, ys, width=bar_w, align="center", color=colors)
                else:
                    plt.bar(xs, ys, color=colors)
            plt.title(args.title or f"Sample radial profile ({args.y})")
            plt.xlabel(x_label)
            plt.ylabel(args.y)
            if args.logy:
                plt.yscale("log")
            plt.tight_layout()
            plt.savefig(out, dpi=200)
            print(f"wrote {out}")
    else:
        print("# r_mid_m, value")
        for x, y in zip(xs, ys):
            print(f"{x:.6e}, {y:.6e}")


def _find_latest_density(search_dir: Path) -> Path:
    if not search_dir.exists():
        raise SystemExit(f"search dir not found: {search_dir}")
    best_path = None
    best_mtime = None
    for p in search_dir.rglob("ion_density_grid.dat"):
        try:
            mtime = p.stat().st_mtime
        except Exception:
            continue
        if best_mtime is None or mtime > best_mtime:
            best_mtime = mtime
            best_path = p
    if not best_path:
        raise SystemExit(f"no ion_density_grid.dat found under {search_dir}")
    return best_path


def _load_density_grid(path: Path) -> Tuple[Dict[str, Any], List[List[float]]]:
    if not path.exists():
        raise SystemExit(f"grid not found: {path}")
    meta: Dict[str, Any] = {}
    rows: List[List[float]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].strip().split()
                if len(parts) >= 2 and parts[0] in ("nx", "x0", "norm"):
                    for i in range(0, len(parts) - 1, 2):
                        key = parts[i]
                        val = parts[i + 1]
                        try:
                            meta[key] = float(val) if "." in val or "e" in val.lower() else int(val)
                        except Exception:
                            meta[key] = val
                continue
            rows.append([float(x) for x in line.split()])
    return meta, rows


def _density_plot(args: argparse.Namespace) -> None:
    if args.grid:
        grid_path = Path(args.grid)
    else:
        grid_path = _find_latest_density(Path(args.search_dir))
        print(f"using {grid_path}")

    meta_path = None
    if args.meta:
        meta_path = Path(args.meta)
    else:
        cand = grid_path.parent / "meta.json"
        if cand.exists():
            meta_path = cand

    meta, rows = _load_density_grid(grid_path)
    if not rows:
        raise SystemExit("density grid is empty")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise SystemExit(f"matplotlib required for density plotting: {exc}") from exc

    nx_out = len(rows[0])
    ny_out = len(rows)
    x0 = float(meta.get("x0", 0.0))
    y0 = float(meta.get("y0", 0.0))
    h = float(meta.get("h", 1.0))
    stride = int(meta.get("stride", 1))
    x1 = x0 + (nx_out - 1) * h * stride
    y1 = y0 + (ny_out - 1) * h * stride

    scale = 1000.0 if args.units == "mm" else 1.0
    extent = [x0 * scale, x1 * scale, y0 * scale, y1 * scale]
    xlabel = "x (mm)" if args.units == "mm" else "x (m)"
    ylabel = "y (mm)" if args.units == "mm" else "y (m)"

    def draw_electrodes(ax, m):
        try:
            import matplotlib.patches as patches
        except Exception as exc:
            raise SystemExit(f"matplotlib patches required for electrodes: {exc}") from exc

        grid_t = _as_float(m.get("GRID_T_M")) or 0.0
        gap = _as_float(m.get("GAP_M")) or 0.0
        ybox = _as_float(m.get("YBOX_M")) or 0.0
        tube = m.get("tube") or {}
        tube_x0 = _as_float(tube.get("x_start")) or 0.0
        tube_x1 = _as_float(tube.get("x_end")) or 0.0
        wall_t = _as_float(tube.get("wall_t")) or 0.0
        holder_t = _as_float(tube.get("holder_t")) or 0.0
        holder = m.get("SAMPLE_HOLDER") or {}
        holder_en = bool(holder.get("enabled", False))
        holder_h = _as_float(holder.get("height_m")) or 0.0
        if holder_t <= 0.0:
            holder_t = _as_float(holder.get("thickness_m")) or 0.0

        color = args.electrode_color
        lw = 1.0

        if grid_t > 0.0 and ybox > 0.0:
            xs0 = 0.0
            xs1 = xs0 + grid_t
            xa0 = xs1 + gap
            xa1 = xa0 + grid_t
            ax.add_patch(
                patches.Rectangle(
                    (xs0 * scale, -ybox * scale),
                    (xs1 - xs0) * scale,
                    2.0 * ybox * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )
            ax.add_patch(
                patches.Rectangle(
                    (xa0 * scale, -ybox * scale),
                    (xa1 - xa0) * scale,
                    2.0 * ybox * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )

        if tube_x1 > tube_x0 and ybox > 0.0 and wall_t > 0.0:
            ax.add_patch(
                patches.Rectangle(
                    (tube_x0 * scale, (ybox - wall_t) * scale),
                    (tube_x1 - tube_x0) * scale,
                    wall_t * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )
            ax.add_patch(
                patches.Rectangle(
                    (tube_x0 * scale, (-ybox) * scale),
                    (tube_x1 - tube_x0) * scale,
                    wall_t * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )

        if holder_en and holder_t > 0.0 and holder_h > 0.0 and tube_x1 > 0.0:
            ax.add_patch(
                patches.Rectangle(
                    ((tube_x1 - holder_t) * scale, (-0.5 * holder_h) * scale),
                    holder_t * scale,
                    holder_h * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )

    plt.figure(figsize=(6, 6))
    plt.imshow(rows, origin="lower", aspect="auto", extent=extent, cmap=args.cmap)

    if meta_path and not args.no_electrodes:
        m = _load_json(meta_path) or {}
        draw_electrodes(plt.gca(), m)
        try:
            import matplotlib.patches as patches
        except Exception as exc:
            raise SystemExit(f"matplotlib patches required for electrodes: {exc}") from exc

        m = _load_json(meta_path) or {}
        grid_t = _as_float(m.get("GRID_T_M")) or 0.0
        gap = _as_float(m.get("GAP_M")) or 0.0
        ybox = _as_float(m.get("YBOX_M")) or 0.0
        tube = m.get("tube") or {}
        tube_x0 = _as_float(tube.get("x_start")) or 0.0
        tube_x1 = _as_float(tube.get("x_end")) or 0.0
        wall_t = _as_float(tube.get("wall_t")) or 0.0
        holder_t = _as_float(tube.get("holder_t")) or 0.0
        holder = m.get("SAMPLE_HOLDER") or {}
        holder_en = bool(holder.get("enabled", False))
        holder_h = _as_float(holder.get("height_m")) or 0.0
        if holder_t <= 0.0:
            holder_t = _as_float(holder.get("thickness_m")) or 0.0

        color = args.electrode_color
        lw = 1.0

        if grid_t > 0.0 and ybox > 0.0:
            xs0 = 0.0
            xs1 = xs0 + grid_t
            xa0 = xs1 + gap
            xa1 = xa0 + grid_t
            plt.gca().add_patch(
                patches.Rectangle(
                    (xs0 * scale, -ybox * scale),
                    (xs1 - xs0) * scale,
                    2.0 * ybox * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )
            plt.gca().add_patch(
                patches.Rectangle(
                    (xa0 * scale, -ybox * scale),
                    (xa1 - xa0) * scale,
                    2.0 * ybox * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )

        if tube_x1 > tube_x0 and ybox > 0.0 and wall_t > 0.0:
            plt.gca().add_patch(
                patches.Rectangle(
                    (tube_x0 * scale, (ybox - wall_t) * scale),
                    (tube_x1 - tube_x0) * scale,
                    wall_t * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )
            plt.gca().add_patch(
                patches.Rectangle(
                    (tube_x0 * scale, (-ybox) * scale),
                    (tube_x1 - tube_x0) * scale,
                    wall_t * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )

        if holder_en and holder_t > 0.0 and holder_h > 0.0 and tube_x1 > 0.0:
            plt.gca().add_patch(
                patches.Rectangle(
                    ((tube_x1 - holder_t) * scale, (-0.5 * holder_h) * scale),
                    holder_t * scale,
                    holder_h * scale,
                    fill=False,
                    edgecolor=color,
                    linewidth=lw,
                )
            )
    plt.colorbar(label=args.zlabel)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(args.title or "Ion number density")
    plt.tight_layout()

    out = Path(args.out) if args.out else Path("ion_density.png")
    plt.savefig(out, dpi=200)
    print(f"wrote {out} (stride={stride})")

    if args.closeups and meta_path:
        m = _load_json(meta_path) or {}
        scr_list = m.get("SCREEN_OFF_LIST_M") or []
        acc_list = m.get("ACCEL_OFF_LIST_M") or []
        pairs = []
        if isinstance(scr_list, list) and isinstance(acc_list, list) and len(scr_list) == len(acc_list) and scr_list:
            for s, a in zip(scr_list, acc_list):
                try:
                    pairs.append((float(s), float(a)))
                except Exception:
                    continue
        if not pairs:
            pairs_raw = m.get("APERTURE_PAIRS_M") or ""
            toks = [t.strip() for t in pairs_raw.replace(";", ",").split(",") if t.strip()]
            for tok in toks:
                parts = [p.strip() for p in tok.split(":") if p.strip()]
                if len(parts) < 2:
                    continue
                try:
                    pairs.append((float(parts[0]), float(parts[1])))
                except Exception:
                    continue

        if pairs:
            a_scr = _as_float(m.get("SCREEN_AP_RAD_M")) or _as_float(m.get("AP_RAD_M")) or 0.0
            a_acc = _as_float(m.get("ACCEL_AP_RAD_M")) or a_scr
            scr_rads = m.get("SCREEN_RAD_LIST_M") or []
            acc_rads = m.get("ACCEL_RAD_LIST_M") or []
            try:
                scr_max = max(float(v) for v in scr_rads) if scr_rads else a_scr
            except Exception:
                scr_max = a_scr
            try:
                acc_max = max(float(v) for v in acc_rads) if acc_rads else a_acc
            except Exception:
                acc_max = a_acc
            a_max = max(scr_max, acc_max)
            scr = m.get("screen_chamfer") or {}
            acc = m.get("accel_chamfer") or {}
            def tan_deg(v):
                try:
                    return math.tan(math.radians(float(v)))
                except Exception:
                    return 0.0
            max_delta = max(
                (float(scr.get("up_depth", 0.0)) * tan_deg(scr.get("up_angle_deg", 0.0))),
                (float(scr.get("dn_depth", 0.0)) * tan_deg(scr.get("dn_angle_deg", 0.0))),
                (float(acc.get("up_depth", 0.0)) * tan_deg(acc.get("up_angle_deg", 0.0))),
                (float(acc.get("dn_depth", 0.0)) * tan_deg(acc.get("dn_angle_deg", 0.0))),
            )
            pad = args.close_pad

            grid_t = _as_float(m.get("GRID_T_M")) or 0.0
            gap = _as_float(m.get("GAP_M")) or 0.0
            xs0 = 0.0
            xs1 = xs0 + grid_t
            xa0 = xs1 + gap
            xa1 = xa0 + grid_t
            x_lo = xs0 - pad
            x_hi = xa1 + pad

            close_pairs = []
            eps = 1.0e-12
            for s, a in pairs:
                if 0.5 * (s + a) >= -eps:
                    close_pairs.append((s, a))
            if not close_pairs:
                close_pairs.append(pairs[0])

            for i, (s, a) in enumerate(close_pairs):
                y_min = min(s, a)
                y_max = max(s, a)
                y_center = 0.5 * (y_min + y_max)
                y_span = max(abs(y_min - y_center), abs(y_max - y_center))
                y_extent = y_span + a_max + max_delta + pad
                y_lo = y_center - y_extent
                y_hi = y_center + y_extent

                plt.figure(figsize=(6, 6))
                plt.imshow(rows, origin="lower", aspect="auto", extent=extent, cmap=args.cmap)
                ax = plt.gca()
                if not args.no_electrodes:
                    draw_electrodes(ax, m)
                ax.set_xlim(x_lo * scale, x_hi * scale)
                ax.set_ylim(y_lo * scale, y_hi * scale)
                plt.colorbar(label=args.zlabel)
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
                plt.title(args.title or "Ion number density")
                plt.tight_layout()

                out_close = _add_suffix_before_ext(out, "_closeup" if i == 0 else f"_closeup_p{i}")
                plt.savefig(out_close, dpi=200)
                print(f"wrote {out_close}")


def main(argv: Optional[Iterable[str]] = None) -> None:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--csv", default="results/runlog_compact.csv", help="Path to runlog_compact.csv")
    common.add_argument("--filter", default=None, help="pandas query() expression to filter rows")

    ap = argparse.ArgumentParser(description="Plot runlog_compact.csv (scatter + heatmap).", parents=[common])

    sp = ap.add_subparsers(dest="cmd", required=True)

    sp_summary = sp.add_parser("summary", help="Print column overview", parents=[common])
    sp_summary.add_argument("--top", type=int, default=20, help="Max columns to show per section")

    sp_scatter = sp.add_parser("scatter", help="Scatter plot", parents=[common])
    sp_scatter.add_argument("--x", required=True, help="x column")
    sp_scatter.add_argument("--y", required=True, help="y column")
    sp_scatter.add_argument("--color", default=None, help="color column")
    sp_scatter.add_argument("--size", default=None, help="size column")
    sp_scatter.add_argument("--hover", default="RUN_DIR,runname", help="comma-separated hover columns")
    sp_scatter.add_argument("--facet-row", default=None, help="facet row column")
    sp_scatter.add_argument("--facet-col", default=None, help="facet col column")
    sp_scatter.add_argument("--logx", action="store_true", help="log-scale x axis")
    sp_scatter.add_argument("--logy", action="store_true", help="log-scale y axis")
    sp_scatter.add_argument("--title", default=None, help="plot title")
    sp_scatter.add_argument("--out", default=None, help="output html path")
    sp_scatter.add_argument("--show", action="store_true", help="open interactive window")

    sp_heat = sp.add_parser("heatmap", help="2D heatmap", parents=[common])
    sp_heat.add_argument("--x", required=True, help="x column")
    sp_heat.add_argument("--y", required=True, help="y column")
    sp_heat.add_argument("--z", required=True, help="value column")
    sp_heat.add_argument("--agg", default="mean", choices=["mean", "median", "min", "max", "count"], help="aggregation")
    sp_heat.add_argument("--xbins", type=int, default=None, help="bin count for x (auto if omitted)")
    sp_heat.add_argument("--ybins", type=int, default=None, help="bin count for y (auto if omitted)")
    sp_heat.add_argument("--smooth", type=float, default=None, help="gaussian smoothing sigma (in bins)")
    sp_heat.add_argument("--zsmooth", default="best", choices=["best", "fast", "none"], help="plotly heatmap smoothing")
    sp_heat.add_argument("--logz", action="store_true", help="log-scale color axis")
    sp_heat.add_argument("--title", default=None, help="plot title")
    sp_heat.add_argument("--out", default=None, help="output html path")
    sp_heat.add_argument("--show", action="store_true", help="open interactive window")

    sp_steer = sp.add_parser("steer-curves", help="Grouped steering curves with loss-aware line styling", parents=[common])
    sp_steer.add_argument("--x", default="ACCEL_OFF_Y_M", help="x column")
    sp_steer.add_argument("--x-units", default="m", choices=["m", "mm"], help="x-axis units")
    sp_steer.add_argument("--x-label", default=None, help="custom x-axis label")
    sp_steer.add_argument("--y", default="STEER_ANGLE_DEG", help="y column")
    sp_steer.add_argument("--y-label", default=None, help="custom y-axis label")
    sp_steer.add_argument("--group", default="AP_RAD_M", help="group/line column (e.g. AP_RAD_M or SCREEN_AP_RAD_M)")
    sp_steer.add_argument("--group-units", default="mm", choices=["m", "mm"], help="units for legend labels when grouping by numeric radius")
    sp_steer.add_argument("--loss-col", default="GRID_LOSS_FRAC", help="loss-fraction column controlling line style")
    sp_steer.add_argument("--loss-abs-col", default="GRID_LOSS_A", help="absolute loss column shown in hover")
    sp_steer.add_argument("--loss-threshold", type=float, default=0.10, help="switch to dashed styling above this loss fraction")
    sp_steer.add_argument("--hover", default="RUN_DIR,RUN_TAG,ACC_DN_DEPTH_M,ACC_DN_ANGLE_DEG", help="comma-separated hover columns")
    sp_steer.add_argument("--abs-y", action="store_true", help="plot abs(y) instead of signed steering")
    sp_steer.add_argument("--require-no-sidewalls", action="store_true", help="drop rows where LOST_TO_SIDEWALLS is true")
    sp_steer.add_argument("--show-threshold", action="store_true", help="draw a vertical marker where each line first exceeds the loss threshold")
    sp_steer.add_argument("--debug-summary", action="store_true", help="print a summary of the filtered rows used for the plot")
    sp_steer.add_argument("--logx", action="store_true", help="log-scale x axis")
    sp_steer.add_argument("--logy", action="store_true", help="log-scale y axis")
    sp_steer.add_argument("--title", default=None, help="plot title")
    sp_steer.add_argument("--out", default=None, help="output html/image path")
    sp_steer.add_argument("--show", action="store_true", help="open interactive window")

    sp_maxsafe = sp.add_parser("max-safe-offset-curves", help="Maximum safe offset vs x-column, grouped by radius", parents=[common])
    sp_maxsafe.add_argument("--x", default="GAP_M", help="x column")
    sp_maxsafe.add_argument("--x-units", default="mm", choices=["m", "mm"], help="x-axis units")
    sp_maxsafe.add_argument("--x-label", default=None, help="custom x-axis label")
    sp_maxsafe.add_argument("--offset-col", default="ACCEL_OFF_Y_M", help="offset column to maximize under the loss threshold")
    sp_maxsafe.add_argument("--y-units", default="mm", choices=["m", "mm"], help="y-axis units")
    sp_maxsafe.add_argument("--y-label", default=None, help="custom y-axis label")
    sp_maxsafe.add_argument("--group", default="AP_RAD_M", help="group/line column (e.g. AP_RAD_M or SCREEN_AP_RAD_M)")
    sp_maxsafe.add_argument("--group-units", default="mm", choices=["m", "mm"], help="units for legend labels when grouping by numeric radius")
    sp_maxsafe.add_argument("--loss-col", default="GRID_LOSS_FRAC", help="loss-fraction column defining safe runs")
    sp_maxsafe.add_argument("--loss-threshold", type=float, default=0.10, help="maximum allowed loss fraction for a run to count as safe")
    sp_maxsafe.add_argument("--require-no-sidewalls", action="store_true", help="drop rows where LOST_TO_SIDEWALLS is true")
    sp_maxsafe.add_argument("--logx", action="store_true", help="log-scale x axis")
    sp_maxsafe.add_argument("--logy", action="store_true", help="log-scale y axis")
    sp_maxsafe.add_argument("--title", default=None, help="plot title")
    sp_maxsafe.add_argument("--out", default=None, help="output html/image path")
    sp_maxsafe.add_argument("--show", action="store_true", help="open interactive window")

    sp_prof = sp.add_parser("sample-profile", help="Plot sample_diameter_profile.json")
    sp_prof.add_argument("--json", default=None, help="path to sample profile json")
    sp_prof.add_argument("--search-dir", default="results", help="search dir if --json omitted")
    sp_prof.add_argument("--kind", default="auto", choices=["auto", "diameter"], help="profile kind to auto-find")
    sp_prof.add_argument("--y", default="I_Apm", choices=["I_Apm", "I_A", "fraction", "count"], help="value field")
    sp_prof.add_argument("--csv", default=None, help="optional csv output")
    sp_prof.add_argument("--title", default=None, help="plot title")
    sp_prof.add_argument("--out", default=None, help="output html path")
    sp_prof.add_argument("--logy", action="store_true", help="log-scale y axis")
    sp_prof.add_argument("--line", action="store_true", help="line+markers instead of histogram")
    sp_prof.add_argument("--mirror", action="store_true", help="mirror radial bins to negative r")
    sp_prof.add_argument("--x-units", default="m", choices=["m", "mm"], help="x-axis units")
    sp_prof.add_argument("--show", action="store_true", help="open interactive window")

    sp_den = sp.add_parser("density-plot", help="Plot ion_density_grid.dat")
    sp_den.add_argument("--grid", default=None, help="path to ion_density_grid.dat")
    sp_den.add_argument("--meta", default=None, help="path to meta.json (optional)")
    sp_den.add_argument("--search-dir", default="results", help="search dir if --grid omitted")
    sp_den.add_argument("--out", default=None, help="output image path")
    sp_den.add_argument("--title", default=None, help="plot title")
    sp_den.add_argument("--units", default="m", choices=["m", "mm"], help="axis units")
    sp_den.add_argument("--cmap", default="viridis", help="matplotlib colormap")
    sp_den.add_argument("--zlabel", default="n_i (a.u.)", help="colorbar label")
    sp_den.add_argument("--electrode-color", default="black", help="outline color for electrodes")
    sp_den.add_argument("--no-electrodes", action="store_true", help="disable electrode overlay")
    sp_den.add_argument("--closeups", action="store_true", help="also write closeup density images")
    sp_den.add_argument("--close-pad", type=float, default=5e-4, help="closeup padding in meters")

    args = ap.parse_args(argv)
    if args.cmd in ("summary", "scatter", "heatmap", "steer-curves", "max-safe-offset-curves"):
        df = _load_csv(Path(args.csv))
        df = _apply_filter(df, args.filter)

    if args.cmd == "summary":
        _summary(df, args.top)
    elif args.cmd == "scatter":
        _scatter(df, args)
    elif args.cmd == "heatmap":
        _heatmap(df, args)
    elif args.cmd == "steer-curves":
        _steer_curves(df, args)
    elif args.cmd == "max-safe-offset-curves":
        _max_safe_offset_curves(df, args)
    elif args.cmd == "sample-profile":
        _sample_profile(args)
    elif args.cmd == "density-plot":
        _density_plot(args)
    else:
        raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    main()
