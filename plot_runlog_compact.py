#!/usr/bin/env python3
"""
plot_runlog_compact.py

Interactive plots for runlog_compact.csv (scatter + heatmap).
Uses pandas + plotly to generate self-contained HTML files.
"""
import argparse
import json
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

    out = Path(args.out) if args.out else Path(f"scatter_{_safe_name(args.x)}_vs_{_safe_name(args.y)}.html")
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
        xaxis_title=args.x,
        yaxis_title=args.y,
        template="plotly_white",
    )
    if args.logz:
        fig.update_coloraxes(type="log")

    out = Path(args.out) if args.out else Path(f"heatmap_{_safe_name(args.z)}_by_{_safe_name(args.x)}_{_safe_name(args.y)}.html")
    fig.write_html(out)
    print(f"wrote {out}")
    if args.show:
        fig.show()


def _find_latest_profile(search_dir: Path, kind: str) -> Path:
    if not search_dir.exists():
        raise SystemExit(f"search dir not found: {search_dir}")

    patterns = []
    if kind == "radial":
        patterns = ["sample_radial_profile.json"]
    elif kind == "diameter":
        patterns = ["sample_diameter_profile.json"]
    else:
        patterns = ["sample_diameter_profile.json", "sample_radial_profile.json"]

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
    for b in bins:
        if mode == "diameter":
            x_lo = float(b.get("y_lo_m", 0.0))
            x_hi = float(b.get("y_hi_m", 0.0))
        else:
            x_lo = float(b.get("r_lo_m", 0.0))
            x_hi = float(b.get("r_hi_m", 0.0))
        rows.append(
            {
                "x_lo_m": x_lo,
                "x_hi_m": x_hi,
                "x_mid_m": 0.5 * (x_lo + x_hi),
                "I_Apm": float(b.get("I_Apm", 0.0)),
                "I_A": float(b.get("I_A", 0.0)),
                "fraction": float(b.get("fraction", 0.0)),
                "count": int(b.get("count", 0)),
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
        if _HAS_PLOTLY:
            if args.line:
                fig = go.Figure(data=go.Scatter(x=xs, y=ys, mode="lines+markers"))
            else:
                fig = go.Figure(data=go.Bar(x=xs, y=ys))

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
                    plt.bar(xs, ys, width=bar_w, align="center")
                else:
                    plt.bar(xs, ys)
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

    sp_prof = sp.add_parser("sample-profile", help="Plot sample_radial_profile.json")
    sp_prof.add_argument("--json", default=None, help="path to sample profile json")
    sp_prof.add_argument("--search-dir", default="results", help="search dir if --json omitted")
    sp_prof.add_argument("--kind", default="auto", choices=["auto", "radial", "diameter"], help="profile kind to auto-find")
    sp_prof.add_argument("--y", default="I_Apm", choices=["I_Apm", "I_A", "fraction", "count"], help="value field")
    sp_prof.add_argument("--csv", default=None, help="optional csv output")
    sp_prof.add_argument("--title", default=None, help="plot title")
    sp_prof.add_argument("--out", default=None, help="output html path")
    sp_prof.add_argument("--logy", action="store_true", help="log-scale y axis")
    sp_prof.add_argument("--line", action="store_true", help="line+markers instead of histogram")
    sp_prof.add_argument("--mirror", action="store_true", help="mirror radial bins to negative r")
    sp_prof.add_argument("--x-units", default="m", choices=["m", "mm"], help="x-axis units")
    sp_prof.add_argument("--show", action="store_true", help="open interactive window")

    args = ap.parse_args(argv)
    if args.cmd in ("summary", "scatter", "heatmap"):
        df = _load_csv(Path(args.csv))
        df = _apply_filter(df, args.filter)

    if args.cmd == "summary":
        _summary(df, args.top)
    elif args.cmd == "scatter":
        _scatter(df, args)
    elif args.cmd == "heatmap":
        _heatmap(df, args)
    elif args.cmd == "sample-profile":
        _sample_profile(args)
    else:
        raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    main()
