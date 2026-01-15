#!/usr/bin/env python3
"""
plot_runlog_compact.py

Interactive plots for runlog_compact.csv (scatter + heatmap).
Uses pandas + plotly to generate self-contained HTML files.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import plotly.express as px


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


def _split_list(text: Optional[str]) -> list[str]:
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

    fig = px.imshow(
        pivot,
        labels={"x": args.x, "y": args.y, "color": args.z},
        aspect="auto",
        origin="lower",
        title=args.title,
    )
    if args.logz:
        fig.update_coloraxes(type="log")
    fig.update_layout(template="plotly_white")

    out = Path(args.out) if args.out else Path(f"heatmap_{_safe_name(args.z)}_by_{_safe_name(args.x)}_{_safe_name(args.y)}.html")
    fig.write_html(out)
    print(f"wrote {out}")
    if args.show:
        fig.show()


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
    sp_heat.add_argument("--logz", action="store_true", help="log-scale color axis")
    sp_heat.add_argument("--title", default=None, help="plot title")
    sp_heat.add_argument("--out", default=None, help="output html path")
    sp_heat.add_argument("--show", action="store_true", help="open interactive window")

    args = ap.parse_args(argv)
    df = _load_csv(Path(args.csv))
    df = _apply_filter(df, args.filter)

    if args.cmd == "summary":
        _summary(df, args.top)
    elif args.cmd == "scatter":
        _scatter(df, args)
    elif args.cmd == "heatmap":
        _heatmap(df, args)
    else:
        raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    main()
