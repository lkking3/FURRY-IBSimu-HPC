#!/usr/bin/env python3
"""
beamlet_layout_solver.py

Solves beamlet coordinates for the current two-grid workflow using:
- steering relationships extracted from runlog_compact.csv
- a chosen aperture radius and gap
- an acceptable current-loss threshold
- configurable packing methodology for screen coordinates
- target modes for either ideal center convergence or packed target references
- interactive HTML plots for screen, accelerator, and overlapped layouts

Packing modes:
- linear: legacy 1D row in y
- radial: concentric rings with angular symmetry
- hcp: hexagonal close-packed lattice clipped to a disk

Target modes:
- center: ideal objective, all beamlets aim for the same target center; target radius is only an acceptance region
- packed: target coordinates are packed inside the target radius

The solver emits a 2D layout JSON plus a legacy APERTURE_PAIRS_M projection for the
current 2D cross-section solver. The legacy projection uses the y-coordinate pair list.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    _HAS_PLOTLY = True
except Exception:
    go = None
    _HAS_PLOTLY = False


AP_PACK = [
    (1.0e-3, 2221.0),
    (1.5e-3, 1063.0),
    (2.0e-3,  769.0),
    (3.0e-3,  379.0),
    (4.0e-3,  223.0),
    (5.0e-3,  151.0),
    (6.0e-3,  109.0),
    (7.0e-3,   73.0),
    (8.0e-3,   55.0),
]


def estimate_beamlet_count(ap_rad_m: float) -> float:
    if ap_rad_m is None or ap_rad_m <= 0.0:
        return AP_PACK[0][1]
    if ap_rad_m <= AP_PACK[0][0]:
        return AP_PACK[0][1]
    if ap_rad_m >= AP_PACK[-1][0]:
        return AP_PACK[-1][1]
    for (r0, n0), (r1, n1) in zip(AP_PACK[:-1], AP_PACK[1:]):
        if r0 <= ap_rad_m <= r1:
            t = (ap_rad_m - r0) / (r1 - r0) if r1 > r0 else 0.0
            return n0 + t * (n1 - n0)
    return AP_PACK[-1][1]


@dataclass
class BeamletAssignment:
    index: int
    screen_x_m: float
    screen_y_m: float
    accel_x_m: float
    accel_y_m: float
    accel_offset_x_m: float
    accel_offset_y_m: float
    accel_offset_mag_m: float
    desired_target_x_m: float
    desired_target_y_m: float
    target_x_m: float
    target_y_m: float
    steer_angle_req_deg: float
    steer_angle_used_deg: float
    clipped_to_safe_limit: bool
    hits_target: bool
    target_constraint_feasible: bool
    target_shortfall_m: float
    landing_radius_from_target_center_m: float


@dataclass
class LayoutSummary:
    aperture_radius_m: float
    gap_m: float
    max_loss_frac: float
    packing_mode: str
    target_mode: str
    screen_pitch_m: float
    plate_radius_m: float
    target_radius_m: float
    target_plane_x_m: float
    steer_origin_x_m: float
    screen_center_x_m: float
    screen_center_y_m: float
    target_center_x_m: float
    target_center_y_m: float
    beamlet_count: int
    beamlets_hitting_target: int
    target_constraint_satisfied: bool
    infeasible_beamlet_count: int
    clipped_beamlet_count: int
    max_target_shortfall_m: float
    estimated_full_plate_capacity: float
    safe_curve_points: int
    safe_offset_max_m: float
    safe_steer_max_deg: float
    aperture_pairs_m: str


@dataclass
class LayoutSolution:
    summary: LayoutSummary
    beamlets: List[BeamletAssignment]
    safe_curve: List[Tuple[float, float]]


def _coerce_numeric_like(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        coerced = pd.to_numeric(out[col], errors="coerce")
        if coerced.notna().any():
            out[col] = coerced
    return out


def _load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"csv not found: {path}")
    return _coerce_numeric_like(pd.read_csv(path))


def _apply_filter(df: pd.DataFrame, expr: Optional[str]) -> pd.DataFrame:
    if not expr:
        return df
    try:
        return df.query(expr, engine="python")
    except Exception as exc:
        raise SystemExit(f"bad --filter expression: {exc}") from exc


def _fit_count_for_span(span_m: float, pitch_m: float) -> int:
    if span_m < 0.0:
        raise ValueError("target span must be non-negative")
    if pitch_m <= 0.0:
        raise ValueError("pitch must be positive")
    return max(1, int(math.floor(span_m / pitch_m)) + 1)


def _sort_points(points: Sequence[Tuple[float, float]]) -> List[Tuple[float, float]]:
    return sorted(points, key=lambda p: (round(math.hypot(p[0], p[1]), 12), round(math.atan2(p[1], p[0]), 12), p[0], p[1]))


def _generate_linear_points(count: int, pitch_m: float, center_x_m: float = 0.0, center_y_m: float = 0.0) -> List[Tuple[float, float]]:
    if count < 1:
        return []
    mid = 0.5 * (count - 1)
    return [(center_x_m, center_y_m + (i - mid) * pitch_m) for i in range(count)]


def _generate_radial_points(radius_m: float, pitch_m: float) -> List[Tuple[float, float]]:
    if radius_m < 0.0 or pitch_m <= 0.0:
        return []
    pts: List[Tuple[float, float]] = [(0.0, 0.0)]
    if radius_m <= 1.0e-15:
        return pts
    max_ring = int(math.floor(radius_m / pitch_m + 1.0e-12))
    for ir in range(1, max_ring + 1):
        rr = ir * pitch_m
        ntheta = max(6, int(round(2.0 * math.pi * rr / pitch_m)))
        for it in range(ntheta):
            th = 2.0 * math.pi * float(it) / float(ntheta)
            x = rr * math.cos(th)
            y = rr * math.sin(th)
            if math.hypot(x, y) <= radius_m + 1.0e-12:
                pts.append((x, y))
    return _sort_points(pts)


def _generate_hcp_points(radius_m: float, pitch_m: float) -> List[Tuple[float, float]]:
    if radius_m < 0.0 or pitch_m <= 0.0:
        return []
    pts: List[Tuple[float, float]] = []
    dy = pitch_m * math.sqrt(3.0) / 2.0
    ny = int(math.floor(radius_m / dy)) if dy > 0.0 else 0
    for j in range(-ny, ny + 1):
        y = j * dy
        xoff = 0.5 * pitch_m if (j & 1) else 0.0
        xmax = math.sqrt(max(radius_m * radius_m - y * y, 0.0))
        nx = int(math.floor((xmax - xoff) / pitch_m)) if pitch_m > 0.0 else 0
        for i in range(-nx - 1, nx + 2):
            x = i * pitch_m + xoff
            if math.hypot(x, y) <= radius_m + 1.0e-12:
                pts.append((x, y))
    return _sort_points(pts)


def _generate_points(mode: str, count: Optional[int], pitch_m: float, radius_m: float) -> List[Tuple[float, float]]:
    mode = mode.strip().lower()
    if mode == "linear":
        if count is None:
            count = _fit_count_for_span(2.0 * radius_m, pitch_m)
        return _generate_linear_points(count, pitch_m)
    if mode == "radial":
        pts = _generate_radial_points(radius_m, pitch_m)
    elif mode == "hcp":
        pts = _generate_hcp_points(radius_m, pitch_m)
    else:
        raise SystemExit(f"unknown packing mode: {mode}")
    if count is not None:
        if count < 1:
            raise SystemExit("count must be at least 1")
        if count > len(pts):
            raise SystemExit(f"requested count={count} exceeds available packed positions={len(pts)}")
        return pts[:count]
    return pts


def _generate_target_points(mode: str, count: int, radius_m: float) -> List[Tuple[float, float]]:
    mode = mode.strip().lower()
    if count < 1:
        return []
    if radius_m <= 1.0e-15:
        return [(0.0, 0.0) for _ in range(count)]
    if mode == "linear":
        if count == 1:
            return [(0.0, 0.0)]
        pitch = (2.0 * radius_m) / float(count - 1)
        return _generate_linear_points(count, pitch)

    gen = _generate_radial_points if mode == "radial" else _generate_hcp_points if mode == "hcp" else None
    if gen is None:
        raise SystemExit(f"unknown packing mode: {mode}")

    hi = max(2.0 * radius_m, 1.0e-12)
    lo = hi
    pts_lo = gen(radius_m, lo)
    while len(pts_lo) < count:
        lo *= 0.5
        if lo < 1.0e-12:
            break
        pts_lo = gen(radius_m, lo)
    if len(pts_lo) < count:
        return pts_lo + [(0.0, 0.0)] * (count - len(pts_lo))

    left = lo
    right = hi
    best = pts_lo
    for _ in range(50):
        mid = 0.5 * (left + right)
        pts_mid = gen(radius_m, mid)
        if len(pts_mid) >= count:
            best = pts_mid
            left = mid
        else:
            right = mid
    return best[:count]


def _translate_points(points: Sequence[Tuple[float, float]], dx: float, dy: float) -> List[Tuple[float, float]]:
    return [(x + dx, y + dy) for x, y in points]


def _build_safe_curve(
    df: pd.DataFrame,
    radius_m: float,
    gap_m: float,
    max_loss_frac: float,
    extra_filter: Optional[str],
    radius_tol: float,
    gap_tol: float,
) -> pd.DataFrame:
    work = _apply_filter(df, extra_filter)

    required = ["AP_RAD_M", "GAP_M", "ACCEL_OFF_Y_M", "STEER_ANGLE_DEG", "GRID_LOSS_FRAC"]
    missing = [c for c in required if c not in work.columns]
    if missing:
        raise SystemExit(f"missing required columns: {', '.join(missing)}")

    work = work.copy()
    for col in required:
        work[col] = pd.to_numeric(work[col], errors="coerce")

    mask = (
        work["AP_RAD_M"].sub(radius_m).abs() <= radius_tol
    ) & (
        work["GAP_M"].sub(gap_m).abs() <= gap_tol
    ) & (
        work["GRID_LOSS_FRAC"] <= max_loss_frac
    )
    if "LOST_TO_SIDEWALLS" in work.columns:
        side = work["LOST_TO_SIDEWALLS"].astype(str).str.lower().isin(["true", "1"])
        mask &= ~side

    safe = work.loc[mask, ["ACCEL_OFF_Y_M", "STEER_ANGLE_DEG"]].dropna().copy()
    if safe.empty:
        raise SystemExit("no safe steering data left after radius/gap/loss filtering")

    safe["OFFSET_ABS_M"] = safe["ACCEL_OFF_Y_M"].abs()
    safe["STEER_ABS_DEG"] = safe["STEER_ANGLE_DEG"].abs()

    # Physical constraint: offset cannot exceed the aperture radius
    safe = safe[safe["OFFSET_ABS_M"] <= radius_m].copy()

    grouped = (
        safe.groupby("OFFSET_ABS_M", dropna=False)["STEER_ABS_DEG"]
        .max()
        .reset_index()
        .sort_values("OFFSET_ABS_M", kind="mergesort")
    )
    grouped = grouped[grouped["OFFSET_ABS_M"] >= 0.0].copy()
    if grouped.empty:
        raise SystemExit("safe steering curve is empty after grouping by offset")

    grouped["STEER_ABS_DEG"] = np.maximum.accumulate(grouped["STEER_ABS_DEG"].to_numpy())
    grouped = grouped.drop_duplicates(subset=["STEER_ABS_DEG"], keep="last").reset_index(drop=True)
    return grouped


def _format_pairs(beamlets: Sequence[BeamletAssignment], radius_m: float) -> str:
    return ", ".join(
        f"{b.screen_y_m:.9g}:{b.accel_y_m:.9g}:{radius_m:.9g}" for b in beamlets
    )


def solve_layout(
    *,
    df: pd.DataFrame,
    radius_m: float,
    gap_m: float,
    max_loss_frac: float,
    packing_mode: str,
    target_mode: str,
    screen_pitch_m: float,
    plate_radius_m: float,
    target_radius_m: float,
    target_plane_x_m: float,
    steer_origin_x_m: float,
    screen_center_x_m: float,
    screen_center_y_m: float,
    target_center_x_m: float,
    target_center_y_m: float,
    count: Optional[int],
    extra_filter: Optional[str],
    radius_tol: float,
    gap_tol: float,
    clip_to_safe_limit: bool,
) -> LayoutSolution:
    if target_plane_x_m <= steer_origin_x_m:
        raise SystemExit("target plane must be downstream of the steering origin")
    if screen_pitch_m <= 0.0:
        raise SystemExit("screen pitch must be positive")
    if plate_radius_m < 0.0 or target_radius_m < 0.0:
        raise SystemExit("plate and target radii must be non-negative")

    curve = _build_safe_curve(df, radius_m, gap_m, max_loss_frac, extra_filter, radius_tol, gap_tol)
    safe_offset_max = float(curve["OFFSET_ABS_M"].max())
    safe_steer_max = float(curve["STEER_ABS_DEG"].max())

    if count is None:
        screen_pts0 = _generate_points(packing_mode, None, screen_pitch_m, plate_radius_m)
        if len(screen_pts0) < 1:
            raise SystemExit("no packed positions available for the chosen grid geometry")
    else:
        screen_pts0 = _generate_points(packing_mode, count, screen_pitch_m, plate_radius_m)

    if target_mode == "center":
        target_pts0 = [(0.0, 0.0) for _ in range(len(screen_pts0))]
    elif target_mode == "packed":
        target_pts0 = _generate_target_points(packing_mode, len(screen_pts0), target_radius_m)
    else:
        raise SystemExit(f"unknown target mode: {target_mode}")

    screen_pts = _translate_points(screen_pts0, screen_center_x_m, screen_center_y_m)
    target_pts = _translate_points(target_pts0[:len(screen_pts0)], target_center_x_m, target_center_y_m)

    xspan = target_plane_x_m - steer_origin_x_m
    curve_theta = curve["STEER_ABS_DEG"].to_numpy(dtype=float)
    curve_offset = curve["OFFSET_ABS_M"].to_numpy(dtype=float)
    safe_reach_max = math.tan(math.radians(safe_steer_max)) * xspan if safe_steer_max > 0.0 else 0.0

    beamlets: List[BeamletAssignment] = []
    hit_count = 0
    infeasible_count = 0
    clipped_count = 0
    max_shortfall = 0.0
    eps = 1.0e-12
    for i, ((sx, sy), (tx_desired, ty_desired)) in enumerate(zip(screen_pts, target_pts)):
        dx_des = tx_desired - sx
        dy_des = ty_desired - sy
        transverse_des = math.hypot(dx_des, dy_des)
        theta_req = math.degrees(math.atan2(transverse_des, xspan)) if transverse_des > 0.0 else 0.0

        if target_mode == "center":
            if transverse_des > eps and safe_reach_max > 0.0:
                reach_used = min(transverse_des, safe_reach_max)
                ux_des = dx_des / transverse_des
                uy_des = dy_des / transverse_des
                tx_used = sx + ux_des * reach_used
                ty_used = sy + uy_des * reach_used
            else:
                tx_used = tx_desired
                ty_used = ty_desired
        else:
            if transverse_des > safe_reach_max + eps and transverse_des > eps and safe_reach_max > 0.0:
                scale = safe_reach_max / transverse_des
                tx_used = sx + dx_des * scale
                ty_used = sy + dy_des * scale
            elif transverse_des > eps or safe_reach_max == 0.0:
                tx_used = tx_desired
                ty_used = ty_desired
            else:
                tx_used = tx_desired
                ty_used = ty_desired

        actual_dx = tx_used - sx
        actual_dy = ty_used - sy
        transverse_used = math.hypot(actual_dx, actual_dy)
        theta_used = math.degrees(math.atan2(transverse_used, xspan)) if transverse_used > 0.0 else 0.0
        clipped = transverse_used + eps < transverse_des or theta_used + eps < theta_req
        if theta_req > safe_steer_max + 1.0e-9 and not clip_to_safe_limit:
            raise SystemExit(
                f"beamlet {i} requires {theta_req:.6g} deg but safe limit is only {safe_steer_max:.6g} deg"
            )
        if clipped:
            clipped_count += 1

        off_mag = float(np.interp(theta_used, curve_theta, curve_offset)) if theta_used > 0.0 else 0.0
        off_mag = min(off_mag, radius_m)  # offset cannot exceed aperture radius
        if transverse_used > eps:
            ux = actual_dx / transverse_used
            uy = actual_dy / transverse_used
        else:
            ux = 0.0
            uy = 0.0
        # Electrostatic steering is counter-steered in geometry: move the accelerator aperture
        # outward so the closer inner edge bends the beam inward toward the target.
        ax = sx - ux * off_mag
        ay = sy - uy * off_mag

        landing_radius = math.hypot(tx_used - target_center_x_m, ty_used - target_center_y_m)
        target_shortfall = max(landing_radius - target_radius_m, 0.0)
        hits_target = target_shortfall <= 1.0e-12
        target_feasible = hits_target
        if hits_target:
            hit_count += 1
        else:
            infeasible_count += 1
            max_shortfall = max(max_shortfall, target_shortfall)

        beamlets.append(
            BeamletAssignment(
                index=i,
                screen_x_m=sx,
                screen_y_m=sy,
                accel_x_m=ax,
                accel_y_m=ay,
                accel_offset_x_m=ax - sx,
                accel_offset_y_m=ay - sy,
                accel_offset_mag_m=off_mag,
                desired_target_x_m=tx_desired,
                desired_target_y_m=ty_desired,
                target_x_m=tx_used,
                target_y_m=ty_used,
                steer_angle_req_deg=theta_req,
                steer_angle_used_deg=theta_used,
                clipped_to_safe_limit=clipped,
                hits_target=hits_target,
                target_constraint_feasible=target_feasible,
                target_shortfall_m=target_shortfall,
                landing_radius_from_target_center_m=landing_radius,
            )
        )

    pairs = _format_pairs(beamlets, radius_m)
    summary = LayoutSummary(
        aperture_radius_m=radius_m,
        gap_m=gap_m,
        max_loss_frac=max_loss_frac,
        packing_mode=packing_mode,
        target_mode=target_mode,
        screen_pitch_m=screen_pitch_m,
        plate_radius_m=plate_radius_m,
        target_radius_m=target_radius_m,
        target_plane_x_m=target_plane_x_m,
        steer_origin_x_m=steer_origin_x_m,
        screen_center_x_m=screen_center_x_m,
        screen_center_y_m=screen_center_y_m,
        target_center_x_m=target_center_x_m,
        target_center_y_m=target_center_y_m,
        beamlet_count=len(beamlets),
        beamlets_hitting_target=hit_count,
        target_constraint_satisfied=(infeasible_count == 0),
        infeasible_beamlet_count=infeasible_count,
        clipped_beamlet_count=clipped_count,
        max_target_shortfall_m=max_shortfall,
        estimated_full_plate_capacity=estimate_beamlet_count(radius_m),
        safe_curve_points=len(curve),
        safe_offset_max_m=safe_offset_max,
        safe_steer_max_deg=safe_steer_max,
        aperture_pairs_m=pairs,
    )
    return LayoutSolution(
        summary=summary,
        beamlets=beamlets,
        safe_curve=[(float(r["OFFSET_ABS_M"]), float(r["STEER_ABS_DEG"])) for _, r in curve.iterrows()],
    )


def _presentation_fonts(fig: "go.Figure") -> None:
    fig.update_layout(
        template="plotly_white",
        font=dict(size=18),
        title_font=dict(size=28),
        legend=dict(font=dict(size=16), title_font=dict(size=18)),
    )
    fig.update_xaxes(title_font=dict(size=22), tickfont=dict(size=16), scaleanchor="y", scaleratio=1.0)
    fig.update_yaxes(title_font=dict(size=22), tickfont=dict(size=16))


def _hover_payload(b: BeamletAssignment) -> Dict[str, float]:
    return {
        "screen_x_mm": b.screen_x_m * 1.0e3,
        "screen_y_mm": b.screen_y_m * 1.0e3,
        "accel_x_mm": b.accel_x_m * 1.0e3,
        "accel_y_mm": b.accel_y_m * 1.0e3,
        "desired_target_x_mm": b.desired_target_x_m * 1.0e3,
        "desired_target_y_mm": b.desired_target_y_m * 1.0e3,
        "target_x_mm": b.target_x_m * 1.0e3,
        "target_y_mm": b.target_y_m * 1.0e3,
        "off_x_mm": b.accel_offset_x_m * 1.0e3,
        "off_y_mm": b.accel_offset_y_m * 1.0e3,
        "off_mag_mm": b.accel_offset_mag_m * 1.0e3,
        "theta_req_deg": b.steer_angle_req_deg,
        "theta_used_deg": b.steer_angle_used_deg,
        "target_shortfall_mm": b.target_shortfall_m * 1.0e3,
        "landing_radius_mm": b.landing_radius_from_target_center_m * 1.0e3,
    }


def _trace_for_points(beamlets: Sequence[BeamletAssignment], which: str, name: str, color: str, marker_symbol: str = "circle") -> "go.Scatter":
    xs: List[float] = []
    ys: List[float] = []
    custom: List[List[Any]] = []
    for b in beamlets:
        hp = _hover_payload(b)
        if which == "screen":
            xs.append(hp["screen_x_mm"])
            ys.append(hp["screen_y_mm"])
        elif which == "accel":
            xs.append(hp["accel_x_mm"])
            ys.append(hp["accel_y_mm"])
        else:
            raise ValueError(which)
        custom.append([
            b.index,
            hp["screen_x_mm"], hp["screen_y_mm"],
            hp["accel_x_mm"], hp["accel_y_mm"],
            hp["desired_target_x_mm"], hp["desired_target_y_mm"],
            hp["target_x_mm"], hp["target_y_mm"],
            hp["off_x_mm"], hp["off_y_mm"], hp["off_mag_mm"],
            hp["theta_req_deg"], hp["theta_used_deg"],
            hp["target_shortfall_mm"], hp["landing_radius_mm"],
            1 if b.clipped_to_safe_limit else 0,
            1 if b.hits_target else 0,
            1 if b.target_constraint_feasible else 0,
        ])
    return go.Scatter(
        x=xs,
        y=ys,
        mode="markers+text",
        text=[str(b.index) for b in beamlets],
        textposition="top center",
        name=name,
        marker=dict(size=12, color=color, symbol=marker_symbol, line=dict(width=1, color="black")),
        customdata=np.array(custom, dtype=object),
        hovertemplate=(
            "index=%{customdata[0]}<br>"
            "screen=(%{customdata[1]:.4g}, %{customdata[2]:.4g}) mm<br>"
            "accel=(%{customdata[3]:.4g}, %{customdata[4]:.4g}) mm<br>"
            "desired_target=(%{customdata[5]:.4g}, %{customdata[6]:.4g}) mm<br>"
            "landing=(%{customdata[7]:.4g}, %{customdata[8]:.4g}) mm<br>"
            "accel offset=(%{customdata[9]:.4g}, %{customdata[10]:.4g}) mm<br>"
            "|offset|=%{customdata[11]:.4g} mm<br>"
            "theta req=%{customdata[12]:.4g} deg<br>"
            "theta used=%{customdata[13]:.4g} deg<br>"
            "target shortfall=%{customdata[14]:.4g} mm<br>"
            "landing radius=%{customdata[15]:.4g} mm<br>"
            "clipped=%{customdata[16]}<br>"
            "hits_target=%{customdata[17]}<br>"
            "target_feasible=%{customdata[18]}<extra></extra>"
        ),
    )


def _axis_bounds(values: Iterable[float], pad: float) -> Tuple[float, float]:
    vals = list(values)
    if not vals:
        return (-1.0, 1.0)
    vmin = min(vals)
    vmax = max(vals)
    if abs(vmax - vmin) < 1.0e-12:
        return (vmin - pad, vmax + pad)
    return (vmin - pad, vmax + pad)


def write_layout_plots(sol: LayoutSolution, out_prefix: Path) -> List[Path]:
    if not _HAS_PLOTLY:
        raise SystemExit("plotly is required for layout plotting")
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    beamlets = sol.beamlets
    allx = [b.screen_x_m * 1.0e3 for b in beamlets] + [b.accel_x_m * 1.0e3 for b in beamlets]
    ally = [b.screen_y_m * 1.0e3 for b in beamlets] + [b.accel_y_m * 1.0e3 for b in beamlets]
    pad = max(sol.summary.aperture_radius_m * 1.0e3 * 3.0, sol.summary.screen_pitch_m * 1.0e3 * 0.6, 1.0)
    xr = _axis_bounds(allx, pad)
    yr = _axis_bounds(ally, pad)

    outputs: List[Path] = []
    specs = [
        ("screen", "Screen Grid Layout", [_trace_for_points(beamlets, "screen", "screen", "royalblue")]),
        ("accel", "Accelerator Grid Layout", [_trace_for_points(beamlets, "accel", "accelerator", "crimson", marker_symbol="diamond")]),
    ]

    overlap_fig = go.Figure()
    overlap_fig.add_trace(_trace_for_points(beamlets, "screen", "screen", "royalblue"))
    overlap_fig.add_trace(_trace_for_points(beamlets, "accel", "accelerator", "crimson", marker_symbol="diamond"))
    for b in beamlets:
        overlap_fig.add_trace(go.Scatter(
            x=[b.screen_x_m * 1.0e3, b.accel_x_m * 1.0e3],
            y=[b.screen_y_m * 1.0e3, b.accel_y_m * 1.0e3],
            mode="lines",
            line=dict(color="rgba(80,80,80,0.4)", width=1),
            hoverinfo="skip",
            showlegend=False,
        ))
    overlap_fig.update_layout(title="Overlapped Screen and Accelerator Layout")
    overlap_fig.update_xaxes(title_text="x (mm)", range=list(xr))
    overlap_fig.update_yaxes(title_text="y (mm)", range=list(yr))
    _presentation_fonts(overlap_fig)
    overlap_path = out_prefix.with_name(out_prefix.name + "_overlap.html")
    overlap_fig.write_html(overlap_path)
    outputs.append(overlap_path)

    for key, title, traces in specs:
        fig = go.Figure()
        for tr in traces:
            fig.add_trace(tr)
        fig.update_layout(title=title)
        fig.update_xaxes(title_text="x (mm)", range=list(xr))
        fig.update_yaxes(title_text="y (mm)", range=list(yr))
        _presentation_fonts(fig)
        out = out_prefix.with_name(out_prefix.name + f"_{key}.html")
        fig.write_html(out)
        outputs.append(out)

    side_fig = go.Figure()
    gap_mm = sol.summary.gap_m * 1.0e3
    target_plane_mm = sol.summary.target_plane_x_m * 1.0e3
    target_radius_mm = sol.summary.target_radius_m * 1.0e3
    target_cx = sol.summary.target_center_x_m * 1.0e3
    target_cy = sol.summary.target_center_y_m * 1.0e3
    target_centered = [
        math.hypot(b.screen_x_m * 1.0e3 - target_cx, b.screen_y_m * 1.0e3 - target_cy)
        for b in beamlets
    ] + [
        math.hypot(b.accel_x_m * 1.0e3 - target_cx, b.accel_y_m * 1.0e3 - target_cy)
        for b in beamlets
    ] + [
        math.hypot(b.target_x_m * 1.0e3 - target_cx, b.target_y_m * 1.0e3 - target_cy)
        for b in beamlets
    ]
    side_ymax = max(target_centered + [target_radius_mm, 1.0])
    side_pad = max(side_ymax * 0.08, 1.0)

    for i, b in enumerate(beamlets):
        screen_r = math.hypot(b.screen_x_m * 1.0e3 - target_cx, b.screen_y_m * 1.0e3 - target_cy)
        accel_r = math.hypot(b.accel_x_m * 1.0e3 - target_cx, b.accel_y_m * 1.0e3 - target_cy)
        target_r = math.hypot(b.target_x_m * 1.0e3 - target_cx, b.target_y_m * 1.0e3 - target_cy)
        desired_r = math.hypot(b.desired_target_x_m * 1.0e3 - target_cx, b.desired_target_y_m * 1.0e3 - target_cy)
        feasible = b.target_constraint_feasible
        color = "#2e8b57" if feasible else "#c0392b"
        side_fig.add_trace(go.Scatter(
            x=[0.0, gap_mm, target_plane_mm],
            y=[screen_r, accel_r, target_r],
            mode="lines+markers+text",
            text=["", "", str(b.index)],
            textposition="top center",
            name="within target" if feasible else "outside target",
            legendgroup="within target" if feasible else "outside target",
            showlegend=(i == 0 and feasible) or (i == 0 and not feasible) or (feasible and not any(bb.target_constraint_feasible for bb in beamlets[:i])) or ((not feasible) and not any((not bb.target_constraint_feasible) for bb in beamlets[:i])),
            line=dict(color=color, width=2),
            marker=dict(size=[8, 8, 10], color=["royalblue", "crimson", color], symbol=["circle", "diamond", "circle"]),
            customdata=np.array([[b.index, screen_r, accel_r, target_r, desired_r, b.steer_angle_used_deg, b.target_shortfall_m * 1.0e3, 1 if feasible else 0]], dtype=object),
            hovertemplate=(
                "index=%{customdata[0]}<br>"
                "screen radius=%{customdata[1]:.4g} mm<br>"
                "accel radius=%{customdata[2]:.4g} mm<br>"
                "landing radius=%{customdata[3]:.4g} mm<br>"
                "requested target radius=%{customdata[4]:.4g} mm<br>"
                "theta used=%{customdata[5]:.4g} deg<br>"
                "target shortfall=%{customdata[6]:.4g} mm<br>"
                "target_feasible=%{customdata[7]}<extra></extra>"
            ),
        ))
        if abs(desired_r - target_r) > 1.0e-9:
            side_fig.add_trace(go.Scatter(
                x=[target_plane_mm],
                y=[desired_r],
                mode="markers",
                marker=dict(size=9, color="#cc8f00", symbol="circle-open"),
                hovertemplate=(
                    f"index={b.index}<br>requested target radius={desired_r:.4g} mm<extra></extra>"
                ),
                showlegend=False,
            ))

    for xpos, label, color in [(0.0, "screen", "#1f3fb5"), (gap_mm, "accelerator", "#b53020"), (target_plane_mm, "target", "#1f6f45")]:
        side_fig.add_vline(x=xpos, line_color=color, line_dash="dot", line_width=1.5)
        side_fig.add_annotation(x=xpos, y=side_ymax + side_pad * 0.25, text=label, showarrow=False, font=dict(size=14, color=color))
    side_fig.add_hline(y=target_radius_mm, line_color="#2e8b57", line_dash="dash", line_width=2)
    side_fig.add_annotation(x=target_plane_mm, y=target_radius_mm, text="acceptable target radius", showarrow=False, xanchor="right", yanchor="bottom", font=dict(size=14, color="#2e8b57"))
    side_fig.update_layout(title="Side View of Beam Steering to Target")
    side_fig.update_xaxes(title_text="Axial Distance (mm)", range=[-pad, target_plane_mm + pad])
    side_fig.update_yaxes(title_text="Radial Offset From Target Center (mm)", range=[-0.5, side_ymax + side_pad])
    _presentation_fonts(side_fig)
    side_path = out_prefix.with_name(out_prefix.name + "_side.html")
    side_fig.write_html(side_path)
    outputs.append(side_path)
    return outputs


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Solve beamlet coordinates from safe steering data")
    ap.add_argument("--csv", required=True, help="path to runlog_compact.csv")
    ap.add_argument("--radius-mm", type=float, required=True, help="aperture radius in mm")
    ap.add_argument("--gap-mm", type=float, required=True, help="grid gap in mm")
    ap.add_argument("--max-loss-frac", type=float, default=0.10, help="maximum acceptable grid loss fraction")
    ap.add_argument("--packing-mode", default="linear", choices=["linear", "radial", "hcp"], help="beamlet packing methodology")
    ap.add_argument("--target-mode", default="center", choices=["center", "packed"], help="targeting objective; center is the ideal-convergence mode")
    ap.add_argument("--screen-pitch-mm", type=float, required=True, help="center-to-center beamlet pitch in mm")
    ap.add_argument("--plate-radius-mm", type=float, required=True, help="screen plate packing radius in mm")
    ap.add_argument("--target-radius-mm", type=float, default=None, help="target packing radius in mm; defaults to plate radius")
    ap.add_argument("--target-plane-mm", type=float, required=True, help="downstream target plane x position in mm")
    ap.add_argument("--steer-origin-mm", type=float, default=0.0, help="x origin for steering estimate in mm")
    ap.add_argument("--screen-center-x-mm", type=float, default=0.0, help="screen pattern center x in mm")
    ap.add_argument("--screen-center-y-mm", type=float, default=0.0, help="screen pattern center y in mm")
    ap.add_argument("--target-center-x-mm", type=float, default=0.0, help="target pattern center x in mm")
    ap.add_argument("--target-center-y-mm", type=float, default=0.0, help="target pattern center y in mm")
    ap.add_argument("--count", type=int, default=None, help="optional explicit beamlet count; otherwise use all packed positions")
    ap.add_argument("--extra-filter", default=None, help="optional pandas query to pin geometry settings before solving")
    ap.add_argument("--radius-tol-mm", type=float, default=1.0e-6, help="radius match tolerance in mm")
    ap.add_argument("--gap-tol-mm", type=float, default=1.0e-6, help="gap match tolerance in mm")
    ap.add_argument("--clip-to-safe-limit", action="store_true", help="clip unreachable target points to the safe steering limit instead of failing")
    ap.add_argument("--out-json", default=None, help="optional JSON output path")
    ap.add_argument("--plot-prefix", default=None, help="optional output prefix for screen/accel/overlap HTML plots")
    args = ap.parse_args(argv)

    df = _load_csv(Path(args.csv))
    target_radius_mm = args.plate_radius_mm if args.target_radius_mm is None else args.target_radius_mm
    sol = solve_layout(
        df=df,
        radius_m=args.radius_mm * 1.0e-3,
        gap_m=args.gap_mm * 1.0e-3,
        max_loss_frac=args.max_loss_frac,
        packing_mode=args.packing_mode,
        target_mode=args.target_mode,
        screen_pitch_m=args.screen_pitch_mm * 1.0e-3,
        plate_radius_m=args.plate_radius_mm * 1.0e-3,
        target_radius_m=target_radius_mm * 1.0e-3,
        target_plane_x_m=args.target_plane_mm * 1.0e-3,
        steer_origin_x_m=args.steer_origin_mm * 1.0e-3,
        screen_center_x_m=args.screen_center_x_mm * 1.0e-3,
        screen_center_y_m=args.screen_center_y_mm * 1.0e-3,
        target_center_x_m=args.target_center_x_mm * 1.0e-3,
        target_center_y_m=args.target_center_y_mm * 1.0e-3,
        count=args.count,
        extra_filter=args.extra_filter,
        radius_tol=args.radius_tol_mm * 1.0e-3,
        gap_tol=args.gap_tol_mm * 1.0e-3,
        clip_to_safe_limit=args.clip_to_safe_limit,
    )

    payload: Dict[str, Any] = {
        "summary": asdict(sol.summary),
        "safe_curve": [{"offset_m": x, "steer_deg": y} for x, y in sol.safe_curve],
        "beamlets": [asdict(b) for b in sol.beamlets],
    }

    print(json.dumps(payload["summary"], indent=2))
    print()
    print("APERTURE_PAIRS_M=")
    print(sol.summary.aperture_pairs_m)

    if args.out_json:
        out = Path(args.out_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {out}")

    if args.plot_prefix:
        outs = write_layout_plots(sol, Path(args.plot_prefix))
        for p in outs:
            print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
