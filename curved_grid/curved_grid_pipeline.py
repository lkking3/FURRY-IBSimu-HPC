"""
curved_grid_pipeline.py

Extension of multi_grid_pipeline that adds spherically curved grid support.

Physics background
------------------
On a spherical grid of radius R (concave-upstream, vertex pointing toward the ion
source), an aperture at radial offset r from the grid axis:

  * sits at axial displacement  Δz(r) = R(1 − cos(arcsin(r/R)))  from the vertex
  * has its bore axis tilted inward by  α(r) = arcsin(r/R)  from the beam axis

This means outer apertures are both shifted downstream and pre-angled toward the
axis, producing a distributed focusing effect across the grid without any per-
aperture steering offset.  A convex-upstream (concave_upstream=False) grid produces
the opposite diverging effect.

Relationship to multi_grid_pipeline
------------------------------------
* Aperture, Chamfer, GridDefinition, serialize_grid_stack, and _format_num are
  re-exported from multi_grid_pipeline for convenience.
* CurvedGridDefinition mirrors GridDefinition but adds a curvature field.
* serialize_curved_grid_stack() appends a 9th pipe-delimited field
  (curv_radius:concave_flag:align_flag) to each grid token in GRID_STACK.  The companion
  multi_grid_2d binary must be built with curved-grid support to use this field;
  older builds that enforce ≤8 fields will reject it.
* flatten_curved_grid() strips the curvature field, allowing the output to be
  passed to the unchanged binary for geometry-only runs.

C++ binary requirements
-----------------------
multi_grid_2d.cpp must be modified (see plan) to:
  1. Accept ≤9 pipe fields per grid token.
  2. Parse field[8] as curv_radius:concave_flag via parse_curvature_token().
  3. Store curv_radius / curv_concave_up in GridSpec.
  4. Use curved_surface_x() in grid_fn_idx() for the slab test.
  5. Use in_any_aperture_curved() for tilted aperture checking.
  6. Use per-aperture tilted emission lines in the Bohm injection loop.
"""
from __future__ import annotations

import datetime
import math
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Locate sim/ relative to this file so grid_types is importable regardless
# of the working directory or how this script is invoked.
_HERE = Path(__file__).resolve().parent
_SIM_DIR = _HERE.parent / "sim"
if _SIM_DIR.exists() and str(_SIM_DIR) not in sys.path:
    sys.path.insert(0, str(_SIM_DIR))

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

from grid_types import (
    Aperture,
    Chamfer,
    GridDefinition,
    SimulationCase,
    _format_num,
    serialize_grid_stack,
)

try:
    from sweep_types import SweepAxis, SweepSpec, linspace, arange
    _HAS_SWEEP = True
except ImportError:
    _HAS_SWEEP = False


# ---------------------------------------------------------------------------
# Convenience: beamlet pair builder
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BeamletPair:
    """One matched screen+accel aperture pair for a single ion beamlet.

    Parameters
    ----------
    screen_offset_m:
        Radial position of the screen grid aperture centre from the beam axis [m].
    steering_arc_m:
        Arc-length displacement of the accel aperture from the concentric-sphere
        (zero-steering) position, measured along the accel sphere surface [m].

        Zero means the radial line from the shared focal point (center of
        curvature) passes exactly through both aperture centres — perfect bore
        alignment with no net steering.  Positive values shift the accel aperture
        outward (toward larger |r|); negative values shift it inward.

        The absolute accel radial offset is computed from this via
        :func:`accel_offset_from_steering` using the R values passed to
        :func:`make_apertures_from_pairs`.
    screen_radius_m:
        Bore radius of the screen aperture [m].
    accel_radius_m:
        Bore radius of the accel aperture [m].  Defaults to ``screen_radius_m``.
    """
    screen_offset_m: float
    steering_arc_m: float = 0.0
    screen_radius_m: float = 0.0015
    accel_radius_m: float = -1.0   # sentinel → copies screen_radius_m

    def __post_init__(self) -> None:
        if self.accel_radius_m < 0.0:
            object.__setattr__(self, 'accel_radius_m', self.screen_radius_m)


def accel_offset_from_steering(
    r_scr_m: float,
    R_scr_m: float,
    R_acc_m: float,
    steering_arc_m: float = 0.0,
) -> float:
    """Return the absolute accel aperture radial offset for a given steering arc-length.

    The concentric-sphere zero-steering position of the accel aperture is:

        r_acc_concentric = R_acc × sin(arcsin(r_scr / R_scr))

    A steering arc-length ``s`` shifts the aperture by Δα = s / R_acc along the
    accel sphere, giving:

        r_acc = R_acc × sin(arcsin(r_scr / R_scr) + s / R_acc)

    When either grid is flat (radius = inf) the arc-length is treated as a
    simple linear radial offset added to the concentric position.

    Parameters
    ----------
    r_scr_m:
        Signed radial offset of the screen aperture [m].
    R_scr_m:
        Screen radius of curvature [m]  (math.inf for flat).
    R_acc_m:
        Accel radius of curvature [m]   (math.inf for flat).
    steering_arc_m:
        Arc-length steering offset [m].  Zero → concentric alignment.

    Returns
    -------
    float
        Absolute radial offset of the accel aperture [m], with the same sign
        as *r_scr_m*.
    """
    sign = 1.0 if r_scr_m >= 0.0 else -1.0
    abs_r = abs(r_scr_m)
    if math.isinf(R_scr_m) or math.isinf(R_acc_m):
        # Flat grid(s): concentric scale + linear steering offset
        scale = 1.0 if math.isinf(R_scr_m) else R_acc_m / R_scr_m
        return sign * (abs_r * scale + steering_arc_m)
    alpha_scr = math.asin(min(1.0, abs_r / R_scr_m))
    delta_alpha = steering_arc_m / R_acc_m
    alpha_acc = max(-math.pi / 2.0, min(math.pi / 2.0, alpha_scr + delta_alpha))
    return sign * R_acc_m * math.sin(alpha_acc)


def make_apertures_from_pairs(
    pairs: List["BeamletPair"],
    *,
    mirror: bool = True,
    R_scr_m: float = math.inf,
    R_acc_m: float = math.inf,
) -> Tuple[List[Aperture], List[Aperture]]:
    """Build screen and accel aperture lists from a list of BeamletPair objects.

    Equivalent to ``aperture_pairs_m`` in the flat-grid system.

    Parameters
    ----------
    pairs:
        Beamlet definitions ordered from on-axis to outermost.  On-axis pairs
        (screen_offset_m == 0) are included once regardless of *mirror*.
    mirror:
        When True, each off-axis pair is reflected to produce a symmetric ±y
        aperture set.  Set this to False if your grid definition already
        handles mirroring via ``CurvedGridDefinition(mirror=True)``.
    R_scr_m:
        Screen radius of curvature [m].  Required to compute the absolute accel
        offset from ``BeamletPair.steering_arc_m``.  Defaults to ``math.inf``
        (flat grid).
    R_acc_m:
        Accel radius of curvature [m].  Required to convert steering arc-length
        to an angular increment on the accel sphere.  Defaults to ``math.inf``.

    Returns
    -------
    (screen_apertures, accel_apertures)
        Two lists of :class:`Aperture` objects ready to pass to
        :class:`CurvedGridDefinition`.
    """
    scr: List[Aperture] = []
    acc: List[Aperture] = []
    for p in pairs:
        r_acc = accel_offset_from_steering(
            p.screen_offset_m, R_scr_m, R_acc_m, p.steering_arc_m
        )
        scr.append(Aperture(offset_m=p.screen_offset_m, radius_m=p.screen_radius_m))
        acc.append(Aperture(offset_m=r_acc,             radius_m=p.accel_radius_m))
        if mirror and abs(p.screen_offset_m) > 1e-12:
            r_acc_neg = accel_offset_from_steering(
                -p.screen_offset_m, R_scr_m, R_acc_m, -p.steering_arc_m
            )
            scr.append(Aperture(offset_m=-p.screen_offset_m, radius_m=p.screen_radius_m))
            acc.append(Aperture(offset_m=r_acc_neg,          radius_m=p.accel_radius_m))
    return scr, acc


def concentric_accel_radius(
    R_scr_m: float,
    screen_thickness_m: float,
    gap_after_m: float,
    screen_offsets_m: List[float],
) -> float:
    """Radius of curvature for an accel grid concentric with the screen grid.

    In a concentric sphere design both grids share a single center of curvature
    so that every radial line from that center passes through the aligned
    screen *and* accel aperture centres.  The accel sphere radius is:

        R_acc = R_scr − Δx

    where Δx = screen_thickness + gap_after_m (the vertex gap).

    ``gap_after_m`` is the on-axis distance from the screen downstream vertex
    to the accel upstream vertex.  The C++ solver positions each curved
    aperture via ``curved_surface_x()``, so no sagitta correction is applied.

    Parameters
    ----------
    R_scr_m:
        Screen grid radius of curvature in metres.
    screen_thickness_m:
        Axial thickness of the screen grid electrode in metres.
    gap_after_m:
        On-axis (vertex) gap between screen downstream face and accel upstream
        face in metres.
    screen_offsets_m:
        Unused. Retained for API compatibility.

    Returns
    -------
    float
        Accel radius of curvature in metres (R_scr − Δx).
    """
    delta_x = screen_thickness_m + gap_after_m
    return R_scr_m - delta_x


def concentric_accel_offsets(
    screen_offsets_m: List[float],
    R_scr_m: float,
    R_acc_m: float,
) -> List[float]:
    """Compute accel aperture y-offsets for perfect concentric-sphere bore alignment.

    For a screen aperture at radial offset r_scr the corresponding accel aperture
    must sit at:

        r_acc = (R_acc / R_scr) × r_scr

    so that the bore axis — a radial line from the shared center of curvature —
    passes exactly through both aperture centres.  Because r_acc/R_acc = r_scr/R_scr,
    the tilt angle α is identical on both grids, so the C++ normal-alignment check
    uses the same bore direction for screen and accel apertures.

    Parameters
    ----------
    screen_offsets_m:
        List of screen radial offsets (signed or unsigned).
    R_scr_m:
        Screen radius of curvature [m].
    R_acc_m:
        Accel radius of curvature [m], typically from :func:`concentric_accel_radius`.

    Returns
    -------
    List[float]
        Accel offsets in the same order as *screen_offsets_m*.
    """
    scale = R_acc_m / R_scr_m
    return [r * scale for r in screen_offsets_m]


def uniform_screen_offsets(
    n: int,
    pitch_m: float,
    *,
    offset_m: float = 0.0,
) -> List[float]:
    """Return *n* radial offsets spaced uniformly by *pitch_m* starting at *offset_m*.

    Parameters
    ----------
    n:
        Number of unique aperture positions.  Must be >= 1.
    pitch_m:
        Centre-to-centre spacing between adjacent apertures [m].
    offset_m:
        Radial position of the first aperture [m].  0.0 places the first
        aperture on-axis; positive values shift the entire array outward.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    return [offset_m + i * pitch_m for i in range(n)]


# ---------------------------------------------------------------------------
# Curvature dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GridCurvature:
    """Spherical curvature of a grid electrode.

    radius_m: radius of curvature in metres.  Use math.inf for a flat grid.
        Must be strictly positive or inf.
    concave_upstream: when True the vertex of the sphere points toward the ion
        source (upstream).  Apertures at larger radii tilt *inward*, producing
        a net focusing effect.  When False the vertex points downstream and
        apertures tilt outward (diverging).
    aperture_alignment: controls the orientation of each aperture bore through
        the curved grid material.

        ``'normal'`` — bores are drilled perpendicular to the curved surface
        (each bore axis follows the local surface normal, tilted by
        alpha = arcsin(r/R) from the beam axis).  This is the physically
        correct choice for a precision-machined curved grid and is the default.

        ``'axial'`` — bores are drilled parallel to the beam axis regardless
        of the surface curvature (equivalent to the flat-grid aperture
        orientation).  This models a grid where the holes were punched straight
        through without accounting for the dome shape.
    """
    radius_m: float = math.inf
    concave_upstream: bool = True
    aperture_alignment: str = 'normal'

    def __post_init__(self) -> None:
        if not (math.isinf(self.radius_m) or self.radius_m > 0.0):
            raise ValueError(
                f"GridCurvature.radius_m must be positive or math.inf, got {self.radius_m!r}"
            )
        if self.aperture_alignment not in ('normal', 'axial'):
            raise ValueError(
                f"GridCurvature.aperture_alignment must be 'normal' or 'axial', "
                f"got {self.aperture_alignment!r}"
            )

    @property
    def is_flat(self) -> bool:
        return math.isinf(self.radius_m)


# ---------------------------------------------------------------------------
# Extended grid definition
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CurvedGridDefinition:
    """GridDefinition with an optional spherical curvature.

    All fields from GridDefinition are replicated here (field-duplication avoids
    frozen dataclass inheritance issues in Python).

    A CurvedGridDefinition with curvature=GridCurvature() (default radius=inf)
    is functionally identical to GridDefinition.
    """
    name: str
    voltage_v: float
    thickness_m: float
    gap_after_m: float = 0.0
    apertures: List[Aperture] = field(default_factory=list)
    upstream_chamfer: Chamfer = field(default_factory=Chamfer)
    downstream_chamfer: Chamfer = field(default_factory=Chamfer)
    mirror: bool = False
    curvature: GridCurvature = field(default_factory=GridCurvature)


@dataclass(frozen=True)
class CurvedSimulationCase:
    """SimulationCase that supports CurvedGridDefinition entries."""
    grids: List[CurvedGridDefinition]
    env: Dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Per-aperture curved geometry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CurvedAperturePosition:
    """Computed geometry for a single aperture on a spherically curved grid.

    offset_radial_m : radial distance of the aperture centre from the grid
        symmetry axis.  Equal to the signed Aperture.offset_m from the input.
    offset_axial_m  : axial displacement of the aperture centre from the grid
        vertex (Δz).  Always 0 for flat grids.  For concave_upstream=True the
        displacement is positive (aperture moved downstream); for
        concave_upstream=False it is negative (aperture moved upstream).
    normal_tilt_deg : angle (degrees) between the aperture's bore axis and the
        beam axis.  Positive means tilted inward (focusing).
    aperture_radius_m : bore radius, copied from the source Aperture.
    """
    offset_radial_m: float
    offset_axial_m: float
    normal_tilt_deg: float
    aperture_radius_m: float


def compute_curved_aperture_positions(
    apertures: List[Aperture],
    curvature: GridCurvature,
) -> List[CurvedAperturePosition]:
    """Compute per-aperture geometry on a spherically curved grid.

    Parameters
    ----------
    apertures:
        List of Aperture objects defining radial offsets and bore radii.
    curvature:
        GridCurvature describing the spherical curvature of the grid.

    Returns
    -------
    One CurvedAperturePosition per input aperture (same order).

    Raises
    ------
    ValueError
        If any aperture's |offset_m| >= curvature.radius_m.
    """
    R = curvature.radius_m
    result: List[CurvedAperturePosition] = []
    for ap in apertures:
        r = abs(ap.offset_m)
        if curvature.is_flat:
            axial_offset = 0.0
            tilt_deg = 0.0
        else:
            if r >= R:
                raise ValueError(
                    f"Aperture at |offset|={r} m is outside the radius of "
                    f"curvature R={R} m (arcsin domain exceeded)."
                )
            # Full spherical formula — no paraxial approximation
            axial_offset = R * (1.0 - math.cos(math.asin(r / R)))
            tilt_deg = math.degrees(math.asin(r / R))
            if not curvature.concave_upstream:
                axial_offset = -axial_offset
                tilt_deg = -tilt_deg

        result.append(CurvedAperturePosition(
            offset_radial_m=ap.offset_m,
            offset_axial_m=axial_offset,
            normal_tilt_deg=tilt_deg,
            aperture_radius_m=ap.radius_m,
        ))
    return result


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

def _serialize_curvature(curvature: GridCurvature) -> str:
    """Return the 9th pipe field: 'R:concave_flag:align_flag'.

    Three colon-separated parts:
      R             — radius in metres, or ``inf`` for flat
      concave_flag  — ``1`` = concave_upstream, ``0`` = convex_upstream
      align_flag    — ``1`` = normal-aligned bores, ``0`` = axial-aligned bores
    """
    r_str = "inf" if curvature.is_flat else _format_num(curvature.radius_m)
    concave = "1" if curvature.concave_upstream else "0"
    align = "1" if curvature.aperture_alignment == 'normal' else "0"
    return f"{r_str}:{concave}:{align}"


def flatten_curved_grid(grid: CurvedGridDefinition) -> GridDefinition:
    """Strip curvature from a CurvedGridDefinition, returning a flat GridDefinition.

    The aperture radial offsets are preserved as-is.  Use this when the C++
    binary does not yet support the 9th curvature field.
    """
    return GridDefinition(
        name=grid.name,
        voltage_v=grid.voltage_v,
        thickness_m=grid.thickness_m,
        gap_after_m=grid.gap_after_m,
        apertures=list(grid.apertures),
        upstream_chamfer=grid.upstream_chamfer,
        downstream_chamfer=grid.downstream_chamfer,
        mirror=grid.mirror,
    )


def serialize_curved_grid_stack(
    grids: List[CurvedGridDefinition],
) -> str:
    """Serialize a list of CurvedGridDefinition objects to a GRID_STACK string.

    Curved grids append a 9th pipe-delimited field ``curv_r:concave_flag`` to
    each grid token.  Flat grids (radius=inf) still append this field so the
    binary can always rely on field[8] being present when operating in curved
    mode.

    ``gap_after_m`` is passed directly to the C++ binary as the vertex gap —
    the on-axis distance from the downstream face of this grid to the upstream
    face of the next.  The C++ solver uses ``curved_surface_x()`` to compute
    each aperture's actual axial position on the sphere, so no sagitta
    correction is needed here.

    Parameters
    ----------
    grids:
        List of 2–4 CurvedGridDefinition objects.

    Returns
    -------
    str
        The GRID_STACK string for the multi_grid_2d binary.

    Raises
    ------
    ValueError
        If len(grids) < 2 or > 4.
    """
    if not 2 <= len(grids) <= 4:
        raise ValueError("serialize_curved_grid_stack expects 2 to 4 grids")

    encoded: List[str] = []
    for grid in grids:
        curved_positions = compute_curved_aperture_positions(
            grid.apertures, grid.curvature
        )
        # Pass gap_after_m directly — the C++ binary positions each aperture
        # on the curved surface via curved_surface_x(), so gap_after_m is the
        # vertex gap and needs no sagitta correction here.
        effective_gap = grid.gap_after_m

        # Build aperture token (offset:radius pairs)
        if not grid.apertures:
            raise ValueError(f"grid '{grid.name}' has no apertures")
        ap_tok = ",".join(
            f"{_format_num(ap.offset_m)}:{_format_num(ap.radius_m)}"
            for ap in grid.apertures
        )
        up_ch = f"{_format_num(grid.upstream_chamfer.depth_m)}:{_format_num(grid.upstream_chamfer.angle_deg)}"
        dn_ch = f"{_format_num(grid.downstream_chamfer.depth_m)}:{_format_num(grid.downstream_chamfer.angle_deg)}"
        mirror_flag = "1" if grid.mirror else "0"
        curv_tok = _serialize_curvature(grid.curvature)

        encoded.append("|".join([
            grid.name,
            _format_num(grid.voltage_v),
            _format_num(grid.thickness_m),
            _format_num(effective_gap),
            ap_tok,
            up_ch,
            dn_ch,
            mirror_flag,
            curv_tok,
        ]))
    return ";".join(encoded)


# ---------------------------------------------------------------------------
# Case loading and env building
# ---------------------------------------------------------------------------

def load_curved_case(config_path: Path) -> CurvedSimulationCase:
    """Load a CurvedSimulationCase from a Python config file.

    The file must define ``CASE = CurvedSimulationCase(...)``.  Falls back
    gracefully if CASE is a plain SimulationCase (wraps grids in
    CurvedGridDefinition with flat curvature).

    Uses attribute inspection (not isinstance) to reconstruct the case from
    its fields, which avoids the ``__main__`` vs module identity issue that
    arises when the pipeline is run as a script.
    """
    namespace: Dict = {
        "Aperture": Aperture,
        "Chamfer": Chamfer,
        "GridDefinition": GridDefinition,
        "BeamletPair": BeamletPair,
        "make_apertures_from_pairs": make_apertures_from_pairs,
        "concentric_accel_radius": concentric_accel_radius,
        "concentric_accel_offsets": concentric_accel_offsets,
        "GridCurvature": GridCurvature,
        "CurvedGridDefinition": CurvedGridDefinition,
        "CurvedSimulationCase": CurvedSimulationCase,
        "SimulationCase": SimulationCase,
    }
    code = compile(config_path.read_text(encoding="utf-8"), str(config_path), "exec")
    exec(code, namespace)
    case = namespace.get("CASE")
    if case is None or not hasattr(case, "grids") or not hasattr(case, "env"):
        raise TypeError(
            f"{config_path} must define CASE = CurvedSimulationCase(...)"
        )
    curved_grids: List[CurvedGridDefinition] = []
    for g in case.grids:
        if hasattr(g, "curvature"):
            # Already a CurvedGridDefinition (or duck-type equivalent)
            curved_grids.append(CurvedGridDefinition(
                name=g.name,
                voltage_v=g.voltage_v,
                thickness_m=g.thickness_m,
                gap_after_m=g.gap_after_m,
                apertures=list(g.apertures),
                upstream_chamfer=g.upstream_chamfer,
                downstream_chamfer=g.downstream_chamfer,
                mirror=g.mirror,
                curvature=GridCurvature(
                    radius_m=g.curvature.radius_m,
                    concave_upstream=g.curvature.concave_upstream,
                    aperture_alignment=getattr(
                        g.curvature, 'aperture_alignment', 'normal'
                    ),
                ),
            ))
        else:
            # Plain GridDefinition — wrap with flat curvature
            curved_grids.append(CurvedGridDefinition(
                name=g.name,
                voltage_v=g.voltage_v,
                thickness_m=g.thickness_m,
                gap_after_m=g.gap_after_m,
                apertures=list(g.apertures),
                upstream_chamfer=g.upstream_chamfer,
                downstream_chamfer=g.downstream_chamfer,
                mirror=g.mirror,
            ))
    return CurvedSimulationCase(grids=curved_grids, env=dict(case.env))


def _run_profile_plot(env: Dict[str, str], repo_root: Path) -> None:
    """Generate sample_profile_I_Apm.html next to the run's sample_diameter_profile.json.

    Reconstructs the output directory using the same naming logic as the C++ binary
    (RESULTS_DIR / RUN_PREFIX[_RUN_TAG]_RUN_STAMP[_jSLURM_JOB_ID]).

    Self-contained: reads the JSON with stdlib only and writes a standalone HTML
    using Plotly via CDN — no local plotly/pandas install required.
    """
    import json as _json

    results_base = env.get("RESULTS_DIR", "results")
    stamp        = env.get("RUN_STAMP", "")
    prefix       = env.get("RUN_PREFIX", "run")
    tag          = env.get("RUN_TAG", "")
    jobid        = env.get("SLURM_JOB_ID", "")

    runname = prefix
    if tag:
        runname += f"_{tag}"
    runname += f"_{stamp}"
    if jobid:
        runname += f"_j{jobid}"

    base_path = Path(results_base)
    if not base_path.is_absolute():
        base_path = repo_root / base_path
    outdir = base_path / runname

    json_file = outdir / "sample_diameter_profile.json"
    if not json_file.exists():
        print(f"[profile-plot] sample_diameter_profile.json not found at {outdir} — skipping.", flush=True)
        return

    out_html = outdir / "sample_profile_I_Apm.html"
    print(f"[profile-plot] generating {out_html.name} ...", flush=True)

    try:
        data = _json.loads(json_file.read_text())
        bins = data.get("bins", [])
        if not bins:
            print("[profile-plot] no bins in profile JSON — skipping.", flush=True)
            return

        xs     = [(b["y_lo_m"] + b["y_hi_m"]) * 0.5 * 1000.0 for b in bins]  # mm
        ys     = [b["I_Apm"] for b in bins]
        colors = ["royalblue" if b.get("in_sample", True) else "crimson" for b in bins]
        bw_mm  = data.get("bin_width_m", (xs[1] - xs[0]) / 1000.0 if len(xs) > 1 else 1.0)
        if isinstance(bw_mm, float):
            bw_mm = bw_mm * 1000.0

        stats  = data.get("stats", {})
        p2a    = stats.get("peak_to_avg", float("nan"))
        rms_nu = stats.get("rms_nonuniform", float("nan"))
        title  = (f"Sample Diameter Profile — I_Apm<br>"
                  f"<sub>peak/avg={p2a:.3f}  rms_nonuniform={rms_nu:.3f}  "
                  f"total_I={data.get('total_I_A', 0):.4f} A</sub>")

        # Build a self-contained HTML using Plotly CDN — no local install needed.
        bar_traces = []
        # Group consecutive same-colour bins into single traces for a compact legend.
        for label, colour in (("in sample", "royalblue"), ("outside sample", "crimson")):
            xi = [x for x, c in zip(xs, colors) if c == colour]
            yi = [y for y, c in zip(ys, colors) if c == colour]
            if xi:
                bar_traces.append(
                    f'{{"type":"bar","name":"{label}","x":{xi},"y":{yi},'
                    f'"marker":{{"color":"{colour}"}},"width":{bw_mm},"showlegend":true}}'
                )

        traces_js = ",".join(bar_traces)
        html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Sample Diameter Profile</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
</head>
<body>
<div id="plot" style="width:100%;height:600px;"></div>
<script>
var traces = [{traces_js}];
var layout = {{
  title: {{text: "{title}", font: {{size: 18}}}},
  xaxis: {{title: "y (mm)", tickfont: {{size: 16}}, titlefont: {{size: 18}}}},
  yaxis: {{title: "I (A/m)", tickfont: {{size: 16}}, titlefont: {{size: 18}}}},
  barmode: "overlay",
  bargap: 0,
  template: "plotly_white",
  legend: {{font: {{size: 14}}}}
}};
Plotly.newPlot("plot", traces, layout, {{responsive: true}});
</script>
</body>
</html>"""

        out_html.write_text(html, encoding="utf-8")
        print(f"[profile-plot] wrote {out_html}", flush=True)

    except Exception as exc:
        print(f"[profile-plot] ERROR: {exc}", flush=True)


def build_curved_env(case: CurvedSimulationCase) -> Dict[str, str]:
    """Build the environment dict for multi_grid_2d using curved serialization."""
    env = os.environ.copy()
    env["GRID_STACK"] = serialize_curved_grid_stack(case.grids)
    for key, value in case.env.items():
        env[key] = str(value)
    return env


# ---------------------------------------------------------------------------
# Curved case utilities: round-trip import / export
# ---------------------------------------------------------------------------

def pairs_from_case(case: CurvedSimulationCase) -> List[BeamletPair]:
    """Extract a :class:`BeamletPair` list from a loaded :class:`CurvedSimulationCase`.

    The first grid is treated as the screen electrode and the second as the
    accel electrode.  Apertures are matched by position index; the result is
    truncated to ``min(len(screen.apertures), len(accel.apertures))``.

    Parameters
    ----------
    case:
        A loaded :class:`CurvedSimulationCase` (output of
        :func:`load_curved_case`).

    Returns
    -------
    List[BeamletPair]
        One pair per matched aperture, in the same order as the aperture
        lists.  If the grids were defined with ``mirror=True`` the list will
        include the mirrored (negative-offset) entries.
    """
    if len(case.grids) < 2:
        raise ValueError(
            "pairs_from_case requires at least 2 grids (screen + accel)"
        )
    screen = case.grids[0]
    accel  = case.grids[1]
    R_scr = screen.curvature.radius_m
    R_acc = accel.curvature.radius_m
    n = min(len(screen.apertures), len(accel.apertures))
    pairs: List[BeamletPair] = []
    for i in range(n):
        r_scr = screen.apertures[i].offset_m
        r_acc = accel.apertures[i].offset_m
        # Compute the steering arc-length: deviation of the accel aperture from
        # the concentric-sphere (zero-steering) position along the accel sphere.
        if math.isinf(R_scr) or math.isinf(R_acc):
            scale = 1.0 if math.isinf(R_scr) else R_acc / R_scr
            steering_arc = r_acc - r_scr * scale
        else:
            abs_r_scr = abs(r_scr)
            abs_r_acc = abs(r_acc)
            alpha_scr = math.asin(min(1.0, abs_r_scr / R_scr))
            alpha_acc = math.asin(min(1.0, abs_r_acc / R_acc))
            sign = 1.0 if r_scr >= 0.0 else -1.0
            # Preserve sign: inward shift = negative arc for positive r
            delta_alpha = (alpha_acc - alpha_scr) * sign
            steering_arc = delta_alpha * R_acc
        pairs.append(BeamletPair(
            screen_offset_m=r_scr,
            steering_arc_m=round(steering_arc, 12),
            screen_radius_m=screen.apertures[i].radius_m,
            accel_radius_m=accel.apertures[i].radius_m,
        ))
    return pairs


def bore_clearance_report(case: CurvedSimulationCase) -> List[Dict]:
    """Compute bore clearance for every aperture in the case.

    For an aperture at radial offset *r* on a curved grid of radius *R* and
    axial thickness *t*, the maximum lateral displacement of the bore axis
    through the material is:

        bore_disp = t × tan(arcsin(r / R))

    The aperture passes the clearance check when its bore radius exceeds this
    displacement (i.e. ``margin = bore_radius − bore_disp > 0``).

    Parameters
    ----------
    case:
        A loaded :class:`CurvedSimulationCase`.

    Returns
    -------
    List[dict]
        One dict per aperture with keys:
        ``grid``, ``offset_m``, ``alpha_deg``, ``bore_disp_m``,
        ``bore_radius_m``, ``margin_m``, ``ok``.
    """
    rows: List[Dict] = []
    for grid in case.grids:
        R = grid.curvature.radius_m
        t = grid.thickness_m
        for ap in grid.apertures:
            r = abs(ap.offset_m)
            if math.isinf(R) or R <= 0.0 or r == 0.0:
                alpha_deg = 0.0
                bore_disp = 0.0
            elif r >= R:
                alpha_deg = 90.0
                bore_disp = float("inf")
            else:
                alpha_rad = math.asin(r / R)
                alpha_deg = math.degrees(alpha_rad)
                bore_disp = t * math.tan(alpha_rad)
            margin = ap.radius_m - bore_disp
            rows.append({
                "grid":          grid.name,
                "offset_m":      ap.offset_m,
                "alpha_deg":     alpha_deg,
                "bore_disp_m":   bore_disp,
                "bore_radius_m": ap.radius_m,
                "margin_m":      margin,
                "ok":            math.isfinite(margin) and margin > 0.0,
            })
    return rows


def write_case_file(
    path: Path,
    *,
    screen_offsets_m: List[float],
    accel_offsets_m: Optional[List[float]] = None,
    ap_radius_m: float = 0.0015,
    R_scr_m: float = 0.100,
    R_acc_m: Optional[float] = None,
    screen_voltage_v: float = 0.0,
    accel_voltage_v: float = -15000.0,
    screen_thickness_m: float = 0.005,
    accel_thickness_m: float = 0.005,
    gap_after_m: float = 0.006,
    env: Optional[Dict[str, str]] = None,
) -> None:
    """Write a curved grid case file (``.py``) pre-populated with aperture pairs.

    The output file is valid Python that can be loaded by
    :func:`load_curved_case`.  It uses named constants and is intended to be
    human-readable and directly editable before running the solver.

    Parameters
    ----------
    path:
        Output ``.py`` file path.
    screen_offsets_m:
        Radial offsets of screen apertures in metres (positive values only;
        mirroring is handled by ``mirror=True`` in the generated
        ``CurvedGridDefinition``).
    accel_offsets_m:
        Radial offsets of accel apertures in metres.  When *None*, the
        concentric-sphere formula is applied using *R_scr_m*, *R_acc_m*,
        *screen_thickness_m*, and *gap_after_m*.
    ap_radius_m:
        Bore radius in metres (shared by screen and accel apertures).
    R_scr_m:
        Screen radius of curvature in metres.
    R_acc_m:
        Accel radius of curvature in metres.  When *None* it is computed via
        :func:`concentric_accel_radius` from *R_scr_m*, *screen_thickness_m*,
        *gap_after_m*, and *screen_offsets_m*.
    screen_voltage_v, accel_voltage_v:
        Electrode voltages in volts.
    screen_thickness_m, accel_thickness_m:
        Electrode axial thicknesses in metres.
    gap_after_m:
        Nominal inter-grid gap at the screen vertex in metres.
    env:
        Optional dict of additional environment variables to embed in the
        ``CASE.env`` block (e.g. ``ION_M_AMU``, ``PLASMA_NI_M3``).
    """
    import datetime

    if R_acc_m is None:
        R_acc_m = concentric_accel_radius(
            R_scr_m, screen_thickness_m, gap_after_m, screen_offsets_m
        )

    if accel_offsets_m is None:
        accel_offsets_m = concentric_accel_offsets(
            screen_offsets_m, R_scr_m, R_acc_m
        )

    # Convert absolute accel offsets to steering arc-lengths
    steering_arcs_m: List[float] = []
    for r_scr, r_acc in zip(screen_offsets_m, accel_offsets_m):
        abs_r_scr = abs(r_scr)
        abs_r_acc = abs(r_acc)
        if math.isinf(R_scr_m) or math.isinf(R_acc_m):
            scale = 1.0 if math.isinf(R_scr_m) else R_acc_m / R_scr_m
            steering_arcs_m.append(r_acc - r_scr * scale)
        else:
            alpha_scr = math.asin(min(1.0, abs_r_scr / R_scr_m))
            alpha_acc = math.asin(min(1.0, abs_r_acc / R_acc_m))
            sign = 1.0 if r_scr >= 0.0 else -1.0
            steering_arcs_m.append((alpha_acc - alpha_scr) * sign * R_acc_m)

    def _fmt_list(vals: List[float]) -> str:
        return "[" + ", ".join(f"{v:.8g}" for v in vals) + "]"

    env_lines = ""
    if env:
        for k, v in sorted(env.items()):
            try:
                float(v)
                env_lines += f'        "{k}": {v},\n'
            except (ValueError, TypeError):
                env_lines += f'        "{k}": "{v}",\n'

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    filename  = Path(path).name

    content = (
        f'"""\n'
        f'Curved grid case file — generated by beamlet_layout_gui  [{timestamp}]\n'
        f'\n'
        f'Load with:  python curved_grid_pipeline.py --config {filename}\n'
        f'Dry-run:    python curved_grid_pipeline.py --config {filename} --dry-run\n'
        f'Profile:    python curved_grid_pipeline.py --config {filename} --plot profile.png\n'
        f'"""\n'
        f'import math\n'
        f'from curved_grid_pipeline import (\n'
        f'    Aperture, Chamfer, BeamletPair, GridCurvature,\n'
        f'    CurvedGridDefinition, CurvedSimulationCase,\n'
        f'    make_apertures_from_pairs,\n'
        f')\n'
        f'\n'
        f'# \u2500\u2500 Electrode geometry \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n'
        f'R_SCR_M         = {R_scr_m:.8g}   # screen radius of curvature (m)\n'
        f'R_ACC_M         = {R_acc_m:.8g}   # accel  radius of curvature (m)  [concentric]\n'
        f'SCREEN_T_M      = {screen_thickness_m:.8g}\n'
        f'ACCEL_T_M       = {accel_thickness_m:.8g}\n'
        f'GAP_AFTER_M     = {gap_after_m:.8g}\n'
        f'VS_V            = {screen_voltage_v:.1f}\n'
        f'VA_V            = {accel_voltage_v:.1f}\n'
        f'AP_RADIUS_M     = {ap_radius_m:.8g}\n'
        f'\n'
        f'# \u2500\u2500 Aperture offsets (m) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n'
        f'SCREEN_OFFSETS_M = {_fmt_list(screen_offsets_m)}\n'
        f'# Arc-length on the accel sphere from the concentric (zero-steering) position.\n'
        f'# Zero = radial line from focal point passes through both aperture centres.\n'
        f'# Positive = accel aperture displaced outward; negative = inward.\n'
        f'STEERING_ARCS_M  = {_fmt_list(steering_arcs_m)}\n'
        f'\n'
        f'PAIRS = [\n'
        f'    BeamletPair(screen_offset_m=s, steering_arc_m=arc, screen_radius_m=AP_RADIUS_M)\n'
        f'    for s, arc in zip(SCREEN_OFFSETS_M, STEERING_ARCS_M)\n'
        f']\n'
        f'screen_apertures, accel_apertures = make_apertures_from_pairs(\n'
        f'    PAIRS, mirror=False, R_scr_m=R_SCR_M, R_acc_m=R_ACC_M\n'
        f')\n'
        f'\n'
        f'# \u2500\u2500 Simulation case \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n'
        f'CASE = CurvedSimulationCase(\n'
        f'    grids=[\n'
        f'        CurvedGridDefinition(\n'
        f'            name="screen",\n'
        f'            voltage_v=VS_V,\n'
        f'            thickness_m=SCREEN_T_M,\n'
        f'            gap_after_m=GAP_AFTER_M,\n'
        f'            apertures=screen_apertures,\n'
        f'            curvature=GridCurvature(radius_m=R_SCR_M, concave_upstream=True),\n'
        f'            mirror=True,\n'
        f'        ),\n'
        f'        CurvedGridDefinition(\n'
        f'            name="accel",\n'
        f'            voltage_v=VA_V,\n'
        f'            thickness_m=ACCEL_T_M,\n'
        f'            gap_after_m=0.0,\n'
        f'            apertures=accel_apertures,\n'
        f'            curvature=GridCurvature(radius_m=R_ACC_M, concave_upstream=True),\n'
        f'            mirror=True,\n'
        f'        ),\n'
        f'    ],\n'
        f'    env={{\n'
        f'{env_lines}'
        f'    }},\n'
        f')\n'
    )
    Path(path).write_text(content, encoding="utf-8")


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_curved_grid_profile(
    grids: List[CurvedGridDefinition],
    *,
    units: str = "mm",
    show_normals: bool = True,
    show_bore_alignments: bool = True,
    normal_length_m: float = 0.003,
    figsize: Tuple[float, float] = (12.0, 7.0),
    ax=None,
    title: Optional[str] = None,
):
    """Render a 2-D cross-section of the curved grid stack.

    For each grid, draws the curved upstream and downstream faces, masks out
    aperture holes, fills the body, and optionally overlays aperture normal
    vectors.

    Parameters
    ----------
    grids:
        List of CurvedGridDefinition objects (2–4 grids).
    units:
        ``'mm'`` (default) or ``'m'`` — controls axis scaling and labels.
    show_normals:
        Whether to draw aperture normal direction arrows.
    show_bore_alignments:
        When True (default) and at least two grids are present, draws:

        * A dashed line through each screen aperture centre along the surface
          normal (the radial line from the focal point), extending to the edges
          of the plot — showing the ideal zero-steering bore axis.
        * For each paired accel aperture, a perpendicular indicator from the
          point on the bore axis to the actual aperture centre — showing how
          far the accel aperture deviates from the concentric position.  The
          indicator is omitted when the steering arc is negligibly small.
    normal_length_m:
        Length of normal arrows in metres.
    figsize:
        Matplotlib figure size in inches.
    ax:
        Existing Axes to draw into.  Creates a new figure when None.
    title:
        Optional figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if not _HAS_MPL:
        raise ImportError(
            "matplotlib and numpy are required for plot_curved_grid_profile"
        )

    scale = 1e3 if units == "mm" else 1.0
    unit_label = "mm" if units == "mm" else "m"

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Compute axial cursor (x position of each grid's vertex)
    x_cursor = 0.0
    grid_x_starts: List[float] = []
    for grid in grids:
        grid_x_starts.append(x_cursor)
        curved_positions = compute_curved_aperture_positions(
            grid.apertures, grid.curvature
        )
        max_sagitta = max(
            (abs(cp.offset_axial_m) for cp in curved_positions), default=0.0
        )
        x_cursor += grid.thickness_m + max_sagitta + grid.gap_after_m

    # Y range for drawing grid bodies
    all_offsets = [abs(ap.offset_m) for g in grids for ap in g.apertures]
    y_max_ap = max(all_offsets, default=0.010)
    y_draw = y_max_ap * 1.35
    n_pts = 400
    y_dense = np.linspace(-y_draw, y_draw, n_pts)

    colors = ["#607D8B", "#455A64", "#37474F", "#263238"]

    for gi, grid in enumerate(grids):
        x0 = grid_x_starts[gi]
        R = grid.curvature.radius_m
        concave_up = grid.curvature.concave_upstream
        color = colors[gi % len(colors)]

        # --- Curved surface x as function of y (dense array) ---
        if grid.curvature.is_flat:
            x_upstream = np.full_like(y_dense, x0)
        else:
            ry = np.clip(y_dense / R, -1.0, 1.0)
            dz = R * (1.0 - np.sqrt(1.0 - ry ** 2))
            x_upstream = x0 + dz if concave_up else x0 - dz

        x_downstream = x_upstream + grid.thickness_m

        # Compute aperture positions on this grid
        curved_positions = compute_curved_aperture_positions(
            grid.apertures, grid.curvature
        )
        # Also add mirrored positions for display when mirror=True
        display_positions = list(curved_positions)
        if grid.mirror:
            for cp in curved_positions:
                if abs(cp.offset_radial_m) > 1e-12:
                    display_positions.append(CurvedAperturePosition(
                        offset_radial_m=-cp.offset_radial_m,
                        offset_axial_m=cp.offset_axial_m,
                        normal_tilt_deg=-cp.normal_tilt_deg,
                        aperture_radius_m=cp.aperture_radius_m,
                    ))

        # --- Fill full grid body (no masking) ---
        ax.fill_betweenx(
            y_dense * scale,
            x_upstream * scale,
            x_downstream * scale,
            color=color,
            alpha=0.85,
            linewidth=0,
            zorder=2,
        )
        # Upstream and downstream face outlines
        ax.plot(x_upstream * scale, y_dense * scale,
                color=color, lw=0.8, alpha=0.6, zorder=3)
        ax.plot(x_downstream * scale, y_dense * scale,
                color=color, lw=0.8, alpha=0.6, zorder=3)

        # --- Cut aperture holes as white polygons ---
        # For 'normal' alignment: the bore is perpendicular to the curved surface.
        #   Each aperture is drawn as a PARALLELOGRAM — the upstream and downstream
        #   faces are offset along the tilted bore axis, so the slot leans inward.
        # For 'axial' alignment: the bore is parallel to the beam axis (x).
        #   Each aperture is drawn as an upright RECTANGLE — the slot cuts
        #   straight through parallel to x regardless of surface tilt.
        use_normal_align = grid.curvature.aperture_alignment == 'normal'
        t = grid.thickness_m

        bg = ax.get_facecolor()   # match axes background for "cut-out" effect

        # Extend aperture polygons beyond both grid faces by this margin so
        # the white cut-out fully covers the fill with no edge slivers.
        ext = t * 0.5

        for cp in display_positions:
            y_ap = cp.offset_radial_m
            r_ap = cp.aperture_radius_m
            x_ap = x0 + cp.offset_axial_m   # aperture centre on upstream curved face

            if use_normal_align and not grid.curvature.is_flat:
                # Bore axis direction:  d = (cos_a, -sin_a)
                # Width direction:      n = (sin_a,  cos_a)
                # Extend ±ext along the bore axis beyond both faces so the
                # parallelogram fully clears the curved fill on both sides.
                alpha_rad = math.radians(cp.normal_tilt_deg)
                sin_a = math.sin(alpha_rad)
                cos_a = math.cos(alpha_rad)
                # Upstream extension: step back ext along bore axis (upstream)
                x0e = x_ap - ext * cos_a
                y0e = y_ap + ext * sin_a
                total = t + 2 * ext   # full depth of cut-out along bore axis
                corners = [
                    (x0e - r_ap * sin_a,              y0e - r_ap * cos_a),
                    (x0e + r_ap * sin_a,              y0e + r_ap * cos_a),
                    (x0e + r_ap * sin_a + total * cos_a, y0e + r_ap * cos_a - total * sin_a),
                    (x0e - r_ap * sin_a + total * cos_a, y0e - r_ap * cos_a - total * sin_a),
                ]
            else:
                # Axial bore: upright rectangle extended ±ext along x
                corners = [
                    (x_ap - ext, y_ap - r_ap),
                    (x_ap - ext, y_ap + r_ap),
                    (x_ap + t + ext, y_ap + r_ap),
                    (x_ap + t + ext, y_ap - r_ap),
                ]

            scaled = [(cx * scale, cy * scale) for cx, cy in corners]
            patch = mpatches.Polygon(
                scaled, closed=True,
                facecolor=bg, edgecolor=color,
                linewidth=0.6, zorder=4,
            )
            ax.add_patch(patch)

        # --- Draw aperture normal arrows ---
        if show_normals:
            for cp in display_positions:
                r = cp.offset_radial_m
                dz = cp.offset_axial_m
                x_ap = (x0 + dz) * scale
                y_ap = r * scale
                alpha_rad = math.radians(cp.normal_tilt_deg)
                # Aperture axis direction (downstream-inward for concave-up at r>0)
                dx_arr = math.cos(alpha_rad) * normal_length_m * scale
                dy_arr = -math.sin(alpha_rad) * normal_length_m * scale
                # Arrow starts just upstream of aperture, points downstream
                ax.annotate(
                    "",
                    xy=(x_ap + dx_arr, y_ap + dy_arr),
                    xytext=(x_ap, y_ap),
                    arrowprops=dict(
                        arrowstyle="->",
                        color="royalblue",
                        lw=1.5,
                    ),
                )

        # --- Label grid ---
        x_label = (np.mean(x_upstream) + grid.thickness_m / 2) * scale
        ax.text(
            x_label,
            (y_draw + 0.002) * scale,
            grid.name,
            ha="center",
            va="bottom",
            fontsize=8,
            color=color,
        )

    ax.set_xlabel(f"Axial position ({unit_label})")
    ax.set_ylabel(f"Radial position ({unit_label})")
    ax.set_aspect("equal")
    ax.axhline(0, color="gray", lw=0.5, linestyle=":")
    ax.grid(True, alpha=0.3)

    # --- Bore axis alignment lines + steering offset perpendiculars ---
    _bore_axis_drawn = False
    _steer_drawn = False

    if show_bore_alignments and len(grids) >= 2:
        scr_grid = grids[0]
        acc_grid = grids[1]
        x0_scr = grid_x_starts[0]
        x0_acc = grid_x_starts[1]

        # Build display positions for screen and accel (mirrors included)
        def _disp(grid):
            pos = compute_curved_aperture_positions(grid.apertures, grid.curvature)
            out = list(pos)
            if grid.mirror:
                for cp in pos:
                    if abs(cp.offset_radial_m) > 1e-12:
                        out.append(CurvedAperturePosition(
                            offset_radial_m=-cp.offset_radial_m,
                            offset_axial_m=cp.offset_axial_m,
                            normal_tilt_deg=-cp.normal_tilt_deg,
                            aperture_radius_m=cp.aperture_radius_m,
                        ))
            return out

        scr_disp = _disp(scr_grid)
        acc_disp = _disp(acc_grid)
        n_pairs = min(len(scr_disp), len(acc_disp))

        # x span for line clipping — extend well past both sides and let
        # matplotlib's axes clip the line to the visible area.
        x_span = x_cursor + y_draw
        x_margin = y_draw * 0.15

        for i in range(n_pairs):
            cp_scr = scr_disp[i]
            cp_acc = acc_disp[i]

            # Screen aperture centre (metres)
            x_scr_ap = x0_scr + cp_scr.offset_axial_m
            y_scr_ap = cp_scr.offset_radial_m

            # Downstream unit normal of this screen aperture:
            #   d = (cos_α,  −sign(r) · sin_α)
            alpha_rad = math.radians(abs(cp_scr.normal_tilt_deg))
            sign_r = 1.0 if cp_scr.offset_radial_m >= 0.0 else -1.0
            cos_a = math.cos(alpha_rad)
            sin_a = math.sin(alpha_rad)
            nx, ny = cos_a, -sign_r * sin_a   # unit normal downstream

            # Dashed bore-axis line: parametric, clipped to the plot x range
            if cos_a > 1e-9:
                t_lo = (-x_margin - x_scr_ap) / cos_a
                t_hi = (x_span + x_margin - x_scr_ap) / cos_a
            else:
                t_lo, t_hi = -x_span, x_span
            lx = np.array([x_scr_ap + t_lo * nx, x_scr_ap + t_hi * nx])
            ly = np.array([y_scr_ap + t_lo * ny, y_scr_ap + t_hi * ny])

            ax.plot(
                lx * scale, ly * scale,
                color="#43A047", lw=0.9, linestyle="--", alpha=0.55, zorder=1,
                label="_nolegend_",
            )
            _bore_axis_drawn = True

            # Accel aperture centre
            x_acc_ap = x0_acc + cp_acc.offset_axial_m
            y_acc_ap = cp_acc.offset_radial_m

            # Project accel aperture onto the bore axis line
            dx = x_acc_ap - x_scr_ap
            dy = y_acc_ap - y_scr_ap
            t_proj = dx * nx + dy * ny
            x_foot = x_scr_ap + t_proj * nx
            y_foot = y_scr_ap + t_proj * ny

            # Perpendicular offset magnitude
            perp_x = x_acc_ap - x_foot
            perp_y = y_acc_ap - y_foot
            offset_mag = math.hypot(perp_x, perp_y)

            if offset_mag > 1e-9:
                # Draw a double-headed line from foot → actual accel position
                ax.plot(
                    [x_foot * scale, x_acc_ap * scale],
                    [y_foot * scale, y_acc_ap * scale],
                    color="#E64A19", lw=1.6, zorder=5, label="_nolegend_",
                )
                # Dot at the concentric (foot) position
                ax.plot(
                    x_foot * scale, y_foot * scale,
                    "o", color="#E64A19", markersize=3.5, zorder=5,
                )
                _steer_drawn = True

    legend_patches = [
        mpatches.Patch(color=colors[gi % len(colors)], alpha=0.7,
                       label=f"{g.name}  ({g.voltage_v:+.0f} V)")
        for gi, g in enumerate(grids)
    ]
    if show_normals:
        legend_patches.append(
            mpatches.Patch(color="royalblue", label="Aperture normal")
        )
    if _bore_axis_drawn:
        legend_patches.append(
            plt.Line2D([0], [0], color="#43A047", lw=0.9, linestyle="--",
                       alpha=0.55, label="Bore axis (concentric)")
        )
    if _steer_drawn:
        legend_patches.append(
            plt.Line2D([0], [0], color="#E64A19", lw=1.6,
                       label="Steering offset from concentric")
        )
    ax.legend(handles=legend_patches, fontsize=8, loc="upper right")

    if title:
        ax.set_title(title)

    return fig


# ---------------------------------------------------------------------------
# Parameter sweep support
# ---------------------------------------------------------------------------

def load_sweep_spec(config_path: Path):
    """Execute a sweep config file and return ``(make_case_fn, sweep_spec)``.

    The config file must define:

    * ``make_case(**overrides)`` — a callable that accepts keyword overrides for
      any of the module-level constants and returns a
      :class:`CurvedSimulationCase`.
    * ``SWEEP`` — a :class:`SweepSpec` instance that describes the axes and
      reduction strategy.

    The execution namespace is populated with all curved-grid types *plus*
    the sweep helpers so that the config file can reference them without an
    explicit import.

    Parameters
    ----------
    config_path:
        Path to the ``.py`` sweep config file.

    Returns
    -------
    tuple[callable, SweepSpec]
        ``(make_case_fn, sweep_spec)``

    Raises
    ------
    ImportError
        If ``sweep_types`` could not be imported (``_HAS_SWEEP`` is False).
    ValueError
        If the config file does not define both ``make_case`` and ``SWEEP``.
    """
    if not _HAS_SWEEP:
        raise ImportError(
            "sweep_types is not available — ensure sim/sweep_types.py exists "
            "and is importable."
        )

    namespace: Dict = {
        # curved-grid types
        "math": __import__("math"),
        "Aperture": Aperture,
        "Chamfer": Chamfer,
        "GridDefinition": GridDefinition,
        "SimulationCase": SimulationCase,
        "BeamletPair": BeamletPair,
        "make_apertures_from_pairs": make_apertures_from_pairs,
        "concentric_accel_radius": concentric_accel_radius,
        "concentric_accel_offsets": concentric_accel_offsets,
        "GridCurvature": GridCurvature,
        "CurvedGridDefinition": CurvedGridDefinition,
        "CurvedSimulationCase": CurvedSimulationCase,
        # sweep helpers
        "SweepAxis": SweepAxis,
        "SweepSpec": SweepSpec,
        "linspace": linspace,
        "arange": arange,
    }
    code = compile(
        Path(config_path).read_text(encoding="utf-8"),
        str(config_path),
        "exec",
    )
    exec(code, namespace)

    make_case_fn = namespace.get("make_case")
    sweep_spec = namespace.get("SWEEP")
    if make_case_fn is None:
        raise ValueError(
            f"{config_path} must define a callable 'make_case(**overrides)'."
        )
    if sweep_spec is None:
        raise ValueError(
            f"{config_path} must define 'SWEEP = SweepSpec(...)'."
        )
    return make_case_fn, sweep_spec


def write_sweep_case_file(
    path: Path,
    *,
    sweep_axes,                          # List[SweepAxis] or List[Tuple[str, str]]
    sweep_mode: str = "product",
    sweep_reduce: bool = False,
    sweep_tag: str = "",
    # Geometry — can be supplied directly or derived from source_case_path
    source_case_path: Optional[Path] = None,
    gui_overrides: Optional[Dict[str, str]] = None,
    screen_offsets_m: Optional[List[float]] = None,
    accel_offsets_m: Optional[List[float]] = None,
    ap_radius_m: float = 0.0015,
    R_scr_m: float = 0.100,
    R_acc_m: Optional[float] = None,
    screen_voltage_v: float = 0.0,
    accel_voltage_v: float = -15000.0,
    screen_thickness_m: float = 0.005,
    accel_thickness_m: float = 0.005,
    gap_after_m: float = 0.006,
    env: Optional[Dict[str, str]] = None,
) -> None:
    """Write a parameterized curved-grid sweep case file (``.py``).

    The generated file is valid Python that can be loaded by both
    :func:`load_curved_case` (single run via ``CASE = make_case()``) and
    :func:`load_sweep_spec` (parameter sweep via ``make_case`` + ``SWEEP``).

    Layout of the generated file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    1. Module-level constants (``R_SCR_M``, ``AP_RADIUS_M``, …).
    2. ``def make_case(**overrides)`` — each constant is retrieved via
       ``overrides.get("CONST_NAME", CONST_NAME)``, the aperture pairs are
       recomputed, and a ``CurvedSimulationCase`` is returned.
    3. ``CASE = make_case()`` for backward compatibility with
       :func:`load_curved_case`.
    4. ``from sweep_types import …`` and ``SWEEP = SweepSpec(…)`` block.

    Parameters
    ----------
    path:
        Output ``.py`` file path.
    sweep_axes:
        List of ``(param_name, values_repr)`` tuples.  *param_name* must match
        one of the module-level constants (e.g. ``"R_SCR_M"``).
        *values_repr* is the raw Python expression to embed verbatim in the
        ``SweepAxis`` call (e.g. ``"[0.08, 0.10, 0.12]"`` or
        ``"linspace(1e16, 5e16, 5)"``).
    sweep_mode:
        Passed to ``SweepSpec(mode=…)``.  Typical values: ``"grid"`` (full
        Cartesian product), ``"zip"`` (paired iteration).
    sweep_reduce:
        Passed to ``SweepSpec(reduce=…)``.  When True the SLURM orchestrator
        will submit a reduction job after the array completes.
    screen_offsets_m, accel_offsets_m, ap_radius_m, R_scr_m, R_acc_m,
    screen_voltage_v, accel_voltage_v, screen_thickness_m, accel_thickness_m,
    gap_after_m, env:
        Same semantics as :func:`write_case_file`.
    """
    import datetime as _dt

    # ------------------------------------------------------------------
    # Normalise sweep_axes: accept both List[SweepAxis] and List[Tuple]
    # ------------------------------------------------------------------
    _raw_axes = []
    for ax in sweep_axes:
        if hasattr(ax, "name") and hasattr(ax, "values"):
            # SweepAxis object from sweep_types
            _raw_axes.append((ax.name, repr(list(ax.values))))
        else:
            # Legacy (param_name, values_repr_string) tuple
            _raw_axes.append(tuple(ax))
    sweep_axes_tuples: List[Tuple[str, str]] = _raw_axes

    # ------------------------------------------------------------------
    # If a source case file is given, load geometry from it and apply
    # any GUI overrides on top.
    # ------------------------------------------------------------------
    if source_case_path is not None:
        _src = load_curved_case(source_case_path)
        _g0  = _src.grids[0]
        _g1  = _src.grids[1] if len(_src.grids) > 1 else None

        # Base values from the loaded case
        R_scr_m           = _g0.curvature.radius_m if not _g0.curvature.is_flat else float("inf")
        screen_thickness_m = _g0.thickness_m
        gap_after_m        = _g0.gap_after_m
        screen_voltage_v   = _g0.voltage_v
        accel_voltage_v    = _g1.voltage_v  if _g1 else accel_voltage_v
        accel_thickness_m  = _g1.thickness_m if _g1 else accel_thickness_m

        _bore_radii = [ap.radius_m for ap in _g0.apertures if ap.radius_m > 0]
        ap_radius_m = _bore_radii[0] if _bore_radii else ap_radius_m

        screen_offsets_m = sorted({abs(ap.offset_m) for ap in _g0.apertures})
        accel_offsets_m  = None   # will be recomputed below
        R_acc_m          = None

        if env is None:
            env = dict(_src.env)
        else:
            merged = dict(_src.env)
            merged.update(env)
            env = merged

        # Apply GUI field overrides (string values from the editor)
        _ov = gui_overrides or {}
        _GEOM_KEYS_MM  = {"R_SCR", "T_SCR", "T_ACC", "GAP", "AP_RAD"}
        _GEOM_KEY_MAP  = {
            "R_SCR":  "R_scr_m",
            "T_SCR":  "screen_thickness_m",
            "T_ACC":  "accel_thickness_m",
            "GAP":    "gap_after_m",
            "VS":     "screen_voltage_v",
            "VA":     "accel_voltage_v",
            "AP_RAD": "ap_radius_m",
        }
        for gui_key, local_name in _GEOM_KEY_MAP.items():
            if gui_key in _ov:
                try:
                    val_mm = float(_ov[gui_key])
                    val_m  = val_mm * 1e-3 if gui_key in _GEOM_KEYS_MM else val_mm
                    locals()[local_name]   # ensure name exists (suppress linter)
                    # Assign into local scope via explicit reassignment
                    if local_name == "R_scr_m":            R_scr_m            = val_m
                    elif local_name == "screen_thickness_m": screen_thickness_m = val_m
                    elif local_name == "accel_thickness_m":  accel_thickness_m  = val_m
                    elif local_name == "gap_after_m":        gap_after_m        = val_m
                    elif local_name == "screen_voltage_v":   screen_voltage_v   = val_m
                    elif local_name == "accel_voltage_v":    accel_voltage_v    = val_m
                    elif local_name == "ap_radius_m":        ap_radius_m        = val_m
                except (ValueError, KeyError):
                    pass
        # ENV overrides from GUI
        for gui_key in _ov:
            if gui_key not in _GEOM_KEY_MAP:
                if env is None:
                    env = {}
                env[gui_key] = _ov[gui_key]

    if screen_offsets_m is None:
        screen_offsets_m = [0.0]

    if R_acc_m is None:
        R_acc_m = concentric_accel_radius(
            R_scr_m, screen_thickness_m, gap_after_m, screen_offsets_m
        )

    if accel_offsets_m is None:
        accel_offsets_m = concentric_accel_offsets(
            screen_offsets_m, R_scr_m, R_acc_m
        )

    # Convert absolute accel offsets to steering arc-lengths
    steering_arcs_m: List[float] = []
    for r_scr, r_acc in zip(screen_offsets_m, accel_offsets_m):
        abs_r_scr = abs(r_scr)
        abs_r_acc = abs(r_acc)
        if math.isinf(R_scr_m) or math.isinf(R_acc_m):
            scale = 1.0 if math.isinf(R_scr_m) else R_acc_m / R_scr_m
            steering_arcs_m.append(r_acc - r_scr * scale)
        else:
            alpha_scr = math.asin(min(1.0, abs_r_scr / R_scr_m))
            alpha_acc = math.asin(min(1.0, abs_r_acc / R_acc_m))
            sign = 1.0 if r_scr >= 0.0 else -1.0
            steering_arcs_m.append((alpha_acc - alpha_scr) * sign * R_acc_m)

    def _fmt_list(vals: List[float]) -> str:
        return "[" + ", ".join(f"{v:.8g}" for v in vals) + "]"

    env_lines = ""
    if env:
        for k, v in sorted(env.items()):
            try:
                float(v)
                env_lines += f'        "{k}": overrides.get("{k}", {v}),\n'
            except (ValueError, TypeError):
                env_lines += f'        "{k}": overrides.get("{k}", "{v}"),\n'

    # Build the SWEEP block
    axes_lines = ""
    for param_name, values_repr in sweep_axes_tuples:
        axes_lines += f'    SweepAxis(name="{param_name}", values={values_repr}),\n'

    timestamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    filename = Path(path).name

    content = (
        f'"""\n'
        f'Curved grid sweep case — generated by curved_grid_pipeline  [{timestamp}]\n'
        f'\n'
        f'Single run:  python curved_grid_pipeline.py --config {filename}\n'
        f'Dry-run:     python curved_grid_pipeline.py --config {filename} --dry-run\n'
        f'Sweep:       python curved_grid_pipeline.py --config {filename} --sweep\n'
        f'Array task:  python curved_grid_pipeline.py --config {filename} \\\n'
        f'                 --sweep-task $SLURM_ARRAY_TASK_ID\n'
        f'"""\n'
        f'import math\n'
        f'from curved_grid_pipeline import (\n'
        f'    Aperture, Chamfer, BeamletPair, GridCurvature,\n'
        f'    CurvedGridDefinition, CurvedSimulationCase,\n'
        f'    make_apertures_from_pairs,\n'
        f'    concentric_accel_radius, concentric_accel_offsets,\n'
        f')\n'
        f'\n'
        f'# ── Electrode geometry ──────────────────────────────────────────\n'
        f'R_SCR_M         = {R_scr_m:.8g}   # screen radius of curvature (m)\n'
        f'R_ACC_M         = {R_acc_m:.8g}   # accel  radius of curvature (m)  [concentric]\n'
        f'SCREEN_T_M      = {screen_thickness_m:.8g}\n'
        f'ACCEL_T_M       = {accel_thickness_m:.8g}\n'
        f'GAP_AFTER_M     = {gap_after_m:.8g}\n'
        f'VS_V            = {screen_voltage_v:.1f}\n'
        f'VA_V            = {accel_voltage_v:.1f}\n'
        f'AP_RADIUS_M     = {ap_radius_m:.8g}\n'
        f'\n'
        f'# ── Aperture offsets (m) ─────────────────────────────────────────\n'
        f'SCREEN_OFFSETS_M = {_fmt_list(screen_offsets_m)}\n'
        f'# Arc-length on the accel sphere from the concentric (zero-steering) position.\n'
        f'STEERING_ARCS_M  = {_fmt_list(steering_arcs_m)}\n'
        f'\n'
        f'\n'
        f'# ── Parameterised case factory ────────────────────────────────────\n'
        f'def make_case(**overrides):\n'
        f'    """Return a CurvedSimulationCase built from module constants.\n'
        f'\n'
        f'    Each constant is replaced by the corresponding value in *overrides*\n'
        f'    when supplied, allowing the sweep runner to vary individual\n'
        f'    parameters without mutating module state.\n'
        f'    """\n'
        f'    r_scr_m      = overrides.get("R_SCR_M",     R_SCR_M)\n'
        f'    r_acc_m_base = overrides.get("R_ACC_M",     R_ACC_M)\n'
        f'    screen_t_m   = overrides.get("SCREEN_T_M",  SCREEN_T_M)\n'
        f'    accel_t_m    = overrides.get("ACCEL_T_M",   ACCEL_T_M)\n'
        f'    gap_after_m  = overrides.get("GAP_AFTER_M", GAP_AFTER_M)\n'
        f'    vs_v         = overrides.get("VS_V",         VS_V)\n'
        f'    va_v         = overrides.get("VA_V",         VA_V)\n'
        f'    ap_radius_m  = overrides.get("AP_RADIUS_M", AP_RADIUS_M)\n'
        f'\n'
        f'    # Re-derive R_ACC_M when R_SCR_M or geometry changes\n'
        f'    if "R_ACC_M" not in overrides:\n'
        f'        r_acc_m = concentric_accel_radius(\n'
        f'            r_scr_m, screen_t_m, gap_after_m, SCREEN_OFFSETS_M\n'
        f'        )\n'
        f'    else:\n'
        f'        r_acc_m = r_acc_m_base\n'
        f'\n'
        f'    pairs = [\n'
        f'        BeamletPair(\n'
        f'            screen_offset_m=s,\n'
        f'            steering_arc_m=arc,\n'
        f'            screen_radius_m=ap_radius_m,\n'
        f'        )\n'
        f'        for s, arc in zip(SCREEN_OFFSETS_M, STEERING_ARCS_M)\n'
        f'    ]\n'
        f'    screen_apertures, accel_apertures = make_apertures_from_pairs(\n'
        f'        pairs, mirror=False, R_scr_m=r_scr_m, R_acc_m=r_acc_m\n'
        f'    )\n'
        f'\n'
        f'    return CurvedSimulationCase(\n'
        f'        grids=[\n'
        f'            CurvedGridDefinition(\n'
        f'                name="screen",\n'
        f'                voltage_v=vs_v,\n'
        f'                thickness_m=screen_t_m,\n'
        f'                gap_after_m=gap_after_m,\n'
        f'                apertures=screen_apertures,\n'
        f'                curvature=GridCurvature(radius_m=r_scr_m, concave_upstream=True),\n'
        f'                mirror=True,\n'
        f'            ),\n'
        f'            CurvedGridDefinition(\n'
        f'                name="accel",\n'
        f'                voltage_v=va_v,\n'
        f'                thickness_m=accel_t_m,\n'
        f'                gap_after_m=0.0,\n'
        f'                apertures=accel_apertures,\n'
        f'                curvature=GridCurvature(radius_m=r_acc_m, concave_upstream=True),\n'
        f'                mirror=True,\n'
        f'            ),\n'
        f'        ],\n'
        f'        env={{\n'
        f'{env_lines}'
        f'        }},\n'
        f'    )\n'
        f'\n'
        f'\n'
        f'# ── Backward-compatible single-run case ────────────────────────────\n'
        f'CASE = make_case()\n'
        f'\n'
        f'\n'
        f'# ── Sweep specification ──────────────────────────────────────────\n'
        f'from sweep_types import SweepAxis, SweepSpec, linspace, arange  # noqa: E402\n'
        f'\n'
        f'SWEEP = SweepSpec(\n'
        f'    axes=[\n'
        f'{axes_lines}'
        f'    ],\n'
        f'    mode="{sweep_mode}",\n'
        f'    reduce={sweep_reduce!r},\n'
        f'    tag="{sweep_tag}",\n'
        f')\n'
    )
    Path(path).write_text(content, encoding="utf-8")


def _make_run_tag(idx: int, params: Dict) -> str:
    """Build a compact RUN_TAG encoding the task index and swept values.

    Example: t0003_AP_RADIUS_M-0.0015_R_SCR_M-0.1
    """
    parts = [f"t{idx:04d}"]
    for k, v in params.items():
        val_str = f"{v:.6g}" if isinstance(v, float) else str(v)
        parts.append(f"{k}-{val_str}")
    return "_".join(parts)


def run_sweep(
    config_path: Path,
    binary_path: Path,
    pipeline_dir: Path,
) -> int:
    """Generate a sweep manifest and submit a SLURM array job.

    Steps
    ~~~~~
    1. Load the sweep spec from *config_path* via :func:`load_sweep_spec`.
    2. Expand all sweep axes into a flat list of parameter dicts.
    3. Write ``sweep_manifest.json`` to *pipeline_dir*.
    4. Submit the SLURM array via ``sbatch``.
    5. If ``sweep_spec.reduce`` is True, submit a dependent reduction job.

    Parameters
    ----------
    config_path:
        Path to the ``.py`` sweep config file.
    binary_path:
        Path to the multi-grid solver binary that each array task will run.
    pipeline_dir:
        Directory containing ``orchestrate_curved_grid_sweep.slurm`` and
        where ``sweep_manifest.json`` will be written.

    Returns
    -------
    int
        Exit code (0 = success).
    """
    import json
    import itertools

    make_case_fn, sweep_spec = load_sweep_spec(config_path)

    # Expand axes into a list of parameter dicts.
    # mode="grid" → Cartesian product; mode="zip" → paired iteration.
    axes = sweep_spec.axes
    if sweep_spec.mode == "zip":
        combos = list(zip(*[ax.values for ax in axes]))
        runs_params = [
            {ax.name: val for ax, val in zip(axes, combo)}
            for combo in combos
        ]
    else:
        # Default: full Cartesian product
        value_lists = [ax.values for ax in axes]
        runs_params = [
            {ax.name: val for ax, val in zip(axes, combo)}
            for combo in itertools.product(*value_lists)
        ]

    n_runs = len(runs_params)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    manifest: Dict = {
        "meta": {
            "config": str(config_path),
            "binary": str(binary_path),
            "mode": sweep_spec.mode,
            "reduce": sweep_spec.reduce,
            "n_runs": n_runs,
            "created": timestamp,
        },
        "runs": [
            {
                "idx": i,
                "tag": _make_run_tag(i, params),
                "params": params,
            }
            for i, params in enumerate(runs_params)
        ],
    }

    # Manifest lives alongside the config: <stem>_manifest.json
    # This matches the SLURM script's derivation: ${SWEEP_CONFIG%.py}_manifest.json
    manifest_path = config_path.with_name(config_path.stem + "_manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(
        f"[sweep] wrote manifest with {n_runs} runs → {manifest_path.name}",
        flush=True,
    )

    # Print per-axis summary
    for ax in axes:
        print(
            f"[sweep]   axis '{ax.name}': {len(ax.values)} values "
            f"({ax.values[0]} … {ax.values[-1]})",
            flush=True,
        )
    print(f"[sweep] total runs: {n_runs}", flush=True)

    # Use relative paths — submit from pipeline_dir so the SLURM script and
    # config are referenced by filename only, binary by ../sim/ relative path.
    sbatch_cmd = [
        "sbatch",
        f"--array=0-{n_runs - 1}",
        f"--export=SWEEP_CONFIG={config_path.name},SWEEP_BINARY=../sim/{binary_path.name}",
        "orchestrate_curved_grid_sweep.slurm",
    ]
    print(f"[sweep] submitting: {' '.join(sbatch_cmd)}", flush=True)
    result = subprocess.run(sbatch_cmd, capture_output=True, text=True, cwd=pipeline_dir)
    if result.returncode != 0:
        print(f"[sweep] sbatch error:\n{result.stderr}", flush=True)
        return result.returncode

    array_job_id = result.stdout.strip().split()[-1]
    print(f"[sweep] submitted array job {array_job_id}", flush=True)

    # Optionally submit a reduction job dependent on the array completing.
    if sweep_spec.reduce:
        reduce_script = pipeline_dir.parent / "tools" / "reduce_best.py"
        if not reduce_script.exists():
            # Fall back: search relative to pipeline_dir
            reduce_script = pipeline_dir / "reduce_best.py"

        if not reduce_script.exists():
            print(
                f"[sweep] WARNING: reduce_best.py not found — skipping "
                f"reduction job submission.",
                flush=True,
            )
        else:
            dep_cmd = [
                "sbatch",
                f"--dependency=afterok:{array_job_id}",
                "--export=ALL,"
                f"SWEEP_MANIFEST={manifest_path}",
                "--wrap",
                f"python {reduce_script} --manifest {manifest_path}",
            ]
            print(
                f"[sweep] submitting reduction job (dep=afterok:{array_job_id})",
                flush=True,
            )
            dep_result = subprocess.run(dep_cmd, capture_output=True, text=True)
            if dep_result.returncode != 0:
                print(
                    f"[sweep] reduction sbatch error:\n{dep_result.stderr}",
                    flush=True,
                )
            else:
                reduce_job_id = dep_result.stdout.strip().split()[-1]
                print(
                    f"[sweep] submitted reduction job {reduce_job_id}",
                    flush=True,
                )

    return 0


def run_sweep_task(
    config_path: Path,
    binary_path: Path,
    task_id: int,
    pipeline_dir: Path,
) -> int:
    """Run one task from a sweep SLURM array.

    Derives the parameter set for *task_id* directly from the sweep config's
    SWEEP spec — no manifest file required.  This ensures tasks always use
    the current config, eliminating stale-manifest bugs.

    Parameters
    ----------
    config_path:
        Path to the ``.py`` sweep config file (must define ``make_case`` and
        ``SWEEP``).
    binary_path:
        Path to the multi-grid solver binary.
    task_id:
        Zero-based combination index — typically ``$SLURM_ARRAY_TASK_ID``.
    pipeline_dir:
        Working directory passed to the solver subprocess.

    Returns
    -------
    int
        Exit code (0 = success).
    """
    make_case_fn, sweep_spec = load_sweep_spec(config_path)
    combos = sweep_spec.combinations()

    if task_id < 0 or task_id >= len(combos):
        raise IndexError(
            f"task_id {task_id} out of range — config defines "
            f"{len(combos)} combinations."
        )

    params: Dict = combos[task_id]
    tag = _make_run_tag(task_id, params)

    print(
        f"[sweep-task {task_id}] tag={tag}  params={params}",
        flush=True,
    )

    case = make_case_fn(**params)

    env = build_curved_env(case)
    env["RUN_TAG"] = tag
    if not env.get("RUN_STAMP"):
        env["RUN_STAMP"] = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    print(
        f"[sweep-task {task_id}] launching {binary_path} …", flush=True
    )
    subprocess.run([str(binary_path)], check=True, env=env, cwd=pipeline_dir)
    print(f"[sweep-task {task_id}] solver finished.", flush=True)

    # Auto-generate sample diameter profile plot if enabled.
    if int(env.get("WRITE_PNG", "1")):
        _run_profile_plot(env, pipeline_dir)

    return 0


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Run the curved-grid IBSimu pipeline."
    )
    parser.add_argument(
        "--config",
        default="curved_grid_case_example.py",
        help="Python config file defining CASE = CurvedSimulationCase(...)",
    )
    parser.add_argument(
        "--binary",
        default="multi_grid_2d",
        help="Path to the multi-grid solver binary",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the generated GRID_STACK and aperture geometry, then exit",
    )
    parser.add_argument(
        "--plot",
        metavar="FILE",
        default=None,
        help="Save a cross-section profile PNG to FILE and exit",
    )
    parser.add_argument(
        "--sweep",
        action="store_true",
        help=(
            "Generate sweep_manifest.json and submit a SLURM array job. "
            "Requires the config file to define make_case() and SWEEP."
        ),
    )
    parser.add_argument(
        "--sweep-task",
        metavar="INT",
        type=int,
        default=None,
        help=(
            "Run a single sweep task by its zero-based index. "
            "The manifest path is taken from --manifest or the SWEEP_MANIFEST "
            "environment variable."
        ),
    )
    parser.add_argument(
        "--manifest",
        metavar="PATH",
        default=None,
        help=(
            "Path to sweep_manifest.json used by --sweep-task. "
            "Overrides the SWEEP_MANIFEST environment variable."
        ),
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = repo_root / config_path
    binary_path = Path(args.binary)
    if not binary_path.is_absolute():
        binary_path = repo_root / binary_path

    # ── Sweep dispatch ──────────────────────────────────────────────────────
    if args.sweep:
        return run_sweep(config_path, binary_path, repo_root)

    if args.sweep_task is not None:
        return run_sweep_task(
            config_path, binary_path, args.sweep_task, repo_root
        )

    case = load_curved_case(config_path)

    if args.dry_run:
        stack = serialize_curved_grid_stack(case.grids)
        print("GRID_STACK =", stack)
        print()
        for grid in case.grids:
            positions = compute_curved_aperture_positions(
                grid.apertures, grid.curvature
            )
            print(f"Grid '{grid.name}'  R={grid.curvature.radius_m} m  "
                  f"concave_upstream={grid.curvature.concave_upstream}")
            for cp in positions:
                print(
                    f"  r={cp.offset_radial_m*1e3:+7.3f} mm  "
                    f"dz={cp.offset_axial_m*1e3:+6.3f} mm  "
                    f"tilt={cp.normal_tilt_deg:+7.3f} deg  "
                    f"bore_r={cp.aperture_radius_m*1e3:.3f} mm"
                )
        return 0

    if args.plot:
        fig = plot_curved_grid_profile(
            case.grids,
            title=f"Curved grid profile — {config_path.stem}",
        )
        out = Path(args.plot)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Profile saved to {out}")
        return 0

    if not binary_path.exists():
        raise FileNotFoundError(
            f"multi-grid solver not found: {binary_path}"
        )

    env = build_curved_env(case)

    # Pin RUN_STAMP now so Python and C++ agree on the output directory name.
    if not env.get("RUN_STAMP"):
        env["RUN_STAMP"] = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    stack = env.get("GRID_STACK", "")
    print(f"[curved-grid] binary    : {binary_path}", flush=True)
    print(f"[curved-grid] config    : {config_path}", flush=True)
    print(f"[curved-grid] cwd       : {repo_root}", flush=True)
    print(f"[curved-grid] GRID_STACK: {stack[:120]}{'...' if len(stack) > 120 else ''}", flush=True)
    print(f"[curved-grid] RUN_SOLVE : {env.get('RUN_SOLVE', '(not set)')}", flush=True)
    print(f"[curved-grid] YBOX_M    : {env.get('YBOX_M', '(not set)')}", flush=True)
    print(f"[curved-grid] launching solver...", flush=True)
    subprocess.run([str(binary_path)], check=True, cwd=repo_root, env=env)
    print(f"[curved-grid] solver finished.", flush=True)

    # Auto-generate sample diameter profile plot if WRITE_PNG is active.
    if int(env.get("WRITE_PNG", "1")):
        _run_profile_plot(env, repo_root)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
