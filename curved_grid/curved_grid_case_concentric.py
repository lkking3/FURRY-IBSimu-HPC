"""
curved_grid_case_example.py

Two-grid ion accelerator with spherically curved electrodes using a
CONCENTRIC SPHERE design: both grids share a single center of curvature so
that every radial line from that centre passes through the aligned screen and
accel aperture centres.

Concentric geometry
-------------------
Screen:  vertex at x = 0,       radius R_scr = 100 mm → centre at x = +100 mm
Accel:   vertex at x = Δx,      radius R_acc = R_scr − Δx → same centre x = +100 mm

where Δx = screen_thickness + effective_gap
         = 5 mm + (6 mm − sagitta(25 mm, 100 mm))
         = 5 mm + (6 mm − 3.18 mm)
         = 5 mm + 2.82 mm  ≈ 7.82 mm

→ R_acc ≈ 92.18 mm, scale = R_acc/R_scr ≈ 0.9218

Aperture layout (screen → accel, same radial angle α)
───────────────────────────────────────────────────────
r_scr = 0 mm    → r_acc = 0 mm       α = 0°
r_scr = 15 mm   → r_acc ≈ 13.83 mm  α = 8.63°
r_scr = 25 mm   → r_acc ≈ 23.04 mm  α = 14.48°

Bore clearance check (normal-aligned bore, 5 mm thick grids, α = 14.48°):
  bore_displacement = 5 × tan(14.48°) = 1.29 mm  <  1.5 mm bore radius  ✓
"""
import math

from curved_grid_pipeline import (
    Aperture,
    Chamfer,
    BeamletPair,
    make_apertures_from_pairs,
    concentric_accel_radius,
    concentric_accel_offsets,
    GridCurvature,
    CurvedGridDefinition,
    CurvedSimulationCase,
)

# ---- geometry parameters ----
R_SCR_M      = 0.100   # screen radius of curvature [m]
AP_RADIUS_M  = 0.0015  # bore radius [m]
SCREEN_T_M   = 0.005   # screen thickness [m]
GAP_AFTER_M  = 0.006   # nominal gap after screen vertex [m]

SCREEN_OFFSETS_M = [0.000, 0.015, 0.025]

# Compute accel radius and offsets for concentric-sphere bore alignment.
# concentric_accel_radius accounts for the sagitta-adjusted effective gap,
# matching exactly what serialize_curved_grid_stack will apply.
R_ACC_M = concentric_accel_radius(
    R_scr_m=R_SCR_M,
    screen_thickness_m=SCREEN_T_M,
    gap_after_m=GAP_AFTER_M,
    screen_offsets_m=SCREEN_OFFSETS_M,
)
ACCEL_OFFSETS_M = concentric_accel_offsets(SCREEN_OFFSETS_M, R_SCR_M, R_ACC_M)

# Build paired aperture lists.  mirror=False because the grids use mirror=True.
PAIRS = [
    BeamletPair(
        screen_offset_m=s,
        accel_offset_m=a,
        screen_radius_m=AP_RADIUS_M,
    )
    for s, a in zip(SCREEN_OFFSETS_M, ACCEL_OFFSETS_M)
]
screen_apertures, accel_apertures = make_apertures_from_pairs(PAIRS, mirror=False)

CASE = CurvedSimulationCase(
    grids=[
        CurvedGridDefinition(
            name="screen",
            voltage_v=0.0,
            thickness_m=SCREEN_T_M,
            gap_after_m=GAP_AFTER_M,
            apertures=screen_apertures,
            upstream_chamfer=Chamfer(0.0, 0.0),
            downstream_chamfer=Chamfer(0.0, 0.0),
            mirror=True,
            curvature=GridCurvature(radius_m=R_SCR_M, concave_upstream=True),
        ),
        CurvedGridDefinition(
            name="accel",
            voltage_v=-15000.0,
            thickness_m=0.005,
            gap_after_m=0.0,
            apertures=accel_apertures,
            upstream_chamfer=Chamfer(0.001, 5.0),
            downstream_chamfer=Chamfer(0.0013, 5.0),
            mirror=True,
            # Accel uses R_ACC_M so its bore directions also point to the
            # shared centre of curvature — same α as the screen apertures.
            curvature=GridCurvature(radius_m=R_ACC_M, concave_upstream=True),
        ),
    ],
    env={
        "RUN_SOLVE":      "1",
        "WRITE_PNG":      "1",
        "PNG_NAME":       "curved_grid_geom.png",
        "RESULTS_DIR":    "results",
        "RUN_PREFIX":     "curved_grid_test",
        "ENABLE_IONS":    "1",
        "X_RIGHT_M":      "0.025",
        "X_RIGHT_PHYS_M": "0.55",
        "YBOX_M":         "0.090",
        "TUBE_WALL_T_M":  "0.0002",
        "H":              "0.0001",
        "SC_FACTOR":      "0.0005",
    },
)
