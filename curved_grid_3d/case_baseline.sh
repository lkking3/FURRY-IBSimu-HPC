#!/usr/bin/env bash
# Baseline 3-D curved two-grid extraction case.
# Physics carried over from the validated 2-D curved module; geometry from the
# assembled half-stack STLs (4 mm gap, concentric, optical axis +z).
# Source this before launching the binary, or let run_grid_3d.slurm read it.

# ---- geometry / domain (SI metres) ----
export H=3.0e-4                 # cell size (0.3 mm) -> ~114 M nodes, ~0.8 B matrix nonzeros.
                                # DO NOT drop to 0.2 mm on this build: ~481 M nodes -> ~3.4 B
                                # nonzeros overflows IBSimu's 32-bit CRowMatrix index -> SIGSEGV
                                # in construct_add(). 0.3 mm is the validated resolution.
export LX=0.095                 # half radial extent (x>=0, mirror at x=0)
export LY=0.095                 # +/- transverse extent
export LZ=0.2134                # axial length (stack + 170 mm drift)
export Z_PLASMA=0.007           # plasma fills z < this (just upstream of screen)

# Full-grid STLs (clean 2-manifold). The HALF model is realised by the x>=0
# mesh domain + mirror BC at x=0, not by cutting the STL.
export SCREEN_STL=geometry/screen_pos.stl
export ACCEL_STL=geometry/accel_pos.stl
export APERTURE_FILE=geometry/apertures_screen.dat

# ---- electrodes ----
export VS_V=0.0                 # screen / plasma potential
export VA_V=-10000.0            # accel  (-10 kV)

# ---- plasma + ion species ----
# DESIGN operating point: n_i = 1.5e18, q = +1, m = 2 amu, 10 kV bias.
#   Bohm-over-bore (r=1.5 mm) at this density ~= 13.6 mA/bore; the meniscus
#   micro-model sets the actual per-bore current (see case_meniscus_cell.sh).
#   I scales linearly with PLASMA_NI_M3 and 1/sqrt(ION_M_AMU).
export PLASMA_NI_M3=1.5e18
export PLASMA_TE_EV=3.7
export ION_Q_E=1.0             # q = +1
export ION_M_AMU=2.0           # m = 2 amu
export ION_TP_EV=0.0
export ION_TT_EV=0.2
export ION_J_SCALE=1.0

# ---- injection mode ----
#   bohm     : Bohm flux emitted at each aperture (this model does extraction).
#   meniscus : inject the pre-extracted beamlet from meniscus_cell (run that
#              first to make MENISCUS_FILE); extraction physics + current come
#              from the fine micro-model, this model does drift + space charge.
export INJECTION_MODE=${INJECTION_MODE:-meniscus}   # design workflow: inject the vetted meniscus beamlet
export MENISCUS_FILE=${MENISCUS_FILE:-meniscus_cache/meniscus_beamlet.dat}
export MENISCUS_L_EXTRACT=0.017   # screen surface -> accel exit, along the normal [m]
export MENISCUS_NPER=500          # macroparticles injected per aperture (meniscus mode)

# ---- solver / iteration ----
export ION_NPART=60000          # total trajectories (~250 per aperture)
export ION_ITER_MAX=25          # ceiling; early-stop ends sooner when converged
export ION_MAX_STEPS=4000
export SC_ALPHA=0.25            # under-relaxation; 0.25 damps the high-current meniscus oscillation (0.7 sawtooths)
export CONV_TOL=0.01            # early stop when current AND half-angle change < 1%/iter
export CONV_MIN_ITER=4          # but run at least this many iterations first
export ION_THREADS=${ION_THREADS:-16}   # IBSimu threads mesh build + solver (match SLURM cpus)

# ---- geometry cache: build the mesh once (slow), reuse it across runs/sweep ----
# Name is tagged by H AND the domain dims (LX/LY/LZ) so a mesh built for a shorter
# domain can never be silently reused for a longer one -- that mismatch put the
# diagnostic/target plane outside the mesh and gave n@plane=0 / transmission=0%.
# Still delete the file if you change the STLs or VS/VA (those aren't in the name).
export GEOM_CACHE=${GEOM_CACHE:-geom_cache/geom_h${H}_x${LX}_y${LY}_z${LZ}.dat}

# ---- drift space-charge compensation (physical; carried over from 2-D) ----
export SC_FACTOR=${SC_FACTOR:-0.0005}        # 1=full SC, ->0=fully neutralised drift
export SC_RAMP_START_Z=${SC_RAMP_START_Z:-0.050}  # start of compensation (m, past accel)
export SC_RAMP_LEN_Z=${SC_RAMP_LEN_Z:-0.0}        # ramp length (m); 0 = step

# ---- visualisation outputs ----
export WRITE_PNG=1               # GeomPlotter cross-section PNGs (headless)
export WRITE_VTK=1               # epot.vtk + beam_density.vtk + beam_trajectories.vtk
export VTK_MIRROR=1             # reflect the x>=0 half across x=0 -> full beam in the VTK
export VTK_STRIDE=2             # field-volume downsample (2 -> 1/8 nodes); raise for full run
export TRAJ_STRIDE=4            # write every Nth beamlet line to beam_trajectories.vtk
export WRITE_ENVELOPE=1         # envelope.csv: beam radius + divergence vs z (waist finder)
export ENV_NZ=60                # number of z-planes in the envelope scan

# ---- diagnostic planes ----
# PIN the target to the physical 2.5" target plane. Do NOT leave this unset: the
# code default is z_diag (= LZ-5h), which now sits at the domain end (~0.2119 m),
# not the target -- that would measure the on-target current 42 mm too far
# downstream and understate it.
export Z_TARGET=0.17            # physical 2.5" target plane [m]

export RESULTS_DIR=${RESULTS_DIR:-results_3d}
