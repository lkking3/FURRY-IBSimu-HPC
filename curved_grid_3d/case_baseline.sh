#!/usr/bin/env bash
# Baseline 3-D curved two-grid extraction case.
# Physics carried over from the validated 2-D curved module; geometry from the
# assembled half-stack STLs (4 mm gap, concentric, optical axis +z).
# Source this before launching the binary, or let run_grid_3d.slurm read it.

# ---- geometry / domain (SI metres) ----
export H=2.0e-4                 # cell size (0.2 mm)  -> ~481 M nodes, ~39 GB
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

# ---- plasma + ion species (Ar+, carried over) ----
export PLASMA_NI_M3=1.0e16
export PLASMA_TE_EV=3.7
export ION_Q_E=1.0
export ION_M_AMU=40.0           # Argon
export ION_TP_EV=0.0
export ION_TT_EV=0.2
export ION_J_SCALE=1.0

# ---- solver / iteration ----
export ION_NPART=60000          # total trajectories (~250 per aperture)
export ION_ITER_MAX=15          # self-consistent Vlasov iterations (divergence settles slowly)
export ION_MAX_STEPS=4000
export SC_ALPHA=0.5             # space-charge averaging factor
export ION_THREADS=1            # IBSimu Poisson solve is single-core

# ---- visualisation outputs ----
export WRITE_PNG=1               # GeomPlotter cross-section PNGs (headless)
export WRITE_VTK=1               # epot.vtk + beam_density.vtk + beam_trajectories.vtk
export VTK_MIRROR=1             # reflect the x>=0 half across x=0 -> full beam in the VTK
export VTK_STRIDE=2             # field-volume downsample (2 -> 1/8 nodes); raise for full run
export TRAJ_STRIDE=4            # write every Nth beamlet line to beam_trajectories.vtk
export WRITE_ENVELOPE=1         # envelope.csv: beam radius + divergence vs z (waist finder)
export ENV_NZ=60                # number of z-planes in the envelope scan

export RESULTS_DIR=results_3d
