#!/usr/bin/env bash
# Parameters for the axisymmetric meniscus micro-model (sim/meniscus_cell).
# Keep species/voltages identical to the 3-D run so the beamlet is consistent.
export VS_V=0.0                   # screen / plasma potential
export VA_V=-10000.0              # accel (-10 kV)  -> 10 kV bias
export PLASMA_NI_M3=1.5e18        # DESIGN ion density (matches the 3-D run); sets Bohm supply
export PLASMA_TE_EV=3.7
export ION_Q_E=1.0                # q = +1
export ION_M_AMU=2.0             # m = 2 amu
export ION_TP_EV=0.0
export ION_TT_EV=0.2
export ION_J_SCALE=1.0

# single-bore geometry (SI metres) -- set to match the real grid bore + gap
export MEN_H=5.0e-6               # mesh 5 um; Debye ~12 um at 1.5e18/3.7eV -> ~2.4 cells/lambda_D (resolved)
export MEN_T_SCR=6.0e-3
export MEN_T_ACC=6.0e-3
export MEN_GAP=4.0e-3
export MEN_RBORE_SCR=1.5e-3       # CONFIRM against the CAD bore radius
export MEN_RBORE_ACC=1.5e-3
export MEN_RMAX=4.0e-3
export MEN_NPART=4000
export MEN_ITER_MAX=30            # headroom: the stiffer 1.5e18 sheath needs more iters to settle (single bore is cheap)
export MEN_RESULTS_DIR=results_men     # aux outputs (men_geom.dat, men_epot.dat)

# beamlet output (shared with the 3-D solver) + reuse cache
export MENISCUS_FILE=meniscus_cache/meniscus_beamlet.dat
export MENISCUS_REUSE=0                 # 0 = FORCE re-solve (do not trust the stale cache); 1 = reuse if present
