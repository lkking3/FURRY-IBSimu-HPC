# 3-D Curved Two-Grid Extraction Module (`grid_3d_curved`)

Full 3-D, self-consistent IBSimu model of the curved two-grid ion-optics
system, built from the SolidWorks CAD grids. Companion to the validated 2-D
curved module (`sim/multi_grid_2d_curved.cpp`); physics and species are carried
over for continuity and cross-checking.

## What this models
- **Geometry**: screen (`V_s = 0`) and accel (`V_a = -10 kV`) spherical-cap grids
  imported as STL solids, assembled concentric with a **4 mm face-to-face gap**
  (10 mm midplane-to-midplane), optical axis = `+z`, screen upstream.
  Screen R_c = 200 mm, accel R_c = 190 mm; ~3 mm bores.
- **Symmetry**: a **HALF** model with a Cartesian mirror on the `x = 0` plane.
  The aperture pattern is hexagonal (6-fold); its mirror planes lie at
  30 deg/90 deg/150 deg, so **no two are orthogonal** and a quarter model is
  *not* valid on IBSimu's Cartesian mesh. The `x = 0` plane is a verified clean
  mirror (reflection error p99 ~ 0.26 mm). Half the disk -> factor 2.
- **Physics**: self-consistent (Vlasov) space charge. Ar+ ions emitted at the
  Bohm velocity from each of 241 apertures along its local surface normal
  (`add_cylindrical_beam_with_energy`); Boltzmann electrons via the
  positive-exponential plasma model (`set_pexp_plasma`), as in 2-D.

## Files
```
sim/grid_3d_curved.cpp            # the module (MODE_3D + STLSolid)
sim/Makefile.3d                   # build:  make -f Makefile.3d
curved_grid_3d/
  case_baseline.sh                # carried-over physics + domain (env vars)
  run_grid_3d.slurm               # SLURM launcher (1 CPU, 48 GB)
  geometry/
    screen_half_pos.stl           # positioned half screen (mm)
    accel_half_pos.stl            # positioned half accel  (mm)
    apertures_screen.dat          # 241 emission sites (x y z nx ny nz r, metres)
tools/extract_apertures_3d.py     # regenerates the aperture file from an STL
```

## Build & run (on the cluster, via SLURM only)
```
cd sim && make -f Makefile.3d            # needs pkg-config ibsimu-1.0.6dev
cd ../curved_grid_3d
sbatch run_grid_3d.slurm case_baseline.sh
```
Outputs land in `results_3d/`: `geom.dat`, `epot.dat`, `pdb.dat` (for the IBSimu
analysis/plotter), and `diagnostics.json` (merged-beam half-angle, currents).

## Domain & cost (h = 0.2 mm)
476 x 951 x 1068 = **483 M nodes ~ 39 GB**. The SLURM script requests 48 GB on
1 CPU (IBSimu's Poisson solve is single-threaded). For a cheap pipeline shake-out
set `H=4.0e-4` (~6 GB) in the case file first.

## Verification plan (do before trusting results)
1. **Compile/link on the cluster.** IBSimu is absent from the prep sandbox; the
   module follows the IBSimu 1.0.6 RADIS STL example + the 2-D module. Reconcile
   any diagnostics-block signature drift (`trajectories_at_plane`) if it warns.
2. **Coarse shake-out** (`H=4e-4`, `ION_ITER_MAX=4`): confirm it solves, particles
   transmit through bores, `diagnostics.json` is written.
3. **2-D vs 3-D on-axis beamlet**: a single central aperture should reproduce the
   2-D curved result for the same bore/voltages.
4. **Mesh convergence**: rerun the merged half-angle at `H=0.3 mm` and `0.2 mm`;
   the 0.2 mm value is the production estimate of discretization error.
5. **Symmetry sanity**: the half-model field/trajectories must be mirror-flat at
   `x = 0` (no kink across the plane).

## Known limitations / refinements
- **Flat plasma boundary on a curved grid**: `InitialPlasma(AXIS_Z, Z_PLASMA)`
  is a flat plane; for the dished screen it conforms only near the apex
  (same approximation as the 2-D module). A conformal/meniscus plasma surface is
  the eventual upgrade (the `set_pexp_plasma` meniscus you flagged for later).
- **Resolution**: h = 0.2 mm resolves the 3 mm bore by ~15 cells -- coarser than
  the 2-D 0.1 mm; quantify via step 4.
- **Aperture set** is extracted from the CAD by ray-casting (241/half, ~2.6 mm
  effective emission radius). Regenerate with `extract_apertures_3d.py` if the
  CAD changes; cross-check the count against the SolidWorks layout.
- **STL closedness**: both half STLs have 0 open edges (closed); trimesh reports
  `watertight=False` only due to non-manifold winding around the many bores --
  not a leak. Re-validate if the STLs are re-exported.
