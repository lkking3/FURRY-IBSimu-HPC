# Upgrade Plan: Self-Consistent Meniscus → Informed 3-D Injection

Status: design / not yet implemented. Author note: companion to `sim/grid_3d_curved.cpp`.

## 1. Problem with the present injection

The 3-D model currently injects ions from each aperture with a **fixed Bohm
current density** over a flat disk (`add_cylindrical_beam_with_energy`), 2h
upstream of the screen, along the local surface normal, with a flat
`InitialPlasma(AXIS_Z, Z_PLASMA)` plane. Consequences:

- The plasma **meniscus** (the curved sheath edge where ions are born) is never
  solved — beamlet birth angle, emittance and current-density profile are
  imposed, not physical.
- The total current is **hand-tuned** via `PLASMA_NI_M3` (the 20 A → 4 A fix)
  instead of being self-limited by the extraction physics.
- The flat plasma plane does not conform to the **dished** grid; peripheral
  apertures (screen surface at z≈33 mm) are not even in contact with the flat
  plasma region (z < ~7 mm).

## 2. Why it must be an external (multi-scale) solve

Debye length at n_i = 4e17 m⁻³, Te = 3.7 eV:  λ_D ≈ **23 µm**.
The meniscus sheath is a few λ_D (≈ 50–100 µm). The full-grid mesh is
**200–300 µm** — and cannot be refined locally (IBSimu uses a uniform Cartesian
mesh). Resolving the meniscus in the full model would need a µm-scale mesh over
the whole 190 mm disk: physically impossible (~10¹²–10¹³ nodes).

Therefore the meniscus must be solved in a **small, high-resolution micro-model**
of a single aperture, and its emitted beamlet used as the injection boundary
condition for the full transport model. This is the standard scale-separation
approach for ion-optics extraction.

## 3. Key simplification: one meniscus serves all apertures

The grids are **concentric with a uniform 4 mm gap and fixed voltages**, so the
*local* extraction field at every bore is essentially identical — only the
beamlet **aim direction** (tilt toward the common centre of curvature) varies
with radius. So:

- the **meniscus shape and emitted phase space** can come from a *single*
  high-res axisymmetric solve;
- the **per-aperture variation** is a geometric transform (position + tilt +
  azimuth) applied when placing that beamlet in 3-D — which we already compute
  (the aperture normals in `apertures_screen.dat`).

This makes the MVP cheap: one micro-solve, reused 210×.

## 4. Architecture

```
[A] Meniscus micro-model            [B] 3-D transport model (grid_3d_curved)
  single aperture, MODE_CYL            full half-grid, MODE_3D
  µm mesh, self-consistent plasma  →   loads beamlet phase space,
  -> emitted beamlet phase space       places it per-aperture (pos/tilt/azimuth),
     (r, r', I, E) at bore exit        transports with mutual space charge
```

### Component A — meniscus / extraction micro-model (`sim/meniscus_cell.cpp`)
- **Geometry:** one bore in the local screen+accel gap (4 mm, VS/VA). Axisymmetric
  (`MODE_CYL`, r–z) for the meniscus; cheap (r ≤ ~4 mm, z ≤ ~16 mm).
- **Mesh:** ~5 µm near the meniscus (resolves λ_D ≈ 23 µm comfortably). A few
  10⁶ nodes — seconds-to-minutes.
- **Physics:** `set_initial_plasma` + `set_pexp_plasma` (already used in 2-D) with
  Bohm injection from a generous upstream plasma volume; iterate to convergence.
  The meniscus forms self-consistently and the extracted current is
  **space-charge/Bohm self-limited** → physical magnitude, no n_i fudge.
- **Output:** phase space at the bore-exit plane — `(r, r', I_weight, E)` for the
  axisymmetric case (revolve to (x,y,vx,vy,vz) when placing in 3-D).
- **Reuse:** the existing 2-D curved solver already does self-consistent plasma
  extraction; adapt it (MODE_CYL, single bore, fine mesh, phase-space dump).

### Component B — 3-D injection from meniscus data (edit `grid_3d_curved.cpp`)
- Add `INJECTION_MODE = bohm | meniscus` (env). `bohm` = current behaviour.
- `meniscus`: load the micro-model phase space once; for each aperture:
  1. sample N macroparticles from the (r, r', I) distribution (revolve in azimuth);
  2. transform to the aperture's 3-D centre, rotate so the beamlet axis = the
     aperture's inward normal (the tilt-toward-focus), random azimuth;
  3. inject via `pdb.add_particle(Particle3D)` with per-particle current weight.
- Total current = sum of per-beamlet currents from the micro-model → physical,
  not set by `PLASMA_NI_M3`.

## 5. Data interface (micro-model → 3-D)
Plain-text `meniscus_beamlet.dat`:
```
# E_eV  I_total_A_per_bore  n_points
# r_m   rp_rad   weight     (one row per phase-space sample)
...
```
Loader in the 3-D module mirrors `read_apertures`. Keep it text + simple so the
micro-model and transport model stay decoupled.

## 6. Coupling levels (phasing)
- **Phase 1 (MVP, one-way):** single axisymmetric meniscus → injected at all
  apertures with geometric tilt. Captures real beamlet emittance + self-limited
  current. Assumes the local field is the same at every bore (true to 1st order).
- **Phase 2 (radial stations):** if the gap/field actually varies with radius
  (check against the assembled geometry), run the micro-model at 3–5 radial
  stations and interpolate by aperture radius.
- **Phase 3 (co-iterated, the "simultaneous" version):** after a 3-D solve, sample
  the *actual* local field along each bore axis (includes neighbour space charge),
  feed it back as the micro-model's boundary field, re-solve the meniscus, update
  injection, repeat until the emitted distribution stabilises. Only needed if
  neighbour space charge perturbs the local extraction appreciably (likely small,
  since the meniscus is dominated by the 4 mm gap field).

## 7. Concrete task list
1. `sim/meniscus_cell.cpp` — MODE_CYL single-bore self-consistent extraction;
   dump `meniscus_beamlet.dat`. (Adapt from 2-D curved solver.) + Makefile target.
2. Driver to run it at the design point (and later a radial sweep).
3. `grid_3d_curved.cpp`: add `INJECTION_MODE`, a meniscus loader, the
   per-aperture sample+transform+`add_particle` injector.
4. Aperture→beamlet placement: reuse `apertures_screen.dat` normals; add azimuthal
   sampling and (Phase 2) radial interpolation.
5. Case-file knobs: `INJECTION_MODE`, `MENISCUS_FILE`, particles-per-aperture.

## 8. Validation
- Micro-model vs the existing 2-D curved single-aperture result (same gap/V): same
  meniscus shape, perveance, beamlet divergence.
- Debye-resolution check: halve the micro-model mesh, confirm meniscus unchanged.
- 3-D `meniscus` vs `bohm` injection: compare total current, merged half-angle,
  focal spot. Expect the self-limited current to land physically without the n_i
  fudge, and a tighter/realistic emittance.
- Confirm total current is now insensitive to `PLASMA_NI_M3` above the
  space-charge-limited regime (the meniscus self-limits).

## 9. Risks / open questions
- **Axisymmetry off-axis:** tilted bores at large radius aren't perfectly
  axisymmetric; Phase 1 treats meniscus as axisym + rigid tilt. If edge beamlets
  matter, a small 3-D micro-box at the edge station validates the approximation.
- **`add_particle` throughput:** 210 apertures × (say) 500 macroparticles = ~10⁵
  particles; fine, but check injection time and per-particle current weighting API.
- **Number of radial stations** for Phase 2 convergence.
- **Co-iteration stability/cost** (Phase 3) — only pursue if Phase 1/2 show the
  local field is materially perturbed by neighbours.
- The micro-model energy/zero-point must match the 3-D model's reference so the
  beamlet is injected at the correct potential/energy.

---

## 10. IMPLEMENTED (Phase 1) — build & run order

Files added/changed:
- `sim/meniscus_cell.cpp` — axisymmetric (MODE_CYL) micro-model + `Makefile.3d` target.
- `sim/grid_3d_curved.cpp` — `INJECTION_MODE=bohm|meniscus`; meniscus loader
  (`read_beamlet`); per-aperture `add_particle` injection at the accel exit
  (placement = aperture surface + normal*`MENISCUS_L_EXTRACT`, axisymmetric
  beamlet revolved with the golden angle); in meniscus mode the plasma model is
  off (drift transport only); injected-current + transmission diagnostics.
- `curved_grid_3d/case_meniscus_cell.sh`, `run_meniscus.slurm` — micro-model run.
- `curved_grid_3d/case_baseline.sh` — `INJECTION_MODE` + `MENISCUS_*` knobs.

Run order on the cluster:
```
cd sim && make -f Makefile.3d meniscus_cell grid_3d_curved
cd ../curved_grid_3d
# 1) solve the meniscus once (cheap) -> results_men/meniscus_beamlet.dat
sbatch run_meniscus.slurm case_meniscus_cell.sh
# 2) run the 3-D transport in meniscus mode (reuses the geometry cache)
cp case_baseline.sh case_men.sh
sed -i 's/^export INJECTION_MODE=.*/export INJECTION_MODE=meniscus/' case_men.sh
sbatch run_grid_3d.slurm case_men.sh
```

Verification when it lands (what to check):
- `results_men/...`: meniscus printed `I_per_bore`, `r_rms`, `rp_rms`, `emittance`,
  and converged (dI<0.5%). Halve `MEN_H` once to confirm the meniscus is mesh-
  converged.
- 3-D `diagnostics.json`: `transmission_pct` sensible (high), `I_full_grid_A` ~=
  2 x I_per_bore x apertures x transmission, and now *self-limited* (insensitive
  to PLASMA_NI_M3 in the space-charge-limited regime).
- Compare meniscus-mode vs bohm-mode `merged_half_angle_deg` / envelope: expect a
  more physical emittance and a tighter, current-consistent beam.

## 11. Compile-time items to reconcile on the cluster (IBSimu not in prep env)
- `ParticleDataBaseCyl`, `MODE_CYL`, `add_2d_beam_with_energy` (cyl) — mirrored
  from two_grid_2d.cpp; should match.
- `pdb.add_particle( IQ, q, m, ParticleP3D(t,x,vx,y,vy,z,vz) )` — verify the exact
  signature/argument order against the installed `particledatabase.hpp`.
- `FuncSolid(&fn)` with `bool fn(double,double,double)` — matches the 2-D usage.
