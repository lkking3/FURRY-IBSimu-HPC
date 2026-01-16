# FURRY-IBSimu-HPC

HPC SLURM submitted workflow for sweeping and optimizing a two-grid electrostatic ion accelerator model using **IBSimu** on **Omega cluster @ DIII-D**.

This repo is built around a single orchestrator job that:
1) Builds the solver.
2) Runs a parameter sweep as a SLURM array.
3) Outputs data to results folder, itemized by job name, certain key geometric parameters, and rundate.
   3.1) Data is delivered in the form of .json folders and .png images
   3.2) Pertinent simulation data, external calculations, and title information are stored in meta.json
   3.3) Beam informatics and statistics are stored in Beam_Metrics.json
   3.4) Information about simulation runtime & breakdown of highest utilization is presented in timing.json
   3.5) Predictions about beam propagation down a tube of length "X_RIGHT_PHYS_M" (default 0.55 m) can be found in Beam_Envelope.json
5) Update_runlog_compact.py compiles individual result files into a single master CSV file. Default conversion does not include all information from .json files for readability.
6) Plot_runlog_compact.py can be run by CLI to generate plotted information from the compact csv for both scatter plots and heatmaps (Further utility is yet to be implemented).

> !!!! **Do not run parameter sweeps on the login node.**  !!!!
> Use `sbatch` and compute partitions only.

---

## Contents

- `orchestrate_pipeline.slurm` — **main entrypoint**; prepares workspace, submits sweep array.
- `two_grid_2d.cpp` — IBSimu C++ simulation (compiled to `two_grid_2d`),  processes .json and .png output.
- `reduce_best.py` (or embedded reducer) — selects top candidates from `results/`. (**Depreciated**, still has functional architecture but needs a new heuristic by which to choose the best output.)
- `update_runlog_compact.py` — optional tooling to build/refresh a master CSV across runs.
- `Makefile` — builds `two_grid_2d`.
- `plot_runlog_compact.py'
Your run outputs will live under:
- `results/<RUN_PREFIX>_<stamp>_<tag>_j<jobid>/...`

## Simulation Control Variables

Below are the “big levers” you’ll most likely touch when doing parameter studies.

1) Sweep resources / scheduling
These control cost and throughput (per array element):
   **SWEEP_CPUS**, **SWEEP_MEM**, **SWEEP_TIME**, **SWEEP_PART**
   **ARRAY_CONC** (how many array tasks run at once)
   **MAX_ARRAY_CHUNK** (array chunking limit for scheduler friendliness)

3) Results naming / location
   **RESULTS_DIR** (defaults to ./results)
   **RUN_PREFIX** (prefix for run directories)

4) Core geometry sweep lists (Stage-1)
These define the parameter sweep grid:
   **AP_RAD_LIST** — aperture radius
   **GRID_T_LIST** — grid thickness
   **GAP_LIST** — grid gap
Chamfers (usually left at 0 unless you’re explicitly studying them):
   **SCR_*** and **ACC_*** lists: depths + angles

5) Electrical potentials
   **VS_V** — screen potential (acceleration voltage reference)
   **VA_V** — accel/tube potential
   **SAMPLE_V** — sample/endplate potential (often equals VA_V)

6) Physics toggles (field-only vs ions)
   **ENABLE_IONS_STAGE1** — set to 0 for fast field-only scans, 1 for plasma/beam iterations
   **ENABLE_IONS_STAGE2** — usually 1 if Stage-2 is meant to evaluate real beam performance
Plasma / source controls:
   **PLASMA_NI_M3** (ion density)
   **PLASMA_TE_EV** (electron temperature)
   **ION_J_SCALE** (current density scale factor)
   **SC_FACTOR** (space charge scaling / neutralization proxy)

7) Drift truncation + footprint prediction (compute saver)
If you truncate the drift (short X_RIGHT_M) but want performance at a longer physical drift:
   **X_RIGHT_M** — simulated drift length
   **X_RIGHT_PHYS_M** (or equivalent in your build) — “real” drift length for prediction
The code should switch the far boundary to Neumann and disable the physical sample plate when truncated, while still writing envelope_pred.json.

## Requirements (Omega)

This workflow assumes Omega module environment.

Typical modules:
```bash
module load gcc/11.x
module load libibsimu/1.0.6
**Libibsimu might not be natively installed on your HPC, in my case I had to install css-omega-modules to my home directory and then type the command: "module use ~/css-omega-modules" before loading libibsimu.

Credit for the creation of IBSimu goes to T. Kalvas et. al.
T. Kalvas, et. al., "IBSIMU: A three-dimensional simulation software for charged particle optics", Rev. Sci. Instrum. 81, 02B703, (2010). [Paper-PDF]
https://ibsimu.sourceforge.net/publications.html
