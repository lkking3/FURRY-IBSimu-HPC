#!/usr/bin/env python3
"""
layout_optimizer.py

Iterative design optimization loop for beamlet_layout_solver.

Provides three building blocks:
  - IterationRecord / IterationLog  : persistent per-iteration audit trail
  - SurrogateModel                  : bootstrap PchipInterpolator on the safe
                                      steering curve with uncertainty estimates
  - LayoutOptimizer                 : orchestrates record_solution(),
                                      ingest_results(), generate_slurm_sweep(),
                                      and get_convergence_data()

Typical workflow
----------------
1. GUI calls optimizer.record_solution(sol, ap_rad_m, gap_m)
   after each "Generate Layout" run.
2. User submits the SLURM sweep exported by optimizer.generate_slurm_sweep().
3. New results land in the results directory.
4. GUI calls optimizer.ingest_results(path) → CSV grows, surrogate refits,
   solver is re-run, new iteration is recorded.
5. Repeat until convergence (layout_score → 0, uncertainty shrinks).
"""
from __future__ import annotations

import datetime
import json
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class IterationRecord:
    iteration_id: int
    timestamp: str                        # ISO-8601
    ap_rad_m: float
    gap_m: float
    proposed_offsets_m: List[float]       # unique |ACCEL_OFF_Y_M| values used
    slurm_sweep_path: Optional[str]       # path of last exported .slurm file
    ingested_run_dirs: List[str]          # RUN_DIRs added during this iteration
    surrogate_r2: Optional[float]         # surrogate R² after last ingest
    uncertainty_mean_deg: Optional[float] # mean bootstrap std (degrees)
    layout_score: Optional[float]         # clipped_count + shortfall_m * 1000


@dataclass
class IterationLog:
    records: List[IterationRecord] = field(default_factory=list)
    csv_path: str = ""                    # master runlog CSV being enriched

    # ------------------------------------------------------------------
    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "csv_path": self.csv_path,
            "records": [asdict(r) for r in self.records],
        }
        path.write_text(json.dumps(data, indent=2))

    @classmethod
    def load(cls, path: Path) -> "IterationLog":
        if not path.exists():
            return cls()
        raw = json.loads(path.read_text())
        records = [IterationRecord(**r) for r in raw.get("records", [])]
        return cls(records=records, csv_path=raw.get("csv_path", ""))

    def append(self, record: IterationRecord) -> None:
        self.records.append(record)

    def latest(self) -> Optional[IterationRecord]:
        return self.records[-1] if self.records else None

    def next_id(self) -> int:
        return (self.records[-1].iteration_id + 1) if self.records else 1


# ---------------------------------------------------------------------------
# Surrogate model
# ---------------------------------------------------------------------------

class SurrogateModel:
    """
    Lightweight 1-D surrogate: STEER_ANGLE_DEG vs |ACCEL_OFF_Y_M| for a
    fixed (AP_RAD_M, GAP_M) slice of the CSV.

    Uncertainty is estimated via bootstrap resampling of the data points,
    each bootstrap sample fitted with a PchipInterpolator.
    """

    N_BOOTSTRAP = 200
    _RNG_SEED = 42

    def __init__(self) -> None:
        self._interps: List[PchipInterpolator] = []
        self._x_data: Optional[np.ndarray] = None
        self._y_data: Optional[np.ndarray] = None
        self._r2: Optional[float] = None
        self.fitted = False

    # ------------------------------------------------------------------
    def fit(
        self,
        df: pd.DataFrame,
        ap_rad_m: float,
        gap_m: float,
        max_loss_frac: float = 0.10,
        radius_tol: float = 1e-6,
        gap_tol: float = 1e-6,
    ) -> None:
        """Filter CSV to the requested geometry slice and fit bootstrap ensemble."""
        self.fitted = False
        self._interps = []
        self._x_data = None
        self._y_data = None
        self._r2 = None

        needed = {"AP_RAD_M", "GAP_M", "ACCEL_OFF_Y_M", "STEER_ANGLE_DEG", "GRID_LOSS_FRAC"}
        missing = needed - set(df.columns)
        if missing:
            print(f"[SurrogateModel] FAIL — missing columns: {missing}")
            print(f"[SurrogateModel] CSV columns present: {list(df.columns)}")
            return

        n_total = len(df)
        mask_rad = np.abs(df["AP_RAD_M"] - ap_rad_m) <= radius_tol
        mask_gap = np.abs(df["GAP_M"] - gap_m) <= gap_tol
        mask_loss = df["GRID_LOSS_FRAC"] <= max_loss_frac

        print(f"[SurrogateModel] CSV rows total       : {n_total}")
        print(f"[SurrogateModel] Rows pass radius tol : {mask_rad.sum()}  "
              f"(ap_rad_m={ap_rad_m:.6g}, tol=±{radius_tol:.3g})")
        print(f"[SurrogateModel] Rows pass gap tol    : {mask_gap.sum()}  "
              f"(gap_m={gap_m:.6g}, tol=±{gap_tol:.3g})")
        print(f"[SurrogateModel] Rows pass loss frac  : {mask_loss.sum()}  "
              f"(max_loss_frac={max_loss_frac})")

        sub = df[mask_rad & mask_gap & mask_loss].dropna(
            subset=["ACCEL_OFF_Y_M", "STEER_ANGLE_DEG"]
        ).copy()
        print(f"[SurrogateModel] Rows after combined filter + dropna: {len(sub)}")

        if "LOST_TO_SIDEWALLS" in sub.columns:
            before = len(sub)
            sub = sub[sub["LOST_TO_SIDEWALLS"].isin([False, 0, "False", "0", "false"])]
            print(f"[SurrogateModel] Rows after sidewall filter: {len(sub)}  "
                  f"(dropped {before - len(sub)} sidewall hits)")

        if len(sub) < 2:
            print(f"[SurrogateModel] FAIL — only {len(sub)} row(s) survived filters, need ≥2")
            return

        # Group by |offset|, average steering angle
        sub["_off_abs"] = sub["ACCEL_OFF_Y_M"].abs()
        grouped = (
            sub.groupby("_off_abs")["STEER_ANGLE_DEG"]
            .mean()
            .reset_index()
            .sort_values("_off_abs")
        )
        self._x_data = grouped["_off_abs"].values.astype(float)
        self._y_data = grouped["STEER_ANGLE_DEG"].values.astype(float)
        print(f"[SurrogateModel] Unique offset groups  : {len(self._x_data)}")

        if len(self._x_data) < 2:
            print("[SurrogateModel] FAIL — need ≥2 unique offset groups after groupby")
            return

        rng = np.random.default_rng(self._RNG_SEED)
        n = len(self._x_data)

        for _ in range(self.N_BOOTSTRAP):
            idx = rng.integers(0, n, size=n)
            xs = self._x_data[idx]
            ys = self._y_data[idx]
            order = np.argsort(xs)
            xs, ys = xs[order], ys[order]
            _, uniq = np.unique(xs, return_index=True)
            xs, ys = xs[uniq], ys[uniq]
            if len(xs) < 2:
                continue
            try:
                self._interps.append(PchipInterpolator(xs, ys, extrapolate=False))
            except Exception:
                pass

        if not self._interps:
            return

        # R² on the grouped data points
        pred = np.nanmean(
            [[float(f(x)) for x in self._x_data] for f in self._interps], axis=0
        )
        ss_res = float(np.nansum((self._y_data - pred) ** 2))
        ss_tot = float(np.nansum((self._y_data - np.nanmean(self._y_data)) ** 2))
        self._r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else float("nan")
        self.fitted = True

    # ------------------------------------------------------------------
    def predict(self, offsets_m) -> Tuple[np.ndarray, np.ndarray]:
        """Return (mean_deg, std_deg) arrays. NaN where extrapolation would occur."""
        offsets_m = np.asarray(offsets_m, dtype=float)
        nan_arr = np.full(offsets_m.shape, float("nan"))
        if not self.fitted or not self._interps:
            return nan_arr.copy(), nan_arr.copy()
        preds = np.array([[float(f(x)) for x in offsets_m] for f in self._interps])
        mean = np.nanmean(preds, axis=0)
        std = np.nanstd(preds, axis=0)
        return mean, std

    # ------------------------------------------------------------------
    def suggest_next_offsets(
        self,
        n: int,
        existing_offsets_m: Optional[List[float]] = None,
        offset_range_m: Optional[Tuple[float, float]] = None,
    ) -> List[float]:
        """
        Return up to n offsets with the highest prediction uncertainty,
        avoiding already-tested values (active learning / exploration).
        """
        if not self.fitted or not self._interps or n < 1:
            return []

        if offset_range_m is None and self._x_data is not None and len(self._x_data) >= 2:
            offset_range_m = (float(self._x_data[0]), float(self._x_data[-1]))
        if offset_range_m is None:
            return []

        grid = np.linspace(offset_range_m[0], offset_range_m[1], 500)
        _, std = self.predict(grid)
        std = np.nan_to_num(std, nan=0.0)

        if existing_offsets_m:
            spacing = (offset_range_m[1] - offset_range_m[0]) / 500.0
            for ex in existing_offsets_m:
                mask = np.abs(grid - abs(ex)) < 3.0 * spacing
                std[mask] = 0.0

        top_idx = np.argsort(std)[::-1][:n]
        return sorted(float(grid[i]) for i in top_idx)

    # ------------------------------------------------------------------
    def coverage_report(self) -> Dict:
        if not self.fitted or self._x_data is None or len(self._x_data) == 0:
            return {
                "n_data_points": 0,
                "offset_range_m": None,
                "mean_uncertainty_deg": None,
                "r2_score": None,
            }
        _, std = self.predict(self._x_data)
        return {
            "n_data_points": int(len(self._x_data)),
            "offset_range_m": [float(self._x_data.min()), float(self._x_data.max())],
            "mean_uncertainty_deg": float(np.nanmean(std)),
            "r2_score": float(self._r2) if self._r2 is not None else None,
        }


# ---------------------------------------------------------------------------
# Main optimizer class
# ---------------------------------------------------------------------------

class LayoutOptimizer:
    """
    Wraps beamlet_layout_solver with an iterative feedback loop.

    Parameters
    ----------
    log_path : Path
        Path to the iteration log JSON file (created if absent).
    csv_path : Path
        Path to the master runlog_compact CSV.
    """

    def __init__(self, log_path: Path, csv_path: Path) -> None:
        self.log_path = Path(log_path)
        self.csv_path = Path(csv_path)
        self.log = IterationLog.load(self.log_path)
        self.log.csv_path = str(self.csv_path)
        self.surrogate = SurrogateModel()
        self._df: Optional[pd.DataFrame] = None
        self._reload_csv()

    # ------------------------------------------------------------------
    # CSV helpers
    # ------------------------------------------------------------------

    def _reload_csv(self) -> None:
        if self.csv_path.exists():
            try:
                self._df = pd.read_csv(self.csv_path)
                for col in self._df.columns:
                    coerced = pd.to_numeric(self._df[col], errors="coerce")
                    if coerced.notna().any():
                        self._df[col] = coerced
            except Exception:
                self._df = None
        else:
            self._df = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def refit_surrogate(
        self,
        ap_rad_m: float,
        gap_m: float,
        max_loss_frac: float = 0.10,
        radius_tol: float = 1e-6,
        gap_tol: float = 1e-6,
    ) -> Dict:
        """
        Refit the surrogate on the current CSV for (ap_rad_m, gap_m).
        Returns coverage_report dict.
        """
        if self._df is not None:
            self.surrogate.fit(
                self._df, ap_rad_m, gap_m, max_loss_frac, radius_tol, gap_tol
            )
        return self.surrogate.coverage_report()

    # ------------------------------------------------------------------
    def record_solution(self, sol, ap_rad_m: float, gap_m: float) -> IterationRecord:
        """
        Record a new iteration from a LayoutSolution.
        Called automatically by the GUI after each successful solve.

        Returns the created IterationRecord.
        """
        proposed = sorted(
            {round(abs(b.accel_offset_mag_m), 12) for b in sol.beamlets}
        )
        layout_score = (
            float(sol.summary.clipped_beamlet_count)
            + float(sol.summary.max_target_shortfall_m) * 1000.0
        )
        cov = self.surrogate.coverage_report()
        rec = IterationRecord(
            iteration_id=self.log.next_id(),
            timestamp=datetime.datetime.now().isoformat(timespec="seconds"),
            ap_rad_m=ap_rad_m,
            gap_m=gap_m,
            proposed_offsets_m=proposed,
            slurm_sweep_path=None,
            ingested_run_dirs=[],
            surrogate_r2=cov["r2_score"],
            uncertainty_mean_deg=cov["mean_uncertainty_deg"],
            layout_score=layout_score,
        )
        self.log.append(rec)
        self.log.save(self.log_path)
        return rec

    # ------------------------------------------------------------------
    def ingest_results(
        self,
        source: Path,
        ap_rad_m: Optional[float] = None,
        gap_m: Optional[float] = None,
        max_loss_frac: float = 0.10,
        radius_tol: float = 1e-6,
        gap_tol: float = 1e-6,
    ) -> Dict:
        """
        Ingest new simulation results into the master CSV and refit surrogate.

        Parameters
        ----------
        source : Path
            Either a CSV file (merged directly) or a results directory
            (update_runlog_compact.py is invoked as a subprocess).
        ap_rad_m, gap_m : floats or None
            If provided, the surrogate is refit for this geometry slice.
            If None, uses the latest iteration's values (if any).

        Returns
        -------
        dict with keys: n_rows_before, n_rows_after, n_new_rows,
                        surrogate_r2, uncertainty_mean_deg, new_run_dirs
        """
        source = Path(source)
        n_before = len(self._df) if self._df is not None else 0
        new_run_dirs: List[str] = []

        if source.is_dir():
            # Delegate to update_runlog_compact.py
            cmd = [
                sys.executable,
                str(Path(__file__).parent / "update_runlog_compact.py"),
                "--results-dir", str(source),
                "--csv", str(self.csv_path),
                "--scan",
                "--refresh",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(
                    f"update_runlog_compact.py failed:\n{result.stderr}"
                )
            self._reload_csv()
        elif source.suffix.lower() == ".csv":
            # Merge CSV rows not already in master
            incoming = pd.read_csv(source)
            if self._df is not None and "RUN_DIR" in self._df.columns and "RUN_DIR" in incoming.columns:
                # Primary dedup: by RUN_DIR
                existing_dirs = set(self._df["RUN_DIR"].dropna().astype(str))
                new_rows = incoming[~incoming["RUN_DIR"].astype(str).isin(existing_dirs)]
                new_run_dirs = list(new_rows["RUN_DIR"].dropna().astype(str))
                if len(new_rows) > 0:
                    merged = pd.concat([self._df, new_rows], ignore_index=True)
                    merged.to_csv(self.csv_path, index=False)
            elif self._df is not None:
                # Fallback dedup: composite key (AP_RAD_M, GAP_M, ACCEL_OFF_Y_M)
                _COMPOSITE_KEY = ["AP_RAD_M", "GAP_M", "ACCEL_OFF_Y_M"]
                key_cols = [c for c in _COMPOSITE_KEY if c in self._df.columns and c in incoming.columns]
                if key_cols:
                    existing_keys = set(
                        self._df[key_cols].dropna().round(12).apply(tuple, axis=1)
                    )
                    mask = ~incoming[key_cols].round(12).apply(tuple, axis=1).isin(existing_keys)
                    new_rows = incoming[mask]
                else:
                    new_rows = incoming
                new_run_dirs = list(new_rows.get("RUN_DIR", pd.Series(dtype=str)).dropna().astype(str))
                if len(new_rows) > 0:
                    merged = pd.concat([self._df, new_rows], ignore_index=True)
                    merged.to_csv(self.csv_path, index=False)
            else:
                # No existing CSV — just write incoming
                incoming.to_csv(self.csv_path, index=False)
                new_run_dirs = list(incoming.get("RUN_DIR", pd.Series(dtype=str)).dropna().astype(str))
            self._reload_csv()
        else:
            raise ValueError(f"source must be a directory or .csv file, got: {source}")

        n_after = len(self._df) if self._df is not None else 0
        n_new = n_after - n_before

        # Determine geometry for surrogate refit
        if ap_rad_m is None or gap_m is None:
            latest = self.log.latest()
            if latest is not None:
                ap_rad_m = ap_rad_m or latest.ap_rad_m
                gap_m = gap_m or latest.gap_m

        cov: Dict = {"r2_score": None, "mean_uncertainty_deg": None}
        if ap_rad_m is not None and gap_m is not None:
            cov = self.refit_surrogate(ap_rad_m, gap_m, max_loss_frac, radius_tol, gap_tol)

        # Update the latest iteration record with ingest info
        latest = self.log.latest()
        if latest is not None:
            latest.ingested_run_dirs = new_run_dirs
            latest.surrogate_r2 = cov.get("r2_score")
            latest.uncertainty_mean_deg = cov.get("mean_uncertainty_deg")
            self.log.save(self.log_path)

        return {
            "n_rows_before": n_before,
            "n_rows_after": n_after,
            "n_new_rows": n_new,
            "surrogate_r2": cov.get("r2_score"),
            "uncertainty_mean_deg": cov.get("mean_uncertainty_deg"),
            "new_run_dirs": new_run_dirs,
        }

    # ------------------------------------------------------------------
    def generate_slurm_sweep(
        self,
        out_path: Path,
        orig_slurm_path: Path,
        proposed_offsets_m: List[float],
        ap_rad_m: float,
        gap_m: float,
        aperture_pairs_m: str,
        delta_m: float = 5e-5,
        n_steps: int = 3,
        run_prefix: str = "optimizer_sweep",
        iteration_id: Optional[int] = None,
        fast_mode: bool = False,
    ) -> List[float]:
        """
        Generate a bash wrapper that patches only the sweep-list lines in
        ``orig_slurm_path`` via sed and submits the result to sbatch.

        The original orchestrate_pipeline.slurm is never modified.  All
        resource settings, physics parameters, chamfer values, etc. stay
        exactly as they are in that file.

        Lines patched by sed (regex anchored to start of line):
          AP_RAD_LIST, SCREEN_AP_RAD_LIST, ACCEL_AP_RAD_LIST
          GAP_LIST
          ACCEL_OFF_LIST  (the make_seq line is replaced with a literal list)
          RUN_PREFIX

        APERTURE_PAIRS_M uses ``${VAR:-default}`` in the original script, so
        it is passed as an exported env var rather than patched.

        Returns the sorted list of new ACCEL_OFF_Y_M values to be swept.
        """
        out_path = Path(out_path)
        orig_slurm_path = Path(orig_slurm_path)

        # ------------------------------------------------------------------
        # Build neighborhood around each proposed offset
        # ------------------------------------------------------------------
        sweep_set: set = set()
        for off in proposed_offsets_m:
            for k in range(-n_steps, n_steps + 1):
                val = round(abs(off) + k * delta_m, 12)
                if val >= 0.0:
                    sweep_set.add(val)
                if abs(off) > 1e-12:
                    val_neg = round(-abs(off) + k * delta_m, 12)
                    sweep_set.add(val_neg)

        # Remove offsets already tested for this geometry slice
        if self._df is not None:
            needed = {"AP_RAD_M", "GAP_M", "ACCEL_OFF_Y_M"}
            if needed.issubset(self._df.columns):
                slice_df = self._df[
                    (np.abs(self._df["AP_RAD_M"] - ap_rad_m) <= 1e-9)
                    & (np.abs(self._df["GAP_M"] - gap_m) <= 1e-9)
                ]
                existing_off = set(slice_df["ACCEL_OFF_Y_M"].dropna().round(12))
                sweep_set -= existing_off

        # Physical constraint: offset cannot exceed aperture radius
        sweep_set = {v for v in sweep_set if abs(v) <= ap_rad_m}

        new_offsets = sorted(sweep_set)
        if not new_offsets:
            # Fallback: full neighborhood ignoring dedup (prevents empty sweep)
            new_offsets = sorted({
                round(abs(off) + k * delta_m, 12)
                for off in proposed_offsets_m
                for k in range(-n_steps, n_steps + 1)
                if 0.0 <= round(abs(off) + k * delta_m, 12) <= ap_rad_m
            })

        iter_tag = f"iter{iteration_id}" if iteration_id is not None else "iterN"
        off_list_str = " ".join(f"{v:.12g}" for v in new_offsets)
        n_total = len(new_offsets)

        # ------------------------------------------------------------------
        # Build sed expressions — use | as delimiter to avoid / conflicts
        # ------------------------------------------------------------------
        ap_patch  = f"( {ap_rad_m:.12g} )"
        gap_patch = f"( {gap_m:.12g} )"
        off_patch = f"( {off_list_str} )"

        # Escape | and & for sed replacement strings
        def _sed_esc(s: str) -> str:
            return s.replace("\\", "\\\\").replace("|", "\\|").replace("&", "\\&")

        sed_lines = [
            f"s|^AP_RAD_LIST=.*|AP_RAD_LIST={_sed_esc(ap_patch)}|",
            f"s|^SCREEN_AP_RAD_LIST=.*|SCREEN_AP_RAD_LIST={_sed_esc(ap_patch)}|",
            f"s|^ACCEL_AP_RAD_LIST=.*|ACCEL_AP_RAD_LIST={_sed_esc(ap_patch)}|",
            f"s|^GAP_LIST=.*|GAP_LIST={_sed_esc(gap_patch)}|",
            # Replace the make_seq call line with a literal list
            f"s|^ACCEL_OFF_LIST=.*|ACCEL_OFF_LIST={_sed_esc(off_patch)}|",
            f's|^RUN_PREFIX=.*|RUN_PREFIX="{_sed_esc(run_prefix)}"|',
        ]

        if fast_mode:
            sed_lines += [
                "s|^WRITE_PNG=.*|WRITE_PNG=0|",
                "s|^PLOT_PROFILE=.*|PLOT_PROFILE=0|",
                "s|^WRITE_DENSITY_PNG=.*|WRITE_DENSITY_PNG=0|",
                "s|^WRITE_DENSITY_GRID=.*|WRITE_DENSITY_GRID=0|",
                "s|^PLOT_DENSITY=.*|PLOT_DENSITY=0|",
                "s|^X_RIGHT_M=.*|X_RIGHT_M=0.02|",
            ]
        sed_expr = " \\\n  ".join(f"-e '{e}'" for e in sed_lines)

        # Escape APERTURE_PAIRS_M for shell single-quote (replace ' with '"'"')
        pairs_escaped = aperture_pairs_m.replace("'", "'\"'\"'")

        fast_note = "FAST MODE: vis outputs disabled, X_RIGHT_M pinned to 0.02 m" if fast_mode else "full mode"
        script = f"""\
#!/bin/bash
# ============================================================
# layout_optimizer — targeted sweep wrapper — {iter_tag}
# Generated : {datetime.datetime.now().isoformat(timespec="seconds")}
# AP_RAD_M  : {ap_rad_m:.6g} m
# GAP_M     : {gap_m:.6g} m
# delta_m   : {delta_m:.3g} m   n_steps : {n_steps}
# New offsets: {n_total}
# Mode      : {fast_note}
#
# Usage: bash {out_path.name}
#   (run from the same directory as orchestrate_pipeline.slurm)
# ============================================================
set -euo pipefail

ORIG_SLURM="{orig_slurm_path}"

if [[ ! -f "$ORIG_SLURM" ]]; then
  echo "ERROR: orchestrate_pipeline.slurm not found at: $ORIG_SLURM" >&2
  exit 1
fi

# APERTURE_PAIRS_M uses ${{VAR:-default}} in the original script — export it
export APERTURE_PAIRS_M='{pairs_escaped}'

# Patch sweep lists and pipe to sbatch
TMP=$(mktemp /tmp/optimizer_{iter_tag}_XXXX.slurm)
trap 'rm -f "$TMP"' EXIT

sed \\
  {sed_expr} \\
  "$ORIG_SLURM" > "$TMP"

echo "================================================"
echo " layout_optimizer sweep — {iter_tag}"
echo " Mode             : {fast_note}"
echo " AP_RAD_LIST      : {ap_patch}"
echo " GAP_LIST         : {gap_patch}"
echo " ACCEL_OFF_LIST   : {off_list_str}"
echo " APERTURE_PAIRS_M : $APERTURE_PAIRS_M"
echo " Submitting       : $TMP"
echo "================================================"
sbatch "$TMP"
"""

        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(script)

        # Record the sweep path in the latest iteration
        latest = self.log.latest()
        if latest is not None:
            latest.slurm_sweep_path = str(out_path)
            self.log.save(self.log_path)

        return new_offsets

    # ------------------------------------------------------------------
    def generate_grid_sweep(
        self,
        out_path: Path,
        orig_slurm_path: Path,
        gap_values_m: List[float],
        ap_rad_m: float,
        offset_step_m: float = 2e-4,
        run_prefix: str = "grid_sweep",
        fast_mode: bool = False,
    ) -> List[Tuple[float, float]]:
        """
        Generate a bash wrapper for a 2D (gap × offset) characterisation sweep
        in single-aperture mode.

        Unlike generate_slurm_sweep(), this method:
          - Sweeps multiple GAP values simultaneously
          - Caps offsets at ap_rad_m (physical constraint)
          - Exports APERTURE_PAIRS_M="" to force two_grid_2d.cpp into
            single-aperture mode, so ACCEL_OFF_Y_M env var is used directly
            rather than acc_off_list.front() derived from APERTURE_PAIRS_M.

        Returns sorted list of (gap_m, offset_m) tuples that will be run.
        """
        out_path = Path(out_path)
        orig_slurm_path = Path(orig_slurm_path)

        # Build offset grid capped at aperture radius
        offset_arr = np.arange(0.0, ap_rad_m + 1e-12, offset_step_m)
        offset_list = sorted({round(float(v), 12) for v in offset_arr})

        # Deduplicate against existing runs on (AP_RAD_M, GAP_M, ACCEL_OFF_Y_M)
        planned: List[Tuple[float, float]] = []
        if self._df is not None and {"AP_RAD_M", "GAP_M", "ACCEL_OFF_Y_M"}.issubset(self._df.columns):
            for g in gap_values_m:
                slice_df = self._df[
                    (np.abs(self._df["AP_RAD_M"] - ap_rad_m) <= 1e-9)
                    & (np.abs(self._df["GAP_M"] - g) <= 1e-9)
                ]
                existing_off = set(slice_df["ACCEL_OFF_Y_M"].dropna().round(12))
                for o in offset_list:
                    if round(o, 12) not in existing_off:
                        planned.append((g, o))
        else:
            planned = [(g, o) for g in gap_values_m for o in offset_list]

        if not planned:
            planned = [(g, o) for g in gap_values_m for o in offset_list]

        # Build list strings for patching
        gap_list_str = " ".join(f"{g:.12g}" for g in sorted({g for g, _ in planned}))
        off_list_str  = " ".join(f"{o:.12g}" for o in sorted({o for _, o in planned}))
        ap_patch  = f"( {ap_rad_m:.12g} )"
        gap_patch = f"( {gap_list_str} )"
        off_patch = f"( {off_list_str} )"
        n_gaps    = len({g for g, _ in planned})
        n_offsets = len({o for _, o in planned})

        def _sed_esc(s: str) -> str:
            return s.replace("\\", "\\\\").replace("|", "\\|").replace("&", "\\&")

        sed_lines = [
            f"s|^AP_RAD_LIST=.*|AP_RAD_LIST={_sed_esc(ap_patch)}|",
            f"s|^SCREEN_AP_RAD_LIST=.*|SCREEN_AP_RAD_LIST={_sed_esc(ap_patch)}|",
            f"s|^ACCEL_AP_RAD_LIST=.*|ACCEL_AP_RAD_LIST={_sed_esc(ap_patch)}|",
            f"s|^GAP_LIST=.*|GAP_LIST={_sed_esc(gap_patch)}|",
            f"s|^ACCEL_OFF_LIST=.*|ACCEL_OFF_LIST={_sed_esc(off_patch)}|",
            f's|^RUN_PREFIX=.*|RUN_PREFIX="{_sed_esc(run_prefix)}"|',
        ]
        if fast_mode:
            sed_lines += [
                "s|^WRITE_PNG=.*|WRITE_PNG=0|",
                "s|^PLOT_PROFILE=.*|PLOT_PROFILE=0|",
                "s|^WRITE_DENSITY_PNG=.*|WRITE_DENSITY_PNG=0|",
                "s|^WRITE_DENSITY_GRID=.*|WRITE_DENSITY_GRID=0|",
                "s|^PLOT_DENSITY=.*|PLOT_DENSITY=0|",
                "s|^X_RIGHT_M=.*|X_RIGHT_M=0.02|",
            ]
        sed_expr = " \\\n  ".join(f"-e '{e}'" for e in sed_lines)
        fast_note = "FAST MODE" if fast_mode else "full mode"

        gap_labels = [f"{g*1e3:.4g}" for g in sorted({g for g, _ in planned})]

        script = f"""\
#!/bin/bash
# ============================================================
# layout_optimizer — 2D grid characterisation sweep
# Generated : {datetime.datetime.now().isoformat(timespec="seconds")}
# AP_RAD_M  : {ap_rad_m:.6g} m  (max offset = {ap_rad_m*1e3:.4g} mm)
# GAP values : {gap_labels} mm
# Offsets    : 0 to {ap_rad_m*1e3:.4g} mm in {offset_step_m*1e3:.4g} mm steps
# Total runs : {n_gaps} gaps x {n_offsets} offsets = {n_gaps * n_offsets} runs
# Mode       : {fast_note}
#
# APERTURE_PAIRS_M is forced to "" so two_grid_2d.cpp uses
# ACCEL_OFF_Y_M directly (single-aperture characterisation mode).
#
# Usage: bash {out_path.name}
#   (run from the same directory as orchestrate_pipeline.slurm)
# ============================================================
set -euo pipefail

ORIG_SLURM="{orig_slurm_path}"

if [[ ! -f "$ORIG_SLURM" ]]; then
  echo "ERROR: orchestrate_pipeline.slurm not found at: $ORIG_SLURM" >&2
  exit 1
fi

# Force single-aperture mode — must be empty so two_grid_2d.cpp falls back
# to the scalar ACCEL_OFF_Y_M env var set from ACCEL_OFF_LIST each iteration.
export APERTURE_PAIRS_M=""

# Patch sweep lists and submit
TMP=$(mktemp /tmp/grid_sweep_XXXX.slurm)
trap 'rm -f "$TMP"' EXIT

sed \\
  {sed_expr} \\
  "$ORIG_SLURM" > "$TMP"

echo "================================================"
echo " layout_optimizer 2D grid sweep"
echo " Mode             : {fast_note}"
echo " AP_RAD_LIST      : {ap_patch}"
echo " GAP_LIST         : {gap_patch}"
echo " ACCEL_OFF_LIST   : {off_patch}"
echo " APERTURE_PAIRS_M : (empty — single-aperture mode)"
echo " Total runs       : {n_gaps} x {n_offsets} = {n_gaps * n_offsets}"
echo " Submitting       : $TMP"
echo "================================================"
sbatch "$TMP"
"""
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(script)
        return sorted(planned)

    # ------------------------------------------------------------------
    def get_convergence_data(self) -> Dict:
        """
        Return lists suitable for plotting convergence over iterations.

        Returns
        -------
        dict with keys: iteration_ids, layout_scores, surrogate_r2s,
                        uncertainty_mean_degs, timestamps
        """
        ids, scores, r2s, uncerts, stamps = [], [], [], [], []
        for rec in self.log.records:
            ids.append(rec.iteration_id)
            scores.append(rec.layout_score)
            r2s.append(rec.surrogate_r2)
            uncerts.append(rec.uncertainty_mean_deg)
            stamps.append(rec.timestamp)
        return {
            "iteration_ids": ids,
            "layout_scores": scores,
            "surrogate_r2s": r2s,
            "uncertainty_mean_degs": uncerts,
            "timestamps": stamps,
        }

    # ------------------------------------------------------------------
    def suggest_next_offsets(
        self,
        n: int = 3,
        ap_rad_m: Optional[float] = None,
        gap_m: Optional[float] = None,
        max_loss_frac: float = 0.10,
        radius_tol: float = 1e-6,
        gap_tol: float = 1e-6,
    ) -> List[float]:
        """
        Convenience wrapper: refit (if needed) and return active-learning
        offset suggestions for the given geometry.
        """
        if ap_rad_m is None or gap_m is None:
            latest = self.log.latest()
            if latest is None:
                return []
            ap_rad_m = ap_rad_m or latest.ap_rad_m
            gap_m = gap_m or latest.gap_m

        if not self.surrogate.fitted:
            self.refit_surrogate(ap_rad_m, gap_m, max_loss_frac, radius_tol, gap_tol)

        existing = self.log.latest().proposed_offsets_m if self.log.latest() else []
        return self.surrogate.suggest_next_offsets(n, existing_offsets_m=existing)
