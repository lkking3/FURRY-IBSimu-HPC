#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import traceback
import webbrowser
from dataclasses import asdict
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import beamlet_layout_solver as solver

try:
    import layout_optimizer as opt_mod
    _HAS_OPTIMIZER = True
except ImportError:
    _HAS_OPTIMIZER = False

def _find_curved_pipeline() -> bool:
    """Search well-known sibling directories for curved_grid_pipeline.py and
    add the first match to sys.path so it can be imported."""
    import sys
    _THIS_DIR = Path(__file__).resolve().parent
    _CANDIDATES = [
        _THIS_DIR,
        _THIS_DIR.parent / "curved_grid",
        _THIS_DIR / "curved_grid",
        _THIS_DIR / "two_grid_curved",
        _THIS_DIR / "Multi_Grid",
        _THIS_DIR.parent / "two_grid_curved",
        _THIS_DIR.parent / "Multi_Grid",
    ]
    for d in _CANDIDATES:
        if (d / "curved_grid_pipeline.py").exists():
            p = str(d)
            if p not in sys.path:
                sys.path.insert(0, p)
            return True
    return False

_HAS_CURVED = False
if _find_curved_pipeline():
    try:
        from curved_grid_pipeline import (
            load_curved_case,
            bore_clearance_report,
            pairs_from_case,
            write_case_file,
            plot_curved_grid_profile,
            uniform_screen_offsets,
            accel_offset_from_steering,
            concentric_accel_radius,
        )
        _HAS_CURVED = True
    except Exception:
        _HAS_CURVED = False


def _find_sim_dir() -> bool:
    """Add the sim/ directory alongside the repo root to sys.path so
    sweep_types (and grid_types) can be imported by the GUI."""
    import sys
    _THIS_DIR = Path(__file__).resolve().parent
    _CANDIDATES = [
        _THIS_DIR.parent / "sim",
        _THIS_DIR / "sim",
        _THIS_DIR,
    ]
    for d in _CANDIDATES:
        if (d / "sweep_types.py").exists():
            p = str(d)
            if p not in sys.path:
                sys.path.insert(0, p)
            return True
    return False

_HAS_SWEEP = False
if _find_sim_dir():
    try:
        from sweep_types import SweepAxis, SweepSpec, linspace, arange
        _HAS_SWEEP = True
    except Exception:
        _HAS_SWEEP = False


class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Beamlet Layout Solver")
        self.geometry("1500x900")
        self.minsize(1200, 760)

        self.last_plot_paths = {}
        self.last_solution = None
        self._build_vars()
        self._build_ui()
        self._set_defaults()

    def _build_vars(self) -> None:
        self.csv_var = tk.StringVar()
        self.radius_var = tk.StringVar()
        self.gap_var = tk.StringVar()
        self.max_loss_var = tk.StringVar()
        self.packing_var = tk.StringVar()
        self.screen_pitch_var = tk.StringVar()
        self.target_mode_var = tk.StringVar()
        self.plate_radius_var = tk.StringVar()
        self.target_radius_var = tk.StringVar()
        self.target_plane_var = tk.StringVar()
        self.steer_origin_var = tk.StringVar()
        self.screen_center_x_var = tk.StringVar()
        self.screen_center_y_var = tk.StringVar()
        self.target_center_x_var = tk.StringVar()
        self.target_center_y_var = tk.StringVar()
        self.count_var = tk.StringVar()
        self.radius_tol_var = tk.StringVar()
        self.gap_tol_var = tk.StringVar()
        self.extra_filter_var = tk.StringVar()
        self.out_json_var = tk.StringVar()
        self.plot_prefix_var = tk.StringVar()
        self.clip_var = tk.BooleanVar(value=True)
        self.auto_open_var = tk.BooleanVar(value=False)
        self.capacity_var = tk.StringVar()
        self.preview_status_var = tk.StringVar(value="Hover over a beamlet in the preview to inspect coordinates.")
        self.pairs_status_var = tk.StringVar(value="APERTURE_PAIRS_M not generated yet.")
        # Optimizer state
        self._optimizer: "opt_mod.LayoutOptimizer | None" = None
        self.opt_log_var = tk.StringVar()
        self.opt_status_var = tk.StringVar(value="No optimizer log loaded.")
        self.opt_surrogate_var = tk.StringVar(value="")
        self.opt_suggestions_var = tk.StringVar(value="")
        # Curved grid import/export state
        self.curved_case_var = tk.StringVar()
        self._loaded_curved_case = None
        self.curved_status_var = tk.StringVar(value="No case loaded.")
        self.curved_action_status_var = tk.StringVar(value="")
        # Aperture layout generation (curved grid tab)
        self.ap_n_var = tk.StringVar(value="3")
        self.ap_pitch_var = tk.StringVar(value="10.0")
        self.ap_offset_var = tk.StringVar(value="0.0")
        self.ap_steering_var = tk.StringVar(value="0.0")
        self.ap_axial_var = tk.BooleanVar(value=True)
        self.ap_offsets_label_var = tk.StringVar(value="Not applied — using loaded case.")
        self._aperture_offsets_override = None  # List[float] or None
        self._steering_arcs_override = None     # List[float] or None
        # Parameter sweep state
        self._sweep_axes: list = []             # List[SweepAxis]
        self._sweep_mode_var  = tk.StringVar(value="product")
        self._sweep_reduce_var = tk.BooleanVar(value=False)
        self._sweep_tag_var   = tk.StringVar(value="sweep")

    def _set_defaults(self) -> None:
        self.csv_var.set(str(Path("results") / "runlog_compact_Aprad_offset.csv"))
        self.radius_var.set("1.5")
        self.gap_var.set("2")
        self.max_loss_var.set("0.10")
        self.packing_var.set("radial")
        self.screen_pitch_var.set("4")
        self.target_mode_var.set("center")
        self.plate_radius_var.set("86")
        self.target_radius_var.set("33")
        self.target_plane_var.set("550")
        self.steer_origin_var.set("0")
        self.screen_center_x_var.set("0")
        self.screen_center_y_var.set("0")
        self.target_center_x_var.set("0")
        self.target_center_y_var.set("0")
        self.count_var.set("1076")
        self.radius_tol_var.set("0.01")   # 0.01 mm = 10 µm
        self.gap_tol_var.set("0.05")       # 0.05 mm = 50 µm
        self.extra_filter_var.set("")
        self.out_json_var.set(str(Path("results") / "beamlet_layout_gui.json"))
        self.plot_prefix_var.set(str(Path("results") / "beamlet_layout_gui"))
        self._bind_preview_updates()
        self._update_capacity_preview()

    def _build_ui(self) -> None:
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        left = ttk.Frame(self, padding=10)
        left.grid(row=0, column=0, sticky="nsw")

        right = ttk.Frame(self, padding=10)
        right.grid(row=0, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=0)
        right.rowconfigure(2, weight=1)

        frm = ttk.LabelFrame(left, text="Inputs", padding=10)
        frm.grid(row=0, column=0, sticky="new")

        rows = []
        rows.append(("CSV", self.csv_var, self._browse_csv))
        rows.append(("Radius (mm)", self.radius_var, None))
        rows.append(("Gap (mm)", self.gap_var, None))
        rows.append(("Max Loss Frac", self.max_loss_var, None))
        rows.append(("Packing", self.packing_var, "combo"))
        rows.append(("Screen Pitch (mm)", self.screen_pitch_var, None))
        rows.append(("Target Mode", self.target_mode_var, "target_combo"))
        rows.append(("Plate Radius (mm)", self.plate_radius_var, None))
        rows.append(("Target Radius (mm)", self.target_radius_var, None))
        rows.append(("Target Plane (mm)", self.target_plane_var, None))
        rows.append(("Steer Origin (mm)", self.steer_origin_var, None))
        rows.append(("Screen Center X (mm)", self.screen_center_x_var, None))
        rows.append(("Screen Center Y (mm)", self.screen_center_y_var, None))
        rows.append(("Target Center X (mm)", self.target_center_x_var, None))
        rows.append(("Target Center Y (mm)", self.target_center_y_var, None))
        rows.append(("Count", self.count_var, None))
        rows.append(("Radius Tol (mm)", self.radius_tol_var, None))
        rows.append(("Gap Tol (mm)", self.gap_tol_var, None))
        rows.append(("Extra Filter", self.extra_filter_var, None))
        rows.append(("Out JSON", self.out_json_var, self._browse_json))
        rows.append(("Plot Prefix", self.plot_prefix_var, self._browse_prefix))

        for i, (label, var, action) in enumerate(rows):
            ttk.Label(frm, text=label).grid(row=i, column=0, sticky="w", padx=(0, 8), pady=3)
            if action == "combo":
                cb = ttk.Combobox(frm, textvariable=var, values=["linear", "radial", "hcp"], state="readonly", width=24)
                cb.grid(row=i, column=1, sticky="ew", pady=3)
            elif action == "target_combo":
                cb = ttk.Combobox(frm, textvariable=var, values=["center", "packed"], state="readonly", width=24)
                cb.grid(row=i, column=1, sticky="ew", pady=3)
            else:
                ent = ttk.Entry(frm, textvariable=var, width=30)
                ent.grid(row=i, column=1, sticky="ew", pady=3)
                if callable(action):
                    ttk.Button(frm, text="Browse", command=action).grid(row=i, column=2, sticky="ew", padx=(6, 0), pady=3)

        frm.columnconfigure(1, weight=1)

        flags = ttk.Frame(left, padding=(0, 10, 0, 0))
        flags.grid(row=1, column=0, sticky="ew")
        ttk.Checkbutton(flags, text="Clip To Safe Limit", variable=self.clip_var).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(flags, text="Auto Open Plots", variable=self.auto_open_var).grid(row=1, column=0, sticky="w")
        ttk.Label(flags, textvariable=self.capacity_var, wraplength=320, justify="left").grid(row=2, column=0, sticky="w", pady=(8, 0))

        btns = ttk.Frame(left, padding=(0, 12, 0, 0))
        btns.grid(row=2, column=0, sticky="ew")
        ttk.Button(btns, text="Generate Layout", command=self._run_solver).grid(row=0, column=0, sticky="ew")
        ttk.Button(btns, text="Copy APERTURE_PAIRS_M", command=self._copy_pairs_to_clipboard).grid(row=1, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(btns, text="Open Screen Plot", command=lambda: self._open_plot("screen")).grid(row=2, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(btns, text="Open Accel Plot", command=lambda: self._open_plot("accel")).grid(row=3, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(btns, text="Open Overlap Plot", command=lambda: self._open_plot("overlap")).grid(row=4, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(btns, text="Open Side Plot", command=lambda: self._open_plot("side")).grid(row=5, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(btns, text="Open JSON", command=self._open_json).grid(row=6, column=0, sticky="ew", pady=(6, 0))
        ttk.Label(btns, textvariable=self.pairs_status_var, wraplength=320, justify="left").grid(row=7, column=0, sticky="ew", pady=(8, 0))

        ttk.Label(right, textvariable=self.preview_status_var, justify="left").grid(row=0, column=0, sticky="ew")
        self.main_notebook = ttk.Notebook(right)
        self.main_notebook.grid(row=2, column=0, sticky="nsew")

        summary_tab = ttk.Frame(self.main_notebook, padding=8)
        summary_tab.columnconfigure(0, weight=1)
        summary_tab.rowconfigure(0, weight=1)
        self.summary_text = tk.Text(summary_tab, wrap="word")
        self.summary_text.grid(row=0, column=0, sticky="nsew")
        self.main_notebook.add(summary_tab, text="Summary")

        beamlet_tab = ttk.Frame(self.main_notebook, padding=8)
        beamlet_tab.columnconfigure(0, weight=1)
        beamlet_tab.rowconfigure(0, weight=1)
        cols = [
            "index", "screen_x_mm", "screen_y_mm", "accel_x_mm", "accel_y_mm",
            "off_x_mm", "off_y_mm", "off_mag_mm", "target_x_mm", "target_y_mm",
            "theta_req_deg", "theta_used_deg", "clipped", "hits_target", "target_feasible", "shortfall_mm",
        ]
        self.tree = ttk.Treeview(beamlet_tab, columns=cols, show="headings", height=18)
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=95, anchor="center")
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb = ttk.Scrollbar(beamlet_tab, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(beamlet_tab, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        self.main_notebook.add(beamlet_tab, text="Beamlets")

        self.preview_canvases = {}
        for key, title in [("screen", "Screen"), ("accel", "Accelerator"), ("overlap", "Overlap"), ("target", "Target Impact"), ("side", "Side View")]:
            frame = ttk.Frame(self.main_notebook, padding=8)
            frame.columnconfigure(0, weight=1)
            frame.rowconfigure(0, weight=1)
            canvas = tk.Canvas(frame, background="white", highlightthickness=1, highlightbackground="#cccccc")
            canvas.grid(row=0, column=0, sticky="nsew")
            canvas.bind("<Configure>", lambda _e: self._redraw_previews())
            self.main_notebook.add(frame, text=title)
            self.preview_canvases[key] = canvas

        self._build_optimizer_tab()
        self._build_curved_tab()

    # ------------------------------------------------------------------
    # Optimizer tab
    # ------------------------------------------------------------------

    def _build_optimizer_tab(self) -> None:
        opt_tab = ttk.Frame(self.main_notebook, padding=8)
        opt_tab.columnconfigure(0, weight=1)
        self.main_notebook.add(opt_tab, text="Optimizer")

        if not _HAS_OPTIMIZER:
            ttk.Label(opt_tab, text="layout_optimizer.py not found. Place it alongside this script.",
                      foreground="red").grid(row=0, column=0, sticky="w")
            return

        # -- Log file row --
        log_row = ttk.LabelFrame(opt_tab, text="Optimizer Log", padding=6)
        log_row.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        log_row.columnconfigure(1, weight=1)
        ttk.Label(log_row, text="Log file:").grid(row=0, column=0, sticky="w")
        ttk.Entry(log_row, textvariable=self.opt_log_var, width=55).grid(row=0, column=1, sticky="ew", padx=4)
        ttk.Button(log_row, text="Browse", command=self._browse_opt_log).grid(row=0, column=2)
        ttk.Button(log_row, text="Load / Create", command=self._load_optimizer).grid(row=0, column=3, padx=(4, 0))

        # -- Surrogate status --
        surr_frame = ttk.LabelFrame(opt_tab, text="Surrogate Model Status", padding=6)
        surr_frame.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        surr_frame.columnconfigure(0, weight=1)
        ttk.Label(surr_frame, textvariable=self.opt_surrogate_var, justify="left",
                  wraplength=700).grid(row=0, column=0, sticky="w")
        ttk.Label(surr_frame, text="Suggested next offsets (mm):").grid(row=1, column=0, sticky="w", pady=(4, 0))
        ttk.Label(surr_frame, textvariable=self.opt_suggestions_var, foreground="#0055cc",
                  justify="left", wraplength=700).grid(row=2, column=0, sticky="w")

        # -- Iteration log treeview --
        log_frame = ttk.LabelFrame(opt_tab, text="Iteration Log", padding=6)
        log_frame.grid(row=2, column=0, sticky="nsew", pady=(0, 6))
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        opt_tab.rowconfigure(2, weight=1)
        iter_cols = ["#", "Timestamp", "Score", "R²", "Uncertainty (°)", "Clipped", "AP_RAD (mm)", "GAP (mm)"]
        self.iter_tree = ttk.Treeview(log_frame, columns=iter_cols, show="headings", height=8)
        for c in iter_cols:
            self.iter_tree.heading(c, text=c)
            self.iter_tree.column(c, width=110, anchor="center")
        self.iter_tree.column("#", width=40)
        self.iter_tree.column("Timestamp", width=155)
        self.iter_tree.grid(row=0, column=0, sticky="nsew")
        iter_vsb = ttk.Scrollbar(log_frame, orient="vertical", command=self.iter_tree.yview)
        self.iter_tree.configure(yscrollcommand=iter_vsb.set)
        iter_vsb.grid(row=0, column=1, sticky="ns")

        # -- Convergence canvas --
        conv_frame = ttk.LabelFrame(opt_tab, text="Convergence (layout score per iteration)", padding=6)
        conv_frame.grid(row=3, column=0, sticky="ew", pady=(0, 6))
        conv_frame.columnconfigure(0, weight=1)
        self.conv_canvas = tk.Canvas(conv_frame, height=120, background="white",
                                     highlightthickness=1, highlightbackground="#cccccc")
        self.conv_canvas.grid(row=0, column=0, sticky="ew")
        self.conv_canvas.bind("<Configure>", lambda _e: self._draw_convergence())

        # -- Action buttons --
        act_frame = ttk.Frame(opt_tab, padding=(0, 4, 0, 0))
        act_frame.grid(row=4, column=0, sticky="ew")
        ttk.Button(act_frame, text="Ingest Results", command=self._ingest_results).grid(row=0, column=0, padx=(0, 8))
        ttk.Button(act_frame, text="Export SLURM Sweep", command=self._export_slurm_sweep).grid(row=0, column=1)
        ttk.Button(act_frame, text="Export Grid Sweep", command=self._export_grid_sweep).grid(row=0, column=2, padx=(8, 0))
        ttk.Label(act_frame, textvariable=self.opt_status_var, foreground="#555555",
                  wraplength=600, justify="left").grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 0))

    def _browse_opt_log(self) -> None:
        p = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("All Files", "*.*")],
            title="Select or create optimizer log file",
        )
        if p:
            self.opt_log_var.set(p)

    def _load_optimizer(self) -> None:
        if not _HAS_OPTIMIZER:
            return
        log_path_str = self.opt_log_var.get().strip()
        if not log_path_str:
            messagebox.showinfo("Optimizer", "Enter or browse for a log file path first.")
            return
        csv_path_str = self.csv_var.get().strip()
        if not csv_path_str:
            messagebox.showinfo("Optimizer", "Set the CSV path in Inputs first.")
            return
        try:
            self._optimizer = opt_mod.LayoutOptimizer(
                log_path=Path(log_path_str),
                csv_path=Path(csv_path_str),
            )
            self._refresh_optimizer_tab()
            self.opt_status_var.set(f"Optimizer loaded. {len(self._optimizer.log.records)} iteration(s) on record.")
        except Exception as exc:
            messagebox.showerror("Optimizer Error", f"{exc}\n\n{traceback.format_exc()}")

    def _refresh_optimizer_tab(self) -> None:
        if self._optimizer is None:
            return

        # Surrogate status
        cov = self._optimizer.surrogate.coverage_report()
        if cov["n_data_points"] == 0:
            surr_text = "Surrogate not yet fitted. Generate a layout and ingest results first."
        else:
            r2_str = f"{cov['r2_score']:.4f}" if cov["r2_score"] is not None else "N/A"
            unc_str = f"{cov['mean_uncertainty_deg']:.4f}°" if cov["mean_uncertainty_deg"] is not None else "N/A"
            rng = cov["offset_range_m"]
            rng_str = f"[{rng[0]*1e3:.3g}, {rng[1]*1e3:.3g}] mm" if rng else "N/A"
            surr_text = (
                f"Data points: {cov['n_data_points']}  |  R²: {r2_str}  |  "
                f"Mean uncertainty: ±{unc_str}  |  Offset range: {rng_str}"
            )
        self.opt_surrogate_var.set(surr_text)

        # Suggestions
        suggestions = self._optimizer.suggest_next_offsets(n=5)
        if suggestions:
            self.opt_suggestions_var.set(
                "  ".join(f"{v*1e3:.3g}" for v in suggestions)
            )
        else:
            self.opt_suggestions_var.set("(fit surrogate first)")

        # Iteration treeview
        for item in self.iter_tree.get_children():
            self.iter_tree.delete(item)
        for rec in self._optimizer.log.records:
            score_str = f"{rec.layout_score:.2f}" if rec.layout_score is not None else "N/A"
            r2_s = f"{rec.surrogate_r2:.4f}" if rec.surrogate_r2 is not None else "N/A"
            unc_s = f"±{rec.uncertainty_mean_deg:.4f}" if rec.uncertainty_mean_deg is not None else "N/A"
            clipped = sum(1 for o in rec.proposed_offsets_m if o > 0)  # proxy for non-zero offsets used
            self.iter_tree.insert("", tk.END, values=(
                rec.iteration_id,
                rec.timestamp,
                score_str,
                r2_s,
                unc_s,
                len(rec.proposed_offsets_m),
                f"{rec.ap_rad_m * 1e3:.3g}",
                f"{rec.gap_m * 1e3:.3g}",
            ))

        self._draw_convergence()

    def _draw_convergence(self) -> None:
        canvas = self.conv_canvas
        canvas.delete("all")
        if self._optimizer is None:
            return
        data = self._optimizer.get_convergence_data()
        ids = data["iteration_ids"]
        scores = [s for s in data["layout_scores"] if s is not None]
        if len(ids) < 1 or len(scores) < 1:
            canvas.create_text(10, 10, anchor="nw", text="No iterations yet.", fill="#888888")
            return

        w = canvas.winfo_width() or 400
        h = canvas.winfo_height() or 120
        pad_l, pad_r, pad_t, pad_b = 50, 20, 12, 28

        valid = [(i, s) for i, s in zip(ids, data["layout_scores"]) if s is not None]
        if len(valid) < 1:
            return
        x_vals = [v[0] for v in valid]
        y_vals = [v[1] for v in valid]
        x_min, x_max = min(x_vals), max(x_vals)
        y_min, y_max = 0.0, max(y_vals) * 1.1 if max(y_vals) > 0 else 1.0
        if x_max == x_min:
            x_max = x_min + 1

        def to_cx(x):
            return pad_l + (x - x_min) / (x_max - x_min) * (w - pad_l - pad_r)

        def to_cy(y):
            return pad_t + (1.0 - (y - y_min) / (y_max - y_min)) * (h - pad_t - pad_b)

        # Axes
        canvas.create_line(pad_l, pad_t, pad_l, h - pad_b, fill="#aaaaaa")
        canvas.create_line(pad_l, h - pad_b, w - pad_r, h - pad_b, fill="#aaaaaa")
        canvas.create_text(4, pad_t, anchor="nw", text=f"{y_max:.1f}", font=("TkDefaultFont", 7), fill="#666666")
        canvas.create_text(pad_l, h - pad_b + 4, anchor="n", text=str(x_min), font=("TkDefaultFont", 7), fill="#666666")
        canvas.create_text(w - pad_r, h - pad_b + 4, anchor="n", text=str(x_max), font=("TkDefaultFont", 7), fill="#666666")
        canvas.create_text((w // 2), h - 6, anchor="s", text="Iteration", font=("TkDefaultFont", 8), fill="#444444")
        canvas.create_text(10, h // 2, anchor="center", angle=90, text="Score", font=("TkDefaultFont", 8), fill="#444444")

        # Line + points
        pts = [(to_cx(x), to_cy(y)) for x, y in zip(x_vals, y_vals)]
        if len(pts) >= 2:
            flat = [coord for pt in pts for coord in pt]
            canvas.create_line(*flat, fill="#0055cc", width=2)
        for px, py in pts:
            canvas.create_oval(px - 3, py - 3, px + 3, py + 3, fill="#0055cc", outline="")

    def _ingest_results(self) -> None:
        if not _HAS_OPTIMIZER:
            messagebox.showinfo("Optimizer", "layout_optimizer module not available.")
            return
        if self._optimizer is None:
            messagebox.showinfo("Optimizer", "Load an optimizer log first (Optimizer tab).")
            return
        source_str = filedialog.askopenfilename(
            title="Select runlog_compact.csv to ingest",
            filetypes=[("CSV", "*.csv"), ("All Files", "*.*")],
        )
        if not source_str:
            return
        try:
            ap_rad_m = float(self.radius_var.get()) * 1e-3
            gap_m = float(self.gap_var.get()) * 1e-3
            max_loss = float(self.max_loss_var.get())
            r_tol = float(self.radius_tol_var.get()) * 1e-3
            g_tol = float(self.gap_tol_var.get()) * 1e-3
            summary = self._optimizer.ingest_results(
                source=Path(source_str),
                ap_rad_m=ap_rad_m,
                gap_m=gap_m,
                max_loss_frac=max_loss,
                radius_tol=r_tol,
                gap_tol=g_tol,
            )
            msg = (
                f"Ingest complete: +{summary['n_new_rows']} rows "
                f"({summary['n_rows_before']} → {summary['n_rows_after']})"
            )
            r2 = summary.get("surrogate_r2")
            if r2 is not None:
                msg += f"  |  R²: {r2:.4f}"
            unc = summary.get("uncertainty_mean_deg")
            if unc is not None:
                msg += f"  |  Uncertainty: ±{unc:.4f}°"
            self.opt_status_var.set(msg)
            # Update CSV var to point at the enriched file
            self.csv_var.set(str(self._optimizer.csv_path))
            self._refresh_optimizer_tab()
            # Re-run solver with updated data and record new iteration
            if summary["n_new_rows"] > 0 and self.last_solution is not None:
                self._run_solver()
        except Exception as exc:
            messagebox.showerror("Ingest Error", f"{exc}\n\n{traceback.format_exc()}")

    def _export_slurm_sweep(self) -> None:
        if not _HAS_OPTIMIZER:
            messagebox.showinfo("Optimizer", "layout_optimizer module not available.")
            return
        if self._optimizer is None:
            messagebox.showinfo("Optimizer", "Load an optimizer log first (Optimizer tab).")
            return
        if self.last_solution is None:
            messagebox.showinfo("Optimizer", "Generate a layout first.")
            return

        # Dialog for sweep parameters
        dlg = tk.Toplevel(self)
        dlg.title("Export SLURM Sweep Wrapper")
        dlg.resizable(False, False)
        dlg.grab_set()

        orig_slurm_var = tk.StringVar(value=str(Path("orchestrate_pipeline.slurm")))
        delta_var = tk.StringVar(value="0.05")   # mm
        steps_var = tk.StringVar(value="3")
        out_var = tk.StringVar(value=str(Path("results") / "optimizer_sweep.sh"))
        fast_mode_var = tk.BooleanVar(value=False)

        frm = ttk.Frame(dlg, padding=12)
        frm.grid()
        frm.columnconfigure(1, weight=1)

        warn = ttk.Label(
            frm,
            text=(
                "\u26a0  SLURM Sweep exports a full beamlet layout (APERTURE_PAIRS_M is set).\n"
                "   Output data is useful for validating beam loss/transmission but\n"
                "   CANNOT train the steering surrogate — ACCEL_OFF_Y_M in the CSV\n"
                "   will record the aperture position, not the swept offset value.\n"
                "   Use  Export Grid Sweep  to generate steering-characterisation data."
            ),
            foreground="dark orange",
            justify="left",
            wraplength=500,
        )
        warn.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))

        ttk.Label(frm, text="orchestrate_pipeline.slurm:").grid(row=1, column=0, sticky="w")
        ttk.Entry(frm, textvariable=orig_slurm_var, width=48).grid(row=1, column=1, sticky="ew", padx=4)
        def _browse_orig():
            p = filedialog.askopenfilename(filetypes=[("SLURM", "*.slurm"), ("All Files", "*.*")],
                                           title="Select orchestrate_pipeline.slurm")
            if p:
                orig_slurm_var.set(p)
        ttk.Button(frm, text="Browse", command=_browse_orig).grid(row=1, column=2)

        ttk.Label(frm, text="Neighborhood delta (mm):").grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(frm, textvariable=delta_var, width=12).grid(row=2, column=1, sticky="w", padx=4, pady=(6, 0))
        ttk.Label(frm, text="Steps each side:").grid(row=3, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frm, textvariable=steps_var, width=12).grid(row=3, column=1, sticky="w", padx=4)
        ttk.Label(frm, text="Output wrapper (.sh):").grid(row=4, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frm, textvariable=out_var, width=48).grid(row=4, column=1, sticky="ew", padx=4)

        def _browse_out():
            p = filedialog.asksaveasfilename(defaultextension=".sh",
                                             filetypes=[("Shell script", "*.sh"), ("All Files", "*.*")])
            if p:
                out_var.set(p)
        ttk.Button(frm, text="Browse", command=_browse_out).grid(row=4, column=2)

        ttk.Checkbutton(
            frm,
            text="Fast mode  (kill vis outputs: PNG/density/profile;  pin X_RIGHT_M=0.02 m)",
            variable=fast_mode_var,
        ).grid(row=5, column=0, columnspan=3, sticky="w", pady=(8, 0))

        info = ttk.Label(frm,
            text="The wrapper patches only the sweep-list lines via sed.\n"
                 "The original .slurm file is never modified.",
            foreground="#555555", justify="left")
        info.grid(row=6, column=0, columnspan=3, sticky="w", pady=(2, 0))

        def do_export():
            try:
                orig_slurm = Path(orig_slurm_var.get().strip())
                delta_m = float(delta_var.get()) * 1e-3
                n_steps = int(steps_var.get())
                out_path = Path(out_var.get())
                sol = self.last_solution
                ap_rad_m = float(self.radius_var.get()) * 1e-3
                gap_m = float(self.gap_var.get()) * 1e-3
                latest = self._optimizer.log.latest()
                iter_tag = f"iter{latest.iteration_id}" if latest else "iterN"
                new_offsets = self._optimizer.generate_slurm_sweep(
                    out_path=out_path,
                    orig_slurm_path=orig_slurm,
                    proposed_offsets_m=latest.proposed_offsets_m if latest else [],
                    ap_rad_m=ap_rad_m,
                    gap_m=gap_m,
                    aperture_pairs_m=sol.summary.aperture_pairs_m or "",
                    delta_m=delta_m,
                    n_steps=n_steps,
                    run_prefix=f"optimizer_{iter_tag}",
                    iteration_id=latest.iteration_id if latest else None,
                    fast_mode=fast_mode_var.get(),
                )
                dlg.destroy()
                self.opt_status_var.set(
                    f"Sweep wrapper exported: {len(new_offsets)} new offset(s) → {out_path}\n"
                    f"Run on HPC with:  bash {out_path.name}"
                )
                self._refresh_optimizer_tab()
            except Exception as exc:
                messagebox.showerror("Export Error", f"{exc}\n\n{traceback.format_exc()}")

        ttk.Button(frm, text="Export", command=do_export).grid(row=7, column=0, columnspan=3, pady=(10, 0))

    def _export_grid_sweep(self) -> None:
        if not _HAS_OPTIMIZER:
            messagebox.showinfo("Optimizer", "layout_optimizer module not available.")
            return
        if self._optimizer is None:
            messagebox.showinfo("Optimizer", "Load an optimizer log first (Optimizer tab).")
            return

        dlg = tk.Toplevel(self)
        dlg.title("Export 2D Grid Sweep")
        dlg.resizable(False, False)
        dlg.grab_set()

        orig_slurm_var = tk.StringVar(value=str(Path("orchestrate_pipeline.slurm")))
        gaps_var       = tk.StringVar(value="0.5, 1.0, 1.5, 2.0")   # mm
        step_var       = tk.StringVar(value="0.2")                    # mm
        ap_rad_var     = tk.StringVar(value=self.radius_var.get())    # mm
        out_var        = tk.StringVar(value=str(Path("results") / "grid_sweep.sh"))
        fast_mode_var  = tk.BooleanVar(value=True)

        frm = ttk.Frame(dlg, padding=12)
        frm.grid()
        frm.columnconfigure(1, weight=1)

        ttk.Label(frm, text="orchestrate_pipeline.slurm:").grid(row=0, column=0, sticky="w")
        ttk.Entry(frm, textvariable=orig_slurm_var, width=48).grid(row=0, column=1, sticky="ew", padx=4)
        def _browse_orig():
            p = filedialog.askopenfilename(filetypes=[("SLURM", "*.slurm"), ("All Files", "*.*")],
                                           title="Select orchestrate_pipeline.slurm")
            if p:
                orig_slurm_var.set(p)
        ttk.Button(frm, text="Browse", command=_browse_orig).grid(row=0, column=2)

        ttk.Label(frm, text="Gap values (mm, comma-separated):").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(frm, textvariable=gaps_var, width=32).grid(row=1, column=1, sticky="w", padx=4, pady=(6, 0))

        ttk.Label(frm, text="Offset step (mm):").grid(row=2, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frm, textvariable=step_var, width=12).grid(row=2, column=1, sticky="w", padx=4)

        ttk.Label(frm, text="Aperture radius (mm):").grid(row=3, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frm, textvariable=ap_rad_var, width=12).grid(row=3, column=1, sticky="w", padx=4)

        ttk.Label(frm, text="Output wrapper (.sh):").grid(row=4, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frm, textvariable=out_var, width=48).grid(row=4, column=1, sticky="ew", padx=4)
        def _browse_out():
            p = filedialog.asksaveasfilename(defaultextension=".sh",
                                             filetypes=[("Shell script", "*.sh"), ("All Files", "*.*")])
            if p:
                out_var.set(p)
        ttk.Button(frm, text="Browse", command=_browse_out).grid(row=4, column=2)

        ttk.Checkbutton(
            frm,
            text="Fast mode  (disable vis outputs, pin X_RIGHT_M=0.02 m)",
            variable=fast_mode_var,
        ).grid(row=5, column=0, columnspan=3, sticky="w", pady=(8, 0))

        ttk.Label(
            frm,
            text="Sets APERTURE_PAIRS_M=\"\" so the sim uses ACCEL_OFF_Y_M directly\n"
                 "(single-aperture characterisation mode — prevents the offset=0 bug).",
            foreground="#555555", justify="left",
        ).grid(row=6, column=0, columnspan=3, sticky="w", pady=(2, 0))

        def do_grid_export():
            try:
                orig_slurm  = Path(orig_slurm_var.get().strip())
                gap_values_m = [float(g.strip()) * 1e-3 for g in gaps_var.get().split(",") if g.strip()]
                step_m      = float(step_var.get()) * 1e-3
                ap_rad_m    = float(ap_rad_var.get()) * 1e-3
                out_path    = Path(out_var.get())
                planned = self._optimizer.generate_grid_sweep(
                    out_path=out_path,
                    orig_slurm_path=orig_slurm,
                    gap_values_m=gap_values_m,
                    ap_rad_m=ap_rad_m,
                    offset_step_m=step_m,
                    run_prefix="grid_sweep",
                    fast_mode=fast_mode_var.get(),
                )
                n_gaps    = len({g for g, _ in planned})
                n_offsets = len({o for _, o in planned})
                dlg.destroy()
                self.opt_status_var.set(
                    f"Grid sweep exported: {n_gaps} gaps × {n_offsets} offsets = {len(planned)} runs → {out_path}\n"
                    f"Run on HPC with:  bash {out_path.name}"
                )
            except Exception as exc:
                messagebox.showerror("Grid Sweep Error", f"{exc}\n\n{traceback.format_exc()}")

        ttk.Button(frm, text="Export", command=do_grid_export).grid(row=7, column=0, columnspan=3, pady=(10, 0))

    def _bind_preview_updates(self) -> None:
        for var in [self.radius_var, self.packing_var, self.target_mode_var, self.screen_pitch_var, self.plate_radius_var, self.target_radius_var, self.count_var]:
            var.trace_add("write", lambda *_: self._update_capacity_preview())

    def _update_capacity_preview(self) -> None:
        try:
            radius_mm = float(self.radius_var.get())
            pitch_mm = float(self.screen_pitch_var.get())
            plate_radius_mm = float(self.plate_radius_var.get())
            target_radius_mm = self._float_or_none(self.target_radius_var.get())
            if target_radius_mm is None:
                target_radius_mm = plate_radius_mm
            count = self._int_or_none(self.count_var.get())
            mode = self.packing_var.get().strip().lower()
            screen_all = solver._generate_points(mode, None, pitch_mm * 1.0e-3, plate_radius_mm * 1.0e-3)
            cap = solver.estimate_beamlet_count(radius_mm * 1.0e-3)
            if count is None:
                txt = (
                    f"Available grid positions: {len(screen_all)}\n"
                    f"Auto-count will use all available grid positions.\n"
                    f"Target radius {target_radius_mm:g} mm does not limit count; target points will compress/overlap as needed.\n"
                    f"Estimated full-plate capacity for aperture radius {radius_mm:g} mm: {cap:.0f}"
                )
            else:
                txt = (
                    f"Requested count: {count}\n"
                    f"Available grid positions: {len(screen_all)}\n"
                    f"Target radius {target_radius_mm:g} mm does not limit count; target points will compress/overlap as needed.\n"
                    f"Estimated full-plate capacity for aperture radius {radius_mm:g} mm: {cap:.0f}"
                )
        except Exception:
            txt = "Enter radius, pitch, plate radius, and optional target radius to preview available grid positions."
        self.capacity_var.set(txt)


    def _set_preview_status(self, text: str) -> None:
        self.preview_status_var.set(text)

    def _preview_world_bounds_mm(self, sol, mode: str = "overlap") -> tuple[float, float, float, float]:
        xs = []
        ys = []
        if mode == "target":
            cx_mm = float(sol.summary.target_center_x_m) * 1e3
            cy_mm = float(sol.summary.target_center_y_m) * 1e3
            tr_mm = float(sol.summary.target_radius_m) * 1e3
            xs.extend([cx_mm - tr_mm, cx_mm + tr_mm])
            ys.extend([cy_mm - tr_mm, cy_mm + tr_mm])
            for b in sol.beamlets:
                xs.append(b.target_x_m * 1e3)
                ys.append(b.target_y_m * 1e3)
        elif mode == "side":
            gap_mm = float(sol.summary.gap_m) * 1e3
            target_plane_mm = float(sol.summary.target_plane_x_m) * 1e3
            cx_mm = float(sol.summary.target_center_x_m) * 1e3
            cy_mm = float(sol.summary.target_center_y_m) * 1e3
            xs.extend([0.0, gap_mm, target_plane_mm])
            ys.append(float(sol.summary.target_radius_m) * 1e3)
            for b in sol.beamlets:
                ys.extend([
                    math.hypot(b.screen_x_m * 1e3 - cx_mm, b.screen_y_m * 1e3 - cy_mm),
                    math.hypot(b.accel_x_m * 1e3 - cx_mm, b.accel_y_m * 1e3 - cy_mm),
                    math.hypot(b.target_x_m * 1e3 - cx_mm, b.target_y_m * 1e3 - cy_mm),
                    math.hypot(b.desired_target_x_m * 1e3 - cx_mm, b.desired_target_y_m * 1e3 - cy_mm),
                ])
        else:
            for b in sol.beamlets:
                xs.extend([b.screen_x_m * 1e3, b.accel_x_m * 1e3, b.target_x_m * 1e3, b.desired_target_x_m * 1e3])
                ys.extend([b.screen_y_m * 1e3, b.accel_y_m * 1e3, b.target_y_m * 1e3, b.desired_target_y_m * 1e3])
        pad = max(float(sol.summary.aperture_radius_m) * 1e3 * 2.5, float(sol.summary.screen_pitch_m) * 1e3 * 0.5, 1.0)
        x0, x1 = (min(xs) - pad, max(xs) + pad) if xs else (-10.0, 10.0)
        y0, y1 = (min(ys) - pad, max(ys) + pad) if ys else (-10.0, 10.0)
        if mode == "side":
            y0 = min(-0.5, y0)
        if abs(x1 - x0) < 1e-9:
            x0 -= 1.0
            x1 += 1.0
        if abs(y1 - y0) < 1e-9:
            y0 -= 1.0
            y1 += 1.0
        return x0, x1, y0, y1

    def _preview_margins(self) -> tuple[float, float, float, float]:
        return 64.0, 22.0, 22.0, 50.0

    def _world_to_canvas(self, x_mm: float, y_mm: float, bounds, width: float, height: float):
        x0, x1, y0, y1 = bounds
        left, right, top, bottom = self._preview_margins()
        avail_w = max(width - left - right, 1.0)
        avail_h = max(height - top - bottom, 1.0)
        sx = avail_w / max(x1 - x0, 1e-9)
        sy = avail_h / max(y1 - y0, 1e-9)
        scale = min(sx, sy)
        used_w = scale * (x1 - x0)
        used_h = scale * (y1 - y0)
        ox = left + 0.5 * (avail_w - used_w)
        oy = top + 0.5 * (avail_h - used_h)
        cx = ox + (x_mm - x0) * scale
        cy = oy + (y1 - y_mm) * scale
        return cx, cy, scale

    def _nice_tick_step(self, span_mm: float, target_ticks: int = 6) -> float:
        span_mm = abs(float(span_mm))
        if span_mm <= 1.0e-9:
            return 1.0
        raw = span_mm / max(target_ticks, 1)
        mag = 10.0 ** math.floor(math.log10(raw))
        norm = raw / mag
        if norm <= 1.0:
            step = 1.0
        elif norm <= 2.0:
            step = 2.0
        elif norm <= 5.0:
            step = 5.0
        else:
            step = 10.0
        return step * mag

    def _preview_axis_labels(self, mode: str) -> tuple[str, str]:
        if mode == "side":
            return "Axial Distance (mm)", "Radial Offset From Target Center (mm)"
        return "X (mm)", "Y (mm)"

    def _draw_measurement_axes(self, canvas: tk.Canvas, bounds, width: float, height: float, mode: str = "overlap") -> None:
        x0, x1, y0, y1 = bounds
        left, right, top, bottom = self._preview_margins()
        x_axis_y = height - bottom
        y_axis_x = left
        canvas.create_line(left, x_axis_y, width - right, x_axis_y, fill="#666666", width=1.5)
        canvas.create_line(y_axis_x, top, y_axis_x, height - bottom, fill="#666666", width=1.5)

        x_step = self._nice_tick_step(x1 - x0)
        y_step = self._nice_tick_step(y1 - y0)

        x_tick = math.ceil(x0 / x_step) * x_step
        while x_tick <= x1 + 1.0e-9:
            cx, _, _ = self._world_to_canvas(x_tick, y0, bounds, width, height)
            if left - 1.0 <= cx <= width - right + 1.0:
                canvas.create_line(cx, x_axis_y, cx, x_axis_y + 6, fill="#666666", width=1.0)
                canvas.create_text(cx, x_axis_y + 18, text=f"{x_tick:g}", fill="#444444", font=("Segoe UI", 9))
            x_tick += x_step

        y_tick = math.ceil(y0 / y_step) * y_step
        while y_tick <= y1 + 1.0e-9:
            _, cy, _ = self._world_to_canvas(x0, y_tick, bounds, width, height)
            if top - 1.0 <= cy <= height - bottom + 1.0:
                canvas.create_line(y_axis_x - 6, cy, y_axis_x, cy, fill="#666666", width=1.0)
                canvas.create_text(y_axis_x - 10, cy, text=f"{y_tick:g}", fill="#444444", font=("Segoe UI", 9), anchor="e")
            y_tick += y_step

        x_label, y_label = self._preview_axis_labels(mode)
        canvas.create_text((left + width - right) * 0.5, height - 18, text=x_label, fill="#333333", font=("Segoe UI", 10, "bold"))
        canvas.create_text(20, (top + height - bottom) * 0.5, text=y_label, fill="#333333", font=("Segoe UI", 10, "bold"), angle=90)

        if x0 <= 0.0 <= x1:
            cx0, cy0, _ = self._world_to_canvas(0.0, y0, bounds, width, height)
            cx1, cy1, _ = self._world_to_canvas(0.0, y1, bounds, width, height)
            canvas.create_line(cx0, cy0, cx1, cy1, fill="#dddddd", dash=(3, 2))
        if y0 <= 0.0 <= y1:
            cx0, cy0, _ = self._world_to_canvas(x0, 0.0, bounds, width, height)
            cx1, cy1, _ = self._world_to_canvas(x1, 0.0, bounds, width, height)
            canvas.create_line(cx0, cy0, cx1, cy1, fill="#dddddd", dash=(3, 2))

    def _bind_hover(self, canvas: tk.Canvas, item_ids, text: str) -> None:
        for item in item_ids:
            canvas.tag_bind(item, "<Enter>", lambda _e, t=text: self._set_preview_status(t))
            canvas.tag_bind(item, "<Leave>", lambda _e: self._set_preview_status("Hover over a beamlet in the preview to inspect coordinates."))

    def _draw_side_preview(self, canvas: tk.Canvas, sol, bounds, w: float, h: float) -> None:
        cx_mm = float(sol.summary.target_center_x_m) * 1e3
        cy_mm = float(sol.summary.target_center_y_m) * 1e3
        screen_x_mm = 0.0
        accel_x_mm = float(sol.summary.gap_m) * 1e3
        target_x_mm = float(sol.summary.target_plane_x_m) * 1e3
        for b in sol.beamlets:
            screen_r = math.hypot(b.screen_x_m * 1e3 - cx_mm, b.screen_y_m * 1e3 - cy_mm)
            accel_r = math.hypot(b.accel_x_m * 1e3 - cx_mm, b.accel_y_m * 1e3 - cy_mm)
            target_r = math.hypot(b.target_x_m * 1e3 - cx_mm, b.target_y_m * 1e3 - cy_mm)
            desired_r = math.hypot(b.desired_target_x_m * 1e3 - cx_mm, b.desired_target_y_m * 1e3 - cy_mm)
            sx, sy, _ = self._world_to_canvas(screen_x_mm, screen_r, bounds, w, h)
            ax, ay, _ = self._world_to_canvas(accel_x_mm, accel_r, bounds, w, h)
            tx, ty, _ = self._world_to_canvas(target_x_mm, target_r, bounds, w, h)
            dx, dy, _ = self._world_to_canvas(target_x_mm, desired_r, bounds, w, h)
            feasible = b.target_constraint_feasible
            line_color = "#2e8b57" if feasible else "#c0392b"
            point_outline = "#1f6f45" if feasible else "#7f1d1d"
            hover = (
                f"Beamlet {b.index}: screen radius={screen_r:.4g} mm; "
                f"accel radius={accel_r:.4g} mm; landing radius={target_r:.4g} mm; "
                f"requested radius={desired_r:.4g} mm; theta used={b.steer_angle_used_deg:.4g} deg; "
                f"shortfall={b.target_shortfall_m*1e3:.4g} mm"
            )
            items = [
                canvas.create_line(sx, sy, ax, ay, fill="#666666", width=1.5),
                canvas.create_line(ax, ay, tx, ty, fill=line_color, width=2.0),
                canvas.create_oval(sx-3, sy-3, sx+3, sy+3, fill="#5b7cfa", outline="#1f3fb5"),
                canvas.create_oval(ax-3, ay-3, ax+3, ay+3, fill="#ff8a7a", outline="#b53020"),
                canvas.create_oval(tx-4, ty-4, tx+4, ty+4, fill=line_color, outline=point_outline),
                canvas.create_text(tx + 8, ty - 8, text=str(b.index), fill=point_outline, font=("Segoe UI", 9, "bold"), anchor="sw"),
            ]
            if abs(desired_r - target_r) > 1e-6:
                items.append(canvas.create_oval(dx-3, dy-3, dx+3, dy+3, outline="#cc8f00", width=1.0))
            self._bind_hover(canvas, items, hover)

        for xpos_mm, label, color in [(screen_x_mm, "screen", "#1f3fb5"), (accel_x_mm, "accelerator", "#b53020"), (target_x_mm, "target", "#1f6f45")]:
            cx, cy0, _ = self._world_to_canvas(xpos_mm, bounds[2], bounds, w, h)
            _, cy1, _ = self._world_to_canvas(xpos_mm, bounds[3], bounds, w, h)
            canvas.create_line(cx, cy0, cx, cy1, fill=color, dash=(3, 2))
            canvas.create_text(cx, 12, text=label, fill=color, font=("Segoe UI", 9, "bold"), anchor="n")

        _, target_rad_y, _ = self._world_to_canvas(bounds[0], float(sol.summary.target_radius_m) * 1e3, bounds, w, h)
        x_left, _, _ = self._world_to_canvas(bounds[0], 0.0, bounds, w, h)
        x_right, _, _ = self._world_to_canvas(bounds[1], 0.0, bounds, w, h)
        canvas.create_line(x_left, target_rad_y, x_right, target_rad_y, fill="#2e8b57", width=2, dash=(4, 2))
        canvas.create_text(12, 12, anchor="nw", text="Side View preview", fill="#333333", font=("Segoe UI", 10, "bold"))
        canvas.create_text(12, 30, anchor="nw", text="Green line/dot: landing inside target\nRed line/dot: landing outside target\nAmber ring: requested target radius", fill="#2f2f2f", font=("Segoe UI", 9))

    def _draw_preview(self, canvas: tk.Canvas, mode: str) -> None:
        canvas.delete("all")
        sol = self.last_solution
        if sol is None:
            return
        w = max(canvas.winfo_width(), 200)
        h = max(canvas.winfo_height(), 200)
        bounds = self._preview_world_bounds_mm(sol, mode)
        self._draw_measurement_axes(canvas, bounds, w, h, mode)

        if mode == "side":
            self._draw_side_preview(canvas, sol, bounds, w, h)
            return

        for b in sol.beamlets:
            sx, sy, scale = self._world_to_canvas(b.screen_x_m * 1e3, b.screen_y_m * 1e3, bounds, w, h)
            ax, ay, scale = self._world_to_canvas(b.accel_x_m * 1e3, b.accel_y_m * 1e3, bounds, w, h)
            tx, ty, scale = self._world_to_canvas(b.target_x_m * 1e3, b.target_y_m * 1e3, bounds, w, h)
            dx, dy, scale = self._world_to_canvas(b.desired_target_x_m * 1e3, b.desired_target_y_m * 1e3, bounds, w, h)
            rad_px = max(sol.summary.aperture_radius_m * 1e3 * scale, 2.0)

            hover = (
                f"Beamlet {b.index}: screen=({b.screen_x_m*1e3:.4g}, {b.screen_y_m*1e3:.4g}) mm; "
                f"accel=({b.accel_x_m*1e3:.4g}, {b.accel_y_m*1e3:.4g}) mm; "
                f"landing=({b.target_x_m*1e3:.4g}, {b.target_y_m*1e3:.4g}) mm; "
                f"offset=({b.accel_offset_x_m*1e3:.4g}, {b.accel_offset_y_m*1e3:.4g}) mm; "
                f"theta used={b.steer_angle_used_deg:.4g} deg; hits_target={1 if b.hits_target else 0}; "
                f"shortfall={b.target_shortfall_m*1e3:.4g} mm"
            )

            items = []
            if mode in ("screen", "overlap"):
                items.append(canvas.create_oval(sx-rad_px, sy-rad_px, sx+rad_px, sy+rad_px, fill="#5b7cfa", outline="#1f3fb5", width=1.5))
                items.append(canvas.create_text(sx, sy-rad_px-8, text=str(b.index), fill="#1f3fb5", font=("Segoe UI", 9, "bold")))
            if mode in ("accel", "overlap"):
                items.append(canvas.create_oval(ax-rad_px, ay-rad_px, ax+rad_px, ay+rad_px, fill="#ff8a7a", outline="#b53020", width=1.5))
                items.append(canvas.create_text(ax, ay-rad_px-8, text=str(b.index), fill="#b53020", font=("Segoe UI", 9, "bold")))
            if mode == "overlap":
                items.append(canvas.create_line(sx, sy, ax, ay, fill="#666666", width=1.5))
                items.append(canvas.create_line(ax, ay, tx, ty, fill="#999999", width=1.0, dash=(4, 2)))
                landing_fill = "#2e8b57" if b.target_constraint_feasible else "#c0392b"
                landing_outline = "" if b.target_constraint_feasible else "#7f1d1d"
                items.append(canvas.create_oval(tx-2, ty-2, tx+2, ty+2, fill=landing_fill, outline=landing_outline))
                if abs(dx - tx) > 1e-6 or abs(dy - ty) > 1e-6:
                    items.append(canvas.create_oval(dx-2, dy-2, dx+2, dy+2, outline="#cc8f00", width=1))
            if mode == "target":
                landing_fill = "#2e8b57" if b.target_constraint_feasible else "#c0392b"
                landing_outline = "#1f6f45" if b.target_constraint_feasible else "#7f1d1d"
                items.append(canvas.create_oval(tx-3, ty-3, tx+3, ty+3, fill=landing_fill, outline=landing_outline, width=1.0))
                items.append(canvas.create_text(tx + 8, ty - 8, text=str(b.index), fill=landing_outline, font=("Segoe UI", 9, "bold"), anchor="sw"))
            self._bind_hover(canvas, items, hover)

        if mode in ("overlap", "target"):
            cx, cy, scale = self._world_to_canvas(float(sol.summary.target_center_x_m) * 1e3, float(sol.summary.target_center_y_m) * 1e3, bounds, w, h)
            tr_px = max(float(sol.summary.target_radius_m) * 1e3 * scale, 1.0)
            canvas.create_oval(cx-tr_px, cy-tr_px, cx+tr_px, cy+tr_px, outline="#2e8b57", width=2, dash=(3, 2))
            canvas.create_line(cx-5, cy, cx+5, cy, fill="#2e8b57", width=2)
            canvas.create_line(cx, cy-5, cx, cy+5, fill="#2e8b57", width=2)
            canvas.create_text(12, 12, anchor="nw", text=f"{mode.title()} preview", fill="#333333", font=("Segoe UI", 10, "bold"))
            if mode == "overlap":
                canvas.create_text(12, 30, anchor="nw", text="Green dashed circle: target acceptance radius\nGreen dot: actual landing point\nAmber ring: requested landing point", fill="#2f2f2f", font=("Segoe UI", 9))
            else:
                canvas.create_text(12, 30, anchor="nw", text="Green dashed circle: target acceptance radius\nGreen dots: beamlet landing points\nRed dots: landing points outside target", fill="#2f2f2f", font=("Segoe UI", 9))
        else:
            canvas.create_text(12, 12, anchor="nw", text=f"{mode.title()} preview", fill="#333333", font=("Segoe UI", 10, "bold"))

    def _redraw_previews(self) -> None:
        for key, canvas in getattr(self, "preview_canvases", {}).items():
            self._draw_preview(canvas, key)


    def _browse_csv(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All Files", "*.*")])
        if p:
            self.csv_var.set(p)

    def _browse_json(self) -> None:
        p = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json"), ("All Files", "*.*")])
        if p:
            self.out_json_var.set(p)

    def _browse_prefix(self) -> None:
        p = filedialog.asksaveasfilename(defaultextension="", filetypes=[("All Files", "*.*")])
        if p:
            self.plot_prefix_var.set(p)

    def _float_or_none(self, text: str):
        t = text.strip()
        if not t:
            return None
        return float(t)

    def _int_or_none(self, text: str):
        t = text.strip()
        if not t:
            return None
        return int(float(t))

    def _run_solver(self) -> None:
        try:
            df = solver._load_csv(Path(self.csv_var.get()))
            target_radius_mm = self._float_or_none(self.target_radius_var.get())
            plate_radius_mm = float(self.plate_radius_var.get())
            if target_radius_mm is None:
                target_radius_mm = plate_radius_mm
            sol = solver.solve_layout(
                df=df,
                radius_m=float(self.radius_var.get()) * 1.0e-3,
                gap_m=float(self.gap_var.get()) * 1.0e-3,
                max_loss_frac=float(self.max_loss_var.get()),
                packing_mode=self.packing_var.get().strip().lower(),
                target_mode=self.target_mode_var.get().strip().lower(),
                screen_pitch_m=float(self.screen_pitch_var.get()) * 1.0e-3,
                plate_radius_m=plate_radius_mm * 1.0e-3,
                target_radius_m=target_radius_mm * 1.0e-3,
                target_plane_x_m=float(self.target_plane_var.get()) * 1.0e-3,
                steer_origin_x_m=float(self.steer_origin_var.get()) * 1.0e-3,
                screen_center_x_m=float(self.screen_center_x_var.get()) * 1.0e-3,
                screen_center_y_m=float(self.screen_center_y_var.get()) * 1.0e-3,
                target_center_x_m=float(self.target_center_x_var.get()) * 1.0e-3,
                target_center_y_m=float(self.target_center_y_var.get()) * 1.0e-3,
                count=self._int_or_none(self.count_var.get()),
                extra_filter=self.extra_filter_var.get().strip() or None,
                radius_tol=float(self.radius_tol_var.get()) * 1.0e-3,
                gap_tol=float(self.gap_tol_var.get()) * 1.0e-3,
                clip_to_safe_limit=self.clip_var.get(),
            )

            payload = {
                "summary": asdict(sol.summary),
                "safe_curve": [{"offset_m": x, "steer_deg": y} for x, y in sol.safe_curve],
                "beamlets": [asdict(b) for b in sol.beamlets],
            }
            out_json = Path(self.out_json_var.get())
            out_json.parent.mkdir(parents=True, exist_ok=True)
            out_json.write_text(json.dumps(payload, indent=2))

            plot_prefix = Path(self.plot_prefix_var.get())
            plot_paths = solver.write_layout_plots(sol, plot_prefix)
            self.last_plot_paths = {
                "overlap": plot_paths[0],
                "screen": plot_paths[1],
                "accel": plot_paths[2],
                "side": plot_paths[3],
                "json": out_json,
            }

            self.last_solution = sol
            self._render_summary(sol, payload)
            self._render_beamlets(sol)
            self._set_pairs_ready(sol.summary.aperture_pairs_m)
            # Auto-record iteration if optimizer is active
            if _HAS_OPTIMIZER and self._optimizer is not None:
                try:
                    self._optimizer.refit_surrogate(
                        float(self.radius_var.get()) * 1e-3,
                        float(self.gap_var.get()) * 1e-3,
                        float(self.max_loss_var.get()),
                        float(self.radius_tol_var.get()) * 1e-3,
                        float(self.gap_tol_var.get()) * 1e-3,
                    )
                    self._optimizer.record_solution(
                        sol,
                        float(self.radius_var.get()) * 1e-3,
                        float(self.gap_var.get()) * 1e-3,
                    )
                    self._refresh_optimizer_tab()
                except Exception:
                    pass  # never crash the main solver flow due to optimizer
            if sol.summary.infeasible_beamlet_count > 0:
                self._set_preview_status(
                    f"Warning: {sol.summary.infeasible_beamlet_count} beamlets cannot reach the target area within the safe steering limit."
                )
            else:
                self._set_preview_status("Hover over a beamlet in the preview to inspect coordinates.")
            self._redraw_previews()

            if self.auto_open_var.get():
                self._open_plot("screen")
                self._open_plot("accel")
                self._open_plot("overlap")
                self._open_plot("side")
        except Exception as exc:
            messagebox.showerror("Solver Error", f"{exc}\n\n{traceback.format_exc()}")

    def _render_summary(self, sol, payload) -> None:
        self.summary_text.delete("1.0", tk.END)
        self.summary_text.insert(tk.END, json.dumps(payload["summary"], indent=2))
        self.summary_text.insert(tk.END, "\n\nGenerated Files:\n")
        for key in ["json", "screen", "accel", "overlap", "side"]:
            p = self.last_plot_paths.get(key)
            if p:
                self.summary_text.insert(tk.END, f"{key}: {p}\n")

    def _render_beamlets(self, sol) -> None:
        for item in self.tree.get_children():
            self.tree.delete(item)
        for b in sol.beamlets:
            self.tree.insert("", tk.END, values=(
                b.index,
                f"{b.screen_x_m * 1e3:.4g}",
                f"{b.screen_y_m * 1e3:.4g}",
                f"{b.accel_x_m * 1e3:.4g}",
                f"{b.accel_y_m * 1e3:.4g}",
                f"{b.accel_offset_x_m * 1e3:.4g}",
                f"{b.accel_offset_y_m * 1e3:.4g}",
                f"{b.accel_offset_mag_m * 1e3:.4g}",
                f"{b.target_x_m * 1e3:.4g}",
                f"{b.target_y_m * 1e3:.4g}",
                f"{b.steer_angle_req_deg:.4g}",
                f"{b.steer_angle_used_deg:.4g}",
                "1" if b.clipped_to_safe_limit else "0",
                "1" if b.hits_target else "0",
                "1" if b.target_constraint_feasible else "0",
                f"{b.target_shortfall_m * 1e3:.4g}",
            ))

    def _set_pairs_ready(self, text: str) -> None:
        sol = self.last_solution
        warning = ""
        if sol is not None and getattr(sol.summary, "infeasible_beamlet_count", 0) > 0:
            warning = f" Warning: {sol.summary.infeasible_beamlet_count} beamlets miss the target radius at the safe limit."
        if text:
            self.pairs_status_var.set("APERTURE_PAIRS_M ready. Click the copy button to place it on the clipboard." + warning)
        else:
            self.pairs_status_var.set("APERTURE_PAIRS_M not generated yet.")

    def _copy_pairs_to_clipboard(self) -> None:
        sol = self.last_solution
        if sol is None or not sol.summary.aperture_pairs_m:
            messagebox.showinfo("No APERTURE_PAIRS_M", "Generate a layout first.")
            return
        self.clipboard_clear()
        self.clipboard_append(sol.summary.aperture_pairs_m)
        self.update()
        self.pairs_status_var.set("APERTURE_PAIRS_M copied to clipboard.")

    def _open_plot(self, key: str) -> None:
        p = self.last_plot_paths.get(key)
        if not p:
            messagebox.showinfo("No Plot", f"No {key} plot has been generated yet.")
            return
        webbrowser.open(Path(p).resolve().as_uri())

    def _open_json(self) -> None:
        p = self.last_plot_paths.get("json")
        if not p:
            messagebox.showinfo("No JSON", "No JSON output has been generated yet.")
            return
        webbrowser.open(Path(p).resolve().as_uri())


    # ------------------------------------------------------------------
    # Curved Grid tab — class-level field definitions
    # ------------------------------------------------------------------

    # (key, label, unit_label, default)
    _GEOM_DEFS = [
        ("R_SCR",  "Screen R",         "mm",  "100"),
        ("T_SCR",  "Screen thickness", "mm",  "5"),
        ("T_ACC",  "Accel thickness",  "mm",  "5"),
        ("GAP",    "Gap after screen", "mm",  "6"),
        ("VS",     "Screen voltage",   "V",   "0.0"),
        ("VA",     "Accel voltage",    "V",   "-15000.0"),
        ("AP_RAD", "Bore radius",      "mm",  "1.5"),
    ]

    # (key, default, description)
    _ENV_DEFS = [
        # ── Run config (not sweep targets) ───────────────────────────────
        ("RUN_SOLVE",      "1",                    "Run Poisson/PIC solve (0 = geometry only)"),
        ("WRITE_PNG",      "1",                    "Write PNG outputs"),
        ("PNG_NAME",       "curved_grid_geom.png", "Geometry PNG filename"),
        ("RESULTS_DIR",    "results",              "Output directory name"),
        ("RUN_PREFIX",     "curved_grid_test",     "Run folder prefix"),
        ("ENABLE_IONS",    "1",                    "Enable ion particle tracing"),
        # ── Domain geometry ───────────────────────────────────────────────
        ("X_RIGHT_M",      "0.025",                "Right boundary of close-up view (m)"),
        ("X_RIGHT_PHYS_M", "0.55",                 "Physical right boundary / tube length (m)"),
        ("YBOX_M",         "0.090",                "Half-height of simulation box (m)"),
        ("TUBE_WALL_T_M",  "0.0002",               "Drift tube wall thickness (m)"),
        # ── Physics / solver ─────────────────────────────────────────────
        ("H",              "0.0001",               "Mesh cell size (m)"),
        ("ION_N",          "1e16",                 "Ion number density (m^-3)"),
        ("SC_FACTOR",      "0.0005",               "Space-charge compensation factor"),
    ]

    # Env vars that are run-config knobs, not physics — excluded from the
    # sweep axis dropdown (still editable in the Simulation Inputs tab).
    _SWEEP_ENV_EXCLUDE = frozenset({
        "RUN_SOLVE", "WRITE_PNG", "PNG_NAME",
        "RESULTS_DIR", "RUN_PREFIX", "ENABLE_IONS",
    })

    # ------------------------------------------------------------------
    # Curved Grid tab — builder
    # ------------------------------------------------------------------

    def _build_curved_tab(self) -> None:
        ctab = ttk.Frame(self.main_notebook, padding=8)
        ctab.columnconfigure(0, weight=1)
        ctab.rowconfigure(1, weight=1)
        self.main_notebook.add(ctab, text="Curved Grid")

        if not _HAS_CURVED:
            ttk.Label(
                ctab,
                text="curved_grid_pipeline.py not found. Place it alongside this script.",
                foreground="red",
            ).grid(row=0, column=0, sticky="w")
            return

        # Initialise StringVar dicts with defaults
        self.curved_geom_vars = {k: tk.StringVar(value=d) for k, _, _, d in self._GEOM_DEFS}
        self.curved_env_vars  = {k: tk.StringVar(value=d) for k, d, _ in self._ENV_DEFS}

        # ── Import section ──────────────────────────────────────────────────
        imp_frm = ttk.LabelFrame(ctab, text="Import Curved Case", padding=8)
        imp_frm.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        imp_frm.columnconfigure(1, weight=1)

        ttk.Label(imp_frm, text="Case file:").grid(row=0, column=0, sticky="w")
        ttk.Entry(imp_frm, textvariable=self.curved_case_var, width=52).grid(
            row=0, column=1, sticky="ew", padx=4)
        ttk.Button(imp_frm, text="Browse",
                   command=self._browse_curved_case).grid(row=0, column=2)
        ttk.Button(imp_frm, text="Load Case",
                   command=self._load_curved_case_file).grid(row=0, column=3, padx=(4, 0))
        ttk.Label(imp_frm, textvariable=self.curved_status_var,
                  foreground="#555555", wraplength=680, justify="left").grid(
            row=1, column=0, columnspan=4, sticky="w", pady=(4, 0))

        # ── Sub-notebook ────────────────────────────────────────────────────
        sub_nb = ttk.Notebook(ctab)
        sub_nb.grid(row=1, column=0, sticky="nsew", pady=(0, 6))

        # ── Tab 1 : Grid Geometry ───────────────────────────────────────────
        geom_tab = ttk.Frame(sub_nb, padding=8)
        geom_tab.columnconfigure(0, weight=1)
        sub_nb.add(geom_tab, text="Grid Geometry")

        info_frm = ttk.LabelFrame(geom_tab, text="Loaded Case Summary", padding=6)
        info_frm.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        info_frm.columnconfigure(0, weight=1)
        self.curved_info_text = tk.Text(
            info_frm, height=4, wrap="word", state="disabled",
            font=("Consolas", 9), background="#f8f8f8")
        self.curved_info_text.grid(row=0, column=0, sticky="ew")

        edit_frm = ttk.LabelFrame(geom_tab, text="Edit Parameters", padding=8)
        edit_frm.grid(row=1, column=0, sticky="ew")
        edit_frm.columnconfigure(1, weight=1)
        edit_frm.columnconfigure(4, weight=1)

        for idx, (key, label, unit, _) in enumerate(self._GEOM_DEFS):
            col_off = (idx % 2) * 3
            row = idx // 2
            ttk.Label(edit_frm, text=f"{label} ({unit}):").grid(
                row=row, column=col_off, sticky="w",
                padx=(0 if col_off == 0 else 20, 4), pady=3)
            ttk.Entry(edit_frm, textvariable=self.curved_geom_vars[key], width=14).grid(
                row=row, column=col_off + 1, sticky="w", pady=3)

        layout_frm = ttk.LabelFrame(geom_tab, text="Aperture Layout", padding=8)
        layout_frm.grid(row=2, column=0, sticky="ew", pady=(8, 0))

        ttk.Label(layout_frm, text="N apertures:").grid(row=0, column=0, sticky="w", padx=(0, 4))
        ttk.Spinbox(layout_frm, from_=1, to=50, textvariable=self.ap_n_var,
                    width=5).grid(row=0, column=1, sticky="w")
        ttk.Label(layout_frm, text="Offset (mm):").grid(row=0, column=2, sticky="w", padx=(16, 4))
        ttk.Entry(layout_frm, textvariable=self.ap_offset_var, width=8).grid(
            row=0, column=3, sticky="w")
        ttk.Label(layout_frm, text="Pitch (mm):").grid(row=0, column=4, sticky="w", padx=(12, 4))
        ttk.Entry(layout_frm, textvariable=self.ap_pitch_var, width=8).grid(
            row=0, column=5, sticky="w")
        ttk.Label(layout_frm, text="Steering arc (mm):").grid(
            row=0, column=6, sticky="w", padx=(12, 4))
        ttk.Entry(layout_frm, textvariable=self.ap_steering_var, width=8).grid(
            row=0, column=7, sticky="w")
        ttk.Button(layout_frm, text="Apply Layout",
                   command=self._apply_aperture_layout).grid(
            row=0, column=8, padx=(16, 0))
        ttk.Label(layout_frm, textvariable=self.ap_offsets_label_var,
                  foreground="#555555").grid(
            row=1, column=0, columnspan=9, sticky="w", pady=(4, 0))

        # ── Tab 2 : Simulation Inputs ───────────────────────────────────────
        env_tab = ttk.Frame(sub_nb, padding=8)
        env_tab.columnconfigure(0, weight=1)
        sub_nb.add(env_tab, text="Simulation Inputs")

        env_frm = ttk.LabelFrame(env_tab, text="Environment Variables", padding=8)
        env_frm.grid(row=0, column=0, sticky="ew")
        env_frm.columnconfigure(1, weight=1)
        env_frm.columnconfigure(4, weight=1)

        for idx, (key, _, desc) in enumerate(self._ENV_DEFS):
            col_off = (idx % 2) * 3
            row = idx // 2
            ttk.Label(env_frm, text=f"{key}:", anchor="w").grid(
                row=row, column=col_off, sticky="w",
                padx=(0 if col_off == 0 else 20, 4), pady=2)
            ent = ttk.Entry(env_frm, textvariable=self.curved_env_vars[key], width=22)
            ent.grid(row=row, column=col_off + 1, sticky="w", pady=2)
            ent.bind("<Enter>", lambda e, d=desc, k=key:
                     self._set_preview_status(f"{k}: {d}"))
            ent.bind("<Leave>", lambda e:
                     self._set_preview_status(
                         "Hover over a beamlet in the preview to inspect coordinates."))

        # ── Tab 3 : Bore Clearance ──────────────────────────────────────────
        bc_tab = ttk.Frame(sub_nb, padding=8)
        bc_tab.columnconfigure(0, weight=1)
        bc_tab.rowconfigure(0, weight=1)
        sub_nb.add(bc_tab, text="Bore Clearance")

        bc_cols = [
            "Grid", "Offset (mm)", "\u03b1 (\u00b0)",
            "Bore Disp (mm)", "Bore r (mm)", "Margin (mm)", "OK",
        ]
        self.bc_tree = ttk.Treeview(bc_tab, columns=bc_cols, show="headings", height=10)
        for c in bc_cols:
            self.bc_tree.heading(c, text=c)
            self.bc_tree.column(c, width=110, anchor="center")
        self.bc_tree.column("Grid", width=75)
        self.bc_tree.column("OK",   width=35)
        self.bc_tree.grid(row=0, column=0, sticky="nsew")
        bc_vsb = ttk.Scrollbar(bc_tab, orient="vertical", command=self.bc_tree.yview)
        self.bc_tree.configure(yscrollcommand=bc_vsb.set)
        bc_vsb.grid(row=0, column=1, sticky="ns")
        self.bc_tree.tag_configure("fail", foreground="#c0392b")
        self.bc_tree.tag_configure("pass", foreground="#2e8b57")

        # ── Tab 4 : Parameter Sweep ─────────────────────────────────────────
        self._build_sweep_tab(sub_nb)

        # ── Action buttons ──────────────────────────────────────────────────
        act_frm = ttk.Frame(ctab, padding=(0, 2, 0, 0))
        act_frm.grid(row=2, column=0, sticky="ew")
        ttk.Button(act_frm, text="Copy Pairs from Case",
                   command=self._copy_pairs_from_case).grid(row=0, column=0, padx=(0, 6))
        ttk.Button(act_frm, text="Open Profile Plot",
                   command=self._open_curved_profile_plot).grid(row=0, column=1, padx=(0, 14))
        ttk.Separator(act_frm, orient="vertical").grid(row=0, column=2, sticky="ns", padx=(0, 14))
        ttk.Button(act_frm, text="Save to File",
                   command=self._save_curved_case_file).grid(row=0, column=3, padx=(0, 6))
        ttk.Button(act_frm, text="Save As\u2026",
                   command=self._save_curved_case_as).grid(row=0, column=4)
        ttk.Label(act_frm, textvariable=self.curved_action_status_var,
                  foreground="#555555", wraplength=640, justify="left").grid(
            row=1, column=0, columnspan=5, sticky="w", pady=(6, 0))

    def _apply_aperture_layout(self) -> None:
        """Compute uniform offset + steering lists from UI fields and store as overrides."""
        if not _HAS_CURVED:
            return
        try:
            n        = int(self.ap_n_var.get())
            offset   = float(self.ap_offset_var.get())  * 1e-3   # mm → m
            pitch    = float(self.ap_pitch_var.get())   * 1e-3   # mm → m
            steering = float(self.ap_steering_var.get()) * 1e-3  # mm → m
        except ValueError as exc:
            messagebox.showerror("Aperture Layout", f"Invalid input: {exc}")
            return
        if pitch <= 0:
            messagebox.showerror("Aperture Layout", "Pitch must be > 0.")
            return
        offsets = uniform_screen_offsets(n, pitch, offset_m=offset)
        self._aperture_offsets_override = offsets
        self._steering_arcs_override = [steering] * n
        offsets_mm  = ", ".join(f"{o * 1e3:.4g}" for o in offsets)
        steering_mm = f"{steering * 1e3:.4g}"
        self.ap_offsets_label_var.set(
            f"Applied: offsets [{offsets_mm}] mm  |  steering {steering_mm} mm  "
            f"(will be used on next Save)"
        )

    def _browse_curved_case(self) -> None:
        # Default the dialog to the two_grid_curved folder if it exists
        _gui_dir = Path(__file__).resolve().parent
        _initial = str(
            next(
                (d for d in [
                    _gui_dir / "two_grid_curved",
                    _gui_dir / "Multi_Grid",
                    _gui_dir,
                ] if d.is_dir()),
                _gui_dir,
            )
        )
        p = filedialog.askopenfilename(
            initialdir=_initial,
            filetypes=[("Python case files", "*case*.py"), ("Python", "*.py"), ("All Files", "*.*")],
            title="Select a curved grid CASE file  (e.g. curved_grid_case_example.py)",
        )
        if p:
            self.curved_case_var.set(p)

    # Pipeline / library filenames that are NOT case files
    _PIPELINE_STEMS = frozenset({
        "curved_grid_pipeline", "multi_grid_pipeline",
        "beamlet_layout_solver", "beamlet_layout_gui",
        "layout_optimizer",
    })

    def _load_curved_case_file(self) -> None:
        if not _HAS_CURVED:
            return
        path_str = self.curved_case_var.get().strip()
        if not path_str:
            messagebox.showinfo("Curved Grid", "Enter or browse for a case file path first.")
            return

        # Guard: reject pipeline/library files early with a helpful message
        stem = Path(path_str).stem
        if stem in self._PIPELINE_STEMS:
            messagebox.showerror(
                "Wrong file selected",
                f"'{Path(path_str).name}' is the pipeline library, not a case file.\n\n"
                "Please select a CASE file instead — for example:\n"
                "  curved_grid_case_example.py\n\n"
                "A case file defines:\n"
                "  CASE = CurvedSimulationCase(grids=[...], env={{...}})\n\n"
                "You can create one with the 'Export Case File from Solution' button below.",
            )
            return

        try:
            case = load_curved_case(Path(path_str))
            self._loaded_curved_case = case

            # Populate grid info text
            self.curved_info_text.configure(state="normal")
            self.curved_info_text.delete("1.0", tk.END)
            for g in case.grids:
                R_str = (
                    f"{g.curvature.radius_m * 1e3:.4g} mm"
                    if not g.curvature.is_flat else "flat"
                )
                orient = "concave\u2191" if g.curvature.concave_upstream else "concave\u2193"
                n_ap = len(g.apertures)
                offsets_mm = ", ".join(
                    f"{ap.offset_m * 1e3:.3g}" for ap in g.apertures
                )
                self.curved_info_text.insert(
                    tk.END,
                    f"{g.name:10s}  V={g.voltage_v:+.0f} V"
                    f"  t={g.thickness_m * 1e3:.4g} mm"
                    f"  R={R_str}  {orient}"
                    f"  gap={g.gap_after_m * 1e3:.3g} mm"
                    f"  mirror={'yes' if g.mirror else 'no'}"
                    f"  {n_ap} aperture(s): [{offsets_mm}] mm\n"
                )
            self.curved_info_text.configure(state="disabled")

            # Bore clearance table
            for item in self.bc_tree.get_children():
                self.bc_tree.delete(item)

            def _inf_str(v: float) -> str:
                return "\u221e" if not math.isfinite(v) else f"{v * 1e3:.3f}"

            report = bore_clearance_report(case)
            n_fail = sum(1 for row in report if not row["ok"])
            for row in report:
                tag = "pass" if row["ok"] else "fail"
                self.bc_tree.insert("", tk.END, tags=(tag,), values=(
                    row["grid"],
                    f"{row['offset_m'] * 1e3:.3f}",
                    f"{row['alpha_deg']:.2f}",
                    _inf_str(row["bore_disp_m"]),
                    f"{row['bore_radius_m'] * 1e3:.3f}",
                    _inf_str(row["margin_m"]),
                    "\u2713" if row["ok"] else "\u2717",
                ))

            n_grids = len(case.grids)
            n_ap    = len(report)
            if n_fail == 0:
                status = (
                    f"Loaded {n_grids} grid(s), {n_ap} aperture(s).  "
                    f"All bore clearances OK."
                )
            else:
                status = (
                    f"Loaded {n_grids} grid(s), {n_ap} aperture(s).  "
                    f"\u26a0 {n_fail} aperture(s) fail bore clearance!"
                )
            self.curved_status_var.set(status)
            self.curved_action_status_var.set("")

            # ── Populate editable fields (only if tab has been built) ───────
            if hasattr(self, "curved_geom_vars"):
                g0 = case.grids[0]
                g1 = case.grids[1] if len(case.grids) > 1 else None
                gv = self.curved_geom_vars
                gv["R_SCR"].set(
                    f"{g0.curvature.radius_m * 1e3:.6g}"
                    if not g0.curvature.is_flat else "inf"
                )
                gv["T_SCR"].set(f"{g0.thickness_m * 1e3:.6g}")
                gv["T_ACC"].set(f"{g1.thickness_m * 1e3:.6g}" if g1 else "5")
                gv["GAP"].set(f"{g0.gap_after_m * 1e3:.6g}")
                gv["VS"].set(f"{g0.voltage_v:.1f}")
                gv["VA"].set(f"{g1.voltage_v:.1f}" if g1 else "-15000.0")
                bore_radii = [ap.radius_m for ap in g0.apertures if ap.radius_m > 0]
                gv["AP_RAD"].set(
                    f"{bore_radii[0] * 1e3:.6g}" if bore_radii else "1.5"
                )

            if hasattr(self, "curved_env_vars"):
                for key, default, _ in self._ENV_DEFS:
                    self.curved_env_vars[key].set(case.env.get(key, default))

            # Back-populate aperture layout fields from the loaded case
            if hasattr(self, "ap_n_var"):
                existing = sorted({abs(ap.offset_m) for ap in case.grids[0].apertures})
                if existing:
                    n = len(existing)
                    offset_m = existing[0]
                    pitch_m  = (existing[1] - existing[0]) if n >= 2 else 0.0
                    self.ap_n_var.set(str(n))
                    self.ap_offset_var.set(f"{offset_m * 1e3:.4g}")
                    self.ap_pitch_var.set(f"{pitch_m * 1e3:.4g}" if pitch_m > 0 else "10.0")
                    self._aperture_offsets_override = None
                    self._steering_arcs_override = None
                    # Back-populate steering from first pair (representative)
                    try:
                        pairs = pairs_from_case(case)
                        rep_steering = pairs[0].steering_arc_m if pairs else 0.0
                    except Exception:
                        rep_steering = 0.0
                    self.ap_steering_var.set(f"{rep_steering * 1e3:.4g}")
                    offsets_mm = ", ".join(f"{o * 1e3:.3g}" for o in existing)
                    self.ap_offsets_label_var.set(
                        f"From case: [{offsets_mm}] mm  |  steering {rep_steering * 1e3:.4g} mm"
                        f"  (Apply Layout to override)"
                    )

        except TypeError as exc:
            msg = str(exc)
            if "must define CASE" in msg or "CASE" in msg:
                messagebox.showerror(
                    "Curved Case Error",
                    f"The selected file does not define a valid case.\n\n"
                    f"A case file must contain:\n"
                    f"  CASE = CurvedSimulationCase(grids=[...], env={{...}})\n\n"
                    f"You selected:\n  {path_str}\n\n"
                    f"Make sure you are opening a case file such as\n"
                    f"  curved_grid_case_example.py\n"
                    f"not the pipeline library itself.\n\n"
                    f"Detail: {exc}",
                )
            else:
                messagebox.showerror(
                    "Curved Case Error", f"{exc}\n\n{traceback.format_exc()}")
        except Exception as exc:
            messagebox.showerror(
                "Curved Case Error", f"{exc}\n\n{traceback.format_exc()}")

    def _save_curved_case_file(self, path=None) -> None:
        """Write current geometry + env fields back to the case file."""
        if not _HAS_CURVED:
            return
        lc = self._loaded_curved_case
        if lc is None:
            messagebox.showinfo("Curved Grid", "Load a case file first.")
            return
        if path is None:
            path = Path(self.curved_case_var.get().strip())
        try:
            gv = self.curved_geom_vars
            R_scr_m = float(gv["R_SCR"].get()) * 1e-3
            t_scr_m = float(gv["T_SCR"].get()) * 1e-3
            t_acc_m = float(gv["T_ACC"].get()) * 1e-3
            gap_m   = float(gv["GAP"].get())   * 1e-3
            vs      = float(gv["VS"].get())
            va      = float(gv["VA"].get())
            ap_m    = float(gv["AP_RAD"].get()) * 1e-3
            env_dict = {k: v.get().strip()
                        for k, v in self.curved_env_vars.items()
                        if v.get().strip()}
            if self._aperture_offsets_override is not None:
                offsets = self._aperture_offsets_override
            else:
                offsets = sorted({abs(ap.offset_m) for ap in lc.grids[0].apertures})
            # Compute accel offsets from steering arcs if an override is active
            if self._steering_arcs_override is not None:
                R_acc_m = concentric_accel_radius(R_scr_m, t_scr_m, gap_m, offsets)
                accel_offsets_m = [
                    accel_offset_from_steering(r, R_scr_m, R_acc_m, s)
                    for r, s in zip(offsets, self._steering_arcs_override)
                ]
            else:
                R_acc_m = None
                accel_offsets_m = None
            write_case_file(
                path,
                screen_offsets_m=offsets,
                accel_offsets_m=accel_offsets_m,
                ap_radius_m=ap_m,
                R_scr_m=R_scr_m,
                R_acc_m=R_acc_m,
                screen_voltage_v=vs,
                accel_voltage_v=va,
                screen_thickness_m=t_scr_m,
                accel_thickness_m=t_acc_m,
                gap_after_m=gap_m,
                env=env_dict,
            )
            self.curved_case_var.set(str(path))
            self.curved_action_status_var.set(f"Saved: {Path(path).name}")
            self._load_curved_case_file()
        except Exception as exc:
            messagebox.showerror("Save Error", f"{exc}\n\n{traceback.format_exc()}")

    def _save_curved_case_as(self) -> None:
        """Pick a new path then save all current field values to it."""
        if not _HAS_CURVED:
            return
        if self._loaded_curved_case is None:
            messagebox.showinfo("Curved Grid", "Load a case file first.")
            return
        p = filedialog.asksaveasfilename(
            defaultextension=".py",
            filetypes=[("Python", "*.py"), ("All Files", "*.*")],
            title="Save curved case file as\u2026",
        )
        if p:
            self._save_curved_case_file(Path(p))

    def _copy_pairs_from_case(self) -> None:
        if not _HAS_CURVED:
            return
        case = self._loaded_curved_case
        if case is None:
            messagebox.showinfo("Curved Grid", "Load a case file first.")
            return
        try:
            pairs = pairs_from_case(case)
            lines = [
                f"    BeamletPair("
                f"screen_offset_m={p.screen_offset_m:.8g}, "
                f"steering_arc_m={p.steering_arc_m:.8g}, "
                f"screen_radius_m={p.screen_radius_m:.8g})"
                for p in pairs
            ]
            text = "PAIRS = [\n" + ",\n".join(lines) + "\n]"
            self.clipboard_clear()
            self.clipboard_append(text)
            self.update()
            self.curved_action_status_var.set(
                f"Copied {len(pairs)} BeamletPair(s) to clipboard."
            )
        except Exception as exc:
            messagebox.showerror("Copy Error", f"{exc}\n\n{traceback.format_exc()}")

    def _open_curved_profile_plot(self) -> None:
        if not _HAS_CURVED:
            return
        case = self._loaded_curved_case
        if case is None:
            messagebox.showinfo("Curved Grid", "Load a case file first.")
            return
        try:
            import tempfile
            fig = plot_curved_grid_profile(
                case.grids,
                title=Path(self.curved_case_var.get()).stem,
            )
            tmp = tempfile.NamedTemporaryFile(
                suffix=".png", delete=False, prefix="curved_profile_")
            fig.savefig(tmp.name, dpi=130, bbox_inches="tight")
            tmp.close()
            webbrowser.open(Path(tmp.name).resolve().as_uri())
            self.curved_action_status_var.set(f"Profile plot saved: {tmp.name}")
        except ImportError:
            messagebox.showinfo(
                "No Matplotlib",
                "matplotlib and numpy are required for profile plots.\n"
                "Install with:  pip install matplotlib numpy",
            )
        except Exception as exc:
            messagebox.showerror("Plot Error", f"{exc}\n\n{traceback.format_exc()}")

    def _export_curved_case_dialog(self) -> None:
        """Kept as a thin alias for Save As (dialog replaced by inline fields)."""
        self._save_curved_case_as()

    # ------------------------------------------------------------------
    # Parameter Sweep tab
    # ------------------------------------------------------------------

    def _build_sweep_tab(self, sub_nb: ttk.Notebook) -> None:
        sw_tab = ttk.Frame(sub_nb, padding=8)
        sw_tab.columnconfigure(0, weight=1)
        sub_nb.add(sw_tab, text="Parameter Sweep")

        if not _HAS_SWEEP:
            ttk.Label(
                sw_tab,
                text="sweep_types.py not found in sim/.\n"
                     "Place sim/sweep_types.py in the repository and restart.",
                foreground="red",
            ).grid(row=0, column=0, sticky="w")
            return

        # ── Axis list ───────────────────────────────────────────────────
        axes_frm = ttk.LabelFrame(sw_tab, text="Sweep Axes", padding=6)
        axes_frm.grid(row=0, column=0, sticky="nsew", pady=(0, 6))
        axes_frm.columnconfigure(0, weight=1)
        axes_frm.rowconfigure(0, weight=1)
        sw_tab.rowconfigure(0, weight=1)

        ax_cols = ["Parameter", "Values", "Count"]
        self._sweep_tree = ttk.Treeview(
            axes_frm, columns=ax_cols, show="headings", height=6, selectmode="browse")
        self._sweep_tree.heading("Parameter", text="Parameter")
        self._sweep_tree.heading("Values",    text="Values")
        self._sweep_tree.heading("Count",     text="N")
        self._sweep_tree.column("Parameter", width=160, anchor="w")
        self._sweep_tree.column("Values",    width=380, anchor="w")
        self._sweep_tree.column("Count",     width=42,  anchor="center")
        self._sweep_tree.grid(row=0, column=0, sticky="nsew")
        ax_vsb = ttk.Scrollbar(axes_frm, orient="vertical", command=self._sweep_tree.yview)
        self._sweep_tree.configure(yscrollcommand=ax_vsb.set)
        ax_vsb.grid(row=0, column=1, sticky="ns")

        ax_btns = ttk.Frame(axes_frm)
        ax_btns.grid(row=1, column=0, sticky="w", pady=(4, 0))
        ttk.Button(ax_btns, text="Add Axis…",    command=self._sweep_add_axis).grid(row=0, column=0, padx=(0, 4))
        ttk.Button(ax_btns, text="Edit…",        command=self._sweep_edit_axis).grid(row=0, column=1, padx=(0, 4))
        ttk.Button(ax_btns, text="Remove Axis",  command=self._sweep_remove_axis).grid(row=0, column=2)

        # ── Mode + options ──────────────────────────────────────────────
        opt_frm = ttk.LabelFrame(sw_tab, text="Sweep Options", padding=6)
        opt_frm.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        opt_frm.columnconfigure(3, weight=1)

        ttk.Label(opt_frm, text="Mode:").grid(row=0, column=0, sticky="w", padx=(0, 4))
        ttk.Radiobutton(opt_frm, text="Product  (Cartesian)",
                        variable=self._sweep_mode_var, value="product",
                        command=self._sweep_update_total).grid(row=0, column=1, sticky="w")
        ttk.Radiobutton(opt_frm, text="Zip  (element-wise)",
                        variable=self._sweep_mode_var, value="zip",
                        command=self._sweep_update_total).grid(row=0, column=2, sticky="w", padx=(12, 0))

        ttk.Label(opt_frm, text="Tag:").grid(row=1, column=0, sticky="w", pady=(6, 0), padx=(0, 4))
        ttk.Entry(opt_frm, textvariable=self._sweep_tag_var, width=22).grid(
            row=1, column=1, columnspan=2, sticky="w", pady=(6, 0))

        ttk.Checkbutton(opt_frm, text="Run reduce_best after all tasks finish",
                        variable=self._sweep_reduce_var).grid(
            row=2, column=0, columnspan=4, sticky="w", pady=(6, 0))

        self._sweep_total_var = tk.StringVar(value="Total combinations: —")
        ttk.Label(opt_frm, textvariable=self._sweep_total_var,
                  foreground="#0055cc").grid(row=3, column=0, columnspan=4, sticky="w", pady=(4, 0))

        # ── Generate ────────────────────────────────────────────────────
        gen_frm = ttk.Frame(sw_tab, padding=(0, 2, 0, 0))
        gen_frm.grid(row=2, column=0, sticky="ew")
        gen_frm.columnconfigure(2, weight=1)

        ttk.Button(gen_frm, text="Generate Sweep Case File…",
                   command=self._sweep_generate).grid(row=0, column=0, padx=(0, 8))
        ttk.Button(gen_frm, text="Clear All Axes",
                   command=self._sweep_clear).grid(row=0, column=1)

        self._sweep_status_var = tk.StringVar(value="")
        ttk.Label(gen_frm, textvariable=self._sweep_status_var,
                  foreground="#555555", wraplength=680, justify="left").grid(
            row=1, column=0, columnspan=3, sticky="w", pady=(6, 0))

    # ── All parameter names the sweep axis dialog offers ────────────────
    @property
    def _all_sweep_param_names(self) -> list:
        # All geometry params are valid sweep targets
        geom = [k for k, *_ in self._GEOM_DEFS]
        # Physics / domain env vars only — exclude run-config knobs
        env  = [k for k, *_ in self._ENV_DEFS if k not in self._SWEEP_ENV_EXCLUDE]
        return geom + env + ["custom…"]

    def _sweep_axis_dialog(self, existing: "SweepAxis | None" = None) -> "SweepAxis | None":
        """Modal dialog for adding or editing one SweepAxis.
        Returns the new SweepAxis, or None if cancelled.
        """
        dlg = tk.Toplevel(self)
        dlg.title("Add Sweep Axis" if existing is None else "Edit Sweep Axis")
        dlg.resizable(False, False)
        dlg.grab_set()

        param_var  = tk.StringVar(value=existing.name if existing else "AP_RAD")
        custom_var = tk.StringVar(value="" if existing is None or existing.name in self._all_sweep_param_names else existing.name)
        values_var = tk.StringVar()
        result: list = [None]   # mutable container so inner func can write it

        # Pre-fill values string from existing axis
        if existing is not None:
            vals = existing.values
            if len(vals) <= 20:
                values_var.set(", ".join(str(v) for v in vals))
            else:
                values_var.set(", ".join(str(v) for v in vals))

        frm = ttk.Frame(dlg, padding=14)
        frm.grid()
        frm.columnconfigure(1, weight=1)

        row = 0
        ttk.Label(frm, text="Parameter:").grid(row=row, column=0, sticky="w", padx=(0, 8))
        param_cb = ttk.Combobox(frm, textvariable=param_var,
                                values=self._all_sweep_param_names,
                                state="readonly", width=28)
        param_cb.grid(row=row, column=1, sticky="ew")
        row += 1

        custom_lbl = ttk.Label(frm, text="Custom name:")
        custom_ent = ttk.Entry(frm, textvariable=custom_var, width=28)

        def _on_param_change(*_):
            if param_var.get() == "custom…":
                custom_lbl.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=(4, 0))
                custom_ent.grid(row=1, column=1, sticky="ew", pady=(4, 0))
            else:
                custom_lbl.grid_remove()
                custom_ent.grid_remove()
        param_cb.bind("<<ComboboxSelected>>", _on_param_change)
        _on_param_change()

        row = 2
        ttk.Label(frm, text="Values:").grid(row=row, column=0, sticky="nw", padx=(0, 8), pady=(8, 0))
        ttk.Entry(frm, textvariable=values_var, width=48).grid(
            row=row, column=1, sticky="ew", pady=(8, 0))
        row += 1

        ttk.Label(
            frm,
            text=(
                "Examples:\n"
                "  Comma list:        0.001, 0.0015, 0.002\n"
                "  linspace(a, b, n): linspace(1e16, 5e16, 5)\n"
                "  arange(a, b, s):   arange(0.001, 0.003, 0.0005)"
            ),
            foreground="#666666", justify="left",
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(2, 8))
        row += 1

        def _parse_values(s: str) -> list:
            s = s.strip()
            if s.startswith("linspace("):
                inner = s[len("linspace("):-1]
                a, b, n = [x.strip() for x in inner.split(",")]
                return linspace(float(a), float(b), int(n))
            if s.startswith("arange("):
                inner = s[len("arange("):-1]
                a, b, step = [x.strip() for x in inner.split(",")]
                return arange(float(a), float(b), float(step))
            # plain comma-separated
            return [float(x.strip()) for x in s.split(",") if x.strip()]

        def do_ok():
            pname = param_var.get()
            if pname == "custom…":
                pname = custom_var.get().strip()
                if not pname:
                    messagebox.showerror("Sweep Axis", "Enter a custom parameter name.", parent=dlg)
                    return
            try:
                vals = _parse_values(values_var.get())
                if not vals:
                    raise ValueError("No values parsed — check your input.")
            except Exception as exc:
                messagebox.showerror("Sweep Axis", f"Could not parse values:\n{exc}", parent=dlg)
                return
            result[0] = SweepAxis(name=pname, values=vals)
            dlg.destroy()

        ttk.Button(frm, text="OK",     command=do_ok).grid(row=row, column=0, pady=(4, 0))
        ttk.Button(frm, text="Cancel", command=dlg.destroy).grid(row=row, column=1, sticky="w", pady=(4, 0))

        dlg.wait_window()
        return result[0]

    def _sweep_add_axis(self) -> None:
        ax = self._sweep_axis_dialog()
        if ax is None:
            return
        self._sweep_axes.append(ax)
        self._sweep_tree_refresh()
        self._sweep_update_total()

    def _sweep_edit_axis(self) -> None:
        sel = self._sweep_tree.selection()
        if not sel:
            messagebox.showinfo("Parameter Sweep", "Select an axis to edit.")
            return
        idx = self._sweep_tree.index(sel[0])
        existing = self._sweep_axes[idx]
        ax = self._sweep_axis_dialog(existing=existing)
        if ax is None:
            return
        self._sweep_axes[idx] = ax
        self._sweep_tree_refresh()
        self._sweep_update_total()

    def _sweep_remove_axis(self) -> None:
        sel = self._sweep_tree.selection()
        if not sel:
            messagebox.showinfo("Parameter Sweep", "Select an axis to remove.")
            return
        idx = self._sweep_tree.index(sel[0])
        del self._sweep_axes[idx]
        self._sweep_tree_refresh()
        self._sweep_update_total()

    def _sweep_clear(self) -> None:
        self._sweep_axes.clear()
        self._sweep_tree_refresh()
        self._sweep_update_total()

    def _sweep_tree_refresh(self) -> None:
        for item in self._sweep_tree.get_children():
            self._sweep_tree.delete(item)
        for ax in self._sweep_axes:
            vals = ax.values
            if len(vals) <= 6:
                val_str = ", ".join(str(v) for v in vals)
            else:
                val_str = (
                    f"{vals[0]}, {vals[1]}, … {vals[-2]}, {vals[-1]}"
                    f"  ({len(vals)} values)"
                )
            self._sweep_tree.insert("", tk.END, values=(ax.name, val_str, len(vals)))

    def _sweep_update_total(self) -> None:
        if not self._sweep_axes:
            self._sweep_total_var.set("Total combinations: —")
            return
        try:
            spec = SweepSpec(
                axes=self._sweep_axes,
                mode=self._sweep_mode_var.get(),
                reduce=self._sweep_reduce_var.get(),
                tag=self._sweep_tag_var.get(),
            )
            n = spec.total
            self._sweep_total_var.set(f"Total combinations: {n}  →  SLURM --array=0-{n - 1}")
        except Exception as exc:
            self._sweep_total_var.set(f"⚠  {exc}")

    def _sweep_generate(self) -> None:
        """Write a sweep-enabled case file and print the sbatch command."""
        if not _HAS_SWEEP:
            messagebox.showinfo("Parameter Sweep", "sweep_types not available.")
            return
        if not self._sweep_axes:
            messagebox.showinfo("Parameter Sweep", "Add at least one sweep axis first.")
            return

        # Need a case file to wrap
        case_path_str = self.curved_case_var.get().strip()
        if not case_path_str:
            messagebox.showinfo("Parameter Sweep",
                                "Load a Curved Grid case file first (use Browse / Load Case).")
            return

        # Ask where to save the sweep case file
        out_path_str = filedialog.asksaveasfilename(
            defaultextension=".py",
            initialfile=Path(case_path_str).stem + "_sweep.py",
            initialdir=str(Path(case_path_str).parent),
            filetypes=[("Python", "*.py"), ("All Files", "*.*")],
            title="Save sweep case file as…",
        )
        if not out_path_str:
            return

        try:
            spec = SweepSpec(
                axes=self._sweep_axes,
                mode=self._sweep_mode_var.get(),
                reduce=self._sweep_reduce_var.get(),
                tag=self._sweep_tag_var.get(),
            )
        except Exception as exc:
            messagebox.showerror("Sweep Error", str(exc))
            return

        # Import write_sweep_case_file from the pipeline
        try:
            from curved_grid_pipeline import write_sweep_case_file  # type: ignore
        except ImportError:
            messagebox.showerror(
                "Sweep Error",
                "curved_grid_pipeline does not export write_sweep_case_file.\n"
                "Ensure the pipeline has been updated with sweep support.",
            )
            return

        # Gather current geometry / env overrides from the GUI fields so the
        # generated sweep file has the most up-to-date base values.
        gui_overrides: dict = {}
        if hasattr(self, "curved_geom_vars"):
            for k, var in self.curved_geom_vars.items():
                v = var.get().strip()
                if v:
                    gui_overrides[k] = v
        if hasattr(self, "curved_env_vars"):
            for k, var in self.curved_env_vars.items():
                v = var.get().strip()
                if v:
                    gui_overrides[k] = v

        try:
            out_path = Path(out_path_str)
            write_sweep_case_file(
                out_path,
                sweep_axes=spec.axes,
                sweep_mode=spec.mode,
                sweep_reduce=spec.reduce,
                sweep_tag=spec.tag,
                source_case_path=Path(case_path_str),
                gui_overrides=gui_overrides,
            )
        except Exception as exc:
            messagebox.showerror("Sweep Error", f"{exc}\n\n{traceback.format_exc()}")
            return

        # Locate the sweep SLURM script
        _gui_dir = Path(__file__).resolve().parent
        slurm_candidates = [
            _gui_dir.parent / "curved_grid" / "orchestrate_curved_grid_sweep.slurm",
            _gui_dir / "orchestrate_curved_grid_sweep.slurm",
        ]
        slurm_path = next((p for p in slurm_candidates if p.exists()), None)
        slurm_name = slurm_path.name if slurm_path else "orchestrate_curved_grid_sweep.slurm"

        sbatch_cmd = (
            f"sbatch \\\n"
            f"  --array=0-{spec.total - 1} \\\n"
            f"  --export=SWEEP_CONFIG={out_path.name},SWEEP_BINARY=../sim/multi_grid_2d_curved \\\n"
            f"  {slurm_name}"
        )

        msg = (
            f"Sweep case file written:\n  {out_path}\n\n"
            f"Run on HPC (from curved_grid/):\n\n{sbatch_cmd}\n\n"
            f"Total tasks: {spec.total}  |  Mode: {spec.mode}  |  Reduce: {spec.reduce}"
        )

        dlg = tk.Toplevel(self)
        dlg.title("Sweep Generated")
        dlg.resizable(True, True)
        dlg.grab_set()
        frm = ttk.Frame(dlg, padding=14)
        frm.grid(sticky="nsew")
        dlg.columnconfigure(0, weight=1)
        dlg.rowconfigure(0, weight=1)
        frm.columnconfigure(0, weight=1)
        frm.rowconfigure(0, weight=1)
        txt = tk.Text(frm, wrap="word", width=72, height=16,
                      font=("Consolas", 9), background="#f8f8f8")
        txt.insert("1.0", msg)
        txt.configure(state="disabled")
        txt.grid(row=0, column=0, sticky="nsew")
        vsb2 = ttk.Scrollbar(frm, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=vsb2.set)
        vsb2.grid(row=0, column=1, sticky="ns")

        def _copy():
            self.clipboard_clear()
            self.clipboard_append(sbatch_cmd)
            self.update()

        btn_frm = ttk.Frame(frm)
        btn_frm.grid(row=1, column=0, columnspan=2, pady=(8, 0), sticky="w")
        ttk.Button(btn_frm, text="Copy sbatch command", command=_copy).grid(row=0, column=0, padx=(0, 8))
        ttk.Button(btn_frm, text="Close", command=dlg.destroy).grid(row=0, column=1)

        self._sweep_status_var.set(
            f"Wrote {out_path.name}  ({spec.total} tasks, mode={spec.mode})"
        )


def main() -> int:
    app = App()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
