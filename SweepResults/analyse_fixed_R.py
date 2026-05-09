"""
analyse_fixed_R.py
------------------
Preliminary analysis of the fixed-R (R_s=0.3 m) curved-grid sweep.

Sweep parameters:
    ni  : 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0  x10^18 m^-3   (7 values)
    ra  : 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0  mm             (7 values)
    R   : fixed at 0.300 m (screen), 0.292 m (accel)

N_APERTURES_FULL : 469  (preliminary full-grid count, subject to recount by ra)

Outputs (saved to FIGS_DIR):
    fig1_Idel_heatmap.png         -- 2D heatmap of I_del in (ni, ra) space
    fig2_current_vs_aprad.png     -- per-beamlet I vs ra, lines by ni   [paper fig]
    fig3_current_vs_ni.png        -- per-beamlet I vs ni, lines by ra   [paper fig]
    fig4_loss_budget.png          -- grid-loss % + leakage % vs ra      [paper fig]
    fig5_uniformity.png           -- peak_to_avg and edge_peaking vs ra
    fig6_best_profile.png         -- radial profile at design point      [paper fig]
    fig7_sidewall_map.png         -- diagnostic: which configs hit tube
    sweep_summary.csv             -- full scalar table
"""

import json
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm

# ── Paths ────────────────────────────────────────────────────────────────────
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
FIGS_DIR    = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(FIGS_DIR, exist_ok=True)

N_SIM_AP   = 25       # apertures in each simulation
N_FULL_AP  = 469      # preliminary full-grid aperture count
SCALE      = N_FULL_AP / N_SIM_AP
R_TARGET_A = 4.0      # delivery requirement
LOSS_LIMIT = 0.10     # 10 % grid-loss ceiling
R_SAMPLE_M = 0.03175  # 31.75 mm sample radius

# ── Plot style ───────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":     "serif",
    "font.size":       10,
    "axes.labelsize":  11,
    "axes.titlesize":  11,
    "legend.fontsize": 9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top":       True,
    "ytick.right":     True,
    "figure.dpi":      150,
})

# ── Ingest ───────────────────────────────────────────────────────────────────
records = {}   # keyed by (ni, ra) -- first successful copy wins

for run_dir in sorted(glob.glob(os.path.join(RESULTS_DIR, "**/"), recursive=False)):
    bm_path   = os.path.join(run_dir, "beam_metrics.json")
    meta_path = os.path.join(run_dir, "meta.json")
    prof_path = os.path.join(run_dir, "sample_diameter_profile.json")

    if not all(os.path.exists(p) for p in [bm_path, meta_path, prof_path]):
        continue

    with open(meta_path)  as f: meta = json.load(f)
    with open(bm_path)    as f: bm   = json.load(f)
    with open(prof_path)  as f: prof = json.load(f)

    ni = meta["physics"]["PLASMA_NI_M3"]
    ra = meta["AP_RAD_M"]
    key = (ni, ra)
    if key in records:
        continue   # keep first successful copy

    col   = bm["collimation"]
    cur   = bm["currents"]
    samp  = bm["sample"]
    stats = prof["stats"]

    # ── Per-beamlet currents ──────────────────────────────────────────────────
    # Prefer the new per_beamlet windowed measurements (written by the
    # updated multi_grid_2d_curved.cpp).  These measure current at the
    # screen-exit plane and accel-exit plane through each aperture's own bore,
    # avoiding the cross-aperture contamination in the global planes.
    #
    # Fall back to the accel-exit global plane for I_ag (which is consistent
    # within a run) and NaN for grid_loss when only old data is available.
    pb_list = bm.get("beamlet_currents", {}).get("per_beamlet", [])

    if pb_list:
        # New per-beamlet measurement available
        I_pg_total  = sum(b["I_pg_A"] for b in pb_list)
        I_ag_total  = sum(b["I_ag_A"] for b in pb_list)
        grid_loss_frac = (I_pg_total - I_ag_total) / I_pg_total if I_pg_total > 0 else np.nan

        # Per-beamlet loss as function of radial position (for radial-loss plot)
        pb_radial = [
            {"r_m": abs(b["scr_offset_m"]), "loss_frac": b["loss_frac"]}
            for b in pb_list
        ]
        # Average current exiting accel per simulated beamlet (for delivered current scaling)
        I_beamlet_ag = I_ag_total / N_SIM_AP if N_SIM_AP > 0 else 0.0
    else:
        # Legacy data: use global accel-exit plane (I_ag_out_A is self-consistent
        # within a run but I_pg_in_A is unreliable, so grid_loss is NaN)
        I_beamlet_ag   = cur["I_ag_out_A"] / N_SIM_AP
        grid_loss_frac = np.nan   # broken for curved-grid global planes
        pb_radial      = []

    # per-beamlet current reaching target plane
    I_beamlet_sm  = samp["I_A"]           / N_SIM_AP
    # per-beamlet in-sample current (within R_sample)
    I_beamlet_in  = prof["total_I_A_in"]  / N_SIM_AP

    # total delivered to full grid, in-sample
    I_del = I_beamlet_in * N_FULL_AP

    records[key] = {
        "ni_m3":          ni,
        "ni_e18":         ni / 1e18,
        "ra_mm":          ra * 1e3,
        "I_beamlet_ag_A": I_beamlet_ag,
        "I_beamlet_sm_A": I_beamlet_sm,
        "I_beamlet_in_A": I_beamlet_in,
        "I_del_A":        I_del,
        "grid_loss_frac": grid_loss_frac,
        "has_per_beamlet": bool(pb_list),
        "pb_radial":      pb_radial,          # list[{r_m, loss_frac}] or []
        "leakage_frac":   stats["leakage_frac"],
        "peak_to_avg":    stats["peak_to_avg"],
        "rms_nonuniform": stats["rms_nonuniform"],
        "edge_peaking":   stats["edge_peaking_index"],
        "lost_to_tube":   col["lost_to_sidewalls"],
        "has_sample":     col["has_sample_beam"],
        "y_rms_sm_mm":    samp["y_rms_m"] * 1e3,
        "y_absmax_sm_mm": samp["y_absmax_m"] * 1e3,
        "profile_bins":   prof["bins"],
    }

SCALAR_EXCLUDE = {"profile_bins", "pb_radial"}
df = pd.DataFrame([{k: v for k, v in r.items() if k not in SCALAR_EXCLUDE}
                   for r in records.values()])
df = df.sort_values(["ni_e18", "ra_mm"]).reset_index(drop=True)

# Mask tube-hit runs -- replace I_del with NaN for valid heatmap
df["I_del_clean"] = df.apply(
    lambda r: np.nan if r["lost_to_tube"] else r["I_del_A"], axis=1)

n_with_pb = df["has_per_beamlet"].sum()
print(f"Loaded {len(df)} unique (ni, ra) combinations")
print(f"  {n_with_pb}/{len(df)} have per-beamlet windowed current data (new C++ diagnostic)")
if n_with_pb < len(df):
    print(f"  NOTE: {len(df)-n_with_pb} runs use legacy data — grid_loss_frac set to NaN")
print(f"Lost to sidewalls: {df['lost_to_tube'].sum()}")
print(f"Configurations with I_del >= {R_TARGET_A} A: "
      f"{(df['I_del_clean'] >= R_TARGET_A).sum()}")
has_loss_data = df["grid_loss_frac"].notna()
print(f"... also satisfying grid_loss <= {LOSS_LIMIT*100:.0f}% (of {has_loss_data.sum()} with valid loss): "
      f"{((df['I_del_clean'] >= R_TARGET_A) & (df['grid_loss_frac'] <= LOSS_LIMIT)).sum()}")

# ── Pivot tables ─────────────────────────────────────────────────────────────
nis = sorted(df["ni_e18"].unique())
ras = sorted(df["ra_mm"].unique())

def pivot(col):
    return df.pivot_table(index="ni_e18", columns="ra_mm", values=col)

piv_Idel   = pivot("I_del_clean")
piv_Iblt   = pivot("I_beamlet_ag_A")
piv_gloss  = pivot("grid_loss_frac")
piv_leak   = pivot("leakage_frac")
piv_pta    = pivot("peak_to_avg")
piv_edge   = pivot("edge_peaking")
piv_tube   = pivot("lost_to_tube")

# ── Colour helpers ───────────────────────────────────────────────────────────
NI_CMAP  = plt.cm.plasma
RA_CMAP  = plt.cm.viridis
ni_colors = {n: NI_CMAP(i / (len(nis) - 1)) for i, n in enumerate(nis)}
ra_colors = {r: RA_CMAP(i / (len(ras) - 1)) for i, r in enumerate(ras)}

# ═══════════════════════════════════════════════════════════════════════════
# Fig 1 -- I_del heatmap in (ni, ra) space
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(6.5, 4.5))

Zidel = piv_Idel.values.astype(float)
ni_arr = piv_Idel.index.values
ra_arr = piv_Idel.columns.values

im = ax.pcolormesh(ra_arr, ni_arr, Zidel,
                   cmap="viridis", shading="nearest",
                   vmin=0, vmax=min(Zidel[np.isfinite(Zidel)].max(), 15))
cb = fig.colorbar(im, ax=ax, label=r"$I_\mathrm{del}$ (A, 469 apertures)")

# 4 A contour
try:
    cs4 = ax.contour(ra_arr, ni_arr, Zidel, levels=[R_TARGET_A],
                     colors="white", linewidths=1.5)
    ax.clabel(cs4, fmt="4 A", fontsize=8, inline=True)
except Exception:
    pass

# 10 % grid-loss contour overlay
Zgloss = piv_gloss.values.astype(float)
try:
    csg = ax.contour(ra_arr, ni_arr, Zgloss, levels=[LOSS_LIMIT],
                     colors="red", linewidths=1.2, linestyles="--")
    ax.clabel(csg, fmt="10 % loss", fontsize=8, inline=True)
except Exception:
    pass

# Tube-hit mask (hatch)
tube_mask = piv_tube.values.astype(bool)
for i, ni_v in enumerate(ni_arr):
    for j, ra_v in enumerate(ra_arr):
        if tube_mask[i, j]:
            ax.add_patch(plt.Rectangle(
                (ra_v - 0.25, ni_v - 0.25), 0.5, 0.5,
                fill=False, hatch="////", edgecolor="gray", linewidth=0))

ax.set_xlabel(r"Aperture radius $r_a$ (mm)")
ax.set_ylabel(r"Ion density $n_i$ ($\times 10^{18}$ m$^{-3}$)")
ax.set_title(r"Delivered current $I_\mathrm{del}$ — fixed $R_s = 0.30$ m")
ax.set_xticks(ras)
ax.set_yticks(nis)
fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig1_Idel_heatmap.png"), dpi=200)
plt.close(fig)
print("Saved fig1_Idel_heatmap.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 2 -- Per-beamlet I vs ra, lines by ni  [paper: fig:current_vs_aprad]
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(5.5, 4.0))

for ni_v in nis:
    sub = df[df["ni_e18"] == ni_v].sort_values("ra_mm")
    ax.plot(sub["ra_mm"], sub["I_beamlet_ag_A"] * 1e3,
            marker="o", markersize=4, color=ni_colors[ni_v],
            label=rf"$n_i = {ni_v:.1f}\times10^{{18}}$ m$^{{-3}}$")

ax.set_xlabel(r"Aperture radius $r_a$ (mm)")
ax.set_ylabel(r"Per-beamlet extracted current (mA)")
ax.set_title(r"Per-beamlet current vs aperture radius")
ax.legend(fontsize=7.5, ncol=2, loc="upper left")
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig2_current_vs_aprad.png"), dpi=200)
plt.close(fig)
print("Saved fig2_current_vs_aprad.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 3 -- Per-beamlet I vs ni, lines by ra  [replaces current_vs_curvrad]
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(5.5, 4.0))

for ra_v in ras:
    sub = df[df["ra_mm"] == ra_v].sort_values("ni_e18")
    ax.plot(sub["ni_e18"], sub["I_beamlet_ag_A"] * 1e3,
            marker="s", markersize=4, color=ra_colors[ra_v],
            label=rf"$r_a = {ra_v:.1f}$ mm")

ax.set_xlabel(r"Ion density $n_i$ ($\times 10^{18}$ m$^{-3}$)")
ax.set_ylabel(r"Per-beamlet extracted current (mA)")
ax.set_title(r"Per-beamlet current vs ion density")
ax.legend(fontsize=7.5, ncol=2, loc="upper left")
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig3_current_vs_ni.png"), dpi=200)
plt.close(fig)
print("Saved fig3_current_vs_ni.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 4 -- Loss budget vs ra, lines by ni  [paper: fig:loss_budget]
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), sharey=False)

for ni_v in nis:
    sub = df[df["ni_e18"] == ni_v].sort_values("ra_mm")
    c   = ni_colors[ni_v]
    lbl = rf"$n_i={ni_v:.1f}$"
    axes[0].plot(sub["ra_mm"], sub["grid_loss_frac"] * 100,
                 marker="o", markersize=4, color=c, label=lbl)
    axes[1].plot(sub["ra_mm"], sub["leakage_frac"] * 100,
                 marker="o", markersize=4, color=c, label=lbl)

for ax, title, ylabel in zip(
    axes,
    ["Grid-interception loss fraction", "Target-leakage fraction $\\lambda$"],
    ["Grid loss (%)", "Leakage (%)"]
):
    ax.axhline(LOSS_LIMIT * 100, color="red", linestyle="--",
               linewidth=1.0, label="10 % limit")
    ax.set_xlabel(r"Aperture radius $r_a$ (mm)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

axes[0].legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig4_loss_budget.png"), dpi=200)
plt.close(fig)
print("Saved fig4_loss_budget.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 5 -- Uniformity stats (peak-to-avg, edge peaking) vs ra
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))

for ni_v in nis:
    sub = df[df["ni_e18"] == ni_v].sort_values("ra_mm")
    c   = ni_colors[ni_v]
    lbl = rf"$n_i={ni_v:.1f}$"
    axes[0].plot(sub["ra_mm"], sub["peak_to_avg"],
                 marker="o", markersize=4, color=c, label=lbl)
    axes[1].plot(sub["ra_mm"], sub["edge_peaking"],
                 marker="o", markersize=4, color=c, label=lbl)

axes[0].set_xlabel(r"Aperture radius $r_a$ (mm)")
axes[0].set_ylabel(r"Peak-to-average ratio $\mathcal{P}$")
axes[0].set_title("Peak-to-average uniformity")
axes[0].axhline(1.0, color="gray", linestyle=":", linewidth=0.8)
axes[0].legend(fontsize=7, ncol=2)

axes[1].set_xlabel(r"Aperture radius $r_a$ (mm)")
axes[1].set_ylabel(r"Edge-peaking index $\mathcal{E}$")
axes[1].set_title("Edge-peaking index")
axes[1].axhline(1.0, color="gray", linestyle=":", linewidth=0.8,
                label=r"$\mathcal{E}=1$ (flat)")

for ax in axes:
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig5_uniformity.png"), dpi=200)
plt.close(fig)
print("Saved fig5_uniformity.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 6 -- Radial profile at candidate design points
# ═══════════════════════════════════════════════════════════════════════════

# Identify candidates: I_del >= 4 A, grid_loss <= 10%, no tube loss
candidates = df[
    (df["I_del_clean"] >= R_TARGET_A) &
    (df["grid_loss_frac"] <= LOSS_LIMIT) &
    (~df["lost_to_tube"])
].copy()

print(f"\nDesign-space candidates (I_del>=4A, loss<=10%, no tube hit): {len(candidates)}")
if len(candidates) > 0:
    # Best point: minimum ni then minimum ra (lowest stress on source)
    candidates_sorted = candidates.sort_values(["ni_e18", "ra_mm"])
    print(candidates_sorted[["ni_e18", "ra_mm", "I_del_clean",
                              "grid_loss_frac", "leakage_frac",
                              "peak_to_avg"]].to_string(index=False))
    best = candidates_sorted.iloc[0]
    best_key = (best["ni_m3"], best["ra_mm"] / 1e3)
else:
    print("No candidates meet all constraints -- relaxing grid_loss limit")
    best_key = df.loc[df["I_del_clean"].idxmax(), ["ni_m3"]].values[0], \
               df.loc[df["I_del_clean"].idxmax(), ["ra_mm"]].values[0] / 1e3
    best = df.loc[df["I_del_clean"].idxmax()]

# Load full profile for best point
best_record = records[best_key]
bins   = best_record["profile_bins"]
y_mid  = np.array([(b["y_lo_m"] + b["y_hi_m"]) / 2 * 1e3 for b in bins])
I_Apm  = np.array([b["I_Apm"] for b in bins])
in_smp = np.array([b["in_sample"] for b in bins], dtype=bool)

# Normalise to mean in-sample density
mu_in = I_Apm[in_smp].mean() if in_smp.any() else 1.0
I_norm = I_Apm / mu_in

fig, ax = plt.subplots(figsize=(5.5, 4.0))
ax.bar(y_mid[in_smp],  I_norm[in_smp],  width=0.5, color="steelblue",
       alpha=0.85, label="In-sample")
ax.bar(y_mid[~in_smp], I_norm[~in_smp], width=0.5, color="salmon",
       alpha=0.70, label="Leakage")
ax.axvline(-R_SAMPLE_M * 1e3, color="k", linestyle="--", linewidth=0.9)
ax.axvline( R_SAMPLE_M * 1e3, color="k", linestyle="--", linewidth=0.9,
            label=rf"$R_s = {R_SAMPLE_M*1e3:.1f}$ mm")
ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)

ni_label = best["ni_e18"]
ra_label = best["ra_mm"]
ax.set_xlabel("Radial position $y$ (mm)")
ax.set_ylabel(r"Normalised current density $I / \langle I \rangle_\mathrm{in}$")
ax.set_title(rf"Target-plane profile: $n_i={ni_label:.1f}\times10^{{18}}$ m$^{{-3}}$, "
             rf"$r_a={ra_label:.1f}$ mm")
ax.legend(fontsize=8)
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig6_best_profile.png"), dpi=200)
plt.close(fig)
print("Saved fig6_best_profile.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 7 -- Sidewall / tube-hit diagnostic map
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(5.5, 4.0))
hit_df   = df[df["lost_to_tube"]]
nohit_df = df[~df["lost_to_tube"]]

ax.scatter(nohit_df["ra_mm"], nohit_df["ni_e18"], marker="o",
           s=40, color="steelblue", alpha=0.8, label="Beam contained")
ax.scatter(hit_df["ra_mm"],   hit_df["ni_e18"],   marker="x",
           s=60, color="crimson", linewidths=1.5, label="Hits tube wall")

ax.set_xlabel(r"Aperture radius $r_a$ (mm)")
ax.set_ylabel(r"Ion density $n_i$ ($\times 10^{18}$ m$^{-3}$)")
ax.set_title("Sidewall containment map ($R_s = 0.300$ m)")
ax.legend(fontsize=9)
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
fig.tight_layout()
fig.savefig(os.path.join(FIGS_DIR, "fig7_sidewall_map.png"), dpi=200)
plt.close(fig)
print("Saved fig7_sidewall_map.png")

# ═══════════════════════════════════════════════════════════════════════════
# Fig 8 -- Radial grid-loss profile (per-beamlet; only if new diagnostic data)
# ═══════════════════════════════════════════════════════════════════════════
records_with_pb = {k: v for k, v in records.items() if v["pb_radial"]}

if records_with_pb:
    mid_ra = sorted({k[1] for k in records_with_pb.keys()})[len(records_with_pb) // 2]
    pb_ni_groups = {}
    for (ni, ra), rec in records_with_pb.items():
        if ra == mid_ra:
            pb_ni_groups[ni] = rec["pb_radial"]

    if pb_ni_groups:
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
        ni_vals_sorted = sorted(pb_ni_groups.keys())
        cmap = plt.cm.plasma
        for idx, ni_val in enumerate(ni_vals_sorted):
            pb = pb_ni_groups[ni_val]
            r_mm   = np.array([p["r_m"] * 1e3 for p in pb])
            lf_pct = np.array([p["loss_frac"] * 100 for p in pb])
            order  = np.argsort(r_mm)
            color  = cmap(idx / max(1, len(ni_vals_sorted) - 1))
            ax.plot(r_mm[order], lf_pct[order], "o-", color=color, markersize=4,
                    label=rf"$n_i = {ni_val/1e18:.1f}\times10^{{18}}$ m$^{{-3}}$")

        ax.axhline(LOSS_LIMIT * 100, color="k", linestyle="--", linewidth=0.9,
                   label=f"{LOSS_LIMIT*100:.0f}% loss limit")
        ax.set_xlabel(r"Beamlet radial offset $|y_s|$ (mm)")
        ax.set_ylabel("Per-beamlet accel-grid loss (%)")
        ax.set_title(rf"Radial loss profile ($r_a = {mid_ra*1e3:.1f}$ mm, $R_s = 0.300$ m)")
        ax.legend(fontsize=8, ncol=2)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS_DIR, "fig8_radial_loss.png"), dpi=200)
        plt.close(fig)
        print("Saved fig8_radial_loss.png")
else:
    print("Skipping fig8_radial_loss.png: no per-beamlet data in current results "
          "(re-run sweep after updating multi_grid_2d_curved.cpp)")

# ═══════════════════════════════════════════════════════════════════════════
# CSV summary
# ═══════════════════════════════════════════════════════════════════════════
csv_path = os.path.join(FIGS_DIR, "sweep_summary.csv")
df_out = df.drop(columns=["has_per_beamlet"], errors="ignore")
df_out.to_csv(csv_path, index=False)
print(f"\nSaved sweep_summary.csv ({len(df_out)} rows)")
print("Done.")
