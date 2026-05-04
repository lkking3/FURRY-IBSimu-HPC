from __future__ import annotations

import argparse
import datetime
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict

# Locate sim/ relative to this file so grid_types is importable regardless
# of the working directory or how this script is invoked.
_HERE = Path(__file__).resolve().parent
_SIM_DIR = _HERE.parent / "sim"
if _SIM_DIR.exists() and str(_SIM_DIR) not in sys.path:
    sys.path.insert(0, str(_SIM_DIR))

from grid_types import (
    Aperture,
    Chamfer,
    GridDefinition,
    SimulationCase,
    _format_num,
    serialize_grid_stack,
    load_case,
)


def _run_profile_plot(env: Dict[str, str], repo_root: Path) -> None:
    """Generate sample_profile_I_Apm.html next to the run's sample_diameter_profile.json.

    Reconstructs the output directory using the same naming logic as the C++ binary
    (RESULTS_DIR / RUN_PREFIX[_RUN_TAG]_RUN_STAMP[_jSLURM_JOB_ID]).

    Self-contained: reads the JSON with stdlib only and writes a standalone HTML
    using Plotly via CDN — no local plotly/pandas install required.
    """
    import json as _json

    results_base = env.get("RESULTS_DIR", "results")
    stamp        = env.get("RUN_STAMP", "")
    prefix       = env.get("RUN_PREFIX", "run")
    tag          = env.get("RUN_TAG", "")
    jobid        = env.get("SLURM_JOB_ID", "")

    runname = prefix
    if tag:
        runname += f"_{tag}"
    runname += f"_{stamp}"
    if jobid:
        runname += f"_j{jobid}"

    base_path = Path(results_base)
    if not base_path.is_absolute():
        base_path = repo_root / base_path
    outdir = base_path / runname

    json_file = outdir / "sample_diameter_profile.json"
    if not json_file.exists():
        print(f"[profile-plot] sample_diameter_profile.json not found at {outdir} — skipping.", flush=True)
        return

    out_html = outdir / "sample_profile_I_Apm.html"
    print(f"[profile-plot] generating {out_html.name} ...", flush=True)

    try:
        data = _json.loads(json_file.read_text())
        bins = data.get("bins", [])
        if not bins:
            print("[profile-plot] no bins in profile JSON — skipping.", flush=True)
            return

        xs     = [(b["y_lo_m"] + b["y_hi_m"]) * 0.5 * 1000.0 for b in bins]  # mm
        ys     = [b["I_Apm"] for b in bins]
        colors = ["royalblue" if b.get("in_sample", True) else "crimson" for b in bins]
        bw_mm  = data.get("bin_width_m", (xs[1] - xs[0]) / 1000.0 if len(xs) > 1 else 1.0)
        if isinstance(bw_mm, float):
            bw_mm = bw_mm * 1000.0

        stats  = data.get("stats", {})
        p2a    = stats.get("peak_to_avg", float("nan"))
        rms_nu = stats.get("rms_nonuniform", float("nan"))
        title  = (f"Sample Diameter Profile — I_Apm<br>"
                  f"<sub>peak/avg={p2a:.3f}  rms_nonuniform={rms_nu:.3f}  "
                  f"total_I={data.get('total_I_A', 0):.4f} A</sub>")

        bar_traces = []
        for label, colour in (("in sample", "royalblue"), ("outside sample", "crimson")):
            xi = [x for x, c in zip(xs, colors) if c == colour]
            yi = [y for y, c in zip(ys, colors) if c == colour]
            if xi:
                bar_traces.append(
                    f'{{"type":"bar","name":"{label}","x":{xi},"y":{yi},'
                    f'"marker":{{"color":"{colour}"}},"width":{bw_mm},"showlegend":true}}'
                )

        traces_js = ",".join(bar_traces)
        html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Sample Diameter Profile</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
</head>
<body>
<div id="plot" style="width:100%;height:600px;"></div>
<script>
var traces = [{traces_js}];
var layout = {{
  title: {{text: "{title}", font: {{size: 18}}}},
  xaxis: {{title: "y (mm)", tickfont: {{size: 16}}, titlefont: {{size: 18}}}},
  yaxis: {{title: "I (A/m)", tickfont: {{size: 16}}, titlefont: {{size: 18}}}},
  barmode: "overlay",
  bargap: 0,
  template: "plotly_white",
  legend: {{font: {{size: 14}}}}
}};
Plotly.newPlot("plot", traces, layout, {{responsive: true}});
</script>
</body>
</html>"""

        out_html.write_text(html, encoding="utf-8")
        print(f"[profile-plot] wrote {out_html}", flush=True)

    except Exception as exc:
        print(f"[profile-plot] ERROR: {exc}", flush=True)


def build_env(case: SimulationCase) -> Dict[str, str]:
    env = os.environ.copy()
    env["GRID_STACK"] = serialize_grid_stack(case.grids)
    for key, value in case.env.items():
        env[key] = str(value)
    return env


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the multi-grid IBSimu pipeline.")
    parser.add_argument("--config", default="multi_grid_case_example.py",
                        help="Python config file defining CASE = SimulationCase(...)")
    parser.add_argument("--binary", default="multi_grid_2d",
                        help="Path to the multi-grid solver binary")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the generated GRID_STACK and exit")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = repo_root / config_path
    binary_path = Path(args.binary)
    if not binary_path.is_absolute():
        binary_path = repo_root / binary_path

    case = load_case(config_path)
    env = build_env(case)

    if args.dry_run:
        print(env["GRID_STACK"])
        return 0

    if not binary_path.exists():
        raise FileNotFoundError(f"multi-grid solver not found: {binary_path}")

    # Pin RUN_STAMP now so Python and C++ agree on the output directory name.
    if not env.get("RUN_STAMP"):
        env["RUN_STAMP"] = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    subprocess.run([str(binary_path)], check=True, cwd=repo_root, env=env)

    # Auto-generate sample diameter profile plot if WRITE_PNG is active.
    if int(env.get("WRITE_PNG", "1")):
        _run_profile_plot(env, repo_root)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
