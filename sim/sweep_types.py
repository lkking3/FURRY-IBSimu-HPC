"""
sweep_types.py

Core types for parametric sweeps across IBSimu simulation variables.

Usage in a case file
--------------------
    from sweep_types import SweepAxis, SweepSpec, linspace, arange

    def make_case(**overrides):
        AP_RADIUS_M = overrides.get("AP_RADIUS_M", 0.0015)
        ION_N       = overrides.get("ION_N",       1e16)
        ...
        return SimulationCase(...)

    CASE = make_case()

    SWEEP = SweepSpec(
        axes=[
            SweepAxis("AP_RADIUS_M", linspace(0.001, 0.003, 5)),
            SweepAxis("ION_N",       [1e16, 2e16, 5e16]),
        ],
        mode="product",   # "product" or "zip"
        reduce=True,      # run reduce_best after all tasks finish
    )

Usage in a pipeline / sweep runner
-----------------------------------
    from sweep_types import SweepSpec

    sweep_spec = namespace["SWEEP"]          # loaded via exec()
    combos     = sweep_spec.combinations()   # list[dict[str, float]]
    manifest   = sweep_spec.manifest()       # list[dict] — JSON-serialisable

    # Write manifest to disk, then sbatch --array=0-<N-1>
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Sequence


# ---------------------------------------------------------------------------
# Convenience constructors
# ---------------------------------------------------------------------------

def linspace(start: float, stop: float, n: int) -> List[float]:
    """Return *n* evenly-spaced values from *start* to *stop* (inclusive)."""
    if n < 1:
        raise ValueError(f"linspace: n must be >= 1, got {n}")
    if n == 1:
        return [float(start)]
    step = (stop - start) / (n - 1)
    return [start + i * step for i in range(n)]


def arange(start: float, stop: float, step: float) -> List[float]:
    """Return values from *start* up to (but not including) *stop* spaced by *step*.

    Mirrors numpy.arange semantics.  Uses integer counting internally to avoid
    floating-point drift.
    """
    if step == 0:
        raise ValueError("arange: step must be non-zero")
    n = max(0, math.ceil((stop - start) / step))
    return [start + i * step for i in range(n)]


# ---------------------------------------------------------------------------
# Core types
# ---------------------------------------------------------------------------

@dataclass
class SweepAxis:
    """One dimension of a parameter sweep.

    Parameters
    ----------
    name:
        The keyword argument name that will be passed to ``make_case(**overrides)``.
        This is also the key in the override dict written to the sweep manifest.
    values:
        Sequence of values to sweep over.  All values must be JSON-serialisable
        (floats, ints, or strings).
    label:
        Optional human-readable axis label used in plot titles / filenames.
        Defaults to *name*.
    """
    name: str
    values: List[Any]
    label: str = ""

    def __post_init__(self) -> None:
        self.values = list(self.values)
        if not self.values:
            raise ValueError(f"SweepAxis '{self.name}': values must be non-empty")
        if not self.label:
            self.label = self.name


@dataclass
class SweepSpec:
    """Full specification for a parametric sweep.

    Parameters
    ----------
    axes:
        List of :class:`SweepAxis` objects defining the sweep dimensions.
    mode:
        ``"product"`` — Cartesian product (all combinations, default).
        ``"zip"``     — Element-wise pairing; all axes must have the same length.
    reduce:
        If ``True``, a reduce_best post-processing job is submitted after all
        array tasks complete (via SLURM ``--dependency=afterok``).
    tag:
        Optional string appended to ``RUN_TAG`` in every task to identify this
        sweep family in the results directory.
    """
    axes: List[SweepAxis] = field(default_factory=list)
    mode: Literal["product", "zip"] = "product"
    reduce: bool = False
    tag: str = ""

    def __post_init__(self) -> None:
        if not self.axes:
            raise ValueError("SweepSpec: at least one axis is required")
        if self.mode not in ("product", "zip"):
            raise ValueError(f"SweepSpec: mode must be 'product' or 'zip', got '{self.mode}'")
        if self.mode == "zip":
            lengths = {ax.name: len(ax.values) for ax in self.axes}
            unique  = set(lengths.values())
            if len(unique) > 1:
                raise ValueError(
                    f"SweepSpec mode='zip': all axes must have the same length, "
                    f"got {lengths}"
                )

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def total(self) -> int:
        """Total number of parameter combinations."""
        if self.mode == "product":
            result = 1
            for ax in self.axes:
                result *= len(ax.values)
            return result
        else:  # zip
            return len(self.axes[0].values)

    def combinations(self) -> List[Dict[str, Any]]:
        """Return a list of override dicts, one per combination.

        Each dict maps axis name → value and can be passed directly to
        ``make_case(**combo)``.
        """
        if self.mode == "product":
            keys   = [ax.name   for ax in self.axes]
            combos = list(itertools.product(*[ax.values for ax in self.axes]))
            return [dict(zip(keys, combo)) for combo in combos]
        else:  # zip
            return [
                {ax.name: ax.values[i] for ax in self.axes}
                for i in range(self.total)
            ]

    def manifest(self) -> List[Dict[str, Any]]:
        """Return the sweep manifest: a JSON-serialisable list of task records.

        Each record contains:
        - ``task_id`` (int)      — 0-based SLURM array index
        - ``overrides`` (dict)   — parameter overrides to pass to make_case
        """
        combos = self.combinations()
        return [
            {"task_id": i, "overrides": combo}
            for i, combo in enumerate(combos)
        ]

    # ------------------------------------------------------------------
    # Human-readable summary
    # ------------------------------------------------------------------

    def summary(self) -> str:
        lines = [f"SweepSpec  mode={self.mode}  total={self.total}  reduce={self.reduce}"]
        for ax in self.axes:
            vals = ax.values
            if len(vals) <= 6:
                val_str = str(vals)
            else:
                val_str = f"[{vals[0]}, {vals[1]}, … {vals[-2]}, {vals[-1]}]  ({len(vals)} values)"
            lines.append(f"  {ax.name:30s}  {val_str}")
        return "\n".join(lines)
