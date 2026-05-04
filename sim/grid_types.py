"""
grid_types.py

Shared geometry types and grid serialization for all IBSimu grid pipelines.

Import from here rather than from a specific pipeline module:

    from grid_types import Aperture, Chamfer, GridDefinition, SimulationCase
    from grid_types import serialize_grid_stack, load_case

Both multi_grid_pipeline and curved_grid_pipeline depend on these types.
Keeping them here means each pipeline is self-contained — neither needs to
import from the other.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass(frozen=True)
class Aperture:
    offset_m: float
    radius_m: float


@dataclass(frozen=True)
class Chamfer:
    depth_m: float = 0.0
    angle_deg: float = 0.0


@dataclass(frozen=True)
class GridDefinition:
    name: str
    voltage_v: float
    thickness_m: float
    gap_after_m: float = 0.0
    apertures: List[Aperture] = field(default_factory=list)
    upstream_chamfer: Chamfer = field(default_factory=Chamfer)
    downstream_chamfer: Chamfer = field(default_factory=Chamfer)
    mirror: bool = False


@dataclass(frozen=True)
class SimulationCase:
    grids: List[GridDefinition]
    env: Dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Serialization helpers
# ---------------------------------------------------------------------------

def _format_num(value: float) -> str:
    return f"{value:.12g}"


def _serialize_apertures(apertures: List[Aperture]) -> str:
    if not apertures:
        raise ValueError("each grid needs at least one aperture")
    return ",".join(
        f"{_format_num(ap.offset_m)}:{_format_num(ap.radius_m)}"
        for ap in apertures
    )


def _serialize_chamfer(chamfer: Chamfer) -> str:
    return f"{_format_num(chamfer.depth_m)}:{_format_num(chamfer.angle_deg)}"


def serialize_grid_stack(grids: List[GridDefinition]) -> str:
    """Encode a list of GridDefinitions into the GRID_STACK env-var format
    consumed by the C++ solver binaries."""
    if not 2 <= len(grids) <= 4:
        raise ValueError("multi-grid solver expects 2 to 4 grids")
    encoded = []
    for grid in grids:
        encoded.append("|".join([
            grid.name,
            _format_num(grid.voltage_v),
            _format_num(grid.thickness_m),
            _format_num(grid.gap_after_m),
            _serialize_apertures(grid.apertures),
            _serialize_chamfer(grid.upstream_chamfer),
            _serialize_chamfer(grid.downstream_chamfer),
            "1" if grid.mirror else "0",
        ]))
    return ";".join(encoded)


# ---------------------------------------------------------------------------
# Case loading
# ---------------------------------------------------------------------------

def load_case(
    config_path: Path,
    extra_namespace: Optional[Dict] = None,
) -> SimulationCase:
    """Execute a Python case file and return the SimulationCase it defines.

    The file must assign ``CASE = SimulationCase(...)``.

    Parameters
    ----------
    config_path:
        Path to the ``.py`` case file.
    extra_namespace:
        Additional names to inject before exec — used by pipeline-specific
        loaders that need to expose extra types (e.g. curved-grid types).
    """
    namespace: Dict = {
        "Aperture": Aperture,
        "Chamfer": Chamfer,
        "GridDefinition": GridDefinition,
        "SimulationCase": SimulationCase,
    }
    if extra_namespace:
        namespace.update(extra_namespace)

    code = compile(config_path.read_text(encoding="utf-8"), str(config_path), "exec")
    exec(code, namespace)

    case = namespace.get("CASE")
    if case is None or not hasattr(case, "grids") or not hasattr(case, "env"):
        raise TypeError(f"{config_path} must define CASE = SimulationCase(...)")
    if isinstance(case, SimulationCase):
        return case
    return SimulationCase(grids=list(case.grids), env=dict(case.env))
