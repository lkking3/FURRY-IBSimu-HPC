from __future__ import annotations

import argparse
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


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


def _format_num(value: float) -> str:
    return f"{value:.12g}"


def _serialize_apertures(apertures: List[Aperture]) -> str:
    if not apertures:
        raise ValueError('each grid needs at least one aperture')
    return ','.join(f"{_format_num(ap.offset_m)}:{_format_num(ap.radius_m)}" for ap in apertures)


def _serialize_chamfer(chamfer: Chamfer) -> str:
    return f"{_format_num(chamfer.depth_m)}:{_format_num(chamfer.angle_deg)}"


def serialize_grid_stack(grids: List[GridDefinition]) -> str:
    if not 2 <= len(grids) <= 4:
        raise ValueError('multi-grid solver expects 2 to 4 grids')
    encoded = []
    for grid in grids:
        encoded.append('|'.join([
            grid.name,
            _format_num(grid.voltage_v),
            _format_num(grid.thickness_m),
            _format_num(grid.gap_after_m),
            _serialize_apertures(grid.apertures),
            _serialize_chamfer(grid.upstream_chamfer),
            _serialize_chamfer(grid.downstream_chamfer),
            '1' if grid.mirror else '0',
        ]))
    return ';'.join(encoded)


def load_case(config_path: Path) -> SimulationCase:
    namespace = {
        'Aperture': Aperture,
        'Chamfer': Chamfer,
        'GridDefinition': GridDefinition,
        'SimulationCase': SimulationCase,
    }
    code = compile(config_path.read_text(encoding='utf-8'), str(config_path), 'exec')
    exec(code, namespace)
    case = namespace.get('CASE')
    if case is None or not hasattr(case, 'grids') or not hasattr(case, 'env'):
        raise TypeError(f'{config_path} must define CASE = SimulationCase(...)')
    if isinstance(case, SimulationCase):
        return case
    return SimulationCase(grids=list(case.grids), env=dict(case.env))


def build_env(case: SimulationCase) -> Dict[str, str]:
    env = os.environ.copy()
    env['GRID_STACK'] = serialize_grid_stack(case.grids)
    for key, value in case.env.items():
        env[key] = str(value)
    return env


def main() -> int:
    parser = argparse.ArgumentParser(description='Run the optional multi-grid IBSimu pipeline.')
    parser.add_argument('--config', default='multi_grid_case_example.py', help='Python config file defining CASE = SimulationCase(...)')
    parser.add_argument('--binary', default='multi_grid_2d', help='Path to the multi-grid solver binary')
    parser.add_argument('--dry-run', action='store_true', help='Print the generated GRID_STACK and exit')
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
        print(env['GRID_STACK'])
        return 0

    if not binary_path.exists():
        raise FileNotFoundError(f'multi-grid solver not found: {binary_path}')

    subprocess.run([str(binary_path)], check=True, cwd=repo_root, env=env)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
