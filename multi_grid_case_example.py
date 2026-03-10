from multi_grid_pipeline import Aperture, Chamfer, GridDefinition, SimulationCase

CASE = SimulationCase(
    grids=[
        GridDefinition(
            name='screen',
            voltage_v=0.0,
            thickness_m=0.005,
            gap_after_m=0.002,
            apertures=[Aperture(0.0, 0.002)],
            upstream_chamfer=Chamfer(0.0, 0.0),
            downstream_chamfer=Chamfer(0.0, 0.0),
            mirror=True,
        ),
        GridDefinition(
            name='accel',
            voltage_v=-15000.0,
            thickness_m=0.005,
            gap_after_m=0.0015,
            apertures=[Aperture(0.0, 0.002)],
            upstream_chamfer=Chamfer(0.002, 5.0),
            downstream_chamfer=Chamfer(0.0026, 5.0),
            mirror=True,
        ),
        GridDefinition(
            name='focus',
            voltage_v=5000.0,
            thickness_m=0.005,
            gap_after_m=0.0,
            apertures=[Aperture(0.0, 0.004)],
            upstream_chamfer=Chamfer(5, 0.002),
            downstream_chamfer=Chamfer(0.0026, 5.0),
            mirror=True,
        ),
    ],
    env={
        'RUN_SOLVE': '1',
        'WRITE_PNG': '1',
        'PNG_NAME': 'multi_grid_geom.png',
        'RESULTS_DIR': 'results',
        'RUN_PREFIX': 'multi_grid_test',
        'ENABLE_IONS': '1',
        'X_RIGHT_M': '0.02',
        'X_RIGHT_PHYS_M': '0.55',
        'YBOX_M': '0.086',
        'TUBE_WALL_T_M': '0.0002',
    },
)
