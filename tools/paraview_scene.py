#!/usr/bin/env python
# =============================================================================
# paraview_scene.py  --  polished 3-D ParaView scene for a grid_3d_curved run.
#
# Shows:  screen + accel electrodes (CAD STL, scaled mm->m),
#         the FULL beam (the half model is already mirrored in the VTK writer)
#           as a log-scaled, opacity-modulated density volume,
#         beamlet trajectory lines,
#         an E-potential mid-plane slice,
#         the 2.5-inch-diameter target disk at the drift end.
#
# Run either way:
#   GUI :  ParaView -> View -> Python Shell ->
#              exec(open('tools/paraview_scene.py').read())
#   batch: pvpython tools/paraview_scene.py [RESULTS_DIR] [--screenshot out.png]
#
# Edit the CONFIG block for your paths / target plane, then re-run.
# =============================================================================
import os, sys
from paraview.simple import *

# ----------------------------- CONFIG ----------------------------------------
args      = [a for a in sys.argv[1:] if not a.startswith('-')]
RESULTS   = args[0] if args else 'results_3d'              # run output dir
GEOM_DIR  = 'curved_grid_3d/geometry'                      # electrode STLs (mm)
MM        = 1.0e-3                                         # STL mm -> sim metres
Z_CLIP    = 0.17                                           # beam termination plane (m); None to disable
Z_TARGET  = Z_CLIP if Z_CLIP is not None else 0.2134       # 2.5" target sits at the termination plane
TARGET_R  = 0.5 * 2.5 * 0.0254                             # 2.5 inch dia -> radius (m)
SHOT      = None
if '--screenshot' in sys.argv:
    i = sys.argv.index('--screenshot')
    SHOT = sys.argv[i+1] if i+1 < len(sys.argv) else 'scene.png'

def p(f):  return os.path.join(RESULTS, f)

# Wipe any existing pipeline so re-running in the same session doesn't stack a
# new (clipped) volume on top of the old un-clipped one.
for _s in list(GetSources().values()):
    Delete(_s)

view = GetActiveViewOrCreate('RenderView')
view.OrientationAxesVisibility = 1
try:    view.BackgroundColorMode = 'Gradient'
except Exception: pass
view.Background  = [0.10, 0.11, 0.14]
view.Background2 = [0.02, 0.02, 0.03]

# ----------------------------- electrodes (CAD) ------------------------------
def load_electrode(stl, color):
    if not os.path.exists(stl):
        print('skip (missing):', stl); return None
    r = STLReader(FileNames=[stl])
    t = Transform(Input=r); t.Transform.Scale = [MM, MM, MM]   # mm -> m
    d = Show(t, view); d.Representation = 'Surface'
    d.DiffuseColor = color; d.Opacity = 1.0; d.Specular = 0.3
    return t

load_electrode(os.path.join(GEOM_DIR, 'screen_pos.stl'), [0.75, 0.78, 0.82])  # screen: light grey
load_electrode(os.path.join(GEOM_DIR, 'accel_pos.stl'),  [0.55, 0.40, 0.30])  # accel:  bronze

# ----------------------------- beam density (volume) -------------------------
if os.path.exists(p('beam_density.vtk')):
    bd = LegacyVTKReader(FileNames=[p('beam_density.vtk')])
    src = bd
    if Z_CLIP is not None:                       # terminate the beam at z = Z_CLIP
        clip = Clip(Input=bd)
        clip.ClipType = 'Plane'
        clip.ClipType.Origin = [0.0, 0.0, Z_CLIP]
        clip.ClipType.Normal = [0.0, 0.0, 1.0]
        clip.Invert = 1                          # keep z < Z_CLIP
        src = clip
    dd = Show(src, view); dd.SetRepresentationType('Volume')
    ColorBy(dd, ('POINTS', 'beam_density'))
    ctf = GetColorTransferFunction('beam_density')
    otf = GetOpacityTransferFunction('beam_density')
    dd.RescaleTransferFunctionToDataRange(True, False)
    ctf.ApplyPreset('Inferno (matplotlib)', True)
    # log scale compresses the wide density range so the faint halo shows
    try:
        ctf.UseLogScale = 1; otf.UseLogScale = 1
    except Exception: pass
    rng = bd.PointData['beam_density'].GetRange()
    lo  = max(rng[1] * 1e-3, 1e-30); hi = max(rng[1], lo * 10)
    ctf.RescaleTransferFunction(lo, hi)
    otf.RescaleTransferFunction(lo, hi)
    # opacity ramp: faint at low density, solid core at high density
    otf.Points = [lo, 0.0, 0.5, 0.0,
                  lo*10,  0.04, 0.5, 0.0,
                  hi*0.3, 0.25, 0.5, 0.0,
                  hi,     0.85, 0.5, 0.0]
    try: dd.OpacityUnitDistanceScaling = 1
    except Exception: pass

# ----------------------------- trajectories (lines) --------------------------
if os.path.exists(p('beam_trajectories.vtk')):
    tj = LegacyVTKReader(FileNames=[p('beam_trajectories.vtk')])
    if Z_CLIP is not None:                       # terminate beamlets at z = Z_CLIP
        tjc = Clip(Input=tj)
        tjc.ClipType = 'Plane'
        tjc.ClipType.Origin = [0.0, 0.0, Z_CLIP]
        tjc.ClipType.Normal = [0.0, 0.0, 1.0]
        tjc.Invert = 1                           # keep z < Z_CLIP
        tj = tjc
    td = Show(tj, view); td.Representation = 'Surface'
    td.AmbientColor = [1.0, 0.85, 0.2]; td.DiffuseColor = [1.0, 0.85, 0.2]
    td.LineWidth = 1.0; td.Opacity = 0.35

# ----------------------------- E-pot mid-plane slice -------------------------
if os.path.exists(p('epot.vtk')):
    ep = LegacyVTKReader(FileNames=[p('epot.vtk')])
    sl = Slice(Input=ep); sl.SliceType = 'Plane'
    sl.SliceType.Origin = [0.0, 0.0, 0.05]; sl.SliceType.Normal = [0.0, 1.0, 0.0]
    sd = Show(sl, view); ColorBy(sd, ('POINTS', 'potential'))
    GetColorTransferFunction('potential').ApplyPreset('Cool to Warm', True)
    sd.Opacity = 0.5

# ----------------------------- 2.5" target disk ------------------------------
disk = Disk(); disk.InnerRadius = 0.0; disk.OuterRadius = TARGET_R
disk.CircumferentialResolution = 96
dt = Transform(Input=disk); dt.Transform.Translate = [0.0, 0.0, Z_TARGET]
# Disk normal is +z already (faces the beam); show translucent red + edges
dd2 = Show(dt, view); dd2.Representation = 'Surface With Edges'
dd2.DiffuseColor = [0.9, 0.15, 0.15]; dd2.Opacity = 0.35
dd2.EdgeColor = [1.0, 0.3, 0.3]
print('target: 2.5 in dia (R=%.4f m) at z=%.4f m' % (TARGET_R, Z_TARGET))

# ----------------------------- camera / render -------------------------------
Render()
ResetCamera(view)
cam = GetActiveCamera()
cam.Azimuth(35); cam.Elevation(20)
view.CameraViewUp = [0, 1, 0]
Render()

if SHOT:
    SaveScreenshot(SHOT, view, ImageResolution=[1920, 1080])
    print('saved', SHOT)
