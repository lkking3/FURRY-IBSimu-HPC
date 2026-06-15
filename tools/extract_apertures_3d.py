#!/usr/bin/env python3
"""
extract_apertures_3d.py
Extract per-aperture emission geometry (centre, inward normal, radius) from a
curved grid STL, for Bohm-injection sites in sim/grid_3d_curved.cpp.

Output (one aperture per line, SI metres, simulation coords = STL_mm / 1000):
    x  y  z   nx ny nz   r_ap
(x,y,z) upstream surface point; (nx,ny,nz) unit inward (downstream) normal.
"""
from __future__ import annotations
import argparse, sys
import numpy as np
try:
    import trimesh
except ImportError:
    sys.exit("need: pip install trimesh scipy networkx shapely rtree")
try:
    from sklearn.cluster import DBSCAN
except ImportError:
    sys.exit("need: pip install scikit-learn")


def sphere_fit(V):
    A = np.c_[2 * V, np.ones(len(V))]
    b = (V ** 2).sum(1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    C = sol[:3]
    return C, float(np.sqrt(sol[3] + C @ C))


def extract(stl_path, grid_spacing, ray_back, eps, min_samples):
    m = trimesh.load(stl_path, process=True)
    V = m.vertices
    # Fit curvature centre on UPSTREAM face only (normal points -z), so the flat
    # back / rim / bore walls do not bias the radius. Axis forced to x=y=0.
    fn = m.face_normals
    up = fn[:, 2] < -0.3
    up_vidx = np.unique(m.faces[up].ravel())
    Vup = V[up_vidx] if up.sum() > 100 else V
    C_fit, R = sphere_fit(Vup)
    C = np.array([0.0, 0.0, C_fit[2]])
    rmax = 0.99 * np.hypot(V[:, 0], V[:, 1]).max()

    x_lo, x_hi = V[:, 0].min(), V[:, 0].max()
    y_lo, y_hi = V[:, 1].min(), V[:, 1].max()
    xs = np.arange(x_lo + 0.3, x_hi, grid_spacing)
    ys = np.arange(y_lo + 0.3, y_hi, grid_spacing)
    pts, dirs = [], []
    for x in xs:
        for y in ys:
            rr = (x - C[0]) ** 2 + (y - C[1]) ** 2
            if rr > rmax * rmax:
                continue
            disc = R * R - rr
            if disc <= 0:
                continue
            z = C[2] - np.sqrt(disc)
            p = np.array([x, y, z])
            n = C - p
            n /= np.linalg.norm(n)
            pts.append(p - n * ray_back)
            dirs.append(n)
    pts = np.asarray(pts)
    dirs = np.asarray(dirs)
    if len(pts) == 0:
        sys.exit("No candidates; check STL units/orientation.")

    hit = m.ray.intersects_any(ray_origins=pts, ray_directions=dirs)
    holes = pts[~hit] + dirs[~hit] * ray_back
    if len(holes) == 0:
        sys.exit("No apertures detected.")

    lab = DBSCAN(eps=eps, min_samples=min_samples).fit(holes[:, :2]).labels_
    apertures = []
    for k in sorted(set(lab)):
        if k == -1:
            continue
        pk = holes[lab == k]
        c2 = pk[:, :2].mean(0)
        rr = c2[0] ** 2 + c2[1] ** 2
        disc = R * R - rr
        if disc <= 0:
            continue
        z = C[2] - np.sqrt(disc)
        centre = np.array([c2[0], c2[1], z])
        n = C - centre
        n /= np.linalg.norm(n)
        rad = float(np.sqrt(((pk[:, :2] - c2) ** 2).sum(1).mean())) * 1.3
        apertures.append((centre, n, rad, len(pk)))

    if apertures:
        med = np.median([a[2] for a in apertures])
        apertures = [a for a in apertures if a[2] <= 2.5 * med]
    return apertures, C, R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stl", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--grid-spacing", type=float, default=1.0)
    ap.add_argument("--ray-back", type=float, default=1.5)
    ap.add_argument("--eps", type=float, default=2.0)
    ap.add_argument("--min-samples", type=int, default=2)
    ap.add_argument("--min-cluster", type=int, default=3)
    ap.add_argument("--unit-scale", type=float, default=1.0e-3)
    args = ap.parse_args()

    apertures, C, R = extract(args.stl, args.grid_spacing, args.ray_back,
                              args.eps, args.min_samples)
    apertures = [a for a in apertures if a[3] >= args.min_cluster]
    s = args.unit_scale
    with open(args.out, "w") as f:
        f.write("# apertures from %s\n" % args.stl)
        f.write("# x y z nx ny nz r_ap (SI metres, sim coords)\n")
        f.write("# curv centre mm %.3f %.3f %.3f R %.3f\n" % (C[0], C[1], C[2], R))
        f.write("# n_apertures %d\n" % len(apertures))
        for centre, n, rad, _ in apertures:
            f.write("%.8e %.8e %.8e  %.8e %.8e %.8e  %.8e\n" %
                    (centre[0]*s, centre[1]*s, centre[2]*s, n[0], n[1], n[2], rad*s))
    radii = np.array([a[2] for a in apertures]) * s
    print("wrote %d apertures -> %s" % (len(apertures), args.out))
    print("  r_ap m mean %.4e min %.4e max %.4e" % (radii.mean(), radii.min(), radii.max()))


if __name__ == "__main__":
    main()
