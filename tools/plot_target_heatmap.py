#!/usr/bin/env python3
"""
plot_target_heatmap.py -- 2-D current-density heatmap on the target plane from
the per-trajectory dump (results_dir/target_profile.dat) written by
grid_3d_curved (WRITE_TARGET). Half model (x>=0) is mirrored across x=0 for the
full circular surface.

  python tools/plot_target_heatmap.py <results_dir> [--bin-mm 1.5] [--target-in 2.5]
"""
import sys, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

ap=argparse.ArgumentParser()
ap.add_argument('results_dir')
ap.add_argument('--bin-mm', type=float, default=1.5)
ap.add_argument('--target-in', type=float, default=2.5)   # target diameter, inches
ap.add_argument('--out', default=None)
a=ap.parse_args()

dat=f"{a.results_dir}/target_profile.dat"
d=np.loadtxt(dat)
x,y,I = d[:,0], d[:,1], d[:,2]
# mirror the half model across x=0 -> full disk
x=np.concatenate([x,-x]); y=np.concatenate([y,y]); I=np.concatenate([I,I])

R=max(np.hypot(x,y).max()*1.05, a.target_in*0.0254/2*1.1)
bw=a.bin_mm*1e-3
nb=int(2*R/bw)
edges=np.linspace(-R,R,nb+1)
H,_,_=np.histogram2d(x,y,bins=[edges,edges],weights=I)   # summed current per cell [A]
J=H/(bw*bw)                                              # current density [A/m^2]
Jt=J.T

rt=a.target_in*0.0254/2.0
inside=np.hypot(x,y)<=rt
print(f"total current (full) = {I.sum():.4f} A ; within {a.target_in}in (r={rt*1e3:.1f}mm) = "
      f"{I[inside].sum():.4f} A ({100*I[inside].sum()/max(I.sum(),1e-30):.1f}%)")
print(f"peak current density = {J.max():.3e} A/m^2")

fig,ax=plt.subplots(figsize=(7.5,6.5))
im=ax.imshow(Jt, origin='lower', extent=[-R*1e3,R*1e3,-R*1e3,R*1e3],
             cmap='inferno', aspect='equal')
th=np.linspace(0,2*np.pi,200)
ax.plot(rt*1e3*np.cos(th), rt*1e3*np.sin(th), 'c--', lw=1.6, label=f'{a.target_in}" target')
ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
ax.set_title(f'Current density on target plane  (total {I.sum():.2f} A)')
ax.legend(loc='upper right')
cb=fig.colorbar(im, ax=ax); cb.set_label('J  (A/m$^2$)')
out=a.out or f"{a.results_dir}/target_heatmap.png"
plt.tight_layout(); plt.savefig(out, dpi=130); print("saved", out)
