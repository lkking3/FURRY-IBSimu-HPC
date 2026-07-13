#!/usr/bin/env python
# =============================================================================
# pa_zmap.py -- on-target uniformity (peak-to-average) and absolute peak power
# density vs axial position z, across a set of runs at different space-charge
# compensation levels. Builds the (compensation x z) operating map.
#
# Each run's beam_density.vtk slice is calibrated to the conserved beam current
# (lossless transport -> integral J dA = I at every plane), so peak power
# density = max(J) * V is absolute (W/cm^2), not relative.
#
#   pvpython/python tools/pa_zmap.py I_A V_kV run1:scfac1 run2:scfac2 ...
#   e.g.  python tools/pa_zmap.py 4.326 10 results_3d_test_14:0.0005 \
#                                          results_3d_test_15:0.005
# scfac = SC_FACTOR (residual fraction); compensation% = (1-scfac)*100.
# Writes tools/../outputs and a pa_zmap.csv you can extend.
# =============================================================================
import sys, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
try:    from scipy.ndimage import gaussian_filter as gf
except Exception:  gf=lambda a,s:a

def load_vtk_struct(fn):
    raw=open(fn,'rb').read(); key=b'LOOKUP_TABLE default\n'
    off=raw.find(key)+len(key)
    hdr=raw[:off].decode('latin1')
    dims=[l for l in hdr.splitlines() if l.startswith('DIMENSIONS')][0].split()[1:]
    org =[l for l in hdr.splitlines() if l.startswith('ORIGIN')][0].split()[1:]
    spc =[l for l in hdr.splitlines() if l.startswith('SPACING')][0].split()[1:]
    nx,ny,nz=map(int,dims); ox,oy,oz=map(float,org); sp=float(spc[0])
    return raw,off,nx,ny,nz,ox,oy,oz,sp

def analyse(run, I, V, r_target=0.5*2.5*0.0254, zlo=0.05, zhi=0.212, smooth_mm=1.0):
    raw,off,nx,ny,nz,ox,oy,oz,sp=load_vtk_struct(os.path.join(run,'beam_density.vtk'))
    xs=ox+np.arange(nx)*sp; ys=oy+np.arange(ny)*sp
    X,Y=np.meshgrid(xs,ys); R=np.hypot(X,Y); dA=sp*sp; sig=smooth_mm*1e-3/sp
    def slc(k): return np.frombuffer(raw[off+k*nx*ny*4:off+(k+1)*nx*ny*4],dtype='>f4').astype(np.float64).reshape(ny,nx)
    Z=[];PA=[];PPD=[];FR=[];RM=[]
    for k in range(int(zlo/sp),min(nz,int(zhi/sp))):
        s=slc(k); tot=s.sum()
        if tot<=0: continue
        # calibrate to conserved current: J = density * I / (sum density * dA)
        J=s*(I/(tot*dA)); Js=gf(J,sig)
        inc=R<=r_target
        foot=inc&(Js>0.05*Js[inc].max())
        if foot.sum()<20: continue
        PA.append(np.percentile(Js[foot],99)/Js[foot].mean())
        PPD.append(Js[inc].max()*V/1e4)                 # W/cm^2 (J[A/m^2]*V[V] -> W/m^2 /1e4)
        FR.append(J[inc].sum()*dA/I*100)
        RM.append(R[s>0.01*s.max()].max()*1e3)
        Z.append(k*sp)
    return map(np.array,(Z,PA,PPD,FR,RM))

def main():
    I=float(sys.argv[1]); V=float(sys.argv[2])*1e3
    runs=[a.split(':') for a in sys.argv[3:]]
    outdir=os.path.join(os.path.dirname(__file__),'..','outputs'); os.makedirs(outdir,exist_ok=True)
    import csv
    cw=csv.writer(open(os.path.join(outdir,'pa_zmap.csv'),'w',newline=''))
    cw.writerow(['run','SC_FACTOR','comp_pct','z_m','P2A_99','peak_Wcm2','pct_on_target','r_max_mm'])
    cmap=plt.cm.viridis; sfs=[float(s) for _,s in runs]
    fig,ax=plt.subplots(1,2,figsize=(13,5))
    for (run,sf),c in zip(runs,np.linspace(0,1,max(len(runs),2))):
        sf=float(sf); Z,PA,PPD,FR,RM=analyse(run,I,V)
        comp=(1-sf)*100; col=cmap(c)
        lab='%.3g%% comp (SC_FACTOR %.4g)'%(comp,sf)
        ax[0].plot(Z*1e3,PA,'o-',ms=3,color=col,label=lab)
        ax[1].plot(Z*1e3,PPD,'o-',ms=3,color=col,label=lab)
        for i in range(len(Z)):
            cw.writerow([os.path.basename(run.rstrip('/')),sf,comp,'%.4f'%Z[i],'%.3f'%PA[i],'%.1f'%PPD[i],'%.1f'%FR[i],'%.1f'%RM[i]])
    ax[0].set_xlabel('target z (mm)'); ax[0].set_ylabel('peak(99%)/avg on target'); ax[0].set_title('Uniformity'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
    ax[1].set_xlabel('target z (mm)'); ax[1].set_ylabel('peak power density (W/cm$^2$)'); ax[1].set_title('Thermal load'); ax[1].grid(alpha=.3); ax[1].set_yscale('log'); ax[1].legend(fontsize=8)
    plt.tight_layout(); out=os.path.join(outdir,'pa_zmap.png'); plt.savefig(out,dpi=140)
    print('saved',out,'and pa_zmap.csv')

if __name__=='__main__': main()
