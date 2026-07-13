#!/usr/bin/env python
# =============================================================================
# scc_0d.py -- 0D rate-balance estimate of the space-charge-compensation (SCC)
# onset length for a DC ion beam drifting through its own fill gas.
#
# Question answered: how far does the beam propagate past the accel exit before
# its uncompensated ion space charge is dealt with?  Two processes remove it:
#
#   (1) NEUTRALISATION (charge exchange): fast ion -> fast neutral, which carries
#       no space charge.  Ion population decays as exp(-z / L_cx),
#       L_cx = 1 / (n_gas * sigma_cx).
#   (2) COMPENSATION (impact ionisation): the beam ionises the gas, and the freed
#       electrons are trapped in the -10 kV well and neutralise the remaining
#       ions.  Electrons build over L_ion = 1 / (n_gas * sigma_ion); equivalently
#       the DC buildup time tau_scc = 1 / (n_gas * sigma_ion * v_beam).
#
# Both scale as 1/pressure, so drift pressure sets where compensation kicks in.
#
# CROSS SECTIONS are the dominant uncertainty. Defaults below are representative
# for D+ (or H+) on D2/H2 at ~5 keV/amu; PIN from ORNL/Barnett (ALADDIN) for the
# final number. Refs: C.F. Barnett, ORNL-6086 "Redbook" (1990); Phys. Rev. 149,62
# (1966) H+ + H2 charge transfer, 2-50 keV.
#
#   python tools/scc_0d.py [E_keV] [m_amu] [sigma_cx_cm2] [sigma_ion_cm2]
#   default: 10 keV, m=2, sigma_cx=1.0e-16, sigma_ion=0.5e-16
# =============================================================================
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

E_keV   = float(sys.argv[1]) if len(sys.argv)>1 else 10.0
m_amu   = float(sys.argv[2]) if len(sys.argv)>2 else 2.0
sig_cx  = (float(sys.argv[3]) if len(sys.argv)>3 else 1.0e-16)*1e-4   # cm^2 -> m^2
sig_ion = (float(sys.argv[4]) if len(sys.argv)>4 else 0.5e-16)*1e-4
QE=1.602176634e-19; AMU=1.66053907e-27; kB=1.380649e-23; T=300.0
vb=np.sqrt(2*E_keV*1e3*QE/(m_amu*AMU))
Ldrift=0.170

def ngas(p_mTorr): return (p_mTorr*1e-3*133.322)/(kB*T)

print("beam: %.0f keV m=%g amu  v=%.3e m/s   drift=%.0f mm"%(E_keV,m_amu,vb,Ldrift*1e3))
print("sigma_cx=%.2e cm^2  sigma_ion=%.2e cm^2 (PIN from Barnett/ALADDIN)\n"%(sig_cx*1e4,sig_ion*1e4))
print(" P(mTorr)  n_gas(/m^3)  L_cx(neutralise) L_ion(compensate) tau_scc  ion-frac@drift-end")
for p in [10,20,40,60,100]:
    n=ngas(p); Lcx=1/(n*sig_cx); Lion=1/(n*sig_ion); tau=1/(n*sig_ion*vb)
    ionfrac=np.exp(-Ldrift/Lcx)
    print("  %5d   %.2e   %6.1f mm       %6.1f mm      %5.1f ns   %.3f"%(p,n,Lcx*1e3,Lion*1e3,tau*1e9,ionfrac))

# --- plot onset lengths vs pressure, with the sim's good/marginal offset band ---
P=np.linspace(5,100,120); n=ngas(P)
Lcx=1/(n*sig_cx)*1e3; Lion=1/(n*sig_ion)*1e3
fig,ax=plt.subplots(figsize=(8,5))
ax.plot(P,Lcx,color='C0',lw=2.4,label='L_cx  (ion neutralisation length)')
ax.plot(P,Lion,color='C3',lw=2.4,label='L_ion (compensation buildup length)')
ax.axhspan(0,3,color='green',alpha=.10); ax.axhline(3,color='green',ls=':',label='sim: filled (offset ≤3 mm)')
ax.axhline(10,color='orange',ls=':',label='sim: degraded (offset 10 mm)')
ax.axhline(Ldrift*1e3,color='gray',ls='--',lw=1,label='drift length (170 mm)')
ax.axvspan(20,60,color='gray',alpha=.06)
ax.set_xlabel('drift pressure (mTorr)'); ax.set_ylabel('onset length (mm)')
ax.set_yscale('log'); ax.set_ylim(1,1000); ax.grid(alpha=.3,which='both')
ax.set_title('SCC onset length vs drift pressure (D$_2$, %.0f keV D$^+$)'%E_keV); ax.legend(fontsize=8)
plt.tight_layout(); import os
out=os.path.join(os.path.dirname(__file__),'..','outputs','scc_onset_length.png')
plt.savefig(out,dpi=140); print("\nsaved",out)
