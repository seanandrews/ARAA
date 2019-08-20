import numpy as np
import os
import sys
from astropy.io import ascii
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# set up plot
fig = plt.figure(figsize=(6.33, 2.5))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0, 0])

# set up axes, labels
Slims = [0.002, 500.]
Rlims = [2.0, 500.]

ax.set_xlim(Rlims)
ax.set_xscale('log')
ax.set_xticks([10, 100])
ax.set_xticklabels(['10', '100'])
ax.set_xlabel('$r \;$ (au)')

ax.set_ylim(Slims)
ax.set_yscale('log')
ax.set_yticks([0.01, 0.1, 1, 10, 100])
ax.set_yticklabels(['0.01', '0.1', '1', '10', '100'])
ax.set_ylabel('$\Sigma \;$ (g cm$^{-2}$)')


### Set some constants
sigSB, G = 5.67051e-5, 6.67259e-8
au, Lsun, Msun = 1.4960e13, 3.826e33, 1.989e33
mu, mH, kB = 2.37, 1.673534e-24, 1.380658e-16


### Set the TW Hya properties
logL = -0.47
logM = -0.09


### Radius grid
rau = np.logspace(-1, 3.5, 1024)


### Toomre unstable zone (Q < 1) 
Tmid = (0.02 * 10.**(logL) * Lsun / (8.*np.pi*(rau * au)**2 * sigSB))**0.25
cs = np.sqrt(kB * Tmid / (mu*mH))
om = np.sqrt(G * 10.**logM * Msun / (rau*au)**3)
sig_GI = cs * om / (np.pi * G)
ax.fill_between(rau, sig_GI, Slims[1], interpolate=True, color='y', alpha=0.2)

### Hayashi 1981 MMSN
sigs_mmsn = 7.1 * rau**(-1.5)
sigs_mmsn[(rau >= 2.7)] = 30 * rau[(rau >= 2.7)]**(-1.5)
sigg_mmsn = 1700. * rau**(-1.5)
zone = (rau >= 0.35) & (rau <= 36.)
ax.plot(rau[zone], sigg_mmsn[zone], '--', color='gray', linewidth=3, alpha=0.5)
ax.plot(rau[zone], sigs_mmsn[zone], '--m', linewidth=3, alpha=0.5)


### Andrews 2012 CO (model sA) and dust (model pC)
sigd10, gam, rc, zeta = 0.29, 1.0, 35., 0.014
sigc = (10./rc)**gam * np.exp((10./rc))
sig_A12_g = sigc * (rau/rc)**(-gam) * np.exp(-(rau/rc)**(2.-gam))
ax.plot(rau, sig_A12_g/zeta, '-', color='gray', linewidth=2, alpha=1.0)

sigd10, gam, rc = 0.39, 0.75, 60.
sig_A12_d = sigd10 * (rau/10.)**(-gam) 
sig_A12_d[rau > rc] = 1e-10
ax.plot(rau, sig_A12_d, '-m', linewidth=2, alpha=1.0)


### Flaherty 2018 turbulence (variable \gamma)
Mg, gam, Rc = 0.05, 0.94, 10.**1.67
sig_F18_g = (2.-gam)*(Mg*Msun/(2.*np.pi*(Rc*au)**2)) * \
            (rau/Rc)**(-gam) * np.exp(-(rau/Rc)**(2.-gam))
ax.plot(rau, sig_F18_g, '-', color='gray', linewidth=2, alpha=1.0)


### Zhang 2017 13C18O
sig_Z17_g = 13. * (rau/20.5)**(-0.9)
zone = (rau >= 5.) & (rau <= 21.)
ax.plot(rau[zone], sig_Z17_g[zone], '-', color='gray', linewidth=2, alpha=1.0)


### van Boekel 17 models
rcm, sigg, sigs, sigl = np.loadtxt('data/twhya_vb17_struct.txt').T
gint = interp1d(rcm/au, sigg, fill_value='extrapolate')
sig_vB17_g = gint(rau)
dint = interp1d(rcm/au, sigs+sigl, fill_value='extrapolate')
sig_vB17_d = dint(rau)
ax.plot(rau, sig_vB17_g, '-', color='gray', linewidth=2, alpha=1.0)
ax.plot(rau, sig_vB17_d, '-m', linewidth=2, alpha=1.0)


### Hogerheijde 16 models
sig10, Rin, Rout, p0, p1, Rbr = 0.39, 4.07, 200., 0.53, 8.0, 47.1
sigbr = sig10 * (Rbr/10.)**(-p0)
sig_in = sigbr * (rau/Rbr)**(-p0)
sig_out = sigbr * (rau/Rbr)**(-p1)
sig_H16_d = sig_in
sig_H16_d[rau < Rin] = 1e-10
sig_H16_d[(rau >= Rbr) & (rau <= Rout)] = sig_out[(rau >= Rbr) & (rau <= Rout)]
sig_H16_d[rau > Rout] = 1e-10
ax.plot(rau, sig_H16_d, '-m', linewidth=2, alpha=1.0)



# annotations
ax.text(300., 90., 'unstable ($Q < 1$)', rotation=-22, fontsize=11, 
        color='goldenrod', horizontalalignment='right')


fig.subplots_adjust(left=0.2, right=0.8, bottom=0.165, top=0.98)
fig.savefig('tw_sigma.pdf')
fig.clf()
