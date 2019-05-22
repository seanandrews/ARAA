import numpy as np
import os
import sys
from astropy.io import ascii
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rc('font', size=9)


# set up plot
fig = plt.figure(figsize=(7., 3.1))
gs = gridspec.GridSpec(1, 1)
Slims = [0.05, 2000.]
Rlims = [2.0, 500.]

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Rlims)
ax0.set_xscale('log')
ax0.set_ylim(Slims)
ax0.set_yscale('log')

# constants
sigSB = 5.67051e-5
au = 1.4960e13
Lsun = 3.826e33
Msun = 1.989e33
mu = 2.37
mH = 1.673534e-24
kB = 1.380658e-16
G = 6.67259e-8

# TW Hya properties
logL = -0.47
logM = -0.09

# Toomre stability (Q = 1 disk) 
rau = np.logspace(-1, 3.5, 1024)
Tmid = (0.02 * 10.**(logL) * Lsun / (8.*np.pi*(rau * au)**2 * sigSB))**0.25
cs = np.sqrt(kB * Tmid / (mu*mH))
om = np.sqrt(G * 10.**logM * Msun / (rau*au)**3)
sig_GI = cs * om / (np.pi * G)
ax0.fill_between(rau, sig_GI, Slims[1], interpolate=True, 
                 color='k', alpha=0.2)

# Andrews 2012 CO (model sA) and dust (model pC) : color = C0
A12_col = 'C0'
sigd10, gam, rc, zeta = 0.29, 1.0, 35., 0.014
sigc = (10./rc)**gam * np.exp((10./rc))
sig_A12_g = sigc * (rau/rc)**(-gam) * np.exp(-(rau/rc)**(2.-gam))
ax0.plot(rau, sig_A12_g/zeta, A12_col)

sigd10, gam, rc = 0.39, 0.75, 60.
sig_A12_d = sigd10 * (rau/10.)**(-gam) 
sig_A12_d[rau > rc] = 1e-10
ax0.plot(rau, sig_A12_d, ':'+A12_col)


# Flaherty 2018 turbulence (variable \gamma): color = C1
F18_col = 'C1'
Mg, gam, Rc = 0.05, 0.94, 10.**1.67
sig_F18_g = (2.-gam)*(Mg*Msun/(2.*np.pi*(Rc*au)**2)) * \
            (rau/Rc)**(-gam) * np.exp(-(rau/Rc)**(2.-gam))
ax0.plot(rau, sig_F18_g, F18_col)

# Zhang 2017 13C18O: color = C2
Z17_col = 'C2'
sig_Z17_g = 13. * (rau/20.5)**(-0.9)
zone = (rau >= 5.) & (rau <= 21.)
ax0.plot(rau[zone], sig_Z17_g[zone], Z17_col)

# van Boekel 17 models: color = C3
vB17_col = 'C3'
rcm, sigg, sigs, sigl = np.loadtxt('data/twhya_vb17_struct.txt').T
gint = interp1d(rcm/au, sigg, fill_value='extrapolate')
sig_vB17_g = gint(rau)
dint = interp1d(rcm/au, sigs+sigl, fill_value='extrapolate')
sig_vB17_d = dint(rau)
ax0.plot(rau, sig_vB17_g, vB17_col)
ax0.plot(rau, sig_vB17_d, ':'+vB17_col)


# Hogerheijde 16 models: color = C4
H16_col = 'C4'
sig10, Rin, Rout, p0, p1, Rbr = 0.39, 4.07, 200., 0.53, 8.0, 47.1
sigbr = sig10 * (Rbr/10.)**(-p0)
sig_in = sigbr * (rau/Rbr)**(-p0)
sig_out = sigbr * (rau/Rbr)**(-p1)
sig_H16_d = sig_in
sig_H16_d[rau < Rin] = 1e-10
sig_H16_d[(rau >= Rbr) & (rau <= Rout)] = sig_out[(rau >= Rbr) & (rau <= Rout)]
sig_H16_d[rau > Rout] = 1e-10
ax0.plot(rau, sig_H16_d, ':'+H16_col)


#ax0.set_xticks([1, 10, 100, 1000])
#ax0.set_xticklabels(['1', '10', '100', '1000'])
#ax0.set_yticks([10, 100])
#ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('$\Sigma$  (g cm$^{-2}$)')
ax0.set_xlabel('$r$  (au)')


fig.subplots_adjust(left=0.2, right=0.8, bottom=0.12, top=0.99)
fig.savefig('tw_sigma.pdf')
fig.clf()

