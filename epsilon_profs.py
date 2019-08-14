import numpy as np
import os
import sys
from astropy.io import ascii
import cloudpickle as cp
from scipy.interpolate import interp1d
import scipy.integrate as sci

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# set up plot
fig = plt.figure(figsize=(6.33, 2.25))
gs = gridspec.GridSpec(1, 2)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])

# set up axes, labels
alims = [0., 200.]
elims = [1.7, 4.7]
Rlims = [0., 220.]
wlims = [0.5, 20.]

# Panel (a) setups  [epsilon(r) profiles]
ax0.set_xlim(alims)
#ax0.set_xscale('log')
#ax0.set_xticks([0.1, 1, 10, 100, 1000])
#ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$r \;$ (au)')

ax0.set_ylim(elims)
ax0.set_ylabel('$\\varepsilon$')

# Panel (b) setups  [Rmm versus wavelength]
ax1.set_xlim(wlims)
ax1.set_xscale('log')
ax1.set_xticks([1, 10])
ax1.set_xticklabels(['1', '10'])
ax1.set_xlabel('$\\lambda \;$ (mm)')

ax1.set_ylim(Rlims)
#ax1.set_yscale('log')
#ax1.set_yticks([10, 100])
#ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('$R_{\\rm mm} \;$ (au)')


### Set some constants
cc = 2.9979e10


### Load the data from Tazzari+ 2016
with open("data/tazzari_profiles.dat", 'rb') as f:
    data = cp.load(f, encoding='latin1')


### index profiles
name = ['DRTau', 'FTTau', 'AS209']
col = ['C0', 'C1', 'C2']
for i in range(len(name)):
    rau, wl = data[name[i]]['gridrad'], data[name[i]]['wle']
    a, b = 0, len(wl)-1
    Ia  = data[name[i]]['intensity'][a][:,1]
    eIa = 0.5*(data[name[i]]['intensity'][a][:,2] - \
               data[name[i]]['intensity'][a][:,0])
    Ib  = data[name[i]]['intensity'][b][:,1]
    eIb = 0.5*(data[name[i]]['intensity'][b][:,2] - \
               data[name[i]]['intensity'][b][:,0])

    eps  = np.log(Ia/Ib) / np.log(wl[b]/wl[a])
    eeps = np.sqrt( (1./(Ia*np.log(wl[b]/wl[a])))**2 * eIa**2 + \
                    (1./(Ib*np.log(wl[b]/wl[a])))**2 * eIb**2 )
    ax0.fill_between(rau, eps+eeps, eps-eeps, facecolor=col[i], alpha=0.5)
    ax0.plot(rau, eps, '-'+col[i])


### calculate effective radii
    reffs = np.zeros(len(wl)) 
    ereffs_lo, ereffs_hi = np.zeros(len(wl)), np.zeros(len(wl))
    for j in range(len(reffs)):
        Fcum = sci.cumtrapz(2.*np.pi*data[name[i]]['intensity'][j][:,1]*rau,
                            rau, initial=0.)
        fint = interp1d(Fcum / Fcum[-1], rau)
        reffs[j] = fint(0.90)

        Fcum = sci.cumtrapz(2.*np.pi*data[name[i]]['intensity'][j][:,0]*rau,
                            rau, initial=0.)
        lint = interp1d(Fcum / Fcum[-1], rau)
        ereffs_lo[j] = reffs[j] - lint(0.90)

        Fcum = sci.cumtrapz(2.*np.pi*data[name[i]]['intensity'][j][:,2]*rau,
                            rau, initial=0.)
        hint = interp1d(Fcum / Fcum[-1], rau)
        ereffs_hi[j] = hint(0.90) - reffs[j]

    ax1.errorbar(wl, reffs, yerr=[5.*ereffs_lo, 7.*ereffs_hi], fmt='o', color=col[i], 
                 markersize=4) 


fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('epsilon_profs.pdf')
fig.clf()
