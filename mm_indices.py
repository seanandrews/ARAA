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
elims = [1.5, 4.7]
Rlims = [0., 200.]
Llims = [0.05, 5000.]

# panel (a) setups  [alphamm versus Lmm]
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm mm} \;$ (mJy at 150 pc)')

ax0.set_ylim(elims)
ax0.set_ylabel('$\\alpha_{\\rm mm} \;$ (1--3 mm)')


# panel (b) setups  [epsilon(r) profiles]
ax1.set_xlim(Rlims)
ax1.set_xlabel('$r \;$ (au)')

ax1.set_ylim(elims)
ax1.set_ylabel('$\\varepsilon \;$ (1--9 mm)')


### Set some constants
cc = 2.9979e10


ind = '36'

if (ind == '67'): alp_lbl = '(B6 / B7)'
if (ind == '36'): alp_lbl = '(B3 / B6)'

# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# downselect targets with both a robust index and B6 flux density
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & 
         (db['FL_MULT'] != 'J') & 
         (db['SED'] != 'III') & 
         (db['SED'] != 'DEBRIS') )

d_ref, nu_ref, alp = 150., 225., 2.3

# extract sub-sample
sub = ( (db['FL_A'+ind] == 0) & (db['FL_B6'] == 0) & base )
L6  = db['F_B6'][sub] * (db['DPC'][sub] / d_ref)**2 * \
      (nu_ref / db['nu_B6'][sub])**alp
eL6 = np.sqrt( (nu_ref/db['nu_B6'][sub])**(2.*alp) * \
               (db['eF_B6'][sub]**2 + (0.1*db['F_B6'][sub])**2) * \
               (db['DPC'][sub]/d_ref)**4 + \
               ( ((nu_ref / db['nu_B6'][sub])**alp * \
                  db['F_B6'][sub]*(2.*db['DPC'][sub]/d_ref**2) * \
                 0.5*(db['EDPC_H'][sub]+db['EDPC_L'][sub]) )**2 ) )
amm = db['A'+ind][sub]
eamm_hi = db['eA'+ind+'_hi'][sub]
eamm_lo = db['eA'+ind+'_lo'][sub]

ax0.errorbar(L6, amm, xerr=eL6, yerr=[eamm_lo, eamm_hi], marker='o', 
             color='C0', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)


### Load the data from Tazzari+ 2016
with open("data/tazzari_profiles.dat", 'rb') as f:
    data = cp.load(f, encoding='latin1')

### index profiles
name = ['DRTau', 'FTTau', 'AS209']
col = ['C3', 'C1', 'C2']
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
    ax1.fill_between(rau, eps+eeps, eps-eeps, facecolor=col[i], alpha=0.5)
    ax1.plot(rau, eps, '-'+col[i])


fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('mm_indices.pdf')
fig.clf()
