import numpy as np
import os
import sys
from astropy.io import ascii
from km_estimator import km_estimator

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# set up plot
fig = plt.figure(figsize=(6.33, 2.6))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0, 0])

# set up axes, labels
plims = [0., 1.]
Mlims = [0.05, 20000.]		# earth masses

ax.set_xlim(Mlims)
ax.set_xscale('log')
ax.set_xticks([0.1, 1, 10, 100, 1000, 10000])
ax.set_xticklabels(['0.1', '1', '10', '100', '1000', '10$^4$'])

ax.set_ylim(plims)
ax.set_ylabel('$p$ ($\ge M$)')
ax.set_xlabel('$M \;$ (M$_{\\boldsymbol{\oplus}}$)')

# show the spectral index test?
show_test = False



### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# baseline selections
base = ( (db['FL_MULT'] != 'J') & (db['SED'] != 'III') & 
         (db['SED'] != 'DEBRIS') & (db['SFR'] != 'Oph') )


### Set some constants
d_ref, nu_ref = 150., 340.
h, c, k = 6.626e-27, 2.9979e10, 1.381e-16
pc, mearth, mjup = 3.0857e18, 5.974e27, 1.898e30
kappa, Tdust, alp = 3.5, 20., 2.3


### Calculate CDFs for a test of the B6-->B7 conversion
if (show_test == True):
    dualband = ( (db['FL_B7'] == 0) & (db['FL_B6'] == 0) & base )
    flags = (db['FL_B7'][dualband] == 1)
    L7 = 1e-3 * db['F_B7'][dualband] * (db['DPC'][dualband] / d_ref)**2
    nu7 = db['nu_B7'][dualband] * 1e9
    L6 = 1e-3 * db['F_B6'][dualband] * (db['DPC'][dualband] / d_ref)**2
    nu6 = db['nu_B6'][dualband] * 1e9

    Bnu = (2 * h * nu7**3 / c**2) / (np.exp(h * nu7 / (k * Tdust)) - 1)
    M7 = (d_ref * pc)**2 * 1e-23 * L7 / (kappa * Bnu) 
    M6 = (d_ref * pc)**2 * 1e-23 * L6 * (nu7 / nu6)**alp / (kappa * Bnu) 

    Ms7, pMs7, epMs7, mukm = km_estimator(M7 / mearth, flags)
    Ms6, pMs6, epMs6, mukm = km_estimator(M6 / mearth, flags)

    ax.fill_between(Ms7, pMs7+epMs7, pMs7-epMs7, 
                    facecolor='C0', alpha=0.5, step='post')
    ax.plot(Ms7, pMs7, 'C0', drawstyle='steps-post')

    ax.fill_between(Ms6, pMs6+epMs6, pMs6-epMs6,
                    facecolor='C1', alpha=0.5, step='post')
    ax.plot(Ms6, pMs6, 'C1', drawstyle='steps-post')


### Selection and combination of mm luminosities
# calculate luminosities and upper limits
L7 = 1e-3 * db['F_B7'] * (nu_ref / db['nu_B7'])**alp * (db['DPC'] / d_ref)**2
L6 = 1e-3 * db['F_B6'] * (nu_ref / db['nu_B6'])**alp * (db['DPC'] / d_ref)**2
limL7 = 1e-3 * db['LIM_B7'] * (nu_ref/db['nu_B7'])**alp * (db['DPC']/d_ref)**2
limL6 = 1e-3 * db['LIM_B6'] * (nu_ref/db['nu_B6'])**alp * (db['DPC']/d_ref)**2

# Planck function at the reference frequency
Bnu = (2*h*(nu_ref*1e9)**3 / c**2) / (np.exp(h*nu_ref*1e9 / (k*Tdust)) - 1)

# targets with B7 detections
detB7 = ( (db['FL_B7'] == 0) & base )
M_detB7 = L7[detB7] * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
d_detB7 = (db['FL_B7'][detB7] == 1)

# targets with **only** B6 detections (i.e., no B7 or B7 limit)
detB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & base )
M_detB6 = L6[detB6] * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
d_detB6 = (db['FL_B7'][detB6] == 0)

# targets with **only** limits or missing data
# (there should be **no entry** without a limit in at least B6 or B7)
lims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & base )
dlims = np.ma.column_stack( (limL7[lims], limL6[lims]) )
M_lims = np.ma.min(dlims, 1) * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
d_lims = np.ones(len(M_lims), dtype=bool)


### Solid mass distribution
# combine all sub-samples
Ms = np.ma.concatenate( (M_detB7, M_detB6, M_lims) )
flags = np.ma.concatenate( (d_detB7, d_detB6, d_lims) )

# calculate the combined CDF
Msolids, pMsolids, epMsolids, mukm = km_estimator(Ms, flags)



### Gas mass distribution
have_Mg = ( (db['FL_Mgas'] == 0) & base )
flagg = (db['FL_Mgas'][have_Mg] == 1)
Mg = db['Mgas'][have_Mg] * mjup / mearth 

# cumulative distribution
Mgas, pMgas, epMgas, mukm = km_estimator(Mg, flagg)



### Plot the distributions 
ax.fill_between(Mgas, pMgas+epMgas, pMgas-epMgas,
                facecolor='gray', alpha=0.3, step='post')
ax.plot(Mgas, pMgas, 'gray', drawstyle='steps-post', linewidth=3)

ax.fill_between(Msolids, pMsolids+epMsolids, pMsolids-epMsolids, 
                facecolor='m', alpha=0.3, step='post')
ax.plot(Msolids, pMsolids, 'm', drawstyle='steps-post', linewidth=3)


### Annotations



fig.subplots_adjust(left=0.2, right=0.8, bottom=0.16, top=0.98)
fig.savefig('mdisk_dist.pdf')
fig.clf()
