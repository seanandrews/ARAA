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
fig = plt.figure(figsize=(6.33, 2.5))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0, 0])

# set up axes, labels
plims = [0., 1.]
Mlims = [0.0125, 8000.]		# earth masses

ax.set_xlim(Mlims)
ax.set_xscale('log')
ax.set_xticks([0.1, 1, 10, 100, 1000])
ax.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax.set_xlabel('$M \;$ (M$_{\\boldsymbol{\oplus}}$)')

ax.set_ylim(plims)
ax.yaxis.get_ticklocs(minor=True)
ax.minorticks_on()
ax.set_ylabel('$p$ ($\ge M$)')

# show the spectral index test?
show_test = False



### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)
print(np.unique(db['SFR']))

# baseline selections
base = ( (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'HJ') & 
         (db['FL_MULT'] != 'HJB') & (db['SFR'] != 'sOri') &  
         (db['SED'] != 'I') & (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') ) 


### Set some constants
d_ref, nu_ref = 150., 340.
h, c, k = 6.626e-27, 2.9979e10, 1.381e-16
pc, mearth, mjup = 3.0857e18, 5.974e27, 1.898e30
kappa, Tdust, alp = 3.5, 20., 2.2


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
print(len(Ms))

# calculate the combined CDF
Msolids, pMsolids, epMsolids, mukm = km_estimator(Ms, flags)



### Gas mass distribution
dets_Mg = ( (db['FL_Mgas'] == 0) & base )
lims_Mg = ( (db['FL_Mgas'] == 1) & base )
Mg_dets = db['Mgas'][dets_Mg] * mjup / mearth
Mg_lims = db['Mgas'][lims_Mg] * mjup / mearth
Mg = np.ma.concatenate((Mg_dets, Mg_lims))
flag_dets = np.zeros_like(Mg_dets)
flag_lims = np.ones_like(Mg_lims)
flaggs = np.ma.concatenate((flag_dets, flag_lims))
flagg = (flaggs == 1)

# cumulative distribution
Mgas, pMgas, epMgas, mukm = km_estimator(Mg, flagg)


### Solids for these gas measurements
# targets with B7 detections
DdetB7 = ( (db['FL_B7'] == 0) & dets_Mg )
DM_detB7 = L7[DdetB7] * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
Dd_detB7 = (db['FL_B7'][DdetB7] == 1)

LdetB7 = ( (db['FL_B7'] == 0) & lims_Mg )
LM_detB7 = L7[LdetB7] * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
Ld_detB7 = (db['FL_B7'][LdetB7] == 1)

# targets with **only** B6 detections (i.e., no B7 or B7 limit)
DdetB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & dets_Mg )
DM_detB6 = L6[DdetB6] * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
Dd_detB6 = (db['FL_B7'][DdetB6] == 0)

LdetB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & lims_Mg )
LM_detB6 = L6[LdetB6] * (d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
Ld_detB6 = (db['FL_B7'][LdetB6] == 0)

# targets with **only** limits or missing data
# (there should be **no entry** without a limit in at least B6 or B7)
Dlims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & dets_Mg )
Ddlims = np.ma.column_stack( (limL7[Dlims], limL6[Dlims]) )
DM_lims = np.ma.min(Ddlims, 1)*(d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
Dd_lims = np.ones(len(DM_lims), dtype=bool)

Llims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & lims_Mg )
Ldlims = np.ma.column_stack( (limL7[Llims], limL6[Llims]) )
LM_lims = np.ma.min(Ldlims, 1)*(d_ref * pc)**2 * 1e-23 / (kappa * Bnu) / mearth
Ld_lims = np.ones(len(LM_lims), dtype=bool)

### Solid mass distribution
# combine all sub-samples
Msg = np.ma.concatenate( (DM_detB7, LM_detB7, DM_detB6, LM_detB6, DM_lims, LM_lims) )
flags = np.ma.concatenate( (Dd_detB7, Ld_detB7, Dd_detB6, Ld_detB6, Dd_lims, Ld_lims) )

Mgs, pMgs, epMgs, mukm = km_estimator(Msg, flagg)


### Plot the distributions 
ax.fill_between(Mgs, pMgs+epMgs, pMgs-epMgs,
                facecolor='orange', alpha=0.3, step='post')
ax.plot(Mgs, pMgs, 'orange', drawstyle='steps-post', alpha=0.6, linewidth=3)

ax.fill_between(Mgas, pMgas+epMgas, pMgas-epMgas,
                facecolor='gray', alpha=0.3, step='post')
ax.plot(Mgas, pMgas, 'gray', drawstyle='steps-post', linewidth=3)

ax.fill_between(Msolids, pMsolids+epMsolids, pMsolids-epMsolids, 
                facecolor='m', alpha=0.3, step='post')
ax.plot(Msolids, pMsolids, 'm', drawstyle='steps-post', linewidth=3)


### Annotations
ax.text(0.07, 0.78, '$\\boldsymbol{M_{\\rm s}}$', color='m', fontsize=13,
        ha='right')
ax.text(0.13, 0.695, '(all disks)', color='m', fontsize=10, 
        horizontalalignment='right', alpha=0.8)
ax.text(0.6, 0.37, '$\\boldsymbol{M_{\\rm s}}$', color='orange', fontsize=13, 
        ha='right')
ax.text(0.63, 0.285, '($M_{\\rm g}$ sample)', color='orange', fontsize=10, 
        horizontalalignment='right', alpha=0.8)
ax.text(10, 0.78, '$\\boldsymbol{M_{\\rm g}}$', color='gray', fontsize=13)

ax.plot([35.,35.],[0., 0.45], ':', color='m', linewidth=2, zorder=0)
ax.plot([3000,3000],[0., 0.45], ':', color='gray', linewidth=2, zorder=0)


fig.subplots_adjust(left=0.2, right=0.8, bottom=0.165, top=0.975)
fig.savefig('mdisk_dist.pdf')
fig.clf()
