import numpy as np
import os
import sys
from astropy.io import ascii

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
Llims = [0.05, 5000.]
Rlims = [1.25, 800.]
Mlims = [0.0125, 8.]

# Panel (a) setups  [mass-luminosity relation]
ax0.set_xlim(Mlims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1])
ax0.set_xticklabels(['0.1', '1'])
ax0.set_xlabel('$M_{\\ast} \;$ (M$_\odot$)')

ax0.set_ylim(Llims)
ax0.set_yscale('log')
ax0.set_yticks([0.1, 1, 10, 100, 1000])
ax0.set_yticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_ylabel('$L_{\\rm mm} \;$ (mJy at 150 pc)')

# Panel (b) setups  [mass-size relation]
ax1.set_xlim(Mlims)
ax1.set_xscale('log')
ax1.set_xticks([0.1, 1])
ax1.set_xticklabels(['0.1', '1'])
ax1.set_xlabel('$M_{\\ast} \;$ (M$_\odot$)')

ax1.set_ylim(Rlims)
ax1.set_yscale('log')
ax1.set_yticks([10, 100])
ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('$R_{\\rm mm} \;$ (au)')


### Load the database

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# baseline selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'HJB') &
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') &
         (db['FL_MULT'] != 'WJ') & (db['FL_MULT'] != 'HJ') &
         (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') &
         (db['FL_logMs'] == 0) )


### Set some constants
d_ref, nu_ref, alp = 150., 340., 2.3


# calculate luminosities, upper limits, masses, and sizes + uncertainties
L7 = db['F_B7'] * (nu_ref / db['nu_B7'])**alp * (db['DPC'] / d_ref)**2
L6 = db['F_B6'] * (nu_ref / db['nu_B6'])**alp * (db['DPC'] / d_ref)**2
eL7 = np.sqrt( (nu_ref/db['nu_B7'])**(2.*alp) * \
               (db['eF_B7']**2 + (0.1*db['F_B7'])**2) * \
               (db['DPC']/d_ref)**4 + \
               ( ((nu_ref / db['nu_B7'])**alp * \
                  db['F_B7']*(2.*db['DPC']/d_ref**2) * \
                 0.5*(db['EDPC_H']+db['EDPC_L']) )**2 ) )
eL6 = np.sqrt( (nu_ref/db['nu_B6'])**(2.*alp) * \
               (db['eF_B6']**2 + (0.1*db['F_B6'])**2) * \
               (db['DPC']/d_ref)**4 + \
               ( ((nu_ref / db['nu_B6'])**alp * \
                  db['F_B6']*(2.*db['DPC']/d_ref**2) * \
                 0.5*(db['EDPC_H']+db['EDPC_L']) )**2 ) )
limL7 = db['LIM_B7'] * (nu_ref/db['nu_B7'])**alp * (db['DPC']/d_ref)**2
limL6 = db['LIM_B6'] * (nu_ref/db['nu_B6'])**alp * (db['DPC']/d_ref)**2
Mstar = 10.**(db['logMs'])
eM_hi = 10.**(db['elogMs_H']+db['logMs']) - Mstar
eM_lo = Mstar - 10.**(db['logMs']-db['elogMs_L'])
Rmm = 10.**db['R7']
lim_Rmm = 10.**db['LIM_R7']
eRmm_hi = 10.**(db['eR7_hi']+db['R7']) - Rmm
eRmm_lo = Rmm - 10.**(db['R7']-db['eR7_lo'])

 

### Lmm versus Mstar
# targets with B7 detections
detB7 = ( (db['FL_B7'] == 0) & base )
L_detB7 = L7[detB7] 
eL_detB7 = eL7[detB7]
M_detB7 = Mstar[detB7]
eMhi_detB7 = eM_hi[detB7]
eMlo_detB7 = eM_lo[detB7]


# targets with **only** B6 detections (i.e., no B7 or B7 limit)
detB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & base )
L_detB6 = L6[detB6]
eL_detB6 = eL6[detB6]
M_detB6 = Mstar[detB6]
eMhi_detB6 = eM_hi[detB6]
eMlo_detB6 = eM_lo[detB6]


# targets with **only** limits or missing data
# (there should be **no entry** without a limit in at least B6 or B7)
lims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & base )
dlims = np.ma.column_stack( (limL7[lims], limL6[lims]) )
L_lims = np.ma.min(dlims, 1)
M_lims = Mstar[lims]
eMhi_lims = eM_hi[lims]
eMlo_lims = eM_lo[lims]


# combine all detections 
Lmm = np.ma.concatenate( (L_detB7, L_detB6) )
eLmm = np.ma.concatenate( (eL_detB7, eL_detB6) )
Ms = np.ma.concatenate( (M_detB7, M_detB6) )
eMs_hi = np.ma.concatenate( (eMhi_detB7, eMhi_detB6) )
eMs_lo = np.ma.concatenate( (eMlo_detB7, eMlo_detB6) )


# plot
ax0.errorbar(M_lims, L_lims, yerr=0.35*L_lims, uplims=True, marker='None', 
             color='gray', alpha=0.5, capsize=1.5, linestyle='None')
ax0.errorbar(M_lims, L_lims, xerr=[eMlo_lims, eMhi_lims], yerr=0, 
             marker='None', color='gray', alpha=0.4, linestyle='None')
ax0.errorbar(Ms, Lmm, xerr=[eMs_lo, eMs_hi], yerr=eLmm, marker='o', 
             color='C0', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

mx = np.logspace(-3, 2, 1024)
my = 100. * mx**1.7
ax0.plot(mx, my, '--r')



### Rmm versus Mstar

# selections
lim = ( (db['FL_B7'] == 0) & (db['FL_R7'] == 1) & base )
det = ( (db['FL_B7'] == 0) & (db['FL_R7'] == 0) & base )

# limits
Rlim = Rmm[lim]
Mlim = Mstar[lim]
eMhi_lim = eM_hi[lim]
eMlo_lim = eM_lo[lim]

# detections
Rdet = Rmm[det]
eRdet_hi = eRmm_hi[det]
eRdet_lo = eRmm_lo[det]
Mdet = Mstar[det]
eMdet_hi = eM_hi[det]
eMdet_lo = eM_lo[det]


# plot
ax1.errorbar(Mlim, Rlim, yerr=0.25*Rlim, uplims=True, marker='None', 
             capsize=1.5, color='gray', alpha=0.5, linestyle='None')
ax1.errorbar(Mlim, Rlim, xerr=[eMlo_lim, eMhi_lim], yerr=0,
             marker='None', color='gray', alpha=0.4, linestyle='None')
ax1.errorbar(Mdet, Rdet, xerr=[eMdet_lo, eMdet_hi], yerr=[eRdet_lo, eRdet_hi], 
             marker='o', color='C0', markersize=3, linestyle='None', 
             elinewidth=1.0, alpha=0.65)
my = 120. * mx**0.9
ax1.plot(mx, my, '--r')

ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.19, top=0.98)
fig.savefig('mstar.pdf')
fig.clf()
