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
Rlims = [2.0, 500.]
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
base = ( (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'T') &
         (db['FL_MULT'] != 'T') & (db['SED'] != 'III') &
         (db['SED'] != 'DEBRIS') & (db['SFR'] != 'Oph') &
         (db['FL_logMs'] == 0) )


### Set some constants
d_ref, nu_ref, alp = 150., 340., 2.3


### Lmm vs Mstar 
# calculate luminosities and upper limits
L7 = db['F_B7'] * (nu_ref / db['nu_B7'])**alp * (db['DPC'] / d_ref)**2
L6 = db['F_B6'] * (nu_ref / db['nu_B6'])**alp * (db['DPC'] / d_ref)**2
limL7 = db['LIM_B7'] * (nu_ref/db['nu_B7'])**alp * (db['DPC']/d_ref)**2
limL6 = db['LIM_B6'] * (nu_ref/db['nu_B6'])**alp * (db['DPC']/d_ref)**2
Mstar = 10.**(db['logMs'])

# targets with B7 detections
detB7 = ( (db['FL_B7'] == 0) & base )
L_detB7 = L7[detB7] 
M_detB7 = Mstar[detB7]

# targets with **only** B6 detections (i.e., no B7 or B7 limit)
detB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & base )
L_detB6 = L6[detB6]
M_detB6 = Mstar[detB6]

# targets with **only** limits or missing data
# (there should be **no entry** without a limit in at least B6 or B7)
lims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & base )
dlims = np.ma.column_stack( (limL7[lims], limL6[lims]) )
L_lims = np.ma.min(dlims, 1)
M_lims = Mstar[lims]

# combine all detections 
Lmm = np.ma.concatenate( (L_detB7, L_detB6) )
Ms = np.ma.concatenate( (M_detB7, M_detB6) )

ax0.errorbar(M_lims, L_lims, yerr=0.35*L_lims, uplims=True, marker='None', 
             capsize=1.5, alpha=0.5, linestyle='None')
ax0.plot(Ms, Lmm, 'oC0', markersize=2)



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.19, top=0.98)
fig.savefig('mstar.pdf')
fig.clf()
