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
fig = plt.figure(figsize=(6.33, 2.2))
gs = gridspec.GridSpec(1, 2)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])

# set up axes, labels
Llims  = [0.05, 5000.]
Rlims  = [2.0, 500.]
COlims = [5, 2000]

# Panel (a) setups  [size-luminosity relation]
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm mm} \;$ (mJy at 150 pc)')

ax0.set_ylim(Rlims)
ax0.set_yscale('log')
ax0.set_yticks([10, 100])
ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('$R_{\\rm mm} \;$ (au)')

# Panel (b) setups  [R_dust versus R_gas]
ax1.set_xlim(COlims)
ax1.set_xscale('log')
ax1.set_xticks([10, 100, 1000])
ax1.set_xticklabels(['10', '100', '1000'])
ax1.set_xlabel('$R_{\\rm CO} \;$ (au)')

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
         (db['SED'] != 'DEBRIS') & (db['SFR'] != 'Oph') )


### Set some constants
d_ref = 150.


### luminosities and sizes
det = ((db['FL_B7'] == 0) & (db['FL_R7'] == 0))
lim = ((db['FL_B7'] == 0) & (db['FL_R7'] == 1))

Ldet = db['F_B7'][det] * (db['DPC'][det] / d_ref)**2
Rdet = db['R7'][det] * db['DPC'][det]
Llim = db['F_B7'][lim] * (db['DPC'][lim] / d_ref)**2
Rlim = db['LIM_R7'][lim] * db['DPC'][lim]

ax0.errorbar(Llim, Rlim, yerr=0.25*Rlim, uplims=True, marker='None', capsize=2, 
             alpha=0.5, linestyle='None')
ax0.plot(Ldet, Rdet, 'oC0', markersize=4)


### continuum sizes and CO sizes
det = ((db['FL_RCO'] == 0) & (db['FL_R7'] == 0))
lim = ((db['FL_RCO'] == 0) & (db['FL_R7'] == 1))

RCO_det = db['R_CO'][det] * db['DPC'][det]
eRCO_det = db['eR_CO'][det] * db['DPC'][det]
Rmm_det = db['R7'][det] * db['DPC'][det]
#RCO_lim = db['R_CO'][det] * db['DPC'][det]
#Rmm_lim = db['R7'][det] * db['DPC'][det]

# simple scaling
mdl_x = np.logspace(0, 4, 128)

ax1.plot(mdl_x, mdl_x, '--', color='gray')
ax1.plot(3.*mdl_x, mdl_x, ':', color='gray')
#ax1.plot(Llim, Rlim, '11', color='C1', markersize=4.0)
ax1.errorbar(RCO_det, Rmm_det, xerr=eRCO_det, marker='o', color='C0', 
             markersize=4.0, linestyle='None')


fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('sizes.pdf')
fig.clf()
