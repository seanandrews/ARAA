import numpy as np
import os
import sys
from astropy.io import ascii

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
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
Llims  = [0.05, 5000.]
Rlims  = [1.25, 800.]
COlims = [4.25, 1500]

# Panel (a) setups  [size-luminosity relation]
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm 0.9 \, mm} \;$ (mJy at 150 pc)')
locmin = mpl.ticker.LogLocator(base=10.,
                               subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                               numticks=5)
ax0.xaxis.set_minor_locator(locmin)
ax0.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax0.set_ylim(Rlims)
ax0.set_yscale('log')
ax0.set_yticks([10, 100])
ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('$R_{\\rm 0.9 \, mm} \;$ (au, to 0.9$L$)')

# Panel (b) setups  [R_dust versus R_gas]
ax1.set_xlim(COlims)
ax1.set_xscale('log')
ax1.set_xticks([10, 100, 1000])
ax1.set_xticklabels(['10', '100', '1000'])
ax1.set_xlabel('$R_{\\rm CO} \;$ (au, to 0.9$L$)')

ax1.set_ylim(Rlims)
ax1.set_yscale('log')
ax1.set_yticks([10, 100])
ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('$R_{\\rm 0.9 \, mm} \;$ (au, to 0.9$L$)')


### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# baseline selections
base = ( (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'B') & 
         (db['FL_MULT'] != 'T') & (db['SED'] != 'III') &
         (db['SED'] != 'DEBRIS') )


### Set some constants
d_ref = 150.


### luminosities and sizes
# simple model
ml = np.logspace(-2, 4, 1024)
mr = 7.*(ml/1.)**0.5
ax0.plot(ml, mr, '--r')

# size upper limits
lim = ((db['FL_B7'] == 0) & (db['FL_R7'] == 1))
Llim = db['F_B7'][lim] * (db['DPC'][lim] / d_ref)**2
Rlim = 10.**db['LIM_R7'][lim]
ax0.errorbar(Llim, Rlim, yerr=0.25*Rlim, uplims=True, marker='None', capsize=2,
             color='gray', alpha=0.65, linestyle='None')

# detections
det = ((db['FL_B7'] == 0) & (db['FL_R7'] == 0))
Ldet = db['F_B7'][det] * (db['DPC'][det] / d_ref)**2
Lerr = np.sqrt( (db['eF_B7'][det]**2 + (0.1*db['F_B7'][det])**2) * \
                (db['DPC'][det]/d_ref)**4 + \
                ( db['F_B7'][det]*(2.*db['DPC'][det]/d_ref**2) * \
                  0.5 * (db['EDPC_H'][det]+db['EDPC_L'][det]) )**2 )
Rdet = 10.**db['R7'][det] 
Rerr_hi = 10.**(db['eR7_hi'][det]+db['R7'][det]) - Rdet
Rerr_lo = Rdet - 10.**(db['R7'][det]-db['eR7_lo'][det])
ax0.errorbar(Ldet, Rdet, xerr=Lerr, yerr=[Rerr_lo, Rerr_hi], marker='o', 
             color='C0', markersize=3, linestyle='None', elinewidth=1.0, 
             alpha=0.65)



### continuum sizes and CO sizes
# simple models
mx = np.logspace(0, 4, 1024)
ax1.plot(mx, mx, ':', color='gray')
ax1.plot(2.5*mx, mx, '--r')

# data
det = ((db['FL_RCO'] == 0) & (db['FL_R7'] == 0))
RCO  = db['R_CO'][det] * db['DPC'][det]
eRCO = np.sqrt((db['eR_CO'][det]*db['DPC'][det])**2 + \
               (db['R_CO'][det]*(0.5*(db['EDPC_H'][det]+db['EDPC_L'][det])))**2)
Rmm = 10.**db['R7'][det] 
Rmm_hi = 10.**(db['eR7_hi'][det]+db['R7'][det]) - Rmm
Rmm_lo = Rmm - 10.**(db['R7'][det]-db['eR7_lo'][det])
ax1.errorbar(RCO, Rmm, xerr=eRCO, yerr=[Rmm_lo, Rmm_hi], marker='o', 
             color='C0', markersize=3, linestyle='None', elinewidth=1.0, 
             alpha=0.65)


ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('sizes.pdf')
fig.clf()
