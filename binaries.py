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
Llims  = [0.05, 5000.]
Rlims  = [1.25, 800.]
Dlims  = [0.2, 5000]

# Panel (a) setups  [luminosity versus projected separation]
ax0.set_ylim(Llims)
ax0.set_yscale('log')
ax0.set_yticks([0.1, 1, 10, 100, 1000])
ax0.set_yticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_ylabel('$\sum_i L_{{\\rm mm}, i} \;$ (mJy at 150 pc)', labelpad=3)

ax0.set_xlim(Dlims)
ax0.set_xscale('log')
ax0.set_xticks([1, 10, 100, 1000])
ax0.set_xticklabels(['1', '10', '100', '1000'])
ax0.set_xlabel('$\Delta \;$ (au)')

# Panel (b) setups  [R_dust versus R_gas]
ax1.set_ylim(Rlims)
ax1.set_yscale('log')
ax1.set_yticks([10, 100])
ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('$R_{\\rm mm} \;$ (au)')

ax1.set_xlim(Llims)
ax1.set_xscale('log')
ax1.set_xticks([0.1, 1, 10, 100, 1000])
ax1.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax1.set_xlabel('$L_{\\rm mm} \;$ (mJy at 150 pc)')


### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# baseline selections
baseJ  = (db['FL_MULT'] == 'J')
baseCB = (db['FL_MULT'] == 'CB')



### Set some constants
d_ref = 150.


### luminosities and projected separations
# simple model
#ml = np.logspace(-2, 4, 1024)
#mr = 7.*(ml/1.)**0.5
#ax0.plot(ml, mr, '--r')

# upper limits
lim = ((db['FL_B6'] == 0) & baseJ)
Llim = db['LIM_B6'][lim] * (db['DPC'][lim] / d_ref)**2
Dlim = db['SEP'][lim] * db['DPC'][lim]
ax0.errorbar(Dlim, Llim, yerr=0.25*Llim, uplims=True, marker='None', capsize=2,
             color='gray', alpha=0.65, linestyle='None')

lim = ((db['FL_B6'] == 1) & baseCB)
Llim = db['LIM_B6'][lim] * (db['DPC'][lim] / d_ref)**2
Dlim = db['SEP'][lim] * db['DPC'][lim]
print(len(Llim))
ax0.errorbar(Dlim, Llim, yerr=0.25*Llim, uplims=True, marker='None', capsize=2,
             color='orange', alpha=0.65, linestyle='None')


# detections (do this in B6?)
detJ = ((db['FL_B6'] == 0) & baseJ) 
LdetJ = db['F_B6'][detJ] * (db['DPC'][detJ] / d_ref)**2
LerrJ = np.sqrt( (db['eF_B6'][detJ]**2 + (0.1*db['F_B6'][detJ])**2) * \
                 (db['DPC'][detJ]/d_ref)**4 + \
                 ( db['F_B6'][detJ]*(2.*db['DPC'][detJ]/d_ref**2) * \
                   0.5 * (db['EDPC_H'][detJ]+db['EDPC_L'][detJ]) )**2 )
DdetJ = db['SEP'][detJ] * db['DPC'][detJ]
DerrJ = np.sqrt((db['DPC'][detJ]*0.1*db['SEP'][detJ])**2 + 
                (db['SEP'][detJ] * \
                 0.5*(db['EDPC_H'][detJ]+db['EDPC_L'][detJ]))**2)
ax0.errorbar(DdetJ, LdetJ, xerr=DerrJ, yerr=LerrJ, marker='o', 
             color='C0', markersize=3, linestyle='None', elinewidth=1.0, 
             alpha=0.75)

detC = ((db['FL_B6'] == 0) & baseCB)
LdetC = db['F_B6'][detC] * (db['DPC'][detC] / d_ref)**2
LerrC = np.sqrt( (db['eF_B6'][detC]**2 + (0.1*db['F_B6'][detC])**2) * \
                (db['DPC'][detC]/d_ref)**4 + \
                ( db['F_B6'][detC]*(2.*db['DPC'][detC]/d_ref**2) * \
                  0.5 * (db['EDPC_H'][detC]+db['EDPC_L'][detC]) )**2 )
print(len(LdetC))
DdetC = db['SEP'][detC] * db['DPC'][detC]
DerrC = np.sqrt((db['DPC'][detC]*0.1*db['SEP'][detC])**2 + 
                (db['SEP'][detC] * \
                 0.5*(db['EDPC_H'][detC]+db['EDPC_L'][detC]))**2)
ax0.errorbar(DdetC, LdetC, xerr=DerrC, yerr=LerrC, marker='o', 
             color='r', markersize=3, linestyle='None', elinewidth=1.0, 
             alpha=0.75)




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
             alpha=0.75)


ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('binaries.pdf')
fig.clf()
