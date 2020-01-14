import numpy as np
import os
import sys
from astropy.io import ascii
from do_regression import do_regression

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
import matplotlib as mpl
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


do_reg = True


# set up plot
fig = plt.figure(figsize=(6.33, 2.25))
gs = gridspec.GridSpec(1, 2)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])

# set up axes, labels
Llims  = [0.05, 5000.]
Rlims  = [1.25, 800.]
Dlims  = [0.125, 8000]

# Panel (a) setups  [luminosity versus projected separation]
ax0.set_ylim(Llims)
ax0.set_yscale('log')
ax0.set_yticks([0.1, 1, 10, 100, 1000])
ax0.set_yticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_ylabel('$\sum L_{{\\rm 1.3 \, mm}} \;$ (mJy at 150 pc)', labelpad=3)

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
ax1.set_ylabel('$R_{\\rm 1.3 \, mm} \;$ (au, at 0.68$L$)')

ax1.set_xlim(Llims)
ax1.set_xscale('log')
ax1.set_xticks([0.1, 1, 10, 100, 1000])
ax1.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax1.set_xlabel('$L_{\\rm 1.3 \, mm} \;$ (mJy at 150 pc)')
locmin = mpl.ticker.LogLocator(base=10.,
                               subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                               numticks=5)
ax1.xaxis.set_minor_locator(locmin)
ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())



### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp2.csv')
db = ascii.read('temp2.csv', format='csv', fast_reader=True)

# baseline selections
baseJ   = ((db['FL_MULT'] != 'CB') & (db['FL_MULT'] != 'HCB') & 
           (db['FL_MULT'] != 'HC') & (db['FL_MULT'] != 'B') & 
           (db['FL_MULT'] != 'W')  & (db['FL_MULT'] != 'S') & 
           (db['FL_MULT'] != 'U') & (db['SED'] != 'III') & 
           (db['SED'] != 'I') & (db['SED'] != 'DEBRIS') ) 

baseCB  = ((db['FL_MULT'] != 'J')  & (db['FL_MULT'] != 'WJ') & 
           (db['FL_MULT'] != 'HJ') & (db['FL_MULT'] != 'HJB') & 
           (db['FL_MULT'] != 'HC') & (db['FL_MULT'] != 'B') &
           (db['FL_MULT'] != 'W')  & (db['FL_MULT'] != 'S') &
           (db['FL_MULT'] != 'U') & (db['SED'] != 'III') & 
           (db['SED'] != 'I') & (db['SED'] != 'DEBRIS') )



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
             color='gray', alpha=0.35, linestyle='None')

lim = ((db['FL_B6'] == 1) & baseCB)
Llim = db['LIM_B6'][lim] * (db['DPC'][lim] / d_ref)**2
Dlim = db['SEP'][lim] * db['DPC'][lim]
ax0.errorbar(Dlim, Llim, yerr=0.25*Llim, uplims=True, marker='None', capsize=2,
             color='C2', alpha=0.35, linestyle='None')


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
             color='m', markersize=3, linestyle='None', elinewidth=1.0, 
             alpha=0.65)

if (do_reg == True):
    flags = np.zeros_like(DdetJ)
    xx = np.log10(DdetJ)
    yy = np.log10(LdetJ)
    ex = 0.4343*DerrJ/xx
    ey = 0.4343*LerrJ/yy
    delta = (flags == 0)
    rposts = do_regression(xx, yy, ex, ey, delta, 'binaries')

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
             color='C2', markersize=3, linestyle='None', elinewidth=1.0, 
             alpha=0.65)


ax0.text(0.5, 0.11, 'circumbinary', ha='left', fontsize=8, color='C2')
ax0.text(100, 0.20, 'binary pairs', ha='left', fontsize=8, color='m')






### luminosities and sizes
# simple models
mx = np.logspace(-3, 4, 1024)
my = 10.**(2.10)*(mx/1e3)**0.49
nrat = (225./340.)
ax1.plot(mx*nrat**2.3, my*nrat**0.3, '--r', alpha=0.6)

# data
det = ((db['FL_B6'] == 0) & (db['FL_R6'] == 0))
Ldet = db['F_B6'][det] * (db['DPC'][det] / d_ref)**2
Lerr = np.sqrt( (db['eF_B6'][det]**2 + (0.1*db['F_B6'][det])**2) * \
                (db['DPC'][det]/d_ref)**4 + \
                ( db['F_B6'][det]*(2.*db['DPC'][det]/d_ref**2) * \
                  0.5 * (db['EDPC_H'][det]+db['EDPC_L'][det]) )**2 )
Rdet = 10.**db['R6'][det]
Rerr_hi = 10.**(db['eR6_hi'][det]+db['R6'][det]) - Rdet
Rerr_lo = Rdet - 10.**(db['R6'][det]-db['eR6_lo'][det])


# singles
mult = ((db['FL_MULT'][det] != 'B') & (db['FL_MULT'][det] != 'W'))
ax1.errorbar(Ldet[mult], Rdet[mult], xerr=Lerr[mult], 
             yerr=[Rerr_lo[mult], Rerr_hi[mult]], marker='o',
             color='C0', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

# binaries
mult = (db['FL_MULT'][det] == 'B')
ax1.errorbar(Ldet[mult], Rdet[mult], xerr=Lerr[mult],
             yerr=[Rerr_lo[mult], Rerr_hi[mult]], marker='o',
             color='m', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

mult = (db['FL_MULT'][det] == 'W')
ax1.errorbar(Ldet[mult], Rdet[mult], xerr=Lerr[mult],
             yerr=[Rerr_lo[mult], Rerr_hi[mult]], marker='o',
             color='m', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

mult = (db['FL_MULT'][det] == 'HJB')
ax1.errorbar(Ldet[mult], Rdet[mult], xerr=Lerr[mult],
             yerr=[Rerr_lo[mult], Rerr_hi[mult]], marker='o',
             color='m', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

mult = (db['FL_MULT'][det] == 'HC')
ax1.errorbar(Ldet[mult], Rdet[mult], xerr=Lerr[mult],
             yerr=[Rerr_lo[mult], Rerr_hi[mult]], marker='o',
             color='m', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

mult = (db['FL_MULT'][det] == 'HCB')
ax1.errorbar(Ldet[mult], Rdet[mult], xerr=Lerr[mult],
             yerr=[Rerr_lo[mult], Rerr_hi[mult]], marker='o',
             color='m', markersize=3, linestyle='None', elinewidth=1.0,
             alpha=0.65)

ax1.text(1.5, 2.5, 'binaries', ha='left', fontsize=8, color='m')
ax1.text(100, 280, 'singles', ha='left', fontsize=8, color='C0')





ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('binaries.pdf')
fig.clf()
