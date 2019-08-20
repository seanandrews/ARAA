import numpy as np
import os
import sys
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from do_regression import do_regression
plt.rc('font', size=9)


sfr  = 'Oph'
band = 7
eFsys = 0.1
do_reg = True

lbl = 'mstar_'+sfr+'_B'+str(band)+'lowcut'


# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# set up plots
fig = plt.figure(figsize=(6.33, 2.6))
gs = gridspec.GridSpec(1, 1)
Llims = [0.125, 8000.]
Mlims = [0.0125, 8.]
dref = 150.

ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(Mlims)
ax.set_xscale('log')
ax.set_ylim(Llims)
ax.set_yscale('log')
ax.set_yticks([1, 10, 100, 1000])
ax.set_yticklabels(['1', '10', '100', '1000'])
ax.set_xticks([0.1, 1])
ax.set_xticklabels(['0.1', '1'])
ax.set_xlabel('$M_\\ast$  ($M_\odot$)')
ax.set_ylabel('$L_{\\rm mm}$  (mJy at 140 pc)')


# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & \
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & \
         (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') & 
         (db['SFR'] == sfr) & (db['FL_logMs'] == 0) & 
         (db['logMs'] >= -1.0) & (db['logMs'] <= 0.3) )


# measurements available
ok = ((db['FL_B'+str(band)] == 0) & base)

L = db['F_B'+str(band)][ok] * (db['DPC'][ok] / dref)**2	
eF = np.sqrt(db['eF_B'+str(band)][ok]**2 + (eFsys * db['F_B'+str(band)][ok])**2)
eD = 0.5*(db['EDPC_H'][ok] + db['EDPC_L'][ok])
eL = np.sqrt((db['DPC'][ok]/dref)**4 * eF**2 + 4.*(L/db['DPC'][ok])**2 * eD**2)
logL = np.log10(L)
elogL = 0.4343 * eL / L

logM = db['logMs'][ok]
elogM = 0.5*(db['elogMs_H'][ok] + db['elogMs_L'][ok])


# upper limits
ok = ((db['FL_B'+str(band)] == 1) & base)

Ll = db['LIM_B'+str(band)][ok] * (db['DPC'][ok] / dref)**2 
logLl = np.log10(Ll)

logMl = db['logMs'][ok]
elogMl = 0.5*(db['elogMs_H'][ok] + db['elogMs_L'][ok])


# plot the data
ax.plot(10.**logMl, Ll, '.C2')
ax.errorbar(10.**logM, L, yerr=eL, fmt='.C0')


# do regression if requested
if (do_reg == True):
    flags_det = np.zeros(len(logL))
    flags_lim = np.ones(len(logLl))
    flags = np.ma.concatenate((flags_det, flags_lim))
    xx = np.ma.concatenate((logM, logMl))
    yy = np.ma.concatenate((logL, logLl))
    ex = np.ma.concatenate((elogM, elogMl))
    ey = np.ma.concatenate((elogL, np.zeros(len(logLl))))
    delta = (flags == 0)
    rposts = do_regression(xx, yy, ex, ey, delta, lbl)


# plot regression draws
xreg = np.linspace(np.log10(Mlims[0]), np.log10(Mlims[1]), num=100)
if (os.path.isfile('outfiles/'+lbl+'.regression.npz') == True):
    chain = (np.load('outfiles/'+lbl+'.regression.npz'))['chain']
    slopes  = chain['beta']
    intcpts = chain['alpha']
    scat    = np.sqrt(chain['sigsqr'])
    yreg_lo = np.zeros(len(xreg))
    yreg_hi = np.zeros(len(xreg))
    ysca_lo = np.zeros(len(xreg))
    ysca_hi = np.zeros(len(xreg))
    for i in range(len(xreg)):
        yreg_lo[i] = np.percentile(intcpts+slopes*xreg[i], 15.865)
        yreg_hi[i] = np.percentile(intcpts+slopes*xreg[i], 84.135)
        devs = np.random.normal(0, scat)
        ysca_lo[i] = np.percentile(intcpts+slopes*xreg[i]+devs, 15.865)
        ysca_hi[i] = np.percentile(intcpts+slopes*xreg[i]+devs, 84.135)
    ax.fill_between(10.**xreg, 10.**ysca_lo, 10.**ysca_hi, facecolor='r', 
                    alpha=0.2, interpolate=True)
    ax.fill_between(10.**xreg, 10.**yreg_lo, 10.**yreg_hi, facecolor='r', 
                    alpha=0.6, interpolate=True)


fig.subplots_adjust(wspace=0.31, hspace=0.0)
fig.subplots_adjust(left=0.2, right=0.8, bottom=0.15, top=0.98)
fig.savefig(lbl+'.pdf')
fig.clf()
