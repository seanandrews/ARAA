import numpy as np
import os
import sys
from astropy.io import ascii
from scipy.interpolate import interp1d
from km_estimator import km_estimator

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# set up plots
fig = plt.figure(figsize=(6.33, 2.25))
gs = gridspec.GridSpec(1, 2)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])

# set up axes, labels
Llims = [0.02, 5000.]
Rlims = [1.25, 800]
Tlims = [0.2, 50.]

ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm mm} \;$  (mJy at 150 pc)')
ax0.set_ylim([0., 1.01])
ax0.set_ylabel('$p(> L_{\\rm mm})$')

ax1.set_xlim(Tlims)
ax1.set_xscale('log')
ax1.set_xticks([1, 10])
ax1.set_xticklabels(['1', '10'])
ax1.set_xlabel('$t_\\ast \;$  (Myr)')
ax1.set_ylim(Llims)
ax1.set_yscale('log')
ax1.set_yticks([0.1, 1, 10, 100, 1000])
ax1.set_yticklabels(['0.1', '1', '10', '100', '1000'])
ax1.set_ylabel('$L_{\\rm mm} \;$  (mJy at 150 pc)')


# for safety, copy over database file
os.system('cp -r DISKS.csv temp2.csv')

dref = 150.

# now load database
db = ascii.read('temp2.csv', format='csv', fast_reader=True)



# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'HJB') & 
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & 
         (db['FL_MULT'] != 'WJ') & (db['FL_MULT'] != 'HJ') & 
         (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') &
         (db['FL_logMs'] == 0) )



### cluster CDF comparisons

sfr = ['Oph', 'Tau', 'Lup', 'Cha', 'Usco']
sfr_name = ['Oph', 'Tau', 'Lup', 'Cha', 'USco']
nsamp = 100


# discretized log(Mstar) bins
binM_lo = [-1.0, -0.7, -0.4, -0.1, 0.2] 
binM_hi = [-0.7, -0.4, -0.1,  0.2, 0.5]

# model of the IMF (Chabrier)
mm = np.logspace(binM_lo[0], binM_hi[-1], 1024)
imf = 0.093 * np.exp(-0.5*(np.log10(mm) - np.log10(0.2))**2/(0.55**2))
imf[mm >= 1.] = 0.041*mm[mm >= 1.]**(-1.35)
int_imf_tot = np.trapz(imf, np.log10(mm))
imf /= int_imf_tot

# figure out how many disks you should draw into each bin, based on IMF
fperbin = np.zeros(len(binM_lo))
for i in range(len(binM_lo)):
    inbin = ((np.log10(mm) >= binM_lo[i]) & (np.log10(mm) < binM_hi[i]))
    fperbin[i] = np.trapz(imf[inbin], np.log10(mm[inbin]))
nb = np.round(fperbin * nsamp)
nperbin = nb.astype(np.int64)

col = ['C3', 'C1', 'C2', 'C0', 'C4']
for i in range(len(sfr)):

    # detections
    det = ((db['FL_B7'] == 0) & base & (db['SFR'] == sfr[i]))
    Mdet = db['logMs'][det]
    Ldet = db['F_B7'][det] * (db['DPC'][det] / dref)**2
    flags_det = np.zeros_like(Mdet)

    # upper limits
    lim = ((db['FL_B7'] == 1) & (db['FL_B6'] != 0) & base & \
           (db['SFR'] == sfr[i]))
    Mlim = db['logMs'][lim]
    Llim = db['LIM_B7'][lim] * (db['DPC'][lim] / dref)**2
    flags_lim = np.ones_like(Mlim)

    # combine
    M = np.ma.concatenate((Mdet, Mlim))
    L = np.ma.concatenate((Ldet, Llim))
    flags = np.ma.concatenate((flags_det, flags_lim))

    # construct a bunch of CDFs for each SFR
    if (sfr[i] == 'Usco'): 
        xmin = 0.58
    else: xmin = 0.8
    dxx = np.logspace(np.log10(xmin), 3, 512)
    dyy = np.zeros((512, 500))
    for ix in range(500):
        # randomly draw a distribution following the IMF sampling rules
        indices = np.array([], dtype=np.int64)
        for j in range(len(binM_lo)):
            inbin = (np.where( (M >= binM_lo[j]) & (M < binM_hi[j]) ))[0]
            draws = np.random.choice(inbin, nperbin[j])
            indices = np.append(indices, draws)
        Ldist = L[indices]
        delta = (flags[indices] == 1)

        # construct the CDF
        xx, yy, eyy, mukm = km_estimator(Ldist, delta)

        # interpolate the CDFs onto a common grid
        cint = interp1d(xx, yy, fill_value='extrapolate')
        tempy = cint(dxx)
        tempy[tempy < 0.] = 0.
        tempy[tempy > 1.] = 1.
        dyy[:,ix] = tempy

        # and plot the CDF
    #    ax0.plot(xx, yy, col[i], alpha=0.05, drawstyle='steps-post')

    ax0.fill_between(dxx, np.percentile(dyy, 2.5, axis=1), 
                     np.percentile(dyy, 97.5, axis=1), facecolor=col[i], 
                     alpha=0.5)

    #ax0.plot(dxx, np.percentile(dyy, 16., axis=1), 'k', drawstyle='steps-post')
    #ax0.plot(dxx, np.percentile(dyy, 50., axis=1), 'k', drawstyle='steps-post')
    #ax0.plot(dxx, np.percentile(dyy, 84., axis=1), 'k', drawstyle='steps-post')

    # annotations
    ax0.text(0.05, 0.33-0.065*i, sfr_name[i], transform=ax0.transAxes, 
             horizontalalignment='left', fontsize=8, fontweight='bold', 
             color=col[i])



### Lmm vs age

# detections
det = ((db['FL_B7'] == 0) & base)
Mdet = db['logMs'][det]
Ldet = db['F_B7'][det] * (db['DPC'][det] / dref)**2
tdet = 10.**(db['logt'][det]) / 1e6
flags_det = np.zeros_like(Mdet)

# upper limits
lim = ((db['FL_B7'] == 1) & (db['FL_B6'] != 0) & base)
Mlim = db['logMs'][lim]
Llim = db['LIM_B7'][lim] * (db['DPC'][lim] / dref)**2
tlim = 10.**(db['logt'][lim]) / 1e6
flags_lim = np.ones_like(Mlim)

# combine
M = np.ma.concatenate((Mdet, Mlim))
L = np.ma.concatenate((Ldet, Llim))
t = np.ma.concatenate((tdet, tlim))
flags = np.ma.concatenate((flags_det, flags_lim))

# sort by Ms
sort_indices = np.argsort(M)
s_M = M[sort_indices]
s_L = L[sort_indices]
s_t = t[sort_indices]
s_flags = flags[sort_indices]
norm = mpl.colors.Normalize(vmin=np.log10(0.1), vmax=np.log10(1.2))
cmap = mpl.cm.get_cmap('coolwarm_r')
rgba = cmap(norm(s_M))
for i in range(len(s_M)):
    lcol = rgba[i]
    lcol[3] = 0.9
    if (s_flags[i] == 0):
        ax1.plot(s_t[i], s_L[i], 'o', markerfacecolor=lcol, 
                 markeredgecolor='None', markersize=4.5)
    else:
        ax1.errorbar(s_t[i], s_L[i], yerr=0.32*s_L[i], uplims=True,
                     marker='None', capsize=2, color=lcol, alpha=0.45, 
                     linestyle='None')


# colorbar
cax = fig.add_axes([0.905, 0.19, 0.02, 0.79])
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, 
                                orientation='vertical') 
cax.yaxis.tick_right()
cax.tick_params(axis='both', which='major', labelsize=7)



# annotations
ax0.text(0.88, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.88, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('age.pdf')
fig.clf()
