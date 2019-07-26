import numpy as np
import os
import sys
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from km_estimator import km_estimator
plt.rc('font', size=9)


# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

dref = 150.

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# set up plots
fig = plt.figure(figsize=(6.33, 2.3))
gs = gridspec.GridSpec(1, 2)
Llims = [0.05, 5000.]
Rlims = [2.0, 500.]
Tlims = [0.2, 50.]

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm mm}$  (mJy at 150 pc)')
ax0.set_ylim([0., 1.])
ax0.set_ylabel('$p(> L_{\\rm mm})$')

ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlim(Tlims)
ax1.set_xscale('log')
ax1.set_xticks([1, 10])
ax1.set_xticklabels(['1', '10'])
ax1.set_xlabel('$t_\\ast$  (Myr)')
ax1.set_ylim(Llims)
ax1.set_yscale('log')
ax1.set_yticks([0.1, 1, 10, 100, 1000])
ax1.set_yticklabels(['0.1', '1', '10', '100', '1000'])
ax1.set_ylabel('$L_{\\rm mm}$  (mJy at 150 pc)')


# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & \
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & \
         (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') &
         (db['FL_logMs'] == 0) )



### cluster CDF comparisons

sfr = ['Tau', 'Lup', 'Cha', 'Usco']

nsamp = 100


# discretized log(Mstar) bins
binM_lo = [-1.0, -0.7, -0.4, -0.1, 0.2] 
binM_hi = [-0.7, -0.4, -0.1,  0.2, 0.5]

# model of the IMF (Chabrier 2005)
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

col = ['C0', 'C1', 'C2', 'C3']
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

        # and plot the CDF
        ax0.plot(xx, yy, col[i], alpha=0.05, drawstyle='steps-post')











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
    mrkr = 'o' if (s_flags[i] == 0) else 'x'
    ax1.plot(s_t[i], s_L[i], mrkr, markerfacecolor=lcol, 
             markeredgecolor='None', markersize=5)


fig.subplots_adjust(wspace=0.35)
fig.subplots_adjust(left=0.09, right=0.91, bottom=0.17, top=0.98)
fig.savefig('age.pdf')
fig.clf()
