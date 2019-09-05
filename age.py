import numpy as np
import os
import sys
from astropy.io import ascii
from scipy.interpolate import interp1d
from km_estimator import km_estimator

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
#import matplotlib.ticker
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
Llims = [0.05, 5000.]
Tlims = [0.3333, 30.]

ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm 0.9 \, mm} \;$  (mJy at 150 pc)')
ax0.set_ylim([0., 1.0])
ax0.set_ylabel('$p(> L_{{\\rm 0.9 \, mm}})$')
locmin = mpl.ticker.LinearLocator(numticks=26)
ax0.yaxis.set_minor_locator(locmin)
ax0.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

locmin = mpl.ticker.LogLocator(base=10., 
                               subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                               numticks=5)
ax0.xaxis.set_minor_locator(locmin)
ax0.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())


ax1.set_xlim(Tlims)
ax1.set_xscale('log')
ax1.set_xticks([1, 10])
ax1.set_xticklabels(['1', '10'])
ax1.set_xlabel('$t_\\ast \;$  (Myr)')
ax1.set_ylim(Llims)
ax1.set_yscale('log')
ax1.set_yticks([0.1, 1, 10, 100, 1000])
ax1.set_yticklabels(['0.1', '1', '10', '100', '1000'])
ax1.set_ylabel('$L_{\\rm 0.9 \, mm} \;$  (mJy at 150 pc)')


# for safety, copy over database file
os.system('cp -r DISKS.csv temp2.csv')

dref = 150.

# now load database
db = ascii.read('temp2.csv', format='csv', fast_reader=True)



# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'HJB') & 
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & 
         (db['FL_MULT'] != 'WJ') & (db['FL_MULT'] != 'HJ') & 
         (db['FL_MULT'] != 'HCB') & (db['FL_MULT'] != 'HC') & 
         (db['SED'] != 'I') & (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') &
         (db['FL_logMs'] == 0) & (db['FL_logt'] != 9) )


### Set some constants
d_ref, nu_ref, alp = 150., 340., 2.3


# calculate luminosities, upper limits, masses, and ages + uncertainties
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
logM = db['logMs']
age = 10.**db['logt']
eage_hi = 10.**(db['elogt_H']+db['logt']) - age
eage_lo = age - 10.**(db['logt']-db['elogt_L'])
reg = db['SFR']
names = db['NAME']


# targets with B7 detections
detB7 = ( (db['FL_B7'] == 0) & base )
L_detB7 = L7[detB7]
eL_detB7 = eL7[detB7]
M_detB7 = logM[detB7]
age_detB7 = age[detB7] / 1e6
eage_hi_detB7 = eage_hi[detB7] / 1e6
eage_lo_detB7 = eage_lo[detB7] / 1e6
sfr_detB7 = reg[detB7]
flags_detB7 = np.zeros_like(L_detB7)


# targets with **only** B6 detections (i.e., no B7 or B7 limit)
detB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & base )
L_detB6 = L6[detB6]
eL_detB6 = eL6[detB6]
M_detB6 = logM[detB6]
age_detB6 = age[detB6] / 1e6
eage_hi_detB6 = eage_hi[detB6] / 1e6
eage_lo_detB6 = eage_lo[detB6] / 1e6
sfr_detB6 = reg[detB6]
flags_detB6 = np.zeros_like(L_detB6)


# targets with **only** limits or missing data
# (there should be **no entry** without a limit in at least B6 or B7)
lims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & base )
dlims = np.ma.column_stack( (limL7[lims], limL6[lims]) )
L_lims = np.ma.min(dlims, 1)
eL_lims = np.zeros_like(L_lims)
M_lims = logM[lims]
age_lims = age[lims] / 1e6
eage_hi_lims = eage_hi[lims] / 1e6
eage_lo_lims = eage_lo[lims] / 1e6
sfr_lims = reg[lims]
flags_lims = np.ones_like(L_lims)
name_lims = names[lims]


# combine everything
Lmm = np.ma.concatenate( (L_detB7, L_detB6, L_lims) )
eLmm = np.ma.concatenate( (eL_detB7, eL_detB6, eL_lims) )
Ms = np.ma.concatenate( (M_detB7, M_detB6, M_lims) )
ages = np.ma.concatenate( (age_detB7, age_detB6, age_lims) )
eages_h = np.ma.concatenate( (eage_hi_detB7, eage_hi_detB6, eage_hi_lims) )
eages_l = np.ma.concatenate( (eage_lo_detB7, eage_lo_detB6, eage_lo_lims) )
sfrs = np.ma.concatenate( (sfr_detB7, sfr_detB6, sfr_lims) )
dflags = np.ma.concatenate( (flags_detB7, flags_detB6, flags_lims) )





### cluster CDF comparisons

sfr = ['Oph', 'Tau', 'Lup', 'Cha', 'IC348', 'Usco']
sfr_name = ['Oph', 'Tau', 'Lup', 'Cha', 'IC 348', 'USco']
med_age = ['1.4', '1.5', '2.1', '2.6', '2.8', '5.2']
nsamp = 100


# discretized log(Mstar) bins
binM_lo = [-1.0, -0.8, -0.6, -0.4, -0.2]	#, 0.0]
binM_hi = [-0.8, -0.6, -0.4, -0.2,  0.0]	#, 0.24]


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

col = ['C3', 'C1', 'C2', 'C9', 'C0', 'C4']
for i in range(len(sfr)):
    print(sfr[i])

    M = Ms[sfrs == sfr[i]]
    L = Lmm[sfrs == sfr[i]]
    flags = dflags[sfrs == sfr[i]]

    # construct a bunch of CDFs for each SFR
    xmin = 0.8
    if (sfr[i] == 'Usco'): xmin = 0.58
    if (sfr[i] == 'IC348'): xmin= 5.
    dxx = np.logspace(np.log10(xmin), 3, 512)
    dyy = np.zeros((512, 500))
    for ix in range(500):
        # randomly draw a distribution following the IMF sampling rules
        indices = np.array([], dtype=np.int64)
        for j in range(len(binM_lo)):
            inbin = (np.where( (M >= binM_lo[j]) & (M < binM_hi[j]) ))[0]
            if (len(inbin) > 0):
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
      #  ax0.plot(xx, yy, col[i], alpha=0.05, drawstyle='steps-post')

    ax0.fill_between(dxx, np.percentile(dyy, 2.5, axis=1), 
                     np.percentile(dyy, 97.5, axis=1), facecolor=col[i], 
                     alpha=0.45)
    ax0.plot(dxx, np.percentile(dyy, 50, axis=1), '-'+col[i], alpha=0.65)

    #ax0.plot(dxx, np.percentile(dyy, 16., axis=1), 'k', drawstyle='steps-post')
    #ax0.plot(dxx, np.percentile(dyy, 50., axis=1), 'k', drawstyle='steps-post')
    #ax0.plot(dxx, np.percentile(dyy, 84., axis=1), 'k', drawstyle='steps-post')

    # annotations
    ax0.text(0.58, 0.89-0.065*i, sfr_name[i]+' -', 
             transform=ax0.transAxes, 
             horizontalalignment='left', fontsize=8, fontweight='bold', 
             color=col[i])
    ax0.text(0.76, 0.89-0.065*i, med_age[i]+' Myr',
             transform=ax0.transAxes,
             horizontalalignment='left', fontsize=8, fontweight='bold', 
             color=col[i])




### Lmm vs age

# add limits on Ms
cut = ( (Ms >= binM_lo[0]) & (Ms <= binM_hi[-1]) )
c_Lmm = Lmm[cut]
c_eLmm = eLmm[cut]
c_Ms = Ms[cut]
c_ages = ages[cut]
c_eages_h = eages_h[cut]
c_eages_l = eages_l[cut]
c_dflags = dflags[cut] 

# sort by Ms ***
sort_indices = np.argsort(c_Ms)
s_M = c_Ms[sort_indices]
s_L = c_Lmm[sort_indices]
s_eL = c_eLmm[sort_indices]
s_t = c_ages[sort_indices]
s_eth = c_eages_h[sort_indices]
s_etl = c_eages_l[sort_indices]
s_flags = c_dflags[sort_indices]
norm = mpl.colors.Normalize(vmin=binM_lo[0], vmax=binM_hi[-1])
cmap = mpl.cm.get_cmap('cool')
rgba = cmap(norm(s_M))
for i in range(len(s_M)):
    lcol = rgba[i]
    lcol[3] = 0.9
    if (s_flags[i] == 1):
        ax1.errorbar(s_t[i], s_L[i], yerr=0.32*s_L[i], uplims=True,
                     marker='None', capsize=2, color=lcol, alpha=0.45,
                     linestyle='None')
    else:
        ax1.plot(s_t[i], s_L[i], 'o', markerfacecolor=lcol,
                 markeredgecolor='None', markersize=4.5)


# colorbar
cax = fig.add_axes([0.905, 0.19, 0.02, 0.785])
nt = [np.log10(0.1), np.log10(0.2), np.log10(0.4), np.log10(0.6), np.log10(0.8), np.log10(1.0), np.log10(1.2)]
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, 
                                orientation='vertical', ticks=nt) 
cax.yaxis.tick_right()
cax.tick_params(axis='both', which='major', labelsize=10)
cax.set_yticklabels(['0.1', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2'])
cax.set_ylabel('$M_\\ast \;$ (M$_\\odot$)', fontsize=10, labelpad=8)




# annotations
ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.09, right=0.89, bottom=0.19, top=0.975)
fig.savefig('age.pdf')
fig.clf()
