import numpy as np
import os
import sys
from astropy.io import ascii
import cloudpickle as cp
from scipy.interpolate import interp1d
import scipy.integrate as sci

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib.ticker import FixedLocator
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
elims = [1.5, 4.7]
Rlims = [0., 150.]
Llims = [0.05, 5000.]
Llims = [np.log10(0.05), np.log10(5000)]

# panel (a) setups  [alphamm versus Lmm]
ax0.set_xlim(Llims)
ax0.set_xticks([-1, 0, 1, 2, 3])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm 1.3 \, mm} \;$ (mJy at 150 pc)')
ax0.yaxis.get_ticklocs(minor=True)
ax0.minorticks_on()
ax0.set_ylim(elims)
ax0.set_ylabel('$\\alpha_{\\rm mm}$')


# panel (b) setups  [epsilon(r) profiles]
ax1.set_xlim(Rlims)
ax1.set_xlabel('$r \;$ (au)')
ax1.xaxis.get_ticklocs(minor=True)
ax1.yaxis.get_ticklocs(minor=True)
ax1.minorticks_on()
ax1.set_ylim(elims)
ax1.set_ylabel('$\\varepsilon \;$ (1--9 mm)')


### Set some constants
cc = 2.9979e10


ind = '67'

if (ind == '67'): alp_lbl = '(B6 / B7)'
if (ind == '36'): alp_lbl = '(B3 / B6)'

# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# downselect targets with both a robust index and B6 flux density
base = ( (db['SED'] != 'I') & (db['SED'] != 'III') & 
         (db['SED'] != 'DEBRIS') )

d_ref, nu_ref, alp = 150., 225., 2.3

# extract sub-sample
sub = ( (db['FL_A'+ind] == 0) & (db['FL_B6'] == 0) & base )
L6  = db['F_B6'][sub] * (db['DPC'][sub] / d_ref)**2 * \
      (nu_ref / db['nu_B6'][sub])**alp
eL6 = np.sqrt( (nu_ref/db['nu_B6'][sub])**(2.*alp) * \
               (db['eF_B6'][sub]**2 + (0.1*db['F_B6'][sub])**2) * \
               (db['DPC'][sub]/d_ref)**4 + \
               ( ((nu_ref / db['nu_B6'][sub])**alp * \
                  db['F_B6'][sub]*(2.*db['DPC'][sub]/d_ref**2) * \
                 0.5*(db['EDPC_H'][sub]+db['EDPC_L'][sub]) )**2 ) )
amm = db['A'+ind][sub]
eamm_hi = db['eA'+ind+'_hi'][sub]
eamm_lo = db['eA'+ind+'_lo'][sub]
names = db['NAME'][sub]

eL6 = [np.log10(L6)-np.log10(L6-eL6), np.log10(L6+eL6)-np.log10(L6)]
L6 = np.log10(L6)
#ax0.errorbar(L6, amm, xerr=eL6, yerr=[eamm_lo, eamm_hi], 
#             marker='o', color='gray', markersize=3, 
#             linestyle='None', elinewidth=0.5, alpha=0.35, zorder=1)
ax0.plot(L6, amm, marker='o', color='gray', markersize=3, alpha=0.45, 
         linestyle='None', zorder=1)

ax0.text(0.68, 0.90, '0.9--1.3 mm', transform=ax0.transAxes, 
         fontsize=9, color='gray')
ax0.text(0.68, 0.83, '1.3--3 mm', transform=ax0.transAxes, 
         fontsize=9, color='C1')




ind = '36'

if (ind == '67'): alp_lbl = '(B6 / B7)'
if (ind == '36'): alp_lbl = '(B3 / B6)'

# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# downselect targets with both a robust index and B6 flux density
base = ( (db['SED'] != 'I') & (db['SED'] != 'III') &
         (db['SED'] != 'DEBRIS') )

d_ref, nu_ref, alp = 150., 225., 2.3

# extract sub-sample
sub = ( (db['FL_A'+ind] == 0) & (db['FL_B6'] == 0) & base )
L6  = db['F_B6'][sub] * (db['DPC'][sub] / d_ref)**2 * \
      (nu_ref / db['nu_B6'][sub])**alp
eL6 = np.sqrt( (nu_ref/db['nu_B6'][sub])**(2.*alp) * \
               (db['eF_B6'][sub]**2 + (0.1*db['F_B6'][sub])**2) * \
               (db['DPC'][sub]/d_ref)**4 + \
               ( ((nu_ref / db['nu_B6'][sub])**alp * \
                  db['F_B6'][sub]*(2.*db['DPC'][sub]/d_ref**2) * \
                 0.5*(db['EDPC_H'][sub]+db['EDPC_L'][sub]) )**2 ) )
amm = db['A'+ind][sub]
eamm_hi = db['eA'+ind+'_hi'][sub]
eamm_lo = db['eA'+ind+'_lo'][sub]
names = db['NAME'][sub]

eL6 = [np.log10(L6)-np.log10(L6-eL6), np.log10(L6+eL6)-np.log10(L6)]
L6 = np.log10(L6)
ax0.errorbar(L6, amm, xerr=eL6, yerr=[eamm_lo, eamm_hi],                
             marker='o', color='C1', markersize=3, 
             linestyle='None', elinewidth=1.0, alpha=0.65, zorder=2)




# highlight individuals
uztau = (names == 'UZ_Tau_E')
ax0.errorbar(L6[uztau], amm[uztau], yerr=[eamm_lo[uztau], eamm_hi[uztau]],
             marker='o', color='C4', markersize=5, linestyle='None',
             elinewidth=2., zorder=3)

drtau = (names == 'DR_Tau')
ax0.errorbar(L6[drtau], amm[drtau], yerr=[eamm_lo[drtau], eamm_hi[drtau]], 
             marker='o', color='C3', markersize=5, linestyle='None', 
             elinewidth=2., zorder=4)

fttau = (names == 'FT_Tau')
ax0.errorbar(L6[fttau], amm[fttau], yerr=[eamm_lo[fttau], eamm_hi[fttau]],
             marker='o', color='C2', markersize=5, linestyle='None',
             elinewidth=2., zorder=5)

as209 = (names == 'AS_209')
ax0.errorbar(L6[as209], amm[as209], yerr=[eamm_lo[as209], eamm_hi[as209]],
             marker='o', color='C0', markersize=5, linestyle='None',
             elinewidth=2., zorder=6)



# B6/B7 indices
ok67 = ((db['FL_B7'] == 0) & (db['FL_B6'] == 0) & (db['FL_A67'] == 0) & base)
names = db['NAME'][ok67]
num67 = len(names)
a67_samples = []
for i in range(num67):
    a67_ = np.load('outputs/'+names[i]+'.alpha67.posterior.npz')['amm']
    a67_samples = np.append(a67_samples, a67_)

N, bins = np.histogram(a67_samples, range=[0, 5], bins=100, density=True)
ax0.plot(1.5*N+Llims[0]-0.10, bins[1:], 'gray')


ok36 = ((db['FL_B6'] == 0) & (db['FL_B3'] == 0) & (db['FL_A36'] == 0) & base)
names = db['NAME'][ok36]
num36 = len(names)
a36_samples = []
for i in range(num36):
    a36_ = np.load('outputs/'+names[i]+'.alpha36.posterior.npz')['amm']
    a36_samples = np.append(a36_samples, a36_)

N, bins = np.histogram(a36_samples, range=[0, 5], bins=100, density=True)
ax0.plot(1.5*N+Llims[0], bins[1:], 'C1')




ax0.fill_between(Llims, [3.7, 3.7], [4.0, 4.0], facecolor='y', alpha=0.2)
ax0.text(np.log10(0.3), 3.83, 'ISM', fontsize=10, color='goldenrod', 
         va='center')




### UZ Tau from Anjali
rau, alpl, alph = np.loadtxt('data/uztau_epsilon.dat').T
ax1.fill_between(rau, alpl, alph, facecolor='C4', alpha=0.25)
ax1.plot(rau, 0.5*(alpl+alph), '-C4')

### Load the data from Tazzari+ 2016
with open("data/tazzari_profiles.dat", 'rb') as f:
    data = cp.load(f, encoding='latin1')

### index profiles
name = ['DRTau', 'FTTau', 'AS209']
col = ['C3', 'C2', 'C0']
for i in range(len(name)):
    rau, wl = data[name[i]]['gridrad'], data[name[i]]['wle']
    a, b = 0, len(wl)-1
    Ia  = data[name[i]]['intensity'][a][:,1]
    eIa = 0.5*(data[name[i]]['intensity'][a][:,2] - \
               data[name[i]]['intensity'][a][:,0])
    Ib  = data[name[i]]['intensity'][b][:,1]
    eIb = 0.5*(data[name[i]]['intensity'][b][:,2] - \
               data[name[i]]['intensity'][b][:,0])

    eps  = np.log(Ia/Ib) / np.log(wl[b]/wl[a])
    eeps = np.sqrt( (1./(Ia*np.log(wl[b]/wl[a])))**2 * eIa**2 + \
                    (1./(Ib*np.log(wl[b]/wl[a])))**2 * eIb**2 )
    ax1.fill_between(rau, eps+eeps, eps-eeps, facecolor=col[i], alpha=0.5)
    ax1.plot(rau, eps, '-'+col[i])

ax1.text(0.93, 0.07, 'UZ Tau E', transform=ax1.transAxes, ha='right', 
         fontsize=8, color='C4')
ax1.text(0.93, 0.13, 'AS 209', transform=ax1.transAxes, ha='right', 
         fontsize=8, color='C0')
ax1.text(0.93, 0.19, 'FT Tau', transform=ax1.transAxes, ha='right',  
         fontsize=8, color='C2')
ax1.text(0.93, 0.25, 'DR Tau', transform=ax1.transAxes, ha='right',  
         fontsize=8, color='C3')




ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')


fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.10, right=0.90, bottom=0.19, top=0.98)
fig.savefig('mm_indices.pdf')
fig.clf()
