import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
from astropy.io import ascii
plt.rc('font', size=9)


ind = '67'


if (ind == '67'): alp_lbl = '(B6 / B7)'
if (ind == '36'): alp_lbl = '(B3 / B6)'


# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)


# downselect targets with both a robust index and B6 flux density
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & \
         (db['FL_MULT'] != 'J') & (db['SED'] != 'III') )

base = ( (db['FL_MULT'] != 'J') & (db['SED'] != 'III') )


# choose a reference distance (in pc)
d0 = 150.

# set up plots
fig = plt.figure(figsize=(6.33, 4.0))
gs = gridspec.GridSpec(2, 2)
Llims = [0.2, 5000]
Rlims = [2, 500]
Mlims = [0.02, 5]
tlims = [0.05, 20]
alims = [1.0, 5.0]


# L_B6 versus alpha
#
# extract sub-sample 
Lalp_wind = ( (db['FL_A'+ind] == 0) & (db['FL_B6'] == 0) )
name = db['NAME'][base & Lalp_wind]
dpc  = db['DPC'][base & Lalp_wind]
F_B6 = db['F_B6'][base & Lalp_wind]
alp  = db['A'+ind][base & Lalp_wind]
print(len(name))

# convert to B6 luminosities (in mJy at some reference distance)
L_B6 = F_B6 * (dpc/d0)**2

# plot
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_ylim(alims)
ax0.plot(L_B6, alp, 'o', markersize=2.0)
ax0.set_xticks([1, 10, 100, 1000])
ax0.set_xticklabels(['1', '10', '100', '1000'])
ax0.set_xlabel('$F_\\nu \\times (d / {\\rm 150 pc})^2$ (mJy)')
ax0.set_ylabel('$\\alpha$ '+alp_lbl)


# R_B7 versus alpha
#
# extract sub-sample 
Ralp_wind = ( (db['FL_A'+ind] == 0) & (db['FL_R7'] == 0) )
name = db['NAME'][base & Ralp_wind]
dpc  = db['DPC'][base & Ralp_wind]
R_B7 = db['R7'][base & Ralp_wind]
alp  = db['A'+ind][base & Ralp_wind]
print(len(name))

# convert to radii in au units
Rau_B7 = R_B7 * dpc

# plot
ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlim(Rlims)
ax1.set_xscale('log')
ax1.set_ylim(alims)
ax1.plot(Rau_B7, alp, 'o', markersize=2.0)
ax1.set_xticks([10, 100])
ax1.set_xticklabels(['10', '100'])
ax1.set_xlabel('$R_{\\rm eff}$ (au)')
ax1.set_ylabel('$\\alpha$ '+alp_lbl)


# Mstar versus alpha
#
# extract sub-sample 
Malp_wind = ( (db['FL_A'+ind] == 0) & (db['FL_logMs'] == 0) )
name = db['NAME'][base & Malp_wind]
dpc  = db['DPC'][base & Malp_wind]
logM = db['logMs'][base & Malp_wind]
alp  = db['A'+ind][base & Malp_wind]
print(len(name))

# plot
ax2 = fig.add_subplot(gs[1, 0])
ax2.set_xlim(Mlims)
ax2.set_xscale('log')
ax2.set_ylim(alims)
ax2.plot(10.**logM, alp, 'o', markersize=2.0)
ax2.set_xticks([0.1, 1])
ax2.set_xticklabels(['0.1', '1'])
ax2.set_xlabel('$M_\\ast$ (M$_\odot$)')
ax2.set_ylabel('$\\alpha$ '+alp_lbl)


# age versus alpha
#
# extract sub-sample 
Talp_wind = ( (db['FL_A'+ind] == 0) & (db['FL_logMs'] == 0) )
name = db['NAME'][base & Talp_wind]
dpc  = db['DPC'][base & Talp_wind]
logt = db['logt'][base & Talp_wind]
alp  = db['A'+ind][base & Talp_wind]
print(len(name))

# plot
ax3 = fig.add_subplot(gs[1, 1])
ax3.set_xlim(tlims)
ax3.set_xscale('log')
ax3.set_ylim(alims)
ax3.plot(1e-6 * 10.**logt, alp, 'o', markersize=2.0)
ax3.set_xticks([0.1, 1, 10])
ax3.set_xticklabels(['0.1', '1', '10'])
ax3.set_xlabel('age (Myr)')
ax3.set_ylabel('$\\alpha$ '+alp_lbl)



fig.subplots_adjust(wspace=0.30, hspace=0.35)
fig.subplots_adjust(left=0.07, right=0.93, bottom=0.09, top=0.98)
fig.savefig('alpha_'+ind+'.pdf')
fig.clf()
