import numpy as np
import os
import sys
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rc('font', size=9)


# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# set up plot
fig = plt.figure(figsize=(7., 2.25))
gs = gridspec.GridSpec(1, 3)
Llims = [0.2, 5000.]
Rlims = [2.0, 500.]
Mlims = [0.02, 5.]
Alims = [0., 6.]

# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & \
         (db['FL_MULT'] != 'J') & (db['SED'] != 'III') )


# Lmm vs Mstar
ok = ((db['FL_B7'] == 0) & (db['FL_logMs'] == 0) & base)
L = db['F_B7'][ok]*(db['DPC'][ok]/140.)**2
M = 10.**(db['logMs'][ok])

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_ylim(Mlims)
ax0.set_yscale('log')
ax0.plot(L, M, '.C0')

ax0.set_xticks([1, 10, 100, 1000])
ax0.set_xticklabels(['1', '10', '100', '1000'])
ax0.set_yticks([0.1, 1])
ax0.set_yticklabels(['0.1', '1'])
ax0.set_ylabel('$M_\\ast$  ($M_\odot$)')
ax0.set_xlabel('$L_{\\rm mm}$  (mJy at 140 pc)')


# Reff(continuum) vs Mstar
ok = ((db['FL_R7'] == 0) & (db['FL_logMs'] == 0) & base)
M = 10.**(db['logMs'][ok])
R = db['R7'][ok]*(db['DPC'][ok])

ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlim(Rlims)
ax1.set_xscale('log')
ax1.set_ylim(Mlims)
ax1.set_yscale('log')
ax1.plot(R, M, '.C0')

ax1.set_xticks([10, 100])
ax1.set_xticklabels(['10', '100'])
ax1.set_yticklabels([])
ax1.set_xlabel('continuum $R_{\\rm eff}$  (au)')


# alpha(continuum) vs Mstar
ok = ((db['FL_A67'] == 0) & (db['FL_logMs'] == 0) & base)
M = 10.**(db['logMs'][ok])
amm = db['A67'][ok]

ax2 = fig.add_subplot(gs[0, 2])
ax2.set_xlim(Alims)
ax2.set_ylim(Mlims)
ax2.set_yscale('log')
ax2.plot(amm, M, '.C0')

#ax2.set_xticks([10, 100])
#ax2.set_xticklabels(['10', '100'])
ax2.set_yticklabels([])
ax2.set_xlabel('$\\alpha_{\\rm mm}$')



fig.subplots_adjust(wspace=0.05, hspace=0.0)
fig.subplots_adjust(left=0.07, right=0.93, bottom=0.17, top=0.99)
fig.savefig('mstar.pdf')
fig.clf()

