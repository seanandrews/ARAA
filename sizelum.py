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

# selection for both B7 flux density and effective size
ok = ((db['FL_B7'] == 0) & (db['FL_R7'] == 0))
L = db['F_B7'][ok]*(db['DPC'][ok]/140.)**2
R = db['R7'][ok]*(db['DPC'][ok])


# set up plot
fig = plt.figure(figsize=(7., 2.5))
gs = gridspec.GridSpec(1, 2)
Llims = [0.5, 2000.]
Rlims = [2.0, 500.]


# dust size versus luminosity
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_ylim(Rlims)
ax0.set_yscale('log')
ax0.plot(L, R, 'oC0')

ax0.set_xticks([1, 10, 100, 1000])
ax0.set_xticklabels(['1', '10', '100', '1000'])
ax0.set_yticks([10, 100])
ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('continuum $R_{\\rm eff}$  (au)')
ax0.set_xlabel('$L_{\\rm mm}$  (mJy at 140 pc)')


# dust size versus gas size
ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlim(Llims)
ax1.set_xscale('log')
ax1.set_ylim(Rlims)
ax1.set_yscale('log')
ax1.plot(L, R, 'oC0')

ax1.set_xticks([1, 10, 100, 1000])
ax1.set_xticklabels(['1', '10', '100', '1000'])
ax1.set_yticks([10, 100])
ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('continuum $R_{\\rm eff}$  (au)')
ax1.set_xlabel('CO $R_{\\rm eff}$  (au)')


fig.subplots_adjust(wspace=0.30, hspace=0.0)
fig.subplots_adjust(left=0.07, right=0.93, bottom=0.15, top=0.99)
fig.savefig('sizes.pdf')
fig.clf()

