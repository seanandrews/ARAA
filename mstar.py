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

# set up plots
fig = plt.figure(figsize=(6.33, 2.6))
gs = gridspec.GridSpec(1, 1)
Llims = [0.125, 8000.]
Mlims = [0.0125, 8.]
dref = 150.

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Mlims)
ax0.set_xscale('log')
ax0.set_ylim(Llims)
ax0.set_yscale('log')
ax0.set_yticks([1, 10, 100, 1000])
ax0.set_yticklabels(['1', '10', '100', '1000'])
ax0.set_xticks([0.1, 1])
ax0.set_xticklabels(['0.1', '1'])
ax0.set_xlabel('$M_\\ast$  ($M_\odot$)')
ax0.set_ylabel('$L_{\\rm mm}$  (mJy at 140 pc)')


# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & \
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & \
         (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') & 
         (db['SFR'] == 'Cha') )


### Lmm vs Mstar (all in B7)

# B7 measurements available
ok = ((db['FL_B7'] == 0) & base)
L7 = db['F_B7'][ok] * (db['DPC'][ok] / dref)**2	
M7 = 10.**(db['logMs'][ok])
id7 = db['NAME'][ok]

# B7 upper limits
ok = ((db['FL_B7'] == 1) & base)
limL7 = db['LIM_B7'][ok] * (db['DPC'][ok] / dref)**2 
limM7 = 10.**(db['logMs'][ok])


#print(len(L6))

#L = np.concatenate((L7, L6))
#M = np.concatenate((M7, M6))

ax0.plot(limM7, limL7, '.C1')
ax0.plot(M7, L7, '.C0')



fig.subplots_adjust(wspace=0.31, hspace=0.0)
fig.subplots_adjust(left=0.2, right=0.8, bottom=0.15, top=0.98)
fig.savefig('mstar.pdf')
fig.clf()

