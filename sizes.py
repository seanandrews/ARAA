import numpy as np
import os
import sys
from astropy.io import ascii

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# set up plot
fig = plt.figure(figsize=(6.33, 2.2))
gs = gridspec.GridSpec(1, 2)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])

# set up axes, labels
Llims = [0.5, 2000.]
Rlims = [2.0, 500.]

# Panel (a) setups  [size-luminosity relation]
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([1, 10, 100, 1000])
ax0.set_xticklabels(['1', '10', '100', '1000'])
ax0.set_xlabel('$L_{\\rm mm}$  (mJy at 140 pc)')

ax0.set_ylim(Rlims)
ax0.set_yscale('log')
ax0.set_yticks([10, 100])
ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('continuum $R_{\\rm eff}$  (au)')

# Panel (b) setups  [R_dust versus R_gas]
ax1.set_xlim(Rlims)
ax1.set_xscale('log')
ax1.set_xticks([1, 10, 100, 1000])
ax1.set_xticklabels(['1', '10', '100', '1000'])
ax1.set_xlabel('CO $R_{\\rm eff}$  (au)')

ax1.set_ylim(Rlims)
ax1.set_yscale('log')
ax1.set_yticks([10, 100])
ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('continuum $R_{\\rm eff}$  (au)')


### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# baseline selections
base = ( (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'T') & 
         (db['FL_MULT'] != 'T') & (db['SED'] != 'III') &
         (db['SED'] != 'DEBRIS') & (db['SFR'] != 'Oph') )



# selection for both B7 flux density and effective size
ok = ((db['FL_B7'] == 0) & (db['FL_R7'] == 0))
L = db['F_B7'][ok]*(db['DPC'][ok]/140.)**2
R = db['R7'][ok]*(db['DPC'][ok])


# dust size versus luminosity
ax0.plot(L, R, 'oC0', markersize=4.0)


# dust size versus gas size
ax1.plot(L, R, 'oC0', markersize=4.0)


fig.subplots_adjust(wspace=0.40)
fig.subplots_adjust(left=0.105, right=0.895, bottom=0.19, top=0.98)
fig.savefig('sizes.pdf')
fig.clf()

