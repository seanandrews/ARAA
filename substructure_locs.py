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
fig = plt.figure(figsize=(6.33, 2.25))
gs = gridspec.GridSpec(1, 2)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])

# set up axes, labels
Llims = [0.05, 80.]
Rlims = [2.0, 500.]
dlims = [0., 2.]
Mlims = [0.0125, 8.]

# Panel (a) setups  [feature location versus Lstar]
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10])
ax0.set_xticklabels(['0.1', '1', '10'])
ax0.set_xlabel('$L_{\\ast} \;$ (L$_\odot$)')

ax0.set_ylim(Rlims)
ax0.set_yscale('log')
ax0.set_yticks([10, 100])
ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('feature location (au)')

# Panel (b) setups  [mass-size relation]
ax1.set_xlim(Mlims)
ax1.set_xscale('log')
ax1.set_xticks([0.1, 1])
ax1.set_xticklabels(['0.1', '1'])
ax1.set_xlabel('$M_{\\ast} \;$ (M$_\odot$)')

ax1.set_ylim(dlims)
#ax1.set_yscale('log')
#ax1.set_yticks([10, 100])
#ax1.set_yticklabels(['10', '100'])
#ax1.set_ylabel('$R_{\\rm mm} \;$ (au)')


### Load the database

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# baseline selections
#base = ( (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'T') &
#         (db['FL_MULT'] != 'T') & (db['SED'] != 'III') &
#         (db['SED'] != 'DEBRIS') & (db['SFR'] != 'Oph') &
#         (db['FL_logMs'] == 0) )


# load substructure locations file
ss = ascii.read('data/substructure_locations.txt')


# grab Lstar and uncertainties for each substructure from database
L, eLhi, eLlo = np.zeros(len(ss)), np.zeros(len(ss)), np.zeros(len(ss))
for i in range(len(ss)):
    ind = np.where(db['NAME'] == ss['name'][i])[0][0]
    L[i] = 10.**(db['logLs'][ind])
    eLhi[i] = 10.**(db['elogLs_H'][ind]+db['logLs'][ind])-L[i]
    eLlo[i] = L[i] - 10.**(db['logLs'][ind]-db['elogLs_L'][ind])


# overlay some models
Lgrid = np.logspace(-2, 2, 1024)
sigSB, au, Lsun = 5.67051e-5, 1.496e13, 3.826e33
T_co = 21.
r_co = np.sqrt(0.02 * Lgrid * Lsun / (8 * np.pi * sigSB * T_co**4)) / au

T_nn = 19.
r_nn = np.sqrt(0.02 * Lgrid * Lsun / (8 * np.pi * sigSB * T_nn**4)) / au


ax0.plot(Lgrid, r_co, 'k')
ax0.plot(Lgrid, r_nn, 'k')




# gaps
gap  = (ss['type'] == 'G')
ring = (ss['type'] == 'R')

ax0.errorbar(L[gap], ss['rau'][gap], xerr=[eLlo[gap], eLhi[gap]], 
             yerr=ss['erau'][gap], fmt='o', color='C0', markersize=4,
             elinewidth=1)

#ax0.errorbar(L[ring], ss['rau'][ring], xerr=[eLlo[ring], eLhi[ring]],
#             yerr=ss['erau'][ring], fmt='o', color='C1', markersize=4,
#             elinewidth=1)

#ax0.plot(Ms, Lmm, 'oC0', markersize=2)



fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.19, top=0.98)
fig.savefig('substructure_locs.pdf')
fig.clf()
