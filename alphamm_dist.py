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


# set up plots
fig = plt.figure(figsize=(3.5, 2.0))
gs = gridspec.GridSpec(1, 1)
alims = [1.0, 5.0]


# extract sub-sample 
alp_ok = ( (db['FL_A'+ind] == 0) )
mu_alp = db['A'+ind][base & alp_ok]

print(len(mu_alp))

# plot
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(alims)
ax0.set_ylim([0, 1])
ax0.plot(L_B6, alp, 'o', markersize=2.0)
ax0.set_xlabel('$\\alpha$ '+alp_lbl)

fig.subplots_adjust(left=0.07, right=0.93, bottom=0.09, top=0.98)
fig.savefig('alpha_'+ind+'_dist.pdf')
fig.clf()
