import numpy as np
import os
import sys
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from km_estimator import km_estimator

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

from matplotlib import rc
rc("font", **{"sans-serif": ["Roboto"]})
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)
rc("mathtext", default="regular")

rc("axes", linewidth=2)
rc("xtick.major", width=2)
rc("xtick.minor", width=1.5)
rc("xtick.major", size=4.5)
rc("xtick.minor", size=2.6)
rc("ytick.major", width=2)
rc("ytick.minor", width=1.5)
rc("ytick.major", size=4.5)
rc("ytick.minor", size=2.6)

rc("axes", labelsize=12)
rc("xtick", labelsize=12)
rc("ytick", labelsize=12)



# set up plot
fig = plt.figure(figsize=(6.33, 2.6))
gs = gridspec.GridSpec(1, 1)
Mlims = [0.05, 20000.]		# earth masses

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Mlims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10, 100, 1000, 10000])
ax0.set_xticklabels(['0.1', '1', '10', '100', '1000', '10$^4$'])
ax0.set_ylim([0, 1])
ax0.set_ylabel('$p$ ($\ge M$)')
ax0.set_xlabel('$M \;$ (M$_\oplus$)')


# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)
print(len(db))

# base selections
base = ( (db['FL_MULT'] != 'J') & (db['SED'] != 'III') & 
         (db['SFR'] != 'Oph') )
print(len(db[base]))

# reference distance
dref = 150.


# DUST DISK MASSES
# B7 data
have_b7 = ( (db['FL_B7'] != -1) & base )
flag_b7 = (db['FL_B7'][have_b7] == 1)
Fb7 = 1e-3 * db['F_B7'][have_b7]
limb7 = 1e-3 * db['LIM_B7'][have_b7]
Fb7[flag_b7] = limb7[flag_b7]
Lb7 = Fb7 * (db['DPC'][have_b7] / dref)**2
names = db['NAME'][have_b7]
print(len(Lb7))

# convert to dust masses
h, c, k = 6.626e-27, 2.9979e10, 1.381e-16
mean_kap7 = 3.5
mean_T = 20.
mean_nu = 340e9
Bnu = (2.*h*mean_nu**3/c**2) / (np.exp(h*mean_nu/(k*mean_T)) + 1)
mdust_simple = (dref*3.09e18)**2 * 1e-23 * Lb7 / (mean_kap7 * Bnu) 	# in g
mdust_simple /= 5.974e27

# cumulative distribution
Mdust, pMdust, epMdust, mukm = km_estimator(mdust_simple, flag_b7)


# GAS DISK MASSES
have_Mg = ( (db['FL_Mgas'] == 0) & base )
flag_Mg = (db['FL_Mgas'][have_Mg] == 1)
Mg = db['Mgas'][have_Mg] * 1.898e30 / 5.974e27

# cumulative distribution
Mgas, pMgas, epMgas, mukm = km_estimator(Mg, flag_Mg)


# plots
ax0.fill_between(Mgas, pMgas+epMgas, pMgas-epMgas,
                 facecolor='C1', alpha=0.5, step='post')
ax0.plot(Mgas, pMgas, 'C1', drawstyle='steps-post')

ax0.fill_between(Mdust, pMdust+epMdust, pMdust-epMdust, 
                 facecolor='C0', alpha=0.5, step='post')
ax0.plot(Mdust, pMdust, 'C0', drawstyle='steps-post')



fig.subplots_adjust(left=0.2, right=0.8, bottom=0.16, top=0.98)
fig.savefig('mdisk_dist.pdf')
fig.clf()

