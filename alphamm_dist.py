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
fig = plt.figure(figsize=(6.33, 2.4))
gs = gridspec.GridSpec(1, 2)
Alims = [0., 4.]
Mlims = [0.00125, 80.]

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(Alims)
ax0.set_ylim([0, 1])
ax0.set_xlabel('$\\alpha_{\\rm mm}$')
ax0.set_ylabel('frac')


# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'T') & \
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & \
         (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') ) 


# B6/B7 indices
ok67 = ((db['FL_B7'] == 0) & (db['FL_B6'] == 0) & (db['FL_A67'] == 0) & base)
names = db['NAME'][ok67]
num67 = len(names)
a67_samples = []
for i in range(num67):
    a67_ = np.load('outputs/'+names[i]+'.alpha67.posterior.npz')['amm']
    a67_samples = np.append(a67_samples, a67_)


ok36 = ((db['FL_B6'] == 0) & (db['FL_B3'] == 0) & (db['FL_A36'] == 0) & base)
names = db['NAME'][ok36]
num36 = len(names)
a36_samples = []
for i in range(num36):
    a36_ = np.load('outputs/'+names[i]+'.alpha36.posterior.npz')['amm']
    a36_samples = np.append(a36_samples, a36_)
print(num36)


# plot the histograms
N, bins, patches = ax0.hist(a67_samples, range=[0, 4], bins=100, density=True,
                            color='C0', align='mid', histtype=u'step', 
                            linewidth=1.5)
xx = np.linspace(0, 4, 1024)
yy = 0.6*np.exp(-0.5*((xx-2.3)/0.6)**2)
ax0.plot(xx, yy, 'r')
#ax0.plot([2.3, 2.3], [0, 1], '--k')

N, bins, patches = ax0.hist(a36_samples, range=[0, 4], bins=100, density=True,
                            color='C1', align='mid', histtype=u'step', 
                            linewidth=1.5)
#ax0.plot([2.5, 2.5], [0, 1], '--k')


fig.subplots_adjust(wspace=0.31, hspace=0.0)
fig.subplots_adjust(left=0.09, right=0.91, bottom=0.17, top=0.98)
fig.savefig('alphamm_dist.pdf')
fig.clf()
