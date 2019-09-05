import numpy as np
import os
import sys
from astropy.io import ascii
from post_summary import post_summary

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# for safety, copy over database file
os.system('cp -r DISKS.csv temp2.csv')

dref = 150.

# now load database
db = ascii.read('temp2.csv', format='csv', fast_reader=True)



# base selections
base = ( (db['FL_MULT'] != 'B') & (db['FL_MULT'] != 'HJB') & 
         (db['FL_MULT'] != 'J') & (db['FL_MULT'] != 'CB') & 
         (db['FL_MULT'] != 'WJ') & (db['FL_MULT'] != 'HJ') & 
         (db['FL_MULT'] != 'HCB') & (db['FL_MULT'] != 'HC') & 
         (db['SED'] != 'I') & (db['SED'] != 'III') & (db['SED'] != 'DEBRIS') &
         (db['FL_logMs'] == 0) & (db['FL_logt'] != 9) )


### Set some constants
d_ref, nu_ref, alp = 150., 340., 2.3


# calculate luminosities, upper limits, masses, and ages + uncertainties
L7 = db['F_B7'] * (nu_ref / db['nu_B7'])**alp * (db['DPC'] / d_ref)**2
L6 = db['F_B6'] * (nu_ref / db['nu_B6'])**alp * (db['DPC'] / d_ref)**2
eL7 = np.sqrt( (nu_ref/db['nu_B7'])**(2.*alp) * \
               (db['eF_B7']**2 + (0.1*db['F_B7'])**2) * \
               (db['DPC']/d_ref)**4 + \
               ( ((nu_ref / db['nu_B7'])**alp * \
                  db['F_B7']*(2.*db['DPC']/d_ref**2) * \
                 0.5*(db['EDPC_H']+db['EDPC_L']) )**2 ) )
eL6 = np.sqrt( (nu_ref/db['nu_B6'])**(2.*alp) * \
               (db['eF_B6']**2 + (0.1*db['F_B6'])**2) * \
               (db['DPC']/d_ref)**4 + \
               ( ((nu_ref / db['nu_B6'])**alp * \
                  db['F_B6']*(2.*db['DPC']/d_ref**2) * \
                 0.5*(db['EDPC_H']+db['EDPC_L']) )**2 ) )
limL7 = db['LIM_B7'] * (nu_ref/db['nu_B7'])**alp * (db['DPC']/d_ref)**2
limL6 = db['LIM_B6'] * (nu_ref/db['nu_B6'])**alp * (db['DPC']/d_ref)**2
logM = db['logMs']
age = db['logt']
eage_hi = 10.**(db['elogt_H']+db['logt']) - age
eage_lo = age - 10.**(db['logt']-db['elogt_L'])
reg = db['SFR']
names = db['NAME']


# targets with B7 detections
detB7 = ( (db['FL_B7'] == 0) & base )
L_detB7 = L7[detB7]
eL_detB7 = eL7[detB7]
M_detB7 = logM[detB7]
age_detB7 = age[detB7] 
eage_hi_detB7 = eage_hi[detB7] / 1e6
eage_lo_detB7 = eage_lo[detB7] / 1e6
sfr_detB7 = reg[detB7]
flags_detB7 = np.zeros_like(L_detB7)
names_detB7 = names[detB7]


# targets with **only** B6 detections (i.e., no B7 or B7 limit)
detB6 = ( (db['FL_B7'] != 0) & (db['FL_B6'] == 0) & base )
L_detB6 = L6[detB6]
eL_detB6 = eL6[detB6]
M_detB6 = logM[detB6]
age_detB6 = age[detB6]
eage_hi_detB6 = eage_hi[detB6] / 1e6
eage_lo_detB6 = eage_lo[detB6] / 1e6
sfr_detB6 = reg[detB6]
flags_detB6 = np.zeros_like(L_detB6)
names_detB6 = names[detB6]


# targets with **only** limits or missing data
# (there should be **no entry** without a limit in at least B6 or B7)
lims = ( (db['FL_B7'] != 0) & (db['FL_B6'] != 0) & base )
dlims = np.ma.column_stack( (limL7[lims], limL6[lims]) )
L_lims = np.ma.min(dlims, 1)
eL_lims = np.zeros_like(L_lims)
M_lims = logM[lims]
age_lims = age[lims]
eage_hi_lims = eage_hi[lims] / 1e6
eage_lo_lims = eage_lo[lims] / 1e6
sfr_lims = reg[lims]
flags_lims = np.ones_like(L_lims)
names_lims = names[lims]


# combine everything
Lmm = np.ma.concatenate( (L_detB7, L_detB6, L_lims) )
eLmm = np.ma.concatenate( (eL_detB7, eL_detB6, eL_lims) )
Ms = np.ma.concatenate( (M_detB7, M_detB6, M_lims) )
ages = np.ma.concatenate( (age_detB7, age_detB6, age_lims) )
eages_h = np.ma.concatenate( (eage_hi_detB7, eage_hi_detB6, eage_hi_lims) )
eages_l = np.ma.concatenate( (eage_lo_detB7, eage_lo_detB6, eage_lo_lims) )
sfrs = np.ma.concatenate( (sfr_detB7, sfr_detB6, sfr_lims) )
dflags = np.ma.concatenate( (flags_detB7, flags_detB6, flags_lims) )
nnames = np.ma.concatenate( (names_detB7, names_detB6, names_lims) )


# discretized log(Mstar) bins
binM_lo = [-1.0, -0.8, -0.6, -0.4, -0.2]	#, 0.0]
binM_hi = [-0.8, -0.6, -0.4, -0.2,  0.0]	#, 0.24]


# add limits on Ms
cut = ( (Ms >= binM_lo[0]) & (Ms <= binM_hi[-1]) )
c_sfrs = sfrs[cut]
c_names = nnames[cut]

region = 'Usco'


# new cut on sfr
in_names = c_names[c_sfrs == region]

nsamples = 100000
lages = []
for i in range(len(in_names)):

    # fetch the age posterior
    iage = np.load('outputs/'+in_names[i]+'.age-mass.posterior.npz')['logAGE']

    # randomly draw nsamples from that posterior
    lages = np.append(lages, np.random.choice(iage, nsamples))


peak, mu0, mu1 = post_summary(lages, prec=0.01)
print(peak)

# plot the histogram
fig = plt.figure(figsize=(6.33, 2.4))
gs = gridspec.GridSpec(1, 1)

ax0 = fig.add_subplot(gs[0, 0])
agelims = [5.0, 7.5]
ax0.set_xlim(agelims)
ax0.set_ylim([0, 1])

N, bins, patches = ax0.hist(lages, range=agelims, bins=100, density=True,
                            color='C0', align='mid', histtype=u'step',
                            linewidth=1.5)

ax0.plot([peak, peak], [0, 1], ':k')
ax0.plot([6.18, 6.18], [0, 1], ':r')


fig.subplots_adjust(wspace=0.31, hspace=0.0)
fig.subplots_adjust(left=0.09, right=0.91, bottom=0.17, top=0.98)
fig.savefig(region+'.age_dist.pdf')
fig.clf()
