import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii


# plot layout
xinch = 3.5
yinch = 2.9
lbound, rbound, bbound, tbound = 0.13, 0.87, 0.11, 0.99
fig = plt.figure(figsize=(xinch, yinch))
gs = gridspec.GridSpec(1, 1)

# plotting stuff for each parameter
Llims = [-2.3, 0.6]     # log Lstar / Lsun
Tlims = [3.55, 3.45]      # log Teff / K


ax = fig.add_subplot(gs[0,0])

os.system('cp -r DISKS.csv temp2.csv')
db = ascii.read('temp2.csv', format='csv', fast_reader=True)
ndb = len(db)

ax.plot(db['logTeff'], db['logLs'], '.k')
ax.plot([3.498, 3.498],[-5, 2], 'k', alpha=0.5)

# MIST isochrones 
iso_dir = '/pool/asha1/COMPLETED/Lupus_sizes/analysis/MIST_isochrones/MIST_v1.1_vvcrit0.4_full_isos/'
iso_file = iso_dir+'MIST_v1.1_feh_p0.00_afe_p0.0_vvcrit0.4_full.iso'
iso_dat = ascii.read(iso_file)
iso_AGE = iso_dat['col2']
iso_M = iso_dat['col4']
iso_L = iso_dat['col9']
iso_T = iso_dat['col14']

#for i in range(len(np.unique(iso_AGE))): print(np.unique(iso_AGE)[i])

ax.plot(iso_T[iso_AGE == 5.7], iso_L[iso_AGE == 5.7], 'C0', alpha=1.0)
ax.plot(iso_T[iso_AGE == 6.0], iso_L[iso_AGE == 6.0], 'C1', alpha=1.0)
ax.plot(iso_T[iso_AGE == 6.300000000000001], iso_L[iso_AGE == 6.300000000000001], 'C2', alpha=1.0)


# MIST mass tracks
mt_dir = 'SP/ScottiePippen/data/MIST/'
mt_file = mt_dir+'MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS/00010M.track.eep'
mt_dat = ascii.read(mt_file)
mt_AGE = np.log10(mt_dat['col1'])
mt_M = mt_dat['col2']
mt_L = mt_dat['col7']
mt_T = mt_dat['col12']
ax.plot(mt_T[mt_AGE <= 7.5], mt_L[mt_AGE <= 7.5], 'C0', alpha=1.0)

mt_file = mt_dir+'MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS/00020M.track.eep'
mt_dat = ascii.read(mt_file)
mt_AGE = np.log10(mt_dat['col1'])
mt_M = mt_dat['col2']
mt_L = mt_dat['col7']
mt_T = mt_dat['col12']
ax.plot(mt_T[mt_AGE <= 7.5], mt_L[mt_AGE <= 7.5], 'C2', alpha=1.0)

mt_file = mt_dir+'MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS/00030M.track.eep'
mt_dat = ascii.read(mt_file)
mt_AGE = np.log10(mt_dat['col1'])
mt_M = mt_dat['col2']
mt_L = mt_dat['col7']
mt_T = mt_dat['col12']
ax.plot(mt_T[mt_AGE <= 7.5], mt_L[mt_AGE <= 7.5], 'C3', alpha=1.0)


# Baraffe models
mt_dir = 'SP/ScottiePippen/data/Baraffe15/'
#dat = ascii.read(mt_dir+'BHAC15_iso.dat')
agestr = np.array(['0.0005', '0.0010', '0.0020', '0.0030', '0.0040', '0.0050', '0.0080'])
t10, l10 = np.zeros(len(agestr)), np.zeros(len(agestr))
t20, l20 = np.zeros(len(agestr)), np.zeros(len(agestr))
t30, l30 = np.zeros(len(agestr)), np.zeros(len(agestr))
for i in range(len(agestr)):
    inp = ascii.read(mt_dir+agestr[i]+'.dat', comment='!')
    t10[i] = np.log10(inp['col2'][inp['col1'] == 0.05])
    l10[i] = inp['col3'][inp['col1'] == 0.05]
    t20[i] = np.log10(inp['col2'][inp['col1'] == 0.2])
    l20[i] = inp['col3'][inp['col1'] == 0.2]
    t30[i] = np.log10(inp['col2'][inp['col1'] == 0.3])
    l30[i] = inp['col3'][inp['col1'] == 0.3]

ax.plot(t10, l10, '--C0', alpha=1.0)
ax.plot(t20, l20, '--C2', alpha=1.0)
ax.plot(t30, l30, '--C3', alpha=1.0)

dat = ascii.read(mt_dir+'0.0005.dat', comment='!')
ax.plot(np.log10(dat['col2']), dat['col3'], ':C0')
dat = ascii.read(mt_dir+'0.0010.dat', comment='!')
ax.plot(np.log10(dat['col2']), dat['col3'], ':C1')
dat = ascii.read(mt_dir+'0.0020.dat', comment='!')
ax.plot(np.log10(dat['col2']), dat['col3'], ':C2')





# plot appearances
ax.set_ylim(Llims)
ax.set_xlim(Tlims)
ax.set_xlabel(r'${\rm log}$'+' '+r'$T_{\rm eff}$'+' / '+r'${\rm K}$',
              fontsize=8)
ax.set_ylabel(r'${\rm log}$'+' '+r'$L_\ast$'+' / '+r'$L_\odot$', fontsize=8)
ax.yaxis.labelpad = 2

ax.tick_params(axis='both', which='major', labelsize=7)
fig.subplots_adjust(left=lbound, right=rbound, bottom=bbound, top=tbound)
fig.savefig('HR_overlap.pdf')
fig.clf()
