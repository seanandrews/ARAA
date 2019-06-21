import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii

reg = 'Cha'


# plot layout
xinch = 5.5
yinch = 3.4
lbound, rbound, bbound, tbound = 0.07, 0.93, 0.10, 0.99
fig = plt.figure(figsize=(xinch, yinch))
gs = gridspec.GridSpec(1, 1)

# plotting stuff for each parameter
Llims = [-2.9, 2.9]     # log Lstar / Lsun
Tlims = [4.0, 3.4]      # log Teff / K


ax = fig.add_subplot(gs[0,0])

# MIST --> Baraffe boundary
ax.plot([3.498, 3.498],[-5, 5], ':k', alpha=0.2)

# MIST isochrones 
iso_dir = '/pool/asha1/COMPLETED/Lupus_sizes/analysis/MIST_isochrones/MIST_v1.1_vvcrit0.4_full_isos/'
iso_dat = ascii.read(iso_dir+'MIST_v1.1_feh_p0.00_afe_p0.0_vvcrit0.4_full.iso')
iso_AGE, iso_L, iso_T = iso_dat['col2'], iso_dat['col9'], iso_dat['col14']

#for i in range(len(np.unique(iso_AGE))): print(np.unique(iso_AGE)[i])

# MIST isochrones @ 0.5, 1, 2, 4, 8, 15 Myr
age = [5.7, 6.0, 6.300000000000001, 6.6000000000000005, 6.9, 7.15]
for i in range(len(age)):
    ax.plot(iso_T[iso_AGE == age[i]], iso_L[iso_AGE == age[i]], 'C0', alpha=0.5)

# MIST mass tracks @ 0.2, 0.4, 0.8, 1.6, 3.2 Msun
ms = ['00020', '00040', '00080', '00160', '00320']
mdir = 'SP/ScottiePippen/data/MIST/MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS/'
for i in range(len(ms)):
    mt_dat = ascii.read(mdir+ms[i]+'M.track.eep')
    mt_AGE, mt_L, mt_T = mt_dat['col1'], mt_dat['col7'], mt_dat['col12']
    ax.plot(mt_T[mt_AGE <= 5e7], mt_L[mt_AGE <= 5e7], 'r', alpha=0.5)


# Baraffe isochrones @ 0.5, 1, 2, 4, 8, 15 Myr
mdir = 'SP/ScottiePippen/data/Baraffe15/'
age = ['0.0005', '0.0010', '0.0020', '0.0040', '0.0080', '0.0150']
for i in range(len(age)):
    dat = ascii.read(mdir+age[i]+'.dat', comment='!')
    logT, logL = np.log10(dat['col2']), dat['col3']
    ax.plot(logT[logT <= 3.498], logL[logT <= 3.498], 'C0', alpha=0.5)
    
# Baraffe mass tracks @ 0.02, 0.05, 0.10 Msun
age = ['0.0005', '0.0010', '0.0020', '0.0030', '0.0040', '0.0050', '0.0080', 
       '0.0100', '0.0150', '0.0200', '0.0250', '0.0300', '0.0400', '0.0500']
t10, l10 = np.zeros(len(age)), np.zeros(len(age))
t05, l05 = np.zeros(len(age)), np.zeros(len(age))
t02, l02 = np.zeros(len(age)), np.zeros(len(age))
for i in range(len(age)):
    inp = ascii.read(mdir+age[i]+'.dat', comment='!')
    t10[i] = np.log10(inp['col2'][inp['col1'] == 0.10])
    l10[i] = inp['col3'][inp['col1'] == 0.10]
    t05[i] = np.log10(inp['col2'][inp['col1'] == 0.05])
    l05[i] = inp['col3'][inp['col1'] == 0.05]
    t02[i] = np.log10(inp['col2'][inp['col1'] == 0.02])
    l02[i] = inp['col3'][inp['col1'] == 0.02]
ax.plot(t10, l10, 'r', alpha=0.5)
ax.plot(t05, l05, 'r', alpha=0.5)
ax.plot(t02, l02, 'r', alpha=0.5)


# datapoints
os.system('cp -r DISKS.csv temp2.csv')
db = ascii.read('temp2.csv', format='csv', fast_reader=True)
base = ((db['SFR'] == reg))
print(len(db), len(db[base]))

ax.errorbar(db['logTeff'][base], db['logLs'][base], xerr=db['elogTeff'][base], 
            yerr=[db['elogLs_L'][base], db['elogLs_H'][base]], fmt='.k',
            ecolor='gray', alpha=0.3)
ax.plot(db['logTeff'][base], db['logLs'][base], '.k')
            


# plot appearances
ax.set_ylim(Llims)
ax.set_xlim(Tlims)
ax.set_xlabel(r'${\rm log}$'+' '+r'$T_{\rm eff}$'+' / '+r'${\rm K}$',
              fontsize=8)
ax.set_ylabel(r'${\rm log}$'+' '+r'$L_\ast$'+' / '+r'$L_\odot$', fontsize=8)
ax.yaxis.labelpad = 2

ax.tick_params(axis='both', which='major', labelsize=7)
fig.subplots_adjust(left=lbound, right=rbound, bottom=bbound, top=tbound)
fig.savefig('HR_'+reg+'.pdf')
fig.clf()
