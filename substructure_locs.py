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
Rlims = [0.0, 170.]
dlims = [0.0, 1.8]
Mlims = [0.0, 2.5]

# Panel (a) setups  [feature location versus Lstar]
ax0.set_xlim(Llims)
ax0.set_xscale('log')
ax0.set_xticks([0.1, 1, 10])
ax0.set_xticklabels(['0.1', '1', '10'])
ax0.set_xlabel('$L_{\\ast} \;$ (L$_\odot$)')

ax0.set_ylim(Rlims)
#ax0.set_yscale('log')
#ax0.set_yticks([10, 100])
#ax0.set_yticklabels(['10', '100'])
ax0.set_ylabel('$R_{\\rm gap} \;$ (au)')

# Panel (b) setups  [mass-size relation]
ax1.set_xlim(Mlims)
#ax1.set_xscale('log')
#ax1.set_xticks([0.1, 1])
#ax1.set_xticklabels(['0.1', '1'])
ax1.set_xlabel('$M_{\\ast} \;$ (M$_\odot$)')

ax1.set_ylim(dlims)
#ax1.set_yscale('log')
#ax1.set_yticks([10, 100])
#ax1.set_yticklabels(['10', '100'])
ax1.set_ylabel('$(R_{\\rm ring} - R_{\\rm gap})/R_{\\rm gap}$')


### Load the database

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)


# load substructure locations file
ss = ascii.read('data/substructure_locations.txt')


# grab Lstar, Mstar, and uncertainties for each substructure from database
L, eLhi, eLlo = np.zeros(len(ss)), np.zeros(len(ss)), np.zeros(len(ss))
M, eMhi, eMlo = np.zeros(len(ss)), np.zeros(len(ss)), np.zeros(len(ss))
for i in range(len(ss)):
    ind = np.where(db['NAME'] == ss['name'][i])[0][0]
    L[i] = 10.**(db['logLs'][ind])
    eLhi[i] = 10.**(db['elogLs_H'][ind]+db['logLs'][ind])-L[i]
    eLlo[i] = L[i] - 10.**(db['logLs'][ind]-db['elogLs_L'][ind])
    M[i] = 10.**(db['logMs'][ind])
    eMhi[i] = 10.**(db['elogMs_H'][ind]+db['logMs'][ind])-M[i]
    eMlo[i] = M[i] - 10.**(db['logMs'][ind]-db['elogMs_L'][ind])


# overlay some snowline models
Lgrid = np.logspace(-2, 2, 1024)
sigSB, au, Lsun = 5.67051e-5, 1.496e13, 3.826e33

Tcond = [19., 25.5, 64.5]
Tlo = [17., 19., 57.]
Thi = [21., 32., 72.]
col = ['C1', 'C2', 'C3']
for i in range(len(Tcond)):
    rcond = np.sqrt(0.02*Lgrid*Lsun / (8*np.pi*sigSB*Tcond[i]**4)) / au
    rlo = np.sqrt(0.02*Lgrid*Lsun / (8*np.pi*sigSB*Tlo[i]**4)) / au
    rhi = np.sqrt(0.02*Lgrid*Lsun / (8*np.pi*sigSB*Thi[i]**4)) / au
    ax0.fill_between(Lgrid, rlo, rhi, facecolor=col[i], alpha=0.3)
    ax0.plot(Lgrid, rcond, '--'+col[i])

# gaps
gap  = (ss['type'] == 'G')
ring = (ss['type'] == 'R')

ax0.errorbar(L[gap], ss['rau'][gap], xerr=[eLlo[gap], eLhi[gap]], 
             yerr=ss['erau'][gap], fmt='o', color='C0', markersize=3,
             elinewidth=1, alpha=0.7)

ax0.text(7., 150., 'N$_2$', color='C1', fontsize=8.5, ha='right')
ax0.text(30., 70., 'CO', color='C2', fontsize=8.5, ha='left')
ax0.text(35., 15., 'CO$_2$', color='C3', fontsize=8.5, ha='left',va='top')





# gap/ring pair separations

# models
xx  = np.linspace(1e-3, Mlims[1], 1024)
Mpl = [15.17, 95.16, 317.907, 3179.07]
col = ['C1', 'C2', 'C3', 'C4']
lbl = ['Neptune', 'Saturn', 'Jupiter', '10 Jupiter']
twk = [-0.02, 0.00, 0.02, 0.0]
for i in range(len(Mpl)):
    yy  = 4 * (Mpl[i]*5.974e27 / (xx * 1.989e33))**(1./3.)
    ax1.plot(xx, yy, '--'+col[i])
    ax1.text(2.55, yy[-1]+twk[i], lbl[i], color=col[i], fontsize=8.5,
             horizontalalignment='left', verticalalignment='center')



disks, dind = np.unique(ss['name'], return_index=True)
for i in range(len(disks)):

    # find all feature pairs in this disk 
    pairs = ss['pair'][ss['name'] == disks[i]]
    locs  = ss['rau'][ss['name'] == disks[i]]
    elocs = ss['erau'][ss['name'] == disks[i]]
    types = ss['type'][ss['name'] == disks[i]]
    Mp = (np.unique(M[ss['name'] == disks[i]]))[0]
    eMp_l = (np.unique(eMlo[ss['name'] == disks[i]]))[0]
    eMp_h = (np.unique(eMhi[ss['name'] == disks[i]]))[0]
    
    # for each pair, calculate the fractional distance
    # and plot it
    npairs = len(np.unique(pairs))
    for j in range(npairs):
        rg = locs[types == 'G'][j]
        erg = elocs[types == 'G'][j]
        rr = locs[types == 'R'][j]
        err = elocs[types == 'R'][j]
        dr = (rr - rg)/rg
        edr = np.sqrt( ( (rg/rr**2)*err )**2 + ( (1./rr)*erg )**2 )

        ax1.errorbar(Mp, dr, xerr=[[eMp_l], [eMp_h]], yerr=edr, fmt='o', 
                     color='C0', markersize=3, elinewidth=1, alpha=0.7)



ax0.text(0.08, 0.86, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.08, 0.86, 'b', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
        

fig.subplots_adjust(wspace=0.37)
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.19, top=0.98)
fig.savefig('substructure_locs.pdf')
fig.clf()
