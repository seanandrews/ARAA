import numpy as np
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)


# set up plot
fig = plt.figure(figsize=(6.33, 3.2))
gs = gridspec.GridSpec(2, 2)
amax_lims = [1.25e-2, 80]
kabs_lims = [0.05, 5]
beta_lims = [0, 3.7]
albp_lims = [-0.1, 1.1]

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(amax_lims)
ax0.set_xscale('log')
ax0.set_ylim(kabs_lims)
ax0.set_yscale('log')
ax0.set_yticks([0.1, 1])
ax0.set_yticklabels(['0.1', '1'])
ax0.set_ylabel('$\kappa_\\nu \;$ (cm$^2$ g$^{-1}$)')
ax0.set_xticklabels([])

ax1 = fig.add_subplot(gs[1, 0])
ax1.set_xlim(amax_lims)
ax1.set_xscale('log')
ax1.set_ylim(beta_lims)
ax1.set_xticks([0.1, 1, 10])
ax1.set_xticklabels(['0.1', '1', '10'])
ax1.set_xlabel('$a_{\\rm max} \;$ (mm)')
ax1.set_yticks([0, 1, 2, 3])
ax1.yaxis.get_ticklocs(minor=True)
ax1.minorticks_on()
ax1.set_ylabel('$\\beta$', labelpad=12)

ax2 = fig.add_subplot(gs[0, 1])
ax2.set_xlim(amax_lims)
ax2.set_xscale('log')
ax2.set_ylim(albp_lims)
ax2.set_xticklabels([])
ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax2.set_ylabel('$\omega_\\nu$')
ax2.yaxis.get_ticklocs(minor=True)
ax2.minorticks_on()
#ax2.plot([0.01, 100], [0.95, 0.95], ':')

ax3 = fig.add_subplot(gs[1, 1])
ax3.set_xlim(amax_lims)
ax3.set_xscale('log')
ax3.set_xticks([0.1, 1, 10])
ax3.set_xticklabels(['0.1', '1', '10'])
ax3.set_ylim(albp_lims)
ax3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax3.set_xlabel('$a_{\\rm max} \;$ (mm)')
ax3.set_ylabel('$\mathcal{P}_\\nu$')
ax3.yaxis.get_ticklocs(minor=True)
ax3.minorticks_on()



# set up loop through data
dfiles = ['q2.5_por0.5', 'q3.5_por0.5', 'q2.5_compact', 'q3.5_compact']
lstys  = ['--C0', '--C1', 'C0', 'C1']


for i in range(len(dfiles)):

    # load the data
    dat = np.load('data/dsharp_'+dfiles[i]+'.npz')
    acm, wl = 10.*dat['acm'], dat['wl']
    ka, ks, alb, pol = dat['kabs'], dat['ksca'], dat['albedo'], dat['pol']

    # mm absorption opacities
    ax0.plot(acm, ka[:, 1], lstys[i], linewidth=2)

    # absorption opacity slope (beta)
    ax1.plot(acm, np.log(ka[:, 1]/ka[:, 2]) / np.log(wl[2]/wl[1]), lstys[i],
             linewidth=2)

    # albedo
    ax2.plot(acm, alb[:, 1], lstys[i], linewidth=2)

    # polarization 
    ax3.plot(acm, np.abs(pol[:, 1]), lstys[i], linewidth=2)


ax0.text(0.09, 0.78, 'a', transform=ax0.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax1.text(0.09, 0.78, 'c', transform=ax1.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax2.text(0.09, 0.78, 'b', transform=ax2.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')
ax3.text(0.09, 0.78, 'd', transform=ax3.transAxes, horizontalalignment='left',
         fontsize=15, color='gray')

ax2.plot([2, 7], [0.34, 0.34], '-k', lw=1.5)
ax2.text(8., 0.315, 'compact', fontsize=8)
ax2.plot([2, 7], [0.23, 0.23], '--k', lw=1.5)
ax2.text(8., 0.205, 'porous', fontsize=8)
ax2.text(2, 0.07, '$q = $', fontsize=8)
ax2.text(5, 0.07, '2.5', fontsize=8, color='C0')
ax2.text(9, 0.07, ',', fontsize=8)
ax2.text(12, 0.07, '3.5', fontsize=8, color='C1')




fig.subplots_adjust(wspace=0.37, hspace=0.00)
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.13, top=0.98)
fig.savefig('opac.pdf')
fig.clf()
