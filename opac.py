import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import os
plt.rc('font', size=9)


# set up the figure
fig = plt.figure(figsize=(6.33, 3.3))
gs = gridspec.GridSpec(2, 2)
amax_lims = [1e-4, 100]
kabs_lims = [0.02, 5]
beta_lims = [0, 4]
albp_lims = [-0.1, 1.1]

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(amax_lims)
ax0.set_xscale('log')
ax0.set_ylim(kabs_lims)
ax0.set_yscale('log')
ax0.set_yticks([0.1, 1])
ax0.set_yticklabels(['0.1', '1'])
ax0.set_ylabel('$\kappa_\\nu$  (cm$^2$ g$^{-1}$)')
ax0.set_xticklabels([])

ax1 = fig.add_subplot(gs[1, 0])
ax1.set_xlim(amax_lims)
ax1.set_xscale('log')
ax1.set_ylim(beta_lims)
ax1.set_xlabel('$a_{\\rm max}$ (cm)')
ax1.set_yticks([0, 1, 2, 3])
ax1.set_ylabel('$\\beta$')

ax2 = fig.add_subplot(gs[0, 1])
ax2.set_xlim(amax_lims)
ax2.set_xscale('log')
ax2.set_ylim(albp_lims)
ax2.set_xticklabels([])
ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax2.set_xlabel('$a_{max}$ (cm)')
ax2.set_ylabel('$\omega_\\nu$')

ax3 = fig.add_subplot(gs[1, 1])
ax3.set_xlim(amax_lims)
ax3.set_xscale('log')
ax3.set_ylim(albp_lims)
ax3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax3.set_xlabel('$a_{\\rm max}$ (cm)')
ax3.set_ylabel('$\mathcal{P}_\\nu$')



# set up loop through data
dfiles = ['q2.5_por0.5', 'q3.5_por0.5', 'q2.5_compact', 'q3.5_compact']
lstys  = ['--C0', '--C3', 'C0', 'C3']


for i in range(len(dfiles)):

    # load the data
    dat = np.load('data/dsharp_'+dfiles[i]+'.npz')
    acm, wl = dat['acm'], dat['wl']
    ka, ks, alb, pol = dat['kabs'], dat['ksca'], dat['albedo'], dat['pol']

    # mm absorption opacities
    ax0.plot(acm, ka[:, 1], lstys[i])

    # absorption opacity slope (beta)
    ax1.plot(acm, np.log(ka[:, 1]/ka[:, 2]) / np.log(wl[2]/wl[1]), lstys[i])

    # albedo
    ax2.plot(acm, alb[:, 1], lstys[i])

    # polarization 
    ax3.plot(acm, pol[:, 1], lstys[i])



fig.subplots_adjust(wspace=0.3, hspace=0.03)
fig.subplots_adjust(left=0.08, right=0.92, bottom=0.12, top=0.99)
fig.savefig('opac.pdf')
fig.clf()
