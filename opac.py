import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import os
plt.rc('font', size=9)

# load and re-package opacity data
d3 = np.load('data/default_q3.5.npz')
acm = d3['acm']
ka_d3, ks_d3, b_d3, w_d3 = d3['kabs0'], d3['ksca0'], d3['beta'], d3['albedo']

d2 = np.load('data/default_q2.5.npz')
ka_d2, ks_d2, b_d2, w_d2 = d2['kabs0'], d2['ksca0'], d2['beta'], d2['albedo']

p3 = np.load('data/porous_q3.5.npz')
ka_p3, ks_p3, b_p3, w_p3 = p3['kabs0'], p3['ksca0'], p3['beta'], p3['albedo']

p2 = np.load('data/porous_q2.5.npz')
ka_p2, ks_p2, b_p2, w_p2 = p2['kabs0'], p2['ksca0'], p2['beta'], p2['albedo']

kabs = np.vstack((ka_d3, ka_d2, ka_p3, ka_p2))
ksca = np.vstack((ks_d3, ks_d2, ks_p3, ks_p2))
bet  = np.vstack((b_d3, b_d2, b_p3, b_p2))
alb  = np.vstack((w_d3, w_d2, w_p3, w_p2))


fig = plt.figure(figsize=(7., 4.2))
gs = gridspec.GridSpec(3, 2)
alims = [1e-4, 100]


# absorption opacities
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(alims)
ax0.set_xscale('log')
ax0.set_ylim([0.02, 5])
ax0.set_yscale('log')
ax0.plot(acm, ka_d3, 'C0')
ax0.plot(acm, ka_p3, '--C0')
ax0.plot(acm, ka_d2, 'C1')
ax0.plot(acm, ka_p2, '--C1')
ax0.set_ylabel('$\kappa_0^{\\rm abs}$  (cm$^2$ g$^{-1}$)')
ax0.set_xticklabels([])
#ax0.set_xlabel('$a_{max}$ (cm)')

# beta
ax1 = fig.add_subplot(gs[1, 0])
ax1.set_xlim(alims)
ax1.set_xscale('log')
ax1.set_ylim([-0.3, 3.9])
ax1.plot(acm, b_d3, 'C0')
ax1.plot(acm, b_p3, '--C0')
ax1.plot(acm, b_d2, 'C1')
ax1.plot(acm, b_p2, '--C1')
ax1.set_xticklabels([])
ax1.set_yticks([0, 1, 2, 3])
ax1.set_ylabel('$\\beta$')

# albedo
ax2 = fig.add_subplot(gs[2, 0])
ax2.set_xlim(alims)
ax2.set_xscale('log')
ax2.set_ylim([-0.1, 1.1])
ax2.plot(acm, w_d3, 'C0')
ax2.plot(acm, w_p3, '--C0')
ax2.plot(acm, w_d2, 'C1')
ax2.plot(acm, w_p2, '--C1')
ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax2.set_xlabel('$a_{max}$ (cm)')
ax2.set_ylabel('$\omega_0$')




fig.subplots_adjust(wspace=0.5, hspace=0.03)
fig.subplots_adjust(left=0.08, right=0.92, bottom=0.09, top=0.99)
fig.savefig('opac.pdf')
fig.clf()
