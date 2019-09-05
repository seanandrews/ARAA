import numpy as np
import os
import sys
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('araa')
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex = True)
from matplotlib.colorbar import Colorbar
from matplotlib import mlab, cm
from astropy.visualization import (AsinhStretch, LinearStretch, ImageNormalize)
from matplotlib.patches import Ellipse
from matplotlib.font_manager import FontProperties


# set up plots
fig = plt.figure(figsize=(6.33, 2.1))
gs  = gridspec.GridSpec(1, 3)
xlims = [1., -1.]
ylims = [-1., 1.]

ax0 = fig.add_subplot(gs[0,0])
ax0.set_xlim(xlims)
ax0.set_ylim(ylims)
ax0.set_xticklabels([])
ax0.set_yticklabels([])
ax1 = fig.add_subplot(gs[0,1])
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax2 = fig.add_subplot(gs[0,2])
ax2.set_xlim(xlims)
ax2.set_ylim(ylims)
ax2.set_xticklabels([])
ax2.set_yticklabels([])


# annotations
ax0.text(0.08, 0.91, 'a', transform=ax0.transAxes, 
         horizontalalignment='right', color='w', fontsize=15, 
         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})

ax1.text(0.08, 0.91, 'b', transform=ax1.transAxes,
         horizontalalignment='right', color='w', fontsize=15,
         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})

ax2.text(0.08, 0.91, 'c', transform=ax2.transAxes,
         horizontalalignment='right', color='w', fontsize=15,
         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})


#ax0.text(0.96, 0.03, '\\textit{scattered light}', transform=ax0.transAxes,
#         horizontalalignment='right', color='w', fontsize=11,
#         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})
#ax1.text(0.96, 0.03, '\\textit{thermal continuum}', transform=ax1.transAxes,
#         horizontalalignment='right', color='w',fontsize=11,
#         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})
#ax2.text(0.96, 0.03, '\\textit{spectral line emission}', 
#         transform=ax2.transAxes,
#         horizontalalignment='right', color='w',fontsize=11,
#         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})


xlims = [0.7, -0.7]
ylims = [-0.7, 0.7]
ax0 = fig.add_subplot(gs[0,0])
ax0.set_xlim(xlims)
ax0.set_ylim(ylims)
ax0.set_xticklabels([])
ax0.set_yticklabels([])
ax0.axis('off')

hdulist = fits.open('data/Fits_for_Sean/IndividualSPHEREpapers/HD143006_Qphi_Jband.fits')
Iscat = hdulist[0].data
hdr = hdulist[0].header
nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
cellsize = 12.25 * 1e-3    # in arcseconds (based on IRDIS plate scale)
RA, DEC = np.meshgrid(cellsize*(np.arange(nx)-0.5*nx+0.5), \
                      cellsize*(np.arange(ny)-0.5*ny+0.5))
ext = (np.max(RA), np.min(RA), np.min(DEC), np.max(DEC))
norm = ImageNormalize(vmin=0, vmax=85, stretch=LinearStretch())
im = ax0.imshow(Iscat, origin='lower', cmap='afmhot', extent=ext, 
                aspect='equal', norm=norm)
beam = Ellipse((xlims[0] + 0.08*np.diff(xlims), xlims[1] - 0.06*np.diff(xlims)),
               0.049, 0.049, 0.)
beam.set_facecolor('w')
ax0.add_artist(beam)


xlims = [0.7, -0.7]
ylims = [-0.7, 0.7]
ax1 = fig.add_subplot(gs[0,1])
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.axis('off')
hdulist = fits.open('data/Fits_for_Sean/WRAY.fits')
Iscat = hdulist[0].data
hdr = hdulist[0].header
nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
cellsize = 12.25 * 1e-3    # in arcseconds (based on IRDIS plate scale)
RA, DEC = np.meshgrid(cellsize*(np.arange(nx)-0.5*nx+0.5), \
                      cellsize*(np.arange(ny)-0.5*ny+0.5))
ext = (np.max(RA), np.min(RA), np.min(DEC), np.max(DEC))
norm = ImageNormalize(vmin=0, vmax=0.001, stretch=LinearStretch())
im = ax1.imshow(Iscat, origin='lower', cmap='afmhot', extent=ext, 
                aspect='equal', norm=norm)
beam = Ellipse((xlims[0] + 0.08*np.diff(xlims), xlims[1] - 0.06*np.diff(xlims)),
               0.049, 0.049, 0.)
beam.set_facecolor('w')
ax0.add_artist(beam)



xlims = [0.7, -0.7]
ylims = [-0.7, 0.7]
ax2 = fig.add_subplot(gs[0,2])
ax2.set_xlim(xlims)
ax2.set_ylim(ylims)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.axis('off')
hdulist = fits.open('data/Fits_for_Sean/IndividualSPHEREpapers/J1604_Epoch2017_08_22_Qphi_Jband.fits')
Iscat = hdulist[0].data
hdr = hdulist[0].header
nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
cellsize = 12.25 * 1e-3    # in arcseconds (based on IRDIS plate scale)
RA, DEC = np.meshgrid(cellsize*(np.arange(nx)-0.5*nx+0.5), \
                      cellsize*(np.arange(ny)-0.5*ny+0.5))
ext = (np.max(RA), np.min(RA), np.min(DEC), np.max(DEC))
norm = ImageNormalize(vmin=0, vmax=250, stretch=LinearStretch())
im = ax2.imshow(Iscat, origin='lower', cmap='afmhot', extent=ext,
                aspect='equal', norm=norm)
beam = Ellipse((xlims[0] + 0.08*np.diff(xlims), xlims[1] - 0.06*np.diff(xlims)),
               0.049, 0.049, 0.)
beam.set_facecolor('w')
ax0.add_artist(beam)



# adjustments for aesthetic purposes
fig.subplots_adjust(wspace=0.05)
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.01, top=0.99)
fig.savefig('shadow_gallery.pdf')
fig.clf()
