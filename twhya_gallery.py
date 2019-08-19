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
xlims = [4.2017, -4.2017]
ylims = [-4.2017, 4.2017]

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
#font0 = FontProperties()
#font1 = font0.copy()
#font1.set_family('sans-serif')
#font1.set_style('italic')
#font1.set_size('small')
#ax0.text(0.96, 0.03, 'scattered light', transform=ax0.transAxes, 
#         horizontalalignment='right', fontproperties=font1, color='w',
#         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})
#ax1.text(0.96, 0.03, 'thermal continuum', transform=ax1.transAxes,
#         horizontalalignment='right', fontproperties=font1, color='w')
#ax2.text(0.96, 0.03, 'spectral line emission', transform=ax2.transAxes,
#         horizontalalignment='right', fontproperties=font1, color='w')
#font2 = font0.copy()
#font2.set_family('sans-serif')
#font2.set_size('large')
#font2.set_weight('semibold')
#ax0.text(0.08, 0.91, 'a', transform=ax0.transAxes, 
#         horizontalalignment='right', fontproperties=font2, color='w',
#         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})
#ax1.text(0.08, 0.91, 'b', transform=ax1.transAxes,
#         horizontalalignment='right', fontproperties=font2, color='w')
#ax2.text(0.08, 0.91, 'c', transform=ax2.transAxes,
#         horizontalalignment='right', fontproperties=font2, color='w')

ax0.text(0.96, 0.03, '\textit{scattered light}', transform=ax0.transAxes,
         horizontalalignment='right', color='w', 
         bbox={'facecolor': 'black', 'alpha': 0.8, 'pad': 2})




# SCATTERED LIGHT

hdulist = fits.open('data/Hband_Qr_cleaned.fits_smooth')
Iscat = hdulist[0].data[0]
hdr = hdulist[0].header
nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
cellsize = 12.25 * 1e-3    # in arcseconds (based on IRDIS plate scale)
RA, DEC = np.meshgrid(cellsize*(np.arange(nx)-0.5*nx+0.5), \
                      cellsize*(np.arange(ny)-0.5*ny+0.5))
ext = (np.max(RA), np.min(RA), np.min(DEC), np.max(DEC))
norm = ImageNormalize(vmin=10, vmax=45, stretch=LinearStretch())
im = ax0.imshow(Iscat, origin='lower', cmap='afmhot', extent=ext, 
                aspect='equal', norm=norm)
beam = Ellipse((xlims[0] + 0.08*np.diff(xlims), xlims[1] - 0.06*np.diff(xlims)),
               0.049, 0.049, 0.)
beam.set_facecolor('w')
ax0.add_artist(beam)
#ax0.plot([-5,5],[0,0],'green')
#ax0.plot([0,0],[-5,5],'green')


# MM CONTINUUM

idir = '/pool/asha0/TALKS/2017/Texas/data/'
hdulist = fits.open(idir+'B6B7_cont.image.tt0.fits')
Imm = np.squeeze(hdulist[0].data)
hdr = hdulist[0].header
beam = (np.pi/180.)**2 * np.pi * hdr['BMAJ'] * hdr['BMIN'] / (4.*np.log(2.))
nu = hdr['CRVAL3']
cc, kk = 2.9979e10, 1.381e-16
Tb = (1e-23 * Imm / beam) * cc**2 / (2.*kk*nu**2)
RA  = 3600. * hdr['CDELT1'] * (np.arange(hdr['NAXIS1'])-(hdr['CRPIX1']-1))
DEC = 3600. * hdr['CDELT2'] * (np.arange(hdr['NAXIS2'])-(hdr['CRPIX2']-1))
ext = (np.max(RA)-0.012, np.min(RA)-0.012, np.min(DEC)-0.015, np.max(DEC)-0.015)
norm = ImageNormalize(vmin=0, vmax=40, stretch=AsinhStretch())
im = ax1.imshow(Tb, origin='lower', cmap='inferno', extent=ext,
                aspect='equal', norm=norm)
beam = Ellipse((xlims[0] + 0.08*np.diff(xlims), xlims[1] - 0.06*np.diff(xlims)),
               hdr['BMAJ']*3600., hdr['BMIN']*3600., 90.-hdr['BPA'])
beam.set_facecolor('w')
ax1.add_artist(beam)
#ax1.plot([-5,5],[0,0],'green')
#ax1.plot([0,0],[-5,5],'green')


# CO

hdulist = fits.open('data/TWHya_CO_highres.pbcor.mom0.clipped.fits')
Ico = np.nan_to_num(np.squeeze(hdulist[0].data))
hdr = hdulist[0].header
RA  = 3600. * hdr['CDELT1'] * (np.arange(hdr['NAXIS1'])-(hdr['CRPIX1']-1))
DEC = 3600. * hdr['CDELT2'] * (np.arange(hdr['NAXIS2'])-(hdr['CRPIX2']-1))
ext = (np.max(RA)+0.05, np.min(RA)+0.05, np.min(DEC)-0.05, np.max(DEC)-0.05)
norm = ImageNormalize(vmin=0, vmax=0.3, stretch=AsinhStretch())
im = ax2.imshow(Ico, origin='lower', cmap='bone', extent=ext,
                aspect='equal', norm=norm)
beam = Ellipse((xlims[0] + 0.08*np.diff(xlims), xlims[1] - 0.06*np.diff(xlims)),
               hdr['BMAJ']*3600., hdr['BMIN']*3600., 90.-hdr['BPA'])
beam.set_facecolor('w')
ax2.add_artist(beam)
#ax2.plot([-5,5],[0,0],'green')
#ax2.plot([0,0],[-5,5],'green')



# adjustments for aesthetic purposes
fig.subplots_adjust(wspace=0.05)
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.01, top=0.99)
fig.savefig('twhya_gallery.pdf')
fig.clf()
