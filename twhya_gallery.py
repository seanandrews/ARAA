import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from matplotlib import mlab, cm
from astropy.visualization import (AsinhStretch, LogStretch, ImageNormalize)
import sys
from matplotlib.patches import Ellipse

plt.rc('font', size=9)


# set up plot
fig = plt.figure(figsize=(7.0, 2.4))
gs  = gridspec.GridSpec(1, 3)
xlims = [3.36, -3.36]
ylims = [-3.36, 3.36]


# SCATTERED LIGHT

# load FITS image
hdulist = fits.open('data/TWHya_Hband_Qr.fits')
Iscat = hdulist[0].data[0]
hdr = hdulist[0].header

# create the image grid
nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
cellsize = 12.25 * 1e-3    # in arcseconds (based on IRDIS plate scale)
RA, DEC = np.meshgrid(cellsize*(np.arange(nx)-0.5*nx+0.5), \
                      cellsize*(np.arange(ny)-0.5*ny+0.5))
ext = (np.max(RA), np.min(RA), np.min(DEC), np.max(DEC))

# plot the image
ax = fig.add_subplot(gs[0,0])
norm = ImageNormalize(vmin=5, vmax=45, stretch=AsinhStretch())
im = ax.imshow(Iscat, origin='lower', cmap='inferno', extent=ext, 
               aspect='equal', norm=norm)

# plot beam
beam = Ellipse((xlims[0] + 0.12*np.diff(xlims), xlims[1] - 0.1*np.diff(xlims)),
               0.048, 0.048, 0.)
beam.set_facecolor('w')
ax.add_artist(beam)

# figure properties
ax.set_xlim(xlims)
ax.set_ylim(ylims)
#ax.set_xticklabels([])
#ax.set_yticklabels([])
ax.set_xlabel('$\Delta \\alpha$  $(^{\prime\prime})$')
ax.set_ylabel('$\Delta \\delta$  $(^{\prime\prime})$')



# MM CONTINUUM

# load FITS image
idir = '/pool/asha0/TALKS/2017/Texas/data/'
hdulist = fits.open(idir+'B6B7_cont.image.tt0.fits')
Imm = np.squeeze(hdulist[0].data)
hdr = hdulist[0].header

# convert to Tb
beam = (np.pi/180.)**2 * np.pi * hdr['BMAJ'] * hdr['BMIN'] / (4.*np.log(2.))
nu = hdr['CRVAL3']
cc, kk = 2.9979e10, 1.381e-16
Tb = (1e-23 * Imm / beam) * cc**2 / (2.*kk*nu**2)

# define coordinate grid
RA  = 3600. * hdr['CDELT1'] * (np.arange(hdr['NAXIS1'])-(hdr['CRPIX1']-1))
DEC = 3600. * hdr['CDELT2'] * (np.arange(hdr['NAXIS2'])-(hdr['CRPIX2']-1))
ext = (np.max(RA)-0.012, np.min(RA)-0.012, np.min(DEC)-0.015, np.max(DEC)-0.015)

# plot the image
ax = fig.add_subplot(gs[0,1])
norm = ImageNormalize(vmin=0, vmax=40, stretch=AsinhStretch())
im = ax.imshow(Tb, origin='lower', cmap='inferno', extent=ext,
               aspect='equal', norm=norm)

# plot beam
beam = Ellipse((xlims[0] + 0.12*np.diff(xlims), xlims[1] - 0.1*np.diff(xlims)),
               hdr['BMAJ']*3600., hdr['BMIN']*3600., 90.-hdr['BPA'])
beam.set_facecolor('w')
ax.add_artist(beam)

# figure properties
ax.set_xlim(xlims)
ax.set_ylim(ylims)
ax.set_xticklabels([])
ax.set_yticklabels([])



# CO

# load FITS image
hdulist = fits.open('data/TWHya_CO_highres.pbcor.mom0.clipped.fits')
Ico = np.nan_to_num(np.squeeze(hdulist[0].data))
hdr = hdulist[0].header

# convert to Tb

# define coordinate grid
RA  = 3600. * hdr['CDELT1'] * (np.arange(hdr['NAXIS1'])-(hdr['CRPIX1']-1))
DEC = 3600. * hdr['CDELT2'] * (np.arange(hdr['NAXIS2'])-(hdr['CRPIX2']-1))
ext = (np.max(RA)-0.012, np.min(RA)-0.012, np.min(DEC)-0.015, np.max(DEC)-0.015)

# plot the image
ax = fig.add_subplot(gs[0,2])
norm = ImageNormalize(vmin=0, vmax=0.3, stretch=AsinhStretch())
im = ax.imshow(Ico, origin='lower', cmap='inferno', extent=ext,
               aspect='equal', norm=norm)

# plot beam
beam = Ellipse((xlims[0] + 0.12*np.diff(xlims), xlims[1] - 0.1*np.diff(xlims)),
               hdr['BMAJ']*3600., hdr['BMIN']*3600., 90.-hdr['BPA'])
beam.set_facecolor('w')
ax.add_artist(beam)

# figure properties
ax.set_xlim(xlims)
ax.set_ylim(ylims)
ax.set_xticklabels([])
ax.set_yticklabels([])



# adjustments for aesthetic purposes
fig.subplots_adjust(hspace=0.03, wspace=0.03)
fig.subplots_adjust(left=0.07, right=0.93, bottom=0.12, top=0.99)
fig.savefig('twhya_gallery.pdf')
fig.clf()
 
