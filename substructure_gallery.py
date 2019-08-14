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
from matplotlib.font_manager import FontProperties
from matplotlib import mlab, cm
from astropy.visualization import (AsinhStretch, LinearStretch, ImageNormalize)
from matplotlib.patches import Ellipse


# set up plot
hr = 0.22
fig = plt.figure(figsize=(6.33, 7.5))
gs  = gridspec.GridSpec(10, 6, width_ratios=(1, 1, 1, 1, 1, 1), 
                        height_ratios=(hr, 1, 1, hr, 1, 1, hr, 1, hr, 1))


# gallery lists
fdir = 'data/'
im_files = 'HD169142.selfcal.concat.GPU-UVMEM.centered_mJyBeam'
offx = 0.
offy = 0.
xlims = [1., -1.]
ylims = [-1., 1.]

vmins = 0.
vmaxs = 15.

cm = 'inferno'


# set some constants
cc, kk = 2.9979e10, 1.381e-16

# loop through gallery images
for i in range(36):

    # load image and header
    hdulist = fits.open(fdir+im_files+'.fits')
    Inu = np.squeeze(hdulist[0].data)
    hdr = hdulist[0].header

    # convert to brightness temperature
    beam = (np.pi/180.)**2 * np.pi * hdr['BMAJ'] * hdr['BMIN'] / (4*np.log(2))
    nu = 232e9	# change
    Inu *= 1e-3	# change
    Tb = (1e-23 * Inu / beam) * cc**2 / (2 * kk * nu**2)

    # define coordinate grid
    RA  = 3600 * hdr['CDELT1'] * (np.arange(hdr['NAXIS1']) - (hdr['CRPIX1']-1))
    DEC = 3600 * hdr['CDELT2'] * (np.arange(hdr['NAXIS2']) - (hdr['CRPIX2']-1))
    ext = (np.max(RA)-offx, np.min(RA)-offx, 
           np.min(DEC)-offy, np.max(DEC)-offy)

    # plot the image
    #ax = fig.add_subplot(gs[np.floor_divide(i, 6), i % 6])
    if (i < 6): rowmap = 1
    if ((i >= 6) & (i < 12)): rowmap = 2
    if ((i >= 12) & (i < 18)): rowmap = 4
    if ((i >= 18) & (i < 24)): rowmap = 5
    if ((i >= 24) & (i < 30)): rowmap = 7
    if (i >= 30): rowmap = 9
    ax = fig.add_subplot(gs[rowmap, i % 6])
    norm = ImageNormalize(vmin=vmins, vmax=vmaxs, stretch=AsinhStretch())
    im = ax.imshow(Tb, origin='lower', cmap=cm, extent=ext, aspect='equal', 
                   norm=norm)

    # panel properties
    ax.set_xlim(xlims)
    ax.set_xticklabels([])
    ax.set_ylim(ylims)
    ax.set_yticklabels([])



    

# adjustments for aesthetic purposes
fig.subplots_adjust(wspace=0.04, hspace=0.04)
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.99)
fig.savefig('substructure_gallery.pdf')
fig.clf()
