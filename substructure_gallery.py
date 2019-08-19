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
hr = 0.27
fig = plt.figure(figsize=(6.33, 7.5))
gs  = gridspec.GridSpec(10, 6, width_ratios=(1, 1, 1, 1, 1, 1), 
                        height_ratios=(0.08, 1, 1, hr, 1, 1, hr, 1, hr, 1))


# gallery lists
fdir = ['data/',	# 0,0
        'data/',	# 0,1
        'data/',	# 0,2
        'data/',	# 0,3
        'data/',	# 0,4
        'data/nienke/',	# 0,5
        'data/paola_tds/',	#1,0
        'data/']		#1,1

im_files = ['CIDA_9_selfcal_cont_image',	# 0,0
	    'sz91.cal.image',			# 0,1
            'SR24S_Band6',			# 0,2
            'HD34282_b7_continuum_selfcal_superuniform.image', #0,3
            'IP_Tau_selfcal_cont_image',	# 0,4
            'SR21_B7_tim',			# 0,5
            'J16042165',			# 1,0
            'lkca15.900.2']			# 1,1

#'HD169142.selfcal.concat.GPU-UVMEM.centered_mJyBeam'

offx = [-0.47, -0.43,  0.00, -0.05, 0.08,  0.10,  0.00, 0.245]
offy = [-0.75, -0.80,  0.00,  0.05, 0.15,  0.05,  0.00, -0.321]
xlims = np.array([ [1., -1.], [1., -1.], [1., -1.], [1., -1.], [1., -1.],
                   [1., -1.], [1., -1.], [1.3, -1.3] ])
ylims = -xlims

vmins = [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.2]
vmaxs = [ 8.0,  3.0, 15.0, 13.0,  7.0, 10.0,  7.0,  8.0]  

cm = ['inferno', 'inferno', 'inferno', 'inferno', 'inferno', 'inferno',
      'inferno', 'inferno']

dinp = ['mm', 'mm', 'mm', 'mm', 'mm', 'mm',
        'mm', 'mm']


# set some constants
cc, kk = 2.9979e10, 1.381e-16

# loop through gallery images
for i in range(len(fdir)):

    # load image and header
    hdulist = fits.open(fdir[i]+im_files[i]+'.fits')
    Inu = np.squeeze(hdulist[0].data)
    hdr = hdulist[0].header

    # mm continuum setups
    if (dinp[i] == 'mm'):
        beam = (np.pi/180)**2 * np.pi*hdr['BMAJ']*hdr['BMIN'] / (4*np.log(2))
        nu = hdr['CRVAL3']
        Tb = (1e-23 * Inu / beam) * cc**2 / (2 * kk * nu**2)

        # define coordinate grid
        RA  = 3600*hdr['CDELT1']*(np.arange(hdr['NAXIS1'])-(hdr['CRPIX1']-1))
        DEC = 3600*hdr['CDELT2']*(np.arange(hdr['NAXIS2'])-(hdr['CRPIX2']-1))

        # color stretch
        norm = ImageNormalize(vmin=vmins[i], vmax=vmaxs[i], 
                              stretch=AsinhStretch())

        # reset image name
        img = Tb

    ext = (np.max(RA)-offx[i], np.min(RA)-offx[i], 
           np.min(DEC)-offy[i], np.max(DEC)-offy[i])

    # set up grid layout
    if (i < 6): rowmap = 1
    if ((i >= 6) & (i < 12)): rowmap = 2
    if ((i >= 12) & (i < 18)): rowmap = 4
    if ((i >= 18) & (i < 24)): rowmap = 5
    if ((i >= 24) & (i < 30)): rowmap = 7
    if (i >= 30): rowmap = 9
    ax = fig.add_subplot(gs[rowmap, i % 6])

    # plot image
    im = ax.imshow(img, origin='lower', cmap=cm[i], extent=ext, aspect='equal', 
                   norm=norm)

    # panel properties
    ax.set_xlim(xlims[i])
    ax.set_ylim(ylims[i])
    ax.axis('off')


# category labels + annotations
ax = fig.add_subplot(gs[0, 0])
ax.axis('off')
ax.text(0.01, 0.80, 'Ring/Cavity', transform=ax.transAxes, fontsize=10,
        verticalalignment='center')
        
ax = fig.add_subplot(gs[3, 0])
ax.axis('off')
ax.text(0.01, 0.13, 'Rings/Gaps', transform=ax.transAxes, fontsize=10)

ax = fig.add_subplot(gs[6, 0])
ax.axis('off')
ax.text(0.01, 0.13, 'Arcs', transform=ax.transAxes, fontsize=10)

ax = fig.add_subplot(gs[8, 0])
ax.axis('off')
ax.text(0.01, 0.13, 'Spirals', transform=ax.transAxes, fontsize=10)


# adjustments for aesthetic purposes
fig.subplots_adjust(wspace=0.04, hspace=0.04)
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.99)
fig.savefig('substructure_gallery.pdf')
fig.clf()
