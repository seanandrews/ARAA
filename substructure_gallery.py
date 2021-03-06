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
        'data/',		#1,1
        'data/paola_tds/',	#1,2
        'data/for_sean/',	#1,3
        'data/Fits_for_Sean/',	#1,4
        'data/Fits_for_Sean/GPIData/',	#1,5
        'data/',	# 2,0
        'data/',	# 2,1
        'data/',	# 2,2
        'data/',	# 2,3
        'data/',	# 2,4
        'data/',	# 2,5
        'data/',	# 3,0
        'data/',	# 3,1
        'data/',	# 3,2
        'data/Fits_for_Sean/DARTTS_Sean/',		# 3,3
        'data/Fits_for_Sean/DARTTS_Sean/',		# 3,4
        'data/Fits_for_Sean/DARTTS_Sean/',		# 3,5
        'data/',	# 4,0
        'data/',	# 4,1
        'data/',	# 4,2
        'data/',	# 4,3
        'data/',	# 4,4
        'data/',	# 4,5
        'data/',	# 5,0
        'data/',	# 5,1
        'data/',	# 5,2
        'data/Fits_for_Sean/IndividualSPHEREpapers/', # 5,3
        'data/Fits_for_Sean/IndividualSPHEREpapers/',	# 5,4
	'data/Fits_for_Sean/']

im_files = ['CIDA_9_selfcal_cont_image',	# 0,0
	    'sz91.cal.image',			# 0,1
            'SR24S_Band6',			# 0,2
            'HD34282_b7_continuum_selfcal_superuniform.image', #0,3
            'IP_Tau_selfcal_cont_image',	# 0,4
            'SR21_B7_tim',			# 0,5
            'J16042165',			# 1,0
            'dmtau.clean_p0.pbcor',		# 1,1
            'DoAr44',		# 1,2
            'IRS48_Hband_PI_halosub',		# 1,3
            'LkCa_15_2015-12-19_Q_phi_star_pol_subtr',		# 1,4
            'Monnier2017_FIG2_GPI_HD_163296_J_band_Stokes_Q_r', # 1,5
            'AS209_continuum',					# 2,0
            'HLTau_B6B7cont_mscale_n2_ap.image.tt0',		# 2,1
            'V1094_B6_v2_fin',					# 2,2
            'DL_Tau_selfcal_cont_image',
            'HD169142.selfcal.concat.GPU-UVMEM.centered_mJyBeam',
            'RULup_continuum',
            'GO_Tau_selfcal_cont_image',
            'Elias24_continuum',
            'HD163296_continuum',
            'J1615_Qphi_Hband',
            'V4046Sgr_Qphi_Hband',
            'PDS66_Qphi_Hband',
            'mwc758_20180301_SL.img',
            'hd135344b.band4.data.selfcal.center',
            'HD143006_continuum',
            'V1247Ori',
            'IRS48_B9_667',
            'uid___A002_X41c27b_X2fc__hd142527_cont_spw0',
            'IMLup_continuum',
            'WaOph6_continuum',
            'Elias27_continuum',
            'SAO206462_Qphi_r2scaled_Jband',
            'MWC758_Qphi_Yband_March2015',
            'HD100456']

dpc = [171.2, 157.8, 138.1, 311.6, 130.1, 137.8,
       149.4, 144.5, 145.3, 134.5, 158.3, 101.3,
       120.6, 135.8, 153.2, 158.2, 114.0, 158.9,
       143.6, 136.0, 101.3, 157.0,  73.7,  98.9,
       159.5, 135.4, 165.3, 398.4, 134.5, 157.3,
       157.8, 123.2, 116.4, 135.4, 159.5, 104.2]

offx = [-0.51, -0.45,  0.00, -0.07,  0.05,  0.07,  
         0.00, -0.01, -0.02,  0.00,  0.00,  0.00,
         0.01, -0.01, -0.15,  0.22,  0.00, -0.04,
        -0.18,  0.10,  0.00,  0.00,  0.00,  0.00,
         0.00,  0.20,  0.00,  0.02,  0.00, -0.20,
         0.00, -0.25,  0.00,  0.00,  0.00,  0.00]

offy = [-0.71, -0.80,  0.00,  0.05,  0.17,  0.05,  
         0.03, -0.03,  0.00,  0.20,  0.00,  0.00,
         0.00,  0.19,  0.08, -0.06,  0.00,  0.09,
        -0.42, -0.38,  0.00,  0.00,  0.00,  0.00,
         0.00, -0.05,  0.02, -0.05,  0.00, -0.20,
         0.00, -0.36,  0.00,  0.00,  0.00,  0.00]

xhw = [0.65, 1.10, 0.90, 1.00, 0.55, 0.95,
       1.20, 0.50, 1.00, 1.50, 0.70, 0.90,
       1.30, 1.05, 2.00, 1.20, 1.00, 0.50,
       1.15, 1.25, 1.20, 1.80, 0.65, 1.40,
       0.85, 1.00, 0.75, 0.75, 1.00, 2.70,
       1.30, 1.10, 1.90, 1.40, 0.70, 0.40]

vmins = [ 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  
          0.0,  1.0,  0.0,   0.0,  0.0,  0.0,
          0.0,  0.0,  0.0,   0.0,  0.0,  0.0,
          0.0,  0.0,  0.0,   0.0,  0.0,  0.0,
          0.0,  0.0,  0.0,   0.5,  0.0,  0.5,
          0.3,  0.25,  0.1,   0.0,  5.0,  0.0]

vmaxs = [ 8.0,  3.0, 15.0,  13.0,  7.0, 12.0,  
          4.0, 20.0, 11.0,   2.0,   60,  1.5,
         25.0, 75.0,  5.0,  13.0, 10.0, 75.0,
          7.0, 30.0,  60.,  1.5e-7, 2.0e-6, 2e-7,
         17.0,  9.0, 13.0,  19.0, 22.0, 20.0,
         10.0, 12.0, 10.0,  1e6,  90, 350]  

dinp = ['mm', 'mm', 'mm', 'mm', 'mm', 'mm',
        'mm', 'mm', 'mm', 'scat', 'scat', 'scat',
        'mm', 'mm', 'mm', 'mm', 'mm', 'mm',
        'mm', 'mm', 'mm', 'scat', 'scat', 'scat',
        'mm', 'mm', 'mm', 'mm', 'mm', 'mm',
        'mm', 'mm', 'mm', 'scat', 'scat', 'scat']

csize = [0, 0, 0, 0, 0, 0,
         0, 0, 9.5, 9.5, 12.25, 14.00,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 12.25, 12.25, 12.25,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 12.25, 12.25, 7.0]

bfwhm = [0, 0, 0, 0, 0, 0,
         0, 0, 0.05, 0.05, 0.049, 0.042,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 0.049, 0.049, 0.049,
         0, 0, 0, 0, 0, 0,
         0, 0, 0, 0.041, 0.027, 0.020]


# set some constants
cc, kk = 2.9979e10, 1.381e-16

# loop through gallery images
for i in range(len(fdir)):

    xlims = [xhw[i], -xhw[i]]
    ylims = [-xhw[i], xhw[i]]

    # mm continuum setups
    if (dinp[i] == 'mm'):
        hdulist = fits.open(fdir[i]+im_files[i]+'.fits')
        Inu = np.squeeze(hdulist[0].data)
        hdr = hdulist[0].header

        print(im_files[i])
        if ((i == 16) or (i == 25)):
            if (i == 16):
                nu = 232e9
                Inu *= 1e-3
            if (i == 25):
                nu = 180e9
        else: nu = hdr['CRVAL3']
        

        beam = (np.pi/180)**2 * np.pi*hdr['BMAJ']*hdr['BMIN'] / (4*np.log(2))
        Tb = (1e-23 * Inu / beam) * cc**2 / (2 * kk * nu**2)

        # define coordinate grid
        RA  = 3600*hdr['CDELT1']*(np.arange(hdr['NAXIS1'])-(hdr['CRPIX1']-1))
        DEC = 3600*hdr['CDELT2']*(np.arange(hdr['NAXIS2'])-(hdr['CRPIX2']-1))

        # color stretch and map
        norm = ImageNormalize(vmin=vmins[i], vmax=vmaxs[i], 
                              stretch=AsinhStretch())
        cm = 'inferno'

        # reset image name
        img = Tb

        # beam
        beam = Ellipse((xlims[0] + 0.09*np.diff(xlims), 
                        xlims[1] - 0.09*np.diff(xlims)),
                        hdr['BMAJ']*3600., hdr['BMIN']*3600., 90.-hdr['BPA'])


    if (dinp[i] == 'scat'):
        hdulist = fits.open(fdir[i]+im_files[i]+'.fits')
        Inu = np.squeeze(hdulist[0].data)
        hdr = hdulist[0].header
 
        # define coordinate grid
        nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
        cellsize = csize[i] * 1e-3
        RA, DEC = np.meshgrid(cellsize*(np.arange(nx)-0.5*nx+0.5), \
                              cellsize*(np.arange(ny)-0.5*ny+0.5))

        # color stretch and map
        norm = ImageNormalize(vmin=vmins[i], vmax=vmaxs[i], 
                              stretch=LinearStretch())
        cm = 'afmhot'

        # reset image name
        img = Inu

        # beam
        beam = Ellipse((xlims[0] + 0.09*np.diff(xlims),
                        xlims[1] - 0.09*np.diff(xlims)),
                        bfwhm[i], bfwhm[i], 0)
        

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
    im = ax.imshow(img, origin='lower', cmap=cm, extent=ext, aspect='equal', 
                   norm=norm)

    beam.set_facecolor('w')
    ax.add_artist(beam)

    #ax.plot([0, 0], [-10, 10], 'g', lw=1)
    #ax.plot([-10, 10], [0, 0], 'g', lw=1)

    # scale bar
    barAU = 10.
    xbarr = xlims[1] - 0.09*np.diff(xlims)
    xbarl = xbarr + barAU/dpc[i]
    ybar  = xlims[1] - 0.09*np.diff(xlims)
    print(xbarr, xbarl, ybar, xlims)
    ax.plot([xbarr, xbarl], [ybar, ybar], '-w', lw=1)

    # panel properties
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
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
