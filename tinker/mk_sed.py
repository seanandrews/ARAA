import numpy as np
import matplotlib.pyplot as plt
import os
import sys


# constants
au = 1.496e13
c  = 2.9979e10
h  = 6.6260755e-27
k  = 1.380658e-16
pc = 3.0857e18
Lsun = 3.826e33
nr = 500
nw = 100

# parameters
r0 = 10.
T0 = 100.
q  = 0.5
sig0 = 100.
p  = 0.0
rout = 1000.
incl = 0.

# radius grid (au)
r = np.logspace(-1, 3, nr)

# wavelength grid (um) + frequency grid (Hz)
wl = np.logspace(-2, 4, nw)
nu = c * 1e4 / wl

# temperature profile
T = T0 * (r/r0)**(-q)

# opacity and surface density profiles
kap = 10. * (nu / 1e12)**1.
sig = sig0 * (r/r0)**(-p)
sig[r > rout] = 0.

# compute emission
Bnu = np.zeros((nr, nw))
tau = np.zeros((nr, nw))
raw_sed = np.zeros(nw)
for j in range(nw):
    # Planck function
    Bnu[:,j] = (2.*h*nu[j]**3/c**2) / (np.exp(h*nu[j]/(k*T))-1.)
    # optical depths
    tau[:,j] = kap[j] * sig / np.cos(np.radians(incl))
    # calculate emission
    raw_sed[j] = np.trapz(Bnu[:,j] * (1. - np.exp(-tau[:,j])) * r * au, r * au)

# calculate the spectrum
fnu = 1e23 * 2. * np.pi * np.cos(np.radians(incl)) * raw_sed / (140.*pc)**2

# convert to an SED
sed = 4.*np.pi*(140.*pc)**2 * nu * (1e-23*fnu) / Lsun

plt.axis([1, 1e4, 1e-3, 10])
plt.loglog(wl, sed)
plt.show()
