import numpy as np
import matplotlib.pyplot as plt

rau = np.logspace(-1, 3, 1000)

# Isella et al. 2009
Rt = [55., 21., 86., 28., 21., 110., 60., 25., 43., 66., 20.]
St = [10., 608., 1.5, 13., 80., 4., 31., 58., 12., 4.7, 50.]
gam = [-0.3, -0.5, 0.8, 0.0, -0.3, 0.7, -0.8, -0.1, 0.8, 0.5, 0.1]
bmaj = np.array([1.05, 0.43, 0.82, 0.80, 0.46, 0.87, 0.83, 0.89, 0.82, 1.42, 1.45])
bmin = np.array([0.72, 0.27, 0.60, 0.58, 0.34, 0.65, 0.70, 0.60, 0.69, 0.85, 0.91])
fwhm = 140.*np.sqrt(bmaj*bmin)


sig = np.zeros((len(Rt), len(rau)))
plt.axis([20, 500, 1e-4, 10])
for i in range(len(Rt)):
    sig[i,:] = St[i] * (Rt[i]/rau)**gam[i] * \
               np.exp( -1. * ((rau/Rt[i])**(2.-gam[i])-1.) / (2.*(2.-gam[i])))
    incond = (rau >= fwhm[i])
    plt.loglog(rau[incond], 0.01*sig[i,incond])
    plt.loglog(rau, 1e-1*(rau/100.)**(-3.), '--')
plt.show()
# very crudely, see things like 1/r**(3/2) to 1/r**10






