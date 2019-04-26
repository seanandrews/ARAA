import numpy as np
import matplotlib.pyplot as plt

# constants
hh, cc, kk = 6.6269655e-27, 2.99792458e10, 1.380658e-16

# rough frequencies
nu0 = cc / 0.1
nu1 = cc / 0.3

# temperatures
T = np.linspace(15., 60., 100)

# alpha_Pl
alp_pl = np.zeros_like(T)
for i in range(len(T)):
    foo = np.log((np.exp(hh*nu1/(kk*T[i]))-1.)/(np.exp(hh*nu0/(kk*T[i]))-1.))
    alp_pl[i] = 3. + foo / np.log(nu0/nu1)

plt.plot(T, alp_pl)
plt.show()


