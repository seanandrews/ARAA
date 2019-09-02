import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import interp1d
import scipy.integrate as sci


disk = 'HTLupA'

x = 0.68

ddir = '/data/sandrews/ALMA_disks/DR/profiles/'
ddir = 'data/'

#rau, ras, Inu, eI, Tb, eTb, Tpl, eTpl = np.loadtxt(ddir+disk+'.profile.txt').T
rau, Inu, eI = np.loadtxt(ddir+disk+'.profile.txt').T
Inu = Inu[rau <= 200]
rau = rau[rau <= 200]


Fcum = sci.cumtrapz(2 * np.pi * Inu * rau, rau, initial=0.)
fint = interp1d(Fcum / Fcum[-1], rau)
reff = fint(x)

print(np.log10(reff))
