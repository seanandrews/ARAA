import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from astropy.io import ascii

# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# select targets with both a robust index and B6 flux density
ok = ((db['FL_A67'] == 0) & (db['FL_B6'] == 0) & (db['FL_MULT'] != 'J'))

# sample
name = db['NAME'][ok]
dpc  = db['DPC'][ok]
F_B6 = db['F_B6'][ok]
alp  = db['A67'][ok]
print(len(name))

# shift to B6 luminosities (in mJy at some reference distance)
d0 = 140.
L_B6 = F_B6 * (dpc/d0)**2

# plot
plt.semilogx(L_B6, alp, 'ok')
plt.axis([0.2, 5000, 0., 7.])
plt.show()


# check indices against one another
ok = ((db['FL_A67'] == 0) & (db['FL_A36'] == 0))

# sample
name = db['NAME'][ok]
a67  = db['A67'][ok]
a36  = db['A36'][ok]
print(len(name))

# plot
plt.plot(a67, a36, 'ok')
plt.axis([0., 7., 0., 7.])
plt.show()



