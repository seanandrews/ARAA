import numpy as np
import os
import sys
from astropy.io import ascii

### Load the database 

# safe copy + load
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)

base = ((db['FL_R7'] == 0) & (db['FL_R6'] == 0))

R7 = 10.**(db['R7'][base])
R6 = 10.**(db['R6'][base])
name = db['NAME'][base]

ratio = R7/R6
for i in range(len(ratio)): print(name[i], ratio[i])
print(np.median(ratio), np.mean(ratio))
