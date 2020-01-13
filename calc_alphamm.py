import numpy as np
import os
import sys
import time
from astropy.io import ascii
from plx_dist import plx_dist
from scipy import stats

nsamples  = 500000
CI_levs = [84.135, 15.865]
esys = 0.1

# for safety, copy over database file
os.system('cp -r DISKS.csv temp2.csv')

# now load database
db = ascii.read('temp2.csv', format='csv', fast_reader=True)

a67, ea67_hi, ea67_lo = np.zeros(len(db)), np.zeros(len(db)), np.zeros(len(db))
fl_a67, lim_a67 = np.zeros(len(db)), np.zeros(len(db))
a36, ea36_hi, ea36_lo = np.zeros(len(db)), np.zeros(len(db)), np.zeros(len(db))
fl_a36, lim_a36 = np.zeros(len(db)), np.zeros(len(db))

t0 = time.time()
indices = [249]
for i in indices:	#range(len(db)):

    # no index information available (no measurement in 1 band, or upper limits
    # in both bands)
    if (db['FL_B7'][i] == -1 or db['FL_B6'][i] == -1):
        fl_a67[i], a67[i], ea67_hi[i], ea67_lo[i], lim_a67[i] = -1, 0, 0, 0, 0
    if (db['FL_B6'][i] == -1 or db['FL_B3'][i] == -1):
        fl_a36[i], a36[i], ea36_hi[i], ea36_lo[i], lim_a36[i] = -1, 0, 0, 0, 0
    if (db['FL_B7'][i] == 1 and db['FL_B6'][i] == 1):
        fl_a67[i], a67[i], ea67_hi[i], ea67_lo[i], lim_a67[i] = -1, 0, 0, 0, 0
    if (db['FL_B6'][i] == 1 and db['FL_B3'][i] == 1):
        fl_a36[i], a36[i], ea36_hi[i], ea36_lo[i], lim_a36[i] = -1, 0, 0, 0, 0


    # dual-band detections 
    if (db['FL_B7'][i] == 0 and db['FL_B6'][i] == 0):

        # simple calculation of posterior mean (for db)
        fl_a67[i], lim_a67[i] = 0, 0
        a67[i] = np.log(db['F_B7'][i]/db['F_B6'][i]) / \
                 np.log(db['nu_B7'][i]/db['nu_B6'][i])

        # full posterior sampling
        sigB7 = np.sqrt(db['eF_B7'][i]**2 + (esys*db['F_B7'][i])**2)
        sigB6 = np.sqrt(db['eF_B6'][i]**2 + (esys*db['F_B6'][i])**2)
        fB7 = np.abs(np.random.normal(db['F_B7'][i], sigB7, nsamples))
        fB6 = np.abs(np.random.normal(db['F_B6'][i], sigB6, nsamples))
        alp67 = np.log(fB7 / fB6) / np.log(db['nu_B7'][i] / db['nu_B6'][i])
        np.savez('outputs/'+db['NAME'][i]+'.alpha67.posterior.npz', amm=alp67)

        # calculate / assign confidence intervals
        CI_dist = np.percentile(alp67, CI_levs)
        ea67_hi[i], ea67_lo[i] = CI_dist[0]-a67[i], a67[i]-CI_dist[1]
       

    if (db['FL_B6'][i] == 0 and db['FL_B3'][i] == 0):

        # simple calculation of posterior mean (for db)
        fl_a36[i], lim_a36[i] = 0, 0
        a36[i] = np.log(db['F_B6'][i]/db['F_B3'][i]) / \
                 np.log(db['nu_B6'][i]/db['nu_B3'][i])

        # full posterior sampling
        sigB6 = np.sqrt(db['eF_B6'][i]**2 + (esys*db['F_B6'][i])**2)
        sigB3 = np.sqrt(db['eF_B3'][i]**2 + (esys*db['F_B3'][i])**2)
        fB6 = np.abs(np.random.normal(db['F_B6'][i], sigB6, nsamples))
        fB3 = np.abs(np.random.normal(db['F_B3'][i], sigB3, nsamples))
        alp36 = np.log(fB6 / fB3) / np.log(db['nu_B6'][i] / db['nu_B3'][i])
        np.savez('outputs/'+db['NAME'][i]+'.alpha36.posterior.npz', amm=alp36)

        # calculate / assign confidence intervals
        CI_dist = np.percentile(alp36, CI_levs)
        ea36_hi[i], ea36_lo[i] = CI_dist[0]-a36[i], a36[i]-CI_dist[1]


    # single-band detections for *lower* limits on alpha
    if (db['FL_B7'][i] == 0 and db['FL_B6'][i] == 1):
        fl_a67[i], a67[i] = 1, -1
        lim_a67[i] = np.log(db['F_B7'][i]/db['LIM_B6'][i]) / \
                     np.log(db['nu_B7'][i]/db['nu_B6'][i])
    if (db['FL_B6'][i] == 0 and db['FL_B3'][i] == 1):
        fl_a36[i], a36[i] = 1, -1
        lim_a36[i] = np.log(db['F_B6'][i]/db['LIM_B3'][i]) / \
                     np.log(db['nu_B6'][i]/db['nu_B3'][i])

    # single-band detections for *upper* limits on alpha
    if (db['FL_B7'][i] == 1 and db['FL_B6'][i] == 0):
        fl_a67[i], a67[i] = 2, -1
        lim_a67[i] = np.log(db['LIM_B7'][i]/db['F_B6'][i]) / \
                     np.log(db['nu_B7'][i]/db['nu_B6'][i])
    if (db['FL_B6'][i] == 1 and db['FL_B3'][i] == 0):
        fl_a36[i], a36[i] = 2, -1
        lim_a36[i] = np.log(db['LIM_B6'][i]/db['F_B3'][i]) / \
                     np.log(db['nu_B6'][i]/db['nu_B3'][i])


print(time.time()-t0)


# add the indices as columns to the database
db['FL_A67'] = fl_a67
db['A67'] = a67
db['eA67_hi'] = ea67_hi
db['eA67_lo'] = ea67_lo
db['LIM_A67'] = lim_a67
db['FL_A36'] = fl_a36
db['A36'] = a36
db['eA36_hi'] = ea36_hi
db['eA36_lo'] = ea36_lo
db['LIM_A36'] = lim_a36

# write out the modified database
ascii.write(db, 'temp2.csv', format='csv', fast_writer=True, overwrite=True)
