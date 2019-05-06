import numpy as np
import os
import sys
from astropy.io import ascii
from plx_dist import plx_dist
from scipy import stats

plx_shift = 0.03
sys_plx   = 0.1
nsamples  = 100000
CI_levs = [84.135, 15.865, 95.45, 4.55]

# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# calculate distance posterior summaries
d, ed_hi, ed_lo = np.zeros(len(db)), np.zeros(len(db)), np.zeros(len(db))
for i in range(len(db)):
    # modified parallax uncertainty
    eplx  = np.sqrt(db['EPI'][i]**2 + sys_plx**2)

    # distance posterior samples
    pdist = plx_dist(db['PI'][i]+plx_shift, eplx, nsamples=nsamples)

    # find the peak of distance posterior
    CI_dist = np.percentile(pdist, CI_levs)
    kde_dist = stats.gaussian_kde(pdist)
    ndisc = np.int(np.round((CI_dist[0] - CI_dist[1]) / 0.1))
    x_dist = np.linspace(CI_dist[1], CI_dist[0], ndisc)
    pk_dist = x_dist[np.argmax(kde_dist.evaluate(x_dist))]

    # mean, upper, lower 1-sigma
    d[i], ed_hi[i], ed_lo[i] = pk_dist, CI_dist[0]-pk_dist, pk_dist-CI_dist[1]

# add these as columns to the database
db['DPC2'] = d
#db['EDPC_HI'] = ed_hi
#db['EDPC_LO'] = ed_lo

# write out the modified database
ascii.write(db, 'temp.csv', format='csv', fast_writer=True, overwrite=True)
