import numpy as np
import os
import sys
from astropy.io import ascii
from plx_dist import plx_dist
from scipy import stats

def calc_d(db, plx_shift=0.03, sys_plx=0.1, nsamples=100000):

    CI_levs = [84.135, 15.865, 95.45, 4.55]

    # calculate distance posteriors and summaries
    d, edhi, edlo = np.zeros(len(db)), np.zeros(len(db)), np.zeros(len(db))
    for i in range(len(db)):
        # modified parallax uncertainty
        eplx  = np.sqrt(db['EPI'][i]**2 + sys_plx**2)

        # distance posterior samples
        pdist = plx_dist(db['PI'][i]+plx_shift, eplx, nsamples=nsamples)

        # store the distance posterior samples
        np.savez('outputs/'+db['NAME'][i]+'.dpc.posterior.npz', dpc=pdist)

        # find the peak of distance posterior
        CI_dist = np.percentile(pdist, CI_levs)
        kde_dist = stats.gaussian_kde(pdist)
        ndisc = np.int(np.round((CI_dist[0] - CI_dist[1]) / 0.1))
        x_dist = np.linspace(CI_dist[1], CI_dist[0], ndisc)
        pk_d = x_dist[np.argmax(kde_dist.evaluate(x_dist))]

        # mean, upper, lower 1-sigma
        d[i], edhi[i], edlo[i] = pk_d, CI_dist[0]-pk_d, pk_d-CI_dist[1]

    # replace these columns in the modified database
    db['DPC'] = d
    db['EDPC_H'] = edhi
    db['EDPC_L'] = edlo

    # return the modified database
    return db
