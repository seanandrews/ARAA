import numpy as np
import sys

def km_estimator(var, flag):
    
    # re-format
    y = var
    delta = np.zeros(len(var), dtype=np.int)
    delta[flag == True] = -1

    # re-arrange to right-censored data (lower limits)
    x = np.max(y)-y
    ntot = len(x)

    # separate out detections only
    xd = x[~flag]
    yd = y[~flag]

    # distinct, ordered values
    isort = np.argsort(xd)
    xds = xd[isort]
    yds = yd[isort]
    u, iuniq = np.unique(xds, return_index=True)
    xds = xds[iuniq]
    yds = yds[iuniq]
    nd = len(xds)

    # K-M variables
    nobs = np.zeros(nd)
    d = np.zeros(nd)
    s = np.zeros(nd)
    for i in range(nd):
        xi = xds[i]
        n1 = np.sum(x >= xi)
        nobs[i] = n1
        n2 = np.sum(xd == xi)
        d[i] = n2
        s[i] = 1-np.float(n2)/n1

    # survivor function and error
    surv = np.zeros(nd)
    var  = np.zeros(nd)
    surv[0] = 1
    var[0] = 0
    for i in range(1,nd):
        s0 = surv[i-1]*s[i-1]
        surv[i] = s0
        var[i] = s0**2*np.sum(d[0:i-1]/(nobs[0:i-1]*(nobs[0:i-1]-1))) 
    kmest = 1.-surv
    kmerr = np.sqrt(var)

    # mean and variance estimates
    s1 = np.zeros(nd)
    d1 = np.zeros(nd)
    for i in range(1,nd):
        s1[i] = surv[i]*(xds[i]-xds[i-1])
        n1 = nobs[i]
        if (n1 > d[i]):
            d1[i] = d[i] / (n1*(n1-d[i]))
    mean_KM = np.max(y) - np.sum(s1)

    km_results = yds, kmest, kmerr, mean_KM

    return km_results
