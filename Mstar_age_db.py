import numpy as np
import os
import sys
import time
from astropy.io import ascii
from lt_taum import lt_taum
from post_summary import post_summary


# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)
ndb = len(db)

# posterior samples
nburn = 5000
nwalk = 10
ntrials = 30000 + nburn

M, eMhi, eMlo = np.zeros(ndb), np.zeros(ndb), np.zeros(ndb)
A, eAhi, eAlo = np.zeros(ndb), np.zeros(ndb), np.zeros(ndb)
fL = np.zeros(ndb)

# loop through database
for i in np.arange(483, 589):    #range(ndb):

    # tracker
    print('start ',i, db['NAME'][i])

    # if there's stellar information, go ahead
    if (np.ma.is_masked(db['logLs'][i])):
        fL[i] = 1
    else:
        # load the (logT, logL) posteriors
        logL = np.load('outputs/'+db['NAME'][i]+'.logL.posterior.npz')['logL']
        logT = np.load('outputs/'+db['NAME'][i]+'.logT.posterior.npz')['logT']
        
        # compute and store log(age/yr) and log(Mstar/Msun) posteriors
        if (db['logTeff'][i] > 3.498):
            pAM = lt_taum(10.**(logT[:ntrials]), logL[:ntrials], 
                          grid_name='MIST', ntrials=ntrials, burn=nburn, 
                          nwalkers=nwalk)
        else:
            pAM = lt_taum(10.**(logT[:ntrials]), logL[:ntrials], 
                          grid_name='Baraffe15', ntrials=ntrials, burn=nburn, 
                          nwalkers=nwalk)
        plogAGE = pAM[:,0]  
        plogMs  = pAM[:,1]  
        np.savez('outputs/'+db['NAME'][i]+'.age-mass.posterior.npz', \
                 logAGE=plogAGE, logM=plogMs)

        # posterior summaries
        A[i], eAhi[i], eAlo[i] = post_summary(plogAGE, prec=0.01)
        M[i], eMhi[i], eMlo[i] = post_summary(plogMs, prec=0.01)

        # replace the new values in the relevant database columns
        db['logMs'][i] = M[i]
        db['elogMs_H'][i] = eMhi[i]
        db['elogMs_L'][i] = eMlo[i]
        db['logt'][i] = A[i]
        db['elogt_H'][i] = eAhi[i]
        db['elogt_L'][i] = eAlo[i]

    # tracker
    print('finish ',i, db['NAME'][i])


# write out the modified database
ascii.write(db, 'temp.csv', format='csv', fast_writer=True, overwrite=True)
