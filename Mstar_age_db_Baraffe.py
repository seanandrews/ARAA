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

fL = np.zeros(ndb)

# loop through database
for i in range(ndb):	

    # if there's stellar information, go ahead
    if (np.ma.is_masked(db['logLs'][i])):
        fL[i] = 1
    else:
        if (db['logTeff'][i] <= 3.498): 
            # tracker
            print('start ',i, db['SPT'][i], db['NAME'][i])

            # load the (logT, logL) posteriors
            L = np.load('outputs/'+db['NAME'][i]+'.logL.posterior.npz')['logL']
            T = np.load('outputs/'+db['NAME'][i]+'.logT.posterior.npz')['logT']
        
            # compute and store log(age/yr) and log(Mstar/Msun) posteriors
            pAM = lt_taum(10.**(T[:ntrials]), L[:ntrials], 
                          grid_name='Baraffe15', ntrials=ntrials, burn=nburn, 
                          nwalkers=nwalk)

            plogAGE = pAM[:,0]  
            plogMs  = pAM[:,1]  
            np.savez('outputs/'+db['NAME'][i]+'.age-mass.posterior.npz', \
                     logAGE=plogAGE, logM=plogMs)

            # posterior summaries
            db['logt'][i], db['elogt_H'][i], db['elogt_L'][i] = post_summary(plogAGE, prec=0.01)
            db['logMs'][i], db['elogMs_H'][i], db['elogMs_L'][i] = post_summary(plogMs, prec=0.01)

            # tracker
            print('finish ',i, db['NAME'][i])

# write out the modified database
ascii.write(db, 'temp.csv', format='csv', fast_writer=True, overwrite=True)
