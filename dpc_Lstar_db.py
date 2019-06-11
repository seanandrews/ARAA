import numpy as np
import os
import sys
from astropy.io import ascii
from plx_dist import plx_dist
from post_summary import post_summary


# load database
os.system('cp -r DISKS.csv temp.csv')
db = ascii.read('temp.csv', format='csv', fast_reader=True)
ndb = len(db)

# parallax systematics
plx_shift, sys_plx = 0.03, 0.1

# posterior samples
ns = 100000

d, edhi, edlo = np.zeros(ndb), np.zeros(ndb), np.zeros(ndb)
L, eLhi, eLlo = np.zeros(ndb), np.zeros(ndb), np.zeros(ndb)
fL = np.zeros(ndb)

# loop through database
for i in range(ndb):

    # tracker
    print(i, db['NAME'][i])

    # modified parallax uncertainty
    eplx = np.sqrt(db['EPI'][i]**2 + sys_plx**2)

    # calculate and store the distance posterior samples
    p_dpc = plx_dist(db['PI'][i]+plx_shift, eplx, nsamples=ns)
    np.savez('outputs/'+db['NAME'][i]+'.dpc.posterior.npz', dpc=p_dpc)

    # distance posterior summaries
    d[i], edhi[i], edlo[i] = post_summary(p_dpc, prec=0.1)   

    # calculate and store the luminosity posterior samples
    if (np.ma.is_masked(db['logLs_in'][i])):
        fL[i] = 1
    else:
        p_logL = np.random.normal(db['logLs_in'][i], db['elogLs_in'][i], ns) + \
                  2.*np.log10(p_dpc / db['d_in'][i])
        np.savez('outputs/'+db['NAME'][i]+'.logL.posterior.npz', logL=p_logL)

        # luminosity posterior summaries
        L[i], eLhi[i], eLlo[i] = post_summary(p_logL, prec=0.001)

        # calculate and store effective temperature posterior samples for (M,t)
        p_logT = np.random.normal(db['logTeff'][i], db['elogTeff'][i], ns)
        np.savez('outputs/'+db['NAME'][i]+'.logT.posterior.npz', logT=p_logT)


# replace the new values in the relevant database columns
db['DPC'] = d
db['EDPC_H'] = edhi
db['EDPC_L'] = edlo
db['logLs'] = np.ma.masked_array(L, fL)
db['elogLs_H'] = np.ma.masked_array(eLhi, fL)
db['elogLs_L'] = np.ma.masked_array(eLlo, fL)

# write out the modified database
ascii.write(db, 'temp.csv', format='csv', fast_writer=True, overwrite=True)
