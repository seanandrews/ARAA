import numpy as np
import os
import sys
from post_summary import post_summary
from astropy.io import ascii


CI_levs = [84.135, 15.865]

# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# load naming key
nkey = ascii.read('ALLDISKS.Reffs.dat')


R7, eR7_hi, eR7_lo = np.zeros(len(db)), np.zeros(len(db)), np.zeros(len(db))
fl_R7, lim_R7 = -np.ones(len(db)), np.zeros(len(db))
for i in range(len(nkey)):
    
    # locate source in database
    ind = np.where(db['NAME'] == nkey['name'][i])[0][0]
    print(db['NAME'][ind])

    # load distance posterior
    dpc = np.load('outputs/'+nkey['name'][i]+'.dpc.posterior.npz')['dpc']
    n_dpc = len(dpc)

    # load rhoeff posterior (90% radius)
    rhoeff = np.load('size_posteriors/'+nkey['file'][i]+'.ppost.npz')['Reff90']
    n_rhoeff = len(rhoeff)

    # draw n_rhoeff random samples from dpc posterior
    d = np.random.choice(dpc, n_rhoeff)

    # log10(Reff) posteriors (au)
    Reff = np.log10(rhoeff * d)

    # save the posteriors
    np.savez('outputs/'+db['NAME'][ind]+'.logReff.posterior.npz', logReff=Reff)

    # posterior summary
    pk_Reff, hi_Reff, lo_Reff = post_summary(Reff, prec=0.001)
    lim_Reff = np.percentile(Reff, 95.45)

    # populate the database arrays appropriately
    fl_R7[ind] = 0
    R7[ind] = pk_Reff
    eR7_hi[ind] = hi_Reff
    eR7_lo[ind] = lo_Reff
    lim_R7[ind] = lim_Reff
    

# add these as columns to the database
db['FL_R7'] = fl_R7
db['logR7'] = R7
db['elogR7_hi'] = eR7_hi
db['elogR7_lo'] = eR7_lo
db['LIM_R7'] = lim_R7

# write out the modified database
ascii.write(db, 'temp.csv', format='csv', fast_writer=True, overwrite=True)
