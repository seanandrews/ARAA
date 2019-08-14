import numpy as np
import os
import sys
from astropy.io import ascii


nsamples  = 500000
CI_levs = [84.135, 15.865]
esys = 0.1


# for safety, copy over database file
os.system('cp -r DISKS.csv temp.csv')

# now load database
db = ascii.read('temp.csv', format='csv', fast_reader=True)

# load naming key
nkey = ascii.read('ALLDISKS.Reffs.dat')


R7, eR7_hi, eR7_lo = np.zeros(len(db)), np.zeros(len(db)), np.zeros_len(db))
fl_R7, lim_R7 = -1*np.ones(len(db)), np.zeros(len(db))
for i in range(len(nkey['file'])):
    
    # locate source in database
    ind = np.where(db['NAME'] == nkey['name'][i])[0][0]

    # load distance posterior
    dpc = np.load('outputs/'+nkey['name'][i]+'.dpc.posterior.npz')['dpc']
    n_dpc = len(dpc)

    # load rhoeff posterior (90% radius)
    rhoeff = np.load('size_posteriors/'+nkey['file'][i]+'.ppost.npz')['Reff90']
    n_rhoeff = len(rhoeff)

    # draw n_rhoeff random samples from dpc posterior
    d = dpc[np.random.random_integers(0, n_dpc, n_rhoeff)]

    # log10(Reff) posteriors (au)
    Reff = np.log10(rhoeff * d)

    # save the posteriors
    np.savez('outputs/'+db['name'][ind]+'.logReff.posterior.npz', logReff=Reff)

    # posterior summary
    CI_Reff = np.percentile(Reff, CI_levs)
    kde_Reff = stats.gaussian_kde(Reff)
    ndisc = np.round((CI_Reff[0] - CI_Reff[1]) / 0.001)
    x_Reff = np.linspace(CI_Reff[1], CI_Reff[0], ndisc)
    pk_Reff = x_Reff[np.argmax(kde_Reff.evaluate(x_Reff))]
    hi_Reff = CI_Reff[0] - pk_Reff
    lo_Reff = pk_Reff - CI_Reff[1]
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
