import numpy as np
import os
import sys
from astropy.io import ascii
from scipy import stats
import linmix
import corner
np.random.seed(2)

def do_regression(XX, YY, EX, EY, DELTA, label):

    # run the linear mixture regression
    lm = linmix.LinMix(XX, YY, EX, EY, delta=DELTA)
    lm.run_mcmc(miniter=10000)

    # save the chains
    os.system('rm -rf outfiles/'+label+'.regression.npz')
    np.savez('outfiles/'+label+'.regression.npz', chain=lm.chain)

    # inference summaries
    alp = lm.chain['alpha']
    bet = lm.chain['beta']
    sca = np.sqrt(lm.chain['sigsqr'])
    cor = lm.chain['corr']

    CI_levs = [15.865, 84.135]
    CI_alp = np.percentile(alp, CI_levs)
    kde_alp = stats.gaussian_kde(alp)
    ndisc = np.round((CI_alp[1] - CI_alp[0]) / 0.01)
    x_alp = np.linspace(CI_alp[0], CI_alp[1], ndisc)
    pk_alp = x_alp[np.argmax(kde_alp.evaluate(x_alp))]

    CI_bet = np.percentile(bet, CI_levs)
    kde_bet = stats.gaussian_kde(bet)
    ndisc = np.round((CI_bet[1] - CI_bet[0]) / 0.01)
    x_bet = np.linspace(CI_bet[0], CI_bet[1], ndisc)
    pk_bet = x_bet[np.argmax(kde_bet.evaluate(x_bet))]

    CI_sca = np.percentile(sca, CI_levs)
    kde_sca = stats.gaussian_kde(sca)
    ndisc = np.round((CI_sca[1] - CI_sca[0]) / 0.01)
    x_sca = np.linspace(CI_sca[0], CI_sca[1], ndisc)
    pk_sca = x_sca[np.argmax(kde_sca.evaluate(x_sca))]

    # dump these to a text file
    alp_str = 'alpha = %5.2f + %5.2f / - %5.2f' % \
              (pk_alp, CI_alp[1]-pk_alp, pk_alp-CI_alp[0])
    bet_str = 'beta = %5.2f + %5.2f / - %5.2f' % \
              (pk_bet, CI_bet[1]-pk_bet, pk_bet-CI_bet[0])
    sca_str = 'dispersion = %5.2f + %5.2f / - %5.2f' % \
              (pk_sca, CI_sca[1]-pk_sca, pk_sca-CI_sca[0])

    os.system('rm -rf outfiles/'+label+'.regression.txt')
    f = open('outfiles/'+label+'.regression.txt', 'w')
    f.write(label + '\n')
    f.write('\n')
    f.write(alp_str + '\n')
    f.write(bet_str + '\n')
    f.write(sca_str)
    f.close()

    # make a covariance plot
    os.system('rm -rf outfiles/'+label+'.corner.png')
    posts = np.column_stack([alp, bet, sca, cor])
    levs = 1.-np.exp(-0.5*(np.arange(3)+1)**2)
    fig = corner.corner(posts, plot_datapoints=False, levels=levs,
              labels=[r'$\alpha$', r'$\beta$', r'$\sigma$', r'$\varrho$'])
    fig.savefig('outfiles/'+label+'.corner.png')
    fig.clf()

    return lm.chain
