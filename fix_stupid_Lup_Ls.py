import numpy as np
from scipy import stats
import os
import sys

def fix_stupid_Lup_Ls(Lstar, eLstar):

    CI_levs = [84.135, 15.865]
    nsamples = 1000000

    plogLstar = np.random.normal(np.log10(Lstar), eLstar / \
                                 Lstar / np.log(10), nsamples) 

    CI_logLstar  = np.percentile(plogLstar, CI_levs)
    pk_logLstar  = np.log10(Lstar)

    return (pk_logLstar, 0.5*(CI_logLstar[0]-CI_logLstar[1]))
