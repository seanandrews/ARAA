import numpy as np
import os
import sys
from post_summary import post_summary


# joint names
oname = 'UZ_Tau_EW'
names = ['UZ_Tau_Ea', 'UZ_Tau_Eb', 'UZ_Tau_Wa', 'UZ_Tau_Wb']

# get the # of posterior samples
npost = len(np.load('outputs/'+names[0]+'.age-mass.posterior.npz')['logM'])
npostL = len(np.load('outputs/'+names[0]+'.logL.posterior.npz')['logL'])

# load the age and mass posteriors for each component
ages, masses = np.zeros((npost, len(names))), np.zeros((npost, len(names)))
lums = np.zeros((npostL, len(names)))
for i in range(len(names)):
    fname = 'outputs/'+names[i]+'.age-mass.posterior.npz'
    ages[:,i], masses[:,i] = np.load(fname)['logAGE'], np.load(fname)['logM']
    lums[:,i] = np.load('outputs/'+names[i]+'.logL.posterior.npz')['logL']

# compute the total mass and average age posteriors
mage, mtot = np.average(ages, 1), np.log10(np.sum(10.**masses,1))
ltot = np.log10(np.sum(10.**lums,1))

# posterior summaries
A, eAhi, eAlo = post_summary(mage, prec=0.01)
M, eMhi, eMlo = post_summary(mtot, prec=0.01)
L, eLhi, eLlo = post_summary(ltot, prec=0.01)

print(L, eLhi, eLlo)
print(M, eMhi, eMlo)
print(A, eAhi, eAlo)

np.savez('outputs/'+oname+'.logL.posterior.npz', logL=ltot)
np.savez('outputs/'+oname+'.age-mass.posterior.npz', logM=mtot, logAGE=mage)
