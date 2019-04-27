import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time
import urllib
from astropy import constants as c
import dsharp_opac as opacity

# load the DSHARP opacity defaults
d = np.load(opacity.get_datafile('default_opacities_smooth.npz'))
a = d['a']
lam = d['lam']
k_abs = d['k_abs']
k_sca = d['k_sca']
gsca = d['g']

# compute the effective scattering opacity
k_sca_eff = (1 - gsca) * k_sca
lam_avg = [0.13, 0.29]
q = [3.5, 2.5]
res = [opacity.size_average_opacity(lam_avg, a, lam, k_abs, k_sca, q=_q, 
       plot=False) for _q in q]
res_eff = [opacity.size_average_opacity(lam_avg, a, lam, k_abs, k_sca_eff, 
           q=_q, plot=False) for _q in q]
eps = [_res['ka'] / (_res['ka'] + _res['ks']) for _res in res]
eps_eff = [_res_eff['ka'] / (_res_eff['ka'] + _res_eff['ks']) 
           for _res_eff in res_eff]

# form: res[q][type][wl]  eps_eff[q][wl]
kabs_qhi = res[0]['ka'][0]
kabs_qlo = res[1]['ka'][0]
ksca_qhi = res_eff[0]['ks'][0]
ksca_qlo = res_eff[1]['ks'][0]
beta_qhi = res[0]['beta']
beta_qlo = res[1]['beta']
walb_qhi = 1.-eps_eff[0][0]
walb_qlo = 1.-eps_eff[1][0]


# spit out things into numpy save files
np.savez('default_q3.5', acm=a, kabs0=kabs_qhi, ksca0=ksca_qhi, beta=beta_qhi, 
         albedo=walb_qhi)
np.savez('default_q2.5', acm=a, kabs0=kabs_qlo, ksca0=kabs_qlo, beta=beta_qlo,
         albedo=walb_qlo)



# make porous versions

density_water = 0.92
density_silicates = 3.30
density_troilite = 4.83
density_organics = 1.50

constants_default = [opacity.diel_warrenbrandt08(),
                     opacity.diel_draine2003('astrosilicates'),
                     opacity.diel_henning('troilite'),
                     opacity.diel_henning('organics'),]

densities_default = np.array([density_water, density_silicates,
                              density_troilite, density_organics,])

fv_default = np.array([0.3642, 0.1670, 0.0258, 0.4430])

porosity = 0.9
constants_porous = [opacity.diel_vacuum()] + constants_default
densities_porous = np.hstack((0.0, densities_default))
fv_porous = np.hstack((porosity, (1-porosity)*fv_default))

rhos_porous = (fv_porous * densities_porous).sum()

t0 = time.time()
d_por = opacity.diel_mixed(constants_porous, fv_porous, rule='Maxwell-Garnett')
res_por = opacity.get_opacities(a, lam, rho_s=rhos_porous, diel_const=d_por,
                                return_all=True, extrapol=True, 
                                extrapolate_large_grains=True)
t1 = time.time()
print((t1-t0)/60.)

r_por = [opacity.size_average_opacity(lam_avg, a, lam, res_por['k_abs'], res_por['k_sca'], q=_q, plot=False) for _q in q]
re_por = [opacity.size_average_opacity(lam_avg, a, lam, res_por['k_abs'], (1.-res_por['g'])*res_por['k_sca'], q=_q, plot=False) for _q in q]

eps = [_res['ka'] / (_res['ka'] + _res['ks']) for _res in r_por]
eps_eff = [_res_eff['ka'] / (_res_eff['ka'] + _res_eff['ks'])
           for _res_eff in re_por]

# form: res[q][type][wl]  eps_eff[q][wl]
p_kabs_qhi = r_por[0]['ka'][0]
p_kabs_qlo = r_por[1]['ka'][0]
p_ksca_qhi = re_por[0]['ks'][0]
p_ksca_qlo = re_por[1]['ks'][0]
p_beta_qhi = r_por[0]['beta']
p_beta_qlo = r_por[1]['beta']
p_walb_qhi = 1.-eps_eff[0][0]
p_walb_qlo = 1.-eps_eff[1][0]


# spit out things into numpy save files
np.savez('porous_q3.5', acm=a, kabs0=p_kabs_qhi, ksca0=p_ksca_qhi, 
         beta=p_beta_qhi, albedo=p_walb_qhi)
np.savez('porous_q2.5', acm=a, kabs0=p_kabs_qlo, ksca0=p_kabs_qlo, 
         beta=p_beta_qlo, albedo=p_walb_qlo)

