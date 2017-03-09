from __future__ import division
import  matplotlib.pyplot as plt
import sys
import os
from random import shuffle
import numpy as np

########### PATHS ##############################################################

mydir = os.path.expanduser("~/GitHub/residence-time")
tools = os.path.expanduser(mydir + "/tools")

sys.path.append(tools + "/DiversityTools/macroeco_distributions")
import macroeco_distributions as md
sys.path.append(tools + "/DiversityTools/distributions")
import distributions as dist
sys.path.append(tools + "/DiversityTools/macroecotools")
import macroecotools as mct
sys.path.append(tools + "/metrics")
import metrics
sys.path.append(tools + "/DiversityTools/mete")
import mete
#sys.path.append(tools + "/pln")
#import pln

from scipy.stats.kde import gaussian_kde
from macroeco_distributions import pln, pln_solver
from numpy import empty


def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D



def get_rad_pln(S, mu, sigma, lower_trunc = True):
    """Obtain the predicted RAD from a Poisson lognormal distribution"""
    abundance = list(empty([S]))
    rank = range(1, int(S) + 1)
    cdf_obs = [(rank[i]-0.5) / S for i in range(0, int(S))]
    j = 0
    cdf_cum = 0
    i = 1
    while j < S:
        cdf_cum += pln.pmf(i, mu, sigma, lower_trunc)
        while cdf_cum >= cdf_obs[j]:
            abundance[j] = i
            j += 1
            if j == S:
                abundance.reverse()
                return abundance
        i += 1



def get_rad_from_obs(ab, dist):
    mu, sigma = pln_solver(ab)
    pred_rad = get_rad_pln(len(ab), mu, sigma)
    return pred_rad



data = mydir + '/results/simulated_data/protected/RAD-Data.csv'

RADs = []

with open(data) as f:

    for d in f:
        d = list(eval(d))
        sim = d.pop(0)
        ct = d.pop(0)

        if len(d) >= 10:
            d = sorted(d, reverse=True)
            RADs.append(d)

print 'Number of RADs:', len(RADs)

mete_r2s = []
zipf_r2s = []
pln_r2s = []

shuffle(RADs)
for i, obs in enumerate(RADs):

    N = int(sum(obs))
    S = int(len(obs))

    print i, N, S, len(pln_r2s)

    if S >= 10 and N > 50:

        if N < 10000:

            result = mete.get_mete_rad(S, N)
            predRAD = result[0]
            mete_r2 = mct.obs_pred_rsquare(np.array(obs), np.array(predRAD))
            mete_r2s.append(mete_r2)

            #zipf_pred = dist.zipf(obs)
            #predRAD = zipf_pred.from_cdf()
            #zipf_r2 = mct.obs_pred_rsquare(np.array(obs), np.array(predRAD))
            #zipf_r2s.append(zipf_r2)

            predRAD = get_rad_from_obs(obs, 'pln')
            pln_r2 = mct.obs_pred_rsquare(np.array(obs), np.array(predRAD))
            pln_r2s.append(pln_r2)

    if len(pln_r2s) > 200: break


fig = plt.figure(111)
kernel = 0.5

D = get_kdens_choose_kernel(mete_r2s, kernel)
plt.plot(D[0],D[1],color = '0.3', lw=3, alpha = 0.99,label= 'METE')

#D = get_kdens_choose_kernel(zipf_r2s, kernel)
#plt.plot(D[0],D[1],color = 'c', lw=3, alpha = 0.99,label= 'Zipf')

D = get_kdens_choose_kernel(pln_r2s, kernel)
plt.plot(D[0],D[1],color = 'm', lw=3, alpha = 0.99, label= 'PLN')

plt.xlim(0.0, 1)
plt.legend(loc=2, fontsize=16)
plt.xlabel('$r$'+r'$^{2}$', fontsize=22)
plt.ylabel('$density$', fontsize=22)
plt.savefig(mydir + '/results/figures/SADfits.png', dpi=600, bbox_inches = "tight")
plt.close()
