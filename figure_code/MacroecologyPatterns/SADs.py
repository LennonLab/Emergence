from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
from random import randint
import numpy as np
import os
from scipy import stats
import sys
from random import shuffle
from math import isnan
from scipy.stats.kde import gaussian_kde
from numpy import empty

mydir = os.path.expanduser('~/GitHub/Emergence')
tools = os.path.expanduser(mydir + "/tools")
data = mydir + '/results/simulated_data/SAR-Data.csv'

sys.path.append(tools + "/DiversityTools/macroecotools")
import macroecotools as mct
sys.path.append(tools + "/DiversityTools/macroeco_distributions")
from macroeco_distributions import pln, pln_solver
sys.path.append(tools + "/DiversityTools/mete")
import mete



def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(0, 1, n)
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



p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 5
fs = 14 # fontsize

#### plot figure ###############################################################
fig = plt.figure(figsize=(3, 2))
fig.add_subplot(1, 1, 1)

data = mydir + '/results/simulated_data/RAD-Data.csv'

RADs = []
with open(data) as f:
    for d in f:
        d = list(eval(d))

        sim = d.pop(0)
        ct = d.pop(0)

        d = sorted(d, reverse=True)
        RADs.append(d)

print 'Number of RADs:', len(RADs)

mete_r2s = []
pln_r2s = []

shuffle(RADs)
for i, obs in enumerate(RADs):

    N = int(sum(obs))
    S = int(len(obs))

    if S > 4 and N > 10 and obs.count(1)/len(obs) < 0.5:

        result = mete.get_mete_rad(S, N)
        pred1 = np.log10(result[0])
        obs1 = np.log10(obs)
        mete_r2 = mct.obs_pred_rsquare(np.array(obs1), np.array(pred1))
        mete_r2s.append(mete_r2)

        pred = get_rad_from_obs(obs, 'pln')
        pred1 = np.log10(pred)
        pln_r2 = mct.obs_pred_rsquare(np.array(obs1), np.array(pred1))
        pln_r2s.append(pln_r2)

        print i, 'N:', N, ' S:', S, ' n:', len(pln_r2s), ' |  mete:', mete_r2, '  pln:', pln_r2

    if len(pln_r2s) > 20: break


kernel = 0.5
D = get_kdens_choose_kernel(mete_r2s, kernel)
plt.plot(D[0],D[1],color = '0.5', lw=3, alpha = 0.99,label= 'log-series')

D = get_kdens_choose_kernel(pln_r2s, kernel)
plt.plot(D[0],D[1],color = '0.1', lw=3, alpha = 0.99, label= 'lognormal')

plt.xlim(0.5, 1)
plt.legend(loc='best', fontsize=fs-5, frameon=False)
plt.xlabel('$r$'+r'$^{2}$', fontsize=fs+3)
plt.ylabel('$density$', fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs-3)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/SADs.png', dpi=200, bbox_inches = "tight")
plt.close()
