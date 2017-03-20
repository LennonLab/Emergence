from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import os
from scipy import stats
from random import shuffle
from math import isnan
from scipy.stats.kde import gaussian_kde


mydir = os.path.expanduser('~/GitHub/simplex')
tools = os.path.expanduser(mydir + "/tools")
data = mydir + '/results/simulated_data/SAR-Data.csv'

def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list), n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


SARs = []
with open(data) as f:
    for d in f:
        d = list(eval(d))
        sim = d.pop(0)
        ct = d.pop(0)
        if max(d) < 10: continue
        SARs.append(d)
        if len(SARs) > 10000: break


z_vals = []
shuffle(SARs)
for sar in SARs:
    area_vector = np.array(range(1, len(sar)+1))**2
    m, b, r, p, std_err = stats.linregress(np.log10(area_vector), np.log10(sar))
    if isnan(m): continue
    z_vals.append(m)
    if len(z_vals) > 10000:
        break

fs = 12
fig = plt.figure()
fig.add_subplot(1, 1, 1)
kernel = 0.2

D = get_kdens_choose_kernel(z_vals, kernel)
plt.plot(D[0],D[1],color = '0.5', lw=3, alpha = 0.99, label= 'SAR '+'$z$'+'-values')
plt.legend(loc=1, fontsize=fs+1)
plt.xlabel('$z$', fontsize=fs+5)
plt.ylabel('density', fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/SAR.png', dpi=200, bbox_inches = "tight")
plt.close()
