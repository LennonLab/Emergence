from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats.kde import gaussian_kde
import sys


mydir = os.path.expanduser('~/GitHub/Emergence')
tools = os.path.expanduser(mydir + "/tools")
data = mydir + '/results/simulated_data/SAR-Data.csv'

def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(0, 1, n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


z_nest = []
z_rand = []

with open(data) as f:
    for d in f:
        d = list(eval(d))
        sim = d.pop(0)
        ct = d.pop(0)
        if ct > 100:
            z1, z2 = d
            z_nest.append(z1)
            z_rand.append(z2)

fs = 14
fig = plt.figure(figsize=(3, 2))
fig.add_subplot(1, 1, 1)
kernel = 0.1

D = get_kdens_choose_kernel(z_nest, kernel)
plt.plot(D[0],D[1],color = 'k', lw=3, alpha = 0.99, label= 'Nested SAR '+'$z$'+'-values')
D = get_kdens_choose_kernel(z_rand, kernel)
plt.plot(D[0],D[1],color = '0.5', lw=3, alpha = 0.99, label= 'R.A. SAR '+'$z$'+'-values')

plt.legend(loc='best', fontsize=fs-5, frameon=False)
plt.xlabel('$z$', fontsize=fs+6)
plt.ylabel('$density$', fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs-3)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/SAR.png', dpi=200, bbox_inches = "tight")
plt.close()
