from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy import stats

mydir = os.path.expanduser('~/GitHub/Emergence')
tools = os.path.expanduser(mydir + "/tools")

_lw = 2
sz = 20

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df2 = pd.DataFrame({'length' : df['length'].groupby(df['sim']).mean()})
df2['NS'] = np.log10(df['avg.pop.size'].groupby(df['sim']).mean())
df2['var'] = np.log10(df['pop.var'].groupby(df['sim']).mean())
df2 = df2[df2['var'] > 1]

#### plot figure ###############################################################
fs = 14
fig = plt.figure(figsize=(3, 2))
fig.add_subplot(1, 1, 1)

Nlist = df2['NS'].tolist()
Vlist = df2['var'].tolist()

plt.scatter(Nlist, Vlist, lw=_lw, color='0.7', s = sz)
m, b, r, p, std_err = stats.linregress(Nlist, Vlist)
Nlist = np.array(Nlist)
plt.plot(Nlist, m*Nlist + b, '-', color='k', label='$z$ = '+str(round(m,2)), lw=_lw)
xlab = r"$log_{10}$"+'(mean)'
ylab = r"$log_{10}$"+'(variance)'
plt.xlabel(xlab, fontsize=fs)
plt.tick_params(axis='both', labelsize=fs-3)
plt.ylabel(ylab, fontsize=fs)
plt.legend(loc='best', fontsize=fs-3, frameon=False)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/TaylorsLaw.png', dpi=200, bbox_inches = "tight")
plt.close()
