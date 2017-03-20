from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy import stats

mydir = os.path.expanduser('~/GitHub/simplex')
tools = os.path.expanduser(mydir + "/tools")

p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 1

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width'].groupby(df['sim']).mean()})
df2['N'] = np.log10(df['total.abundance'].groupby(df['sim']).mean())
df2['NS'] = np.log10(df['avg.pop.size'].groupby(df['sim']).mean())
df2['var'] = np.log10(df['pop.var'].groupby(df['sim']).mean())
df2['Pdens'] = np.log10(df['avg.pop.size'].groupby(df['sim']).mean()/df2['width']**3)

df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

#### plot figure ###############################################################
fs = 12 # fontsize
fig = plt.figure()
fig.add_subplot(1, 1, 1)

df2 = df2[df2['N'] > 0]
Nlist = df2['NS'].tolist()
Vlist = df2['var'].tolist()

plt.scatter(Nlist, Vlist, lw=_lw, color='0.2', s = sz)
m, b, r, p, std_err = stats.linregress(Nlist, Vlist)
Nlist = np.array(Nlist)
plt.plot(Nlist, m*Nlist + b, '-', color='k', label='Taylor\'s Law, '+'$z$ = '+str(round(m,2)))
xlab = r"$log_{10}$"+'(Pop mean)'
ylab = r"$log_{10}$"+'(variance)'
plt.xlabel(xlab, fontsize=fs)
plt.tick_params(axis='both', labelsize=fs-3)
plt.ylabel(ylab, fontsize=fs)
plt.legend(loc=2, fontsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/TaylorsLaw.png', dpi=200, bbox_inches = "tight")
plt.close()
