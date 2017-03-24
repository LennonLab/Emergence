from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import statsmodels.api as sm

mydir = os.path.expanduser('~/GitHub/simplex')

def figplot(x, y, xlab, ylab, fig, n):
    fig.add_subplot(3, 3, n)
    plt.scatter(x, y, lw=0.5, color='0.7', s = 4)
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    plt.plot(x, y, lw=_lw, color='k')
    plt.tick_params(axis='both', labelsize=6)
    plt.xlabel(xlab, fontsize=9)
    plt.ylabel(ylab, fontsize=9)
    plt.xlim(1,5.5)
    return fig

p, fr, _lw, w, sz = 1, 0.1, 1.5, 1, 5

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df = df[df['ct'] > 300]

df2 = pd.DataFrame({'length' : df['length'].groupby(df['sim']).mean()})
df2['sim'] = df['sim'].groupby(df['sim']).mean()
df2['flow'] = df['flow.rate'].groupby(df['sim']).mean()
df2['tau'] = np.log10(df2['length']**p/df2['flow'])

df2['Grow'] = df['active.avg.per.capita.growth'].groupby(df['sim']).max()
df2['Maint'] = np.log10(df['active.avg.per.capita.maint'].groupby(df['sim']).mean())
df2['Disp'] = df['active.avg.per.capita.active.dispersal'].groupby(df['sim']).max()
df2['RPF'] = df['dormant.avg.per.capita.rpf'].groupby(df['sim']).mean()
df2['Eff'] = df['active.avg.per.capita.efficiency'].groupby(df['sim']).min()
df2['MF'] = np.log10(df['all.avg.per.capita.mf'].groupby(df['sim']).mean())

xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fig = plt.figure()

ylab = 'Growth rate'
fig = figplot(df2['tau'], df2['Grow'], xlab, ylab, fig, 1)

ylab = 'Maintenance energy'
fig = figplot(df2['tau'], df2['Maint'], xlab, ylab, fig, 2)

ylab = 'Active disperal rate'
fig = figplot(df2['tau'], df2['Disp'], xlab, ylab, fig, 4)

ylab = 'Random resuscitation\nfrom dormancy'
fig = figplot(df2['tau'], df2['RPF'], xlab, ylab, fig, 5)

ylab = 'Resource specialization'
fig = figplot(df2['tau'], df2['Eff'], xlab, ylab, fig, 7)

ylab = 'Decrease of maintenance\nenergy when dormant'
fig = figplot(df2['tau'], df2['MF'], xlab, ylab, fig, 8)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Trait_vs_Tau.png', dpi=200, bbox_inches = "tight")
