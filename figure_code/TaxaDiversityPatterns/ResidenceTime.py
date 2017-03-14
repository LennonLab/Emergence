from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import statsmodels.api as sm

mydir = os.path.expanduser('~/GitHub/simplex')

def figplot(x, y, xlab, ylab, fig, n):
    fig.add_subplot(3, 3, n)
    #plt.hexbin(df2['tau'], df2['N'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
    #plt.scatter(x, y, lw=_lw, color='0.2', s = sz)
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    plt.plot(x, y, lw=_lw, color='m')
    plt.xlabel(xlab, fontsize=fs+3)
    plt.tick_params(axis='both', labelsize=fs)
    plt.ylabel(ylab, fontsize=fs+3)
    return fig


p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 5

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df = df[df['ct'] > 200]
df = df[df['ct'] < 300]

df2 = pd.DataFrame({'width' : df['width'].groupby(df['sim']).mean()})

df2['sim'] = df['sim'].groupby(df['sim']).mean()
df2['flow'] = df['flow.rate'].groupby(df['sim']).mean()
df2['tau'] = np.log10(df2['width']**p/df2['flow'])
df2['R'] = df['resource.particles'].groupby(df['sim']).mean()
df2['N'] = df['total.abundance'].groupby(df['sim']).mean()
df2['Prod'] = df['ind.production'].groupby(df['sim']).mean()
df2['S'] = df['species.richness'].groupby(df['sim']).mean()
df2['E'] = df['simpson.e'].groupby(df['sim']).mean()
df2['W'] = df['Whittakers.turnover'].groupby(df['sim']).mean()
df2['Dorm'] = df['Percent.Dormant'].groupby(df['sim']).mean()

#df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

#### plot figure ###############################################################
xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 6 # fontsize
fig = plt.figure()

x = df2['tau']
y = df2['N']
ylab = 'N'
fig = figplot(x, y, xlab, ylab, fig, 1)

x = df2['tau'].tolist()
y = df2['Prod'].tolist()
ylab = 'Productivity'
fig = figplot(x, y, xlab, ylab, fig, 2)

x = df2['tau']
y = df2['S']
ylab = 'S'
fig = figplot(x, y, xlab, ylab, fig, 4)

x = df2['tau']
y = df2['E']
ylab = 'Evenness'
fig = figplot(x, y, xlab, ylab, fig, 5)

x = df2['tau']
y = np.log10(df2['W'])
ylab = r"$log_{10}$"+'(' + r"$\beta$" +')'
fig = figplot(x, y, xlab, ylab, fig, 7)

x = df2['tau']
y = df2['Dorm']
ylab = '%Dormant'
fig = figplot(x, y, xlab, ylab, fig, 8)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Taxa_vs_Tau.png', dpi=200, bbox_inches = "tight")
plt.close()
