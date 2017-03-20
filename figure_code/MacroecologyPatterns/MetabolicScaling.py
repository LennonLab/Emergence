from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy import stats

def figplot(x, y, xlab, ylab, fig, n):
    x = np.log10(x)
    y = np.log10(y)
    fig.add_subplot(3, 3, n)
    plt.scatter(x, y, lw=0.5, color='0.2', s = 4)
    m, b, r, p, std_err = stats.linregress(x, y)
    plt.plot(x, m*x + b, '-', color='k')
    plt.xlabel(xlab, fontsize=9)
    plt.tick_params(axis='both', labelsize=6)
    plt.ylabel(ylab, fontsize=9)
    plt.title('$z$ = '+str(round(m,2)), fontsize=8)
    return fig


mydir = os.path.expanduser('~/GitHub/simplex')
tools = os.path.expanduser(mydir + "/tools")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df = df[df['ct'] > 200]

df2 = pd.DataFrame({'length' : df['length'].groupby(df['sim']).mean()})
df2['N'] = df['total.abundance'].groupby(df['sim']).mean()
df2['NS'] = df['avg.pop.size'].groupby(df['sim']).mean()
df2['Pdens'] = df['avg.pop.size'].groupby(df['sim']).min()/df2['length']**1

state = 'all'
df2['G'] = df[state+'.avg.per.capita.growth'].groupby(df['sim']).min()
df2['M'] = df[state+'.avg.per.capita.maint'].groupby(df['sim']).mean()
df2['D'] = df[state+'.avg.per.capita.active.dispersal'].groupby(df['sim']).mean()
df2['E'] = df[state+'.avg.per.capita.efficiency'].groupby(df['sim']).mean()
df2['size'] = df[state+'.size'].groupby(df['sim']).mean()

df2['B'] = df2['G'] * df2['size']
df2['MSB'] = df2['B']/df2['size']

#df2 = df2[np.log10(df2['Pdens']) != -1]
#df2 = df2[np.log10(df2['Pdens']) < -0.7]
#df2 = df2[np.log10(df2['size']) < 1]
#df2 = df2[np.log10(df2['Pdens']) > 0]


df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

fig = plt.figure()
xlab = r"$log_{10}$"+'(body size)'
ylab = r"$log_{10}$"+'(growth rate)'
fig = figplot(df2['size'], df2['B'], xlab, ylab, fig, 1)

xlab = r"$log_{10}$"+'(body size)'
ylab = r"$log_{10}$"+'(m.s. growth rate)'
fig = figplot(df2['size'], df2['MSB'], xlab, ylab, fig, 2)

xlab = r"$log_{10}$"+'(body size)'
ylab = r"$log_{10}$"+'(Pop. density)'
fig = figplot(df2['size'], df2['Pdens'], xlab, ylab, fig, 3)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/MetabolicScaling.png', dpi=200, bbox_inches = "tight")
