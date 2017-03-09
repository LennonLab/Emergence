from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/residence-time')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width'].groupby(df['ct']).mean()})
df2['flow'] = df['flow.rate'].groupby(df['ct']).mean()

p = 1
df2['tau'] = np.log10(df2['width']**p/df2['flow'])
#df2['tau'] = np.log10(df2['flow'])
#df2['tau'] = np.log10(df2['width'])

df2['R'] = df['resource.particles'].groupby(df['ct']).mean()
df2['N'] = df['total.abundance'].groupby(df['ct']).mean()
df2['Prod'] = df['ind.production'].groupby(df['ct']).mean()
df2['S'] = df['species.richness'].groupby(df['ct']).mean()
df2['E'] = np.log10(df['simpson.e'].groupby(df['ct']).mean())
df2['W'] = np.log10(df['Whittakers.turnover'].groupby(df['ct']).mean())
df2['Dorm'] = df['Percent.Dormant'].groupby(df['ct']).mean()

df2['G'] = df['active.avg.per.capita.growth'].groupby(df['ct']).mean()
df2['Maint'] = df['dormant.avg.per.capita.maint'].groupby(df['ct']).mean()
df2['Disp'] = df['active.avg.per.capita.active.dispersal'].groupby(df['ct']).mean()
df2['RPF'] = df['dormant.avg.per.capita.rpf'].groupby(df['ct']).mean()
df2['Eff'] = df['active.avg.per.capita.efficiency'].groupby(df['ct']).mean()
df2['MF'] = df['active.avg.per.capita.mf'].groupby(df['ct']).mean()


E = 0.1
df2['P'] = (E/df2['Maint']) * (E/df2['RPF']) * (df2['MF'])
df2['G'] = (E/df2['G']) * (E/df2['Disp']) * E*df2['Eff']
df2['phi'] = np.log10(df2['G'] * df2['P'])


df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()
df2['x'] = df2['phi'] * df2['tau']
xs = df2['x'].tolist()

#### plot figure ###############################################################

xlab = r"$log_{10}$"+'(' + r"$\tau$" +') * ' + r"$log_{10}$"+'(' + r"$\phi$" +')'
#xlab =  r"$\tau$" +' * ' + r"$\phi$"
fs = 6 # fontsize
fig = plt.figure()

mct = 1
binz = 'log'
ps = 120
ps2 = 5
gd = 30
a = 0.1

xl = -0.5
xh = 2

#### N vs. Tau #################################################################
Vs = df2['N'].tolist()
maxv = max(Vs)
i = Vs.index(maxv)
imax = xs[i]
fig.add_subplot(3, 3, 1)

#plt.hexbin(df2['x'], df2['N'], mincnt=mct, gridsize = gd, bins=binz, cmap=plt.cm.Greys_r)
plt.scatter(df2['x'], df2['N'], lw = 0.0, s = ps2, facecolors='k', alpha=0.8)
plt.scatter(imax, maxv, lw = 1.5, s = ps, facecolors='none', edgecolors='r')
plt.axvline(1, color='r', lw = 2, ls = '-')
plt.ylim(0, 2000)
plt.ylabel(r"$N$", fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### Prod vs. Tau #################################################################
Vs = df2['Prod'].tolist()
maxv = max(Vs)
i = Vs.index(maxv)
imax = xs[i]
fig.add_subplot(3, 3, 2)
plt.scatter(df2['x'], df2['Prod'], lw = 0.0, s = ps2, facecolors='k', alpha=0.8)
plt.scatter(imax, maxv, lw = 1.5, s = ps, facecolors='none', edgecolors='r')
plt.axvline(1, color='r', lw = 2, ls = '-')
plt.ylim(0, 80)
plt.ylabel("Productivity", fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### S vs. Tau #################################################################
Vs = df2['S'].tolist()
maxv = max(Vs)
i = Vs.index(maxv)
imax = xs[i]
fig.add_subplot(3, 3, 4)
plt.scatter(df2['x'], df2['S'], lw = 0.0, s = ps2, facecolors='k', alpha=0.8)
plt.scatter(imax, maxv, lw = 1.5, s = ps, facecolors='none', edgecolors='r')
plt.axvline(1, color='r', lw = 2, ls = '-')
plt.ylim(0, 80)
plt.ylabel(r"$S$", fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### E vs. Tau #################################################################
Vs = df2['E'].tolist()
maxv = min(Vs)
i = Vs.index(maxv)
imax = xs[i]
fig.add_subplot(3, 3, 5)
plt.scatter(df2['x'], df2['E'], lw = 0.0, s = ps2, facecolors='k', alpha=0.8)
plt.scatter(imax, maxv, lw = 1.5, s = ps, facecolors='none', edgecolors='r')
plt.axvline(1, color='r', lw = 2, ls = '-')
plt.ylim(-2, 0)
plt.ylabel(r"$E$", fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### W vs. Tau #################################################################
Vs = df2['W'].tolist()
maxv = max(Vs)
i = Vs.index(maxv)
imax = xs[i]
fig.add_subplot(3, 3, 7)
plt.scatter(df2['x'], df2['W'], lw = 0.0, s = ps2, facecolors='k', alpha=0.8)
plt.scatter(imax, maxv, lw = 1.5, s = ps, facecolors='none', edgecolors='r')
plt.axvline(1, color='r', lw = 2, ls = '-')
#plt.ylim(0, 2)
plt.ylabel(r"$W$", fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### Dorm vs. Tau #################################################################
Vs = df2['Dorm'].tolist()
maxv = min(Vs)
i = Vs.index(maxv)
imax = xs[i]
fig.add_subplot(3, 3, 8)
plt.scatter(df2['x'], df2['Dorm'], lw = 0.0, s = ps2, facecolors='k', alpha=0.8)
plt.scatter(imax, maxv, lw = 1.5, s = ps, facecolors='none', edgecolors='r')
plt.axvline(1, color='r', lw = 2, ls = '-')
plt.ylim(-5, 100)
plt.ylabel("%Dormant", fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Fig3.png', dpi=200, bbox_inches = "tight")
plt.close()
#plt.show()
