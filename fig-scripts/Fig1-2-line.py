from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
from random import randint
import numpy as np
import os
import statsmodels.api as sm
import sys

p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 5


mydir = os.path.expanduser('~/GitHub/residence-time')
df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width'].groupby(df['ct']).mean()})

df2['sim'] = df['sim'].groupby(df['ct']).mean()
df2['flow'] = df['flow.rate'].groupby(df['ct']).mean()
df2['tau'] = np.log10(df2['width']**p/df2['flow'])

df2['R'] = df['resource.particles'].groupby(df['ct']).mean()
df2['N'] = df['total.abundance'].groupby(df['ct']).mean()
df2['Prod'] = df['ind.production'].groupby(df['ct']).mean()
df2['S'] = df['species.richness'].groupby(df['ct']).mean()
df2['E'] = df['simpson.e'].groupby(df['ct']).mean()
df2['W'] = df['Whittakers.turnover'].groupby(df['ct']).mean()
df2['Dorm'] = df['Percent.Dormant'].groupby(df['ct']).mean()

df2['Grow'] = df['active.avg.per.capita.growth'].groupby(df['ct']).mean()
df2['Maint'] = np.log10(df['active.avg.per.capita.maint'].groupby(df['ct']).mean())
df2['Disp'] = df['active.avg.per.capita.active.dispersal'].groupby(df['ct']).mean()
df2['RPF'] = df['dormant.avg.per.capita.rpf'].groupby(df['ct']).mean()
df2['Eff'] = df['active.avg.per.capita.efficiency'].groupby(df['ct']).mean()
df2['MF'] = df['active.avg.per.capita.mf'].groupby(df['ct']).min()

#df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

#### plot figure ###############################################################
xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 6 # fontsize

#### plot figure ###############################################################
xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 6 # fontsize
fig = plt.figure()

#### N vs. Tau #################################################################
ax1 = fig.add_subplot(3, 3, 1)
ax2 = fig.add_subplot(3, 3, 2)
ax3 = fig.add_subplot(3, 3, 4)
ax4 = fig.add_subplot(3, 3, 5)
ax5 = fig.add_subplot(3, 3, 7)
ax6 = fig.add_subplot(3, 3, 8)


sims = list(set(df2['sim'].tolist()))
clrs = []
for sim in sims:
    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)
    clr = '#%02X%02X%02X' % (r1(),r2(),r3())
    clrs.append(clr)


for i, sim in enumerate(sims):
    print 'sim:', sim

    df = df2[df2['sim'] == sim]

    x = df['tau'].tolist()
    y = df['Grow'].tolist()
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax1.plot(x, y, lw=_lw, color=clrs[i])
    ax1.set_xlabel(xlab, fontsize=fs+3)
    ax1.tick_params(axis='both', labelsize=fs)
    ax1.set_ylabel('Growth rate', fontsize=fs+3)
    #ax1.set_xlim(0, 5)

    x = df['tau'].tolist()
    y = df['Maint'].tolist()
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax2.plot(x, y, lw=_lw, color=clrs[i])
    ax2.set_xlabel(xlab, fontsize=fs+3)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.set_ylabel('Maintenance energy', fontsize=fs+3)
    #ax2.set_xlim(0, 5)

    x = df['tau']
    y = df['Disp']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax3.plot(x, y, lw=_lw, color=clrs[i])
    ax3.set_xlabel(xlab, fontsize=fs+3)
    ax3.tick_params(axis='both', labelsize=fs)
    ax3.set_ylabel('Active disperal rate', fontsize=fs+2)
    #ax3.set_xlim(0, 5)

    x = df['tau']
    y = df['RPF']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax4.plot(x, y, lw=_lw, color=clrs[i])
    ax4.set_xlabel(xlab, fontsize=fs+3)
    ax4.tick_params(axis='both', labelsize=fs)
    ax4.set_ylabel('Random resuscitation\nfrom dormancy', fontsize=fs+2)
    #ax4.set_xlim(0, 5)

    x = df['tau']
    y = df['Eff']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax5.plot(x, y, ls='-', lw=_lw, color=clrs[i])

    ax5.set_xlabel(xlab, fontsize=fs+3)
    ax5.tick_params(axis='both', labelsize=fs)
    ax5.set_ylabel('Resource specialization', fontsize=fs+2)
    #ax5.set_xlim(0, 5)

    x = df['tau']
    y = df['MF']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax6.plot(x, y, lw=_lw, color=clrs[i])
    ax6.set_xlabel(xlab, fontsize=fs+3)
    ax6.tick_params(axis='both', labelsize=fs)
    ax6.set_ylabel('Decrease of maintenance\nenergy when dormant', fontsize=fs+2)
    #ax6.set_xlim(0, 5)
    #ax6.set_ylim(0, 40)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Fig2-line.png', dpi=200, bbox_inches = "tight")
plt.close()





fig = plt.figure()
#### N vs. Tau #################################################################
ax1 = fig.add_subplot(3, 3, 1)
ax2 = fig.add_subplot(3, 3, 2)
ax3 = fig.add_subplot(3, 3, 4)
ax4 = fig.add_subplot(3, 3, 5)
ax5 = fig.add_subplot(3, 3, 7)
ax6 = fig.add_subplot(3, 3, 8)


for i, sim in enumerate(sims):
    print 'sim:', sim

    df = df2[df2['sim'] == sim]

    x = df['tau']
    y = df['N']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax1.plot(x, y, lw=_lw, color=clrs[i])
    ax1.set_xlabel(xlab, fontsize=fs+3)
    ax1.tick_params(axis='both', labelsize=fs)
    ax1.set_ylabel('N', fontsize=fs+3)
    #ax1.set_ylim(0, 1300)
    #ax1.set_xlim(0, 5)

    x = df['tau'].tolist()
    y = df['Prod'].tolist()
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax2.plot(x, y, lw=_lw, color=clrs[i])
    ax2.set_xlabel(xlab, fontsize=fs+3)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.set_ylabel('Productivity', fontsize=fs+3)
    #ax2.set_ylim(0, 35)
    #ax2.set_xlim(0, 5)

    x = df['tau']
    y = df['S']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax3.plot(x, y, lw=_lw, color=clrs[i])
    ax3.set_xlabel(xlab, fontsize=fs+3)
    ax3.tick_params(axis='both', labelsize=fs)
    ax3.set_ylabel('S', fontsize=fs+3)
    #ax3.set_ylim(0, 50)
    #ax3.set_xlim(0, 5)

    x = df['tau']
    y = df['E']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax4.plot(x, y, lw=_lw, color=clrs[i])
    ax4.set_xlabel(xlab, fontsize=fs+3)
    ax4.tick_params(axis='both', labelsize=fs)
    ax4.set_ylabel('Evenness', fontsize=fs+3)
    ax4.set_ylim(0, 1)
    #ax4.set_xlim(0, 5)

    x = df['tau']
    y = np.log10(df['W'])
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax5.plot(x, y, lw=_lw, color=clrs[i])
    ax5.set_xlabel(xlab, fontsize=fs+3)
    ax5.tick_params(axis='both', labelsize=fs)
    ax5.set_ylabel(r"$log_{10}$"+'(' + r"$\beta$" +')', fontsize=fs+3)
    ax5.set_ylim(-2, 0.5)
    #ax5.set_xlim(0, 5)

    x = df['tau']
    y = df['Dorm']
    lowess = sm.nonparametric.lowess(y, x, frac=fr)
    x, y = lowess[:, 0], lowess[:, 1]
    ax6.plot(x, y, lw=_lw, color=clrs[i])
    ax6.set_xlabel(xlab, fontsize=fs+3)
    ax6.tick_params(axis='both', labelsize=fs)
    ax6.set_ylabel('%Dormant', fontsize=fs+3)
    ax6.set_ylim(0, 100)
    #ax6.set_xlim(0, 5)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Fig1-line.png', dpi=200, bbox_inches = "tight")
plt.close()
