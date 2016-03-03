from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

from random import randint
import statsmodels
import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table


def get_color(): # FUNCTION TO ASSIGN COLORS TO Sp_

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    return '0.5'


mydir = os.path.expanduser('~/GitHub/ResTime')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

dat = pd.read_csv(mydir + '/results/simulated_data/examples/2015_12_08/SimData.csv')
dat = pd.DataFrame(dat)

dat = dat[np.isfinite(dat['total.abundance'])]
dat = dat[np.isfinite(dat['species.richness'])]
dat = dat[np.isfinite(dat['simpson.e'])]
dat = dat[np.isfinite(dat['Whittakers.turnover'])]
dat = dat[np.isfinite(dat['Jaccards.dissimilarity'])]
dat = dat[dat['total.abundance'] > 0]
dat = dat[dat['species.richness'] > 0]

tau = np.log10((dat['width']*dat['height'])/dat['flow.rate']).tolist()
dat['tau'] = list(tau)
dat = dat[dat['tau'] < 6]


#### plot figure ###############################################################
xlab = r"$log$"+'(' + r"$\tau$" +')'
fs = 8 # fontsize
fig = plt.figure()

#### N vs. Tau #################################################################
fig.add_subplot(2, 2, 1)

dat2 = dat.groupby('RowID')

for group in dat2.groups.keys():

    df = dat2.get_group(group)
    N = np.log10(df['total.abundance']).tolist()
    if len(N) < 10:
        continue
    tau = df['tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(N, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=1, alpha=0.3)

N = np.log10(dat['total.abundance']).tolist()
tau = dat['tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(N, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

plt.ylabel(r"$log$"+'(' + r"$N$" +')', fontsize=fs+6)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.ylim(-2, 5)
#plt.text(2, -1,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2,4)), fontsize=fs+2, color='k')

#### S vs. Tau #################################################################
fig.add_subplot(2, 2, 2)

for group in dat2.groups.keys():

    df = dat2.get_group(group)
    S = df['species.richness']
    if len(S) < 10:
        continue

    S = np.log10(df['species.richness']).tolist()
    tau = df['tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(S, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=1, alpha=0.3)

S = np.log10(dat['species.richness']).tolist()
tau = dat['tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(S, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

#plt.ylim(0, 20)
plt.ylabel(r"$log$"+'(' + r"$S$" +')', fontsize=fs+6)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(2, -1,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2,4)), fontsize=fs+2, color='k')

#### E vs. Tau #################################################################
fig.add_subplot(2, 2, 3)


for group in dat2.groups.keys():

    df = dat2.get_group(group)
    E = df['simpson.e'].tolist()
    if len(E) < 10:
        continue

    tau = df['tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(E, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=1, alpha=0.3)

E = dat['simpson.e'].tolist()
tau = dat['tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(E, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

plt.ylabel(r"$log$"+'(' + r"$Evenness$" +')', fontsize=fs+6)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### W vs. Tau #################################################################

fig.add_subplot(2, 2, 4)

for group in dat2.groups.keys():

    df = dat2.get_group(group)
    W = np.log10(df['Whittakers.turnover']).tolist()
    if len(W) < 10:
        continue

    tau = df['tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(W, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=1, alpha=0.3)

W = np.log10(dat['Whittakers.turnover']).tolist()
tau = dat['tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(W, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

plt.ylabel(r"$log$"+'(' + r"$N$" +')', fontsize=fs+6)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
#plt.savefig(mydir2 + 'Desktop/ResTime/figures/Fig2.png', dpi=600, bbox_inches = "tight")
plt.show()
