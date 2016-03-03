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

dat = dat[np.isfinite(dat['particle.tau'])]
dat = dat[np.isfinite(dat['individual.tau'])]
dat = dat[np.isfinite(dat['resource.concentration'])]
dat = dat[dat['particle.tau'] > 0]
dat = dat[dat['individual.tau'] > 0]

tau = np.log10((dat['width']*dat['height'])/dat['flow.rate']).tolist()
dat['tau'] = list(tau)
dat = dat[dat['tau'] < 6]


#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()

#### N vs. Tau #################################################################
fig.add_subplot(3, 3, 1)

xlab = 'Ecosystem ' + r"$\tau$" + ', ' + r"$log_{10}$"

dat2 = dat.groupby('RowID')

for group in dat2.groups.keys():

    df = dat2.get_group(group)
    Itau = np.log10(df['individual.tau']).tolist()
    if len(Itau) < 10:
        continue
    tau = df['tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(Itau, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=2, alpha=0.3)

Itau = np.log10(dat['individual.tau']).tolist()
tau = dat['tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(Itau, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

ylab = 'Organism ' + r"$\tau$" + ', ' + r"$log_{10}$"
plt.ylabel(ylab, fontsize=fs+6)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.ylim(-2, 5)
#plt.text(2, -1,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2,4)), fontsize=fs+2, color='k')

#### S vs. Tau #################################################################
fig.add_subplot(3, 3, 2)

for group in dat2.groups.keys():

    df = dat2.get_group(group)
    Ptau = df['particle.tau']
    if len(Ptau) < 10:
        continue

    Ptau = np.log10(df['particle.tau']).tolist()
    tau = df['tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(Ptau, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=2, alpha=0.3)

Ptau = np.log10(dat['particle.tau']).tolist()
tau = dat['tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(Ptau, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

#plt.ylim(0, 20)
plt.ylabel(r"$log_{10}$"+'(' + r"$\tau_{P}$" +')', fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(2, -1,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2,4)), fontsize=fs+2, color='k')

#### E vs. Tau #################################################################
fig.add_subplot(3, 3, 3)


for group in dat2.groups.keys():

    df = dat2.get_group(group)
    Itau = np.log10(df['individual.tau']).tolist()
    if len(Itau) < 10:
        continue

    tau = df['particle.tau'].tolist()

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess(Itau, tau)
    c = get_color()
    plt.plot(lowess[:, 0], lowess[:, 1], color = c, ls='-', lw=2, alpha=0.3)

Itau = np.log10(dat['individual.tau']).tolist()
tau = dat['particle.tau'].tolist()
lowess = statsmodels.nonparametric.smoothers_lowess.lowess(Itau, tau)
plt.plot(lowess[:, 0], lowess[:, 1], color = 'r', ls='-', lw=2, alpha=0.8)

ylab = 'Organism ' + r"$\tau$" + ', ' + r"$log_{10}$"
plt.ylabel(ylab, fontsize=fs+6)
plt.xlabel(xlab, fontsize=fs+6)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
#plt.savefig(mydir2 + 'Desktop/ResTime/figures/Fig2.png', dpi=600, bbox_inches = "tight")
plt.show()
