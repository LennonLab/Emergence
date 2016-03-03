from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table


mydir = os.path.expanduser('~/GitHub/ResTime')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

dat = pd.read_csv(mydir + '/results/simulated_data/examples/2015_12_08/SimData.csv')
#dat = pd.read_csv(mydir + '/results/simulated_data/examples/2015_12_27/SimData.csv')
dat = dat[np.isfinite(dat['individual.tau'])]
dat = dat[np.isfinite(dat['particle.tau'])]

dat = dat[dat['particle.tau'] > 0]
dat = dat[dat['individual.tau'] > 0]

dat = dat[dat['barriers'] == 0]

dat['Itau'] = dat['individual.tau']
dat['Ptau'] = np.log10(dat['particle.tau'])
dat['Etau'] = np.log10(dat['width']/dat['flow.rate'])


Etau = dat['Etau'].tolist()
Itau = dat['Itau'].tolist()
Ptau = dat['Ptau'].tolist()

#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()

#### Itau vs. Tau #################################################################
fig.add_subplot(3, 3, 1)

xlab = 'Ecosystem ' + r"$\tau$" + ', ' + r"$log_{10}$"
f2 = smf.ols('Itau ~ Etau + I(Etau ** 2.0)', dat).fit()
#print f2.summary(),'\n\n'
#f2 = smf.ols('Itau ~ tau', d).fit()

a, b, c =  f2.params
p1, p2, p3 = f2.pvalues
r2 = round(f2.rsquared, 2)

st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]
pred_mean_se = data[:,3]
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

tau2, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp = zip(*sorted(zip(Etau, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp)))
plt.scatter(Etau, Itau, color = 'b', alpha = 0.4 , s = 5, linewidths = 0.0, edgecolor = 'k')
plt.fill_between(tau2, pred_ci_low, pred_ci_upp, color='r', lw=0.0, alpha=0.1)
plt.fill_between(tau2, pred_mean_ci_low, pred_mean_ci_upp, color='r', lw=0.0, alpha=0.3)
plt.plot(tau2, fitted, color = 'r', ls='--', lw=0.5, alpha=0.9)

ylab = 'Organism ' + r"$\tau$" + ', ' + r"$log_{10}$"
plt.ylabel(ylab, fontsize=fs+5)
#plt.ylim(0, 7)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

plt.text(4, 0.5,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p3, 3)), fontsize=fs+1, color='k')

#### Particle tau vs. Tau #################################################################

fig.add_subplot(3, 3, 2)
f2 = smf.ols('Ptau ~ Etau', dat).fit()
#print f2.summary(),'\n\n'
#f2 = smf.ols('Ptau ~ tau', d).fit()
#print f2.summary()

a, b  =  f2.params
p1, p2 = f2.pvalues
r2 = round(f2.rsquared, 2)

st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]
pred_mean_se = data[:,3]
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

tau2, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp = zip(*sorted(zip(Etau, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp)))
plt.scatter(Etau, Ptau, color = 'b', alpha = 0.4 , s = 5, linewidths = 0.0, edgecolor = 'k')
plt.fill_between(tau2, pred_ci_low, pred_ci_upp, color='r', lw=0.0, alpha=0.1)
plt.fill_between(tau2, pred_mean_ci_low, pred_mean_ci_upp, color='r', lw=0.0, alpha=0.3)
plt.plot(tau2, fitted, color = 'r', ls='--', lw=0.5, alpha=0.9)

ys = [0,3]
xs = [0,3]
plt.plot(xs, ys, lw=1.0, ls='--', color='0.3')

plt.ylabel(r"$log_{10}$"+'(' + r"$\tau_{P}$" +')', fontsize=fs+5)
#plt.ylim(0, 3)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

plt.text(3.5, 0.5,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p3, 3)), fontsize=fs+1, color='k')

#### Ptau vs. ITau #################################################################
fig.add_subplot(3, 3, 3)

xlab = 'Particle ' + r"$\tau$" + ', ' + r"$log_{10}$"

f2 = smf.ols('Itau ~ Ptau + I(Ptau ** 2.0)', dat).fit()
print f2.summary(),'\n\n'

#f2 = smf.ols('Itau ~ Ptau', d).fit()
#print f2.summary()
#sys.exit()

a, b, c  =  f2.params
p1, p2, p3 = f2.pvalues
r2 = round(f2.rsquared, 2)

st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]
pred_mean_se = data[:,3]
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

tau2, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp = zip(*sorted(zip(Ptau, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp)))
plt.scatter(Ptau, Itau, color = 'b', alpha = 0.4 , s = 5, linewidths = 0.0, edgecolor = 'k')
plt.fill_between(tau2, pred_ci_low, pred_ci_upp, color='r', lw=0.0, alpha=0.1)
plt.fill_between(tau2, pred_mean_ci_low, pred_mean_ci_upp, color='r', lw=0.0, alpha=0.3)
plt.plot(tau2, fitted, color = 'r', ls='--', lw=0.5, alpha=0.9)

ylab = 'Organism ' + r"$\tau$" + ', ' + r"$log_{10}$"
plt.ylabel(ylab, fontsize=fs+5)
#plt.ylim(0.5, 2.5)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)

plt.text(0.1, 2.2,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2,3)), fontsize=fs+1, color='k')


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
#plt.savefig(mydir2 + 'Desktop/ResTime/figures/Fig1.png', dpi=600, bbox_inches = "tight")
plt.show()
