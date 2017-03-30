from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
import scipy as sc
from scipy import stats
import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/Emergence')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


def xfrm(X, _max): return -np.log10(_max - np.array(X))

def figplot(x, y, xlab, ylab, fig, n, binned = 1):

    '''main figure plotting function'''

    fig.add_subplot(2, 2, n)
    y2 = list(y)
    x2 = list(x)

    if binned == 1:
        X, Y = (np.array(t) for t in zip(*sorted(zip(x2, y2))))
        Xi = xfrm(X, max(X)*1.05)
        bins = np.linspace(np.min(Xi), np.max(Xi)+1, 100)
        ii = np.digitize(Xi, bins)
        y2 = np.array([np.mean(Y[ii==i]) for i in range(1, len(bins)) if len(Y[ii==i]) > 0])
        x2 = np.array([np.mean(X[ii==i]) for i in range(1, len(bins)) if len(X[ii==i]) > 0])

    d = pd.DataFrame({'x': list(x2)})
    d['y'] = list(y2)
    f = smf.ols('y ~ x', d).fit()

    m, b, r, p, std_err = stats.linregress(x2, y2)

    st, data, ss2 = summary_table(f, alpha=0.05)
    fitted = data[:,2]
    mean_ci_low, mean_ci_upp = data[:,4:6].T
    ci_low, ci_upp = data[:,6:8].T

    x2, y2, fitted, ci_low, ci_upp = zip(*sorted(zip(x2, y2, fitted, ci_low, ci_upp)))

    if n == 1: lbl = r'$rarity$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'+r'$r^2$' + '=' +str(round(r**2,2))
    elif n == 2: lbl = r'$Nmax$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'+r'$r^2$' + '=' +str(round(r**2,2))
    elif n == 3: lbl = r'$Ev$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'+ r'$r^2$' + '=' + str(round(r**2,2))
    elif n == 4: lbl = r'$S$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'+r'$r^2$' + '=' +   str(round(r**2,2))
    plt.scatter(x2, y2, color = 'SkyBlue', alpha= 1 , s = 12, linewidths=0.5, edgecolor='Steelblue', label=lbl)
    plt.legend(loc='best', fontsize=6)

    plt.fill_between(x2, ci_upp, ci_low, color='b', lw=0.1, alpha=0.15)
    plt.plot(x2, fitted,  color='b', ls='--', lw=1.0, alpha=0.9)
    plt.xlabel(xlab, fontsize=9)
    plt.ylabel(ylab, fontsize=9)
    plt.tick_params(axis='both', labelsize=6)
    plt.xlim(0.9*min(x2), 1.1*max(x2))
    plt.ylim(min(ci_low), max(ci_upp))
    return fig


df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'length' : df['length']})
df2['N'] = np.log10(df['total.abundance'].groupby(df['sim']).max())
df2['D'] = np.log10(df['N.max'].groupby(df['sim']).max())
df2['S'] = np.log10(df['species.richness'].groupby(df['sim']).max())
df2['E'] = np.log10(df['simpson.e'].groupby(df['sim']).min())
df2['R'] = np.log10(df['logmod.skew'].groupby(df['sim']).max())

df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()
fig = plt.figure()

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Rarity, '+r'$log_{10}$'
fig = figplot(df2['N'], df2['R'], xlab, ylab, fig, 1)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Dominance, '+r'$log_{10}$'
fig = figplot(df2['N'], df2['D'], xlab, ylab, fig, 2)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Evenness, ' +r'$log_{10}$'
fig = figplot(df2['N'], df2['E'], xlab, ylab, fig, 3)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Richness, ' +r'$log_{10}$'
fig = figplot(df2['N'], df2['S'], xlab, ylab, fig, 4)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling.png', dpi=600, bbox_inches = "tight")
plt.close()
