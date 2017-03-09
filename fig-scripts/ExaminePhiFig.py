from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
#from scipy import stats

def obs_pred_rsquare(obs, pred):

    """Determines the prop of variability in a data set accounted for by a model

    In other words, this determines the proportion of variation explained by
    the 1:1 line in an observed-predicted plot.

    """
    return 1 - sum((obs - pred) ** 2) / sum((obs - np.mean(obs)) ** 2)


mydir = os.path.expanduser('~/GitHub/residence-time')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width']})
df2['flow'] = df['flow.rate']
df2['tau'] = (df2['width']**3)/df2['flow']
df2['dil'] = np.log10(1/df2['tau'])

df2['N'] = df['total.abundance']
df2['S'] = df['species.richness']
df2['Prod'] = df['ind.production']
df2['Active'] = (100 - df['Percent.Dormant'])/100

df2['AvgG'] = df['active.avg.per.capita.growth']
df2['AvgDisp'] = df['active.avg.per.capita.active.dispersal']
df2['AvgRPF'] = df['dormant.avg.per.capita.RPF']
df2['AvgE'] = 1 - df['active.avg.per.capita.N.efficiency']
df2['AvgMaint'] = df['active.avg.per.capita.maint']
df2['MF'] = df['dormant.avg.per.capita.MF']/np.max(df['dormant.avg.per.capita.MF'])

E = 0.06



fs = 8
gd = 20
mct = 1
binz = 'log'
fig = plt.figure()
fig.add_subplot(3, 3, 1)

df2['P'] = np.log10((df2['AvgRPF']) * df2['AvgMaint'] * 1/df2['MF'] * df2['Active'])
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

r2 = obs_pred_rsquare(df2['dil'], df2['P'])
print 'r2:', r2, '\n'

plt.hexbin(df2['dil'], df2['P'], mincnt=mct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel(r"$persistence$", fontsize=fs+2)
plt.xlabel(r"$1/\tau$", fontsize=fs+2)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.plot([-8,0],[-8,0], 'k-')
plt.ylim(-8,0)
plt.xlim(-8,0)

fig.add_subplot(3, 3, 2)

df2['G'] = np.log10(df2['AvgG'] * df2['AvgDisp'] * df2['Active']) - E
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()
#r2 = obs_pred_rsquare(df2['dil'], df2['G'])
#print 'r2:', r2, '\n'

plt.hexbin(df2['dil'], df2['G'], mincnt=mct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel(r"$growth$", fontsize=fs+2)
plt.xlabel(r"$1/\tau$", fontsize=fs+2)
plt.tick_params(axis='both', which='major', labelsize=fs)

fig.add_subplot(3, 3, 3)

df2['P'] = df2['AvgMaint'] * df2['AvgRPF'] * df2['MF']
df2['G'] = df2['AvgG'] * df2['AvgDisp'] * df2['AvgE']

df2['phi'] = np.log10(df2['P'] / (df2['G'] + E))
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

plt.hexbin(df2['dil'], df2['phi'], mincnt=mct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel(r"$\phi$", fontsize=fs+2)
plt.xlabel(r"$1/\tau$", fontsize=fs+2)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.plot([-8,0],[-8,0], 'k-')
plt.ylim(-8,0)
plt.xlim(-8,0)


r2 = obs_pred_rsquare(df2['dil'], df2['phi'])
print 'r2:', r2


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.5, hspace=0.4)
plt.savefig(mydir + '/results/figures/Phi.png', dpi=200, bbox_inches = "tight")
plt.close()
