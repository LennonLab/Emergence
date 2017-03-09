from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

mydir = os.path.expanduser('~/GitHub/residence-time')
df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width'].groupby(df['ct']).mean()})

df2['sim'] = df['sim'].groupby(df['ct']).mean()
df2['flow'] = df['flow.rate'].groupby(df['ct']).mean()

p = 1
df2['tau'] = np.log10(df2['width']**p/df2['flow'])

df2['R'] = df['resource.particles'].groupby(df['ct']).mean()
df2['N'] = np.log10(df['total.abundance'].groupby(df['ct']).mean())
df2['Prod'] = df['ind.production'].groupby(df['ct']).mean()
df2['S'] = df['species.richness'].groupby(df['ct']).mean()
df2['E'] = df['simpson.e'].groupby(df['ct']).mean()
df2['W'] = np.log10(df['Whittakers.turnover'].groupby(df['ct']).mean())
df2['Dorm'] = df['Percent.Dormant'].groupby(df['ct']).mean()

df2['Grow'] = np.log10(df['active.avg.per.capita.growth'].groupby(df['ct']).mean())
df2['Maint'] = df['active.avg.per.capita.maint'].groupby(df['ct']).mean()
df2['Disp'] = df['active.avg.per.capita.active.dispersal'].groupby(df['ct']).min()
#df2['RPF'] = df['active.avg.per.capita.rpf'].groupby(df['ct']).max()
df2['Eff'] = df['active.avg.per.capita.efficiency'].groupby(df['ct']).mean()
#df2['MF'] = df['active.avg.per.capita.mf'].groupby(df['ct']).min()

#df2['Grow'] = df['dormant.avg.per.capita.growth'].groupby(df['ct']).mean()
#df2['Maint'] = np.log10(df['dormant.avg.per.capita.maint']).groupby(df['ct']).mean()
#df2['Disp'] = df['dormant.avg.per.capita.active.dispersal'].groupby(df['ct']).mean()
df2['RPF'] = df['dormant.avg.per.capita.rpf'].groupby(df['ct']).min()
#df2['Eff'] = df['dormant.avg.per.capita.efficiency'].groupby(df['ct']).mean()
df2['MF'] = df['dormant.avg.per.capita.mf'].groupby(df['ct']).min()

#df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()
#df2 = df2[df2['Prod'] < 50]
#df2 = df2[df2['N'] > 0]


#### plot figure ###############################################################
xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 6 # fontsize
fig = plt.figure()

gd = 25
mnct = 1
binz = 'log'
w = 1

#### N vs. Tau #################################################################
fig.add_subplot(3, 3, 1)

plt.hexbin(df2['tau'], df2['N'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel('N', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
#plt.ylim(1,2000)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.1, 2, 'A', color = 'y', fontweight='bold')

#### production vs. Tau ########################################################
fig.add_subplot(3, 3, 2)

plt.hexbin(df2['tau'], df2['Prod'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel('Productivity', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(4.2, -0.25, 'B', color = 'y', fontweight='bold')

#### S vs. Tau #################################################################
fig.add_subplot(3, 3, 4)

plt.hexbin(df2['tau'], df2['S'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel('S', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.1, 1.6, 'C', color = 'y', fontweight='bold')

#### E vs. Tau #################################################################
fig.add_subplot(3, 3, 5)

plt.hexbin(df2['tau'], df2['E'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel('Evenness', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(4.7, 0.3, 'D', color = 'y', fontweight='bold')

#### W vs. Tau #################################################################
ax5 = fig.add_subplot(3, 3, 7)

plt.hexbin(df2['tau'], df2['W'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel(r"$log_{10}$"+'(' + r"$\beta$" +')', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.1, -3.0, 'E', color = 'y', fontweight='bold')

#### dormancy vs. Tau ########################################################
fig.add_subplot(3, 3, 8)

plt.hexbin(df2['tau'], df2['Dorm'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet)
plt.ylabel('%Dormant', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(4.7, 0.2, 'F', color = 'y', fontweight='bold')

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Fig1-heat.png', dpi=200, bbox_inches = "tight")
plt.close()





#### plot figure ###############################################################
xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 6 # fontsize
fig = plt.figure()


#### AvgGrow vs. Tau #################################################################
fig.add_subplot(3, 3, 1)
plt.hexbin(df2['tau'], df2['Grow'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet, alpha = 1)
plt.ylabel('Specific growth rate', fontsize=fs+2)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(2, 0.15, 'A', color = 'y', fontweight='bold')
#plt.text(1.0, 1.05, 'Growth Syndrome', color = 'Crimson', fontsize = 10, fontweight='bold')


#### AvgActDisp vs. Tau #################################################################
fig.add_subplot(3, 3, 2)
plt.hexbin(df2['tau'], df2['Maint'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet, alpha = 1)
plt.ylabel('Maintenance energy, '+r"$log_{10}$", fontsize=fs+2)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.3, -4.6, 'C', color = 'y', fontweight='bold')
#plt.text(0.5, -0.38, 'Persistence Syndrome', color = 'Steelblue', fontsize = 10, fontweight='bold')

#### E vs. Tau #################################################################
fig.add_subplot(3, 3, 4)
df3 = df2[df2['tau'] > 1.5]
plt.hexbin(df3['tau'], df3['Disp'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet, alpha = 1)
plt.ylabel('Active disperal rate', fontsize=fs+2)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(2, 0.1, 'B', color = 'y', fontweight='bold')


#### AvgEff vs. Tau #################################################################
fig.add_subplot(3, 3, 5)
plt.hexbin(df2['tau'], df2['RPF'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet, alpha = 1)
plt.ylabel('Random resuscitation\nfrom dormancy, ' + r"$log_{10}$", fontsize=fs+2)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.3, -2, 'E', color = 'y', fontweight='bold')


#### AvgRPF vs. Tau #################################################################
fig.add_subplot(3, 3, 7)
df3 = df2[df2['tau'] > 1.5]
plt.hexbin(df3['tau'], df3['Eff'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet, alpha = 1)
plt.ylabel('Resource specialization', fontsize=fs+2)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.3, 0.099, 'D', color = 'y', fontweight='bold')


#### AvgRPF vs. Tau #################################################################
fig.add_subplot(3, 3, 8)
plt.hexbin(df2['tau'], df2['MF'], mincnt=mnct, gridsize = gd, bins=binz, cmap=plt.cm.jet, alpha = 1)
plt.ylabel('Decrease of maintenance\nenergy when dormant, ' + r"$log_{10}$", fontsize=fs+2)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(2.1, 4, 'F', color = 'y', fontweight='bold')

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.savefig(mydir + '/results/figures/Fig2-heat.png', dpi=200, bbox_inches = "tight")
plt.close()
