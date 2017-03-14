from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy import stats

mydir = os.path.expanduser('~/GitHub/simplex')
tools = os.path.expanduser(mydir + "/tools")


p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 1

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width'].groupby(df['sim']).mean()})
#df2 = df2.replace([0], np.nan).dropna()

df2['sim'] = df['sim'].groupby(df['sim']).mean()
df2['flow'] = df['flow.rate'].groupby(df['sim']).mean()

df2['tau'] = np.log10(df2['width']**p/df2['flow'])

df2['N'] = np.log10(df['total.abundance'].groupby(df['sim']).mean())
df2['NS'] = np.log10(df['avg.pop.size'].groupby(df['sim']).mean())
df2['Pdens'] = np.log10(df['avg.pop.size'].groupby(df['sim']).mean()/df2['width']**3)

state = 'all'
df2['G'] = df[state+'.avg.per.capita.growth'].groupby(df['sim']).max()
df2['M'] = df[state+'.avg.per.capita.maint'].groupby(df['sim']).mean()
df2['D'] = df[state+'.avg.per.capita.active.dispersal'].groupby(df['sim']).mean()
df2['E'] = df[state+'.avg.per.capita.efficiency'].groupby(df['sim']).mean()
df2['size'] = np.log10(df[state+'.size'].groupby(df['sim']).mean())

df2['B'] = np.log10(df2['G'])
df2['MSB'] = np.log10(df2['G']/df[state+'.size'].groupby(df['sim']).mean())

df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

#### plot figure ###############################################################
fs = 6 # fontsize
fig = plt.figure()

#### Metabolic rate vs. Body size ##############################################
fig.add_subplot(2, 2, 1)
plt.scatter(df2['size'], df2['B'], lw=_lw, color='0.2', s = sz)
m, b, r, p, std_err = stats.linregress(df2['size'], df2['B'])
print 'MTE:', m
plt.plot(df2['size'], m*df2['size'] + b, '-', color='k')
xlab = r"$log_{10}$"+'(body size)'
ylab = r"$log_{10}$"+'(growth rate)'
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)
plt.ylabel(ylab, fontsize=fs+3)
plt.text(1.0, -2.4, '$z$ = '+str(round(m,2)), fontsize=fs+2)
plt.ylim(-2.6, -0.6)
plt.xlim(0, 2.0)


#### N vs. Tau #################################################################
fig.add_subplot(2, 2, 2)
plt.scatter(df2['size'], df2['MSB'], lw=_lw, color='0.2', s = sz)
m, b, r, p, std_err = stats.linregress(df2['size'], df2['MSB'])
print 'MTE (msb):', m
#Mlist = np.array(Mlist)
plt.plot(df2['size'], m*df2['size'] + b, '-', color='k')
xlab = r"$log_{10}$"+'(body size)'
ylab = r"$log_{10}$"+'(m.s. growth rate)'
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)
plt.ylabel(ylab, fontsize=fs+1)
plt.text(1.0, -1.7, '$z$ = '+str(round(m,2)), fontsize=fs+2)
plt.ylim(-3.2, -1.4)
plt.xlim(0, 2.0)



fig.add_subplot(2, 2, 3)
plt.scatter(df2['size'], df2['Pdens'], lw=_lw, color='0.2', s = sz)
m, b, r, p, std_err = stats.linregress(df2['size'], df2['Pdens'])
print 'Pop density:', m
Mlist = np.array(df2['size'])
plt.plot(df2['size'], m*df2['size'] + b, '-', color='k')
xlab = r"$log_{10}$"+'(Body size)'
ylab = r"$log_{10}$"+'(Pop. density)'
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)
plt.ylabel(ylab, fontsize=fs+3)
plt.text(1.0, -0.4, '$z$ = '+str(round(m,2)), fontsize=fs+2)
plt.ylim(-3.0, 0.0)
plt.xlim(0, 2.0)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/MetabolicScaling.png', dpi=200, bbox_inches = "tight")
plt.close()
