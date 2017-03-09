from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
from random import randint
import numpy as np
import os
from scipy import stats


p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 5


mydir = os.path.expanduser('~/GitHub/residence-time')
df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width'].groupby(df['ct']).mean()})
df2 = df2.replace([0], np.nan).dropna()

df2['sim'] = df['sim'].groupby(df['ct']).mean()
df2['flow'] = df['flow.rate'].groupby(df['ct']).mean()

df2['area'] = df['width'].groupby(df['ct']).mean()**1
df2['S'] = df['species.richness'].groupby(df['ct']).mean()

df2['N'] = df['total.abundance'].groupby(df['ct']).mean()
df2['NS'] = df['avg.pop.size'].groupby(df['ct']).mean()
df2['var'] = df['pop.var'].groupby(df['ct']).mean()

state = 'active'
df2['G'] = df[state+'.avg.per.capita.growth'].groupby(df['ct']).mean()
df2['M'] = df[state+'.avg.per.capita.maint'].groupby(df['ct']).mean()
df2['D'] = df[state+'.avg.per.capita.active.dispersal'].groupby(df['ct']).mean()
df2['E'] = df[state+'.avg.per.capita.efficiency'].groupby(df['ct']).mean()
df2['size'] = df[state+'.size'].groupby(df['ct']).mean()

df2['B'] = df2['M'] + df2['G'] + df2['E']#*df2['G']

df2 = df2[np.log10(df2['NS']) > 0.3]
#df2 = df2[np.log10(df2['area']) > 0]
#df2 = df2[np.log10(df2['var']) > 1]
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()


#### plot figure ###############################################################
fs = 6 # fontsize
fig = plt.figure()

#### N vs. Tau #################################################################


sims = list(set(df2['sim'].tolist()))
clrs = []
for sim in sims:
    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)
    clr = '#%02X%02X%02X' % (r1(),r2(),r3())
    clrs.append(clr)


Nlist = []
Vlist = []

Slist = []
Alist = []

Blist = []
Mlist = []


for i, sim in enumerate(sims):
    print 'sim:', sim

    df = df2[df2['sim'] == sim]

    Nlist.extend(np.log10(df['NS']).tolist())
    Vlist.extend(np.log10(df['var']).tolist())
    
    Alist.extend(np.log10(df['area']).tolist())
    Slist.extend(np.log10(df['S']).tolist())
    
    Mlist.extend(np.log10(df['size']).tolist())
    Blist.extend(np.log10(df['B']).tolist())
    

fig.add_subplot(3, 3, 1)
plt.scatter(Nlist, Vlist, lw=_lw, color=clrs, s = sz)
m, b, r, p, std_err = stats.linregress(Nlist, Vlist)
Nlist = np.array(Nlist)
plt.plot(Nlist, m*Nlist + b, '-')
xlab = r"$log_{10}$"+'(mean)'
ylab = r"$log_{10}$"+'(variance)'
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)
plt.ylabel(ylab, fontsize=fs+3)
plt.text(0.3, 2.5, 'slope = '+str(round(m,2)), fontsize=fs+2)
#plt.ylim(1.4, 2.8)
#plt.xlim(0.25, 0.65)

fig.add_subplot(3, 3, 2)
plt.scatter(Alist, Slist, lw=_lw, color=clrs, s = sz)
m, b, r, p, std_err = stats.linregress(Alist, Slist)
Alist = np.array(Alist)
plt.plot(Alist, m*Alist + b, '-')
xlab = r"$log_{10}$"+'(area)'
ylab = r"$log_{10}$"+'(S)'
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)
plt.ylabel(ylab, fontsize=fs+3)
plt.text(0, 2.3, 'slope = '+str(round(m,2)), fontsize=fs+2)
#plt.ylim(2.15, 2.35)
#plt.xlim(-.2, 1.4)

fig.add_subplot(3, 3, 3)
plt.scatter(Mlist, Blist, lw=_lw, color=clrs, s = sz)
m, b, r, p, std_err = stats.linregress(Mlist, Blist)
Mlist = np.array(Mlist)
plt.plot(Mlist, m*Mlist + b, '-')
xlab = r"$log_{10}$"+'(body size)'
ylab = r"$log_{10}$"+'(metabolic rate)'
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', labelsize=fs)
plt.ylabel(ylab, fontsize=fs+3)
plt.text(-0.38, -1.52, 'slope = '+str(round(m,2)), fontsize=fs+2)
#plt.ylim(-1.55, -1.3)
#plt.xlim(-.4, -0.05)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/MacroEcoPatterns2.png', dpi=200, bbox_inches = "tight")
plt.close()
