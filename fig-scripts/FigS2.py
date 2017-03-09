from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

from random import randint
from scipy import stats


p = 1
fr = 0.2
_lw = 0.5
w = 1
sz = 5


mydir = os.path.expanduser('~/GitHub/residence-time')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width']})
df2['flow'] = df['flow.rate']
df2['tau'] = np.log10((df['width'])/df2['flow'])

df2['N'] = np.log10(df['total.abundance']).groupby(df['ct']).median()
df2['D'] = np.log10(df['N.max']).groupby(df['ct']).median()
df2['S'] = np.log10(df['species.richness']).groupby(df['ct']).median()
df2['E'] = np.log10(df['simpson.e']).groupby(df['ct']).median()
df2['R'] = np.log10(df['logmod.skew']).groupby(df['ct']).median()
df2['sim'] = df['sim'].groupby(df['ct']).mean()

df2 = df2[df2['N'] > 1]
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()


metrics = ['Rarity, '+r'$log_{10}$',
        'Dominance, '+r'$log_{10}$',
        'Evenness, ' +r'$log_{10}$',
        'Richness, ' +r'$log_{10}$']


fig = plt.figure()
fs = 12 # font size used across figures
for index, metric in enumerate(metrics):
    fig.add_subplot(2, 2, index+1)

    Nlist = df2['N'].tolist()
    metlist = []
    if index == 0: metlist = df2['R'].tolist()
    elif index == 1: metlist = df2['D'].tolist()
    elif index == 2: metlist = df2['E'].tolist()
    elif index == 3: metlist = df2['S'].tolist()

    print len(df2['N']), len(metlist)
    df2['y'] = list(metlist)

    m, b, r, p, std_err = stats.linregress(Nlist, metlist)

    r2 = round(r**2,2)
    Int = round(b,2)
    Coef = round(m,2)

    gd = 15
    mct = 0
    plt.scatter(df2['N'], metlist, lw=_lw, color='0.2', s = sz)
    
    if index == 0:
        plt.text(0.7, 0.05, r'$rarity$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.7, -0.1,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(-0.9, 0.2)


    elif index == 1:
        plt.text(0.7, 2.4, r'$Nmax$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.7, 2.0,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')


    elif index == 2:
        plt.text(0.7, -1.3, r'$Ev$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.7, -1.6,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(-1.8, 0.1)

    elif index == 3:
        plt.text(0.7, 2.1, r'$S$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.7, 1.8,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(0.4, 2.5)


    plt.xlabel('$log$'+r'$_{10}$'+'($N$)', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/FigS2.png', dpi=600, bbox_inches = "tight")
#plt.show()
plt.close()
