from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from scipy import stats

mydir = os.path.expanduser('~/GitHub/simplex')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df = df[df['ct'] > 100]


df2 = pd.DataFrame({'length' : df['length']})
df2['N'] = np.log10(df['total.abundance']).groupby(df['sim']).median()
df2['D'] = np.log10(df['N.max']).groupby(df['sim']).median()
df2['S'] = np.log10(df['species.richness']).groupby(df['sim']).median()
df2['E'] = np.log10(df['simpson.e']).groupby(df['sim']).median()
df2['R'] = np.log10(df['logmod.skew']).groupby(df['sim']).median()

df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

metrics = ['Rarity, '+r'$log_{10}$', 'Dominance, '+r'$log_{10}$',
        'Evenness, ' +r'$log_{10}$', 'Richness, ' +r'$log_{10}$']


fig = plt.figure()
p, fr, _lw, w, sz, fs = 1, 0.2, 0.5, 1, 5, 12

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
        plt.text(1.9, -1.4, r'$rarity$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(1.9, -1.6,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(-1.8, 0)
        plt.xlim(0.5, 3.5)

    elif index == 1:
        plt.text(0.75, 3.0, r'$Nmax$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.75, 2.5,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(0.0, 3.5)
        plt.xlim(0.5, 3.5)

    elif index == 2:
        plt.text(0.75, -1.2, r'$Ev$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.75, -1.4,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(-1.6, 0.2)
        plt.xlim(0.5, 3.5)

    elif index == 3:
        plt.text(0.75, 1.6, r'$S$'+ ' = '+str(round(10**Int,2))+'*'+r'$N$'+'$^{'+str(round(Coef,2))+'}$', fontsize=fs-2, color='k')
        plt.text(0.75, 1.4,  r'$r^2$' + '=' +str(r2), fontsize=fs-2, color='k')
        plt.ylim(0.4, 1.8)
        plt.xlim(0.5, 3.5)

    plt.xlabel('$log$'+r'$_{10}$'+'($N$)', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling.png', dpi=600, bbox_inches = "tight")
plt.close()
