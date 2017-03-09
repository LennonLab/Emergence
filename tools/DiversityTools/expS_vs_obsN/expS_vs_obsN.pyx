from __future__ import division
import  matplotlib.pyplot as plt

import numpy as np
import random
from random import randrange, sample
import scipy as sc
from scipy import stats

import os
import sys
from scipy.stats.distributions import t

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

import itertools as it
import pandas as pd
from math import log10
import linecache

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")
mydir2 = os.path.expanduser("~/Desktop/")



""" This code creates a figures of three subplots. First, there is a plot of
observed species richness versus observed total abundance. Second, there is a
plot estimated richness versus number of samples. Third, there is a subplot of
expected richness versus dominance. """


def EstimateS(SiteList):
    m = len(SiteList)
    m_inf = 0
    SpDict = {}

    for site in SiteList:
        if min(site) <= 10: m_inf += 1

        for sp in site:
            if sp in SpDict:
                SpDict[sp] += 1
            else: SpDict[sp] = 1

    IncVals = SpDict.values()
    S = len(IncVals)
    qs = [0]*10

    for i, q in enumerate(qs):
        qs[i] = IncVals.count(i+1)

    # Chao2
    q1 = qs[0]
    q2 = qs[1]
    chao2 = S + (((m-1)/m) * ((q1*(q1-1)) / (2*(q2+1))))

    var = 'und'
    if q1 > 0 and q2 > 0:
        var = q2 * (0.5*(q1/q2)**2 + (q1/q2)**3 + 0.25*(q1/q2)**4)

    # ICE
    num = 0
    n_inf = 0
    for i, qk in enumerate(qs):
        num += (i+1)*i*qk
        n_inf += (i+1)*qk

    ci = 1 - (q1/n_inf)
    gamma = (sum(qs)/ci) * (m_inf/(m_inf-1)) * (num/(n_inf**2)) - 1
    cv = max(0, gamma)
    ice = (S-sum(qs)) + (sum(qs)/ci) + ((q1/ci) * cv)

    return [chao2, var, ice, S]





""" Read in data """

OrC = 'open'
path = mydir2 +'data/micro/EMP' + OrC
IN = path + '/EMP' + OrC + '-SADs.txt'


num_lines = sum(1 for line in open(IN))
lines = sample(range(1, num_lines+1), 1000)

Slist = []
Nlist = []
NofEMP = 0
SofEMP = 0

for i, line in enumerate(lines):

    data = linecache.getline(IN, line)
    SAD = eval(data)
    S = len(SAD)
    N = int(sum(SAD))

    if S < 2 or N < 10: continue

    Slist.append(float(np.log10(S)))
    Nlist.append(float(np.log10(N)))

d = pd.DataFrame({'N': list(Nlist)})
d['S'] = list(Slist)




""" Plot 1 """
fig = plt.figure()
fig.add_subplot(2, 3, 1)

ols = smf.ols('S ~ N', d).fit() # linear regression on log10(S) and log10(N)
ols_ci = ols.conf_int().ix['N'].tolist()
olsd = dict(a = ols.params['Intercept'], b = ols.params['N'], lb = ols_ci[0], ub = ols_ci[1])

print ols.summary()
x = np.arange(float(d.N.min()), float(d.N.max()), 0.1)
get_y = lambda a, b: a + b * x

y = get_y(olsd['a'], olsd['b'])
plt.plot(x, y, color='red', label='OLS')

plt.scatter(d.N, d.S, color='b', alpha=0.3)
plt.xlabel('log10(N)', fontsize=16)
plt.ylabel('log10(S)', fontsize=16)


""" Plot 2 """
fig.add_subplot(2, 3, 2)

ct = 0
numEMP = 0
numEMPopen = 0
RADs = []
OrC = 'open'

path = mydir2 +'data/micro/EMP'+OrC
IN = path + '/EMP' + OrC + '-SbyS.txt'

num_lines = sum(1 for line in open(IN))
print 'number of sites:', num_lines


SampSizes = [16, 24, 32, 48]# , 64, 96, 128, 182, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6142, 8192, 12288, num_lines]

AvgChao = []
AvgICE = []
AvgS = []

ciLowChao = []
ciHiChao = []
ciLowICE = []
ciHiICE = []
ciLowS = []
ciHiS = []

for n in SampSizes:

    Chao2s = []
    ICEs = []
    Ss = []

    for i in range(10):

        SiteData = []
        lines = random.sample(range(1, num_lines+1), n)
        for i, line in enumerate(lines):

            data = linecache.getline(IN, line)
            data = eval(data)

            slist = data[1]
            SiteData.append(slist)

        chao2, var, ice, S = EstimateS(SiteData)

        Chao2s.append(chao2)
        ICEs.append(ice)
        Ss.append(S)

    AvgChao.append(np.mean(Chao2s))
    AvgICE.append(np.mean(ICEs))
    AvgS.append(np.mean(Ss))

    ciLowChao.append(np.percentile(Chao2s, 2.5))
    ciHiChao.append(np.percentile(Chao2s, 97.5))
    ciLowICE.append(np.percentile(ICEs, 2.5))
    ciHiICE.append(np.percentile(ICEs, 97.5))
    ciLowS.append(np.percentile(Ss, 2.5))
    ciHiS.append(np.percentile(Ss, 97.5))

    print'n:',n, ' obsS:', AvgS[-1], ' Avg Chao2:',AvgChao[-1], ' Avg ICE:',AvgICE[-1]


fig = plt.figure()
fig.add_subplot(1, 1, 1)

plt.plot(SampSizes, AvgS, color='grey', ls='-', lw=2, label='obs S')
plt.plot(SampSizes, AvgChao, color='c', ls='-', lw=2, label='Chao')
plt.plot(SampSizes, AvgICE, color='m', ls='-', lw=2, label='ICE')

plt.fill_between(SampSizes, ciLowS, ciHiS, color='grey', alpha=0.3)
plt.fill_between(SampSizes, ciLowChao, ciHiChao, color='c', alpha=0.3)
plt.fill_between(SampSizes, ciLowICE, ciHiICE, color='m', alpha=0.3)

plt.xlabel('Number of samples', fontsize=16)
plt.ylabel('Number of species', fontsize=16)
leg = plt.legend(loc=4,prop={'size':14})
leg.draw_frame(False)


""" Plot 3 """
fig.add_subplot(2, 3, 3)

num_lines = sum(1 for line in open(IN))
lines = sample(range(1, num_lines+1), 1000)

obsS = []
obsN = []
obsD = []
expS = []

def SfromD(sad, D):

    slist = []
    for i in range(2):
        while max(sad) <= D:

            j = randrange(len(sad))
            ab1 = sad.pop(j)
            ab2 = int(round(10**(0.4 + 0.92*np.log10(ab1))))
            ab1 -= ab2
            sad.extend([ab1, ab2])

        slist.append(len(sad))
    return np.mean(slist)


for i, line in enumerate(lines):

    print i, num_lines
    if i > 100: break

    data = linecache.getline(IN, line)
    SAD = eval(data)
    S = len(SAD)
    N = int(sum(SAD))
    D = int(max(SAD))

    if S < 2 or N < 10: continue

    obsS.append(S)
    obsN.append(N)
    obsD.append(D)
    s = SfromD([N], D)
    expS.append(s)


lm = smf.ols('obsS ~ expS', d).fit() # linear regression on log10(S) and log10(N)
st, data, ss2 = summary_table(lm, alpha=0.05)
        # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
        # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
        # Residual, Std Error Residual, Student Residual, Cook's D

print lm.summary()
fittedvalues = data[:,2]
predict_mean_se = data[:,3]
mean_ci_low, mean_ci_upp = data[:,4:6].T

for i in range(len(obsS)):
    plt.scatter(obsS[i], expS[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Steelblue')

plt.fill_between(expS, mean_ci_low, mean_ci_upp, color='b', lw=0.0, alpha=0.3)

plt.scatter(d.N, d.S, color='b', alpha=0.3)
plt.xlabel('Predicted S, ' + 'r$log_{10}$', fontsize=16)
plt.ylabel('Observed S, ' + 'r$log_{10}$', fontsize=16)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig2-Open_NoMicrobe1s.png', dpi=600, bbox_inches = "tight")
plt.close()
