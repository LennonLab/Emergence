from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import os.path
import com1DCORE as c1dc
import cloud
import random
import sys

mydir = '/Users/lisalocey/Desktop/HydroBIDE/tools'
sys.path.append(mydir)
import metrics
import timeSeries

mydir = '/Users/lisalocey/Desktop/HydroBIDE/models/bide'
sys.path.append(mydir)
import bide

# First, some figure functions. These will be moved to their own module.


"""
REStimes = [1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
            3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0,
            6.5, 7.0, 8.0, 9.0, 10, 20.0, 40.0, 60.0, 80.0, 100.0, 200.0, 400.0,
            600.0, 800.0, 1000, 1500.0, 2000.0]
"""
"""
REStimes = [1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
            8.0, 10.0, 16.0, 20.0, 25.0, 32.0, 64.0, 128.0, 200.0,
            600.0, 1000.0, 1500.0, 2000.0]
"""


REStimes = [1.1, 3.0, 4.0, 6.0, 7.0, 8.0, 10.0, 16.0, 64.0, 128.0, 200.0,
            600.0, 1000.0, 1500.0]


                                        
dormancy = False
use_picloud = 'n'
var_to_vary = 'Volume'
varlist = range(1, 4)
           
Nlist = []   # mean total abundance
Plist = []   # mean productivity
AvgAblist = [] # # mean mean abundance (i.e. avg N/S)
Rlist = []  # mean total resources
TOlist = []  # mean biomass Turnover
Qlist = []   # mean cell quota

Slist = []  # richness
Evarlist = [] # evenness
Hlist = [] # Shannon's diversity
CTlist = [] # compositional turnover
RADLISTS = [] # list holding RADs and their info (V, REStime)

tauAtQmaxList = []
tauAtAvgMuMaxList = []
ObstauAtAvgMuMaxList = []

TauList = []
colors = []

for value in varlist:
    
    Qmax = 0.0      
    tauAtQmax = 0.0  
    
    REStimes2 = list(REStimes)
    time = 4000 # number of generations. Needs to be large enough to allow
    # the population to reach a relatively stable state
    
    V = random.choice([10, 20, 30, 40, 50])#, 60, 70, 80 ,90, 100])
    Pcons = random.choice([5, 10, 15, 20])
    Rcons = random.choice([5, 10, 15, 20])
    propQ = random.choice([0.1, 0.2, 0.3, 0.4, 0.5])
    lgp = random.choice([0.85, 0.9, 0.95, 0.99])
    
    optTau = random.choice([1.25, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0])
    a = optTau in REStimes2
    if a is False:
        REStimes2.append(optTau)
        REStimes2.sort()
        
    mean = 1.0/float(optTau)
    std = 0.0#float(random.randrange(0,51))/100 # standard deviation around the mean
    
    maint = float(random.randrange(1,20))/100
    maint = maint*mean      
                            
    random.seed()
    r = lambda: random.randint(0,255)
    color = '#%02X%02X%02X' % (r(),r(),r())
    colors.append(color)
    
    AvgN = [] # mean total abundance
    AvgP = [] # mean productivity
    
    AvgAb = [] # average mean abundance (i.e. avg N/S)
    AvgR = [] # mean total resources
    AvgTO = [] # mean biomass turnover
    AvgQ = [] # mean cell quota
   
    AvgS = [] # mean richness
    AvgEvar = [] # mean evenness
    AvgH = [] # mean diversity
    AvgCT = [] # mean compositional turnover
    
    AvgMaxU = []
    RADlists = []
    REStimes3 = []
    
    for ri, REStime in enumerate(REStimes2): 
        
        print value,'of', len(varlist),' : ', var_to_vary, V,' tau:',REStime,' OptTau:',optTau,
        avg_vals = []
        
        if use_picloud == 'n' and var_to_vary == 'Volume': 
            avg_vals = c1dc.hydrobide(V, REStime, Rcons, Pcons, propQ, time, dormancy, mean, std, lgp, maint)
              
        elif use_picloud == 'y' and var_to_vary == 'Volume':
            job_id = cloud.call(c1dc.hydrobide, V, REStime, Rcons, Pcons, propQ, time, dormancy, mean, std, lgp, maint, _type='m1')
            avg_vals = cloud.result(job_id)
        
        if len(avg_vals) == 1:
            print avg_vals[0]
        
        elif len(avg_vals) > 1:
            REStimes3.append(REStime)
            
            avgN, avgR, avgTO, avgQ, avgS, avgCT, avgEvar, avgH, RADinfo, avgab, avgMaxU, burnIn, lag, avgP = avg_vals
                
            AvgN.append(avgN)
            AvgP.append(avgP)
            
            AvgAb.append(avgab)
            AvgR.append(avgR)
            
            AvgTO.append(avgTO)
            AvgQ.append(avgQ)
        
            AvgS.append(avgS) # mean richness
            AvgEvar.append(avgEvar) # mean evenness
            
            AvgH.append(avgH) # mean diversity
            AvgCT.append(avgCT) # mean compositional turnover
            
            AvgMaxU.append(avgMaxU)
        
            if len(RADinfo) >= 3:
                RADlists.append(RADinfo)
            
            print '[Burn-in time:',burnIn,'] [suitable time lag:',lag,']', 
            print ' N:', round(avgN*V,2),'S:', round(avgS,2),'Q:',round(avgQ,2),
            
	    result1 = isinstance(avgN, float)
            result2 = isinstance(avgR, float)
            if result1 == True and result2 == True:
                print 'R:', round(avgR,2),
                if avgQ >= Qmax:
                    Qmax = avgQ
                    tauAtQmax = REStime
                
                print ''
            else: print ''
            
    tauAtQmaxList.append(tauAtQmax)
    tauAtAvgMuMaxList.append(mean**-1)
    MaxAvgMaxUptake = max(AvgMaxU)
        
    index = AvgMaxU.index(MaxAvgMaxUptake)
    ObstauAtAvgMuMaxList.append(REStimes3[index])
                                            
    Nlist.append(AvgN)
    Plist.append(AvgP)
    
    AvgAblist.append(AvgAb)
    Rlist.append(AvgR)
        
    TOlist.append(AvgTO)
    Qlist.append(AvgQ)
    Slist.append(AvgS)  # richness
        
    Evarlist.append(AvgEvar) # evenness
    Hlist.append(AvgH) # Shannon's diversity
    CTlist.append(AvgCT) # compositional turnover
    
    TauList.append(REStimes3)
    RADLISTS.append(RADlists)
    
    

""" A figure of response curves for random parameter combos """

fig = plt.figure()
fs = 10 # fontsize
depVars = [AvgAblist, Slist, CTlist, Evarlist]
beta = r"$\beta$"
y_labels = ['Avg abundance, N/S','Richness, S', beta+'-diversity' ,'Evenness, Evar']
x_label = r"$\tau$"

for i, depVar in enumerate(depVars):
    fig = figSet(fig, i, var_to_vary, varlist, depVar, TauList, x_label, y_labels[i], colors)

plt.savefig(mydir + 'figs/2x2Rand3.png', dpi=600, bbox_inches = 'tight', pad_inches=0.03)
#plt.show()


""" And now a figure of sweet predictions """

fig = plt.figure()
fs = 10 # fontsize
depVars = [Nlist, Rlist]
y_labels = ['Abundance, N/V','Resources, R/V']
x_label = r"$\tau$"

for i, depVar in enumerate(depVars):
    fig = figPred(fig, i, var_to_vary, varlist, depVar, TauList, x_label, y_labels[i], tauAtQmaxList, tauAtAvgMuMaxList, colors)

fig = ObsPredAvgMuMax(fig, ObstauAtAvgMuMaxList, tauAtAvgMuMaxList, colors)

plt.savefig(mydir + 'figs/sweet_predictions3.png', dpi=600, bbox_inches = 'tight', pad_inches=0.03)
plt.show()
