from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats

mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/Desktop/")

sys.path.append(mydir + "/tools/GenAnalysis/tools/")
import macroecotools
import macroeco_distributions
import predRADs
import mete
import pln

sys.path.append(mydir + "tools/feasible_functions")
import feasible_functions as ff


datasets = []
for name in os.listdir(mydir2 +'data/micro'):
        datasets.append([name, 'micro'])

for name in os.listdir(mydir2 +'data/macro'):
        datasets.append([name, 'macro'])


ct = 0
numMicros = 0
numMacros = 0
numEMP = 0
numEMPopen = 0

models = ['METE', 'MaxEntGeomSeries', 'IntegerPartitions', 'Zipf', 'PoissonLogNormal']

for dataset in datasets:

    name = dataset[0] # name of dataset
    if name == '.DS_Store': continue

    kind = dataset[1] # micro or macro

    for model in models:
        OUT = open(mydir2 + 'data/'+kind+'/'+name+'/'+name+'-'+model+'_SADMetricData_NoMicrobe1s.txt','w+')
        #OUT = open(mydir2 + 'data/'+kind+'/'+name+'/'+name+'-'+model+'_SADMetricData.txt','w+')
        RADs = []

        if kind == 'macro':
            path = kind
            RADs = ff.get_SADs(mydir2 +'/data/'+path, name)
            print 'macro', name, len(RADs)

        if kind == 'micro':
            path = kind
            RADs = ff.GetSADsFromBiom_labeled(mydir2 +'data/'+path+'/EMPopen/', name)
            print 'micro', name, len(RADs)

        ct = 0
        numRADs = len(RADs)
        for RAD in RADs:

            if kind == 'micro':
                RAD = list([x for x in RAD if x > 1])

            elif kind == 'macro':
	            RAD = list([x for x in RAD if x > 0])


            N = sum(RAD)
            S = len(RAD)

            if S < 2 or N < 10: continue

            if model == 'MaxEntGeomSeries':
                # Geometric Series for N and S
                RAD = mete.get_mete_sad_geom(S, N) # False mean no zeros allowed

            elif model == 'METE':
                # Log-series via METE, i.e. mle for N/S
                mete = mete.get_mete_rad(S, N)
                RAD = logSeries[0]

            elif model == 'IntegerPartitionsNS':
                # Feasible set of partitions, based on N and S, sensu Locey and White 2013
	        ObsPred = open(mydir2 + 'data/'+kind+'/'+name+'/'+name+'-'+model+'_SADMetricData_NoMicrobe1s.txt','r')
	        RAD =

            elif model == 'IntegerPartitionsN':
                # Feasible set of partitions, based on N
                RAD =

            elif model == 'Zipf':
                # Zipf
                RAD =

            elif model == 'PoissonLogNormal':
                # Poisson Log-normal
                RAD = pln.get_rad_from_obs(RAD, 'pln')

            elif model == 'NegativeBinomial':
                RAD = pln.get_rad_from_obs(RAD, 'negbin')

            Evar = ff.e_var(RAD)
            ESimp = ff.simpsons_evenness(RAD)
            ENee = ff.NHC_evenness(RAD)
            EPielou = ff.pielous_evenness(RAD)

            EHeip = ff.Heips_evenness(RAD)
            EQ = ff.EQ_evenness(RAD)

            BP = ff.Berger_Parker(RAD)
            SimpDom = ff.simpsons_dom(RAD)
            perc_ones = ff.Singletons(RAD)

            Nmax = max(RAD)

            rareRel = ff.rarityRel(RAD)
            rareOnes = ff.rarityOnes(RAD)

            skew = stats.skew(RAD)

            ct+=1

            print>>OUT, name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew
            print kind, name, numRADs - ct

        OUT.close()
