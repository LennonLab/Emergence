from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats

import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd

########### PATHS ##############################################################

tools = os.path.expanduser("~/tools")
sys.path.append(tools + "/macroeco_distributions")
import macroeco_distributions as md
sys.path.append(tools + "/distributions")
import distributions as dist
sys.path.append(tools + "/macroecotools")
import macroecotools
sys.path.append(tools + "/partitions")
import partitions as parts
sys.path.append(tools + "/FeasibleFunctions")
import feasible_functions as ff
sys.path.append(tools + "/metrics")
import metrics
sys.path.append(tools + "/mete")
import mete
sys.path.append(tools + "/pln")
import pln

data = os.path.expanduser("~/data")
mydir = os.path.expanduser("~/GitHub/rare-bio")
#macroeco_dir = os.path.expanduser("~/GitHub/macroecotools")

########### END PATHS ##########################################################


########### FUNCTIONS ##########################################################

def import_obs_pred_data(input_filename):
    # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "S15, f8, f8",
    names = ['site','obs','pred'], delimiter = " ")
    return data

########### END FUNCTIONS ######################################################


########### GET DATA ###########################################################

datasets = []
for name in os.listdir(data +'/micro'):
        datasets.append([name, 'micro'])

for name in os.listdir(data +'/macro'):
        datasets.append([name, 'macro'])

########### END GET DATA #######################################################



########### DECLARE VARIABLES ##################################################

ct = 0
numMicros = 0
numMacros = 0
numEMP = 0
numEMPopen = 0

models = ['METE', 'MaxEntGeomD']
#models = ['IntPartsNS', 'IntPartsN']
########### END DECLARE VARIABLES ##############################################


########### LOOP THROUGH DATA ##################################################

for dataset in datasets:

    name = dataset[0] # name of dataset
    if name == '.DS_Store': continue

    kind = dataset[1] # micro or macro

    for model in models:
        OUT = open(data + '/'+kind+'/'+name+'/'+name+'-'+model+'_SADMetricData_NoMicrobe1s.txt','w+')
        #OUT = open(mydir2 + 'data/'+kind+'/'+name+'/'+name+'-'+model+'_SADMetricData.txt','w+')
        RADs = []

        if name == 'EMPopen' or name == 'EMPclosed':
            path = data +'/micro/'
            RADs = metrics.GetSADsFromBiom_labeled(path, name)

        else: RADs = metrics.get_SADs(data +'/'+ kind, name)

        if len(RADs) > 10:
            RADs = random.sample(RADs,10)

        if kind == 'micro':

            print 'micro', name, len(RADs)

        ct = 0
        numRADs = len(RADs)
        for RAD in RADs:

            if kind == 'micro':
                RAD = list([int(x) for x in RAD if x > 1])

            elif kind == 'macro':
                RAD = list([int(x) for x in RAD if x > 0])

            N = int(sum(RAD))
            S = int(len(RAD))

            if S < 10 or max(RAD) < 2: continue

            if model == 'MaxEntGeomD':
                # Geometric Series for N and S
                result = mete.get_mete_sad_geom(S, N) # False mean no zeros allowed
                RAD = result[0]

            elif model == 'METE':
                # Log-series via METE, i.e. mle for N/S
                result = mete.get_mete_rad(S, N)
                RAD = result[0]

            elif model == 'IntPartsNS':
                # Feasible set of partitions, based on N and S, sensu Locey and White 2013
                if kind == 'macro' or N > 10**4: continue
                partitions = parts.rand_partitions(N, S, 100)
                unique_RADs = [list(x) for x in set(tuple(x) for x in partitions)]
                RAD = ff.get_hottest_SAD(unique_RADs)

            elif model == 'IntPartsN':
                # Feasible set of partitions, based on N
                partitions = parts.rand_partitions(N, S, 100)
                unique_RADs = [list(x) for x in set(tuple(x) for x in partitions)]
                RAD = ff.get_hottest_SAD(unique_RADs)

            elif model == 'Zipf':
                zipf_pred = dist.zipf(RAD)
                RAD = zipf_pred.from_cdf()

            elif model == 'PLN': # Poisson Log-normal
                RAD = pln.get_rad_from_obs(RAD, 'pln')

            print model, kind, name, numRADs - ct

            Evar = 0.0 #metrics.e_var(RAD)
            ESimp = metrics.simpsons_evenness(RAD)
            ENee = 0.0 #metrics.NHC_evenness(RAD)
            EPielou = 0.0 #metrics.pielous_evenness(RAD)

            EHeip = 0.0 # metrics.Heips_evenness(RAD)
            EQ = 0.0 # metrics.EQ_evenness(RAD)

            BP = metrics.Berger_Parker(RAD)
            SimpDom = 0.0 # metrics.simpsons_dom(RAD)
            perc_ones = 0.0 # metrics.Singletons(RAD)

            Nmax = max(RAD)

            rareRel = 0.0 # metrics.rarityRel(RAD)
            rareOnes = 0.0 # metrics.rarityOnes(RAD)

            skew = stats.skew(RAD)

            ct+=1

            print>>OUT, name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew

        OUT.close()

########### LOOP THROUGH DATA ##################################################
