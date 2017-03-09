from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats


mydir = os.path.expanduser("~/")

sys.path.append(mydir + "tools/metrics")
import metrics
sys.path.append(mydir + "tools/getRADs")
import getRADs

def getMetrics():
    datasets = []
    for name in os.listdir(mydir +'data/micro'):
            datasets.append([name, 'micro'])
    for name in os.listdir(mydir +'data/macro'):
            datasets.append([name, 'macro'])


    for dataset in datasets:

        name = dataset[0] # name of dataset
        kind = dataset[1] # micro or macro

        #print name, kind
        #xnames = ['AGSOIL', '.DS_Store', 'SLUDGE']
        #xnames = ['HMP']

        if name == '.DS_Store': continue
        if name in ['MGRAST'] : pass
        else: continue

        #OUT = open(mydir + 'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt','w+')
        OUT = open(mydir + 'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt','w+')
        RADs = []

        if kind == 'macro':
            RADs = getRADs.get_SADs(mydir +'data/'+kind+'/'+name+'/', name)
            print 'macro', name, len(RADs)


        if kind == 'micro':
            #RADs = getPhyloSADs.get_SADs(seqID, dataset[0], 'genus')

            if name == 'EMPclosed' or name == 'EMPopen':
                RADs = getRADs.EMP_SADs(mydir +'data/'+kind+'/'+name+'/', name)

            else:
                RADs = getRADs.get_SADs(mydir +'data/'+kind+'/'+name+'/', name)

            print 'micro', name, len(RADs)

        ct = 0
        numRADs = len(RADs)
        for RAD in RADs:

            if kind == 'micro':
                RAD = list([x for x in RAD if x > 0]) # greater than 1 means singletons are excluded

            elif kind == 'macro':
                RAD = list([x for x in RAD if x > 0])


            N = sum(RAD)
            S = len(RAD)

            if S < 10: continue
            if max(RAD) == min(RAD): continue

            # Evenness
            Var = np.var(RAD, ddof = 1)
            Evar = metrics.e_var(RAD)
            ESimp = metrics.e_simpson(RAD)
            EQ = metrics.EQ(RAD)
            O = metrics.OE(RAD)
            #Camargo = 0.0 # metrics.camargo(RAD)   # Takes too long
            ENee = metrics.NHC(RAD)
            EPielou = metrics.e_pielou(RAD)
            EHeip = metrics.e_heip(RAD)


            # Dominance
            BP = metrics.Berger_Parker(RAD)
            SimpDom = metrics.simpsons_dom(RAD)
            Nmax = max(RAD)
            McN = metrics.McNaughton(RAD)

            # Rarity
            skew = stats.skew(RAD)
            logskew = metrics.Rlogskew(RAD)
            #p_ones = metrics.r_singletons(RAD)
            #p_zpt1 = metrics.p_ZPtOne(RAD)

            # Preston's alpha and richness, from Curtis and Sloan (2002).
            # Estimating prokaryotic diversity and its limits. PNAS.
            preston_a, preston_S = metrics.Preston(RAD)

            # Richness estimators
            chao1, ace, jknife1, jknife2 = metrics.EstimateS1(RAD)
            margalef = metrics.Margalef(RAD)
            menhinick = metrics.Menhinick(RAD)

            ct+=1

            print>>OUT, name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S
            #print name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S

            print name, numRADs - ct
        OUT.close()

    return


getMetrics()
