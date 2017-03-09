from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats
import time

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

        print name
        xnames = ['AGSOIL', '.DS_Store', 'CATLIN', 'CHU', 'EMPclosed']
        if name in xnames: continue

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

        ct = 0
        numRADs = len(RADs)
        for RAD in RADs:

            if kind == 'micro':
                RAD = list([x for x in RAD if x > 0])

            elif kind == 'macro':
                RAD = list([x for x in RAD if x > 0])


            N = sum(RAD)
            S = len(RAD)

            if S < 10 or N < 11: continue

            times = []
            # Evenness
            t = time.clock()
            Evar = metrics.e_var(RAD)
            t = time.clock() - t
            times.append(['Evar', t])

            t = time.clock()
            ESimp = metrics.e_simpson(RAD)
            t = time.clock() - t
            times.append(['Simpsons Evenness', t])

            t = time.clock()
            EQ = metrics.EQ(RAD)
            t = time.clock() - t
            times.append(['EQ', t])

            t = time.clock()
            O = metrics.OE(RAD)
            t = time.clock() - t
            times.append(['O', t])

            t = time.clock()
            Camargo = 'nan' #metrics.camargo(RAD)
            t = time.clock() - t
            times.append(['Camargo', t])
            #ENee = metrics.NHC(RAD)
            #EPielou = metrics.e_pielou(RAD)
            #EHeip = metrics.e_heip(RAD)


            # Dominance
            t = time.clock()
            BP = metrics.Berger_Parker(RAD)
            t = time.clock() - t
            times.append(['Berger Parker', t])

            t = time.clock()
            SimpDom = metrics.simpsons_dom(RAD)
            t = time.clock() - t
            times.append(['Simpsons Dominance', t])

            t = time.clock()
            Nmax = max(RAD)
            t = time.clock() - t
            times.append(['Nmax', t])

            t = time.clock()
            McN = metrics.McNaughton(RAD)
            t = time.clock() - t
            times.append(['McNaughton', t])

            # Rarity
            t = time.clock()
            skew = stats.skew(RAD)
            t = time.clock() - t
            times.append(['Skew', t])

            t = time.clock()
            logskew = stats.skew(np.log(RAD))
            t = time.clock() - t
            times.append(['log-Skew', t])

            t = time.clock()
            p_ones = metrics.r_singletons(RAD)
            t = time.clock() - t
            times.append(['Percent Singletons', t])

            # Preston's alpha:
            t = time.clock()
            preston = metrics.Prestons_a(RAD)
            t = time.clock() - t
            times.append(['Preston\'s alpha', t])

            max_time = ['none', 0.0]
            for i, t in enumerate(times):
                if t[1] > max_time[1]:
                    max_time = t
            ct+=1

            print name, N, S, times,'\n'

            #print name, kind, N, S, Evar, ESimp, EQ, O, Camargo, BP, SimpDom, Nmax, McN, skew, logskew, p_ones, preston[0]
            #print name, numRADs - ct

    return


getMetrics()
