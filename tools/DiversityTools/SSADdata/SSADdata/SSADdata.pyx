from __future__ import division
import sys
import os
#import random
import numpy as np
from scipy import stats


mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics

#OUT1 = open("/GitHub/rare-bio/output/EMPclosed-SSAD-ResultsTable.txt",'w+')


def get_EMP_SSADs():

    DATA = mydir + "data/micro/EMPopen/EMPopen-SSADdata.txt"

    SSADdict = {}

    with open(DATA) as f:

        for d in f:
            if d.strip():

                d = d.split()
                species = d[0]
                #sample = d[1]
                abundance = float(d[2])

                if abundance > 0:
                    if species in SSADdict:
                        SSADdict[species].append(abundance)
                    else:
                        SSADdict[species] = [abundance]



    SSADs = []
    SSADlist = SSADdict.items()

    for tup in SSADlist:

        SSAD = tup[1]
        if len(SSAD) >= 1:
            SSAD.sort()
            SSAD.reverse()
            SSADs.append(SSAD)

    return SSADs



def get_Macro_SSADs():

    for name in os.listdir('/data/macro'):

        if name == 'BCI': continue
        if name == '.DS_Store': continue

        else:
            print name

            SSADdict = {}

            path ='/data/macro/'+name+'/'+name
            OUT = open(path + '-SSADMetricData.txt','w+')
            DATA = path+'-data.txt'

            with open(DATA) as f:

                for d in f:
                    if d.strip():
                        d = d.split()

                        if name == 'GENTRY':
                            species = name + d[1] # Dataset name plus species identifier
                            abundance = float(d[2])

                        else:
                            #site = d[0]
                            #year = d[1]
                            species = name + d[2] # Dataset name plus species identifier
                            abundance = float(d[3])

                        if abundance > 0:
                            if species in SSADdict:
                                SSADdict[species].append(abundance)
                            else:
                                SSADdict[species] = [abundance]



            SSADs = []
            SSADlist = SSADdict.items()

            for tup in SSADlist:

                SSAD = tup[1]
                if len(SSAD) >= 1:
                    SSAD.sort()
                    SSAD.reverse()
                    SSADs.append(SSAD)


            num = len(SSADs)
            ct = 0

            for SSAD in SSADs:

                SSAD = list([x for x in SSAD if x != 0]) # removes zeros

                N = int(sum(SSAD))
                S = int(len(SSAD))

                if S < 1: continue

                print name, 'N:',N,'S:',S,'  ', num - ct

                Evar = metrics.e_var(SSAD)
                ESimp = metrics.simpsons_evenness(SSAD)
                ENee = metrics.NHC_evenness(SSAD)
                EPielou = metrics.pielous_evenness(SSAD)

                EHeip = metrics.Heips_evenness(SSAD)
                EQ = metrics.EQ_evenness(SSAD)

                BP = metrics.Berger_Parker(SSAD)
                SimpDom = 0 #metrics.simpsons_dom(SSAD)

                rareRel = metrics.rarityRel(SSAD)
                rareOnes = metrics.rarityOnes(SSAD)

                skew = stats.skew(SSAD)

                ct += 1

                print>>OUT, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew

            print 'number of SSADs:', len(SSADs)

            OUT.close()

    return


#get_Macro_SSADs()
SSADs = get_EMP_SSADs()

maxN = 0
for ssad in SSADs:
    if max(ssad) > maxN:
        maxN = max(ssad)

print 'abundance of the most abundant OTU in the EMP open reference data set is:', maxN
