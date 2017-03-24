from __future__ import division
from numpy import mean, log10
from os.path import expanduser
from scipy import stats
import sys
import numpy as np

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")
GenPath = mydir + "GitHub/simplex/results/simulated_data/"

from diversity_metrics import *
from spatial_functions import *


def output(iD, sD, rD, ps, sim, N, R, ct, prod, splist2):

    h, l, r, u = ps

    m, nN, rmax, gmax, maintmax, dmax, amp, freq, phase, pmax, dormlim = 1, 3, 1, 1, 1, 1, 0, 0, 0, 1, 1

    SpIDs, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, DormList, ADList, SizeList, indX, indY = [list([]) for _ in xrange(14)]
    RIDs, Rvals, Rtypes = [list([]) for _ in xrange(3)]

    N, S, R = 0, 0, 0

    for k, v in rD.items():
            RIDs.append(k)
            Rvals.append(v['v'])
            Rtypes.append(v['t'])

    for k, v in iD.items():
            IndIDs.append(k)
            SpIDs.append(v['sp'])

            GrowthList.append(v['gr'])
            MaintList.append(v['mt'])
            MFDList.append(v['mf'])
            RPFList.append(v['rp'])
            N_RList.append(v['ef'])
            DispList.append(v['di'])
            ADList.append(v['st'])
            SizeList.append(v['sz'])
            Qs.append(v['q'])
            indX.append(v['x'])
            indY.append(v['y'])

    N = len(IndIDs)
    S = len(list(set(SpIDs)))
    if N > 0:

        numD = ADList.count('d')
        numA = N - numD

        try: pD = 100*(numD/N)
        except: pD = float('NaN')

        RAD, splist = metrics.GetRAD(SpIDs)
        S, R, N = len(RAD), len(RIDs), len(SpIDs)

        Lists = [SpIDs, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList, SizeList]

        aLists, dLists = metrics.separateCom(Lists)
        a_SpeciesIDs, a_IndIDs, a_Qs, a_GrowthList, a_MaintList, a_MFDList, a_RPFList, a_N_RList, a_DispList, a_SizeList = aLists
        d_SpeciesIDs, d_IndIDs, d_Qs, d_GrowthList, d_MaintList, d_MFDList, d_RPFList, d_N_RList, d_DispList, d_SizeList = dLists

        aRAD, asplist = metrics.GetRAD(a_SpeciesIDs)
        dRAD, dsplist = metrics.GetRAD(d_SpeciesIDs)
        aCOM = []
        dCOM = []

        for sp in splist:

            if sp in asplist:
                i = asplist.index(sp)
                aCOM.append(aRAD[i])
            else: aCOM.append(0)

            if sp in dsplist:
                i = dsplist.index(sp)
                dCOM.append(dRAD[i])
            else: dCOM.append(0)


        ES = metrics.e_simpson(RAD)
        EV = metrics.e_var(RAD)
        Nm = max(RAD)

        skew = stats.skew(RAD)
        lms = log10(abs(float(skew)) + 1)
        if skew < 0:
            lms = lms * -1

        wt = 0
        if N == 1: splist2 = list(splist)
        else:
            wt = metrics.WhittakersTurnover(splist, splist2)
            splist2 = list(splist)

        SA, SD = len(aRAD), len(dRAD)

        all_G = mean(GrowthList)
        all_M = mean(MaintList)
        all_avgMF = mean(MFDList)
        all_avgRPF = mean(RPFList)
        all_Disp = mean(DispList)
        all_Size = mean(SizeList)

        a_G = mean(a_GrowthList)
        a_M = mean(a_MaintList)
        a_avgMF = mean(a_MFDList)
        a_avgRPF = mean(a_RPFList)
        a_Disp = mean(a_DispList)
        a_Size = mean(a_SizeList)

        d_G = mean(d_GrowthList)
        d_M = mean(d_MaintList)
        d_avgMF = mean(d_MFDList)
        d_avgRPF = mean(d_RPFList)
        d_Disp = mean(d_DispList)
        d_Size = mean(d_SizeList)

        List = []
        for nr in N_RList:
            List.append(np.var(nr))

        all_NR = mean(List)
        if sum(List) == 0: all_NR = 0

        a_NR = mean(a_N_RList)
        if sum(a_N_RList) == 0: a_NR = 0

        d_NR = mean(d_N_RList)
        if sum(d_N_RList) == 0: d_NR = 0

        a_Q = np.mean(a_Qs)
        d_Q = np.mean(d_Qs)
        all_Q = np.mean(Qs)

        avgN = 0
        if S > 0: avgN = N/S

        Nvar = np.var(RAD)
        Rdens = R/(h*l)

        OUT = open(GenPath + 'SimData.csv', 'a')
        outlist = [sim, ct, m, r, nN, rmax, gmax, maintmax, dmax, 500, u, h, l, N, prod, prod, R, Rdens,\
        S, ES, EV, avgN, Nvar, Nm, lms, wt, amp, freq, phase, all_Q, a_Q, d_Q, all_G, a_G, d_G, all_M, a_M, d_M,\
        all_NR, a_NR, d_NR, all_Disp, a_Disp, d_Disp, all_avgRPF, a_avgRPF, d_avgRPF, all_avgMF, a_avgMF, d_avgMF,\
        all_Size, a_Size, d_Size, numA, SA, numD, SD, pD, dormlim]

        outlist = str(outlist).strip('[]')
        outlist = outlist.replace(" ", "")
        print>>OUT, outlist
        OUT.close()

        rad = str(RAD).strip('[]')
        rad = rad.replace(" ", "")
        OUT = open(GenPath + 'RAD-Data.csv', 'a')
        print>>OUT, sim, ',', ct,',',  rad
        OUT.close()

        SAR = spatial.SAR(indX, indY, SpIDs, h)
        sar = str(SAR).strip('[]')
        sar = sar.replace(" ", "")
        OUT = open(GenPath + 'SAR-Data.csv', 'a')
        print>>OUT, sim, ',', ct,',',  sar
        OUT.close()

        print 'sim:', '%3s' % sim, 'ct:', '%3s' % ct,'  N:', '%4s' %  N,
        print '  S:', '%4s' %  S, '  R:', '%4s' % R, ' u0:', '%4s' % round(u, 4)

    return splist2
