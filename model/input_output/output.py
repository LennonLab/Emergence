from __future__ import division
import sys
from numpy import mean, log10
from os.path import expanduser
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")

from diversity_metrics import *
from spatial_functions import *



def output(IndDict, SpDict, ResDict, sim, ctr2):
    tau = np.log10((w**1)/u1)

    SpIDs, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, DormList, ADList, SizeList, indX, indY, indZ = [list([]) for _ in xrange(15)]
    for key, value in IndDict.items():
            IndIDs.append(key)
            sp = value['spID']
            SpIDs.append(sp)
    
            GrowthList.append(SpDict[sp]['grow'])
            MaintList.append(SpDict[sp]['maint'])
            MFDList.append(SpDict[sp]['mfd'])
            RPFList.append(SpDict[sp]['rpf'])
            N_RList.append(SpDict[sp]['eff'])
            DispList.append(SpDict[sp]['disp'])
            DormList.append(SpDict[sp]['dlim'])
            ADList.append(value['state'])
            SizeList.append(value['size'])
            Qs.append(value['q'])
            indX.append(value['x'])
            indY.append(value['y'])
            indZ.append(value['z'])
    
    numD = ADList.count('d')
    numA = N - numD
    
    percD = 0
    if N > 0: percD = 100*(numD/N)
    
    RAD, splist = [], []
    if len(SpIDs) >= 1:
        RAD, splist = metrics.GetRAD(SpIDs)
    
    S, R, N = len(RAD), len(RIDs), len(SpIDs)
    
    if N > 0:
        Lists = [SpIDs, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, DormList, ADList, SizeList]
    
        aLists, dLists = metrics.separateCom(Lists)
        a_SpeciesIDs, a_IndIDs, a_Qs, a_GrowthList, a_MaintList, a_MFDList, a_RPFList, a_N_RList, a_DispList, a_DormList, a_SizeList = aLists
        d_SpeciesIDs, d_IndIDs, d_Qs, d_GrowthList, d_MaintList, d_MFDList, d_RPFList, d_N_RList, d_DispList, d_DormList, d_SizeList = dLists
    
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
            # log-modulo transformation of skewnness
            lms = log10(abs(float(skew)) + 1)
            if skew < 0:
                lms = lms * -1
    
                wt = 0
                if len(Ns) == 1:
                    splist2 = list(splist)
                if len(Ns) > 1:
                    wt = metrics.WhittakersTurnover(splist, splist2)
                    splist2 = list(splist)
    
                SA = len(aRAD)
                SD = len(dRAD)
    
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
                if sum(List) == 0:
                    all_NR = 0
    
                a_NR = mean(a_N_RList)
                if sum(a_N_RList) == 0:
                    a_NR = 0
    
                d_NR = mean(d_N_RList)
                if sum(d_N_RList) == 0:
                    d_NR = 0
    
                a_Q = np.mean(a_Qs)
                d_Q = np.mean(d_Qs)
                all_Q = np.mean(Qs)
    
                avgN = 0
                if S > 0:
                    avgN = N/S
    
                Nvar = np.var(RAD)
    
                Nmeans = [np.var(x) for x in zip(*N_RList)]
                NR = mean(Nmeans)
                OUT = open(GenPath + 'SimData.csv', 'a')
                outlist = [sim, ctr2, m, r, nN, rmax, gmax, maintmax, dmax, seedCom,\
                u1, h, l, w, N, PRODI, PRODN, R, RDENS,\
                S, ES, EV, avgN, Nvar, Nm, lms, wt, amp, freq, phase,\
                all_Q, a_Q, d_Q, all_G, a_G, d_G, all_M, a_M, d_M, all_NR, a_NR, d_NR,\
                all_Disp, a_Disp, d_Disp, all_avgRPF, a_avgRPF, d_avgRPF,\
                all_avgMF, a_avgMF, d_avgMF, all_Size, a_Size, d_Size, numA, SA, numD, SD, percD, dormlim]
    
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")
                print>>OUT, outlist
                OUT.close()
    
                rad = str(RAD).strip('[]')
                rad = rad.replace(" ", "")
                OUT = open(GenPath + 'RAD-Data.csv', 'a')
                print>>OUT, sim, ',', ctr2,',',  rad
                OUT.close()
    
                SAR = spatial.SAR(indX, indY, indZ, SpIDs, w)
    
                sar = str(SAR).strip('[]')
                sar = sar.replace(" ", "")
                OUT = open(GenPath + 'SAR-Data.csv', 'a')
                print>>OUT, sim, ',', ctr2,',',  sar
                OUT.close()
                
    return
    