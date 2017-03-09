from __future__ import division
from random import shuffle, seed #, randint
import numpy as np
from numpy import sin, pi, mean
import sys
import os
import time
from scipy import stats

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/residence-time/model/metrics")
import metrics
sys.path.append(mydir + "GitHub/residence-time/model/bide")
import bide
sys.path.append(mydir + "GitHub/residence-time/model/randparams")
import randparams as rp
sys.path.append(mydir + "GitHub/residence-time/model/spatial")
import spatial

GenPath = mydir + 'GitHub/residence-time/results/simulated_data/'


OUT = open(GenPath + 'SimData.csv','w+')
print>>OUT, 'sim,ct,immigration.rate,res.inflow,N.types,max.res.val,max.growth.rate,max.met.maint,max.active.dispersal,starting.seed,\
flow.rate,height,length,width,total.abundance,ind.production,biomass.prod.N,resource.particles,resource.concentration,\
species.richness,simpson.e,e.var,avg.pop.size,pop.var,N.max,logmod.skew,Whittakers.turnover,amplitude,frequency,phase,\
all.biomass,active.biomass,dormant.biomass,all.avg.per.capita.growth,active.avg.per.capita.growth,dormant.avg.per.capita.growth,\
all.avg.per.capita.maint,active.avg.per.capita.maint,dormant.avg.per.capita.maint,all.avg.per.capita.efficiency,active.avg.per.capita.efficiency,dormant.avg.per.capita.efficiency,\
all.avg.per.capita.active.dispersal,active.avg.per.capita.active.dispersal,dormant.avg.per.capita.active.dispersal,all.avg.per.capita.rpf,active.avg.per.capita.rpf,dormant.avg.per.capita.rpf,\
all.avg.per.capita.mf,active.avg.per.capita.mf,dormant.avg.per.capita.mf,all.size,active.size,dormant.size,N.active,S.active,N.dormant,S.Dormant,Percent.Dormant,dorm.limit'
OUT.close()

OUT = open(GenPath + 'RAD-Data.csv', 'w+')
OUT.close()

OUT = open(GenPath + 'SAR-Data.csv', 'w+')
OUT.close()


#######################  COMMUNITY PARAMETERS  #########################

Ns, RAD, splist, splist2, RTypes, RX, RY, RZ, RIDs, RVals = [list([]) for _ in xrange(10)]
nN, u1, numA, numD, RDENS, RDiv, RRich, S, ES, EV, BP, SD, Nm, sk, Mu, Maint, ID, RID, N, ct, T, R, PRODI, PRODN = [0]*24
SpDict, IndDict = {}, {}

################ MODEL INPUTS ##################################
params = rp.get_rand_params()
width, height, length, seedCom, m, r, nN, rmax, gmax, maintmax, dmax, amp, freq, phase, rates, pmax, dormlim, smax = params

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
u0 = rates[0]  # initial in-flow speed

processes = range(1, 12)
t = time.clock()
BurnIn = 'not done'
p, sim, ctr2 = 0.0, 1, 1


while sim < 100000:

    seed()
    ct += 1
    numc = 0
    plot_system = 'no'

    # fluctuate flow according to amplitude, frequency, & phase
    u1 = float(u0) #+ u0*(amp * sin(2*pi * ct * freq + phase))
    if u1 > 1.0: u1 = u0

    shuffle(processes)
    for num in processes:

        if num == 1: # Inflow of resources
            RTypes, RVals, RX, RY, RZ, RIDs, RID = bide.ResIn(RTypes, RVals, RX, RY, RZ, RID, RIDs, params, u1)

        elif num == 2: # Resource flow
            RTypes, RVals, RX, RY, RZ, RIDs, RID = bide.res_flow(RTypes, RVals, RX, RY, RZ, RID, RIDs, params, u1)

        elif num == 3: # Inflow of individuals (immigration)
            SpDict, IndDict, ID = bide.immigration(SpDict, IndDict, params, ID, ct, u1)

        elif num == 4: # flowthrough of individuals
            IndDict = bide.ind_flow(SpDict, IndDict, height, length, width, u1)

        elif num == 5: # Active dispersal
            IndDict = bide.ind_disp(SpDict, IndDict, height, length, width, u1)

        elif num == 11: # Search
            IndDict = bide.search(SpDict, IndDict, height, length, width, u0, RTypes, RVals, RX, RY, RZ, RIDs)

        elif num == 6: # Consume
            numc, SpDict, IndDict, height, length, width, u0, Rtypes, Rvals, RX, RY, RZ, RIDs = bide.consume(SpDict, IndDict, height, length, width, u0, RTypes, RVals, RX, RY, RZ, RIDs)

        elif num == 10: # Grow
            IndDict = bide.grow(SpDict, IndDict)

        elif num == 7: # Transition
            IndDict = bide.transition(SpDict, IndDict)

        elif num == 8: # Maintenance
            IndDict = bide.maintenance(SpDict, IndDict)

        elif num == 9: # Reproduction

            Qs = []
            for key, value in IndDict.items(): Qs.append(value['q'])
            p1, TNQ1 = metrics.getprod(Qs)

            SpDict, IndDict, ID = bide.reproduce(u0, SpDict, IndDict, ID)

            Qs = []
            for key, value in IndDict.items(): Qs.append(value['q'])
            p2, TNQ2 = metrics.getprod(Qs)

            PRODI = p2 - p1
            PRODN = TNQ2 - TNQ1


    R, N = len(RIDs), len(list(IndDict))

    Ns.append(N)
    RDENS = R/(width**3)

    tau = np.log10((width**1)/u1)
    minct = 100 + int(round(4**tau))

    #print 'sim:', '%4s' % sim, 'ct:', '%3s' % ctr2, '  t:', '%6s' % str(round(minct - ct)), '  tau:', '%5s' %  round(tau,3), '  width:', '%4s' %  round(width,1), 'flow:', '%5s' %  round(u1,4), '   N:', '%4s' %  N, '  R:', '%3s' % R, '  C:', '%4s' % numc


    if BurnIn == 'not done':
        if N == 0 or ct >= minct:
            BurnIn = 'done'
            Ns = [Ns[-1]] # only keep the most recent N value

    if BurnIn == 'done':

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
        if N > 0:
            percD = 100*(numD/N)


        RAD, splist = [], []
        if len(SpIDs) >= 1:
            RAD, splist = metrics.GetRAD(SpIDs)

        RTAU, INDTAU, TTAU = 0, 0, 0

        # Number of tracers, resource particles, and individuals
        S, R, N = len(RAD), len(RIDs), len(SpIDs)

        if N >= 1 and ct%10 == 0:

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
            lms = np.log10(np.abs(float(skew)) + 1)
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
            u1, height, length, width, N, PRODI, PRODN, R, RDENS,\
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

            SAR = spatial.SAR(indX, indY, indZ, SpIDs, width)

            sar = str(SAR).strip('[]')
            sar = sar.replace(" ", "")
            OUT = open(GenPath + 'SAR-Data.csv', 'a')
            print>>OUT, sim, ',', ctr2,',',  sar
            OUT.close()


        if len(Ns) > 100:
            ctr2 += 1
            print 'ct:', '%4s' % minct, 'sim:', '%3s' % sim, 'tau:', '%5s' %  round(tau,2), 'V:', '%4s' %  width**3,'  F:', '%6s' %  round(u1,3), '  N:', '%4s' %  N, ' S:', '%4s' % S, ' Simp:', '%4s' % round(ES,2), ' Evar:', '%4s' % round(EV,2), 'R:', '%4s' % R, '%D:', '%4s' % round(percD,2), ' C:', '%4s' % numc

            rates = np.roll(rates, -1, axis=0)
            u0 = rates[0]

            Ns, RAD, splist, splist2, RTypes, RX, RY, RZ, RIDs, RVals = [list([]) for _ in xrange(10)]
            nN, u1, numA, numD, RDENS, RDiv, RRich, S, ES, EV, BP, SD, Nm, sk, Mu, Maint, ID, RID, N, ct, T, R, PRODI, PRODN = [0]*24
            SpDict, IndDict = {}, {}

            p, t, BurnIn = 0, 0, 'not done'

            if u0 == max(rates):
                print '\n'
                ################ MODEL INPUTS ##################################
                params = rp.get_rand_params(width)
                width, height, length, seedCom, m, r, nN, rmax, gmax, maintmax, dmax, amp, freq, phase, rates, pmax, dormlim, smax = params
                SpDict, IndDict = {}, {}
                sim += 1
