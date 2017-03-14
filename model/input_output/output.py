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


def output(IndDict, SpDict, ResDict, params, sim, ls, Ns, splist2):

    ''' A function to quantify aspects of abundance and diversity among
        species, individuals, and resource particles.


        SpDict  :  A python dictionary holding properties of specific species:

                   's_min' : lower bound on individual size
                   's_max' : upper bound on individual size

                   'grow' :  maximum growth rate
                   'disp' :  maximum dispersal rate

                   'rpf' :  probability of random transition to active state
                   'maint' : basal metabolic maintenance cost

                   'dlim' : lower size limit at which individuals go dormant
                   'spec' : speciation rate

                   'mfd' : factor by which dormant decreases costs of
                   maintenance energy


                   'eff' : Resource particles can belong to one of n types for
                   which species have varying abilities to grow on.


        IndDict  : A python dictionary to hold properties of individual
                   organisms

                   'size' : The size of the individual is an abstract
                   quantity, but can vary over two orders of magnitude

                   'q' : level of endogenous resources
                   'spID' : species ID
                   'state' : whether active or dormant

                  'x' : x-coordinate
                  'y' : y-coordinate
                  'z' : z-coordinate


        params  :  General model parameters. Not all are used in every function.

                   w : width of the system
                   h : height of the system
                   l : length of the system

                   seed : Number of starting individuals
                   m : immigration rate and the probability of an individual
                   organism immigrating per time step

                   r : Maximum number of resource particles entering per time step
                   rmax : Maximum size of individual resource particles

                   nN : Number of inflowing resource types
                   gmax : Maximum specific growth rate
                   maintmax : Maximum metabolic maintenance
                   dmax : Maximum dispersal rate
                   pmax : maximum probability of random resuscitation
                   dormlim : level of endogenous resource at which
                   individual go dormant
                   smax : Maximum size of any individual

                   amp : amplitude of environmental flux
                   freq : frequency of environmental flux
                   phase : phase of individual immigration and resource inflow
                   rate : rate of system flow through

        ls  :  A list holding the following pieces of information

               ct : time step number
               numc : number of consumption events
               prodI : individual productivity
               prodN : biomass productivity

        sim  :  Simulation number

    '''

    w, h, l, seed, m, r, nN, rmax, gmax, maintmax, dmax, amp, freq, phase, rate, pmax, dormlim, smax, st = params
    ct, numc, prodI, prodN = ls

    SpIDs, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, DormList, ADList, SizeList, indX, indY, indZ = [list([]) for _ in xrange(15)]
    RIDs, Rvals, Rtypes = [list([]) for _ in xrange(3)]

    N, S, R = 0, 0, 0

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

    N = len(IndIDs)
    if N > 0:

        for key, value in ResDict.items():
            RIDs.append(key)
            Rvals.append(value['size'])
            Rtypes.append(value['type'])

        numD = ADList.count('d')
        numA = N - numD

        try: pD = 100*(numD/N)
        except: pD = 0

        RAD, splist = metrics.GetRAD(SpIDs)
        S, R, N = len(RAD), len(RIDs), len(SpIDs)

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
        else:
            wt = metrics.WhittakersTurnover(splist, splist2)
            splist2 = list(splist)

        SA = len(aRAD)
        SD = len(dRAD)

        all_G = mean(GrowthList)
        all_M = mean(MaintList)
        all_avgMF = max(MFDList)
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
        Rdens = R/(w*h*l)

        OUT = open(GenPath + 'SimData.csv', 'a')
        outlist = [sim, ct, m, r, nN, rmax, gmax, maintmax, dmax, seed,\
        rate, h, l, w, N, prodI, prodN, R, Rdens,\
        S, ES, EV, avgN, Nvar, Nm, lms, wt, amp, freq, phase,\
        all_Q, a_Q, d_Q, all_G, a_G, d_G, all_M, a_M, d_M, all_NR, a_NR, d_NR,\
        all_Disp, a_Disp, d_Disp, all_avgRPF, a_avgRPF, d_avgRPF,\
        all_avgMF, a_avgMF, d_avgMF, all_Size, a_Size, d_Size, numA, SA, numD, SD, pD, dormlim]

        outlist = str(outlist).strip('[]')
        outlist = outlist.replace(" ", "")
        print>>OUT, outlist
        OUT.close()

        rad = str(RAD).strip('[]')
        rad = rad.replace(" ", "")
        OUT = open(GenPath + 'RAD-Data.csv', 'a')
        print>>OUT, sim, ',', ct,',',  rad
        OUT.close()

        SAR = spatial.SAR(indX, indY, indZ, SpIDs, w)

        sar = str(SAR).strip('[]')
        sar = sar.replace(" ", "")
        OUT = open(GenPath + 'SAR-Data.csv', 'a')
        print>>OUT, sim, ',', ct,',',  sar
        OUT.close()

        Ns.append(N)

    return [Ns, N, S, R, splist2]
