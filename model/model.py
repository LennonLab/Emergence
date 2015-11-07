from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt

import statsmodels.tsa.stattools as sta
from math import isnan
from random import choice
from scipy import stats
import numpy as np
from numpy import sin, pi, mean
import sys
import os
import time
import psutil

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/tools/metrics")
import metrics
sys.path.append(mydir + "GitHub/simplex/tools/LBM")
import LBM
sys.path.append(mydir + "GitHub/simplex/tools/bide")
import bide
sys.path.append(mydir + "GitHub/simplex/tools/randparams")
import randparams as rp

""" To generate movies:
1.) uncomment line 'ani.save' on or near line 364
2.) adjust the frames variable on or near line 364, to change movie length
3.) change plot_system = 'no' to 'yes' on or near line 66

Because generating animations requires computing time and memory, doing so can
be computationally demanding. To quicken the process, use plot_system = 'no' on
or near line 66.
"""

# https://www.quantstart.com/articles/Basics-of-Statistical-Mean-Reversion-Testing
# http://statsmodels.sourceforge.net/0.5.0/generated/statsmodels.tsa.stattools.adfuller.html


GenPath = mydir + 'GitHub/simplex/results/simulated_data/'

"""
############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(GenPath + 'examples/SimData.csv','w')
OUT2 = open(GenPath + 'examples/RADs.csv','w')
OUT3 = open(GenPath + 'examples/Species.csv','w')
OUT4 = open(GenPath + 'examples/IndRTD.csv','w')
OUT5 = open(GenPath + 'examples/TracerRTD.csv','w')
OUT6 = open(GenPath + 'examples/ResRTD.csv','w')

# printing physical variables, residence times, community diversity properties
# physiological values, trait values, resource values
print>>OUT1, 'RowID, motion, ind.production, biomass.prod.N, biomass.prod.P, biomass.prod.C, res.inflow, N.types, P.types, C.types, max.res.val, max.growth.rate, max.met.maint, max.active.dispersal, barriers, logseries.a, starting.seed, flow.rate, width, height, viscosity, total.abundance, immigration.rate, resource.tau, particle.tau, individual.tau, resource.concentration, shannons.resource.diversity, resource.richness, species.richness, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, tracer.particles, resource.particles, speciation, Whittakers.turnover, Jaccards.dissimilarity, Sorensens.dissimilarity, avg.per.capita.growth, avg.per.capita.maint, avg.per.capita.N.efficiency, avg.per.capita.P.efficiency, avg.per.capita.C.efficiency, avg.per.capita.active.dispersal, amplitude, flux, frequency, phase, disturbance'

OUT1.close()
OUT2.close()
OUT3.close()
OUT4.close()
OUT5.close()
OUT6.close()
"""

def nextFrame(arg):

    """ Function called for each successive animation frame; arg is the frame number """

    global p, BurnIn, t, num_sims, width, height, Rates, u0, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, SpColorDict, GrowthDict, N_RD, P_RD, C_RD, DispDict, MaintDict, one9th, four9ths, one36th, barrier, gmax, dmax, maintmax, IndIDs, Qs, IndID, IndTimeIn, IndExitAge, IndX, IndY,  Ind_scatImage, SpeciesIDs, EnvD, TY, tracer_scatImage, TTimeIn, TIDs, TExitAge, TX, RTypes, RX, RY, RID, RIDs, RVals, RTimeIn, RExitAge, resource_scatImage, bN, bS, bE, bW, bNE, bNW, bSE, bSW, ct1, Mu, Maint, motion, reproduction, speciation, seedCom, m, r, nNi, nP, nC, rmax, sim, RAD, splist, N, ct, splist2, WTs, Jcs, Sos, RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, T, R, LowerLimit, prod_i, prod_q, viscosity, alpha, Ts, Rs, PRODIs, Ns, TTAUs, INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs, PRODNs, PRODPs, PRODCs, lefts, bottoms, Gs, Ms, NRs, PRs, CRs, Ds, RTAUs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, amp, freq, flux, pulse, phase, disturb, envgrads, barriers

    ct += 1
    plot_system = 'yes'
    # fluctuate flow according to amplitude, frequency, & phase
    u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))
    if u1 > 1: u1 == 1.0

    # Fluid dynamics
    nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier])
    rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

    # Inflow of tracers
    TIDs, TTimeIn, TX, TY = bide.NewTracers(motion,TIDs, TX, TY, TTimeIn, width, height, u0)
    # moving tracer particles
    if len(TIDs) > 0:
        if motion == 'fluid':
            TIDs, TX, TY, TExitAge, TTimeIn = bide.fluid_movement('tracer', TIDs, TTimeIn, TExitAge, TX, TY, ux, uy, width, height, u0)
        else:
            TIDs, TX, TY, TExitAge, TTimeIn = bide.nonfluid_movement('tracer', motion, TIDs, TTimeIn, TExitAge, TX, TY, ux, uy, width, height, u0)

    # Inflow of resources
    RTypes, RVals, RX, RY,  RIDs, RID, RTimeIn = bide.ResIn(motion, RTypes, RVals, RX, RY,  RID, RIDs, RTimeIn, r, rmax, nNi, nP, nC, width, height, u1)

    # resource flow
    Lists = [RTypes, RIDs, RID, RVals]
    if len(RTypes) > 0:
        if motion == 'fluid':
            RTypes, RX, RY,  RExitAge, RIDs, RID, RTimeIn, RVals = bide.fluid_movement('resource', Lists, RTimeIn, RExitAge, RX, RY,  ux, uy, width, height, u0)
        else:
            RTypes, RX, RY,  RExitAge, RIDs, RID, RTimeIn, RVals = bide.nonfluid_movement('resource', motion, Lists, RTimeIn, RExitAge, RX, RY,  ux, uy, width, height, u0)

    # Inflow of individuals (immigration)
    if ct == 1: SpeciesIDs, IndX, IndY,  MaintDict, EnvD, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion, 1000, 1, SpeciesIDs, IndX, IndY,  width, height, MaintDict, EnvD, envgrads, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u1, alpha, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)
    else: SpeciesIDs, IndX, IndY,  MaintDict, EnvD, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion, seedCom, m, SpeciesIDs, IndX, IndY,  width, height, MaintDict, EnvD, envgrads, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u1, alpha, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

    # dispersal
    Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
    if len(SpeciesIDs) > 0:
        if motion == 'fluid':
            SpeciesIDs, IndX, IndY,  IndExitAge, IndIDs, IndID, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.fluid_movement('individual', Lists, IndTimeIn, IndExitAge, IndX, IndY,  ux, uy, width, height, u0)
        else:
            SpeciesIDs, IndX, IndY,  IndExitAge, IndIDs, IndID, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.nonfluid_movement('individual', motion, Lists, IndTimeIn, IndExitAge, IndX, IndY,  ux, uy, width, height, u0)

    # Search for resources
    SpeciesIDs, Qs, IndIDs, ID, TimeIn, X, Y, GrowthDict, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.search(reproduction, speciation, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndX, IndY,  width, height, GrowthDict, DispDict, SpColorDict, N_RD, P_RD, C_RD, MaintDict, EnvD, envgrads, nNi, nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

    PRODI, PRODN, PRODC, PRODP = 0, 0, 0, 0

    p1, TNQ1, TPQ1, TCQ1 = metrics.getprod(Qs)

    # Consume
    RTypes, RVals, RIDs, RID, RTimeIn, RExitAge, RX, RY,  SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndX, IndY,  GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.consume(RTypes, RVals, RIDs, RID, RX, RY,  RTimeIn, RExitAge, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndX, IndY,  width, height, GrowthDict, N_RD, P_RD, C_RD, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

    # Reproduction
    SpeciesIDs, Qs, IndIDs, ID, TimeIn, X, Y, GrowthDict, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.reproduce(reproduction, speciation, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndX, IndY,  width, height, GrowthDict, DispDict, SpColorDict, N_RD, P_RD, C_RD, MaintDict, EnvD, envgrads, nNi, nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

    # maintenance
    SpeciesIDs, X, Y, IndExitAge, IndIDs, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.maintenance(SpeciesIDs, IndX, IndY,  IndExitAge, SpColorDict, MaintDict, EnvD, IndIDs, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

    p2, TNQ2, TPQ2, TCQ2 = metrics.getprod(Qs)

    PRODI = p2 - p1
    PRODN = TNQ2 - TNQ1
    PRODP = TPQ2 - TPQ1
    PRODC = TCQ2 - TCQ1

    # disturbance
    if np.random.binomial(1, disturb) == 1: SpeciesIDs, X, Y, IndExitAge, IndIDs, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.decimate(SpeciesIDs, IndX, IndY,  IndExitAge, SpColorDict, MaintDict, EnvD, IndIDs, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

    ax = fig.add_subplot(111)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    if len(SpeciesIDs) >= 1:  RAD, splist = bide.GetRAD(SpeciesIDs)
    else: RAD, splist, N, S = [], [], 0, 0

    N, S, tt, rr = sum(RAD), len(RAD), len(TIDs), len(RIDs)
    Ns.append(N)

    Title = ['Individuals consume resources, grow, reproduce, and die as they move through the environment. \nAverage speed on the x-axis is '+str(u0)+' units per time step. '+str(len(TExitAge))+' tracers have passed through.\nN: '+str(N)+', S: '+str(S)+', tracers: '+str(tt)+', resources: '+str(rr)+', replicate: '+str(len(Ns)), 'ct: '+str(ct)]

    txt.set_text(' '.join(Title))
    ax.set_ylim(0, height)
    ax.set_xlim(0, width)

    if plot_system == 'no':
        ##### PLOTTING THE SYSTEM ##############################################
        resource_scatImage.remove()
        tracer_scatImage.remove()
        Ind_scatImage.remove()
        colorlist = []
        sizelist = []
        for i, val in enumerate(SpeciesIDs):
            colorlist.append(SpColorDict[val])
            sizelist.append(mean(Qs[i]) * 1000)

        resource_scatImage = ax.scatter(RX, RY, s = RVals, c = 'w', edgecolor = 'SpringGreen', lw = 0.6, alpha=0.7)

        Ind_scatImage = ax.scatter(IndX, IndY, s = sizelist, c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.9)
        tracer_scatImage = ax.scatter(TX, TY, s = 200, c = 'r', marker='*', lw=0.0, alpha=0.6)

    plt.draw()

    if BurnIn == 'not done' and len(Ns) > 99:
        AugmentedDickeyFuller = sta.adfuller(Ns)
        val, p = AugmentedDickeyFuller[0:2]

        if p >= 0.05:
            Ns.pop(0)

        elif p < 0.05 or isnan(p) == True:
            BurnIn = 'done'
            Ns = [Ns[-1]] # only keep the most recent N value

    if BurnIn == 'done':

        PRODIs.append(PRODI)
        PRODNs.append(PRODN)
        PRODPs.append(PRODP)
        PRODCs.append(PRODC)

        if len(RExitAge) > 0:
            RTAUs.append(mean(RExitAge))
        if len(IndExitAge) > 0:
            INDTAUs.append(mean(IndExitAge))
        if len(TExitAge) > 0:
            TTAUs.append(mean(TExitAge))

        # Examining the resource RAD
        if len(RTypes) > 0:
            RRAD, Rlist = bide.GetRAD(RTypes)
            RDens = len(RTypes)/(height*width)
            RDiv = float(metrics.Shannons_H(RRAD))
            RRich = len(Rlist)

        RDENs.append(RDens)
        RDIVs.append(RDiv)
        RRICHs.append(RRich)
        # Number of tracers, resource particles, and individuals
        T, R, N = len(TIDs), len(RIDs), len(SpeciesIDs)

        Ts.append(T)
        Rs.append(R)
        Ss.append(S)

        if N >= 1:

            RAD, splist = bide.GetRAD(SpeciesIDs)
            RAD, splist = zip(*sorted(zip(RAD, splist), reverse=True))
            RAD = list(RAD)

            S = len(RAD)
            Ss.append(S)
            # Evenness, Dominance, and Rarity measures
            Ev = metrics.e_var(RAD)
            EVs.append(Ev)
            ES = metrics.e_simpson(RAD)
            ESs.append(ES)

            if len(Ns) == 1:
                splist2 = list(splist)

            if len(Ns) > 1:
                wt = metrics.WhittakersTurnover(splist, splist2)
                jc = metrics.jaccard(splist, splist2)
                so = metrics.sorensen(splist, splist2)
                splist2 = list(splist)
                WTs.append(wt)
                Jcs.append(jc)
                Sos.append(so)

            Nm, BP = [max(RAD), Nm/N]
            NMAXs.append(Nm)
            BPs.append(BP)

            SD = metrics.simpsons_dom(RAD)
            SDs.append(SD)
            sk = stats.skew(RAD)
            SKs.append(sk)

            Gs.append(mean(GrowthList))
            Ms.append(mean(MaintList))
            Ds.append(mean(DispList))

            means = [sum(x)/len(x) for x in zip(*N_RList)]
            NRs.append(mean(means))
            means = [sum(x)/len(x) for x in zip(*P_RList)]
            PRs.append(mean(means))
            means = [sum(x)/len(x) for x in zip(*C_RList)]
            CRs.append(mean(means))


        process = psutil.Process(os.getpid())
        mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)
        # return the memory usage in MB

        if len(Ns) > 99:
            t = time.clock() - t
            print sim, ' N:', int(round(mean(Ns))), 'S:', int(round(mean(Ss))), 'WT:', round(mean(WTs),2), ':  flow:', u0, 'time:', round(t,1), 'seconds', ' MB:',int(round(mem)), 'p-val =', round(p,3)
            t = time.clock()

            SString = str(splist).strip('()')
            RADString = str(RAD).strip('()')
            RADString = str(RAD).strip('[]')
            IndRTD = str(IndExitAge).strip('[]')
            TRTD = str(TExitAge).strip('[]')
            RRTD = str(RExitAge).strip('[]')

            OUT1 = open(GenPath + 'examples/SimData.csv','a')
            OUT2 = open(GenPath + 'examples/RADs.csv','a')
            OUT3 = open(GenPath + 'examples/Species.csv','a')
            OUT4 = open(GenPath + 'examples/IndRTD.csv','a')
            OUT5 = open(GenPath + 'examples/TracerRTD.csv','a')
            OUT6 = open(GenPath + 'examples/ResRTD.csv','a')

            outlist = [sim, motion, mean(PRODIs), mean(PRODNs), mean(PRODPs), mean(PRODCs), r, nNi, nP, nC, rmax, gmax, maintmax, dmax, barriers, alpha, seedCom, u0, width, height, viscosity, N, m, mean(RTAUs), mean(TTAUs), mean(INDTAUs), mean(RDENs), mean(RDIVs), mean(RRICHs), mean(Ss), mean(ESs), mean(EVs), mean(BPs), mean(SDs), mean(NMAXs), mean(SKs), T, R, speciation, mean(WTs), mean(Jcs), mean(Sos), mean(Gs), mean(Ms), mean(NRs), mean(PRs), mean(CRs), mean(Ds), amp, flux, freq, phase, disturb]
            outlist = str(outlist).strip('[]')

            print>>OUT1, outlist
            print>>OUT2, RADString
            print>>OUT3, SString
            print>>OUT4, ct1,',', sim,',', IndRTD
            print>>OUT5, ct1,',', sim,',', TRTD
            print>>OUT6, ct1,',', sim,',', RRTD

            OUT1.close()
            OUT2.close()
            OUT3.close()
            OUT4.close()
            OUT5.close()
            OUT6.close()

            ct1 += 1
            ct = 0
            SpColorDict, GrowthDict, MaintDict, EnvD, N_RD, P_RD, C_RD, RColorDict, DispDict = {}, {}, {}, {}, {}, {}, {}, {}, {}
            width, height, alpha, motion, reproduction, speciation, seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, amp, freq, flux, pulse, phase, disturb, envgrads, barriers, Rates = rp.get_rand_params()

            sim += 1
            if sim > num_sims:
                print "simplex finished"
                sys.exit()

            alpha = np.random.uniform(0.99, 0.999)

            for i in range(barriers):
                lefts.append(np.random.uniform(0.2, .8))
                bottoms.append(np.random.uniform(0.2, 0.8))

            u0 = choice(Rates)
            p = 0
            BurnIn = 'not done'

            RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, ct, IndID, RID, N, ct1, T, R, PRODI, PRODQ = [0]*21
            RAD, splist, IndTimeIn, SpeciesIDs, IndX, IndY,  IndIDs, Qs, IndExitAge, TX, TY, TExitAge, TIDs, TTimeIn, RX, RY,  RIDs, RTypes, RExitAge, RTimeIn, RVals, Gs, Ms, NRs, PRs, CRs, Ds, Ts, Rs, PRODIs, PRODNs, PRODPs, PRODCs, Ns, RTAUs, TTAUs, INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs,BPs, SDs, NMAXs, SKs, MUs, MAINTs, WTs, Jcs, Sos, splist2 = [list([]) for _ in xrange(53)]
            SpColorDict, GrowthDict, MaintDict, EnvD, N_RD, P_RD, C_RD, RColorDict, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
            SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, RColorList, DispList = [list([]) for _ in xrange(8)]

            #Ts, Rs, PRODIs, PRODNs, PRODPs, PRODCs, Ns, RTAUs, TTAUs, INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs = [list([]) for _ in xrange(22)]
            #RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, WT, Jc, So, ct, IndID, RID, N, ct1, T, R, PRODI, PRODQ = [0]*24

            #SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, RColorList, DispList, IndTimeIn, SpeciesIDs, IndX, IndY,  IndIDs, Qs, IndExitAge, splist2 = [list([]) for _ in xrange(16)]
            #TX, TY, TExitAge, TIDs, TTimeIn, WTs, Jcs, Sos, RX, RY,  RIDs, RTypes, RExitAge, RTimeIn, RVals = [list([]) for _ in xrange(15)]

            n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, bN, bS, bE, bW, bNE, bNW, bSE, bSW = LBM.SetLattice(u0, viscosity, width, height, lefts, bottoms, barriers)

            u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))

            ####################### REPLACE ENVIRONMENT ########################
            ax = fig.add_subplot(111)



################ DIMENSIONAL & MODEL CONSTANTS ##################################
width, height, alpha, motion, reproduction, speciation, seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, amp, freq, flux, pulse, phase, disturb, envgrads, barriers, Rates = rp.get_rand_params()
lefts, bottoms = [], []

for b in range(barriers):
    lefts.append(np.random.uniform(0.05, .95))
    bottoms.append(np.random.uniform(0.05, 0.95))

#######################  Ind COMMUNITY PARAMETERS  #########################
RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, ct, IndID, RID, N, ct1, T, R, PRODI, PRODQ = [0]*21
RAD, splist, IndTimeIn, SpeciesIDs, IndX, IndY,  IndIDs, Qs, IndExitAge, TX, TY, TExitAge, TIDs, TTimeIn, RX, RY,  RIDs, RTypes, RExitAge, RTimeIn, RVals, Gs, Ms, NRs, PRs, CRs, Ds, Ts, Rs, PRODIs, PRODNs, PRODPs, PRODCs, Ns, RTAUs, TTAUs, INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs,BPs, SDs, NMAXs, SKs, MUs, MAINTs, WTs, Jcs, Sos, splist2 = [list([]) for _ in xrange(53)]
SpColorDict, GrowthDict, MaintDict, EnvD, N_RD, P_RD, C_RD, RColorDict, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, RColorList, DispList = [list([]) for _ in xrange(8)]


###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
num_sims = 1000

LowerLimit, sim, left1, bottom1, left2, bottom2 = 30, 1, 0.2, 0.2, 0.6, 0.6
viscosity = 10 # unitless but required by an LBM model
u0 = choice(Rates)  # initial in-flow speed

############### INITIALIZE GRAPHICS ############################################
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111) # initiate first plot
Ind_scatImage = ax.scatter([0],[0], alpha=0)
tracer_scatImage = ax.scatter([0],[0], alpha=0)
resource_scatImage = ax.scatter([0],[0], alpha=0)

#####################  Lattice Boltzmann PARAMETERS  ###################
n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, bN, bS, bE, bW, bNE, bNW, bSE, bSW = LBM.SetLattice(u0, viscosity, width, height, lefts, bottoms, barriers)

Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

t = time.clock()
Ns = []
BurnIn = 'not done'
p = 0.0

ani = animation.FuncAnimation(fig, nextFrame, frames=110, interval=40, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/simplex/results/movies/examples/2015_10_05_1751_hydrobide.avi', bitrate=5000)
