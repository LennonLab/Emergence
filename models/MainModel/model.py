from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d
#from mpl_toolkits.mplot3d import Axes3D

from random import choice
from scipy import stats
import numpy as np
from numpy import sin, pi
import sys
import os
import psutil
#import itertools

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics
sys.path.append(mydir + "GitHub/hydrobide/tools/LBM")
import LBM
sys.path.append(mydir + "GitHub/hydrobide/tools/bide")
import bide



def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.pyx """

    #motion = choice(['fluid', 'random_walk', 'conveyor', 'search'])
    motion = choice(['fluid', 'random_walk'])
    D = choice([2, 2]) # number of spatial dimensions

    width = choice([10, 12, 14, 16, 18, 20])
    height = choice([5, 6, 7, 8, 9, 10])
    length = choice([5, 6, 7, 8, 9, 10])

    barriers = choice([2, 4, 6, 8, 16, 32])
    barriers = 4
    #if motion == 'conveyor': barriers = 0

    pulse = choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    flux = choice(['yes'])

    # Sine wave: y(t) = A * sin(2*pi*f*t + phi)
    # let phi = 0, meaning 0 amplitude at time 0
    amp = choice([0.05, 0.1, 0.2, 0.3, 0.4, 0.5]) # A
    freq = choice([0.1, 0.08, 0.06, 0.04, 0.02, 0.01]) # f
    phase = choice([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    # 0 = in phase; 16 = entirely out of phase

    disturb = choice([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1])

    alpha = np.random.uniform(0.99, 0.999)
    reproduction = choice(['fission', 'sexual'])
    speciation = choice(['yes', 'no'])

    seed = choice([1000]) # size of starting community
    m = choice([1]) # m = individuals immigrating per time step
    r = choice([50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
    # r = resource particles flowing in per time step

    nNi = choice([1, 3, 6, 12, 14, 16, 18, 20]) # max number of Nitrogen types
    nP = choice([1, 3, 6, 12, 14, 16, 18, 20]) # max number of Phosphorus types
    nC = choice([1, 3, 6, 12, 14, 16, 18, 20]) # max number of Carbon types

    envgrads = []
    num_envgrads = choice([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    for i in range(num_envgrads):
        g = choice(['horizontal', 'vertical', 'aggregated'])
        if g == 'aggregated':
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([g, x, y])
        else:
            x = choice([0,1])
            envgrads.append([g, x, 1-x])


    rmax = choice([1000, 5000, 10000, 15000]) # maximum resource particle size

    gmax = choice([0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    dmax = choice([0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9])
    maintmax = choice([0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002])

    reproduction = 'fission'

    return [width, height, length, alpha, motion, D, reproduction, speciation, \
            seed, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, amp, freq, \
            flux, pulse, phase, disturb, envgrads, barriers]



def testlengths(TypeOf, function, Lists):
    vals = []
    for List in Lists:
        vals.append(len(List))

    if min(vals) != max(vals):
        print '\n'+TypeOf+': '+function+', lists have different lengths:', vals
        sys.exit()

    return



######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    plot_system = 'yes'
    logdata = 'no'
    global width, height, length, Rates, u0, rho, ux, uy, n0, nN
    global nS, nE, nW, nNE, nNW, nSE, nSW, SpColorDict, GrowthDict, N_RD
    global P_RD, C_RD, DispDict, MaintDict, one9th, four9ths, one36th, barrier
    global gmax, dmax, maintmax, IndIDs, Qs, IndID, IndTimeIn, IndExitAge
    global IndX, IndY, IndZ, Ind_scatImage, SpeciesIDs, EnvD

    global TY, tracer_scatImage, TTimeIn, TZ
    global TIDs, TExitAge, TX, RTypes, RX
    global RY, RZ, RID, RIDs, RVals, RTimeIn
    global RExitAge, resource_scatImage

    global bN, bS, bE, bW, bNE, bNW, bSE, bSW
    global ct1, Mu, Maint, motion, D, reproduction, speciation
    global seed, m, r, nNi, nP, nC, rmax, sim, RAD, splist, N, ct, splist2, WTs
    global RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, T, R, LowerLimit
    global prod_i, prod_q, viscosity, alpha, Ts, Rs, PRODIs, Ns, TTAUs, INDTAUs
    global RDENs, RDIVs, RRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs

    global PRODNs, PRODPs, PRODCs, lefts, bottoms, Gs, Ms, NRs,PRs,CRs,Ds, RTAUs
    global GrowthList, MaintList, N_RList, P_RList, C_RList, DispList
    global amp, freq, flux, pulse, phase, disturb, envgrads, barriers

    Lists = [GrowthList, MaintList, N_RList, P_RList, C_RList,
                DispList, SpeciesIDs, IndX, IndY, IndZ, Qs, IndIDs]
    testlengths('individuals', 'any', Lists)

    u1 = float(u0)
    if flux == 'yes':
        # fluctuate flow rate according to randomly chosen values for:
        # amplitude, frequency, and phase

        u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))
        if u1 > 1:
            u1 == 1.0

    for step in range(1): # adjust number of steps for smooth animation
        # inflow of tracers
        TIDs, TTimeIn, TX, TY, TZ = bide.NewTracers(motion,\
        TIDs, TX, TY, TZ, TTimeIn, width, height, length, u0, D)

        # inflow of resources
        RTypes, RVals, RX, RY, RZ, RIDs, RID, \
        RTimeIn = bide.ResIn(motion, RTypes, RVals, RX, RY,\
        RZ, RID, RIDs, RTimeIn, r, rmax, nNi, nP, nC, width, height,\
        length, u1, D)

        # immigration
        if ct == 0:
            SpeciesIDs, IndX, IndY, IndZ, MaintDict, EnvD, \
            GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, \
            Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, \
            C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion,\
            seed, SpeciesIDs, IndX, IndY, IndZ, width, \
            height, length, MaintDict, EnvD, GrowthDict, DispDict, SpColorDict, \
            IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u1, \
            alpha, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

        else:

            SpeciesIDs, IndX, IndY, IndZ, MaintDict, EnvD, \
            GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, \
            Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, \
            C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion, \
            m, SpeciesIDs, IndX, IndY, IndZ, width, height, \
            length, MaintDict, EnvD, GrowthDict, DispDict, SpColorDict, IndIDs, \
            IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u1, alpha, \
            D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)
        ct += 1

        if motion == 'fluid' or motion == 'conveyor':
            # a 'conveyor' belt action wherein y-coordinates never
            # change occurs when there is 0 turbulence

            # stream & collide
            nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier = LBM.stream([nN, nS,\
                    nE, nW, nNE, nNW, nSE, nSW, barrier])

            rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, \
            nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE,\
                                nNW, nSE, nSW, u0)

            # dispersal
            Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList, \
                                MaintList, N_RList, P_RList, C_RList, DispList]

            if len(SpeciesIDs) > 0:
                SpeciesIDs, IndX, IndY, IndZ, IndExitAge, \
                IndIDs, IndID, IndTimeIn, Qs, GrowthList, MaintList, N_RList, \
                P_RList, C_RList, DispList = bide.fluid_movement('individual', \
                        Lists, IndTimeIn, IndExitAge, IndX, IndY, \
                        IndZ, ux, uy, width, height, u0)

            # resource flow
            Lists = [RTypes, RIDs, RID, RVals]
            if len(RTypes) > 0:
                RTypes, RX, RY, RZ, RExitAge, \
                RIDs, RID, RTimeIn, \
                RVals  = bide.fluid_movement('resource', Lists, RTimeIn,\
                            RExitAge, RX, RY, RZ, \
                            ux, uy, width, height, u0)

            # moving tracer particles
            if len(TIDs) > 0:
                TIDs, TX, TY, TZ, \
                TExitAge, TTimeIn = bide.fluid_movement('tracer',\
                        TIDs, TTimeIn, TExitAge, TX,\
                        TY, TZ, ux, uy, width, height, u0)


        elif motion == 'random_walk' or motion == 'unidirectional':

            # Moving tracer particles
            if len(TIDs) > 0:
                TIDs, TX, TY, TZ, \
                TExitAge, TTimeIn = bide.nonfluid_movement('tracer',\
                        motion, TIDs, TExitAge, TTimeIn,\
                        TX, TY, TZ, width,\
                        height, length, u0, D, GrowthList, MaintList, N_RList,\
                        P_RList, C_RList, DispList)

            # Moving resource particles
            if len(RTypes) > 0:
                Lists = [RTypes, RIDs, RID, RVals]
                RTypes, RX, RY, RZ, RExitAge, \
                RIDs, RID, RTimeIn, \
                RVals = bide.nonfluid_movement('resource', motion, Lists, \
                        RExitAge, RTimeIn, RX, RY, \
                        RZ, width, height, length, u0, D, GrowthList,\
                        MaintList, N_RList, P_RList, C_RList, DispList)

            # Moving individuals
            if len(SpeciesIDs) > 0:
                Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList,\
                        MaintList, N_RList, P_RList, C_RList, DispList]

                SpeciesIDs, IndX, IndY, IndZ, IndExitAge, \
                IndIDs, IndID, IndTimeIn, Qs, GrowthList, MaintList, N_RList,\
                P_RList, C_RList, DispList = bide.nonfluid_movement('individual',\
                        motion, Lists, IndExitAge, IndTimeIn, IndX, \
                        IndY, IndZ, width, height, length, u0,\
                        D, GrowthList, MaintList, N_RList, P_RList, C_RList,\
                        DispList)


        # consume & reproduce
        p1 = len(IndIDs)
        TNQ1 = 0
        TPQ1 = 0
        TCQ1 = 0

        for q in Qs:
            TNQ1 += q[0]
            TPQ1 += q[1]
            TCQ1 += q[2]

        RTypes, RVals, RIDs, RID, RTimeIn, RExitAge, RX, RY, RZ, SpeciesIDs, \
        Qs, IndIDs, IndID, IndTimeIn, IndX, IndY, IndZ, GrowthList, MaintList,\
        N_RList, P_RList, C_RList, DispList = bide.consume(RTypes, RVals, RIDs,\
        RID, RX, RY, RZ, RTimeIn, RExitAge, SpeciesIDs, Qs, IndIDs, IndID,\
        IndTimeIn, IndX, IndY, IndZ, width, height, length, GrowthDict, N_RD,\
        P_RD, C_RD, DispDict, D, GrowthList, MaintList, N_RList, P_RList,\
        C_RList, DispList)

        SpeciesIDs, Qs, IndIDs, ID, TimeIn, X, Y, Z, GrowthDict, DispDict,\
        GrowthList, MaintList, N_RList, P_RList, C_RList,\
        DispList = bide.reproduce(reproduction, speciation, SpeciesIDs, Qs,\
        IndIDs, IndID, IndTimeIn, IndX, IndY, IndZ, width, height, length,\
        GrowthDict, DispDict, SpColorDict, N_RD, P_RD, C_RD, MaintDict, EnvD, D,\
        nNi, nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

        TNQ2 = 0
        TPQ2 = 0
        TCQ2 = 0

        for q in Qs:
            TNQ2 += q[0]
            TPQ2 += q[1]
            TCQ2 += q[2]

        PRODI = len(IndIDs) - p1
        PRODN = TNQ2 - TNQ1
        PRODP = TPQ2 - TPQ1
        PRODC = TCQ2 - TCQ1

        # maintenance
        SpeciesIDs, X, Y, Z, IndExitAge, IndIDs, IndTimeIn, Qs, GrowthList,\
        MaintList, N_RList, P_RList, C_RList, \
        DispList = bide.maintenance(SpeciesIDs, IndX, IndY, IndZ,\
        IndExitAge, SpColorDict, MaintDict, EnvD, IndIDs, IndTimeIn, Qs, D,\
        GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

        # disturbance
        d = np.random.binomial(1, disturb)
        if d == 1:
            SpeciesIDs, X, Y, Z, IndExitAge, IndIDs, IndTimeIn, Qs,\
            GrowthList, MaintList, N_RList, P_RList, C_RList,\
            DispList = bide.decimate(SpeciesIDs, IndX, IndY, IndZ,\
            IndExitAge, SpColorDict, MaintDict, EnvD, IndIDs, IndTimeIn, Qs,\
            D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)


    if D == 3: ax = fig.add_subplot(111, projection='3d')
    else: ax = fig.add_subplot(111)
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                left='off', right='off', labelbottom='off', labelleft='off')

    if len(SpeciesIDs) >= 1: RAD, splist = bide.GetRAD(SpeciesIDs)
    else: RAD, splist, N, S = [], [], 0, 0

    N, S, tt, rr = sum(RAD), len(RAD), len(TIDs), len(RIDs)

    Title = ['Individuals consume resources, grow, reproduce, and die as they move through the environment.'
    '\nAverage speed on the x-axis is '+str(u0)+' units per time step. '+str(len(TExitAge))+' tracers have passed through.',
           '\nMotion is '+motion+'; N: '+str(N)+', S: '+str(S)+', tracers: '+str(tt)+', resources: '+str(rr)+', replicate: '+str(len(Ns)), 'ct: '+str(ct)]

    txt.set_text(' '.join(Title))
    ax.set_ylim(0, height)
    ax.set_xlim(0, width)
    if D == 3:
        ax.set_zlim(0,length)

    if plot_system == 'yes':

        ##### PLOTTING THE SYSTEM ############################################
        resource_scatImage.remove()
        tracer_scatImage.remove()
        Ind_scatImage.remove()

        colorlist = []
        sizelist = []
        for i, val in enumerate(SpeciesIDs):
            colorlist.append(SpColorDict[val])
            sizelist.append(np.mean(Qs[i]) * 1000)

        if D == 2:
            resource_scatImage = ax.scatter(RX, RY, \
            s = RVals, c = 'w', edgecolor = 'SpringGreen', lw = 0.6, alpha=0.7)
        elif D == 3:
            resource_scatImage = ax.scatter(RX, RY, \
            RZ, s = RVals, c = 'w', edgecolor = 'SpringGreen',\
            lw = 0.6, alpha=0.2)

        if D == 2:
            Ind_scatImage = ax.scatter(IndX, IndY, s = sizelist,\
            c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.9)
        elif D == 3:
            Ind_scatImage = ax.scatter(IndX, IndY, IndZ,\
            s = Qs, c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.99)

        if D == 2:
            tracer_scatImage = ax.scatter(TX, TY,\
            s = 200, c = 'r', marker='*', lw=0.0, alpha=0.6)
        elif D == 3:
            tracer_scatImage = ax.scatter(TX, TY,\
            TZ, s = 200, c = 'r', marker='*', lw=0.0, alpha=0.8)

    plt.draw()

    if u0 >= 0.75: LowerLimit = 80
    elif u0 >= 0.5: LowerLimit = 10
    elif u0 >= 0.1: LowerLimit = 4
    elif u0 > 0.025: LowerLimit = 2
    else: LowerLimit = 2

    # Record model values and reset, or not
    if len(TExitAge) >= LowerLimit or ct > 100:
        ct = 95

        PRODIs.append(PRODI)
        PRODNs.append(PRODN)
        PRODPs.append(PRODP)
        PRODCs.append(PRODC)


        RTAUs.append(np.mean(RExitAge))
        INDTAUs.append(np.mean(IndExitAge))
        TTAUs.append(np.mean(TExitAge))
        RExitAge, IndExitAge, TExitAge = [],[],[]

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

        if N == 0:
            Ns.append(0)
            Ss.append(0)

        if N >= 1:

            RAD, splist = bide.GetRAD(SpeciesIDs)
            RAD, splist = zip(*sorted(zip(RAD, splist), reverse=True))
            S = len(RAD)
            Ss.append(S)
            Ns.append(N)

            # Evenness, Dominance, and Rarity measures
            Ev = metrics.e_var(RAD)
            EVs.append(Ev)

            ES = metrics.e_simpson(RAD)
            ESs.append(ES)

            if len(Ns) == 1:
                splist2 = list(splist)

            if len(Ns) > 1:
                wt = metrics.WhittakersTurnover(splist, splist2)
                splist2 = list(splist)
                WTs.append(wt)

            Nm, BP = [max(RAD), Nm/N]
            NMAXs.append(Nm)
            BPs.append(BP)

            SD = metrics.simpsons_dom(RAD)
            SDs.append(SD)
            sk = stats.skew(RAD)
            SKs.append(sk)

            Gs.append(np.mean(GrowthList))
            Ms.append(np.mean(MaintList))
            Ds.append(np.mean(DispList))

            means = [sum(x)/len(x) for x in zip(*N_RList)]
            NRs.append(np.mean(means))
            means = [sum(x)/len(x) for x in zip(*P_RList)]
            PRs.append(np.mean(means))
            means = [sum(x)/len(x) for x in zip(*C_RList)]
            CRs.append(np.mean(means))

        process = psutil.Process(os.getpid())
        mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)
        # return the memory usage in MB

        if len(Ns) >= 2:

            fG = np.mean(Gs)
            fM = np.mean(Ms)
            fNR = np.mean(NRs)
            fPR = np.mean(PRs)
            fCR = np.mean(CRs)
            fD = np.mean(Ds)

            T = np.mean(Ts)
            R = np.mean(Rs)
            PRODI = np.mean(PRODIs)
            PRODN = np.mean(PRODNs)
            PRODP = np.mean(PRODPs)
            PRODC = np.mean(PRODCs)
            N = np.mean(Ns)
            RTAU = np.mean(RTAUs)
            TTAU = np.mean(TTAUs)
            INDTAU = np.mean(INDTAUs)
            RDENS = np.mean(RDENs)
            RDIV = np.mean(RDIVs)
            RRICH = np.mean(RRICHs)
            S = np.mean(Ss)
            ES = np.mean(ESs)
            EV = np.mean(EVs)
            BP = np.mean(BPs)
            SD = np.mean(SDs)
            NMAX = np.mean(NMAXs)
            SK = np.mean(SKs)
            WT = np.mean(WTs)

            print sim, ' N:', int(round(N)), 'S:', int(round(S)), ' pI:', \
            int(PRODI), 'WT:', round(WT,3), ':  flow:', u0, 'motion:',motion, \
            ' MB:',int(round(mem))

            if logdata == 'yes':

                SString = str(splist).strip('()')
                RADString = str(RAD).strip('()')
                IndRTD = str(IndExitAge).strip('[]')
                TRTD = str(TExitAge).strip('[]')
                RRTD = str(RExitAge).strip('[]')

                OUT1 = open(GenPath + '2015_September/12_Sept/SimData.csv','a')
                OUT2 = open(GenPath + '2015_September/12_Sept/RADs.csv','a')
                OUT3 = open(GenPath + '2015_September/12_Sept/Species.csv','a')
                OUT4 = open(GenPath + '2015_September/12_Sept/IndRTD.csv','a')
                OUT5 = open(GenPath + '2015_September/12_Sept/TracerRTD.csv','a')
                OUT6 = open(GenPath + '2015_September/12_Sept/ResRTD.csv','a')

                # RowID, sim, motion, dimensions, ind.production, biomass.prod.N
                # biomass.prod.P, biomass.prod.C, res.inflow, N.types, P.types,
                # C.types
                outlist = [ct1, sim, motion, D, PRODI, PRODN, PRODP, PRODC, r,\
                                                                    nNi, nP, nC]

                # max.res.val, max.growth.rate, max.met.maint,
                # max.active.dispersal, barrier.width, barrier.height
                outlist += [rmax, gmax, maintmax, dmax, barriers]

                # logseries.a, starting.seed, flow.rate, width, height,viscosity
                outlist += [alpha, seed, u0, width, height, viscosity]

                # total.abundance, immigration.rate, resource.tau, particle.tau,
                # individual.ta, resource.concentration, shannons.res.diversity
                outlist += [N, m, RTAU, TTAU, INDTAU, RDENS, RDIV]

                # resource.richness, species.richness, simpson.e, e.var,
                # berger.parker, inv.simp.D, N.max, skew, tracer.particles,
                # resource.particles
                outlist += [RRICH, S, ES, EV, BP, SD, NMAX, SK, T, R]

                # speciation, species.turnover, avg.per.capita.growth,
                # avg.per.capita.maint'
                outlist += [speciation, WT, fG, fM]

                # avg.per.capita.N.efficiency, avg.per.capita.P.efficiency
                # avg.per.capita.C.efficiency, avg.per.capita.active.dispersal
                outlist += [fNR, fPR, fCR, fD]

                # amplitude, flux, frequency, phase, disturbance
                outlist += [amp, flux, freq, phase, disturb]

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

        if len(Ns) >= 2:
            ct1 += 1
            ct = 0

            if u0 == min(Rates):
                SpColorDict, GrowthDict, MaintDict, EnvD, N_RD, P_RD, C_RD, \
                RColorDict, DispDict = {}, {}, {}, {}, {}, {}, {}, {}, {}

                width, height, length, alpha, motion, D, reproduction, \
                speciation, seed, m, r, nNi, nP, nC, rmax, gmax, maintmax, \
                dmax, amp, freq, flux, pulse, phase, disturb, envgrads, \
                barriers = get_rand_params()


                sim += 1
                alpha = np.random.uniform(0.99, 0.999)
                print '\n'

            for i in range(barriers):
                lefts.append(np.random.uniform(0.05, .95))
                bottoms.append(np.random.uniform(0.05, 0.95))

            Rates = np.roll(Rates, -1, axis=0)
            u0 = Rates[0]  # initial in-flow speed

            Ts, Rs, PRODIs, PRODNs, PRODPs, PRODCs, Ns, RTAUs, TTAUs,\
            INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs,\
                    SKs, MUs, MAINTs = [list([]) for _ in xrange(22)]

            RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint,WT = [0]*13

            SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, \
                    RColorList, DispList = [list([]) for _ in xrange(8)]

            IndTimeIn, SpeciesIDs, IndX, IndY, IndZ, IndIDs,\
                    Qs, IndExitAge, splist2 = [list([]) for _ in xrange(9)]

            TX, TY, TZ, TExitAge,\
            TIDs, TTimeIn, WTs = [list([]) for _ in xrange(7)]

            RX, RY, RZ, RIDs, RTypes, \
            RExitAge, RTimeIn, RVals = [list([]) for _ in xrange(8)]

            if motion == 'fluid' or motion == 'conveyor':
                n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy,\
                bN, bS, bE, bW, bNE, bNW, bSE, bSW, BarrierX1, BarrierY1,\
                BarrierX2, BarrierY2 = LBM.SetLattice(u0, viscosity, width, \
                height, lefts, bottoms, barriers)

            u1 = 0.0
            if flux == 'yes':
                u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))

            # inflow of resources
            RTypes, RVals, RX, RY, RZ, RIDs, RID, RTimeIn = bide.ResIn(motion,\
            RTypes, RVals, RX, RY, RZ, RID, RIDs, RTimeIn, r, rmax, nNi, nP,\
            nC, width, height, length, u0, D)

            # immigration
            SpeciesIDs, IndX, IndY, IndZ, MaintDict, EnvD, GrowthDict, DispDict,\
            SpColorDict, IndIDs, ID, TimeIn, Qs, N_RD, P_RD, C_RD, GrowthList,\
            MaintList, N_RList, P_RList, C_RList,\
            DispList = bide.immigration(dmax, gmax, maintmax, motion, seed,\
            SpeciesIDs, IndX, IndY, IndZ, width, height, length, MaintDict, EnvD, \
            GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs,\
            N_RD, P_RD, C_RD, nNi, nP, nC, u0, alpha, D, GrowthList, \
            MaintList, N_RList, P_RList, C_RList, DispList)

            ####################### REPLACE ENVIRONMENT ########################
            if D == 3: ax = fig.add_subplot(111, projection='3d')
            elif D == 2: ax = fig.add_subplot(111)

GenPath = mydir + '/GitHub/hydrobide/results/simulated_data/'

############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(GenPath + '2015_September/12_Sept/SimData.csv','w')
OUT2 = open(GenPath + '2015_September/12_Sept/RADs.csv','w')
OUT3 = open(GenPath + '2015_September/12_Sept/Species.csv','w')
OUT4 = open(GenPath + '2015_September/12_Sept/IndRTD.csv','w')
OUT5 = open(GenPath + '2015_September/12_Sept/TracerRTD.csv','w')
OUT6 = open(GenPath + '2015_September/12_Sept/ResRTD.csv','w')

# printing physical variables, residence times, community diversity properties
# physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, motion, dimensions, ind.production, biomass.prod.N, biomass.prod.P, biomass.prod.C, res.inflow, N.types, P.types, C.types, max.res.val, max.growth.rate, max.met.maint, max.active.dispersal, barrier.width, barrier.height, logseries.a, starting.seed, flow.rate, width, height, viscosity, total.abundance, immigration.rate, resource.tau, particle.tau, individual.tau, resource.concentration, shannons.resource.diversity, resource.richness, species.richness, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, tracer.particles, resource.particles, speciation, species.turnover, avg.per.capita.growth, avg.per.capita.maint, avg.per.capita.N.efficiency, avg.per.capita.P.efficiency, avg.per.capita.C.efficiency, avg.per.capita.active.dispersal, amplitude, flux, frequency, phase, disturbance'

OUT1.close()
OUT2.close()
OUT3.close()
OUT4.close()
OUT5.close()
OUT6.close()


################ DIMENSIONAL & MODEL CONSTANTS ##################################
width, height, length, alpha, motion, D, reproduction, speciation, seed, m, r,\
nNi, nP, nC, rmax, gmax, maintmax, dmax, amp, freq, flux, pulse, phase, disturb,\
envgrads, barriers = get_rand_params()

lefts = []
bottoms = []

for b in range(barriers):
    lefts.append(np.random.uniform(0.05, .95))
    bottoms.append(np.random.uniform(0.05, 0.95))

#######################  Ind COMMUNITY PARAMETERS  #########################
RDens, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = 0,0,0,0,0,0,0,0,0,0,0,0
ct, IndID, RID, N, ct1, T, R, PRODI, PRODQ = 0,0,0,0,0,0,0,0,0

RAD, splist, IndTimeIn, SpeciesIDs, IndX, IndY, IndZ, IndIDs, Qs, \
IndExitAge, TX, TY, TZ, TExitAge, TIDs, TTimeIn, RX, RY, RZ, RIDs, RTypes,\
RExitAge, RTimeIn, RVals, Gs, Ms, NRs, PRs, CRs, Ds, Ts, Rs, PRODIs, PRODNs,\
PRODPs, PRODCs, Ns, RTAUs, TTAUs, INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs,\
BPs, SDs, NMAXs, SKs, MUs, MAINTs, WTs, splist2 = [list([]) for _ in xrange(54)]


SpColorDict, GrowthDict, MaintDict, EnvD, N_RD, P_RD, C_RD, RColorDict, \
        DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}

SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, RColorList, \
        DispList = [list([]) for _ in xrange(8)]

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
LowerLimit, sim = 30, 2
left1, bottom1, left2, bottom2 = 0.2, 0.2, 0.6, 0.6

BarrierX1, BarrierY1, BarrierX2, BarrierY2 = [],[],[],[]
viscosity = 10 # unitless but required by an LBM model

Rates = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.1, 0.09, 0.08, \
            0.07, 0.06, 0.05, 0.04, 0.02, 0.01, 0.0075, 0.005])  # inflow speeds
u0 = Rates[0]  # initial in-flow speed

############### INITIALIZE GRAPHICS ############################################
fig = plt.figure(figsize=(12, 8))

if D == 2:
    ax = fig.add_subplot(111) # initiate first plot
    Ind_scatImage = ax.scatter([0],[0], alpha=0)
    tracer_scatImage = ax.scatter([0],[0], alpha=0)
    resource_scatImage = ax.scatter([0],[0], alpha=0)

    if motion == 'fluid' or motion == 'conveyor':
        #####################  Lattice Boltzmann PARAMETERS  ###################
        n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, bN, bS, \
        bE, bW, bNE, bNW, bSE, bSW = LBM.SetLattice(u0, viscosity, \
        width, height, lefts, bottoms, barriers)

elif D == 3:
    ax = fig.add_subplot(111, projection='3d')
    Ind_scatImage = ax.scatter([0],[0],[0], alpha=0.0)
    tracer_scatImage = ax.scatter([0],[0],[0], alpha=0.0)
    resource_scatImage = ax.scatter([0],[0],[0], alpha=0.0)
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                left='off', right='off', labelbottom='off', labelleft='off')

Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=110, interval=40, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/hydrobide/results/movies/2015_08_05_1741_hydrobide.avi', bitrate=5000)
