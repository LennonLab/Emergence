from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d
#from mpl_toolkits.mplot3d import Axes3D

from random import choice
from scipy import stats
import numpy as np
import sys
import os
import psutil
import itertools

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics
sys.path.append(mydir + "GitHub/hydrobide/tools/LBM")
import LBM
sys.path.append(mydir + "GitHub/hydrobide/tools/bide")
import bide



def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.pyx """

    motion = choice(['fluid', 'random_walk'])
    D = choice([2, 2]) # number of spatial dimensions

    width = choice([5, 6, 7, 8, 9, 10])
    height = choice([5, 6, 7, 8, 9, 10])
    length = choice([5, 6, 7, 8, 9, 10])

    alpha = np.random.uniform(0.99, 0.999)
    reproduction = choice(['fission', 'sexual'])
    speciation = choice(['yes', 'no'])
    predators = choice([0, 1, 2, 4, 8])
    parasites = choice([0, 1, 2, 4, 8])
    env_gradient = choice(['no', 'yes'])

    seedcom = 1000 # size of starting community
    m = choice([1]) # individuals immigrating per time step
    r = choice([50, 100, 150, 200, 250, 300, 350, 400, 450, 500]) # resource particles flowing in per time step

    nNi = choice([2, 4, 6, 8, 10]) # maximum number of Nitrogen types
    nP = choice([2, 4, 6, 8, 10]) # maximum number of Phosphorus types
    nC = choice([2, 4, 6, 8, 10]) # maximum number of Carbon types

    rmax = choice([1000, 5000, 10000]) # maximum resource particle size

    gmax = choice([0.2, 0.3, 0.4, 0.5])
    dmax = choice([0.01, 0.05, 0.1, 0.5, 1.0])
    maintmax = choice([0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008])

    #width, height = 10, 10
    #motion = 'fluid'
    reproduction = 'fission'
    speciation = 'yes'
    #rmax = 10000
    #r = 400
    #gmax = 0.2
    #maintmax = 0.0001

    return [width, height, length, alpha, motion, D, reproduction, speciation, predators, parasites, env_gradient, seedcom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax]



def testlengths(TypeOf, function, Lists):
    vals = []
    for List in Lists: vals.append(len(List))
    if min(vals) != max(vals):
        print '\n'+TypeOf+': '+function+', list lengths are different sizes:', vals
        sys.exit()
    return



######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    plot_system = 'no'
    logdata = 'yes'
    global width, height, length, Rates, u0, shift, sign, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW
    global SpColorDict, GrowthDict, N_RD, P_RD, C_RD, DispDict, MaintDict, one9th, four9ths, one36th, barrier, gmax, dmax, maintmax
    global IndIDs, Qs, IndID, IndTimeIn, IndExitAge, IndXcoords, IndYcoords, IndZcoords, Ind_scatImage, SpeciesIDs

    global TracerYcoords, tracer_scatImage, TracerTimeIn, TracerZcoords, TracerIDs, TracerExitAge, TracerXcoords
    global ResTypes, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResVals, ResTimeIn, ResExitAge, resource_scatImage
    global avg_Q, avg_maint, avg_disp, avg_res, avgTau, avg_growth, avg_prey, avg_symb, avg_parasite

    global barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW
    global BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2

    global ct1, Mu, Maint, motion, D, reproduction, speciation, predators, parasites, symbionts
    global env_gradient, seedcom, m, r, nNi, nP, nC, rmax, sim, RAD, splist, N, ct, splist2, WTs
    global ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, T, R, LowerLimit, prod_i, prod_q, viscosity, alpha
    global Ts, Rs, PRODIs, Ns, TRACERTAUs, INDTAUs, RESDENs, RESDIVs, RESRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs, RESTAUs

    global PRODNs, PRODPs, PRODCs
    global GrowthList, MaintList, N_RList, P_RList, C_RList, DispList
    global Gs, Ms, NRs, PRs, CRs, Ds

    Lists = [GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, Qs, IndIDs]
    testlengths('individuals', 'any', Lists)

    for step in range(1): # adjust number of steps for smooth animation
        # inflow of tracers
        TracerIDs, TracerTimeIn, TracerXcoords, TracerYcoords, TracerZcoords = bide.NewTracers(motion, TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, TracerTimeIn, width, height, length, u0, D)

        # inflow of resources
        ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(motion, ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResTimeIn, r, rmax, nNi, nP, nC, width, height, length, u0, D)

        # immigration
        if ct == 0:
            SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion, seedcom, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, width, height, length, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u0, alpha, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)
        else:

            SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion, m,       SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, width, height, length, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u0, alpha, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)
        ct += 1

        if motion == 'fluid' or motion == 'conveyor':  # a 'conveyor' belt action wherein y-coordinates never change occurs when there is 0 turbulence

            # stream & collide
            nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])
            rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

            # dispersal
            Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
            if len(SpeciesIDs) > 0:
                SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndExitAge, IndIDs, IndID, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.fluid_movement('individual', Lists, IndTimeIn, IndExitAge, IndXcoords, IndYcoords, IndZcoords, ux, uy, width, height, u0)

            # resource flow
            Lists = [ResTypes, ResIDs, ResID, ResVals]
            if len(ResTypes) > 0:
                ResTypes, ResXcoords, ResYcoords, ResZcoords, ResExitAge, ResIDs, ResID, ResTimeIn, ResVals  = bide.fluid_movement('resource', Lists, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords, ux, uy, width, height, u0)

            # moving tracer particles
            if len(TracerIDs) > 0:
                TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, TracerExitAge, TracerTimeIn = bide.fluid_movement('tracer', TracerIDs, TracerTimeIn, TracerExitAge, TracerXcoords, TracerYcoords, TracerZcoords, ux, uy, width, height, u0)

        elif motion == 'random_walk' or motion == 'unidirectional':

            # Moving tracer particles
            if len(TracerIDs) > 0:
                TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, TracerExitAge, TracerTimeIn = bide.nonfluid_movement('tracer', motion, TracerIDs, TracerExitAge, TracerTimeIn, TracerXcoords, TracerYcoords, TracerZcoords, width, height, length, u0, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

            # Moving resource particles
            if len(ResTypes) > 0:
                Lists = [ResTypes, ResIDs, ResID, ResVals]
                ResTypes, ResXcoords, ResYcoords, ResZcoords, ResExitAge, ResIDs, ResID, ResTimeIn, ResVals = bide.nonfluid_movement('resource', motion, Lists, ResExitAge, ResTimeIn, ResXcoords, ResYcoords, ResZcoords, width, height, length, u0, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

            # Moving individuals
            if len(SpeciesIDs) > 0:
                Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
                SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndExitAge, IndIDs, IndID, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.nonfluid_movement('individual', motion, Lists, IndExitAge, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, u0, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)


        # consume & reproduce
        p1 = len(IndIDs)
        TNQ1 = 0
        TPQ1 = 0
        TCQ1 = 0

        for q in Qs:
            TNQ1 += q[0]
            TPQ1 += q[1]
            TCQ1 += q[2]

        ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.consume(ResTypes, ResVals, ResIDs, ResID, ResXcoords, ResYcoords, ResZcoords, ResTimeIn, ResExitAge, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, GrowthDict, N_RD, P_RD, C_RD, DispDict, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)
        SpeciesIDs, Qs, IndIDs, ID, TimeIn, Xcoords, Ycoords, Zcoords, GrowthDict, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.reproduce(reproduction, speciation, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, GrowthDict, DispDict, SpColorDict, N_RD, P_RD, C_RD, MaintDict, D, nNi, nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

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
        SpeciesIDs, Xcoords, Ycoords, Zcoords, IndExitAge, IndIDs, IndTimeIn, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.maintenance(SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndExitAge, SpColorDict, MaintDict, IndIDs, IndTimeIn, Qs, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)


    if D == 3: ax = fig.add_subplot(111, projection='3d')
    else: ax = fig.add_subplot(111)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    if len(SpeciesIDs) >= 1: RAD, splist = bide.GetRAD(SpeciesIDs)
    else: RAD, splist, N, S = [], [], 0, 0

    N, S, tt, rr = sum(RAD), len(RAD), len(TracerIDs), len(ResIDs)

    Title = ['Individuals consume resources, grow, reproduce, and die as they move through the environment.'
    '\nAverage speed on the x-axis is '+str(u0)+' units per time step. '+str(len(TracerExitAge))+' tracers have passed through.',
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
            sizelist.append(Qs[i] * 1000)

        if D == 2: resource_scatImage = ax.scatter(ResXcoords, ResYcoords, s = ResVals, c = 'w', edgecolor = 'SpringGreen', lw = 0.6, alpha=0.7)
        elif D == 3: resource_scatImage = ax.scatter(ResXcoords, ResYcoords, ResZcoords, s = ResVals, c = 'w', edgecolor = 'SpringGreen', lw = 0.6, alpha=0.2)

        if D == 2: Ind_scatImage = ax.scatter(IndXcoords, IndYcoords, s = sizelist, c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.9)
        elif D == 3: Ind_scatImage = ax.scatter(IndXcoords, IndYcoords, IndZcoords, s = Qs, c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.99)

        if D == 2: tracer_scatImage = ax.scatter(TracerXcoords, TracerYcoords, s = 200, c = 'r', marker='*', lw=0.0, alpha=0.6)
        elif D == 3: tracer_scatImage = ax.scatter(TracerXcoords, TracerYcoords, TracerZcoords, s = 200, c = 'r', marker='*', lw=0.0, alpha=0.8)

    plt.draw()

    if u0 >= 0.75: LowerLimit = 80
    elif u0 >= 0.5: LowerLimit = 10
    elif u0 >= 0.1: LowerLimit = 4
    elif u0 > 0.025: LowerLimit = 2
    else: LowerLimit = 2

    # Record model values and reset, or not
    if len(TracerExitAge) >= LowerLimit or ct > 100:
        ct = 95

        PRODIs.append(PRODI)
        PRODNs.append(PRODN)
        PRODPs.append(PRODP)
        PRODCs.append(PRODC)


        RESTAUs.append(np.mean(ResExitAge))
        INDTAUs.append(np.mean(IndExitAge))
        TRACERTAUs.append(np.mean(TracerExitAge))
        ResExitAge, IndExitAge, TracerExitAge = [],[],[]

        # Examining the resource RAD
        if len(ResTypes) > 0:
            ResRAD, Rlist = bide.GetRAD(ResTypes)
            ResDens = len(ResTypes)/(height*width)
            ResDiv = float(metrics.Shannons_H(ResRAD))
            ResRich = len(Rlist)

        RESDENs.append(ResDens)
        RESDIVs.append(ResDiv)
        RESRICHs.append(ResRich)

        # Number of tracers, resource particles, and individuals
        T, R, N = len(TracerIDs), len(ResIDs), len(SpeciesIDs)

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
        mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)    # return the memory usage in MB

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
            RESTAU = np.mean(RESTAUs)
            TRACERTAU = np.mean(TRACERTAUs)
            INDTAU = np.mean(INDTAUs)
            RESDENS = np.mean(RESDENs)
            RESDIV = np.mean(RESDIVs)
            RESRICH = np.mean(RESRICHs)
            S = np.mean(Ss)
            ES = np.mean(ESs)
            EV = np.mean(EVs)
            BP = np.mean(BPs)
            SD = np.mean(SDs)
            NMAX = np.mean(NMAXs)
            SK = np.mean(SKs)
            WT = np.mean(WTs)

            print sim, ' N:', int(round(N)), 'S:', int(round(S)), ' pI:', int(PRODI), 'WT:', round(WT,3), ':  flow:', u0, 'motion:',motion, ' MB:',int(round(mem))

            if logdata == 'yes':

                SString = str(splist).strip('()')
                RADString = str(RAD).strip('()')
                IndRTD = str(IndExitAge).strip('[]')
                TracerRTD = str(TracerExitAge).strip('[]')
                ResRTD = str(ResExitAge).strip('[]')

                OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/SimData.csv','a+')
                OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/RADs.csv','a+')
                OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/Species.csv','a+')
                OUT4 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/IndRTD.csv','a+')
                OUT5 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/TracerRTD.csv','a+')
                OUT6 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/ResRTD.csv','a+')

                # RowID, sim, motion, dimensions, ind.production, biomass.prod.N
                # biomass.prod.P, biomass.prod.C, res.inflow, N.types, P.types, C.types
                outlist = [ct1, sim, motion, D, PRODI, PRODN, PRODP, PRODC, r, nNi, nP, nC]

                # max.res.val, max.growth.rate, max.met.maint, max.active.dispersal, barrier.width, barrier.height
                outlist += [rmax, gmax, maintmax, dmax, BarrierWidth, BarrierHeight]

                # logseries.a, starting.seed, flow.rate, width, height, viscosity
                outlist += [alpha, seedcom, u0, width, height, viscosity]

                # total.abundance, immigration.rate, resource.tau, particle.tau, individual.tau
                # resource.concentration, shannons.resource.diversity
                outlist += [N, m, RESTAU, TRACERTAU, INDTAU, RESDENS, RESDIV]

                # resource.richness, species.richness, simpson.e, e.var, berger.parker, inv.simp.D
                # N.max, skew, tracer.particles, resource.particles
                outlist += [RESRICH, S, ES, EV, BP, SD, NMAX, SK, T, R]

                # speciation, species.turnover, avg.per.capita.growth, avg.per.capita.maint'
                outlist += [speciation, WT, fG, fM]

                # avg.per.capita.N.efficiency, avg.per.capita.P.efficiency
                # avg.per.capita.C.efficiency, avg.per.capita.active.dispersal
                outlist += [fNR, fPR, fCR, fD]

                outlist = str(outlist).strip('[]')

                print>>OUT1, outlist
                print>>OUT2, RADString
                print>>OUT3, SString
                print>>OUT4, ct1,',', sim,',', IndRTD
                print>>OUT5, ct1,',', sim,',', TracerRTD
                print>>OUT6, ct1,',', sim,',', ResRTD

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
                SpColorDict, GrowthDict, MaintDict, N_RD, P_RD, C_RD, ResColorDict, DispDict = {}, {}, {}, {}, {}, {}, {}, {}
                width, height, length, alpha, motion, D, reproduction, speciation, predators, parasites, env_gradient, seedcom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax = get_rand_params()


                sim += 1
                alpha = np.random.uniform(0.99, 0.999)
                print '\n'

            Rates = np.roll(Rates, -1, axis=0)
            u0 = Rates[0]  # initial in-flow speed

            Ts, Rs, PRODIs, PRODNs, PRODPs, PRODCs, Ns, RESTAUs, TRACERTAUs, INDTAUs, RESDENs, RESDIVs, RESRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint,WT = [0,0,0,0,0,0,0,0,0,0,0,0,0]
            SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, ResColorList, DispList = [[],[],[],[],[],[],[],[]]
            IndTimeIn, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndIDs, Qs, IndExitAge, splist2 = [[], [],[],[],[],[],[],[],[]]
            TracerXcoords, TracerYcoords, TracerZcoords, TracerExitAge, TracerIDs, TracerTimeIn, WTs = [[],[],[],[],[],[],[]]
            ResXcoords, ResYcoords, ResZcoords, ResIDs, ResTypes, ResExitAge, ResTimeIn, ResVals = [[],[],[],[],[],[],[],[]]

            if motion == 'fluid' or motion == 'conveyor':
                n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, left1, bottom1, left2, bottom2, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

            # inflow of resources
            ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(motion, ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResTimeIn, r, rmax, nNi, nP, nC, width, height, length, u0, D)
            # immigration
            SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, ID, TimeIn, Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = bide.immigration(dmax, gmax, maintmax, motion, seedcom, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, width, height, length, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u0, alpha, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList)

            ####################### REPLACE ENVIRONMENT ############################
            if D == 3: ax = fig.add_subplot(111, projection='3d')
            elif D == 2: ax = fig.add_subplot(111)

'''
############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/SimData.csv','w')
OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/RADs.csv','w')
OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/Species.csv','w')
OUT4 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/IndRTD.csv','w')
OUT5 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/TracerRTD.csv','w')
OUT6 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/25_Aug/ResRTD.csv','w')

# printing physical variables, residence times, community diversity properties, physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, motion, dimensions, ind.production, biomass.prod.N, biomass.prod.P, biomass.prod.C, res.inflow, N.types, P.types, C.types, max.res.val, max.growth.rate, max.met.maint, max.active.dispersal, barrier.width, barrier.height, logseries.a, starting.seed, flow.rate, width, height, viscosity, total.abundance, immigration.rate, resource.tau, particle.tau, individual.tau, resource.concentration, shannons.resource.diversity, resource.richness, species.richness, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, tracer.particles, resource.particles, speciation, species.turnover, avg.per.capita.growth, avg.per.capita.maint, avg.per.capita.N.efficiency, avg.per.capita.P.efficiency, avg.per.capita.C.efficiency, avg.per.capita.active.dispersal'

OUT1.close()
OUT2.close()
OUT3.close()
OUT4.close()
OUT5.close()
OUT6.close()
'''

################ DIMENSIONAL & MODEL CONSTANTS ##################################
width, height, length, alpha, motion, D, reproduction, speciation, predators, parasites, env_gradient, seedcom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax = get_rand_params()
motion = 'fluid'
#######################  Ind COMMUNITY PARAMETERS  #########################
ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = 0,0,0,0,0,0,0,0,0,0,0,0
ct, IndID, ResID, N, ct1, T, R, PRODI, PRODQ = 0,0,0,0,0,0,0,0,0

RAD, splist, IndTimeIn, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndIDs, Qs, IndExitAge = [],[],[],[],[],[],[],[],[],[]
TracerXcoords, TracerYcoords, TracerZcoords, TracerExitAge, TracerIDs, TracerTimeIn = [],[],[],[],[],[]
ResXcoords, ResYcoords, ResZcoords, ResIDs, ResTypes, ResExitAge, ResTimeIn, ResVals = [],[],[],[],[],[],[],[]
Gs, Ms, NRs, PRs, CRs, Ds = [],[],[],[],[],[]

Ts, Rs, PRODIs, PRODNs, PRODPs, PRODCs, Ns, RESTAUs, TRACERTAUs, INDTAUs, RESDENs, RESDIVs, RESRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
WTs, splist2 = [],[]

SpColorDict, GrowthDict, MaintDict, N_RD, P_RD, C_RD, ResColorDict, DispDict = {}, {}, {}, {}, {}, {}, {}, {}
SpColorList, GrowthList, MaintList, N_RList, P_RList, C_RList, ResColorList, DispList = [[],[],[],[],[],[],[],[]]

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
LowerLimit, shift, sign, sim = 30, 0.0, 0.1, 75
left1, bottom1, left2, bottom2 = 0.2, 0.2, 0.6, 0.6
BarrierWidth, BarrierHeight = 0.0, 0.0

BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = [],[],[],[]
viscosity = 10 # unitless but required by an LBM model

Rates = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.02, 0.01, 0.0075, 0.005])  # inflow speeds
#Rates = np.array([1.0, 0.75, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005])  # inflow speeds
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
        n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, left1, bottom1, left2, bottom2, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

        #BarrierImage1 = plt.bar(left1*width, BarrierHeight*height, BarrierWidth*width, bottom1*height, color = '0.3', edgecolor = '0.4', alpha=0.2)
        BarrierImage1 = plt.fill_between([left1*width, (left1+BarrierWidth)*width], bottom1*height, (bottom1+BarrierHeight)*height,
                color ='0.3', alpha= '0.4', linewidths=0.5, edgecolor='0.2')
        BarrierImage2 = plt.fill_between([left2*width, (left2+BarrierWidth)*width], bottom2*height, (bottom2+BarrierHeight)*height,
                color ='0.3', alpha= '0.4', linewidths=0.5, edgecolor='0.2')

elif D == 3:
    ax = fig.add_subplot(111, projection='3d')
    Ind_scatImage = ax.scatter([0],[0],[0], alpha=0.0)
    tracer_scatImage = ax.scatter([0],[0],[0], alpha=0.0)
    resource_scatImage = ax.scatter([0],[0],[0], alpha=0.0)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=110, interval=40, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/hydrobide/results/movies/2015_08_05_1741_hydrobide.avi', bitrate=5000)
