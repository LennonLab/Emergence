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
    #reproduction = choice(['fission', 'sexual'])
    speciation = choice(['yes', 'no'])
    predators = choice([0, 1, 2, 4, 8])
    parasites = choice([0, 1, 2, 4, 8])
    env_gradient = choice(['no', 'yes'])

    seedcom = 1000 # size of starting community
    m = choice([1]) # individuals immigrating per time step
    r = choice([50, 100, 150, 200]) # resource particles flowing in per time step
    nr = choice([1, 2, 4, 6, 8, 10]) # maximum number of resources types
    rmax = choice([200, 400, 600, 800, 1000]) # maximum resource particle size

    gmax = choice([0.1, 0.2, 0.3, 0.4])
    maintmax = choice([0.01, 0.02, 0.04, 0.05])
    width, height = 10, 10
    #motion = 'fluid'
    reproduction = 'fission'
    #speciation = 'no'
    #rmax = 200
    #r = 200


    return [width, height, length, alpha, motion, D, reproduction, speciation, predators, parasites, env_gradient, seedcom, m, r, nr, rmax, gmax, maintmax]



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
    global SpColorDict, GrowthDict, ResUseDict, DispDict, MaintDict, one9th, four9ths, one36th, barrier, gmax, maintmax
    global IndIDs, Qs, IndID, IndTimeIn, IndExitAge, IndXcoords, IndYcoords, IndZcoords, Ind_scatImage, SpeciesIDs

    global TracerYcoords, tracer_scatImage, TracerTimeIn, TracerZcoords, TracerIDs, TracerExitAge, TracerXcoords
    global ResTypes, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResVals, ResTimeIn, ResExitAge, resource_scatImage
    global avg_Q, avg_maint, avg_disp, avg_res, avgTau, avg_growth, avg_prey, avg_symb, avg_parasite

    global barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW
    global BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2

    global ct1, Mu, Maint, motion, D, reproduction, speciation, predators, parasites, symbionts
    global env_gradient, seedcom, m, r, nr, rmax, sim, RAD, splist, N, ct, splist2, WTs
    global ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, T, R, LowerLimit, prod_i, prod_q, viscosity, alpha
    global Ts, Rs, PRODIs, PRODQs, Ns, TRACERTAUs, INDTAUs, RESDENs, RESDIVs, RESRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs, RESTAUs

    for step in range(1): # adjust number of steps for smooth animation
        # inflow of tracers
        TracerIDs, TracerTimeIn, TracerXcoords, TracerYcoords, TracerZcoords = bide.NewTracers(motion, TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, TracerTimeIn, width, height, length, u0, D)

        # inflow of resources
        ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(motion, ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResTimeIn, r, rmax, nr, width, height, length, u0, D)

        # immigration
        if ct == 0:
            SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict = bide.immigration(gmax, maintmax, motion, seedcom, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, width, height, length, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict, nr, u0, alpha, D)
        else:

            SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict = bide.immigration(gmax, maintmax, motion, m, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, width, height, length, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict, nr, u0, alpha, D)
        ct += 1

        if motion == 'fluid' or motion == 'conveyor':  # a 'conveyor' belt action wherein y-coordinates never change occurs when there is 0 turbulence

            # stream & collide
            nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])
            rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

            # dispersal
            Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict]
            if len(SpeciesIDs) > 0: SpeciesIDs, IndXcoords, IndYcoords, IndExitAge, IndIDs, IndID, IndTimeIn, Qs = bide.fluid_movement('individual', Lists, IndTimeIn, IndExitAge, IndXcoords, IndYcoords, ux, uy, width, height, u0)

            # resource flow
            Lists = [ResTypes, ResIDs, ResID, ResVals]
            if len(ResTypes) > 0: ResTypes, ResXcoords, ResYcoords, ResExitAge, ResIDs, ResID, ResTimeIn, ResVals = bide.fluid_movement('resource', Lists, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ux, uy, width, height, u0)

            # moving tracer particles
            if len(TracerIDs) > 0: TracerIDs, TracerXcoords, TracerYcoords, TracerExitAge, TracerTimeIn = bide.fluid_movement('tracer', TracerIDs, TracerTimeIn, TracerExitAge, TracerXcoords, TracerYcoords, ux, uy, width, height, u0)

        elif motion == 'random_walk' or motion == 'unidirectional':

            # Moving tracer particles
            if len(TracerIDs) > 0: TracerIDs, TracerExitAge, TracerTimeIn, TracerXcoords, TracerYcoords, TracerZcoords = bide.nonfluid_movement('tracer', motion, TracerIDs, TracerExitAge, TracerTimeIn, TracerXcoords, TracerYcoords, TracerZcoords, width, height, length, u0, D)

            # Moving resource particles
            if len(ResTypes) > 0:
                Lists = [ResTypes, ResIDs, ResVals]
                Lists, ResExitAge, ResTimeIn, Xcoords, Ycoords, Zcoords = bide.nonfluid_movement('resource', motion, Lists, ResExitAge, ResTimeIn, ResXcoords, ResYcoords, ResZcoords, width, height, length, u0, D)
                ResTypes, ResIDs, ResVals = Lists
            # Moving individuals
            if len(SpeciesIDs) > 0:
                Lists = [SpeciesIDs, IndIDs, Qs, DispDict]
                Lists, IndExitAge, IndTimeIn, Xcoords, Ycoords, Zcoords = bide.nonfluid_movement('individual', motion, Lists, IndExitAge, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, u0, D)
                SpeciesIDs, IndIDs, Qs = Lists

        # consume & reproduce
        p1, q1 = [len(IndIDs), sum(Qs)]
        ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords = bide.consume(ResTypes, ResVals, ResIDs, ResID, ResXcoords, ResYcoords, ResZcoords, ResTimeIn, ResExitAge, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, GrowthDict, ResUseDict, DispDict, D)
        SpeciesIDs, Qs, indIDs, ID, TimeIn, Xcoords, Ycoords, Zcoords, GrowthDict, DispDict = bide.reproduce(gmax, maintmax, reproduction, speciation, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, GrowthDict, DispDict, SpColorDict, ResUseDict, MaintDict, D, nr)
        PRODI, PRODQ = [len(IndIDs) - p1, sum(Qs) - q1]

        # maintenance
        SpeciesIDs, Xcoords, Ycoords, Zcoords, IndExitAge, IndIDs, IndTimeIn, Qs = bide.maintenance(SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndExitAge, SpColorDict, MaintDict, IndIDs, IndTimeIn, Qs, D)


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

    #print np.mean(ResExitAge), np.mean(IndExitAge), np.mean(TracerExitAge)
    # Record model values and reset, or not
    if len(TracerExitAge) >= LowerLimit or ct > 100:
        ct = 95

        PRODIs.append(PRODI)
        PRODQs.append(PRODQ)

        RESTAUs.append(np.mean(ResExitAge))
        INDTAUs.append(np.mean(IndExitAge))
        TRACERTAUs.append(np.mean(TracerExitAge))
        ResExitAge, IndExitAge, TracerExitAge = [],[],[]

        # Examining the resource RAD
        if len(ResTypes) > 0:
            ResRAD, Rlist = bide.GetRAD(ResTypes)
            ResDens = sum(ResTypes)/(height*width)
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

            # Specific Growth rate and Maintenance

            mu, maint = [0.0, 0.0]
            for i, sp in enumerate(splist):
                mu += (RAD[i] * GrowthDict[sp])
                maint += (RAD[i] * MaintDict[sp])

            Mu, Maint = [mu/N, maint/N]
            MUs.append(Mu)
            MAINTs.append(Maint)

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

        process = psutil.Process(os.getpid())
        mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)    # return the memory usage in MB

        if len(Ns) >= 2 and len(splist) > 0:
            T, R, PRODI, PRODQ, N, RESTAU, TRACERTAU, INDTAU, RESDENS, RESDIV, RESRICH, S, ES, EV, BP, SD, NMAX, SK, MU, MAINT, WT = [np.mean(Ts), np.mean(Rs), np.mean(PRODIs), np.mean(PRODQs), np.mean(Ns), np.mean(RESTAUs), np.mean(TRACERTAUs), np.mean(INDTAUs), np.mean(RESDENs), np.mean(RESDIVs), np.mean(RESRICHs), np.mean(Ss), np.mean(ESs), np.mean(EVs), np.mean(BPs), np.mean(SDs), np.mean(NMAXs), np.mean(SKs), np.mean(MUs), np.mean(MAINTs), np.mean(WTs)]
            print sim, ' N:', int(round(N)), 'S:', int(round(S)), ' pI:', int(PRODI), 'WT:', round(WT,3), ':  flow:', u0, 'motion:',motion, ' MB:',int(round(mem))

            if logdata == 'yes':

                SString = str(splist).strip('()')
                RADString = str(RAD).strip('()')
                IndRTD = str(IndExitAge).strip('[]')
                TracerRTD = str(TracerExitAge).strip('[]')
                ResRTD = str(ResExitAge).strip('[]')

                OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/SimData.csv','a')
                OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/RADs.csv','a')
                OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/Species.csv','a')
                OUT4 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/IndRTD.csv','a')
                OUT5 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/TracerRTD.csv','a')
                OUT6 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/ResRTD.csv','a')

                print>>OUT1, ct1,',', sim,',', motion,',', D,',', PRODI,',', PRODQ,',', r,',', nr,',', rmax,',', BarrierWidth,',', BarrierHeight,',', alpha,',', seedcom,',', LowerLimit,',', u0,',', width,',', height,',', viscosity,',', N,',', m,',',RESTAU,',',TRACERTAU,',', INDTAU,',', RESDENS,',', RESDIV,',', RESRICH,',', S,',', ES,',', EV,',', BP,',', SD,',', NMAX,',', SK,',', MU,',', MAINT,',',T,',',R,',',speciation,',',WT,',',gmax,',',maintmax
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
                SpColorDict, GrowthDict, MaintDict, ResUseDict, ResColorDict, DispDict = {}, {}, {}, {}, {}, {}
                width, height, length, alpha, motion, D, reproduction, speciation, predators, parasites, env_gradient, seedcom, m, r, nr, rmax, gmax, maintmax = get_rand_params()
                sim += 1
                alpha = np.random.uniform(0.99, 0.999)
                print '\n'

            Rates = np.roll(Rates, -1, axis=0)
            u0 = Rates[0]  # initial in-flow speed

            Ts, Rs, PRODIs, PRODQs, Ns, RESTAUs, TRACERTAUs, INDTAUs, RESDENs, RESDIVs, RESRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs = [[], [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint,WT = [0,0,0,0,0,0,0,0,0,0,0,0,0]
            IndTimeIn, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndIDs, Qs, IndExitAge, splist2 = [[], [],[],[],[],[],[],[],[]]
            TracerXcoords, TracerYcoords, TracerZcoords, TracerExitAge, TracerIDs, TracerTimeIn, WTs = [[],[],[],[],[],[],[]]
            ResXcoords, ResYcoords, ResZcoords, ResIDs, ResTypes, ResExitAge, ResTimeIn, ResVals = [[],[],[],[],[],[],[],[]]

            if motion == 'fluid' or motion == 'conveyor':
                n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, left1, bottom1, left2, bottom2, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

            # inflow of resources
            ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(motion, ResTypes, ResVals, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResTimeIn, r, rmax, nr, width, height, length, u0, D)
            # immigration
            SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, MaintDict, GrowthDict, DispDict, SpColorDict, indIDs, ID, TimeIn, Qs, ResUseDict = bide.immigration(gmax, maintmax, motion, seedcom, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, width, height, length, MaintDict, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict, nr, u0, alpha, D)

            ####################### REPLACE ENVIRONMENT ############################
            if D == 3: ax = fig.add_subplot(111, projection='3d')
            elif D == 2: ax = fig.add_subplot(111)


"""
############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/SimData.csv','w')
OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/RADs.csv','w')
OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/Species.csv','w')
OUT4 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/IndRTD.csv','w')
OUT5 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/TracerRTD.csv','w')
OUT6 = open(mydir + '/GitHub/hydrobide/results/simulated_data/2015_August/11_Aug/ResRTD.csv','w')

# printing physical variables, residence times, community diversity properties, physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, motion, D, ind.prod, biomass.prod, res.inflow, res.types, max.res.val, barrier.width, barrier.height, logseries.a, starting.seed, lower.limit, FlowRate, Width, Height, Viscosity, N, immigration.rate, resource.tau, particle.tau, individual.tau, resource.concentration, shannons.resource.diversity, resource.richness, S, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, avg.per.capita.growth, avg.per.capita.maint, tracer.particles, resource.particles, speciation, species.turnover, growth.max, maint.max'
#             ct,   sim, motion, D, PRODI,    PRODQ,        r,           nr,       rmax,        BarrierWidth,  BarrierHeight,  alpha,       seedcom,       LowerLimit,  u0,       width, height, viscosity, N, m,                RESTAU,       TRACERTAU,    INDTAU,         RESDENS,                RESDIV,                      RESRICH,           S, ES,        EV,    BP,            SD,         NMAX,  SK,   MU,                    MAINT,                T,                R,                 speciation,   species turnover

OUT1.close()
OUT2.close()
OUT3.close()
OUT4.close()
OUT5.close()
OUT6.close()
"""


################ DIMENSIONAL & MODEL CONSTANTS ##################################
width, height, length, alpha, motion, D, reproduction, speciation, predators, parasites, env_gradient, seedcom, m, r, nr, rmax, gmax, maintmax = get_rand_params()
motion = 'fluid'
#######################  Ind COMMUNITY PARAMETERS  #########################
ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = 0,0,0,0,0,0,0,0,0,0,0,0
ct, IndID, ResID, N, ct1, T, R, PRODI, PRODQ = 0,0,0,0,0,0,0,0,0

RAD, splist, IndTimeIn, SpeciesIDs, IndXcoords, IndYcoords, IndZcoords, IndIDs, Qs, IndExitAge = [],[],[],[],[],[],[],[],[],[]
TracerXcoords, TracerYcoords, TracerZcoords, TracerExitAge, TracerIDs, TracerTimeIn = [],[],[],[],[],[]
ResXcoords, ResYcoords, ResZcoords, ResIDs, ResTypes, ResExitAge, ResTimeIn, ResVals = [],[],[],[],[],[],[],[]

Ts, Rs, PRODIs, PRODQs, Ns, RESTAUs, TRACERTAUs, INDTAUs, RESDENs, RESDIVs, RESRICHs, Ss, ESs, EVs, BPs, SDs, NMAXs, SKs, MUs, MAINTs = [], [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
WTs, splist2 = [],[]

SpColorDict, GrowthDict, MaintDict, ResUseDict, ResColorDict, DispDict = {}, {}, {}, {}, {}, {}

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
LowerLimit, shift, sign, sim = 30, 0.0, 0.1, 26
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
