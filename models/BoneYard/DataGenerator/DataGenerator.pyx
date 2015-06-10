from __future__ import division
from random import randint, choice
import numpy as np
import sys
import os
from scipy import stats
import time
import psutil

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "/GitHub/hydrobide/tools/LBM")
import LBM
sys.path.append(mydir + "/GitHub/hydrobide/tools/bide")
import bide
sys.path.append(mydir + "tools/metrics")
import metrics


def get_mrrmax():
    """ Get random model parameter values. Others are chosen in bide.pyx """

    m = 1 #int(choice(range(1, 2))) # individuals immigrating per time step
    r = int(choice(range(10, 100))) # resource particles flowing in per time step
    nr = int(choice(range(1, 10))) # maximum number of resources types
    rmax = int(choice(range(10, 1000))) # maximum value of resource particle size

    return [m, r, nr, rmax]


#######################  MICROBE COMMUNITY PARAMETERS  #########################
m, r, nr, rmax = get_mrrmax()

MicID, ResID, N, S = 0, 0, 0, 0
avgTau = str()

COM, MicXcoords, MicYcoords, AvgTaus, RAD, splist = [], [], [], [], [], []
MicIDs, MicQs, MicExitAge, MicTimeIn = [], [], [], []
ResYcoords, RES, ResIDs, ResType, ResXcoords = [], [], [], [], []
TracerIDs, TracerXcoords, TracerYcoords, TracerExitAge = [], [], [], []

microbe_color_dict, GrowthDict, MaintDict = {}, {}, {}
ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}


###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
width = int(choice([5,6,7,8,9,10]))
height = int(choice([5, 5]))
shift, sign, sim, ct1, BarrierWidth, BarrierHeight = 0.0, 0.1, 0, 0, 0.2, 0.2

Rates = np.array([1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025])  # inflow speeds
u0 = Rates[0]  # initial in-flow speed


#####################  Lattice Boltzmann PARAMETERS  ###########################
viscosity =  0.84   # fluid viscosity of water
n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight)


############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/SimData.csv','w+')
OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/RADs.csv','w+')
OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/Species.csv','w+')
# printing physical variables, residence times, community diversity properties, physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, FlowRate, Width, Height, Viscosity, N, immigration.rate, particle.tau, cell.tau, resource.concentration, shannons.resource.diversity, S, resource.richness, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, avg.per.capita.growth, avg.per.capita.maint'
#             ct1,   sim,   u0,     width, height, viscosity, N, m,                  TracerTau,  MicrobeTau, ResDens,                ResDiv,                    S, ResRich,            ES,        Ev,    BP,            SD,         Nm,    sk,         Mu,               Maint
OUT1.close()
OUT2.close()
OUT3.close()

######################################## RUN SIMULATIONS  ######################
start = 1
#t = time.clock()
while sim < 1000:

    # stream
    nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])

    # collide
    rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

    y = np.random.binomial(1, u0/20)
    if y == 1:
        # immigration
        COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict = bide.immigration(m, COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr)
        # new tracers
        TracerIDs, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height)
        # inflow of resources
        RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height)


    # dispersal
    COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs = bide.dispersal(COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height, u0, MicIDs, MicID, MicTimeIn, MicQs, DispParamsDict)

    # resource flow
    RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.resFlow(RES, ux, uy, u0, ResXcoords, ResYcoords, width, height, ResID, ResIDs, r, rmax, nr, ResType)

    # moving tracer particles
    TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy = bide.MoveTracers(TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy)

    # consume and reproduce
    RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, ResType = bide.ConsumeAndReproduce(RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict, ResType, ResUseDict)

    # maintenance
    #COM, MicXcoords, MicYcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs = bide.maintenance(COM, MicXcoords, MicYcoords, MicExitAge, microbe_color_dict, MaintDict, MicIDs, MicID, MicTimeIn, MicQs)


    if len(TracerExitAge) >= 20:
        ct1 += 1

        N = len(COM)
        T = len(TracerIDs)
        R = len(RES)

        process = psutil.Process(os.getpid())
        mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)    # return the memory usage in MB

        # Physical and general community parameters
        OutList = [ct1, sim, u0, width, height, viscosity, N, m]

        # Residence times for tracers and microbes
        TracerTau = float(np.mean(TracerExitAge))
        MicrobeTau = float(np.mean(MicExitAge))
        OutList.extend([TracerTau, MicrobeTau])

        # Examining the resource RAD
        if len(ResType) > 0:
            ResRAD, Rlist = bide.GetRAD(ResType)
            ResDens = sum(RES)/(height*width)
            ResDiv = float(metrics.Shannons_H(ResRAD))
            ResRich = len(Rlist)
            OutList.extend([ResDens, ResDiv, ResRich])
        else: OutList.extend([0, 0, 0])


        if N == 0:
            S = 0
            print sim, '  N:', N, 'S:', S, ' T:', T,' R:', R, ' : flow rate:', u0, ' memory:',mem
            OutList.extend([S, 0, 0, 0, 0, 0, 0, 0, 0])

        else:
            RAD, splist = bide.GetRAD(COM)
            RAD, splist = zip(*sorted(zip(RAD, splist), reverse=True))

            if N != sum(RAD):
                print 'N != sum(RAD)'
                sys.exit()

            S = len(RAD)
            OutList.extend([S])

            if S == 1: OutList.extend([0, 0, 1.0, 0, N, 0, 0, 0])

            elif max(RAD) == min(RAD):

                SD = float(metrics.simpsons_dom(RAD))
                # Specific Growth rate and Maintenance
                Mu, Maint = 0, 0
                for i, sp in enumerate(splist):
                    Mu += RAD[i] * GrowthDict[sp]
                    Maint += RAD[i] * MaintDict[sp]

                Mu = float(Mu/S)
                Maint = float(Maint/S)
                OutList.extend([Mu, Maint])

                OutList.extend([1, 1, S/N, SD, N/S, 0, Mu, Maint])

            else:
                # Evenness, Dominance, and Rarity measures
                Ev = metrics.e_var(RAD)
                ES = float(metrics.e_simpson(RAD))
                Nm = max(RAD)
                BP = float(Nm/N)
                SD = float(metrics.simpsons_dom(RAD))
                sk = float(stats.skew(RAD))

                OutList.extend([ES, Ev, BP, SD, Nm, sk])

                # Specific Growth rate and Maintenance
                Mu, Maint = 0, 0
                for i, sp in enumerate(splist):
                    Mu += RAD[i] * GrowthDict[sp]
                    Maint += RAD[i] * MaintDict[sp]

                Mu = float(Mu/S)
                Maint = float(Maint/S)
                OutList.extend([Mu, Maint])

                process = psutil.Process(os.getpid())
                mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)    # return the memory usage in MB

        print sim, ' N:', N, 'S:', S, ' T:', T,' R:', R, ' : flow rate:', u0, ' memory:',mem

        OutString = str(OutList).strip('[]')
        SString = str(splist).strip('()')
        RADString = str(RAD).strip('()')
        OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/SimData.csv','a')
        OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/RADs.csv','a')
        OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/Species.csv','a')
        print>>OUT1, OutString
        print>>OUT2, RADString
        print>>OUT3, SString
        OUT1.close()
        OUT2.close()
        OUT3.close()


        if u0 == min(Rates):
            microbe_color_dict, GrowthDict, MaintDict = {}, {}, {}
            ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}
            width = int(choice([5,6,7,8,9,10]))
            height = int(choice([5,6]))
            m, r, nr, rmax = get_mrrmax()
            sim += 1
            print '\n'

        Rates = np.roll(Rates, -1, axis=0)
        u0 = Rates[0]  # initial in-flow speed

        MicTimeIn, COM, MicXcoords, MicYcoords, TracerXcoords, TracerYcoords, RES, ResXcoords, ResYcoords, ResIDs, ResType, MicIDs, MicQs, MicExitAge, TracerExitAge, TracerIDs = [list([]) for _ in xrange(16)]
        # Lattice Boltzmann PARAMETERS
        n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight)
