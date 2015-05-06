from __future__ import division
from random import randint, choice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os
import scipy as sc
from scipy.optimize import curve_fit
from scipy import stats


mydir = os.path.expanduser("~/")
sys.path.append(mydir + "/GitHub/hydrobide/tools/LBM")
import LBM
sys.path.append(mydir + "/GitHub/hydrobide/tools/bide")
import bide

sys.path.append(mydir + "tools/metrics")
import metrics


def get_mrrmax():
    m = int(choice(range(1, 10))) # number of individuals immigrating per time step
    r = int(choice(range(1, 10))) # number of resource particles flowing in per time step
    nr = int(choice(range(1, 10))) # max number of resources types
    rmax = int(choice(range(1, 100))) # max value of resource particle sizes
    return [m, r, nr, rmax]

#######################  MICROBE COMMUNITY PARAMETERS  #########################
COM = []
Xs = []
MicXcoords = []
MicYcoords = []
MicIDs = []
MicQs = []
MicExitAge = []

MicID = 0
MicTimeIn = []

microbe_color_dict = {}
GrowthDict = {}
MaintDict = {}
ResUseDict = {}
ResColorDict = {}
DispParamsDict = {}

m, r, nr, rmax = get_mrrmax()

AvgTaus = []
avgTau = str()

########################################## RESOURCE PARAMETERS #################
RES = []
ResIDs = []
ResType = []
ResID = 0
ResYcoords = []
ResXcoords = []

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
width = int(choice([5,6,7,8,9,10]))
height = int(choice([5,6,7,8,9,10]))
shift = 0.0
sign = 0.1

BarrierXcoords1 = []
BarrierYcoords1 = []
BarrierXcoords2 = []
BarrierYcoords2 = []
Rates = np.array([1.0, 0.5, 0.25, 0.1, 0.075, 0.05])  # inflow speeds
u0 = float(choice(Rates))  # initial in-flow speed

#######################  Inert tracer particles  ###############################
TracerXcoords = []
TracerYcoords = []
TracerIDs = []
TracerExitAge = []
TracerIDs, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0)

#####################  Lattice Boltzmann PARAMETERS  ###########################
BarrierWidth = 0.2
BarrierHeight = 0.2
viscosity =  u0*10   	                                       # fluid viscosity

omega = 1 / (3 * viscosity + 0.5)	                  # relaxation parameter
four9ths = 4.0/9.0	    # abbreviations for lattice-Boltzmann weight factors
one9th   = 1.0/9.0
one36th  = 1.0/36.0

n0 = four9ths * (np.ones((height,width)) - 1.5*u0**2)
nN = one9th * (np.ones((height,width)) - 1.5*u0**2)
nS = one9th * (np.ones((height,width)) - 1.5*u0**2)
nE = one9th * (np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nW = one9th * (np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nNE = one36th * (np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nSE = one36th * (np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nNW = one36th * (np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nSW = one36th * (np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW       # macroscopic density
ux = (nE + nNE + nSE - nW - nNW - nSW) / rho	        # macroscopic x velocity
uy = (nN + nNE + nNW - nS - nSE - nSW) / rho	        # macroscopic y velocity

barrier = np.zeros((height, width), bool)                  # Initialize barriers

barrierN = np.roll(barrier,  1, axis=0)           # sites just north of barriers
barrierS = np.roll(barrier, -1, axis=0)           # sites just south of barriers
barrierE = np.roll(barrier,  1, axis=1)           # etc.
barrierW = np.roll(barrier, -1, axis=1)
barrierNE = np.roll(barrierN,  1, axis=1)
barrierNW = np.roll(barrierN, -1, axis=1)
barrierSE = np.roll(barrierS,  1, axis=1)
barrierSW = np.roll(barrierS, -1, axis=1)


############## OPEN OUTPUT DATA FILE AND RUN SIMULATIONS  ######################
OUT = open(mydir + '/GitHub/hydrobide/results/simulated_data/DataGenerator.csv','w+')

print>>OUT, 'RowID, FlowRate, Width, Height, Viscosity, N, S, TracerParticle_ResTime, MicrobeCell_ResTime, Sp_Evenness, BergerParker, InvSimpDiversity, Nmax, Skew, AvgPerCapita_GrowthRate, AvgPerCapita_Maintenance, Immigration_Rate, Resource_Concentration, Shannons_ResourceDiversity, ResourceRichness'


ct = 0
ct1 = 0
sims = 0
while sims < 10: # adjust number of steps for smooth animation

    # immigration
    COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height, MaintDict, u0, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict = bide.immigration(m, COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height, MaintDict, u0, GrowthDict, DispParamsDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr)

    # stream
    nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])

    # collide
    rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, omega, u0 = LBM.collide(rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, omega, u0, four9ths, one9th, one36th)

    # dispersal
    COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs = bide.dispersal(COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height, u0, MicIDs, MicID, MicTimeIn, MicQs, DispParamsDict)

    # resource flow
    RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.resFlow(RES, ux, uy, u0, ResXcoords, ResYcoords, width, height, ResID, ResIDs, r, rmax, nr, ResType)

    # moving tracer particles
    TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy = bide.MoveTracers(TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy)
    TracerIDs, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0)

    # consume and reproduce
    RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, ResType = bide.ConsumeAndReproduce(RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict, ResType, ResUseDict)

    # maintenance
    COM, MicXcoords, MicYcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs = bide.maintenance(COM, MicXcoords, MicYcoords, MicExitAge, microbe_color_dict, MaintDict, MicIDs, MicID, MicTimeIn, MicQs)


    if len(TracerExitAge) >= 10: # Many if/else statements, but it should ultimately save time
        S = len(COM)

        if S > 3:
            RAD, splist = bide.GetRAD(COM)

            if max(RAD) > 1:
                ct += 1
                N = sum(COM)

                # Physical and general community parameters
                OutList = [ct1, u0, width, height, viscosity, m, N, S]
                ct1 += 1

                # Residence times for tracers and microbes
                TracerTau = float(np.mean(TracerExitAge))
                MicrobeTau = float(np.mean(MicExitAge))
                OutList.extend([TracerTau, MicrobeTau])

                # Evenness, Dominance, and Rarity measures
                #Ev = metrics.e_var(RAD)
                ES = float(metrics.e_simpson(RAD))
                #EQ = metrics.EQ(RAD)
                #O = metrics.OE(RAD)
                #EN = metrics.NHC(RAD)
                #EP = metrics.e_pielou(RAD)
                #EH = metrics.e_heip(RAD)
                Nm = max(RAD)
                BP = float(Nm/N)
                SD = float(metrics.simpsons_dom(RAD))
                #Mc = metrics.McNaughton(RAD)
                sk = float(stats.skew(RAD))
                #lsk = metrics.Rlogskew(RAD)

                OutList.extend([ES, BP, SD, Nm, sk])

                # Specific Growth rate and Maintenance
                Mu, Maint = 0, 0
                for i, sp in enumerate(splist):
                    Mu += RAD[i] * GrowthDict[sp]
                    Maint += RAD[i] * MaintDict[sp]

                Mu = float(Mu/S)
                Maint = float(Maint/S)

                # Examining the resource RAD
                ResRAD, Rlist = bide.GetRAD(ResType)
                ResDens = sum(RES)/(height*width)
                ResDiv = float(metrics.Shannons_H(ResRAD))
                ResRich = len(Rlist)

                OutList.extend([Mu, Maint, m, ResDens, ResDiv, ResRich])
                OutString = str(OutList).strip('[]')

                print>>OUT, OutString
                print 'N:', sum(COM), 'S:',len(COM), '  V:', round(height*width,1),'Q:', u0, ' m:',m,' r:',r,' : ', sims, ct

                Xs = []
                TracerExitAge = []
                MicExitAge = []

                if ct >= 4:

                    width = int(choice([5,6,7,8,9,10]))
                    height = int(choice([5,6,7,8,9,10]))

                    m, r, nr, rmax = get_mrrmax()

                    Rates = np.roll(Rates, -1, axis=0)
                    u0 = Rates[0]  # initial in-flow speed
                    viscosity =  u0*10  # fluid viscosity
                    omega = 1 / (3 * viscosity + 0.5) # relaxation parameter

                    microbe_color_dict = {}
                    GrowthDict = {}
                    MaintDict = {}
                    ResUseDict = {}
                    ResColorDict = {}
                    DispParamsDict = {}

                    MicTimeIn, COM, MicXcoords, MicYcoords, TracerXcoords, TracerYcoords, RES, ResXcoords, ResYcoords, ResIDs, ResType, MicIDs, MicQs, MicExitAge, TracerExitAge, TracerIDs = [list([]) for _ in xrange(16)]
                    args = LBM.SetLattice(n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, one9th, four9ths, one36th, rho, u0, width, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, height, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight)
                    n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, one9th, four9ths, one36th, rho, u0, width, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, height, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight = args

                    sims += 1
                    ct = 0

OUT.close()
