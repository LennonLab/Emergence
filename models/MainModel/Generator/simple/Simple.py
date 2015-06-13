from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d
#from mpl_toolkits.mplot3d import Axes3D

from scipy import stats
import numpy as np
import sys
import os
import psutil

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "/tools/metrics")
import metrics
sys.path.append(mydir + "/GitHub/hydrobide/tools/LBM")
import LBM
sys.path.append(mydir + "/GitHub/hydrobide/tools/bide/bide_test")
import bide_test as bide


######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    global width, height, Rates, u0, shift, sign, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW
    global SpColorDict, GrowthDict, ResUseDict, DispParamsDict, MaintDict, one9th, four9ths, one36th, barrier
    global IndIDs, Qs, IndID, IndTimeIn, IndExitAge, IndXcoords, IndYcoords, Ind_scatImage

    global TracerYcoords, tracer_scatImage, TracerTimeIn, TracerIDs, TracerExitAge, TracerXcoords
    global ResTypes, ResXcoords, ResYcoords, ResID, ResIDs, ResVals, ResTimeIn, ResExitAge, resource_scatImage
    global avg_Q, avg_maint, avg_disp, avg_res, avgTau, avg_growth

    global barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW
    global BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2

    global ct1, Mu, Maint, motion, reproduction, mutation, predators, parasites, symbionts
    global env_gradient, seedcom, m, r, nr, rmax, sim, RAD, splist, N, TracerTau, IndTau
    global ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, T, R, stop, prod_i, prod_q, SpeciesIDs, viscosity, alpha

    for step in range(1): # adjust number of steps for smooth animation

        # inflow of tracers
        TracerIDs, TracerTimeIn, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, TracerTimeIn, width, height, u0)

        # inflow of resources
        ResTypes, ResVals, ResXcoords, ResYcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(ResTypes, ResVals, ResXcoords, ResYcoords, ResID, ResIDs, ResTimeIn, r, rmax, nr, width, height, u0)

	# immigration
        SpeciesIDs, IndXcoords, IndYcoords, MaintDict, GrowthDict, DispParamDict, SpColorDict, IDs, ID, TimeIn, Qs, ResUseDict = bide.immigration(m, SpeciesIDs, IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict, nr, u0, alpha)

        # stream
        nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])

        # collide
        rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

        # dispersal
        List = SpeciesIDs, IndIDs, ID, Qs
        if len(SpeciesIDs) > 0: SpeciesIDs, IndXcoords, IndYcoords, IndExitAge, IndIDs, IndID, IndTimeIn, Qs = bide.fluid_movement('individual', List, IndTimeIn, IndExitAge, IndXcoords, IndYcoords, ux, uy, width, height, u0)

        # resource flow
        List = ResTypes, ResIDs, ResID, ResVals
        if len(ResTypes) > 0: ResTypes, ResXcoords, ResYcoords, ResExitAge, ResIDs, ResID, ResTimeIn, ResVals = bide.fluid_movement('resource', List, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ux, uy, width, height, u0)

        # moving tracer particles
        List = TracerIDs
        if len(TracerIDs) > 0: TracerIDs, TracerXcoords, TracerYcoords, TracerExitAge, TracerTimeIn = bide.fluid_movement('tracer', List, TracerTimeIn, TracerExitAge, TracerXcoords, TracerYcoords, ux, uy, width, height, u0)

        # consume and reproduce
        prod_i = len(IndIDs)
        prod_q = sum(Qs)

        ResTypes, ResIDs, ResXcoords, ResYcoords, SpeciesIDs, IndIDs, IndID, IndTimeIn, Qs, IndXcoords, IndYcoords, ResVals = bide.ConsumeAndReproduce(ResTypes, ResIDs, ResXcoords, ResYcoords, SpeciesIDs, IndIDs, IndID, IndTimeIn, Qs, IndXcoords, IndYcoords, ResVals, width, height, GrowthDict, ResUseDict)

        #ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords = ResLists
        #SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords = IndLists

        prod_i = len(IndIDs) - prod_i
        prod_q = sum(Qs) - prod_q

        # maintenance
        SpeciesIDs, IndXcoords, IndYcoords, IndExitAge, IndIDs, IndTimeIn, Qs = bide.maintenance(SpeciesIDs, IndXcoords, IndYcoords, IndExitAge, SpColorDict, MaintDict, IndIDs, IndTimeIn, Qs)


    ########## plot the system #################################################
    fig.add_subplot(1,1,1)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    if len(SpeciesIDs) >= 1:
        RAD, splist = bide.GetRAD(SpeciesIDs)
        N, S = len(SpeciesIDs), len(RAD)

    else:
        RAD, splist, N, S = [], [], 0, 0

    Title = ['Inds consume resources, grow, reproduce, and die as they flow through a fluid environment. Average speed',
           '\non the x-axis is '+str(u0)+' units per time step. '+str(len(TracerExitAge))+' tracers have passed through.',
           'N = '+str(N)+', S = '+str(S)+'.'
           '\nOpen circles are resource particles. Semi-impermeable barriers (grey bars) produce turbulence.']

    txt.set_text(' '.join(Title))
    plt.ylim(0, height)
    plt.xlim(0, width)
    plt.draw()

    ##### PLOTTING THE INDIVIDUALS ############################################
    resource_scatImage.remove()
    resource_scatImage = plt.scatter(ResXcoords, ResYcoords, c = 'w', edgecolor = 'SpringGreen', s = ResVals, lw = 0.6, alpha=0.7)

    tracer_scatImage.remove()
    tracer_scatImage = plt.scatter(TracerXcoords, TracerYcoords, c = 'r', marker='*', lw=0.0, s = 200, alpha=0.6)

    Ind_scatImage.remove()
    colorlist = []
    for i, val in enumerate(SpeciesIDs): colorlist.append(SpColorDict[val])
    Ind_scatImage = plt.scatter(IndXcoords, IndYcoords, c = colorlist, edgecolor = '0.2', s = Qs, lw = 0.2, alpha=0.9)

    # Record model values and reset, or not
    if len(TracerExitAge) >= stop:

        # Examining the resource RAD
        if len(ResTypes) > 0:
            ResRAD, Rlist = bide.GetRAD(ResTypes)
            ResDens = sum(ResTypes)/(height*width)
            ResDiv = float(metrics.Shannons_H(ResRAD))
            ResRich = len(Rlist)

        # Residence times for tracers and Inds
        TracerTau = float(np.mean(TracerExitAge))
        IndTau = float(np.mean(IndExitAge))

        T = len(TracerIDs)
        R = len(ResTypes)
        N = len(SpeciesIDs)

        if N >= 1:
            RAD, splist = bide.GetRAD(SpeciesIDs)
            RAD, splist = zip(*sorted(zip(RAD, splist), reverse=True))
            S = len(RAD)

            # Specific Growth rate and Maintenance

            mu, maint = [0.0, 0.0]
            for i, sp in enumerate(splist):
                mu = RAD[i] * GrowthDict[sp]
                maint = RAD[i] * MaintDict[sp]

            Mu = mu/N
            Maint = maint/N

            # Evenness, Dominance, and Rarity measures
            Ev = float(metrics.e_var(RAD))
            ES = float(metrics.e_simpson(RAD))
            Nm = max(RAD)
            BP = float(Nm/N)
            SD = float(metrics.simpsons_dom(RAD))
            sk = float(stats.skew(RAD))

        process = psutil.Process(os.getpid())
        mem = round(process.get_memory_info()[0] / float(2 ** 20), 1)    # return the memory usage in MB

        print sim, ' N:', N, 'S:', S, ' pI:', round(prod_i,1), 'pQ:', round(prod_q,2), ': flow:', u0, ' MB:',mem

        SString = str(splist).strip('()')
        RADString = str(RAD).strip('()')
        OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/SimData.csv','a')
        OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/RADs.csv','a')
        OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/Species.csv','a')
        print>>OUT1, ct1,',', sim,',', prod_i,',', prod_q,',', r,',', nr,',', rmax,',', BarrierWidth,',', BarrierHeight,',', alpha,',', seedcom,',', stop,',', u0,',', width,',', height,',', viscosity,',', N,',', m,',', TracerTau,',', IndTau,',', ResDens,',', ResDiv,',', ResRich,',', S,',', ES,',', Ev,',', BP,',', SD,',', Nm,',', sk,',', Mu,',', Maint
        print>>OUT2, RADString
        print>>OUT3, SString
        OUT1.close()
        OUT2.close()
        OUT3.close()

        if u0 == min(Rates):
            SpColorDict, GrowthDict, MaintDict = {}, {}, {}
            ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}
            width, height, alpha, motion, seedcom, m, r, nr, rmax = bide.get_rand_params()
            sim += 1
            alpha = np.random.uniform(0.9, 0.999)
            print '\n'

        Rates = np.roll(Rates, -1, axis=0)
        u0 = Rates[0]  # initial in-flow speed

        TracerTau, IndTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        IndTimeIn, SpeciesIDs, IndXcoords, IndYcoords, IndIDs, Qs, IndExitAge = [],[],[],[],[],[],[]
        TracerXcoords, TracerYcoords, TracerExitAge, TracerIDs, TracerTimeIn = [],[],[],[],[]
        ResXcoords, ResYcoords, ResIDs, ResTypes, ResExitAge, ResTimeIn, ResVals = [],[],[],[],[],[],[]

        # Lattice Boltzmann PARAMETERS
        n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

        # Seed or do not seed the community ############################################
        if seedcom > 0:
            # inflow of resources
            ResTypes, ResVals, ResXcoords, ResYcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(ResTypes, ResVals, ResXcoords, ResYcoords, ResID, ResIDs, ResTimeIn, r, rmax, nr, width, height, u0)

            # immigration
            SpeciesIDs, IndXcoords, IndYcoords, MaintDict, GrowthDict, DispParamDict, SpColorDict, IDs, ID, TimeIn, Qs, ResUseDict = bide.immigration(m, SpeciesIDs, IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict, nr, u0, alpha)

        ####################### REPLACE ENVIRONMENT
        tracer_scatImage.remove()
        resource_scatImage.remove()
        Ind_scatImage.remove()

        fig.add_subplot(1,1,1)
        plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

        tracer_scatImage = plt.scatter([0],[0], alpha=0)
        resource_scatImage = plt.scatter([0],[0], alpha=0)
        Ind_scatImage = plt.scatter([0],[0], alpha=0)


############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/SimData.csv','w+')
OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/RADs.csv','w+')
OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/Species.csv','w+')
# printing physical variables, residence times, community diversity properties, physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, ind.prod, biomass.prod, res.inflow, res.types, max.res, barrier.width, barrier.height, logseries.a, starting.seed, stop.point, FlowRate, Width, Height, Viscosity, N, immigration.rate, particle.tau, cell.tau, resource.concentration, shannons.resource.diversity, resource.richness, S, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, avg.per.capita.growth, avg.per.capita.maint'
#             ct1,   sim, prod_i,   prod_q,       r,          nr,        rmax,    BarrierWidth,  BarrierHeight,  alpha,       seedcom,          stop,       u0,       width, height, viscosity, N, m,                TracerTau,    IndTau,ResDens,               ResDiv,                      ResRich,           S, ES,        Ev,    BP,            SD,         Nm,    sk,   Mu,                    Maint
OUT1.close()
OUT2.close()
OUT3.close()

################ DIMENSIONAL & MODEL CONSTANTS ##################################
width, height, alpha, motion, seedcom, m, r, nr, rmax = bide.get_rand_params()

#######################  Ind COMMUNITY PARAMETERS  #########################
TracerTau, IndTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = 0,0,0,0,0,0,0,0,0,0,0,0,0,0
IndID, ResID, N, ct1, T, R, prod_i, prod_q = 0,0,0,0,0,0,0,0

RAD, splist = [], []
IndTimeIn, SpeciesIDs, IndXcoords, IndYcoords, IndIDs, Qs, IndExitAge = [],[],[],[],[],[],[]
TracerXcoords, TracerYcoords, TracerExitAge, TracerIDs, TracerTimeIn = [],[],[],[],[]
ResXcoords, ResYcoords, ResIDs, ResTypes, ResExitAge, ResTimeIn, ResVals = [],[],[],[],[],[],[]

SpColorDict, GrowthDict, MaintDict = {}, {}, {}
ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
stop, shift, sign, sim, BarrierWidth, BarrierHeight = 10, 0.0, 0.1, 0, 0.1, 0.1
BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = [],[],[],[]
viscosity = 1 # unitless but required by an LBM model

Rates = np.array([1.0, 0.75, 0.5, 0.1, 0.075, 0.05, 0.025, 0.01])  # inflow speeds
u0 = Rates[0]  # initial in-flow speed

############### INITIALIZE GRAPHICS ############################################
fig = plt.figure(figsize=(12, 8))

fig.add_subplot(1,1,1) # initiate first plot

Ind_scatImage = plt.scatter([0],[0], alpha=0)
tracer_scatImage = plt.scatter([0],[0], alpha=0)
resource_scatImage = plt.scatter([0],[0], alpha=0)

#####################  Lattice Boltzmann PARAMETERS  ###########################
n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

left = BarrierXcoords1[0]
bottom = BarrierYcoords1[0]
BHeight = BarrierYcoords1[1] - bottom
BWidth = BarrierXcoords1[1] - left
BarrierImage1 = plt.bar(left-0.3, BHeight, BWidth-0.3, bottom, color = '0.3', edgecolor = '0.4', alpha=0.2)

left = BarrierXcoords2[0]
bottom = BarrierYcoords2[0]
BHeight = BarrierYcoords2[1] - bottom
BWidth = BarrierXcoords2[1] - left
BarrierImage2 = plt.bar(left-0.3, BHeight, BWidth-0.3, bottom, color = '0.3', edgecolor = '0.4', alpha=0.2)


# Seed or do not seed the community ############################################
if seedcom > 0:
    # inflow of resources
    ResTypes, ResVals, ResXcoords, ResYcoords, ResIDs, ResID, ResTimeIn = bide.ResIn(ResTypes, ResVals, ResXcoords, ResYcoords, ResID, ResIDs, ResTimeIn, r, rmax, nr, width, height, u0)
    # immigration
    SpeciesIDs, IndXcoords, IndYcoords, MaintDict, GrowthDict, DispParamDict, SpColorDict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict = bide.immigration(m, SpeciesIDs, IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, SpColorDict, IndIDs, IndID, IndTimeIn, Qs, ResUseDict, nr, u0, alpha)


Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=5000, interval=100, blit=False) # 20000 frames is a long movie
plt.show()
