from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from random import randint, choice
from scipy import stats
import numpy as np
import sys
import os
import time
import psutil

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics
sys.path.append(mydir + "/GitHub/hydrobide/tools/LBM")
import lbmMovie as LBM
sys.path.append(mydir + "/GitHub/hydrobide/tools/bide/bideMovie")
import bideMovie as bide



def get_mrrmax():
    """ Get random model parameter values. Others are chosen in bide.pyx """

    m = int(choice(range(1, 10))) # individuals immigrating per time step
    r = int(choice(range(10, 100))) # resource particles flowing in per time step
    nr = int(choice(range(1, 10))) # maximum number of resources types
    rmax = int(choice(range(100, 1000))) # maximum value of resource particle size

    return [m, r, nr, rmax]


######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    global width, height, Rates, u0, shift, sign, Ss, Ns, EVs, NSs, barrier, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, COM
    global MicXcoords, MicYcoords, microbe_scatImage, microbe_color_dict, GrowthDict, MaintDict, AvgTaus, MicIDs, MicQs, MicID, MicIDs

    global MicTimeIn, MicExitAge, avgTau, TracerIDs, TracerExitAge, TracerXcoords, TracerYcoords, tracer_scatImage, resource_scatImage
    global ResXcoords, ResYcoords, ResID, ResIDs, RES, LogSeriesAlpha, omega, MUs, VarMUs, Maints, VarMaints, ResType, ResUseDict, DispParamsDict

    global one9th, four9ths, one36th, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, sim, RAD, splist
    global BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight, ct1, m, r, nr, rmax

    for step in range(1): # adjust number of steps for smooth animation

        # new tracers
        TracerIDs, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0)

        # inflow of resources
        RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, u0)

	# immigration
        COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict = bide.immigration(m, COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr, u0, LogSeriesAlpha)

        # stream
        nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])

        # collide
        rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

        # dispersal
        args1 = [COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height]
        args2 = [u0, MicIDs, MicID, MicTimeIn, MicQs]
        args1, args2 = bide.dispersal(args1, args2)
        COM, MicXcoords, MicYcoords, MicExitAge = args1
        MicIDs, MicID, MicTimeIn, MicQs = args2

        # resource flow
        args = bide.resFlow([RES, ux, uy, u0, ResXcoords, ResYcoords, width, height, ResID, ResIDs])
        RES, ResXcoords, ResYcoords, ResID, ResIDs = args

        # moving tracer particles
        TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords = bide.MoveTracers(TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy)

        # consume and reproduce
        RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, ResType = bide.ConsumeAndReproduce(RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict, ResType, ResUseDict)

        # maintenance
        args = [COM, MicXcoords, MicYcoords, MicExitAge, microbe_color_dict, MaintDict, MicIDs, MicID, MicTimeIn, MicQs, height, width]
        args = bide.maintenance(args)
        COM, MicXcoords, MicYcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs = args

    ########## GENERATE FIGURES ############################################
    fig.add_subplot(1, 1, 1)  # Plot 1: plot of the system

    plt.tick_params(\
    axis='both',       # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    left='off',
    right='off',
    labelbottom='off',
    labelleft='off')   # labels along the bottom edge are off

    T = ['Microbes consume resources, grow, reproduce, and die as they flow through a fluid environment.',
        '\nAverage speed on the x-axis is '+str(u0)+' units per time step. Residence time is estimated from',
        'inert tracers (red stars).\nOpen circles are resource particles. Semi-impermeable barriers (grey bars) produce turbulence.']

    txt.set_text(' '.join(T))
    plt.draw()
    plt.ylim(0,height)
    plt.xlim(0,width)

    ##### PLOTTING THE MICROBES ############################################
    tracer_scatImage.remove()
    tracer_scatImage = plt.scatter(TracerXcoords, TracerYcoords, c = 'r', marker='*', lw=0.0, s = 200, alpha=0.6)

    resource_scatImage.remove()
    resource_scatImage = plt.scatter(ResXcoords, ResYcoords, c = 'w', edgecolor = 'SpringGreen', s = RES, lw = 0.6, alpha=0.7)

    microbe_scatImage.remove()
    colorlist = []

    for i, val in enumerate(COM): colorlist.append(microbe_color_dict[val])
    microbe_scatImage = plt.scatter(MicXcoords, MicYcoords, c = colorlist, edgecolor = 'k', s = MicQs, lw = 0.2, alpha=1.0)


    if len(TracerExitAge) >= 10:
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
            m, r, nr, rmax = get_mrrmax()
            sim += 1
            LogSeriesAlpha = np.random.uniform(0.9, 0.999)
            print '\n'

        Rates = np.roll(Rates, -1, axis=0)
        u0 = Rates[0]  # initial in-flow speed

        MicTimeIn, COM, MicXcoords, MicYcoords, TracerXcoords, TracerYcoords, RES, ResXcoords, ResYcoords, ResIDs, ResType, MicIDs, MicQs, MicExitAge, TracerExitAge, TracerIDs = [list([]) for _ in xrange(16)]
        # Lattice Boltzmann PARAMETERS
        n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

        ####################### REPLACE ENVIRONMENT
        fig.add_subplot(1, 1, 1)

        tracer_scatImage.remove()
        tracer_scatImage = plt.scatter([0],[0], alpha=0.0)

        resource_scatImage.remove()
        resource_scatImage = plt.scatter([0],[0], alpha=0.0)

        microbe_scatImage.remove()
        microbe_scatImage = plt.scatter([0],[0], alpha=0.0)

        plt.tick_params(\
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        left='off',
        right='off',
        labelbottom='off',
        labelleft='off') # labels along the bottom edge are off
        plt.subplots_adjust(wspace=0.35, hspace=0.35)




############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/SimData.csv','w+')
OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/RADs.csv','w+')
OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/Species.csv','w+')
# printing physical variables, residence times, community diversity properties, physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, FlowRate, Width, Height, Viscosity, N, S, immigration.rate, particle.tau, cell.tau, resource.concentration, shannons.resource.diversity, resource.richness, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, avg.per.capita.growth, avg.per.capita.maint'
#             ct1,   sim,   u0,     width, height, viscosity, N, S, m,                  TracerTau,  MicrobeTau, ResDens,                ResDiv,                    ResRich,            ES,        Ev,    BP,            SD,         Nm,    sk,         Mu,               Maint
OUT1.close()
OUT2.close()
OUT3.close()

################ DIMENSIONAL & MODEL CONSTANTS ##################################
m, r, nr, rmax = get_mrrmax()
#######################  MICROBE COMMUNITY PARAMETERS  #########################
MicID, ResID, N, S, ct1 = 0, 0, 0, 0, 0

COM, MicXcoords, MicYcoords, AvgTaus, RAD, splist = [], [], [], [], [], []
MicIDs, MicQs, MicExitAge, MicTimeIn = [], [], [], []
ResYcoords, RES, ResIDs, ResType, ResXcoords = [], [], [], [], []
TracerIDs, TracerXcoords, TracerYcoords, TracerExitAge = [], [], [], []

microbe_color_dict, GrowthDict, MaintDict = {}, {}, {}
ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
shift, sign, sim, ct1, BarrierWidth, BarrierHeight = 0.0, 0.1, 0, 0, 0.2, 0.2
Rates = np.array([1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01])  # inflow speeds
u0 = Rates[0]  # initial in-flow speed

BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = [[],[],[],[]]

width, height = 10, 10
LogSeriesAlpha = np.random.uniform(0.9, 0.999)

#####################  Lattice Boltzmann PARAMETERS  ###########################
viscosity =  0.84   # fluid viscosity of water
n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)


###############  GRAPHICS AND ANIMATION ##########################################################################################################
fig = plt.figure(figsize=(12, 8))
fig.add_subplot(1, 1, 1) # initiate first plot

microbe_scatImage = plt.scatter([0],[0], alpha=0.0)
tracer_scatImage = plt.scatter([0],[0], alpha=0.0)
resource_scatImage = plt.scatter([0],[0], alpha=0.0)

left = BarrierXcoords1[0]
bottom = BarrierYcoords1[0]
bheight = BarrierYcoords1[1] - bottom
bwidth = BarrierXcoords1[1] - left
BarrierImage1 = plt.bar(left-0.3, bheight, bwidth-0.3, bottom, color = '0.3', edgecolor = '0.4', alpha=0.2)

left = BarrierXcoords2[0]
bottom = BarrierYcoords2[0]
bheight = BarrierYcoords2[1] - bottom
bwidth = BarrierXcoords2[1] - left
BarrierImage2 = plt.bar(left-0.3, bheight, bwidth-0.3, bottom, color = '0.3', edgecolor = '0.4', alpha=0.2)

T = ['Microbes consume resource particles and grow, reproduce, and die as they flow through a complex fluid environment.',
        '\nCurrent speed on the x-axis is '+str(round(u0,3)), 'units per time step.']

txt = fig.suptitle(' '.join(T), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=5000, interval=50, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/Hydro-bide/results/movies/HydrobideVideoTest.avi', metadata={'artist':'Guido'}, bitrate=5000)
