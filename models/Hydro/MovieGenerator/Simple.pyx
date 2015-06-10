from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from random import randint, choice
from scipy import stats
import numpy as np
import sys
import os
import psutil

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics
sys.path.append(mydir + "/GitHub/hydrobide/tools/LBM/lbmMovie")
import lbmMovie as LBM
sys.path.append(mydir + "/GitHub/hydrobide/tools/bide/bideMovie")
import bide



def get_mrrmax():
    """ Get random model parameter values. Others are chosen in bide.pyx """

    seed = choice([1000])
    m = choice([0]) # individuals immigrating per time step
    r = choice([200]) # resource particles flowing in per time step
    nr = choice([10]) # maximum number of resources types
    rmax = choice([5000]) # maximum value of resource particle size

    return [seed, m, r, nr, rmax]


######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    global width, height, Rates, u0, shift, sign, Ss, Ns, EVs, NSs, barrier, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, COM
    global MicXcoords, MicYcoords, microbe_scatImage, microbe_color_dict, GrowthDict, MaintDict, AvgTaus, MicIDs, MicQs, MicID, MicIDs

    global MicTimeIn, MicExitAge, avgTau, TracerIDs, TracerExitAge, TracerXcoords, TracerYcoords, tracer_scatImage, resource_scatImage
    global ResXcoords, ResYcoords, ResID, ResIDs, RES, alpha, omega, MUs, VarMUs, Maints, VarMaints, ResType, ResUseDict, DispParamsDict

    global one9th, four9ths, one36th, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, sim, RAD, splist
    global BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight, ct1, m, r, nr, rmax, Mu, Maint

    global N, TracerTau, MicrobeTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, T, R, seed, stop, prod_i, prod_q

    for step in range(1): # adjust number of steps for smooth animation

        # new tracers
        TracerIDs, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0)

        # inflow of resources
        RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, u0)

	# immigration
        COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict = bide.immigration(m, COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr, u0, alpha)

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
        TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, T, R, ResDens, ResDiv, ResRich, TracerTau, MicrobeTau, N, S, Mu, Maint, Ev, ES, Nm , BP , SD , sk = bide.MoveTracers(TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy, T, R, RES, ResType, ResDens, ResDiv, ResRich, TracerTau,MicrobeTau, MicExitAge, COM, N, S, Mu, Maint, GrowthDict, MaintDict,Ev, ES, Nm , BP , SD , sk)

        # consume and reproduce
        p1 = len(COM)
        q1 = sum(MicQs)
        RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, ResType = bide.ConsumeAndReproduce(RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict, ResType, ResUseDict)
        prod_i = len(COM) - p1
        prod_q = sum(MicQs) - q1

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

    N = len(COM)
    RAD, splist = bide.GetRAD(COM)
    S = len(RAD)

    Title = ['Microbes consume resources, grow, reproduce, and die as they flow through a fluid environment. Average speed',
           '\non the x-axis is '+str(u0)+' units per time step. '+str(len(TracerExitAge))+' tracers have passed through.',
           'N = '+str(N)+', S = '+str(S)+'.'
           '\nOpen circles are resource particles. Semi-impermeable barriers (grey bars) produce turbulence.']

    txt.set_text(' '.join(Title))
    plt.draw()
    plt.ylim(0,height)
    plt.xlim(0,width)

    ##### PLOTTING THE MICROBES ############################################
    resource_scatImage.remove()
    resource_scatImage = plt.scatter(ResXcoords, ResYcoords, c = 'w', edgecolor = 'SpringGreen', s = RES, lw = 0.6, alpha=0.7)

    tracer_scatImage.remove()
    tracer_scatImage = plt.scatter(TracerXcoords, TracerYcoords, c = 'r', marker='*', lw=0.0, s = 200, alpha=0.6)

    microbe_scatImage.remove()
    colorlist = []

    for i, val in enumerate(COM): colorlist.append(microbe_color_dict[val])
    microbe_scatImage = plt.scatter(MicXcoords, MicYcoords, c = colorlist, edgecolor = '0.2', s = MicQs, lw = 0.2, alpha=0.9)

    # Record model values and reset, or not
    if len(TracerExitAge) >= stop:

        # Examining the resource RAD
        if len(ResType) > 0:
            ResRAD, Rlist = bide.GetRAD(ResType)
            ResDens = sum(RES)/(height*width)
            ResDiv = float(metrics.Shannons_H(ResRAD))
            ResRich = len(Rlist)

        # Residence times for tracers and microbes
        TracerTau = float(np.mean(TracerExitAge))
        MicrobeTau = float(np.mean(MicExitAge))

        T = len(TracerIDs)
        R = len(RES)

        N = len(COM)
        if N >= 1:

            RAD, splist = bide.GetRAD(COM)
            RAD, splist = zip(*sorted(zip(RAD, splist), reverse=True))
            S = len(RAD)

            # Specific Growth rate and Maintenance

            mu, maint = 0, 0
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
        print>>OUT1, ct1, sim, prod_i, prod_q, r, nr, rmax, BarrierWidth, BarrierHeight, alpha, seed, stop, u0, width, height, viscosity, N, m, TracerTau, MicrobeTau,ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint
        print>>OUT2, RADString
        print>>OUT3, SString
        OUT1.close()
        OUT2.close()
        OUT3.close()


        if u0 == min(Rates):
            microbe_color_dict, GrowthDict, MaintDict = {}, {}, {}
            ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}
            seed, m, r, nr, rmax = get_mrrmax()
            sim += 1
            alpha = np.random.uniform(0.9, 0.999)
            print '\n'

        Rates = np.roll(Rates, -1, axis=0)
        u0 = Rates[0]  # initial in-flow speed

        TracerTau, MicrobeTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        MicTimeIn, COM, MicXcoords, MicYcoords, TracerXcoords, TracerYcoords, RES, ResXcoords, ResYcoords, ResIDs, ResType, MicIDs, MicQs, MicExitAge, TracerExitAge, TracerIDs = [list([]) for _ in xrange(16)]
        # Lattice Boltzmann PARAMETERS
        n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)


        # Seed or do not seed the community ############################################
        if seed > 0:
            # inflow of resources
            RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, 1)
            # immigration
            COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict = bide.immigration(seed, COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr, 1, alpha)

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
print>>OUT1, 'RowID, sim, ind.prod, biomass.prod, res.inflow, res.types, max.res, barrier.width, barrier.height, logseries.a, starting.seed, stop.point, FlowRate, Width, Height, Viscosity, N, immigration.rate, particle.tau, cell.tau, resource.concentration, shannons.resource.diversity, resource.richness, S, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, avg.per.capita.growth, avg.per.capita.maint'
#             ct1,   sim, prod_i,   prod_q,       r,          nr,        rmax,    BarrierWidth,  BarrierHeight,  alpha,       seed,          stop,       u0,       width, height, viscosity, N, m,                TracerTau,    MicrobeTau,ResDens,               ResDiv,                      ResRich,           S, ES,        Ev,    BP,            SD,         Nm,    sk,   Mu,                    Maint
OUT1.close()
OUT2.close()
OUT3.close()

################ DIMENSIONAL & MODEL CONSTANTS ##################################
seed, m, r, nr, rmax = get_mrrmax()
#######################  MICROBE COMMUNITY PARAMETERS  #########################
MicID, ResID, N, S, ct1, Mu, Maint, T, R = 0, 0, 0, 0, 0, 0, 0, 0, 0
prod_i, prod_q = 0, 0
N, TracerTau, MicrobeTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

COM, MicXcoords, MicYcoords, AvgTaus, RAD, splist = [], [], [], [], [], []
MicIDs, MicQs, MicExitAge, MicTimeIn = [], [], [], []
ResYcoords, RES, ResIDs, ResType, ResXcoords = [], [], [], [], []
TracerIDs, TracerXcoords, TracerYcoords, TracerExitAge = [], [], [], []

microbe_color_dict, GrowthDict, MaintDict = {}, {}, {}
ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
stop, shift, sign, sim, BarrierWidth, BarrierHeight = 10, 0.0, 0.1, 0, 0.0, 0.0

#Rates = np.array([1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01])  # inflow speeds
Rates = np.array([1.0, 0.75, 0.5, 0.1, 0.075, 0.05, 0.025, 0.01])  # inflow speeds

u0 = Rates[0]  # initial in-flow speed

BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = [[],[],[],[]]

width, height = 10, 10
alpha = 0.99

#####################  Lattice Boltzmann PARAMETERS  ###########################
viscosity =  1
n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)


###############  GRAPHICS AND ANIMATION ########################################
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

Title = ['Microbes consume resource particles and grow, reproduce, and die as they flow through a complex fluid environment.',
        '\nCurrent speed on the x-axis is '+str(round(u0,3)), 'units per time step.']

txt = fig.suptitle(' '.join(Title), fontsize = 12)


# Seed or do not seed the community ############################################
if seed > 0:
    # inflow of resources
    RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, 1)
    # immigration
    COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict = bide.immigration(seed, COM, MicXcoords, MicYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr, 1, alpha)


ani = animation.FuncAnimation(fig, nextFrame, frames=5000, interval=100, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/Hydro-bide/results/movies/HydrobideVideoTest.avi', metadata={'artist':'Guido'}, bitrate=5000)
