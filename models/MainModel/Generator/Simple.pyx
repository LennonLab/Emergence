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
sys.path.append(mydir + "/GitHub/hydrobide/tools/LBM")
import LBM
sys.path.append(mydir + "/GitHub/hydrobide/tools/bide")
import bide


def checklists(Qs2, Slist2, inds2, Tlist, xcoords, ycoords, zcoords):

    Blist = list(set([len(Qs2), len(Slist2), len(inds2), len(Tlist), len(xcoords), len(ycoords), len(zcoords)]))
    if len(Blist) > 1:
        print 'lists of different sizes:', Blist
        sys.exit()

    N = len(Qs2)
    if N == 0:
        return ['N = 0']

    return ['pass']


def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.pyx """


    motion = choice(['fluid', 'conveyor', 'random_walk', 'uncorrelated'])
    D = str()
    if motion == 'uncorrelated' or motion == 'random_walk':
        D = choice([2, 3])
    else:
        D = 2

    reproduction = choice(['clonal', 'sexual'])
    mutation = choice(['yes', 'no'])
    predators = choice(['yes', 'no'])
    parasites = choice(['yes', 'no'])
    symbionts = choice(['yes', 'no'])
    env_gradient = choice(['no', 'yes'])

    # richness of the metacommunity
    J = choice([100, 1000, 10000])

    # size of starting community
    seed = choice([0, 10, 100, 1000])

    # individuals immigrating per time step
    m = choice([0, 2, 4, 8])

    # resource particles flowing in per time step
    r = choice([0, 10, 50, 100])

    # maximum number of resources types
    nr = choice([1, 2, 4, 8, 16, 32])

    # maximum resource particle size
    rmax = choice([500, 1000, 2000, 4000, 8000])

    # mean and standard deviation for number of prey
    avg_prey = [np.random.uniform(0, 10), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for number of symbionts
    avg_symb = [np.random.uniform(0, 10), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for number of parasites
    avg_parasite = [np.random.uniform(0, 10), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for specific growth rate
    avg_growth = [np.random.uniform(0.1, 1.0), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for propagule cell quota
    avg_Q = [np.random.uniform(0.1, 1.0), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for specific maintenance
    avg_maint = [np.random.uniform(0.01, 0.1), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for specific active dispersal
    avg_disp = [np.random.uniform(0.01, 1.0), np.random.uniform(0.01, 0.1)]

    # mean and standard deviation for specific resource use efficiency
    avg_res = [np.random.uniform(0.01, 1.0), np.random.uniform(0.01, 0.1)]

    return [motion, D, reproduction, mutation, predators, parasites, symbionts, env_gradient, J, seed, m, r, nr, rmax, avg_prey, avg_symb, avg_parasite, avg_growth, avg_Q, avg_maint, avg_disp, avg_res]


######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    global width, height, Rates, u0, shift, sign, Ss, Ns, EVs, NSs, barrier, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, COM
    global IndXcoords, IndYcoords, Ind_scatImage, Sp_color_dict, GrowthDict, MaintDict, AvgTaus, IndIDs, IndQs, IndID, IndIDs

    global IndTimeIn, IndExitAge, avgTau, TracerIDs, TracerExitAge, TracerXcoords, TracerYcoords, tracer_scatImage, resource_scatImage
    global ResXcoords, ResYcoords, ResID, ResIDs, RES, alpha, omega, MUs, VarMUs, Maints, VarMaints, ResType, ResUseDict, DispParamsDict

    global one9th, four9ths, one36th, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, sim, RAD, splist
    global BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight, ct1, Mu, Maint

    global motion, D, reproduction, mutation, predators, parasites, symbionts, env_gradient, J, seed, m, r, nr, rmax, avg_prey, avg_symb, avg_parasite, avg_growth, avg_Q, avg_maint, avg_disp, avg_res
    global N, TracerTau, IndTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, T, R, seed, stop, prod_i, prod_q

    for step in range(1): # adjust number of steps for smooth animation

        random.seed() # use current system time to initiate a random seed (ensures that autocorrelation doesn't crop-up, as is possible when using pseudorandom number generators)

        # new tracers
        TracerIDs, TracerXcoords, TracerYcoords = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0, D)

        # inflow of resources
        RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, u0, D)

	# immigration
       SpeciesIDs,IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamDict, Sp_color_dict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict = bide.immigration(m, SpeciesIDs, IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, Sp_color_dict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict, nr, u0, alpha, D)

        if motion == 'fluid' or motion == 'conveyor':  # a 'conveyor' belt action wherein y-coordinates never change will occur when there is no turbulence in a fluid dynamics model, most analogous to an infinitely viscous fluid

            # stream
            nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])

            # collide
            rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW = LBM.collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0)

            # dispersal
           SpeciesIDs,IndXcoords, IndYcoords, IndExitAge, IndIDs, IndID, IndTimeIn, IndQs = bide.fluid_movement(SpeciesIDs, IndIDs, IndID, IndTimeIn, IndQs, IndExitAge, ux, uy, IndXcoords, IndYcoords, IndZcoords, width, height, u0)

            # resource flow
            ResType, ResXcoords, ResYcoords, ResExitAge, ResIDs, ResID, ResTimeIn, ResVals = bide.fluid_movement(ResType, ResIDs, ResID, ResTimeIn, ResVals, ResExitAge, ux, uy, ResXcoords, ResYcoords, ResZcoords, width, height, u0)

            # moving tracer particles
            TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, T, R, ResDens, ResDiv, ResRich, TracerTau, IndTau, N, S, Mu, Maint, Ev, ES, Nm , BP , SD , sk = bide.fluid_MoveTracers(TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, T, R, RES, ResType, ResDens, ResDiv, ResRich, TracerTau,IndTau, IndExitAge,SpeciesIDs,N, S, Mu, Maint, GrowthDict, MaintDict,Ev, ES, Nm , BP , SD , sk, D)

        elif motion == 'random_walk' or motion == 'uncorrelated':

            # Moving tracer particles
            Lists = [TracerIDs]
            coords = [TracerXcoords, TracerYcoords]
            if D == 3:
                coords.append(TracerZcoords)

            nonfluid_movement('tracer', Lists, TracerExitAge, TracerTimeIn, coords, width, height, length, u0, D)

            # Moving resource particles
            coords = [ResXcoords, ResYcoords]
            if D == 3:
                coords.append(ResZcoords)

            Lists = [ResType, ResIDs, ResVals]
            nonfluid_movement('resource', Lists, ResExitAge, ResTimeIn, coords, width, height, length, u0, D)


            # Moving individuals
            coords = [IndXcoords, IndYcoords]
            if D == 3:
                coords.append(IndZcoords)

            Lists = [SpeciesIDs, IndIDs, Qs, DispParamDict]
            nonfluid_movement('individual', Lists, IndExitAge, IndTimeIn, coords, width, height, length, u0, D)


        # consume and reproduce
        p1 = len(COM)
        q1 = sum(IndQs)

        ResLists, IndLists = bide.ConsumeAndReproduce(ResType, ResVals, ResIDs, ResID, ResCoords, ResTimeIn, ResExitAge, IndType, IndQs, IndIDs, IndID, IndTimeIn, IndCoords, width, height, length, GrowthDict, ResUseDict)
        ResType, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords = ResLists
        IndType, IndQs,   IndIDs, IndID, IndTimeIn,             IndXcoords, IndYcoords, IndZcoords = IndLists

        prod_i = len(COM) - p1
        prod_q = sum(IndQs) - q1

        # maintenance
       SpeciesIDs,IndXcoords, IndYcoords, IndExitAge, IndIDs, IndID, IndTimeIn, IndQs = bide.maintenance(SpeciesIDs, IndXcoords, IndYcoords, IndExitAge, Sp_color_dict, MaintDict, IndIDs, IndID, IndTimeIn, IndQs, height, width, D)


    ########## GENERATE FIGURES ############################################
    if D = 2 or motion == 'fluid' or motion == 'conveyor':
        fig.add_subplot(1, 1, 1)  # Plot 1: plot of the system
        plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    elif D = 3:
        fig.add_subplot(111, projection='3d')
        plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    N = len(COM)
    RAD, splist = bide.GetRAD(COM)
    S = len(RAD)

    Title = ['Inds consume resources, grow, reproduce, and die as they flow through a fluid environment. Average speed',
           '\non the x-axis is '+str(u0)+' units per time step. '+str(len(TracerExitAge))+' tracers have passed through.',
           'N = '+str(N)+', S = '+str(S)+'.'
           '\nOpen circles are resource particles. Semi-impermeable barriers (grey bars) produce turbulence.']

    txt.set_text(' '.join(Title))
    plt.draw()
    plt.ylim(0,height)
    plt.xlim(0,width)

    ##### PLOTTING THE IndS ############################################
    resource_scatImage.remove()
    resource_scatImage = plt.scatter(ResXcoords, ResYcoords, c = 'w', edgecolor = 'SpringGreen', s = RES, lw = 0.6, alpha=0.7)

    tracer_scatImage.remove()
    tracer_scatImage = plt.scatter(TracerXcoords, TracerYcoords, c = 'r', marker='*', lw=0.0, s = 200, alpha=0.6)

    Ind_scatImage.remove()
    colorlist = []

    for i, val in enumerate(COM): colorlist.append(Sp_color_dict[val])
    Ind_scatImage = plt.scatter(IndXcoords, IndYcoords, c = colorlist, edgecolor = '0.2', s = IndQs, lw = 0.2, alpha=0.9)

    # Record model values and reset, or not
    if len(TracerExitAge) >= stop:

        # Examining the resource RAD
        if len(ResType) > 0:
            ResRAD, Rlist = bide.GetRAD(ResType)
            ResDens = sum(RES)/(height*width)
            ResDiv = float(metrics.Shannons_H(ResRAD))
            ResRich = len(Rlist)

        # Residence times for tracers and Inds
        TracerTau = float(np.mean(TracerExitAge))
        IndTau = float(np.mean(IndExitAge))

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
        print>>OUT1, ct1,',', sim,',', prod_i,',', prod_q,',', r,',', nr,',', rmax,',', BarrierWidth,',', BarrierHeight,',', alpha,',', seed,',', stop,',', u0,',', width,',', height,',', viscosity,',', N,',', m,',', TracerTau,',', IndTau,',', ResDens,',', ResDiv,',', ResRich,',', S,',', ES,',', Ev,',', BP,',', SD,',', Nm,',', sk,',', Mu,',', Maint
        print>>OUT2, RADString
        print>>OUT3, SString
        OUT1.close()
        OUT2.close()
        OUT3.close()


        if u0 == min(Rates):
            Sp_color_dict, GrowthDict, MaintDict = {}, {}, {}
            ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}
            motion, D, reproduction, mutation, predators, parasites, symbionts, env_gradient, J, seed, m, r, nr, rmax, avg_prey, avg_symb, avg_parasite, avg_growth, avg_Q, avg_maint, avg_disp, avg_res = get_rand_params()
            sim += 1
            alpha = np.random.uniform(0.9, 0.999)
            print '\n'

        Rates = np.roll(Rates, -1, axis=0)
        u0 = Rates[0]  # initial in-flow speed

        TracerTau, IndTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        IndTimeIn,SpeciesIDs,IndXcoords, IndYcoords, TracerXcoords, TracerYcoords, RES, ResXcoords, ResYcoords, ResIDs, ResType, IndIDs, IndQs, IndExitAge, TracerExitAge, TracerIDs = [list([]) for _ in xrange(16)]

        if motion == 'fluid' or motion == 'conveyor':
            #####################  Lattice Boltzmann PARAMETERS  ###########################
            viscosity =  1 # unitless but required by an LBM model
            n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)

        elif motion == 'random_walk' or motion == 'uncorrelated':


        # Seed or do not seed the community ############################################
        if seed > 0:
            # inflow of resources
            RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, 1)
            # immigration
           SpeciesIDs,IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamDict, Sp_color_dict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict = bide.immigration(seed,SpeciesIDs,IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, Sp_color_dict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict, nr, 1, alpha)

        ####################### REPLACE ENVIRONMENT
        fig.add_subplot(1, 1, 1)

        tracer_scatImage.remove()
        tracer_scatImage = plt.scatter([0],[0], alpha=0.0)

        resource_scatImage.remove()
        resource_scatImage = plt.scatter([0],[0], alpha=0.0)

        Ind_scatImage.remove()
        Ind_scatImage = plt.scatter([0],[0], alpha=0.0)

        ########## GENERATE FIGURES ############################################
        if D = 2 or motion == 'fluid' or motion == 'conveyor':
            fig.add_subplot(1, 1, 1)  # Plot 1: plot of the system
            plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

        elif D = 3:
            fig.add_subplot(111, projection='3d')
            plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')



############## OPEN OUTPUT DATA FILE ###########################################
OUT1 = open(mydir + '/GitHub/hydrobide/results/simulated_data/SimData.csv','w+')
OUT2 = open(mydir + '/GitHub/hydrobide/results/simulated_data/RADs.csv','w+')
OUT3 = open(mydir + '/GitHub/hydrobide/results/simulated_data/Species.csv','w+')
# printing physical variables, residence times, community diversity properties, physiological values, trait values, resource values
print>>OUT1, 'RowID, sim, ind.prod, biomass.prod, res.inflow, res.types, max.res, barrier.width, barrier.height, logseries.a, starting.seed, stop.point, FlowRate, Width, Height, Viscosity, N, immigration.rate, particle.tau, cell.tau, resource.concentration, shannons.resource.diversity, resource.richness, S, simpson.e, e.var, berger.parker, inv.simp.D, N.max, skew, avg.per.capita.growth, avg.per.capita.maint'
#             ct1,   sim, prod_i,   prod_q,       r,          nr,        rmax,    BarrierWidth,  BarrierHeight,  alpha,       seed,          stop,       u0,       width, height, viscosity, N, m,                TracerTau,    IndTau,ResDens,               ResDiv,                      ResRich,           S, ES,        Ev,    BP,            SD,         Nm,    sk,   Mu,                    Maint
OUT1.close()
OUT2.close()
OUT3.close()

################ DIMENSIONAL & MODEL CONSTANTS ##################################
motion, D, reproduction, mutation, predators, parasites, symbionts, env_gradient, J, seed, m, r, nr, rmax, avg_prey, avg_symb, avg_parasite, avg_growth, avg_Q, avg_maint, avg_disp, avg_res = get_rand_params()
#######################  Ind COMMUNITY PARAMETERS  #########################
IndID, ResID, N, S, ct1, Mu, Maint, T, R = 0, 0, 0, 0, 0, 0, 0, 0, 0
prod_i, prod_q = 0, 0
N, TracerTau, IndTau, ResDens, ResDiv, ResRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

SpeciesIDs, IndXcoords, IndYcoords, AvgTaus, RAD, splist = [], [], [], [], [], []
IndIDs, IndQs, IndExitAge, IndTimeIn = [], [], [], []
ResYcoords, RES, ResIDs, ResType, ResXcoords = [], [], [], [], []
TracerIDs, TracerXcoords, TracerYcoords, TracerExitAge = [], [], [], []

Sp_color_dict, GrowthDict, MaintDict = {}, {}, {}
ResUseDict, ResColorDict, DispParamsDict = {}, {}, {}

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
stop, shift, sign, sim, BarrierWidth, BarrierHeight = 10, 0.0, 0.1, 0, 0.0, 0.0

#Rates = np.array([1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01])  # inflow speeds
Rates = np.array([1.0, 0.75, 0.5, 0.1, 0.075, 0.05, 0.025, 0.01])  # inflow speeds

u0 = Rates[0]  # initial in-flow speed

BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = [[],[],[],[]]

width, height = 10, 10
alpha = 0.99


if motion == 'fluid' or motion == 'conveyor':
    #####################  Lattice Boltzmann PARAMETERS  ###########################
    viscosity =  1 # unitless but required by an LBM model
    n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2 = LBM.SetLattice(u0, viscosity, width, height, BarrierWidth, BarrierHeight, BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2)


###############  GRAPHICS AND ANIMATION ########################################
fig = plt.figure(figsize=(12, 8))
fig.add_subplot(1, 1, 1) # initiate first plot

Ind_scatImage = plt.scatter([0],[0], alpha=0.0)
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

Title = ['Inds consume resource particles and grow, reproduce, and die as they flow through a complex fluid environment.',
        '\nCurrent speed on the x-axis is '+str(round(u0,3)), 'units per time step.']

txt = fig.suptitle(' '.join(Title), fontsize = 12)


# Seed or do not seed the community ############################################
if seed > 0:
    # inflow of resources
    RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType = bide.ResIn(RES, ResXcoords, ResYcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, 1)
    # immigration
   SpeciesIDs, IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamDict, Sp_color_dict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict = bide.immigration(seed, SpeciesIDs, IndXcoords, IndYcoords, width, height, MaintDict, GrowthDict, DispParamsDict, Sp_color_dict, IndIDs, IndID, IndTimeIn, IndQs, ResUseDict, nr, 1, alpha)


ani = animation.FuncAnimation(fig, nextFrame, frames=5000, interval=100, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/Hydro-bide/results/movies/HydrobideVideoTest.avi', metadata={'artist':'Guido'}, bitrate=5000)
