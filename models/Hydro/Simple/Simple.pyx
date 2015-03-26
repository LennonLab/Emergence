from __future__ import division
from random import randint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os


mydir = os.path.expanduser("~/GitHub")
sys.path.append(mydir + "/HYDRO-BIDE/tools/LBM")
import LBM
sys.path.append(mydir + "/HYDRO-BIDE/tools/bide")
import bide


def Lattice():

    global n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, one9th, four9ths, one36th, rho, u0, width
    global ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, height
    global BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight

    args1 = [n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, one9th, four9ths, one36th, rho, u0, width]
    args2 = [ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, height]
    args3 = [BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight]

    args1, args2, args3 = LBM.SetLattice(args1, args2, args3)

    n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, one9th, four9ths, one36th, rho, u0, width = args1
    ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW, height = args2
    BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight = args3

    return


######### Function called for each successive animation frame ##################

def nextFrame(arg):	# arg is the frame number

    global width, height, Rates, u0, shift, sign, Ss, Ns, EVs, NSs, barrier, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, COM
    global MicXcoords, MicYcoords, microbe_scatImage, microbe_color_dict, GrowthDict, MaintDict, AvgTaus, MicIDs, MicQs, MicID, MicIDs

    global MicTimeIn, MicExitAge, avgTau, TracerIDs, TracerExitAge, TracerXcoords, TracerYcoords, tracer_scatImage, resource_scatImage
    global ResXcoords, ResYcoords, ResID, ResIDs, RES, DataColor, LogSeriesAlpha, omega, MUs, VarMUs, Maints, VarMaints

    global one9th, four9ths, one36th, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW
    global BarrierXcoords1, BarrierYcoords1, BarrierXcoords2, BarrierYcoords2, BarrierWidth, BarrierHeight

    print 'Frame:', arg, 'list number:', len(AvgTaus)+1, 'list length:', len(TracerExitAge),
    print ' u0:',u0,' N:', len(COM), ' S:', len(list(set(COM)))

    args = bide.NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0)
    TracerIDs, TracerXcoords, TracerYcoords = args

    for step in range(1): # adjust number of steps for smooth animation

	# immigration
        args1 = [COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height, MaintDict]
        args2 = [u0, GrowthDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, LogSeriesAlpha]
        args1, args2 = bide.immigration(args1, args2)
        args1 = [COM, MicXcoords, MicYcoords, MicExitAge, MaintDict] = args1
        args2 = [GrowthDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs] = args2

        # stream
        args = LBM.stream([nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign])
        nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, shift, sign = args

        # collide
        args = LBM.collide(rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, omega, u0, four9ths, one9th, one36th)
        rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, omega, u0 = args

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
        args = [TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy, COM, Ns, Ss, EVs, NSs, GrowthDict, MUs, VarMUs, Maints, VarMaints, MaintDict]
        args = bide.MoveTracers(args)
        TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy, COM, Ns, Ss, EVs, NSs, GrowthDict, MUs, VarMUs, Maints, VarMaints, MaintDict = args

        # consume and reproduce
        args = [RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict]
        args = bide.ConsumeAndReproduce(args)
        RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords = args

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
    tracer_scatImage = plt.scatter(TracerXcoords, TracerYcoords, c = 'r', marker='*', lw=0.0, s = 100, alpha=0.5)

    resource_scatImage.remove()
    resource_scatImage = plt.scatter(ResXcoords, ResYcoords, c = 'w', edgecolor = 'SpringGreen', s = RES, lw = 0.5, alpha=0.4)

    microbe_scatImage.remove()
    colorlist = []

    for i, val in enumerate(COM): colorlist.append(microbe_color_dict[val])
    microbe_scatImage = plt.scatter(MicXcoords, MicYcoords, c = colorlist, edgecolor = 'k', s = MicQs, lw = 0.2, alpha=1.0)

    if len(TracerExitAge) >= 4 and u0 <= 0.01:
        avgTau = np.mean(TracerExitAge)
        AvgTaus.append(avgTau)
        TracerExitAge = []
        MicExitAge = []

    elif len(TracerExitAge) >= 40 and u0 <= 0.1:
        avgTau = np.mean(TracerExitAge)
        AvgTaus.append(avgTau)
        TracerExitAge = []
        MicExitAge = []
    elif len(TracerExitAge) >= 80 and u0 <= 1:
        avgTau = np.mean(TracerExitAge)
        AvgTaus.append(avgTau)
        TracerExitAge = []
        MicExitAge = []
    if len(AvgTaus) >= 1:

        ######################## REPLACE ENVIRONMENT
        fig.add_subplot(1, 1, 1)

        tracer_scatImage.remove()
        tracer_scatImage = plt.scatter([0],[0], alpha=0.0)

        resource_scatImage.remove()
        resource_scatImage = plt.scatter([0],[0], alpha=0.0)

        microbe_scatImage.remove()
        microbe_scatImage = plt.scatter([0],[0], alpha=0.0)

        Ns, Ss, NSs, EVs, MUs, VarMUs, AvgTaus, MicTimeIn, COM, MicXcoords, MicYcoords, TracerXcoords, TracerYcoords = [list([]) for _ in xrange(13)]
        RES, ResXcoords, ResYcoords, ResIDs, MicIDs, MicQs, MicExitAge, TracerExitAge, TracerIDs, Maints, VarMaints = [list([]) for _ in xrange(11)]

        Rates = np.roll(Rates, -1, axis=0)
        u0 = Rates[0]  # initial in-flow speed

        LogSeriesAlpha = float(np.random.uniform(0.99, 0.9999))

        viscosity =  u0*10   	                       # fluid viscosity
        omega = 1 / (3 * viscosity + 0.5)	          # relaxation parameter

        Lattice()
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



################ DIMENSIONAL & MODEL CONSTANTS ##################################
Evar = []
EVs = list()
Ss = list()
Ns = list()
NSs = list()

MUs = list()
VarMUs = list()
Maints = list()
VarMaints = list()

AvgTaus = list()
avgTau = str()

#######################  MICROBE COMMUNITY PARAMETERS  #########################
COM = []
Xs = []

r = lambda: randint(0,255)
DataColor = '#%02X%02X%02X' % (r(),r(),r())

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

########################################## RESOURCE PARAMETERS #################
RES = []
ResIDs = []
ResID = 0

ResYcoords = []
ResXcoords = []
ResColorDict = {}

######################  SIMULATION VARIABLES  ##################################
shift = 0.0
sign = 0.1

BarrierXcoords1 = []
BarrierYcoords1 = []
BarrierXcoords2 = []
BarrierYcoords2 = []
#######################  Inert tracer particles  ###############################
TracerXcoords = []
TracerYcoords = []
TracerIDs = []
TracerExitAge = []

#####################  Lattice Boltzmann PARAMETERS  ###########################
width = 5
height = 10
LogSeriesAlpha = float(np.random.uniform(0.99, 0.9999))

BarrierWidth = 0.2
BarrierHeight = 0.2

Rates = np.array([1, 0.25, 0.1, 0.05]) # inflow speeds

u0 = Rates[0]  # initial in-flow speed
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

rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW   # macroscopic density
ux = (nE + nNE + nSE - nW - nNW - nSW) / rho	# macroscopic x velocity
uy = (nN + nNE + nNW - nS - nSE - nSW) / rho	# macroscopic y velocity

barrier = np.zeros((height, width), bool)  # Initialize barriers

barrierN = np.roll(barrier,  1, axis=0)       # sites just north of barriers
barrierS = np.roll(barrier, -1, axis=0)       # sites just south of barriers
barrierE = np.roll(barrier,  1, axis=1)       # etc.
barrierW = np.roll(barrier, -1, axis=1)
barrierNE = np.roll(barrierN,  1, axis=1)
barrierNW = np.roll(barrierN, -1, axis=1)
barrierSE = np.roll(barrierS,  1, axis=1)
barrierSW = np.roll(barrierS, -1, axis=1)

###############  GRAPHICS AND ANIMATION ##########################################################################################################
fig = plt.figure(figsize=(12, 8))
fig.add_subplot(1, 1, 1) # initiate first plot

Lattice()

microbe_scatImage = plt.scatter([0],[0], alpha=0.0)
tracer_scatImage = plt.scatter([0],[0], alpha=0.0)
resource_scatImage = plt.scatter([0],[0], alpha=0.0)

left = BarrierXcoords1[0]
bottom = BarrierYcoords1[0]
bheight = BarrierYcoords1[1] - bottom
bwidth = BarrierXcoords1[1] - left
BarrierImage1 = plt.bar(left-0.3, bheight, bwidth-0.3, bottom, color = '0.3', edgecolor = '0.3', alpha=0.2)

left = BarrierXcoords2[0]
bottom = BarrierYcoords2[0]
bheight = BarrierYcoords2[1] - bottom
bwidth = BarrierXcoords2[1] - left
BarrierImage2 = plt.bar(left-0.3, bheight, bwidth-0.3, bottom, color = '0.3', edgecolor = '0.3', alpha=0.2)

T = ['Microbes consume resource particles and grow, reproduce, and die as they flow through a complex fluid environment.',
        '\nCurrent speed on the x-axis is '+str(round(u0,3)), 'units per time step.']

txt = fig.suptitle(' '.join(T), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=5000, interval=50, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/Hydro-bide/results/movies/HydrobideVideoTest.avi', metadata={'artist':'Guido'}, bitrate=5000)
