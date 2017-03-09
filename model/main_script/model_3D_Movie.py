from __future__ import division

import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from random import choice, shuffle #, randint
import numpy as np
from numpy import sin, pi, mean
import sys
import os
import time

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/residence-time/model/metrics")
import metrics
sys.path.append(mydir + "GitHub/residence-time/model/bide")
import bide2 as bide
sys.path.append(mydir + "GitHub/residence-time/model/randparams")
import randparams as rp
sys.path.append(mydir + "GitHub/residence-time/model/spatial")

""" To generate movies:
1.) uncomment line 'ani.save' on or near line 364
2.) adjust the frames variable on or near line 364, to change movie length
3.) change plot_system = 'no' to 'yes' on or near line 66

Because generating animations requires computing time and memory, doing so can
be computationally demanding. To quicken the process, use plot_system = 'no' on
or near line 66.
"""

# https://www.quantstart.com/articles/Basics-of-Statistical-Mean-Reversion-Testing
# http://statsmodels.sourceforge.net/0.5.0/generated/statsmodels.tsa.stattools.adfuller.html


GenPath = mydir + 'GitHub/residence-time/results/simulated_data/'


OUT = open(GenPath + 'SimData.csv','w+')
print>>OUT, 'sim,ct,dormancy,immigration,res.inflow,N.types,max.res.val,max.growth.rate,max.met.maint,max.active.dispersal,starting.seed,flow.rate,height,length,width,viscosity,total.abundance,ind.production,biomass.prod.N,resource.tau,particle.tau,individual.tau,resource.particles,resource.concentration,species.richness,simpson.e,N.max,tracer.particles,Whittakers.turnover,avg.per.capita.growth,avg.per.capita.maint,avg.per.capita.N.efficiency,avg.per.capita.active.dispersal,amplitude,flux,frequency,phase,spec.growth,spec.disp,spec.maint,N.active,N.dormant,Percent.Dormant,avg.per.capita.RPF,avg.per.capita.MF'
OUT.close()

OUT = open(GenPath + 'RAD-Data.csv', 'w+')
OUT.close()

#######################  COMMUNITY PARAMETERS  #########################

# Lists
CRList, Ns, SpColorList, RColorList, RAD, splist, splist2, TIDs, TX, TY, RTypes, RX, RY, RZ, RIDs, RVals, SpeciesIDs, indX, indY, indZ, IndIDs, Qs, N_RD, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList = [list([]) for _ in xrange(30)]
# Scalars
nNi, u1, numA, numD, RDENS, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, IndID, RID, N, ct, T, R, PRODI, PRODN = [0]*24
# Dictionaries
SpColorDict, GrowthDict, MaintDict, MainFactorDict, RPFDict, N_RD, RColorDict, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}


################ MODEL INPUTS ##################################
width, height, length, seedCom, m, r, nNi, rmax, gmax, maintmax, dmax, amp, freq, flux, pulse, phase, envgrads, Rates, pmax, mmax, dorm, imm = rp.get_rand_params()

###############  SIMULATION VARIABLES, DIMENSIONAL & MODEL CONSTANTS  ##########
sim = 1
viscosity = 10 # viscosity is unitless but required by LBM model
u0 = Rates[0]  # initial in-flow speed


processes = range(1, 10)
t = time.clock()
BurnIn = 'not done'
p = 0.0

def nextFrame(arg):

    """ Function called for each successive animation frame; arg is the frame number """

    global mmax, pmax, ADList, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth
    global p, BurnIn, t, num_sims, width, height, Rates, u0
    global SpColorDict, GrowthDict, N_RD, CRList, RZ, length
    global DispDict, MaintDict, one9th, four9ths, one36th, barrier
    global gmax, dmax, maintmax, IndIDs, Qs, IndID, indX, indZ
    global indY,  Ind_scatImage, SpeciesIDs, EnvD, TY
    global TIDs, TTX, RTypes, RX, RY, RID, RIDs, RVals
    global resource_scatImage, ct1, Mu, Maint
    global reproduction, speciation, seedCom, m, r, nNi, rmax, sim
    global RAD, splist, N, ct, splist2, WT, RDens, RDiv, RRich, S, ES
    global Ev, BP, SD, Nm, sk, T, R, prod_i, prod_q, viscosity, dorm, imm
    global Ts, Rs, PRODIs, Ns, TTAUs, INDTAUs, RDENs, RDIVs, RRICHs, Ss, ESs, EVs
    global BPs, SDs, NMAXs, SKs, MUs, MAINTs, PRODNs, PRODPs, PRODCs, lefts, bottoms
    global Gs, Ms, NRs, PRs, CRs, Ds, RTAUs, GrowthList, MaintList, N_RList, MFDList, RPFList
    global DispList, amp, freq, flux, pulse, phase, envgrads
    global MainFactorDict, RPFDict, SpecRPF, SpecMF,t

    ct += 1
    plot_system = 'no'

    # fluctuate flow according to amplitude, frequency, & phase
    u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))
    if u1 > 1.0: u1 = 1.0

    shuffle(processes)
    for num in processes:


        if num == 1: # Inflow of resources
            RTypes, RVals, RX, RY, RZ, RIDs, RID = bide.ResIn(RTypes, RVals, RX, RY, RZ, RID, RIDs, r, rmax, nNi, height, length, width, u1)

        elif num == 2: # Resource flow
            Lists = [RTypes, RIDs, RID, RVals]
            RTypes, RX, RY, RZ, RIDs, RID, RVals = bide.flow('resource', Lists, RX, RY, RZ, height, length, width, u0)

        elif num == 3: # Inflow of individuals (immigration)
            CRList, mmax, pmax, dmax, gmax, maintmax, seedCom, SpeciesIDs, indX, indY, indZ, height, length, width, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, Qs, N_RD, nNi, u1, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList, ct, m = bide.immigration(CRList, mmax, pmax, dmax, gmax, maintmax, seedCom, SpeciesIDs, indX, indY, indZ, height, length, width, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, GrowthDict, DispDict, SpColorDict, IndIDs, IndID, Qs, N_RD, nNi, u1, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList, ct, m)

        elif num == 4: # Dispersal
            Lists = [CRList, SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList]
            CRList, SpeciesIDs, indX, indY, indZ, IndIDs, IndID, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList = bide.flow('individual', Lists, indX, indY, indZ, height, length, width, u0)

        elif num == 5: # Forage
            CRList, RVals, RIDs, RX, RY, RZ, Sp_IDs, indX, indY, indZ, Qs, IDs, ID, height, length, width, GD, DispD, colorD, N_RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, nN, GList, MList, MFDList, RPDList, NList, DList, ADList = bide.nearest_forage(CRList, RVals, RIDs, RX, RY, RZ, SpeciesIDs, indX, indY, indZ, Qs, IndIDs, IndID, height, length, width, GrowthDict, DispDict, SpColorDict, N_RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, nNi, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList)

        elif num == 6: # Consume
            numc, CRList, RPFDict, RTypes, RVals, RIDs, RID, RX, RY, RZ, SpeciesIDs, Qs, IndIDs, IndID, indX, indY, indZ, height, length, width, GrowthDict, N_RD, DispDict, GrowthList, MaintList, MFDList, RPFList, MainFactorDict, N_RList, DispList, ADList = bide.consume(CRList, RPFDict, RTypes, RVals, RIDs, RID, RX, RY, RZ, SpeciesIDs, Qs, IndIDs, IndID, indX, indY, indZ, height, length, width, GrowthDict, N_RD, DispDict, GrowthList, MaintList, MFDList, RPFList, MainFactorDict, N_RList, DispList, ADList)

        elif num == 7: # Transition to or from dormancy
            CRList, SpeciesIDs, indX, indY, indZ, IndIDs, Qs, DispList, GrowthList, MaintList, MFDList, RPFList, N_RList, MainFactorDict, RPFDict, ADList = bide.transition(CRList, SpeciesIDs, indX, indY, indZ, IndIDs, Qs, DispList, GrowthList, MaintList, MFDList, RPFList, N_RList, MainFactorDict, RPFDict, ADList)

        elif num == 8: # Maintenance
            CRList, SpeciesIDs, indX, indY, indZ, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList = bide.maintenance(CRList, SpeciesIDs, indX, indY, indZ, SpColorDict, MaintDict, MainFactorDict, RPFDict, EnvD, IndIDs, Qs, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList)

        elif num == 9: # Reproduction
            p1, TNQ1 = metrics.getprod(Qs)

            CRList, SpeciesIDs, indX, indY, indZ, Qs, IndIDs, IndID, height, length, width, GrowthDict, DispDict, SpcolorDict, N_RD, MD, MFD, RPD, EnvD, envGs, nNi, GList, MList, MFDList, RPDList, NList, DList, ADList = bide.reproduce(CRList, SpeciesIDs, indX, indY, indZ, Qs, IndIDs, IndID, height, length, width, GrowthDict, DispDict, SpColorDict, N_RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, nNi, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList)

            p2, TNQ2 = metrics.getprod(Qs)
            PRODI = p2 - p1
            PRODN = TNQ2 - TNQ1


    ax = fig.add_subplot(111, projection='3d')
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    R, N = len(RIDs), len(SpeciesIDs)

    RDENS = R/(height*width)
    numD = ADList.count('d')
    numA = N - numD

    percD = 0
    if N > 0: percD = 100*(numD/N)

    Title = ['Individuals consume resources, grow, reproduce, and die as they move through the environment. \nAverage speed on the x-axis is '+str(u0)+' units per time step. \nActive N: '+str(numA)+', S: '+str(S)+', resources: '+str(R)+', ct: '+str(ct)+', %dormant: '+str(round(percD, 2)) + ', ' + str(len(Ns))]

    txt.set_text(' '.join(Title))
    ax.set_ylim(0, height)
    ax.set_xlim(0, width)

    if plot_system == 'yes':

        resource_scatImage.remove()
        Ind_scatImage.remove()

        colorlist = []
        sizelist = []
        for i, val in enumerate(SpeciesIDs):

            if ADList[i] == 'a':
                colorlist.append('r')
            else: colorlist.append('0.3')

            #colorlist.append(SpColorDict[val])
            sizelist.append(Qs[i][0] * 1000)

        resource_scatImage = ax.scatter(RX, RY, RZ,  c = 'w', edgecolor = 'SpringGreen', lw = 0.6, alpha=0.3, s = RVals*1)
        Ind_scatImage = ax.scatter(indX, indY, indZ, c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.9, s = sizelist)


    Ns.append(N)
    minct = 0
    tau = np.log10((width*length*height)/u1)

    #print 'sim:', '%4s' % sim, 'ct:', '%4s' % ct, 'tau:', '%5s' %  round(tau,2), ' volume:', '%4s' %  width**3,  '  N:', '%4s' %  N, 'R:', '%4s' % R, 'C:', '%4s' % numc, 'D:', '%4s' % percD

    if tau <= 2:     minct = 300
    elif tau <= 2.5: minct = 400
    elif tau <= 3:   minct = 500
    elif tau <= 3.5: minct = 600
    elif tau <= 4:   minct = 700
    elif tau <= 4.5: minct = 800
    elif tau <= 5:   minct = 900
    elif tau <= 5.5: minct = 1000
    else: minct = 1100

    if BurnIn == 'not done':
        if len(indX) == 0 or ct >= minct:
            BurnIn = 'done'
            Ns = [Ns[-1]] # only keep the most recent N value

    if BurnIn == 'done':

        RAD, splist, N, S = [], [], 0, 0
        if len(SpeciesIDs) >= 1:
            RAD, splist = bide.GetRAD(SpeciesIDs)
            S = len(RAD)

        RAD, splist = bide.GetRAD(SpeciesIDs)
        if len(RAD) > 1:
            RAD, splist = zip(*sorted(zip(RAD, splist), reverse=True))
        RAD = list(RAD)

        RTAU, INDTAU, TTAU = 0, 0, 0

        # Number of tracers, resource particles, and individuals
        T, R, N = len(TIDs), len(RIDs), len(SpeciesIDs)

        if N >= 1 and ct%10 == 0:

            spD = DispDict.values()
            spM = MaintDict.values()
            spG = GrowthDict.values()

            if len(spD) > 0: SpecDisp = mean(spD)
            if len(spM) > 0: SpecMaint = mean(spM)
            if len(spG) > 0: SpecGrowth = mean(spG)

            ES = metrics.e_simpson(RAD)

            wt = 0
            if len(Ns) == 1:
                splist2 = list(splist)
            if len(Ns) > 1:
                wt = metrics.WhittakersTurnover(splist, splist2)
                splist2 = list(splist)

            G = mean(GrowthList)
            M = mean(MaintList)
            avgMF = mean(MFDList)
            avgRPF = mean(RPFList)
            Disp = mean(DispList)

            Nmeans = [np.var(x) for x in zip(*N_RList)]
            NR = mean(Nmeans)

            OUT = open(GenPath + 'SimData.csv', 'a')

            outlist = [sim, ct, dorm, m, r, nNi, rmax, gmax, maintmax, dmax, \
            seedCom, u0, height, length, width, viscosity, N, PRODI, PRODN, RTAU, \
            TTAU, INDTAU, R, RDENS, S, ES, Nm, T, wt, G, M, NR, Disp, \
            amp, flux, freq, phase, SpecGrowth, SpecDisp, SpecMaint, numA, \
            numD, percD, avgRPF, avgMF]

            outlist = str(outlist).strip('[]')
            outlist = outlist.replace(" ", "")
            print>>OUT, outlist
            OUT.close()

            rad = str(RAD).strip('[]')
            rad = rad.replace(" ", "")
            OUT = open(GenPath + 'RAD-Data.csv', 'a')
            print>>OUT, sim, ',', ct,',',  rad
            OUT.close()


        if len(Ns) > 100:

            t = time.clock() - t
            print 'sim:', '%4s' % sim, 'tau:', '%5s' %  round(tau,2), 'volume:', '%4s' %  width**3,'  N:', '%4s' %  N, 'S:', '%4s' % S, 'R:', '%4s' % R, '%D:', '%4s' % round(percD,1)
            t = 0

            Rates = np.roll(Rates, -1, axis=0)
            u0 = Rates[0]
            u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))

            # Lists
            CRList, Ns, SpColorList, RColorList, RAD, splist, splist2, TIDs, TX, TY, RTypes, RX, RY, RZ, RIDs, RVals, SpeciesIDs, indX, indY, indZ, IndIDs, Qs, N_RD, GrowthList, MaintList, MFDList, RPFList, N_RList, DispList, ADList = [list([]) for _ in xrange(30)]
            # Scalars
            u1, numA, numD, RDENS, RDiv, RRich, S, ES, Ev, BP, SD, Nm, sk, Mu, Maint, IndID, RID, N, ct, T, R, PRODI, PRODN = [0]*23
            # Dictionaries
            SpColorDict, GrowthDict, MaintDict, MainFactorDict, RPFDict, N_RD, RColorDict, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}

            p = 0
            sim += 1
            BurnIn = 'not done'

            SpColorDict, GrowthDict, MaintDict, MainFactorDict, RPFDict, N_RD, RColorDict, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}

            if u0 == max(Rates):
                if len(Rates) > 1: print '\n'

                width, height, length, seedCom, m, r, nNi, rmax, gmax, maintmax, dmax, amp, freq, flux, pulse, phase, envgrads, Rates, pmax, mmax, dorm, imm = rp.get_rand_params()
                lefts, bottoms = [], []

                u1 = u0 + u0*(amp * sin(2*pi * ct * freq + phase))

            ####################### REPLACE ENVIRONMENT ########################
            ax = fig.add_subplot(111, projection='3d')



############### INITIALIZE GRAPHICS ############################################
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Ind_scatImage = ax.scatter([0], [0], [0])
resource_scatImage = ax.scatter([0], [0], [0])
Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=110, interval=40, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/residence-time/results/movies/examples/2015_10_05_1751_hydrobide.avi', bitrate=5000)
