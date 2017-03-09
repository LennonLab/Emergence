from __future__ import division
from random import choice, randrange
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from scipy import stats


######################### COMMUNITY SIMULATION FUNCTION ########################
def simCOM(gens, dorm, disp, env, lgp = 0.95):
    COM = list(np.random.logseries(lgp, 1000)) # populate the community

    for i, microbe in enumerate(COM):
        x, y = np.random.uniform(0, 1, 2)
        COM[i] = [microbe, choice([0, 1]), x, y] # 0 is active; 1 is dormant

    for g in range(gens):

        i = randrange(len(COM))          # Choose an individual at random to die
	x = np.random.binomial(1, dorm)
	if x == 0 or COM[i][1] == 0: COM.pop(i)

        i = randrange(len(COM))             # Choose an individual to transition
        diff = np.abs(min([COM[i][2], COM[i][3]]) - env)
        COM[i][1] = np.random.binomial(1, diff)

        i = randrange(len(COM))               # Choose an individual to disperse
        x = np.random.normal(loc=COM[i][2], scale=disp)
	y = np.random.normal(loc=COM[i][3], scale=disp)

	if x > 1: x = 1 # keeping the individual within bounds
	elif x < 0: x = 0
	if y > 1: y = 1
	elif y < 1: y = 1
	COM[i][2], COM[i][3] = [x, y]

	for i, val in enumerate(COM):
            if val[1] == 0: # Choose an individual at random to reproduce
                COM.append(val)
                break

        print g, len(COM)
    return COM

################ GET SLOPES AND INTERCEPTS ####### Regression Approach #########
def getVals(AlistsA, SlistsA, ZvalsA, YintsA, AlistsT,
                            SlistsT, ZvalsT, YintsT, COM):


    SlistA, SlistT, AlistA, AlistT, grain = [ [], [], [],[], 2]
    ZvalsA, YintsA, ZvalsT, YintsT = [ [], [], [], [] ]

    for i in range(1, 10):
        slistA, slistT, grain = [[], [], grain/2] # grain on the log2 scale

        for ind in COM:
            if min([ind[2],ind[3]])<= grain:
                slistT.append(ind[0])
                if ind[1]==0: slistA.append(ind[0])

        if len(list(set(slistT))) > 0:
            SlistT.append(len(list(set(slistT))))
            AlistT.append(grain)

            if len(list(set(slistA))) > 0:
                SlistA.append(len(list(set(slistA))))
                AlistA.append(grain)

    ZvalA, YintA, r, p, s = stats.linregress(np.log(AlistA), np.log(SlistA))
    ZvalsA.append(ZvalA)
    YintsA.append(YintA)
    AlistsA.append(AlistA)
    SlistsA.append(SlistA)

    ZvalT, YintT, r, p, s = stats.linregress(np.log(AlistT), np.log(SlistT))
    ZvalsT.append(ZvalT)
    YintsT.append(YintT)
    AlistsT.append(AlistT)
    SlistsT.append(SlistT)


    return [AlistsA, SlistsA, ZvalsA, YintsA, AlistsT, SlistsT, ZvalsT, YintsT]

####################  GENERATE DATA FOR FIGURES  ###############################
AlistsA, SlistsA, ZvalsA, YintsA, AlistsT, SlistsT, ZvalsT, YintsT, gens = [ [],
                                                [], [], [], [], [], [],[], 5000 ]

Dorm, Disp, Env = [ 0.9, 0.1, 0.9]
for i in range(10):
    COM = simCOM(gens, Dorm, Disp, Env)
    AlistsA, SlistsA, ZvalsA, YintsA, AlistsT, SlistsT, ZvalsT, YintsT = getVals(AlistsA, SlistsA, ZvalsA, YintsA, AlistsT, SlistsT, ZvalsT, YintsT, COM)

##############################  FIGURES  #######################################
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

plt.xlabel('log(A)', fontsize = 16)
plt.ylabel('log(S)', fontsize = 16)

for i, slistA in enumerate(SlistsA):
    if i != 0:
        plt.plot(AlistsA[i], slistA, color='gray')
    elif i == 0:
        plt.plot(AlistsA[i], slistA, color='gray', label='Active')

    ax.tick_params(axis='both', labelsize=10)
    plt.yscale('log')
    plt.xscale('log')

for i, slistT in enumerate(SlistsT):
    if i !=0 :
        plt.plot(AlistsT[i], slistT, color='k')
    elif i == 0:
        plt.plot(AlistsT[i], slistT, color='k', label='All')
    ax.tick_params(axis='both', labelsize=14)
    plt.yscale('log')
    plt.xscale('log')

leg = plt.legend(loc=2,prop={'size':16})
leg.draw_frame(False)

plt.savefig('/Users/lisalocey/Desktop/TARs_simple.png',dpi=600,bbox_inches='tight',pad_inches=0.1)
plt.show()
