from __future__ import division
from random import choice, randrange
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from scipy import stats


######################### COMMUNITY SIMULATION FUNCTION ########################
def simCOM(gens, dorm, disp, env, lgp = 0.92):
    COM = list(np.random.logseries(lgp, 2000)) # populate the community

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
def getVals(Zvals, Yints, COM):
    Slist, Alist, grain = [ [], [], 2]

    for i in range(1, 10):
        slist, grain = [[], grain/2] # grain on the log2 scale

        for ind in COM:
            if min([ind[2], ind[3]]) <= grain: slist.append(ind[0])

        if len(list(set(slist))) > 0:
            Slist.append(len(list(set(slist))))
            Alist.append(grain)

    Zval, Yint, r, p, s = stats.linregress(np.log(Alist), np.log(Slist))
    Zvals.append(Zval)
    Yints.append(Yint)
    return [Zvals, Yints]

####################  GENERATE DATA FOR FIGURES  ###############################
inds, num = [ [0, 0, 1, 1, 2, 2],  10]
DormVals, DispVals, EnvVals = [ np.array(range(1, num+1)) * 0.1 ] * 3
Zvals, Yints, gens = [ [],  [],  20000]

for i in DormVals: # rows starting from bottom
    for j in DispVals: # columns starting from left
        COM = simCOM(gens, i, j, 0.5)  # dormancy & dispersal
        Zvals, Yints = getVals(Zvals, Yints, COM)

for i in DormVals:  # rows starting from bottom
    for j in EnvVals:  # columns starting from left
        COM = simCOM(gens, i, 0.5, j) # dormancy & environment
        Zvals, Yints = getVals(Zvals, Yints, COM)

for i in DispVals: # rows starting from bottom
    for j in EnvVals:  # columns starting from left
        COM = simCOM(gens, 0.5, i, j) # dispersal & environment
        Zvals, Yints = getVals(Zvals, Yints, COM)

##############################  FIGURES  #######################################
fig = plt.figure()
Zvals, Yints =  [np.reshape(Zvals, (3, num**2)), np.reshape(Yints, (3, num**2))]

ylabels = ['Dormant\ncapacity']*4 + ['Dispersal\ncapacity']*2
xlabels = ['Dispersal capacity']*2 + ['Environmmental\nHeterogeneity']*4
heatcolors = ['Reds', 'Reds', 'Blues', 'Blues', 'jet', 'jet']

for i in range(6):
    ax = fig.add_subplot(3, 2, i+1)

    if i==0: plt.title('Z-values', fontsize=14)
    elif i==1: plt.title('Y-ints', fontsize=14)

    if i%2 == 0: vals = np.reshape(Zvals[inds[i]], (num, num))
    elif i%2 > 0: vals = np.reshape(Yints[inds[i]], (num, num))

    plt.imshow(vals, origin='lower', cmap = plt.get_cmap(heatcolors[i]))
    plt.xlabel(xlabels[i], fontsize = 10)
    plt.ylabel(ylabels[i], fontsize = 10)

    ax.tick_params(axis='both', labelsize=10)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=8)

plt.subplots_adjust(wspace=0.00, hspace=0.55)
plt.savefig('/Users/lisalocey/Desktop/heat.png',dpi=600,bbox_inches='tight',pad_inches=0.1)
