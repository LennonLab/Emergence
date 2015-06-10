from __future__ import division
import  matplotlib.pyplot as plt
import datetime
import numpy as np
import os.path
import Hydrobide
import sys

mydir = '/Users/lisalocey/Desktop/repos/HYDRO-BIDE/'
sys.path.append(mydir + '/tools/figurefunctions')
import figs

fig = plt.figure()

r = 2 # number rows
c = 2 # number columns
wsp = 0.39 # width space
hsp = 0.39 # height space

use_picloud = 'n'
basic = True   # if true, will not return info on richness, evenness, etc.

Tau = 'random' # can also be 'list'
replicates = 1

systems = [1]
time = 5000 # number of generations. Needs to be large enough to allow
            # the population to reach a relatively stable state


i=1
if 1 in systems:
    # IDEAL CHEMOSTAT
    system = 'ideal'
    Glist, Nlist, Plist, Rlist, MCTlist, Qlist, TauList, mean, tauAtQmaxList = Hydrobide.run_hydrobide(use_picloud, system, Tau, basic, replicates, time)
    
    fig = figs.MCTvsTau(fig, r, c, i, wsp, hsp, MCTlist, TauList, tauAtQmaxList)
    i+=1
    #fig = figs.MUvsD(fig, r, c, i, wsp, hsp, Plist, TauList, tauAtQmaxList)
    fig = figs.BTvsTau(fig, r, c, i, wsp, hsp, Plist, TauList, tauAtQmaxList)


if 2 in systems:
    # IMMIGRATION AND SPECIES DIFFERENCES
    system = 'imNonNeutral'
    Glist, Nlist, Plist, Rlist, MCTlist, Qlist, TauList, mean, tauAtQmaxList = Hydrobide.run_hydrobide(use_picloud, system, Tau, basic, replicates, time)
    
    i+=1
    fig = figs.MCTvsTau(fig, r, c, i, wsp, hsp, MCTlist, TauList, tauAtQmaxList)
    i+=1
    #fig = figs.MUvsD(fig, r, c, i, wsp, hsp, Plist, TauList, tauAtQmaxList)
    fig = figs.BTvsTau(fig, r, c, i, wsp, hsp, Plist, TauList, tauAtQmaxList)


if 3 in systems:
    # COMPLEX ECOSYSTEM
    system = 'complex_ecosystem'
    Glist, Nlist, Plist, Rlist, MCTlist, Qlist, TauList, mean, tauAtQmaxList = Hydrobide.run_hydrobide(use_picloud, system, Tau, basic, replicates, time)
    
    i+=1
    fig = figs.MCTvsTau(fig, r, c, i, wsp, hsp, MCTlist, TauList, tauAtQmaxList)
    i+=1
    fig = figs.MUvsD(fig, r, c, i, wsp, hsp, Plist, TauList, tauAtQmaxList)
    

# save figure
d = datetime.datetime.now()
label = str()
for attr in [ 'year', 'month', 'day', 'hour', 'minute']:
   label += '_'+str(getattr(d, attr))
   
plt.savefig(mydir + 'results/figs/systems'+label+'2.png', dpi=600, bbox_inches = 'tight', pad_inches=0.03)
plt.show()