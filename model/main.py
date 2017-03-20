from __future__ import division
from random import shuffle, seed, randint
from os.path import expanduser
import sys
import numpy as np

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")

from processes import *
from diversity_metrics import *
from spatial_functions import *
from input_output import *

labels.clear()

def iter_procs(procs, iD, sD, rD, ps, ct, ceil = 2000):
    pr = 0

    shuffle(procs)
    for p in procs:

        if p is 'resource_inflow': # Inflow of resources
            rD = resource_inflow.ResIn(rD, ps)

        elif p is 'resource_flow': # Resource flow
            rD = resource_flow.res_flow(rD, ps)

        elif p is 'immigration' and ps[2] > 0: # Inflow of individuals (immigration)
            sD, iD = immigration.immigration(sD, iD, ps)

        elif p is 'passive_dispersal': # flowthrough of individuals
            iD = passive_dispersal.ind_flow(iD, ps)

        elif p is 'active_dispersal': # Active dispersal
            iD = active_dispersal.ind_disp(iD, ps)

        elif p is 'consume': # Consume
            iD, rD = consume.consume(iD, rD, ps)

        elif p is 'growth': # Grow
            iD = growth.grow(iD)

        elif p is 'transition': # Transition
            iD = transition.transition(iD)

        elif p is 'maintenance': # Maintenance
            iD = maintenance.maintenance(iD)

        elif p is 'reproduction': # Reproduction
            sD, iD, pr = reproduction.reproduce(sD, iD, ps)

        elif p is 'disturb' and len(list(iD)) > ceil:
            iD = disturb.disturb(iD, ceil)

    N, R = len(list(iD)), len(list(rD))
    return [iD, sD, rD, N, R, ct+1, pr]



def run_model(procs, sim, rD = {}, sD = {}, iD = {}, ct = 0, splist2 = []):
    seed()
    print '\n'
    h = randint(1, 101)
    l = int(h)
    r = randint(1, 101)
    ps = h, l, r, 10**np.random.uniform(-3, 0)
    sD, iD = immigration.immigration(sD, iD, ps, 1000)

    while ct < 400:
        iD, sD, rD, N, R, ct, prod = iter_procs(procs, iD, sD, rD, ps, ct)
        N, S, R, splist2 = output.output(iD, sD, rD, ps, sim, N, R, ct, prod, splist2)
        if N == 0: break

procs = labels.processes()
for sim in range(10**5): run_model(procs, sim)
