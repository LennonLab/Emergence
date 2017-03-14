from random import shuffle
from os.path import expanduser
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")

from processes import *
from diversity_metrics import *
from randomized_parameters import *
from spatial_functions import *
from input_output import *

labels.clear()


def iter_procs(procs, iD, sD, rD, ps, ct, ceiling = 1000):

    w, h, l, sd, m, r, nN, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = ps
    prodI, prodN, numc = 0, 0, 0

    shuffle(procs)
    for p in procs:

        if p is 'resource_inflow': # Inflow of resources
            rD = resource_inflow.ResIn(rD, ps)

        elif p is 'resource_flow': # Resource flow
            rD = resource_flow.res_flow(rD, ps)

        elif p is 'immigration': # Inflow of individuals (immigration)
            sD, iD = immigration.immigration(sD, iD, ps, ct)

        elif p is 'passive_dispersal': # flowthrough of individuals
            iD = passive_dispersal.ind_flow(sD, iD, ps)

        elif p is 'active_dispersal': # Active dispersal
            iD = active_dispersal.ind_disp(sD, iD, ps)

        elif p is 'search': # Search
            iD = active_dispersal.search(sD, iD, rD, ps)

        elif p is 'consume': # Consume
            numc, iD, rD = consume.consume(sD, iD, rD, ps)

        elif p is 'growth': # Grow
            iD = growth.grow(sD, iD)

        elif p is 'transition': # Transition
            iD = transition.transition(sD, iD)

        elif p is 'maintenance': # Maintenance
            iD = maintenance.maintenance(sD, iD)

        elif p is 'reproduction': # Reproduction
            sD, iD, pI, pN = reproduction.reproduce(sD, iD, ps)

        elif p is 'disturb' and len(list(iD)) > ceiling:
            iD = disturb.disturb(iD, ceiling)

    N = len(list(iD))
    R = len(list(rD))

    return [iD, sD, rD, ct+1, N, R, [ct+1, numc, pI, pN]]



def run_model(procs, sim):
    ps = randparams.get_rand_params()
    sD, iD, rD, splist2, Ns, N, S, R, ct = {}, {}, {}, [], [], 0, 0, 0, 0

    while ct < 600:
        iD, sD, rD, ct, N, R, ls = iter_procs(procs, iD, sD, rD, ps, ct)
        if ct > 200 and ct%10 == 0: Ns, N, S, R, splist2 = output.output(iD, sD, rD, ps, sim, ls, Ns, splist2)
        print 'sim:', '%3s' % sim, 'ct:', '%3s' % ct,'  N:', '%4s' %  N, '  S:', '%4s' %  S, '  R:', '%4s' %  R, ' u0:', '%4s' % round(ps[14], 6)

        if N == 0: break

procs = labels.processes()
for sim in range(10**5): run_model(procs, sim)
