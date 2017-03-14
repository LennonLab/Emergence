# -*- coding: utf-8 -*-
from __future__ import division
from random import choice
import numpy as np

def immigration(sD, iD, ps, ct):

    ''' A function to bring individual organisms into a system

        sD  :  A python dictionary holding properties of specific species:

                   's_min' : lower bound on individual size
                   's_max' : upper bound on individual size

                   'grow' :  maximum growth rate
                   'disp' :  maximum dispersal rate

                   'rpf' :  probability of random transition to active state
                   'maint' : basal metabolic maintenance cost

                   'dlim' : lower size limit at which individuals go dormant
                   'spec' : speciation rate

                   'mfd' : factor by which dormant decreases costs of
                   maintenance energy


                   'eff' : Resource particles can belong to one of n types for
                   which species have varying abilities to grow on.


        iD  : A python dictionary to hold properties of individual
                   organisms

                   'size' : The size of the individual is an abstract
                   quantity, but can vary over two orders of magnitude

                   'q' : level of endogenous resources
                   'spID' : species ID
                   'state' : whether active or dormant

                  'x' : x-coordinate
                  'y' : y-coordinate
                  'z' : z-coordinate


        ps  :  General model parameters. Not all are used in every function.

                   w : width of the system
                   h : height of the system
                   l : length of the system

                   sd : Number of starting individuals
                   m : immigration rate and the probability of an individual
                       organism immigrating per time step

                   r : Maximum number of resource particles entering per time step
                   rm : Maximum size of individual resource particles

                   n : Number of inflowing resource types
                   gm : Maximum specific growth rate
                   mm : Maximum metabolic maintenance
                   dm : Maximum dispersal rate
                   pm : maximum probability of random resuscitation
                   dl : level of endogenous resource at which
                        individual go dormant
                   sm : Maximum size of any individual

                   a : amplitude of environmental flux
                   f : frequency of environmental flux
                   p : phase of individual immigration and resource inflow
                   rate : rate of system flow through

        ct  :  indicates whether the model is on its first time step
    '''

    w, h, l, sd, m, r, n, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = ps

    if ct > 0: sd = 1
    for j in range(sd):

        if m == 0 and ct > 0: break
        if sd == 1 and np.random.binomial(1, u) == 0: continue

        p = np.random.randint(1, 1000)
        if p not in sD:
            # species max & min body size
            smax = np.random.uniform(sm, sm)
            smin = np.random.uniform(1, 1)
            sD[p] = {'s_max' : smax}
            sD[p]['s_min'] = smin

            # species growth rate
            sD[p]['grow'] = np.random.uniform(gm/100, gm)

            # species active dispersal rate
            sD[p]['disp'] = np.random.uniform(dm/100, dm)

            # species RPF factor
            sD[p]['rpf'] = np.random.uniform(pm/100, pm)

            # species maintenance
            sD[p]['maint'] = np.random.uniform(mm/1000, mm)

            # dormancy limit
            sD[p]['dlim'] = np.random.uniform(dl/100, dl)

            # species speciation rate
            sD[p]['spec'] = np.random.uniform(0.01, 0.01)

            # trophic interactions
            #sD[p]['produces'] = np.random.randint(1, 1000)
            #sD[p]['consumes'] = np.random.randint(1, 1000)
            #sD[p]['parasitizes'] = np.random.randint(1, 1000)

            # species maintenance factor
            mfd = np.random.logseries(0.95, 1)[0]
            sD[p]['mfd'] = mfd + 1

            # A set of specific growth rates for three major types of resources
            sD[p]['eff'] = choice([[1/3, 1/3, 1/3], [0.5, 0.3, 0.2], [0.6, 0.2, 0.2],
                            [0.8, 0.1, 0.1], [0.99, 0.005, 0.005], [0.9, 0.05, 0.05]])

        i = id(p)
        iD[i] = {'spID' : p}
        iD[i]['x'] = float(np.random.uniform(0, 0.9*h))
        iD[i]['y'] = float(np.random.uniform(0, 0.9*l))
        iD[i]['z'] = float(np.random.uniform(0, 0.9*w))

        smin = sD[p]['s_min']
        smax = sD[p]['s_max']
        size = np.random.uniform(smin, smax)

        iD[i]['size'] = size # every individual starts at minimum size
        iD[i]['q'] = size - smin # but with some amount of endog resources
        iD[i]['state'] = 'a'

    return [sD, iD]
