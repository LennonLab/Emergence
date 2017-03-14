# -*- coding: utf-8 -*-
from __future__ import division
from random import randint
import numpy as np

def ResIn(rD, ps):

    ''' A function to bring resource particles into a system

        rD  :  A python dictionary holding properties of individual
                  resource particles:

                  'type' : Resource particles can belong to one of 3 types for
                  which species have varying abilities to grow on.

                  'size' : The size of the resource particle is an abstract
                  quantity, but can vary over two orders of magnitude

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
                u : rate of system flow through
    '''


    w, h, l, sd, m, r, n, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = ps

    try: rid = max(list(rD))
    except: rid = 1

    for i in range(r):
        x = np.random.binomial(1, u)

        if x == 1:
            rD[rid] = {'type' : randint(0, n-1)}
            rD[rid]['size'] = np.random.uniform(1, rm)
            rD[rid]['x'] = float(np.random.uniform(0, 0.9*h))
            rD[rid]['y'] = float(np.random.uniform(0, 0.9*h))
            rD[rid]['z'] = float(np.random.uniform(0, 0.9*h))
            rid += 1

    return rD
