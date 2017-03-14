# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np


def ind_flow(sD, iD, ps):

    ''' A function to passively move individuals through a system

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

                   seed : Number of starting individuals
                   m : immigration rate and the probability of an individual
                   organism immigrating per time step

                   r : Maximum number of resource particles entering per time step
                   rmax : Maximum size of individual resource particles

                   nN : Number of inflowing resource types
                   gmax : Maximum specific growth rate
                   maintmax : Maximum metabolic maintenance
                   dmax : Maximum dispersal rate
                   pmax : maximum probability of random resuscitation
                   dormlim : level of endogenous resource at which
                   individual go dormant
                   smax : Maximum size of any individual

                   amp : amplitude of environmental flux
                   freq : frequency of environmental flux
                   phase : phase of individual immigration and resource inflow
                   rate : rate of system flow through

    '''

    w, h, l, sd, m, r, n, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = ps

    for k, val in iD.items():
        x = val['x']
        y = val['y']
        z = val['z']
        q = val['q']
        sp = val['spID']
        sz = val['size']
        d = sD[sp]['disp']

        x1, y1, z1 = 1, 1, 1
        if val['state'] == 'a':
            x1 = np.random.binomial(1, 1-d)
            y1 = np.random.binomial(1, 1-d)
            z1 = np.random.binomial(1, 1-d)
            q -= d * sz * u
            iD[k]['q'] = q

        x += u*x1
        y += u*y1
        z += u*z1

        if x > h or y > l or z > w:
            del iD[k]

        else:
            iD[k]['x'] = x
            iD[k]['y'] = y
            iD[k]['z'] = z

    return iD
