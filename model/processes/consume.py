# -*- coding: utf-8 -*-
from __future__ import division
from random import choice
import numpy as np


def get_closest(rD, cs):

    rIDs = list(rD)
    cl = 0
    if len(rIDs) == 0: return cl

    x1, y1, z1 = cs
    t = min([20, len(rIDs)])
    mD, ct = 10**10, 0

    for ct in range(t):
        j = choice(rIDs)
        x = rD[j]['x']
        y = rD[j]['y']
        z = rD[j]['z']

        d = np.sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)

        if d < mD:
            mD = d
            c = j

        return c



def consume(sD, iD, rD, params):

    ''' A function to simulate consumption of resource particles

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

        ct  :  indicates whether the model is on its first time step
    '''

    w, h, l, sd, m, r, n, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = params

    nc = 0

    for k, val in iD.items():
        sp = val['spID']
        st = val['state']
        x1 = val['x']
        y1 = val['y']
        z1 = val['z']
        sp = val['spID']
        q = val['q']
        sz = val['size']
        mf = sD[sp]['mfd']
        ef = sD[sp]['eff']
        mt = sD[sp]['maint']
        rp = sD[sp]['rpf']

        if q <= 0:
            del iD[k]
            continue

        cd = [x1, y1, z1]
        cl = get_closest(rD, cd)
        if cl == 0: continue

        x2 = rD[cl]['x']
        y2 = rD[cl]['y']
        z2 = rD[cl]['z']
        Rv = rD[cl]['size']
        rt = rD[cl]['type']

        d = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        ir = ((0.75*sz)/np.pi)**(1.0/3)
        rr = ((0.75*Rv*1)/np.pi)**(1.0/3)

        if d <= ir + rr:
            if st == 'd':
                x = np.random.binomial(1, rp)
                if x == 0:
                    continue
                elif x == 1:
                    iD[k]['state'] == 'a'
                    iD[k]['maint'] = mt*mf
                    iD[k]['q'] -= q * rp
        else: continue

        nc += 1
        e = ef[rt]
        if Rv > e * q: # Increase cell quota
            Rv -= e * q
            q += e * q
        else:
            q += Rv
            Rv = 0.0

        if Rv <= 0.0:
            del rD[cl]

        iD[k]['q'] = q

    return [nc, iD, rD]
