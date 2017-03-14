# -*- coding: utf-8 -*-
from __future__ import division
from random import choice
import numpy as np


''' INPUTS TO FUNCTIONS BELOW:

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
                   u : rate of system flow through

        ct  :  indicates whether the model is on its first time step '''




def get_closest(rD, cs):

    rIDs = list(rD)
    cl = 0
    if len(rIDs) == 0: return cl

    x1, y1, z1 = cs
    t = min([20, len(rIDs)])
    mD = 10**10

    for ct in range(t):
        j = choice(rIDs)
        x = rD[j]['x']
        y = rD[j]['y']
        z = rD[j]['z']

        d = np.sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)

        if d < mD:
            mD = d
            cl = j

        return cl




def ind_disp(sD, iD, ps):

    ''' Simulate active resistence of individuals to a flowing system '''

    w, h, l, sd, m, r, n, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = ps

    for k, val in iD.items():
        x = val['x']
        y = val['y']
        z = val['z']
        q = val['q']
        sp = val['spID']
        sz = val['size']
        d = sD[sp]['disp']

        x1, y1, z1 = 0, 0, 0
        if val['state'] == 'a':
            x1 = np.random.binomial(1, d)
            y1 = np.random.binomial(1, d)
            z1 = np.random.binomial(1, d)
            q -= d * sz * u
            iD[k]['q'] = q

        x -= x1*d
        y -= y1*d
        z -= z1*d

        if x > h or y > l or z > w:
            del iD[k]

        else:
            iD[k]['x'] = x
            iD[k]['y'] = y
            iD[k]['z'] = z

    return iD




def search(sD, iD, rD, ps):

    ''' Simulate active searching of individuals for resources '''

    w, h, l, sd, m, r, n, rm, gm, mm, dm, a, f, p, u, pm, dl, sm, st = ps

    for k, val in iD.items():

        x = val['x']
        y = val['y']
        z = val['z']
        sp = val['spID']
        q = val['q']
        d = sD[sp]['disp']
        sz = val['size']

        if val['state'] == 'd': continue

        cs = [x, y, z]
        c = get_closest(rD, cs)

        if c == 0:
            x1 = np.random.uniform(0, h)
            y1 = np.random.uniform(0, l)
            z1 = np.random.uniform(0, w)

        else:
            x1 = rD[c]['x']
            y1 = rD[c]['y']
            z1 = rD[c]['z']

        x2 = abs(x - x1)
        y2 = abs(y - y1)
        z2 = abs(z - z1)

        q -= d * sz
        iD[k]['q'] = q

        if x > x1:
            x -= np.random.uniform(0, d*x2)
        elif x < x1:
            x += np.random.uniform(0, d*x2)

        if y > y1:
            y -= np.random.uniform(0, d*y2)
        elif y < y1:
            y += np.random.uniform(0, d*y2)

        if z > z1:
            z -= np.random.uniform(0, d*z2)
        elif z < z1:
            z += np.random.uniform(0, d*z2)

        iD[k]['x'] = x
        iD[k]['y'] = y
        iD[k]['z'] = z

    return iD
