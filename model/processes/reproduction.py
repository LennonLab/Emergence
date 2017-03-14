# -*- coding: utf-8 -*-
import numpy as np

def getprod(Qs):

    N = 0
    if len(Qs) == 0:
        return [0, N]

    if len(Qs) > 0:
        for q in Qs:
            N += q

        p1 = len(Qs)
        return [p1, N]



def reproduce(sD, iD, ps):

    ''' A function to simulate reproduction

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

    w, h, l, sd, m, r, n, rm, gm, mm, d, a, f, p, u, pm, dl, sm, st = ps

    try: ID = max(list(iD))
    except: ID = 1

    p1, p2, b1, b2 = 0, 0, 0, 0

    Qs = []
    for key, value in iD.items():
        Qs.append(value['q'])
        p1, b1 = getprod(Qs)

    for key, value in iD.items():
        state = value['state']
        sp = value['spID']
        q = value['q']
        smax = sD[sp]['s_max']
        sz = value['size']
        smin = sD[sp]['s_min']

        if state == 'd': continue

        if q < smin:
            del iD[key]
            continue

        if sz/2 < smin: continue

        if q > smax or np.random.binomial(1, q/smax) == 1: # individual is large enough to reproduce
            x = value['x']
            y = value['y']
            z = value['z']

            f1 = np.random.uniform(0.5, 0.5)
            f2 = 1 - f1
            iD[key]['q'] = q*f1
            iD[ID] = {'q' : q*f2}
            iD[key]['size'] = sz*f1
            iD[ID]['size'] = sz*f2
            iD[ID]['state'] = 'a'
            iD[ID]['spID'] = sp

            mu = sD[sp]['spec']
            if mu > 0:
                p = np.random.binomial(1, mu)

                if p == 1: # speciate
                    new_sp = max(list(sD))+1
                    iD[ID]['spID'] = new_sp

                    # species growth rate
                    sD[new_sp] = {'grow' : sD[sp]['grow']}

                    # species active dispersal rate
                    sD[new_sp]['disp'] = sD[sp]['disp']

                    # species RPF factor
                    sD[new_sp]['rpf'] = sD[sp]['rpf']

                    # species maintenance
                    sD[new_sp]['maint'] = sD[sp]['maint']

                    # dormancy limit
                    sD[new_sp]['dlim'] = sD[sp]['dlim']

                    # species maintenance factor
                    sD[new_sp]['mfd'] = sD[sp]['mfd']

                    # species speciation rate
                    sD[new_sp]['spec'] = sD[sp]['spec']

                    sD[new_sp]['s_min'] = sD[sp]['s_min']
                    sD[new_sp]['s_max'] = sD[sp]['s_max']

                    # A set of specific growth rates for three major types of resources
                    sD[new_sp]['eff'] = sD[sp]['eff']

            iD[ID]['x'] = float(x)
            iD[ID]['y'] = float(y)
            iD[ID]['z'] = float(z)
            ID += 1

    for key, value in iD.items(): Qs.append(value['q'])
    p2, b2 = getprod(Qs)

    pI = p2 - p1
    pN = b2 - b1

    return [sD, iD, pI, pN]
