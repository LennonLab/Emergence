# -*- coding: utf-8 -*-

def res_flow(rD, ps):

    ''' A function to move resource particles through a system

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

    for k, val in rD.items():
        x = val['x']+u
        y = val['y']+u
        z = val['z']+u

        if x > h or y > l or z > w:
            del rD[k]
        else:
            rD[k]['x'] = x
            rD[k]['y'] = y
            rD[k]['z'] = z

    return rD
