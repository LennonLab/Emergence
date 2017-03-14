# -*- coding: utf-8 -*-
def maintenance(sD, iD):

    ''' A function to bring resource particles into a system

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
    '''


    for k, val in iD.items():

        sp = val['spID']
        q = val['q']
        sz = val['size']
        smin = sD[sp]['s_min']
        maint = sD[sp]['maint']

        q -= maint
        if q <= 0 or sz < smin:
            del iD[k]
        else:
            iD[k]['q'] = q

    return iD
