# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import sys
import numpy as np
from os.path import expanduser
import time


mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")

from diversity_metrics import *
from spatial_functions import * 

def consume(SpDict, IndDict, h, l, w, u0, Rtypes, Rvals, RX, RY, RZ, RIDs):
    
        ''' A function to bring resource particles into a system
    
        SpDict  :  A python dictionary holding properties of specific species:
                   
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
                     
                     
        IndDict  : A python dictionary to hold properties of individual
                   organisms            
            
                   'size' : The size of the individual is an abstract 
                   quantity, but can vary over two orders of magnitude
                  
                   'q' : level of endogenous resources
                   'spID' : species ID
                   'state' : whether active or dormant
                  
                  'x' : x-coordinate
                  'y' : y-coordinate
                  'z' : z-coordinate
        
        
        params  :  General model parameters. Not all are used in every function.
        
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
                   
        ct  :  indicates whether the model is on its first time step
    '''
    
    numc = 0

    for key, value in IndDict.items():
        sp = value['spID']
        state = value['state']
        x1 = value['x']
        y1 = value['y']
        z1 = value['z']
        sp = value['spID']
        q = value['q']
        sz = value['size']
        mfd = SpDict[sp]['mfd']
        eff = SpDict[sp]['eff']
        maint = SpDict[sp]['maint']
        rpf = SpDict[sp]['rpf']

        if q <= 0:
            del IndDict[key]
            continue

        coords = [x1, y1, z1]
        if len(RIDs):
            closest = get_closest(RIDs, RX, RY, RZ, Rtypes, coords)
            ri = RIDs.index(closest)
            x2 = RX[ri]
            y2 = RY[ri]
            z2 = RZ[ri]
            Rval = Rvals[ri]
            rtype = Rtypes[ri]

        else: continue

        dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        i_radius = ((0.75*sz)/np.pi)**(1.0/3)
        r_radius = ((0.75*Rval*1)/np.pi)**(1.0/3)

        if dist <= i_radius + r_radius:
            if state == 'd':
                x = np.random.binomial(1, rpf)
                if x == 0:
                    continue
                elif x == 1:
                    IndDict[key]['state'] == 'a'
                    IndDict[key]['maint'] = maint*mfd
                    IndDict[key]['q'] -= q * rpf
        else: continue

        numc += 1
        e = eff[rtype]
        if Rval > e * q: # Increase cell quota
            Rval -= e * q
            q += e * q
        else:
            q += Rval
            Rval = 0.0

        if Rval <= 0.0:
            Rvals.pop(ri)
            Rtypes.pop(ri)
            RIDs.pop(ri)
            RX.pop(ri)
            RY.pop(ri)
            RZ.pop(ri)
        else:
            Rvals[ri] = Rval

        IndDict[key]['q'] = q

    return [numc, SpDict, IndDict, h, l, w, u0, Rtypes, Rvals, RX, RY, RZ, RIDs]
