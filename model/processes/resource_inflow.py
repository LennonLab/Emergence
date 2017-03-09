# -*- coding: utf-8 -*-
from __future__ import division
from random import randint
import numpy as np

def ResIn(RDict, params):
    
    ''' A function to bring resource particles into a system
    
        RDict  :  A python dictionary holding properties of individual
                  resource particles:
                   
                  'type' : Resource particles can belong to one of 3 types for
                  which species have varying abilities to grow on.
                      
                  'size' : The size of the resource particle is an abstract 
                  quantity, but can vary over two orders of magnitude
                  
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
    '''
    
    
    w, h, l, seed, m, r, nN, rmax, gmax, maintmax, dmax, amp, freq, phase, rate, pmax, dormlim, smax = params
    
    rid = max(list(RDict))
    
    for i in range(r):
        x = np.random.binomial(1, rate)

        if x == 1:
            RDict[rid] = {'type' : randint(0, nN-1)}
            RDict[rid]['size'] = np.random.uniform(1, rmax)
            RDict[rid]['x'] = float(np.random.uniform(0, 0.9*h))
            RDict[rid]['y'] = float(np.random.uniform(0, 0.9*h))
            RDict[rid]['z'] = float(np.random.uniform(0, 0.9*h))
            rid += 1

    return RDict