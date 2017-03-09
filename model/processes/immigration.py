# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import sys
import numpy as np



def immigration(SpDict, IndDict, params, ct):

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
    
    
    w, h, l, seed, m, r, nN, rmax, g_max, m_max, d_max, amp, freq, phase, rates, p_max, dlim, s_max = params
    
    sd = int(seed)
    if ct > 1:
        sd = 1

    for j in range(sd):

        if m == 0 and ct > 1: break
        if sd == 1 and np.random.binomial(1, u1) == 0: continue

        prop = np.random.randint(1, 1000)
        if prop not in SpDict:

            # species max & min body size
            smax = np.random.uniform(100, 100)
            smin = np.random.uniform(1, 1)
            SpDict[prop] = {'s_max' : smax}
            SpDict[prop]['s_min'] = smin

            # species growth rate
            SpDict[prop]['grow'] = np.random.uniform(g_max/100, g_max)

            # species active dispersal rate
            SpDict[prop]['disp'] = np.random.uniform(d_max/100, d_max)

            # species RPF factor
            SpDict[prop]['rpf'] = np.random.uniform(p_max/100, p_max)

            # species maintenance
            SpDict[prop]['maint'] = np.random.uniform(m_max/100, m_max)

            # dormancy limit
            SpDict[prop]['dlim'] = np.random.uniform(dlim/100, dlim)

            # species speciation rate
            SpDict[prop]['spec'] = np.random.uniform(0.01, 0.01)

            # species maintenance factor
            mfd = np.random.logseries(0.95, 1)[0]
            if mfd > 40: mfd = 40
            SpDict[prop]['mfd'] = mfd + 1

            # A set of specific growth rates for three major types of resources
            SpDict[prop]['eff'] = choice([[1/3, 1/3, 1/3], [0.5, 0.3, 0.2], [0.6, 0.2, 0.2], [0.8, 0.1, 0.1], [0.99, 0.005, 0.005], [0.9, 0.05, 0.05]])
            
            
        IndDict[ID] = {'spID' : prop}
        IndDict[ID]['x'] = float(np.random.uniform(0, 0.9*h))
        IndDict[ID]['y'] = float(np.random.uniform(0, 0.9*l))
        IndDict[ID]['z'] = float(np.random.uniform(0, 0.9*w))

        smin = SpDict[prop]['s_min']
        smax = SpDict[prop]['s_max']

        size = float(np.random.uniform(smin, smax))
        IndDict[ID]['size'] = size
        q = size - smin
        IndDict[ID]['q'] = q

        IndDict[ID]['state'] = 'a'
        ID += 1

    return [SpDict, IndDict]