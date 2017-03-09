# -*- coding: utf-8 -*-
from __future__ import division
import sys
import numpy as np
from os.path import expanduser


mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")


def ind_flow(SpDict, IndDict, h, l, w, u0):
    
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
    

    for key, value in IndDict.items():

        x = value['x']
        y = value['y']
        z = value['z']
        q = value['q']
        sp = value['spID']
        sz = value['size']
        d = SpDict[sp]['disp']

        x1, y1, z1 = 1, 1, 1
        #if value['state'] == 'a':
        #    x1 = np.random.binomial(1, 1-d)
        #    y1 = np.random.binomial(1, 1-d)
        #    z1 = np.random.binomial(1, 1-d)
        #    q -= d * sz * u0
        #    IndDict[key]['q'] = q

        x += u0*x1
        y += u0*y1
        z += u0*z1

        if x > h or y > l or z > w:
            del IndDict[key]

        else:
            IndDict[key]['x'] = x
            IndDict[key]['y'] = y
            IndDict[key]['z'] = z

    return IndDict
