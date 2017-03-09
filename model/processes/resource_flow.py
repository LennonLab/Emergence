# -*- coding: utf-8 -*-
import numpy as np

def res_flow(RTypes, RVals, RX, RY, RZ, RID, RIDs, params, u0):

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

    RX = np.array(RX) + u0
    RY = np.array(RY) + u0
    RZ = np.array(RZ) + u0

    i1 = np.where(RX > h)[0].tolist()
    i2 = np.where(RY > l)[0].tolist()
    i3 = np.where(RZ > w)[0].tolist()
    index = np.array(list(set(i1 + i2 + i3)))

    RTypes = np.delete(RTypes, index).tolist()
    RX = np.delete(RX, index).tolist()
    RY = np.delete(RY, index).tolist()
    RZ = np.delete(RZ, index).tolist()
    RVals = np.delete(RVals, index).tolist()
    RIDs = np.delete(RIDs, index).tolist()

    return [RTypes, RVals, RX, RY, RZ, RIDs, RID]