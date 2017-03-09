# -*- coding: utf-8 -*-
from random import shuffle
from os.path import expanduser
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")

from diversity_metrics import *

def reproduce(u0, SpDict, IndDict, ID):

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
    
    Qs = []
    for key, value in IndDict.items(): 
        Qs.append(value['q'])
        p1, TNQ1 = diversity_metrics.getprod(Qs)
    
    for key, value in IndDict.items():
        state = value['state']
        sp = value['spID']
        q = value['q']
        smax = SpDict[sp]['s_max']
        sz = value['size']
        smin = SpDict[sp]['s_min']

        if state == 'd': continue

        if q < 1.0:
            del IndDict[key]
            continue

        if sz/2 < smin: continue
        if q/smax >= 1.0 or np.random.binomial(1, q/smax) == 1: # individual is large enough to reproduce
            x = value['x']
            y = value['y']
            z = value['z']


            f1 = np.random.uniform(0.5, 0.5)
            f2 = 1 - f1
            IndDict[key]['q'] = q*f1
            IndDict[ID] = {'q' : q*f2}
            IndDict[key]['size'] = sz*f1
            IndDict[ID]['size'] = sz*f2
            IndDict[ID]['state'] = 'a'
            IndDict[ID]['spID'] = sp

            mu = SpDict[sp]['spec']
            if mu > 0:
                p = np.random.binomial(1, mu)

                if p == 1: # speciate
                    new_sp = max(list(SpDict))+1
                    IndDict[ID]['spID'] = new_sp

                    # species growth rate
                    SpDict[new_sp] = {'grow' : SpDict[sp]['grow']}

                    # species active dispersal rate
                    SpDict[new_sp]['disp'] = SpDict[sp]['disp']

                    # species RPF factor
                    SpDict[new_sp]['rpf'] = SpDict[sp]['rpf']

                    # species maintenance
                    SpDict[new_sp]['maint'] = SpDict[sp]['maint']

                    # dormancy limit
                    SpDict[new_sp]['dlim'] = SpDict[sp]['dlim']

                    # species maintenance factor
                    SpDict[new_sp]['mfd'] = SpDict[sp]['mfd']

                    # species speciation rate
                    SpDict[new_sp]['spec'] = SpDict[sp]['spec']

                    SpDict[new_sp]['s_min'] = SpDict[sp]['s_min']
                    SpDict[new_sp]['s_max'] = SpDict[sp]['s_max']

                    # A set of specific growth rates for three major types of resources
                    SpDict[new_sp]['eff'] = SpDict[sp]['eff']

            IndDict[ID]['x'] = float(x)
            IndDict[ID]['y'] = float(y)
            IndDict[ID]['z'] = float(z)
            ID += 1
            
        Qs = []
        for key, value in IndDict.items(): Qs.append(value['q'])
        p2, TNQ2 = diversity_metrics.getprod(Qs)
    
        prodI = p2 - p1
        prodN = TNQ2 - TNQ1

    return [SpDict, IndDict, ID, prodI, prodN]
