from __future__ import division
import sys                                            
import numpy as np
from random import randrange
import math
import random



def get_maxUptake(maxUptake_dict, prop, mean, std, system):
    
    maxUptake = mean
    
    #if system != 'ideal': # uncomment to make 'ideal' systems neutral
    while maxUptake <= 0 or maxUptake > 1.0:
        maxUptake = np.random.normal(loc=mean, scale=std, size=1)
    
    maxUptake_dict[prop] = maxUptake
    
    return [maxUptake, maxUptake_dict]


def get_SpeciesDispParams(system, DispParamsDict, prop, mean, std):
    
    sp_mean = mean
    sp_std = std
        
    sp_std = (randrange(1, 10**5)/10**5) / randrange(10, 100)    
    DispParamsDict[prop] = [sp_mean, sp_std]
    
    return [[sp_mean, sp_std], DispParamsDict]
    
 

def immigration(system, V, Qs2, Slist2, inds2, Tlist, ID, maxUptake_dict, maxUptake_list, DispParamsDict, mean, std, lgp, num_in, PropQ, xcoords, ycoords, zcoords, disp_mu, disp_std):

    """ Immigration """
    propagules = np.random.logseries(lgp, num_in)  # list of propagules
    for prop in propagules:
        #PropQ = randrange(100,1000)/1000
        Qs2.append(PropQ) # add individual's cell quota to cell quota list
        Slist2.append(prop) # add species identity to species list
        inds2.append(ID)
        Tlist.append(0)
        
        scale = 10**5
        x = (float(round(random.randrange(scale)))/scale)*V
        y = (float(round(random.randrange(scale)))/scale)*V
        z = (float(round(random.randrange(scale)))/scale)*V
        
        xcoords.append(x)
        ycoords.append(y)
        zcoords.append(z)
        
        if prop in maxUptake_dict:
            maxUptake = maxUptake_dict[prop]
        else: 
            maxUptake, maxUptake_dict = get_maxUptake(maxUptake_dict, prop, mean, std, system) # maxUptake_method as a 3rd arg 
    
        maxUptake_list.append(maxUptake)
        
        if prop in DispParamsDict:
            disp_params = DispParamsDict[prop]
        else: 
            disp_params, DispParamsDict = get_SpeciesDispParams(system, DispParamsDict, prop, disp_mu, disp_std) # maxUptake_method as a 3rd arg 
            
                    
        ID += 1
            
    return [Qs2, Slist2, inds2, Tlist, ID, maxUptake_dict, maxUptake_list, DispParamsDict, xcoords, ycoords, zcoords]

            
            
def growth(Qs2, Slist2, inds2, maxUptake_list, Tlist, TR, Kq, maint, xcoords, ycoords, zcoords):
    
    """ Growth """
    Glist = []
    
    for i, val in enumerate(Qs2):
        
        Q = Qs2[i] - maint # individual cell quota minus maint
        
        maxUptake = maxUptake_list[i]
        
        if Q < Kq:
            
            Qs2.pop(i)
            Slist2.pop(i)
            
            inds2.pop(i)
            xcoords.pop(i)
            
            ycoords.pop(i)
            zcoords.pop(i)
            
            maxUptake_list.pop(i)
            Tlist.pop(i)
            continue
            
        if TR > 0:        
            g = maxUptake*(1.0 - Kq/Q) # Droop form of g(s) with maintenance costs
            #if Q + g > 1.0:
            #    g = 1.0 - Q
            
            if TR >= g:
                Q += g  # increase cell quota by uptake
                TR -= g # decrease total resources 
                Glist.append(g) # list of growth rates 
                
            elif TR < g:
                Q += TR
                Glist.append(TR) # list of growth rates 
                TR = 0
                
            if TR < 0:
                print 'Error: total resources is in the negative!'
                sys.exit() # kill the simulation for debugging
        
            if Q >= Kq:
                Qs2[i] = Q
            
            elif Q < Kq:
                
                Qs2.pop(i)
                Slist2.pop(i)
                
                inds2.pop(i)
                xcoords.pop(i)
                
                ycoords.pop(i)
                zcoords.pop(i)
                
                maxUptake_list.pop(i)
                Tlist.pop(i)
                
    return [Qs2, Slist2, inds2, maxUptake_list, Tlist, TR, Glist, xcoords, ycoords, zcoords] 
    
    


def reproduction(Qs, Slist, inds, Tlist, maxUptake_list, Kq, ID, maint, xcoords, ycoords, zcoords):

    """ Reproduction """
    new = 0
    N = len(Qs)
    
    for n in range(N):
        i = randrange(len(Qs))
        Q = Qs[i]
    
        if Q > Kq*2:
            
            if Q > 1.0:
                x = [1]
            else:
                x = np.random.binomial(1, Q, 1)
            
            if x[0] == 1:
                
                Q = Q/2.0
                Qs[i] = Q
                
                Qs.append(Q)     # individuals produce cells with half Q
                Slist.append(Slist[i])
                
                inds.append(ID)
                maxUptake_list.append(maxUptake_list[i])
                
                Tlist.append(0)
                xcoords.append(xcoords[i]) 
                
                ycoords.append(ycoords[i])
                zcoords.append(zcoords[i])
                
                new += 1
                ID += 1
    
    if N > 0:
        gofS = new/float(N)
      
    return [Qs, Slist, inds, Tlist, maxUptake_list, Kq, ID, gofS, xcoords, ycoords, zcoords]          
                
    
   
   
def emigration(Qs2, Slist2, inds2, Tlist, maxUptake_list, r, V, xcoords, ycoords, zcoords, system):   
    
    """ Emigration """
    Qs3 = []
    Slist3 = []
    
    inds3 = [] 
    xcoords3 = []
    
    ycoords3 = []
    zcoords3 = []
    
    maxUptake_list3 = [] 
    Tlist3 = []
    
    p_out = r/V # per capita chance of flowing out
    N = len(Qs2)
    
    for i in range(N):
        
        x = np.random.binomial(1, p_out, 1)
        if x[0] == 0: 
            
            Qs3.append(Qs2[i]) # individual emigrates, i.e., washed out
            Slist3.append(Slist2[i]) # remove species identity from species list
            
            inds3.append(inds2[i])
            xcoords3.append(xcoords[i])
            
            ycoords3.append(ycoords[i])
            zcoords3.append(zcoords[i])
            
            maxUptake_list3.append(maxUptake_list[i])
            Tlist3.append(Tlist[i] + 1)
            
    return [Qs3, Slist3, inds3, Tlist3, maxUptake_list3, xcoords3, ycoords3, zcoords3]
   
   
     
def outflow(TR, r, V):
    
    TR -= (r/V) * TR # Total Resources decrease due to outflow
    if TR < 0:
        print 'Error: There is a resource debt: TR =',TR,' D =',r/V
        sys.exit()
    if V <= 0:
        print 'Error: There is no volume: V =',V
        sys.exit()
    
    return [TR, V]
    
    
    
def dispersal(V, Qs2, Slist2, inds2, Tlist, maxUptake_list, DispParamsDict, Kq, ID, xcoords, ycoords, zcoords, system):
    
    """ Envision the local environment as a 3-dimensional space"""
    
    N = len(inds2)
    
    for i in range(N): # Simulate over individuals
        
        spID = Slist2[i] # unique individual ID
        mean, std = DispParamsDict[spID]
        
        go = 'no'
        while go == 'no':
            
            if system != 'non-ideal':
                scale = 10**5
                xcoords[i] = (float(round(randrange(scale)))/scale)*V
                ycoords[i] = (float(round(randrange(scale)))/scale)*V
                zcoords[i] = (float(round(randrange(scale)))/scale)*V
                
                go = 'yes'
                
            else:
                
                x = np.random.normal(loc=mean, scale=std, size=1)
                
                if xcoords[i] + x > V:
                    if xcoords[i] - x < 0:
                        print 'out of bounds for x-coord'
                        sys.exit()
                    else:
                        xcoords[i] = xcoords[i] - x
                else:
                    xcoords[i] = xcoords[i] + x
                
                y = np.random.normal(loc=mean, scale=std, size=1)
                
                if ycoords[i] + y > V:
                    if ycoords[i] - y < 0:
                        print 'out of bounds for y-coord'
                        sys.exit()
                    else:
                        ycoords[i] = ycoords[i] - y
                else:
                    ycoords[i] = ycoords[i] + y
                
                z = np.random.normal(loc=mean, scale=std, size=1)
                
                if zcoords[i] + z > V:
                    if zcoords[i] - z < 0:
                        print 'out of bounds for z-coord'
                        sys.exit()
                    else:
                        zcoords[i] = zcoords[i] - z
                else:
                    zcoords[i] = zcoords[i] + z
                
                go = 'yes'        
                    
    return [Qs2, Slist2, inds2, Tlist, maxUptake_list, DispParamsDict, Kq, ID, xcoords, ycoords, zcoords]
        
        
       