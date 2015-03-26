from __future__ import division
import sys
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import re

########################################################################################################
######   A Section devoted to evenness indices and descriptive statistical functions ###################
minS = 7


def Bray_Curtis(ABV1, ABV2):
    
    """ A function that calculates the Bray-Curtis dissimilarity index, coded
    by Ken Locey.

    Because similarity is more intuitive than dissimilarity, this function
    actually return % similarity.

    A value of 0 means that the two sets being compared have no elements
    in common. A value of 1 means the two sets have identical members.
        
    ABV1 & ABV2  :  Lists that contain species id#s and the species abundance
                    e.g. ABV1 = [[1, 251], [2, 187], [3, 122], ...]
                    where the first element in each sublist is the species
                    id# and the second element is the abundance.
                    
    This function expects the lists to come sorted smallest to largest
    according to species id#  
    """
    list1, list2 = getlists(ABV1, ABV2)
    
    C = float()
    for i, ab in enumerate(list1):
        C += np.abs(ab - list2[i])
    
    C = 100 * (1 - C/(sum(list1) + sum(list2)))
                               
    return C
    

def Berger_Parker(sad):
    return max(sad)/sum(sad)

def Singletons(sad):
    singletons = sad.count(1)
    return 100*(singletons/len(sad))


def e_var(SAD): # Smith and Wilson's evenness metric (Evar) 
    P = np.log(RAD)
    S = len(RAD)
    X = 0
    for x in P:
        X += (x - np.mean(P))**2/S
    evar = 1 - 2/math.pi*np.arctan(X) 
    return(evar)


def Shannons_H(sad):
    H = 0
    for i in sad:
        p = i/sum(sad)
        H += p*np.log(p)
    return H*-1.0
    
            
def Shannons_even(sad):
    H = Shannons_H(sad)
    S = len(sad)
    return H/np.log(S)
    


def simplest_gini(x):
        """Return computed Gini coefficient of inequality. This function was found at http://econpy.googlecode.com/svn/trunk/pytrix/utilities.py """

        #note: follows basic formula
        #see: `calc_gini2`
        #contact: aisaac AT american.edu
        	
        x = sorted(x)  # increasing order
        n = len(x)
        G = sum(xi * (i+1) for i,xi in enumerate(x))
        G = 2.0*G/(n*sum(x)) #2*B
        return G - 1 - (1./n)

def gini_sample(SADs):
    """ Compute Gini's coefficient for each macrostate in a random sample """
    Gs = []
    for sad in SADs:
        G = simplest_gini(sad)
        Gs.append(G)
    return Gs


def Mcintosh_evenness(SAD):
    S = len(SAD)
    N = sum(SAD)
    sum_n = 0
    for n in SAD: sum_n += n**2
    U = np.sqrt(sum_n)    
    E = (N - U)/(N - (N/np.sqrt(S)))
    return E


def pielous_evenness(SAD):
    S = len(SAD)
    N = float(sum(SAD))
    H = 0
    for p in SAD:
        H += -(p/N)*np.log(p/N)
    J = H/np.log(S)
    return J


def NHC_evenness(SAD):
    SAD.sort()
    SAD.reverse()
    x_list = range(1,len(SAD)+1)
    y_list = np.log(SAD)
    slope,intercept,r_value,p_value,std_err = stats.linregress(x_list, y_list)
    
    if slope > 0.0:
        evar = e_var(SAD)
        print slope, p_value, evar
    return slope


def Heips_evenness(SAD):
    S = len(SAD)
    N = float(sum(SAD))
    H = 0.0
    for p in SAD:
        H += -(p/N)*np.log(p/N) 
    H = (np.exp(H) - 1)/(S - 1)
    return H


def simpsons_dom(SAD):
    D = 0.0
    N = sum(SAD)
    S = len(SAD)
    
    for x in SAD:
        D += x*(x-1)
    D = 1 - (D/(N*(N-1)))
    
    return D
    

def simpsons_evenness(SAD):
    D = 0.0
    N = sum(SAD)
    S = len(SAD)
    
    for x in SAD:
        D += (x*x) / (N*N)
    
    E = (1/D)/S
    if E > 1.0:
        print 'Simpsons violation',E
        print N,S, SAD
        sys.exit()
    
    return E


def berger_parker(SAD):
    bp = float(max(SAD))/sum(SAD)
    return bp

    
def EQ_evenness(SAD):
    
    SAD.sort()
    SAD.reverse()
    
    S = len(SAD)
    y_list = list(np.log(SAD))
    x_list = []
    for rank in range(1,S+1):
        x_list.append(rank/S)
    slope, intercept, rval, pval, std_err = stats.linregress(x_list, y_list)
    
    Eq = -2/np.pi*np.arctan(slope)
    return Eq


