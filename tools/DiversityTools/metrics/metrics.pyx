from __future__ import division
import sys
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import mpmath as mpm
from scipy.optimize import fsolve
import re
import random
import os
import math

mydir = os.path.expanduser("~/GitHub/rare-bio")
data = os.path.expanduser("~/data")

################################################################################
#### Section of evenness indices and descriptive statistical functions #########


############### RARITY #########################################################
def r_singletons(sad): # relative abundance of the least abundant taxon

    """ percent of species represented by a single individual """

    return 100 * sad.count(1)/len(sad)


def p_ZPtOne(sad):

    """ percent taxa with less than 0.1% N """
    N = sum(sad)
    S = len(sad)

    sad = np.array(sad)/N
    sad = sad*100
    sad = sad.tolist()

    numR = 0
    for sp in sad:
        if sp < 0.1:
            numR += 1

    return numR/S


def Rlogskew(sad):
    S = len(sad)

    if S <= 2.0:
        print 'S < 2, cannot compute log-skew'
        sys.exit()

    sad = np.log10(sad)
    mu = np.mean(sad)

    num = 0
    denom = 0
    for ni in sad:
        num += ((ni - mu)**3.0)/S
        denom += ((ni - mu)**2.0)/S

    t1 = num/(denom**(3.0/2.0))
    t2 = (S/(S - 2.0)) * np.sqrt((S - 1.0)/S)

    return t1 * t2


############### LOGNORMAL VARIABLES ############################################

def Preston(sad):
    N = sum(sad)
    Nmax = max(sad)

    left = (2 * N)/(np.sqrt(np.pi) * Nmax)

    func = lambda a : left - (math.erf(np.log(2)/a) / a)

    guess = 0.1 # alpha is often ~0.2, but appears to be lower for larger N
    a = fsolve(func, guess)

    expS = (np.sqrt(np.pi) / a) * np.exp( (np.log(2)/(2*a))**2 )

    return a[0], expS[0]



############### DOMINANCE ######################################################


def Berger_Parker(sad):
    return max(sad)/sum(sad)


def McNaughton(sad):
    sad.sort(reverse=True)
    return 100 * (sad[0] + sad[1])/sum(sad)



############ DIVERSITY #########################################################


def Shannons_H(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'

    H = 0
    for i in sad:
        p = i/sum(sad)
        H += p*np.log(p)
    return H*-1.0


def simpsons_dom(sad): # ALSO CONSIDERED A DOMINANCE MEASURE

    sad = filter(lambda a: a != 0, sad)

    D = 0.0
    N = sum(sad)

    for x in sad:
        D += x*x
    D = 1 - (D/(N*N))

    return D


######### EVENNESS #############################################################


def e_shannon(sad):

    sad = filter(lambda a: a != 0, sad)

    H = Shannons_H(sad)
    S = len(sad)
    return H/np.log(S)



def simplest_gini(x):
        """Return computed Gini coefficient of inequality.
        This function was found at http://econpy.googlecode.com/svn/trunk/pytrix/utilities.py """

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


def e_Mcintosh(SAD):

    SAD = filter(lambda a: a != 0, SAD)

    S = len(SAD)
    N = sum(SAD)
    sum_n = 0
    for n in SAD: sum_n += n**2
    U = np.sqrt(sum_n)
    E = (N - U)/(N - (N/np.sqrt(S)))
    return E


def e_pielou(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    S = len(SAD)
    N = float(sum(SAD))
    H = 0
    for p in SAD:
        H += -(p/N)*np.log(p/N)
    J = H/np.log(S)
    return J



def EQ(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    SAD.reverse()
    S = len(SAD)

    y_list = np.log(SAD).tolist()
    x_list = []
    for rank in range(1, S+1):
        x_list.append((rank)/S)

    slope, intercept, rval, pval, std_err = stats.linregress(x_list, y_list)

    Eq = 1 + (-2/np.pi) * np.arctan(slope)
    return Eq



def NHC(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    SAD.sort()
    SAD.reverse()
    x_list = range(1,len(SAD)+1)
    y_list = np.log(SAD)
    slope,intercept,r_value,p_value,std_err = stats.linregress(x_list, y_list)

    return slope


def e_heip(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    S = len(SAD)
    N = float(sum(SAD))
    H = 0.0
    for p in SAD:
        if p < 1.0:
            print 'p < 1.0', p
            sys.exit()

        H += -(p/N)*np.log(p/N)
    H = (np.exp(H) - 1)/(S - 1)
    return H



def e_simpson(SAD): # based on 1/D, not 1 - D
    SAD = filter(lambda a: a != 0, SAD)

    D = 0.0
    N = sum(SAD)
    S = len(SAD)

    for x in SAD:
        D += (x*x) / (N*N)

    E = round((1.0/D)/S, 4)

    if E < 0.0 or E > 1.0:
        print 'Simpsons Evenness =',E
    return E


def berger_parker(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    bp = max(SAD)/sum(SAD)
    return bp



def e_var(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    P = np.log(SAD)
    S = len(SAD)
    mean = np.mean(P)
    X = 0
    for x in P:
        X += (x - mean)**2/S
    evar = 1.0 - 2/np.pi*np.arctan(X)

    if evar < 0.0 or evar > 1.0:
        print 'Evar =',evar
    return evar



def OE(SAD):
    SAD = filter(lambda a: a != 0, SAD)

    S = len(SAD)
    N = sum(SAD)
    o = 0

    for ab in SAD:
        o += min(ab/N, 1/S)

    return o



def camargo(SAD): # function to calculate Camargo's eveness:
    SAD = filter(lambda a: a != 0, SAD)

    S = len(SAD)
    SAD = np.array(SAD)/sum(SAD)
    SAD = SAD.tolist()

    E = 1
    for i in range(0, S-1):
        for j in range(i+1, S-1):

            pi = SAD[i]
            pj = SAD[j]
            E -= abs(pi - pj)/S

    return E



######## RICHNESS ESTIMATORS ###################################################


def Margalef(sad):
    sad = filter(lambda a: a != 0, sad)
    return (len(sad) - 1)/np.log(sum(sad))

def Menhinick(sad):
    sad = filter(lambda a: a != 0, sad)
    return len(sad)/np.sqrt(sum(sad))


def EstimateS2(SiteList):

    """ Chao and ICE estimators of S for two or more samples. These metrics
    account for the occurrence (presence/absence) of taxa in a sample, but not
    the observed abundance """

    m = len(SiteList)
    m_inf = 0
    SpDict = {}

    for site in SiteList:

        if min(site) <= 10: m_inf += 1

        for sp in site:
            if sp in SpDict:
                SpDict[sp] += 1
            else: SpDict[sp] = 1

    IncVals = SpDict.values()
    S = len(IncVals)

    qs = [0]*10

    for i, q in enumerate(qs):

        qs[i] = IncVals.count(i+1)

    # Chao2
    q1 = qs[0]
    q2 = qs[1]
    chao2 = S + (((m-1)/m) * ((q1*(q1-1)) / (2*(q2+1))))

    var = 'und'
    if q1 > 0 and q2 > 0:
        var = q2 * (0.5*(q1/q2)**2 + (q1/q2)**3 + 0.25*(q1/q2)**4)

    # ICE
    num = 0
    n_inf = 0

    for i, qk in enumerate(qs):
        num += (i+1)*i*qk
        n_inf += (i+1)*qk

    ci = 1 - (q1/n_inf)

    gamma = (sum(qs)/ci) * (m_inf/(m_inf-1)) * (num/(n_inf**2)) - 1
    cv = max(0, gamma)

    ice = (S-sum(qs)) + (sum(qs)/ci) + ((q1/ci) * cv)

    return [chao2, var, ice, S]


def EstimateS1(sad):

    """ Chao and ACE estimators of S for a single sample. These metrics account
    for the abundance of taxa in a sample, but inherently include the assumption
    that all species co-occurred or, at least, do not preserve any spatial
    structure.

    This code is based on equations in [Magurran & McGill (2013). Biological
    Diversity: Frontiers in measurement and assessment.] and on equations from
    the EstimateS user guide: http://viceroy.eeb.uconn.edu/estimates/
    EstimateSPages/EstSUsersGuide/EstimateSUsersGuide.htm#AppendixB obtained on
    10 April 2015.
    """

    # Chao1 estimator
    s_obs = len(sad)
    n = sum(sad)
    m_inf = 0

    ones = sad.count(1)
    twos = sad.count(2)

    jknife1 = s_obs + ones
    jknife2 = s_obs + 2*ones - twos

    chao = s_obs + (n/(n-1)) * (ones*(ones-1))/(2*twos+1)

    # ACE estimator
    Srare = 0
    Sabund = 0
    Nrare = 0
    for ab in sad:

        if ab < 11:
            Srare += 1
            Nrare += ab
        elif ab >= 11:
            Sabund += 1

    if Nrare < 2:
        ace = chao # using chao when ace is undefined
        return [chao, ace, jknife1, jknife2]

    else:
        Cace = 1 - (ones/Nrare)

        num = 0
        for i in range(1, 11):
            num += i*(i-1)*sad.count(i)

        denom = Nrare*(Nrare-1)

        if Cace == 0 or denom == 0:
            ace = chao # using chao when ace is undefined
            return [chao, ace, jknife1, jknife2]

        Cvar = max((Srare/Cace)*(num/denom) - 1, 0)
        ace = Sabund + (Srare/Cace) + (ones/Cace)*Cvar

    return chao, ace, jknife1, jknife2





######## ISLAND OF MISFIT FUNCTION #############################################

def get_skews(sad):
    sad = filter(lambda a: a != 0, sad)

    skews = []
    for i in sad:
        skews.append(stats.skew(i))

    return skews


def get_modal(_list):

    """ Finds the mode from a kernel density function across a sample """
    exp_mode = 0.0
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    density.covariance_factor = lambda : .001
    density._compute_covariance()
    D = [xs,density(xs)]
    d = 0
    maxd = 0.0
    while d < len(D[1]):
        if D[1][d] > maxd:
            maxd = D[1][d]
            exp_mode = D[0][d]
        d += 1
    return exp_mode



def get_kdens_choose_kernel(xlist,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(xlist)
    n = len(xlist)
    xs = np.linspace(min(xlist),max(xlist),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D



def get_kdens(xlist):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(xlist)
    #xs = np.linspace(min(xlist),max(xlist),n)
    xs = np.linspace(0.0,1.0,len(xlist))
    density.covariance_factor = lambda : 0.5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


def WhittakersTurnover(site1, site2):

  if len(site1) == 0 and len(site2) == 0:
      return 0
  elif len(site1) == 0 or len(site2) == 0:
      return 1.0

  set1 = set(site1)
  set2 = set(site2)

  gamma = set1.union(set2)         # Gamma species pool
  s     = len(gamma)                                   # Gamma richness
  abar = np.mean([len(set1), len(set2)])   # Mean sample richness
  bw   = s/abar - 1
  return bw
