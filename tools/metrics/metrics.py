from __future__ import division
import  matplotlib.pyplot as plt
import sys
import numpy as np
from numpy import log10
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.optimize import fsolve
import math



def count_pts_within_radius(x, y, radius, logscale=0):
    """Count the number of points within a fixed radius in 2D space"""
    #TODO: see if we can improve performance using KDTree.query_ball_point
    #http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html
    #instead of doing the subset based on the circle
    raw_data = np.array([x, y])
    x = np.array(x)
    y = np.array(y)
    raw_data = raw_data.transpose()

    # Get unique data points by adding each pair of points to a set
    unique_points = set()
    for xval, yval in raw_data:
        unique_points.add((xval, yval))

    count_data = []
    for a, b in unique_points:
        if logscale == 1:
            num_neighbors = len(x[((log10(x) - log10(a)) ** 2 +
                                   (log10(y) - log10(b)) ** 2) <= log10(radius) ** 2])
        else:
            num_neighbors = len(x[((x - a) ** 2 + (y - b) ** 2) <= radius ** 2])
        count_data.append((a, b, num_neighbors))
    return count_data



def plot_color_by_pt_dens(x, y, radius, loglog=0, plot_obj=None):
    """Plot bivariate relationships with large n using color for point density

    Inputs:
    x & y -- variables to be plotted
    radius -- the linear distance within which to count points as neighbors
    loglog -- a flag to indicate the use of a loglog plot (loglog = 1)

    The color of each point in the plot is determined by the logarithm (base 10)
    of the number of points that occur with a given radius of the focal point,
    with hotter colors indicating more points. The number of neighboring points
    is determined in linear space regardless of whether a loglog plot is
    presented.

    """
    plot_data = count_pts_within_radius(x, y, radius, loglog)
    sorted_plot_data = np.array(sorted(plot_data, key=lambda point: point[2]))

    if plot_obj == None:
        plot_obj = plt.axes()

    if loglog == 1:
        plot_obj.set_xscale('log')
        plot_obj.set_yscale('log')
        plot_obj.scatter(sorted_plot_data[:, 0], sorted_plot_data[:, 1],
                         c = np.sqrt(sorted_plot_data[:, 2]), edgecolors='none')
        plot_obj.set_xlim(min(x) * 0.5, max(x) * 2)
        plot_obj.set_ylim(min(y) * 0.5, max(y) * 2)
    else:
        plot_obj.scatter(sorted_plot_data[:, 0], sorted_plot_data[:, 1],
                    c = log10(sorted_plot_data[:, 2]), edgecolors='none')
    return plot_obj

############### RARITY #########################################################
def percent_ones(sad):
    """ percent of species represented by a single individual """

    sad = filter(lambda a: a != 0, sad)

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    return 100 * sad.count(1)/len(sad)


def percent_pt_one(sad):
    """ percent taxa with less than 0.1% N """

    sad = filter(lambda a: a != 0, sad)

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    N = sum(sad)
    S = len(sad)

    sad = np.array(sad)/N
    sad = sad*100
    sad = sad.tolist()

    numR = 0
    for sp in sad:
        if sp < 0.1:
            numR += 1

    return 100 * numR/S


def Rlogskew(sad):

    sad = filter(lambda a: a != 0, sad)

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'


    S = len(sad)

    if S <= 2.0:
        return 'NaN'

    if max(sad) == min(sad):
        return 'NaN'

    sad = np.log10(sad)
    mu = np.mean(sad)

    num = 0
    denom = 0
    for ni in sad:
        num += ((ni - mu)**3.0)/S
        denom += ((ni - mu)**2.0)/S

    t1 = num/(denom**(3.0/2.0))
    t2 = (S/(S - 2.0)) * np.sqrt((S - 1.0)/S)

    return round(t1 * t2, 4)


############### LOGNORMAL VARIABLES ############################################

def Preston(sad):

    sad = filter(lambda a: a != 0, sad)

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    N = sum(sad)

    if N <= 0:
        return 'NaN'

    Nmax = max(sad)

    left = (2 * N)/(np.sqrt(np.pi) * Nmax)

    func = lambda a : left - (math.erf(np.log(2)/a) / a)

    guess = 0.1 # alpha is often ~0.2, but appears to be lower for larger N
    a = fsolve(func, guess)

    expS = (np.sqrt(np.pi) / a) * np.exp( (np.log(2)/(2*a))**2 )

    return a[0], expS[0]



############### DOMINANCE ######################################################


def Berger_Parker(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    return max(sad)/sum(sad)


def McNaughton(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    if len(sad) == 1:
        return 'NaN'

    sad.sort(reverse=True)
    return 100 * (sad[0] + sad[1])/sum(sad)



############ DIVERSITY #########################################################


def Shannons_H(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'

    H = 0
    for i in sad:
        p = i/sum(sad)
        H += p*np.log(p)
    return round(H*-1.0, 6)


def simpsons_dom(sad): # ALSO CONSIDERED A DOMINANCE MEASURE

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'


    D = 0.0
    N = sum(sad)

    for x in sad:
        D += x*x
    D = 1 - (D/(N*N))

    return D


######### EVENNESS #############################################################


def e_shannon(sad):

    sad = filter(lambda a: a != 0, sad)

    if len(sad) <= 1:
        return 'NaN'

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'


    H = Shannons_H(sad)
    S = len(sad)
    return round(H/np.log(S), 6)



def simplest_gini(sad):
    """Return computed Gini coefficient of inequality.
    This function was found at:
    http://econpy.googlecode.com/svn/trunk/pytrix/utilities.py """

    #note: follows basic formula
    #see: `calc_gini2`
    #contact: aisaac AT american.edu

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'

    sad = sorted(sad)  # increasing order
    n = len(sad)
    G = sum(xi * (i+1) for i,xi in enumerate(sad))
    G = 2.0*G/(n*sum(sad)) #2*B
    return round(G - 1 - (1.0/n), 6)



def gini_sample(sads):
    """ Compute Gini's coefficient for each macrostate in a random sample """

    Gs = []
    for sad in sads:
        G = simplest_gini(sad)
        Gs.append(G)
    return Gs


def e_Mcintosh(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'

    S = len(sad)
    N = sum(sad)
    sum_n = 0
    for n in sad: sum_n += n**2
    U = np.sqrt(sum_n)
    E = (N - U)/(N - (N/np.sqrt(S)))
    return round(E, 6)



def EQ(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'
    sad.reverse()
    S = len(sad)

    y_list = np.log(sad).tolist()
    x_list = []

    for rank in range(1, S+1):
        x_list.append((rank)/S)

    slope, intercept, rval, pval, std_err = stats.linregress(x_list, y_list)

    Eq = 1 + (-2/np.pi) * np.arctan(slope)
    return Eq



def NHC(sad):

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'

    sad.sort()
    sad.reverse()
    x_list = range(1,len(sad)+1)
    y_list = np.log(sad)
    slope,intercept,r_value,p_value,std_err = stats.linregress(x_list, y_list)

    return slope


def e_heip(sad):
    sad = filter(lambda a: a != 0, sad)

    if sum(sad) <= 0:
        return 'NaN'

    x = sum(1 for n in sad if n < 0)
    if x >= 1:
        return 'NaN'

    sad = filter(lambda a: a != 0, sad)

    if sum(sad) == 0:
        return 'NaN'

    S = len(sad)
    N = float(sum(sad))
    H = 0.0
    for p in sad:
        if p < 1.0:
            print 'p < 1.0', p
            sys.exit()

        H += -(p/N)*np.log(p/N)
    H = (np.exp(H) - 1)/(S - 1)
    return H



def e_simpson(sad): # based on 1/D, not 1 - D
    sad = filter(lambda a: a != 0, sad)

    D = 0.0
    N = sum(sad)
    S = len(sad)

    for x in sad:
        D += (x*x) / (N*N)

    E = round((1.0/D)/S, 4)

    if E < 0.0 or E > 1.0:
        print 'Simpsons Evenness =',E
    return E



def e_var(sad):
    sad = filter(lambda a: a != 0, sad)

    P = np.log(sad)
    S = len(sad)
    mean = np.mean(P)
    X = 0
    for x in P:
        X += (x - mean)**2/S
    evar = 1.0 - 2/np.pi*np.arctan(X)

    if evar < 0.0 or evar > 1.0:
        print 'Evar =',evar
    return evar



def OE(sad):
    sad = filter(lambda a: a != 0, sad)

    S = len(sad)
    N = sum(sad)
    o = 0

    for ab in sad:
        o += min(ab/N, 1/S)

    return o



def camargo(sad): # function to calculate Camargo's eveness:
    sad = filter(lambda a: a != 0, sad)

    S = len(sad)
    sad = np.array(sad)/sum(sad)
    sad = sad.tolist()

    E = 1
    for i in range(0, S-1):
        for j in range(i+1, S-1):

            pi = sad[i]
            pj = sad[j]
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
    #m_inf = 0

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
    """ Finds the kernel density function across a sample of sads """
    density = gaussian_kde(xlist)
    n = len(xlist)
    xs = np.linspace(min(xlist),max(xlist),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D



def get_kdens(xlist):
    """ Finds the kernel density function across a sample of sads """
    density = gaussian_kde(xlist)
    #xs = np.linspace(min(xlist),max(xlist),n)
    xs = np.linspace(0.0,1.0,len(xlist))
    density.covariance_factor = lambda : 0.5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


def jaccard(seq1, seq2):

    """ Obtained from: https://github.com/doukremt/distance/blob/master/distance/_simpledists.py
        on Sept 8 2015

    Compute the Jaccard distance between the two sequences `seq1` and `seq2`.
    They should contain hashable items.
    The return value is a float between 0 and 1, where 0 means equal, and 1 totally different.
    """

    set1, set2 = set(seq1), set(seq2)
    return 1 - len(set1 & set2) / float(len(set1 | set2))


def sorensen(seq1, seq2):

    if len(seq1) == 0 and len(seq2) == 0:
      return 0
    elif len(seq1) == 0 or len(seq2) == 0:
      return 1.0
    """ Obtained from: https://github.com/doukremt/distance/blob/master/distance/_simpledists.py
        on Sept 8 2015

    Compute the Sorensen distance between the two sequences `seq1` and `seq2`.
    They should contain hashable items.
    The return value is a float between 0 and 1, where 0 means equal, and 1 totally different.
    """

    set1, set2 = set(seq1), set(seq2)
    return 1 - (2 * len(set1 & set2) / float(len(set1) + len(set2)))



def WhittakersTurnover(site1, site2):

  """ citation: """
  if len(site1) == 0 and len(site2) == 0:
      return 0
  elif len(site1) == 0 or len(site2) == 0:
      return 1.0

  set1 = set(site1)
  set2 = set(site2)

  gamma = set1.intersection(set2)         # Gamma species pool
  s = len(gamma)                                   # Gamma richness
  abar = np.mean([len(set1), len(set2)])   # Mean sample richness
  bw   = ((len(set1) - s) + (len(set2) - s))/abar

  return bw



def getprod(Qs):

    N, P, C = 0, 0, 0

    if len(Qs) == 0:
        return [0, N, P, C]

    if len(Qs) > 0:
        for q in Qs:
            N += q[0]
            P += q[1]
            C += q[2]

        p1 = len(Qs)

        return [p1, N, P, C]
