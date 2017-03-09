from __future__ import division
import os, collections, math
from scipy import stats, optimize
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import nbinom
import numpy as np
import pandas as pd
from macroeco_distributions import pln, pln_solver, negbin_solver, trunc_geom
from scipy.optimize import fsolve
import scipy.optimize as opt
from numpy import log, log2, exp, sqrt, log10
from math import erf, pi
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity


"""This code was written using MIT liscenced code from the following Weecology
repos: METE (https://github.com/weecology/METE) and macroecotools
(https://github.com/weecology/macroecotools). """



class predictS:
    def __init__(self, N, Nmax, predictNmax = True):
        self.N = N
        self.Nmax = Nmax
        self.predictNmax = predictNmax

    def alpha(self, a, Nmax, Nmin=1):

        """Numerically solve for Preston's a. Needed to estimate S using the lognormal"""

        y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
        y = y * exp((log(2.0)/(2.0*a))**2.0)
        y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a))
        y += erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
        y -= self.N

        return y # find alpha

    def s(self, a, Nmax, Nmin=1):

        """Predict S from the lognormal using Nmax as the only empirical input"""

        return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10


    def getNmax(self, b=0.6148, slope=0.942904468437):

        """Predict Nmax using N and the scaling law of Nmax with N predicted by the lognormal"""

        #NmaxCalc = 10 ** (b + slope*(log10(self.N)))
        NmaxCalc = (self.N **  slope) * b
        return int(round(NmaxCalc))


    def getS(self, predictNmax=True):

        guess = 0.1 # initial guess for Pre ston's alpha
        Nmin = 1

        if self.predictNmax == True:
            Nmax = self.getNmax()
        else:
            Nmax = self.Nmax

        a = opt.fsolve(self.alpha, guess, (Nmax, Nmin))[0]
        S2 = self.s(a, Nmax, 1)

        return int(round(S2))

class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

class zipf:

    """ A class to obtain a zipf object with inherited mle shape parameter,
    mle form of the rank-abundance distribution, and a rank-abundance curve
    based on fitting the zipf to the observed via a generalized linear model."""

    def __init__(self, obs, estimator):
        self.obs = obs
        self.estimator = estimator

    def zipf_ll(self, ab, a):
        """Log-likelihood of the Zipf distribution with x_min = 1."""
        return sum(stats.zipf.logpmf(ab, a))

    def zipf_solver(self, ab):
        #ab = self.obs
        """Obtain the MLE parameter for a Zipf distribution with x_min = 1."""
        par0 = 1 + len(ab) / (sum(np.log(2 * np.array(ab))))
        def zipf_func(x):
            return -self.zipf_ll(ab, x)
        #par = optimize.fmin(zipf_func, x0 = par0, disp=False)[0]
        estimator = str(self.estimator)
        par = getattr(optimize, estimator)(zipf_func, x0 = par0, disp=False)[0]
        return par

    def from_cdf(self):
        """ Obtain the maximum likelihood form of the Zipf distribution, given
        the mle value for the Zipf shape parameter (a). Using a, this code
        generates a rank-abundance distribution (RAD) from the cumulative
        density function (cdf) using the percent point function (ppf) also known
        as the quantile function.
        see: http://www.esapubs.org/archive/ecol/E093/155/appendix-B.htm
        This is an actual form of the Zipf distribution, obtained from getting
        the mle for the shape parameter.
        """
        p = self.zipf_solver(self.obs)
        S = len(self.obs)
        rv = stats.zipf(a=p)
        rad = []
        for i in range(1, S+1):
            val = (S - i + 0.5)/S
            x = rv.ppf(val)
            rad.append(int(x))
        point = collections.namedtuple('Rad_and_p', ['x', 'y'])
        point_return = point(rad, y = p)
        return point_return

    def from_glm(self):

        """ Fit the Zipf distribution to the observed vector of integer values
        using a generalized linear model.
        Note: This is a fitted curve; not an actual form of the Zipf distribution
        This method was inspired by the vegan
        package's open source code on vegan's public GitHub repository:
        https://github.com/vegandevs/vegan/blob/master/R/rad.zipf.R
        on Thursday, 19 Marth 2015 """

        ranks = np.log(range(1, len(self.obs)+1))
        off = [np.log(sum(self.obs))] * len(self.obs)

        d = pd.DataFrame({'ranks': ranks, 'off': off, 'x':self.obs})

        lm = smf.glm(formula='x ~ ranks', data = d, family = sm.families.Poisson()).fit()
        pred = lm.predict()

        return pred

    def get_pred_iterative(self, pmf, S):
        """Function to get predicted abundances (reverse-sorted) for distributions with no analytical ppf."""

        cdf = [(S - i + 0.5) / S for i in range(1, S + 1)]
        cdf = np.sort(cdf)
        rad = []

        j = 0
        i = 1
        cdf_cum = 0

        # print statements below were used as checks in writing the code
        while j < len(cdf):
            print cdf_cum, len(pmf), i
            cdf_cum += pmf[i-1]

            while cdf_cum >= cdf[j]:
                print cdf_cum, cdf[j], rad.count(0)

                rad.append(i)
                #print 'N:', N, 'S:', S, 'species:', len(rad), 'abundance:',i, 'cum_abs:', sum(rad)#, 'kmin:', min(rad), 'kmax:', max(rad)

                j += 1
                if j == len(cdf):
                    rad.reverse()
                    return np.array(rad)
            i += 1


    def get_emp_cdf(self, pmf):
        """Compute the empirical cdf given a list or an array"""
        pmf = np.array(pmf)
        cdf = []

        for point in pmf:
            point_cdf = len(pmf[pmf <= point]) / len(pmf)
            cdf.append(point_cdf)
        return np.array(cdf)



    """ Note: In their paper, the authors M is our N. Their N is our S. Their kmax is our Nmax."""

    def zipf_rgf_params(self, obs_rad):

        N = sum(obs_rad)
        S = len(obs_rad)

        kmin = min(obs_rad)
        Nmax = max(obs_rad)
        avg_ab = N/S

        #Set gamma and upper bound of b.
        gamma = 1.99999
        b = 1

        _sum = 0
        for k in range(kmin, N):
            pk = math.exp(-b*k)/k**gamma
            if pk > 0:
                _sum += pk
            else:
                break

        A = 1/_sum

        #Find the best b.
        Nmaxtmp = N # initialize Nmaxtmp as N
        b0 = 2*b
        b1 = 0

        while(abs(Nmaxtmp - Nmax) > 1):
            b = (b0 + b1)/2
            sum1 = 0
            kc = 0

            for k in range(N, kmin,-1):
                sum1 += A * math.exp(-b*k)/k**gamma
                if sum1 > 1/S:
                    kc = k
                    break

            sum1 = 0
            sum2 = 0

            for k in range(kc, N):
                s1=k*math.exp(-b*k)/k**gamma
                s2=math.exp(-b*k)/k**gamma

                if s1 <= 0 or s2 <= 0:
                    break

                sum1 += s1
                sum2 += s2

            Nmaxtmp = sum1/sum2

            print Nmaxtmp, Nmax,'b =', b

            if Nmaxtmp > Nmax:
                b1 = b

            else:
                b0 = b

        sum1 = 0
        sum2 = 0

        for k in range(kmin, N):
            sum1 += math.exp(-b*k)/k**gamma
            sum2 += k * math.exp(-b*k)/k**gamma

        A = 1/sum1
        kavm = sum2/sum1

        print gamma,'\t', b,'\t',A,'\t',kavm,'\t',avg_ab  #gamma, b, A, modeling
        #Compare modelling <k> and real <k>. If they are different, guess another gamma (increasing gamma will decrease <k>).

        return [gamma, b, A, N]




    def zipf_pmf(self, gamma, b, A, N):

        pmf = []
        k = 1
        while k <= N:
            pK = A * np.exp(-b*k) / (k**gamma)
            pmf.append(pK)
            k += 1

        return pmf



    def zipf_rgf(self):

        S = len(self.obs)
        gamma, b, A, N = self.zipf_rgf_params(self.obs)
        pmf = self.zipf_pmf(gamma, b, A, N) # pmf for non-built-in function

        rad = self.get_pred_iterative(pmf, S)

        return (rad, gamma)

class lognorm:

    def __init__(self, obs, dist):
        self.obs = obs
        self.dist = dist

    def ppoints(self, n):
        """ numpy analogue or `R`'s `ppoints` function
            see details at http://stat.ethz.ch/R-manual/R-patched/library/stats/html/ppoints.html
            :param n: array type or number
            Obtained from: http://stackoverflow.com/questions/20292216/imitating-ppoints-r-function-in-python
            on 5 April 2016
            """
        if n < 10:
            a = 3/8
        else:
            a = 1/2

        try:
            n = np.float(len(n))
        except TypeError:
            n = np.float(n)
        return (np.arange(n) + 1 - a)/(n + 1 - 2*a)


    def lognorm_glm(self):

        """ Fit the lognormal distribution to the observed vector of integer values
        using a generalized linear model.
        Note: This is a fitted curve; not an actual form of the lognormal distribution
        This method was inspired by the vegan package's open source code on vegan's public
        GitHub repository: https://github.com/vegandevs/vegan/R/rad.lognormal.R
        on Thursday, 5 April 2016
        """

        ranks = np.log(range(1, len(self.obs)+1))
        ranks = -norm.ppf(self.ppoints(len(ranks)))

        d = pd.DataFrame({'rnks': ranks, 'x': self.obs})
        lm = smf.glm(formula='x ~ rnks', data = d, family = sm.genmod.families.family.Poisson(link=sm.genmod.families.links.log)).fit()
        pred = lm.predict()

        return pred


    def get_rad_pln(self, S, mu, sigma, lower_trunc = True):
        """Obtain the predicted RAD from a Poisson lognormal distribution"""
        abundance = list(np.empty([S]))
        rank = range(1, int(S) + 1)
        cdf_obs = [(rank[i]-0.5) / S for i in range(0, int(S))]
        j = 0
        cdf_cum = 0
        i = 1
        while j < S:
            cdf_cum += pln.pmf(i, mu, sigma, lower_trunc)
            while cdf_cum >= cdf_obs[j]:
                abundance[j] = i
                j += 1
                if j == S:
                    abundance.reverse()
                    return abundance
            i += 1

    def get_rad_negbin(self, S, n, p):
        """Obtain the predicted RAD from a negative binomial distribution"""
        abundance = list(np.empty([S]))
        rank = range(1, int(S) + 1)
        cdf_obs = [(rank[i]-0.5) / S for i in range(0, int(S))]
        j = 0
        cdf_cum = 0
        i = 1
        while j < S:
            cdf_cum += nbinom.pmf(i, n, p) / (1 - nbinom.pmf(0, n, p))
            while cdf_cum >= cdf_obs[j]:
                abundance[j] = i
                j += 1
                if j == S:
                    abundance.reverse()
                    return abundance
            i += 1

    def get_rad_from_obs(self):
        if self.dist == 'negbin':
            n, p = negbin_solver(self.obs)
            pred_rad = self.get_rad_negbin(len(self.obs), n, p)
        elif self.dist == 'pln':
            mu, sigma = pln_solver(self.obs)
            pred_rad = self.get_rad_pln(len(self.obs), mu, sigma)
        return (pred_rad, mu, sigma)

def get_Geom(N,S,zeros):

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = trunc_geom.ppf(np.array(cdf), SNratio, N)
    return abd

def e_simpson(SAD): # based on 1/D, not 1 - D

    " Simpson's evenness "
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

def skewness(RAD):
    skew = stats.skew(RAD)
    # log-modulo skewness
    lms = np.log10(np.abs(float(skew)) + 1)
    if skew < 0:
        lms = lms * -1
    return lms

def CV_KDE(oneD_array, select_bandwidth = True):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple