from __future__ import division
import sys
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import re
import random
import os


def camargo(SAD): # function to calculate Camargo's eveness:
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


def EQ_evenness(SAD):

    SAD.reverse()
    S = len(SAD)

    y_list = np.log(SAD).tolist()
    x_list = []
    for rank in range(1, S+1):
        x_list.append((rank)/S)

    slope, intercept, rval, pval, std_err = stats.linregress(x_list, y_list)

    Eq = 1 + (-2/np.pi) * np.arctan(slope)
    return Eq


def simpsons_dom(SAD):
    D = 0.0
    N = sum(SAD)

    for x in SAD:
        D += x*(x-1)
    D = 1 - (D/(N*(N-1)))

    return D


def simpsons_evenness(SAD): # based on 1/D, not 1 - D
    D = 0.0
    N = sum(SAD)
    S = len(SAD)

    for x in SAD:
        D += (x*x) / (N*N)

    E = (1/D)/S

    return E


def e_var(SAD):
    P = np.log(SAD)
    S = len(SAD)
    X = 0
    for x in P:
        X += (x - np.mean(P))**2/S
    evar = 1 - 2/np.pi*np.arctan(X)
    return(evar)


def OE(SAD):
    S = len(SAD)
    N = sum(SAD)
    o = 0

    for ab in SAD:
        o += min(ab/N, 1/S)

    return o


def McNaughton(sad):
    sad.sort(reverse=True)
    return 100 * (sad[0] + sad[1])/sum(sad)



def Rlogskew(sad):
    S = len(sad)
    sad = np.log10(sad)
    mu = np.mean(sad)

    num = 0
    denom = 0
    for ni in sad:
        num += ((ni - mu)**3)/S
        denom += ((ni - mu)**2)/S

    t1 = num/(denom**(3/2))
    t2 = (S/(S - 2)) * np.sqrt((S - 1)/S)
    return t1 * t2



SAD1 = [10]*10
print SAD1
E = camargo(SAD1)
Eq = EQ_evenness(SAD1)
Esimp = simpsons_evenness(SAD1)
Evar = e_var(SAD1)
O = OE(SAD1)
print E,',', Eq, ',', Esimp,',',Evar,',',O,'\n'


SAD2 = [91] + [1]*9

print SAD2
E = camargo(SAD2)
Eq = EQ_evenness(SAD2)
Esimp = simpsons_evenness(SAD2)
Evar = e_var(SAD2)
O = OE(SAD2)
print E,',', Eq, ',', Esimp,',',Evar,',',O,'\n'

print McNaughton(SAD1), McNaughton(SAD2)

print "---------RARITY---------"
SAD1 = [100, 90, 80, 70, 60, 50, 40, 30, 30, 10, 1]
SAD2 = [700, 700, 750, 750, 75, 75, 55, 55, 30, 20, 50]
r1 = Rlogskew(SAD1)
r2 = Rlogskew(SAD2)
print r1, r2
