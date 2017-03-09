from __future__ import division
#import sys
from random import choice, randint
import numpy as np
#from scipy import stats
#from scipy.stats import gaussian_kde
#from scipy.optimize import fsolve
import math


def distance(p0, p1):

    """ take two (x, y) tuples as parameters
    http://stackoverflow.com/questions/5407969/distance-formula-between-two-points-in-a-list"""

    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)



def nearest_neighbor(xlist1, xlist2, ylist1, ylist2, q=1):

    """ xlist1 and ylist1 are assumed to be the lists of individual organisms
    xlist2 and ylist2 are assumed to be the lists of whatever the individual
    organisms are being measured with respect to their distance to """

    n = len(xlist1)
    r = len(xlist2)
    r = min([10, r])

    refpoints = min([10, n])
    DistList = []

    for ref in range(refpoints):

        i = randint(0, n-1)
        x1 = xlist1[i]
        y1 = ylist1[i]
        MinDist = 10000

        for j in range(r):

            x2 = xlist2[j]
            y2 = ylist2[j]
            dist = distance((x1, y1), (x2, y2))

            if dist < MinDist:
                MinDist = dist

        DistList.append(MinDist)

    return np.mean(DistList)





def avg_dist(xlist1, xlist2, ylist1, ylist2, q=1):

    """ xlist1 and ylist1 are assumed to be the lists of individual organisms
    xlist2 and ylist2 are assumed to be the lists of whatever the individual
    organisms are being measured with respect to their distance to """

    nmax = len(xlist1)
    rmax = len(xlist2)

    refpoints = min([100, nmax])
    dist = []

    for n in range(refpoints):
        for j in range(q):

            i1 = choice(range(nmax))
            x1 = xlist1[i1]
            y1 = ylist1[i1]

            i2 = choice(range(rmax))
            x2 = xlist2[i2]
            y2 = ylist2[i2]

            dist.append(distance((x1, y1), (x2, y2)))

    return np.mean(dist)


"""
def nearest_neighbor_dist(xlist1, xlist2, ylist1, ylist2, q=1):

    nmax = len(xlist1)
    refpoints = min([100, nmax])

    for n in range(refpoints):
        for j in range(q):

            i = choice(range(nmax))
            x1 = xlist1[i]
            x2 = xlist2[i]
            y1 = ylist1[i]
            y2 = ylist2[i]


    return
"""
