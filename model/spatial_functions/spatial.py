from __future__ import division
import  matplotlib.pyplot as plt
from random import choice, randint, shuffle
from scipy import stats
from math import isnan
import numpy as np
import math
import sys



def SARt1(X1s, Y1s, indC, SpID1s, h):

    ''' strictly nested '''
    
    newX, newY, newS = [], [], []
    
    Xs, Ys, SpIDs = list(X1s), list(Y1s), list(SpID1s)
    xh, xl, yh, yl = float(h), 0, float(h), 0
    
    shuffle(SpIDs)
    species = []
    areas = []
    while ((xh - xl) * (yh - yl)) > 1:
        for i, x in enumerate(Xs):
            y = Ys[i]
            if x >= xl and x <= xh and y >= yl and y <= yh:
                newX.append(x)
                newY.append(y)
                newS.append(SpIDs[i])
                      
        s = len(list(set(newS)))
        if s > 0:
            species.append(s)
            a = (xh - xl) * (yh - yl)
            areas.append(a)  
            
        xl += 1
        yl += 1
        yh -= 1
        xh -= 1
                
        Xs = list(newX)
        Ys = list(newY)
        SpIDs = list(newS)
        newX, newY, newS = [], [], []
        
    areas.reverse()
    areas = np.array(areas)
    species.reverse()
    species = np.array(species)
        
    m, b, r, p, s = stats.linregress(np.log10(areas), np.log10(species))
    
    '''
    print species
    print areas
    print m
    
    fig = plt.figure()
    fig.add_subplot(2, 2, 1)
    plt.plot(np.log10(areas), np.log10(species), label='z = '+str(round(m,2)))
    plt.plot(np.log10(areas), b + np.log10(areas)*m, 'r')
    plt.legend(loc='best', fontsize=10)
    fig.add_subplot(2, 2, 2)
    
    plt.scatter(X1s, Y1s, color=indC)
    plt.show()
    sys.exit()
    '''
          
    return m
        


def SARt2(Xs, Ys, indC, SpIDs, h):

    ''' random accumulation '''
    
    boxes = [list([]) for _ in xrange(h**2)]

    index = 0
    for i, val in enumerate(SpIDs):
        x = int(round(Xs[i]))
        y = int(round(Ys[i]))

        index = int(round(x + (y * h)))

        if index > len(boxes) - 1:
            index = len(boxes) - 1
        elif index < 0:
            index = 0

        boxes[index].append(val)

    q = []
    species = []
    areas = []
    boxes2 = list(boxes)
    while boxes2:
        i = randint(0, len(boxes2)-1)
        box = boxes2.pop(i)
        q.extend(box)
        s = len(list(set(q)))
        if s > 0: 
            species.append(s)
            areas.append(len(q))
       
    m, b, r, p, s = stats.linregress(np.log10(areas), np.log10(species))
        
    '''
    print m
    fig = plt.figure()
    fig.add_subplot(2, 2, 1)
    plt.plot(np.log10(areas), np.log10(species), label='z = '+str(round(m,2)))
    plt.plot(np.log10(areas), b + np.log10(areas)*m, 'r')
    plt.legend(loc='best', fontsize=10)
    fig.add_subplot(2, 2, 2)

    plt.scatter(Xs, Ys, color=indC)
    plt.show()
    sys.exit()
    '''
        
    return m
        
        
        
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
