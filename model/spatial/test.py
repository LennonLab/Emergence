from __future__ import division
from random import randint
import numpy as np
import math
import sys


def SAR(Xs, Ys, Zs, SpIDs, w):
    
    boxes = [list([]) for _ in xrange(w**2)]
    
    index = 0
    for i, val in enumerate(SpIDs):
        x = int(round(Xs[i]))
        y = int(round(Ys[i]))

        index = int(round(x + (y * w)))

        if index > len(boxes) - 1:
            index = len(boxes) - 1
        elif index < 0:
            index = 0
            
        boxes[index].append(val)
    
    sar = [] 
    i = 1
    while i < len(boxes):
        nums = [randint(0, len(boxes)-1) for p in range(0, i)]
        q = []
        for num in nums:
            q.extend(boxes[num])
        
        sar.append(len(list(set(q))))
        i = i*2
        
    return sar


Xs = np.array([randint(0, 50) for p in range(0, 100)])
Ys = np.array(Xs)
Zs = np.array(Xs)

SpIDs = np.array([randint(21, 67) for p in range(0, 100)])

Xs = Xs.tolist()
Ys = Ys.tolist()
Zs = Zs.tolist()

w = max([max(Xs), max(Ys), max(Zs)])

sar = SAR(Xs, Ys, Zs, SpIDs, w)

print sar