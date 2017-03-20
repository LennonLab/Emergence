# -*- coding: utf-8 -*-
from random import randint
import numpy as np
import time

def ResIn(rD, ps):
    h, l, r, u = ps

    for i in range(r):
        p = np.random.binomial(1, u)
        ID = time.time()
        if p == 1:
            rD[ID] = {'t' : randint(0, 2)}
            rD[ID]['v'] = np.random.uniform(1, 1000)
            rD[ID]['x'] = float(np.random.uniform(0, 0.1*h))
            rD[ID]['y'] = float(np.random.uniform(0, 0.1*l))

    return rD
