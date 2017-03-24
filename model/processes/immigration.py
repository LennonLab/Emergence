# -*- coding: utf-8 -*-
from __future__ import division
from random import choice
import numpy as np
import time
import copy

def immigration(sD, iD, ps, sd=1):
    h, l, r, u = ps

    for j in range(sd):
        if np.random.binomial(1, u) == 0: continue

        p = np.random.randint(1, 1000)
        if p not in sD:
            sD[p] = {'gr' : 10**np.random.uniform(-4, 0)} # growth rate
            sD[p]['di'] = 10**np.random.uniform(-4, 0) # active dispersal rate
            sD[p]['rp'] = 10**np.random.uniform(-4, 0) # RPF factor
            sD[p]['mt'] = 10**np.random.uniform(-4, 0) # maintenance
            sD[p]['mf'] = 10**np.random.uniform(-4, 0)
            es = np.random.uniform(1, 100, 3)
            sD[p]['ef'] = es/sum(es) # growth efficiencies

        ID = time.time()
        iD[ID] = copy.copy(sD[p])
        iD[ID]['sp'] = p
        iD[ID]['x'] = np.random.uniform(0, h)
        iD[ID]['y'] = np.random.uniform(0, l)
        iD[ID]['sz'] = np.random.uniform(1, 100)
        iD[ID]['q'] = np.random.uniform(1, 100)
        iD[ID]['st'] = choice(['a', 'a'])

    return [sD, iD]
