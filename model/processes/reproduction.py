# -*- coding: utf-8 -*-
import numpy as np
import time
import copy

def reproduce(sD, iD, ps, p = 0):

    for k, v in iD.items():
        if v['st'] == 'd' or v['q'] <= 0 or np.isnan(v['sz']): continue

        if np.random.binomial(1, v['gr']) == 1:
            p += 1
            iD[k]['q'] = v['q']/2.0
            iD[k]['sz'] = v['sz']/2.0

            i = time.time()
            iD[i] = copy.copy(iD[k])

            if np.random.binomial(1, 0.01) == 1:
                iD[i]['sp'] = i
                sD[i] = copy.copy(sD[iD[k]['sp']])
                sD[iD[k]['sp']] = copy.copy(sD[v['sp']])

    return [sD, iD, p]