# -*- coding: utf-8 -*-
import numpy as np

def maintenance(iD):
    for k, v in iD.items():

        m = v['mt']
        if v['st'] == 'd': m = v['mt'] * v['mf']

        iD[k]['q'] -= m * v['q']
        if iD[k]['q'] == 0:
            iD[k]['q'] += m * v['q']
            iD[k]['sz'] -= m * v['sz']

        if iD[k]['sz'] < m or v['q'] <= m or np.isnan(v['sz']): del iD[k]

    return iD
