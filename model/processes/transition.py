# -*- coding: utf-8 -*-
import numpy as np

def transition(iD):

    for k, v in iD.items():
        if v['st'] == 'a' and np.random.binomial(1, v['rp']) == 1:
            iD[k]['st'] = 'd'

        elif v['st'] == 'd' and np.random.binomial(1, v['rp']) == 1:
            iD[k]['q'] -= v['rp'] * v['q']
            iD[k]['st'] = 'a'

    return iD
