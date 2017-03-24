# -*- coding: utf-8 -*-

def grow(iD):
    ''' increase overall size and decrease endogenous resources'''

    for k, v in iD.items():
        if v['st'] == 'a':
            iD[k]['sz'] += v['gr'] * v['sz']
            iD[k]['q'] -= v['gr'] * v['q']

    return iD
