# -*- coding: utf-8 -*-
from random import choice, shuffle

def consume(iD, rD, ps):
    ''' increase endogenous resources but not overall size '''
    h, l, r, u = ps

    keys = list(iD)
    shuffle(keys)
    for k in keys:

        if iD[k]['st'] == 'd': continue
        if len(list(rD)) == 0: return [iD, rD]

        c = choice(list(rD))
        e = iD[k]['ef'][rD[c]['t']] * iD[k]['q']

        if  rD[c]['v'] > e:
            rD[c]['v'] -= e
            iD[k]['q'] += e

        else:
            iD[k]['q'] += rD[c]['v']
            del rD[c]

    return [iD, rD]
