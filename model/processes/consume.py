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

        iD[k]['q'] += min([rD[c]['v'], e])
        rD[c]['v'] -= min([rD[c]['v'], e])
        if rD[c]['v'] <= 0: del rD[c]

    return [iD, rD]
