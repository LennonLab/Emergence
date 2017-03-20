# -*- coding: utf-8 -*-

def ind_disp(iD, ps):
    h, l, r, u = ps

    for k, v in iD.items():
        if v['st'] == 'a' and v['q'] > 0:

            iD[k]['q'] -= v['di'] * v['q']
            iD[k]['x'] -= v['di'] * u
            iD[k]['y'] -= v['di'] * u

            if iD[k]['q'] < 0: del iD[k]

    return iD
