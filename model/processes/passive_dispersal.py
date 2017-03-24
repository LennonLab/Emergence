# -*- coding: utf-8 -*-
def ind_flow(iD, ps):
    h, l, r, u = ps

    for k, val in iD.items():
        iD[k]['x'] += u
        iD[k]['y'] += u

        if iD[k]['x'] > h or iD[k]['y'] > l: del iD[k]

    return iD
