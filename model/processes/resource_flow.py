# -*- coding: utf-8 -*-
def res_flow(rD, ps):
    h, l, r, u = ps

    for k, v in rD.items():
        rD[k]['x'] += u
        rD[k]['y'] += u

        if rD[k]['x'] > h or rD[k]['y'] > l or rD[k]['v'] <= 0: del rD[k]

    return rD
