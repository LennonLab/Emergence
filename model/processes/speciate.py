# -*- coding: utf-8 -*-
import numpy as np


def speciate(sD, iD, ps):

    w, h, l, sd, m, r, n, rm, gm, mm, d, a, f, p, u, pm, dl, sm, st = ps

    Qs = []
    for key, value in iD.items():
        q = value['q']

        if q < 0: del iD[key]
        else: Qs.append(value['q'])

    b1 = sum(Qs)
    n1 = len(Qs)

    for key, value in iD.items():
        state = value['state']
        sp = value['spID']
        q = value['q']
        smax = sD[sp]['s_max']
        sz = value['size']
        smin = sD[sp]['s_min']

        if state == 'd' or sz/2 <= smin: continue

        if q > smax or np.random.binomial(1, q/smax) == 1: # individual is large enough to reproduce
            x = value['x']
            y = value['y']
            z = value['z']

            f1 = np.random.uniform(0.5, 0.5)
            iD[key]['q'] = q*f1
            iD[key]['size'] = sz*f1


            i = id(f1)
            iD[i] = {'q' : q*f1}
            iD[i]['size'] = sz*f1
            iD[i]['state'] = 'a'
            iD[i]['spID'] = sp

            mu = sD[sp]['spec']
            if mu > 0:
                p = np.random.binomial(1, mu)
                if p == 1: # speciate
                    new_sp = id(i)
                    iD[i]['spID'] = new_sp

                    # species growth rate
                    sD[new_sp] = {'grow' : sD[sp]['grow']}

                    # species active dispersal rate
                    sD[new_sp]['disp'] = sD[sp]['disp']

                    # species RPF factor
                    sD[new_sp]['rpf'] = sD[sp]['rpf']

                    # species maintenance
                    sD[new_sp]['maint'] = sD[sp]['maint']

                    # dormancy limit
                    sD[new_sp]['dlim'] = sD[sp]['dlim']

                    # species maintenance factor
                    sD[new_sp]['mfd'] = sD[sp]['mfd']

                    # species speciation rate
                    sD[new_sp]['spec'] = sD[sp]['spec']

                    sD[new_sp]['s_min'] = sD[sp]['s_min']
                    sD[new_sp]['s_max'] = sD[sp]['s_max']

                    # A set of specific growth rates for three major types of resources
                    sD[new_sp]['eff'] = sD[sp]['eff']

            iD[i]['x'] = float(x)
            iD[i]['y'] = float(y)
            iD[i]['z'] = float(z)

    Qs = []
    for key, value in iD.items(): Qs.append(value['q'])
    b2 = sum(Qs)
    n2 = len(Qs)

    bB = b2 - b1
    nN = n2 - n1

    return [sD, iD, bB, nN]
