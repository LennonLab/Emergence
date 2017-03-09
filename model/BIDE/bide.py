# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import sys
from math import sqrt, pi


def get_closest(RIDs, RX, RY, RZ, Rtypes, coords):
    closest = 0
    if len(RIDs) == 0:
        return closest

    x1, y1, z1 = coords
    Try = min([20, len(RIDs)])
    minDist, ct = 10**10, 0

    while ct < Try:
        ct += 1
        j = randint(0, len(RIDs)-1)
        x = RX[j]
        y = RY[j]
        z = RZ[j]

        dist = sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)

        if dist < minDist:
            minDist = dist
            closest = RIDs[j]

        return closest


def decomposition(i, n):
    gn = choice([[1/3, 1/3, 1/3], [0.5, 0.3, 0.2], [0.6, 0.2, 0.2],
                [0.8, 0.1, 0.1], [0.99, 0.005, 0.005], [0.9, 0.05, 0.05]])
    return gn


def get_color(ID, colorD): # FUNCTION TO ASSIGN COLORS TO Sp_

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    colorD[ID] = color

    return colorD



def ResIn(RTypes, RVals, RX, RY, RZ, RID, RIDs, params, u1):

    w, h, l, seed, m, r, nN, rmax, gmax, maintmax, dmax, amp, freq, phase, rates, pmax, dormlim, smax = params

    for i in range(r):
        x = np.random.binomial(1, u1)

        if x == 1:
            rval = np.random.uniform(1, rmax)
            rtype = randint(0, nN-1)
            RTypes.append(rtype)
            RVals.append(rval)
            RIDs.append(RID)
            RID += 1

            RX.append(float(np.random.uniform(0, 0.9*h)))
            RY.append(float(np.random.uniform(0, 0.9*l)))
            RZ.append(float(np.random.uniform(0, 0.9*w)))
    return [RTypes, RVals, RX, RY, RZ, RIDs, RID]



def immigration(SpDict, IndDict, params, ID, ct, u1):

    w, h, l, seed, m, r, nN, rmax, g_max, m_max, d_max, amp, freq, phase, rates, p_max, dlim, s_max = params
    # SpDict = 'disp', 's_max', 's_min', 'mfd', 'rpf', 'grow', 'maint, 'dlim', 'eff'
    # IndDict: key = ID
    # IndDict: values: 'spID', 'size', 'q', 'x', 'y', 'z', 'state'

    sd = int(seed)
    if ct > 1:
        sd = 1

    for j in range(sd):

        if m == 0 and ct > 1: break
        if sd == 1 and np.random.binomial(1, u1) == 0: continue

        prop = np.random.randint(1, 1000)
        if prop not in SpDict:

            # species max & min body size
            smax = np.random.uniform(100, 100)
            smin = np.random.uniform(1, 1)
            SpDict[prop] = {'s_max' : smax}
            SpDict[prop]['s_min'] = smin

            # species growth rate
            SpDict[prop]['grow'] = np.random.uniform(g_max/100, g_max)

            # species active dispersal rate
            SpDict[prop]['disp'] = np.random.uniform(d_max/100, d_max)

            # species RPF factor
            SpDict[prop]['rpf'] = np.random.uniform(p_max/100, p_max)

            # species maintenance
            SpDict[prop]['maint'] = np.random.uniform(m_max/100, m_max)

            # dormancy limit
            SpDict[prop]['dlim'] = np.random.uniform(dlim/100, dlim)

            # species speciation rate
            SpDict[prop]['spec'] = np.random.uniform(0.01, 0.01)

            # species maintenance factor
            mfd = np.random.logseries(0.95, 1)[0]
            if mfd > 40: mfd = 40
            SpDict[prop]['mfd'] = mfd + 1

            # A set of specific growth rates for three major types of resources
            SpDict[prop]['eff'] = list(decomposition(1, nN))

        IndDict[ID] = {'spID' : prop}
        IndDict[ID]['x'] = float(np.random.uniform(0, 0.9*h))
        IndDict[ID]['y'] = float(np.random.uniform(0, 0.9*l))
        IndDict[ID]['z'] = float(np.random.uniform(0, 0.9*w))

        smin = SpDict[prop]['s_min']
        smax = SpDict[prop]['s_max']

        size = float(np.random.uniform(smin, smax))
        IndDict[ID]['size'] = size
        q = size - smin
        IndDict[ID]['q'] = q

        IndDict[ID]['state'] = 'a'
        ID += 1

    return [SpDict, IndDict, ID]




def ind_flow(SpDict, IndDict, h, l, w, u0):

    for key, value in IndDict.items():

        x = value['x']
        y = value['y']
        z = value['z']
        q = value['q']
        sp = value['spID']
        sz = value['size']
        d = SpDict[sp]['disp']

        x1, y1, z1 = 1, 1, 1
        #if value['state'] == 'a':
        #    x1 = np.random.binomial(1, 1-d)
        #    y1 = np.random.binomial(1, 1-d)
        #    z1 = np.random.binomial(1, 1-d)
        #    q -= d * sz * u0
        #    IndDict[key]['q'] = q

        x += u0*x1
        y += u0*y1
        z += u0*z1

        if x > h or y > l or z > w:
            del IndDict[key]

        else:
            IndDict[key]['x'] = x
            IndDict[key]['y'] = y
            IndDict[key]['z'] = z

    return IndDict



def ind_disp(SpDict, IndDict, h, l, w, u0):

    for key, value in IndDict.items():
        x = value['x']
        y = value['y']
        z = value['z']
        q = value['q']
        sp = value['spID']
        sz = value['size']
        d = SpDict[sp]['disp']

        x1, y1, z1 = 0, 0, 0
        if value['state'] == 'a':
            x1 = np.random.binomial(1, d)
            y1 = np.random.binomial(1, d)
            z1 = np.random.binomial(1, d)
            q -= d * sz * u0
            IndDict[key]['q'] = q

        x -= x1*d
        y -= y1*d
        z -= z1*d

        if x > h or y > l or z > w:
            del IndDict[key]

        else:
            IndDict[key]['x'] = x
            IndDict[key]['y'] = y
            IndDict[key]['z'] = z

    return IndDict



def search(SpDict, IndDict, h, l, w, u0, RTypes, RVals, RXs, RYs, RZs, RIDs):

    for key, value in IndDict.items():

        x1 = value['x']
        y1 = value['y']
        z1 = value['z']
        sp = value['spID']
        q = value['q']
        d = SpDict[sp]['disp']
        sz = value['size']

        if value['state'] == 'd': continue

        coords = [x1, y1, z1]
        if len(RIDs):
            closest = get_closest(RIDs, RXs, RYs, RZs, RTypes, coords)
            ri = RIDs.index(closest)
            x2 = RXs[ri]
            y2 = RYs[ri]
            z2 = RZs[ri]

        else:
            x2 = np.random.uniform(0, h)
            y2 = np.random.uniform(0, l)
            z2 = np.random.uniform(0, w)

        x = np.abs(x1 - x2)
        y = np.abs(y1 - y2)
        z = np.abs(z1 - z2)

        q -= d * sz
        IndDict[key]['q'] = q

        if x1 > x2:
            x1 -= np.random.uniform(0, d*x)
        elif x1 < x2:
            x1 += np.random.uniform(0, d*x)

        if y1 > y2:
            y1 -= np.random.uniform(0, d*y)
        elif y1 < y2:
            y1 += np.random.uniform(0, d*y)

        if z1 > z2:
            z1 -= np.random.uniform(0, d*z)
        elif z1 < z2:
            z1 += np.random.uniform(0, d*z)

        IndDict[key]['x'] = x1
        IndDict[key]['y'] = y1
        IndDict[key]['z'] = z1

    return IndDict




def res_flow(RTypes, RVals, RX, RY, RZ, RID, RIDs, params, u0):

    w, h, l, seed, m, r, nN, rmax, g_max, m_max, d_max, amp, freq, phase, rates, p_max, dlim, s_max = params

    RX = np.array(RX) + u0
    RY = np.array(RY) + u0
    RZ = np.array(RZ) + u0

    i1 = np.where(RX > h)[0].tolist()
    i2 = np.where(RY > l)[0].tolist()
    i3 = np.where(RZ > w)[0].tolist()
    index = np.array(list(set(i1 + i2 + i3)))

    RTypes = np.delete(RTypes, index).tolist()
    RX = np.delete(RX, index).tolist()
    RY = np.delete(RY, index).tolist()
    RZ = np.delete(RZ, index).tolist()
    RVals = np.delete(RVals, index).tolist()
    RIDs = np.delete(RIDs, index).tolist()

    return [RTypes, RVals, RX, RY, RZ, RIDs, RID]




def grow(SpDict, IndDict):

    for key, value in IndDict.items():

        if value['state'] == 'a':
            sp = value['spID']
            q = value['q']
            sz = value['size']
            g = SpDict[sp]['grow']

            sz += sz * g
            q -= q * g

            if q <= 0:
                del IndDict[key]
                continue

            IndDict[key]['size'] = sz
            IndDict[key]['q'] = q

    return IndDict



def maintenance(SpDict, IndDict):

    for key, value in IndDict.items():

        sp = value['spID']
        q = value['q']
        sz = value['size']
        smin = SpDict[sp]['s_min']
        maint = SpDict[sp]['maint']

        q -= maint
        if q <= 0 or sz < smin:
            del IndDict[key]
        else:
            IndDict[key]['q'] = q

    return IndDict



def transition(SpDict, IndDict):
    for key, value in IndDict.items():

        sp = value['spID']
        q = value['q']
        state = value['state']

        dlim = SpDict[sp]['dlim']
        mfd = SpDict[sp]['mfd']
        maint = SpDict[sp]['maint']
        rpf = SpDict[sp]['rpf']
        smax = SpDict[sp]['s_max']

        if state is 'a' and q <= dlim*smax:
            IndDict[key]['maint'] = maint/mfd
            IndDict[key]['state'] = 'd'

        elif state is 'd' and np.random.binomial(1, rpf) == 1:
            IndDict[key]['maint'] = maint*mfd
            IndDict[key]['q'] -= q * rpf
            IndDict[key]['state'] = 'a'

    return IndDict




def consume(SpDict, IndDict, h, l, w, u0, Rtypes, Rvals, RX, RY, RZ, RIDs):
    numc = 0

    for key, value in IndDict.items():
        sp = value['spID']
        state = value['state']
        x1 = value['x']
        y1 = value['y']
        z1 = value['z']
        sp = value['spID']
        q = value['q']
        sz = value['size']
        mfd = SpDict[sp]['mfd']
        eff = SpDict[sp]['eff']
        maint = SpDict[sp]['maint']
        rpf = SpDict[sp]['rpf']

        if q <= 0:
            del IndDict[key]
            continue

        coords = [x1, y1, z1]
        if len(RIDs):
            closest = get_closest(RIDs, RX, RY, RZ, Rtypes, coords)
            ri = RIDs.index(closest)
            x2 = RX[ri]
            y2 = RY[ri]
            z2 = RZ[ri]
            Rval = Rvals[ri]
            rtype = Rtypes[ri]

        else: continue

        dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        i_radius = ((0.75*sz)/pi)**(1.0/3)
        r_radius = ((0.75*Rval*1)/pi)**(1.0/3)

        if dist <= i_radius + r_radius:
            if state == 'd':
                x = np.random.binomial(1, rpf)
                if x == 0:
                    continue
                elif x == 1:
                    IndDict[key]['state'] == 'a'
                    IndDict[key]['maint'] = maint*mfd
                    IndDict[key]['q'] -= q * rpf
        else: continue

        numc += 1
        e = eff[rtype]
        if Rval > e * q: # Increase cell quota
            Rval -= e * q
            q += e * q
        else:
            q += Rval
            Rval = 0.0

        if Rval <= 0.0:
            Rvals.pop(ri)
            Rtypes.pop(ri)
            RIDs.pop(ri)
            RX.pop(ri)
            RY.pop(ri)
            RZ.pop(ri)
        else:
            Rvals[ri] = Rval

        IndDict[key]['q'] = q

    return [numc, SpDict, IndDict, h, l, w, u0, Rtypes, Rvals, RX, RY, RZ, RIDs]




def reproduce(u0, SpDict, IndDict, ID):

    for key, value in IndDict.items():
        state = value['state']
        sp = value['spID']
        q = value['q']
        smax = SpDict[sp]['s_max']
        sz = value['size']
        smin = SpDict[sp]['s_min']

        if state == 'd': continue

        if q < 1.0:
            del IndDict[key]
            continue

        if sz/2 < smin: continue
        if q/smax >= 1.0 or np.random.binomial(1, q/smax) == 1: # individual is large enough to reproduce
            x = value['x']
            y = value['y']
            z = value['z']


            f1 = np.random.uniform(0.5, 0.5)
            f2 = 1 - f1
            IndDict[key]['q'] = q*f1
            IndDict[ID] = {'q' : q*f2}
            IndDict[key]['size'] = sz*f1
            IndDict[ID]['size'] = sz*f2
            IndDict[ID]['state'] = 'a'
            IndDict[ID]['spID'] = sp

            mu = SpDict[sp]['spec']
            if mu > 0:
                p = np.random.binomial(1, mu)

                if p == 1: # speciate
                    new_sp = max(list(SpDict))+1
                    IndDict[ID]['spID'] = new_sp

                    # species growth rate
                    SpDict[new_sp] = {'grow' : SpDict[sp]['grow']}

                    # species active dispersal rate
                    SpDict[new_sp]['disp'] = SpDict[sp]['disp']

                    # species RPF factor
                    SpDict[new_sp]['rpf'] = SpDict[sp]['rpf']

                    # species maintenance
                    SpDict[new_sp]['maint'] = SpDict[sp]['maint']

                    # dormancy limit
                    SpDict[new_sp]['dlim'] = SpDict[sp]['dlim']

                    # species maintenance factor
                    SpDict[new_sp]['mfd'] = SpDict[sp]['mfd']

                    # species speciation rate
                    SpDict[new_sp]['spec'] = SpDict[sp]['spec']

                    SpDict[new_sp]['s_min'] = SpDict[sp]['s_min']
                    SpDict[new_sp]['s_max'] = SpDict[sp]['s_max']

                    # A set of specific growth rates for three major types of resources
                    SpDict[new_sp]['eff'] = SpDict[sp]['eff']

            IndDict[ID]['x'] = float(x)
            IndDict[ID]['y'] = float(y)
            IndDict[ID]['z'] = float(z)
            ID += 1

    return [SpDict, IndDict, ID]
