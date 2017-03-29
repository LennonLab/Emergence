# -*- coding: utf-8 -*-
from __future__ import division
from random import choice, shuffle, randint, sample
import numpy as np
import time
import copy



def immigration(sD, iD, ps, sd=1):
    h, l, r, u = ps

    for j in range(sd):
        if sd == 1 and np.random.binomial(1, 0.01) == 0: continue

        p = np.random.randint(1, 10)
        if p not in sD:
            sD[p] = {'gr' : 10**np.random.uniform(-3, -1)} # growth rate
            sD[p]['di'] = 10**np.random.uniform(-2, 0) # active dispersal rate
            sD[p]['rp'] = 10**np.random.uniform(-2, 0) # RPF factor
            sD[p]['mt'] = 10**np.random.uniform(-3, -1) # maintenance
            sD[p]['mf'] = 10**np.random.uniform(-4, -1)
            es = np.random.uniform(1, 100, 3)
            sD[p]['ef'] = es/sum(es) # growth efficiencies
            r1 = lambda: randint(0,255)
            r2 = lambda: randint(0,255)
            r3 = lambda: randint(0,255)
            clr = '#%02X%02X%02X' % (r1(),r2(),r3())
            sD[p]['color'] = clr

        ID = time.time()
        iD[ID] = copy.copy(sD[p])
        iD[ID]['sp'] = p
        iD[ID]['x'] = np.random.uniform(0, h)
        iD[ID]['y'] = np.random.uniform(0, l)
        iD[ID]['sz'] = np.random.uniform(1, 100)
        iD[ID]['q'] = np.random.uniform(1, 100)
        iD[ID]['st'] = 'a'

    return [sD, iD]
    
    
    
def disturb(iD, ceiling):
    iDs = list(iD)
    iDs = sample(iDs, len(iDs) - int(ceiling/2))
    
    f_dict = {key: iD[key] for key in iD if key not in iDs}

    return f_dict
    
    
def ind_disp(iD, ps):
    h, l, r, u = ps

    for k, v in iD.items():
        if v['st'] == 'a' and v['q'] > 0:

            iD[k]['q'] -= v['di'] * v['q']
            iD[k]['x'] -= v['di'] * u
            iD[k]['y'] -= v['di'] * u

            if iD[k]['q'] < 0: del iD[k]

    return iD



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
    
    
    
def grow(iD):
    ''' increase overall size and decrease endogenous resources'''

    for k, v in iD.items():
        if v['st'] == 'a':
            iD[k]['sz'] += v['gr'] * v['sz']
            iD[k]['q'] -= v['gr'] * v['q']

    return iD
    


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
    
    

def ind_flow(iD, ps):
    h, l, r, u = ps

    for k, val in iD.items():
        iD[k]['x'] += u
        iD[k]['y'] += u

        if iD[k]['x'] > h or iD[k]['y'] > l: del iD[k]

    return iD
    
    
    
def reproduce(sD, iD, ps, p = 0):

    for k, v in iD.items():
        if v['st'] == 'd' or v['q'] <= 0 or np.isnan(v['sz']): continue

        if np.random.binomial(1, v['gr']) == 1:
            p += 1
            iD[k]['q'] = v['q']/2.0
            iD[k]['sz'] = v['sz']/2.0

            i = time.time()
            iD[i] = copy.copy(iD[k])

            if np.random.binomial(1, 0.0001) == 1:
                iD[i]['sp'] = i
                sD[i] = copy.copy(sD[iD[k]['sp']])
                sD[iD[k]['sp']] = copy.copy(sD[v['sp']])
            
            sD[iD[k]['x']] = iD[i]['x']
            sD[iD[k]['y']] = iD[i]['y']
            
    return [sD, iD, p]
    
    
    
def res_flow(rD, ps):
    h, l, r, u = ps

    for k, v in rD.items():
        rD[k]['x'] += u
        rD[k]['y'] += u

        if rD[k]['x'] > h or rD[k]['y'] > l or rD[k]['v'] <= 0: del rD[k]

    return rD
    
    

def ResIn(rD, ps):
    h, l, r, u = ps

    for i in range(r):
        p = np.random.binomial(1, u)
        ID = time.time()
        if p == 1:
            rD[ID] = {'t' : randint(0, 2)}
            rD[ID]['v'] = np.random.uniform(1, 10)
            rD[ID]['x'] = float(np.random.uniform(0, h))
            rD[ID]['y'] = float(np.random.uniform(0, l))

    return rD
    
    
    
def transition(iD):

    for k, v in iD.items():
        if v['st'] == 'a' and np.random.binomial(1, v['rp']) == 1:
            iD[k]['st'] = 'd'

        elif v['st'] == 'd' and np.random.binomial(1, v['rp']) == 1:
            iD[k]['q'] -= v['rp'] * v['q']
            iD[k]['st'] = 'a'

    return iD