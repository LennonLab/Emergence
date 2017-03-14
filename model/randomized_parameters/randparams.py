from __future__ import division
from random import randint, choice
import numpy as np

def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.py """

    sd = 500 # size of starting community
    w = randint(10, 10)
    h = float(w)
    l = float(h)

    u = 10**np.random.uniform(-6, 0) # rate of flow
    n = randint(3, 3) # number of resource types

    r = randint(10, 10) # number of inflowing resource particles
    rm = np.random.uniform(100, 100) # max size of resource particles

    sm = np.random.uniform(100, 100) # max size of individuals

    a = np.random.uniform(0.1, 0.1) # amplitude
    f = np.random.uniform(0.1, 0.1) # frequency
    p = np.random.uniform(0.1, 0.1) # phase
    m = np.random.uniform(0., 0.) # immigration rate
    g = np.random.uniform(0.1, 0.1) # max growth rate
    dl = np.random.uniform(0.1, 0.1) # max dormancy limit
    dm = np.random.uniform(0.1, 0.1) # dispersal max
    pm = np.random.uniform(0.1, 0.1) # resuscitation max
    mm = np.random.uniform(0.1, 0.1) # maintenance max
    s = choice(['y', 'n']) # stochastic or not

    return [w, h, l, sd, m, r, n, rm, g, mm, dm, a, f, p, u, pm, dl, sm, s]
