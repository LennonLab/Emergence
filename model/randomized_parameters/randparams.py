from __future__ import division
from random import randint
import numpy as np
import sys

def get_rand_params(width=0):
    """ Get random model parameter values. Others are chosen in bide.py """

    seedCom = 1000 # size of starting community

    if width == 0:
        width = 4
    elif width < 10:
        width += 1
    elif width >= 10:
        width = 4

    height = float(width)
    length = float(height)

    low = np.random.uniform(4, 4)
    rates = np.linspace(0, -low, 40)
    rates = 10**rates
    rates = rates.tolist()

    nN = randint(3, 3)
    amp = np.random.uniform(10**-1, 10**-1)
    freq = np.random.uniform(10**-1, 10**-1)
    phase = np.random.uniform(10**-1, 10**-1)
    m = np.random.uniform(0.1, 0.1)

    r = randint(20, 20)
    rmax = np.random.uniform(1000, 1000)

    dormlim = np.random.uniform(0.1, 0.1)
    gmax = np.random.uniform(0.2, 0.2)
    dmax = np.random.uniform(0.2, 0.2)
    pmax = np.random.uniform(0.1, 0.1)
    mmax = np.random.uniform(0.1, 0.1)
    smax = np.random.uniform(100, 100)

    return [width, height, length, seedCom, m, r, nN, rmax, gmax, mmax, dmax, amp, freq, phase, rates, pmax, dormlim, smax]
