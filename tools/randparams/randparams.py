from __future__ import division
from random import choice
import numpy as np
import sys
import os

def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.py """

    motion = choice(['fluid', 'random_walk']) # 'fluid', 'unidirectional'
    D = choice([2, 2]) # number of spatial dimensions

    width = choice([10, 20, 40, 80, 160])
    height = choice([5, 10, 20, 40, 80])
    barriers = choice([0, 0])

    pulse = choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    flux = choice(['yes'])

    # Sine wave: y(t) = A * sin(2*pi*f*t + phi)
    # let phi = 0, meaning 0 amplitude at time 0
    amp = choice([0.05, 0.1, 0.2, 0.3, 0.4, 0.5]) # A
    freq = choice([0.1, 0.08, 0.06, 0.04, 0.02, 0.01]) # f
    phase = choice([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    # 0 = in phase; 16 = entirely out of phase

    disturb = choice([0.00001])#, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1])
    rates = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.02, 0.01, 0.0075, 0.005, 0.001, 0.0005, 0.0001])  # inflow speeds

    alpha = np.random.uniform(0.95, 0.999)
    reproduction = choice(['fission', 'sexual'])
    speciation = choice(['yes', 'no'])

    seedCom = 100 # size of starting community
    m = choice([0.0, 0.001, 0.005, 0.01, 0.05]) # m = probability of immigration
    r = choice([200, 400, 600, 800, 1000])
    r = choice([200]) #resource particles flowing in per time step
    nNi = choice([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) # max number of Nitrogen types
    nP = choice([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) # max number of Phosphorus types
    nC = choice([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) # max number of Carbon types

    envgrads = []
    num_envgrads = choice([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    for i in range(num_envgrads):

        x = np.random.uniform(0, width)
        y = np.random.uniform(0, height)
        envgrads.append([x, y])

    #rmax = choice([10000, 20000, 40000, 80000, 100000]) # maximum resource particle size
    rmax = choice([1000, 2000, 4000, 8000, 10000]) # maximum resource particle size

    gmax = choice([0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9])
    maintmax = choice([0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.002])
    maintmax = maintmax*gmax
    dmax = choice([0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5])

    # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS
    #m = 0
    speciation = 'yes'
    #maintmax = 0.0000001
    reproduction = 'fission'
    #width = 20
    #height = 5
    #num_envgrads = 10
    #barriers = 1
    #r = 100
    #gmax = 0.1
    #rmax = 1000
    #dmax = 0.0000000001

    return [width, height, alpha, motion, reproduction, speciation, \
            seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, amp, freq, \
            flux, pulse, phase, disturb, envgrads, barriers, rates]
