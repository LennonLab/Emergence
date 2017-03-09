from __future__ import division
from random import choice, randint
import numpy as np
import sys
import os

def get_rand_params(fixed):
    """ Get random model parameter values. Others are chosen in bide.py """

    envgrads = []
    seedCom = 100 # size of starting community
    rates = []

    if fixed is True:

        rates = np.array([1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001])

        #motion = 'white_noise'
        #motion = 'brown_noise'
        motion = 'fluid'

        #width = choice([5, 20, 5, 20])
        width = 10
        #height = int(width)
        height = choice([5])

        num_envgrads = 2
        for i in range(num_envgrads):
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([x, y])

        nNi = 3 # max number of Nitrogen types
        nP = 3 # max number of Phosphorus types
        nC = 3 # max number of Carbon types

        amp = 0.001
        freq = 0.001
        phase = 0.0
        pulse = 0.001
        flux = 'yes'

        disturb = 0.00001
        m = 0.01
        speciation = 0.01
        maintmax = 0.001

        reproduction = 'fission'
        alpha = 0.99
        barriers = 0

        #r = choice([1, 10, 100])#resource particles flowing in per time step
        r = 100
        gmax = 0.9
        rmax = 1000
        dmax = 0.1
        #p_of_going_active = np.random.uniform(0.0001, 0.01)
        #maint_factor = np.random.uniform(1, 100)


    elif fixed is False:

        #motion = choice(['fluid', 'white_noise']) # 'fluid', 'unidirectional'
        motion = 'fluid'
        #motion = 'white_noise'

        if motion == 'white_noise':
            rates = np.array([0.00001])
        else:
            rates = [choice([1.0, 0.8, 0.6, 0.4, 0.2, 0.08, 0.06, 0.04, 0.02, 0.008, 0.006, 0.002])]

            #rates = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.2, 0.08, 0.06, 0.04, 0.02, 0.008, 0.006, 0.002])
            #rates = np.array([1.0, 0.7, 0.4, 0.1, 0.07, 0.04, 0.01, 0.007])
            #rates = np.array([1.0, 0.1, 0.01, 0.001])

            #rate = 1.001
            #rates = np.array([rate])


        width = randint(5, 10, 20, 30, 40, 50, 60, 80, 70, 90, 100)
        #width = 5

        #height = randint(10, 100)
        height = 20

        barriers = randint(0, 3)
        #barriers = 2

        pulse = np.random.uniform(0.01, 1.0)
        flux = choice(['yes'])

        # Sine wave: y(t) = amplitude * sin(2 * pi * frequency * t + phase)
        # if phi = 0, then there will be 0 amplitude at time 0
        amp = np.random.uniform(0.05, 0.5) # A
        freq = np.random.uniform(0.01, 0.1) # f
        phase = randint(0, 16) # 0 = in phase; 16 = entirely out of phase

        disturb = np.random.uniform(0.001, 0.1)
        alpha = np.random.uniform(0.95, 0.99)
        reproduction = choice(['fission'])

        speciation = np.random.uniform(0.01, 0.1)
        m = np.random.uniform(0.01, 0.1) # m = probability of immigration

        r = randint(10, 100) #resource particles flowing in per time step
        rmax = randint(100, 1000) # maximum resource particle size

        nNi = randint(1, 10) # max number of Nitrogen types
        nP = randint(1, 10) # max number of Phosphorus types
        nC = randint(1, 10) # max number of Carbon types

        num_envgrads = randint(1, 10)
        for i in range(num_envgrads):
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([x, y])

        #tp_max = np.random.uniform(0.001, 0.1) # maximum probability of transitioning at random from dormant to active (i.e., Scout hypothesis)

        gmax = np.random.uniform(0.1, 0.5)
        dmax = np.random.uniform(0.01, 0.1) # probability of dispersing in a given time step
        maintmax = np.random.uniform(0.0005, 0.005) # maximum metabolic maintanence cost


        # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS

    return [width, height, alpha, motion, reproduction, speciation, \
            seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, amp, freq, \
            flux, pulse, phase, disturb, envgrads, barriers, rates]
