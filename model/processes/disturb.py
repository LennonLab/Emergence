# -*- coding: utf-8 -*-
from random import sample

def disturb(iD, ceiling):

    ''' A function to bring resource particles into a system

        iD  : A python dictionary to hold properties of individual
                   organisms

                   'size' : The size of the individual is an abstract
                   quantity, but can vary over two orders of magnitude

                   'q' : level of endogenous resources
                   'spID' : species ID
                   'state' : whether active or dormant

                  'x' : x-coordinate
                  'y' : y-coordinate
                  'z' : z-coordinate
    '''

    iDs = list(iD)
    iDs = sample(iDs, len(iDs) - ceiling)

    f_dict = {key: iD[key] for key in iD if key not in iDs}

    return f_dict
