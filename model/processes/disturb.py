# -*- coding: utf-8 -*-
from random import sample

def disturb(iD, ceiling):
    iDs = list(iD)
    iDs = sample(iDs, len(iDs) - ceiling/2)

    f_dict = {key: iD[key] for key in iD if key not in iDs}

    return f_dict
