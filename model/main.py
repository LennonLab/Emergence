from random import shuffle
from os.path import expanduser
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/model")

from processes import *
from diversity_metrics import *
from randomized_parameters import *
from spatial_functions import * 
from input_output import *


labels.clear()


def iter_processes(processes, IndDict, SpDict, ResDict, params):

    for process in shuffle(processes):
    
        if process is 'resource_inflow': # Inflow of resources
            ResDict, RID = bide.ResIn(ResDict, params)
    
        elif process is 'resource_flow': # Resource flow
            ResDict, RID = bide.res_flow(ResDict, params)
    
        elif process is 'immigration': # Inflow of individuals (immigration)
            SpDict, IndDict, ID = bide.immigration(SpDict, IndDict, params)
    
        elif process is 'passive dispersal': # flowthrough of individuals
            IndDict = bide.ind_flow(SpDict, IndDict, params)
    
        elif process is 'active disperal': # Active dispersal
            IndDict = bide.ind_disp(SpDict, IndDict, params)
    
        elif process is 'searching': # Search
            IndDict = bide.search(SpDict, IndDict, ResDict, params)
    
        elif process is 'consume': # Consume
            numc, SpDict, IndDict, ResDict = bide.consume(SpDict, IndDict, ResDict, params)
    
        elif process is 'growth': # Grow
            IndDict = bide.grow(SpDict, IndDict)
    
        elif process is 'transition': # Transition
            IndDict = bide.transition(SpDict, IndDict)
    
        elif process is 'maintenance': # Maintenance
            IndDict = bide.maintenance(SpDict, IndDict)
    
        elif process is 'reproduce': # Reproduction
            SpDict, IndDict, ID, prodI, prodN = bide.reproduce(SpDict, IndDict, params)
            
    return [IndDict, SpDict, ResDict]
    
    
    
def run_model(sim, ct2):
    params = rp.get_rand_params()
    SpDict, IndDict, ResDict, Ns = {}, {}, {}, []
    
    while len(Ns) < 1000 or max(Ns[-10]) == 0:
        IndDict, SpDict, ResDict = iter_processes(processes, IndDict, SpDict, ResDict, params)
        output.output(IndDict, SpDict, ResDict, sim, ct2)

    print 'ct:', '%4s' % ct2, 'sim:', '%3s' % sim,'  N:', '%4s' %  Ns[-1], '\n'
    
    
processes = labels.processes()
for sim in range(10**5): run_model(processes, sim)