from __future__ import division
from random import randint, choice
import numpy as np
import sys
import os
import math

mydir = os.path.expanduser("~/Desktop/Repos/HYDRO-BIDE/results/movies")

limit = 0.5

def NewTracers(TracerIDs, TracerXcoords, TracerYcoords, width, height, u0):

    x = np.random.binomial(1, u0/4)

    if x == 1:
        TracerYcoords.append(float(np.random.uniform(0.2*height, 0.8*height)))
        TracerXcoords.append(float(np.random.uniform(0.01, 0.02)))
        TracerIDs.append(0)

    return [TracerIDs, TracerXcoords, TracerYcoords]


def GetRAD(COM):

    RAD = []
    tList = list(set(COM))

    for i, sp in enumerate(tList):
        RAD.append(COM.count(sp))

    return RAD


def e_var(RAD): # Smith and Wilson's evenness metric (Evar)
    P = np.log(RAD)
    S = len(RAD)
    X = 0
    for x in P:
        X += (x - np.mean(P))**2/S
    evar = 1 - 2/np.pi*np.arctan(X)
    return(evar)


def get_color(ID, color_dict): # FUNCTION TO ASSIGN COLORS TO MICROBE SPECIES

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    color_dict[ID] = color

    return color_dict


def immigration(args1, args2):

    COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height, MaintDict = args1
    u0, GrowthDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, LogSeriesAlpha = args2

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity

    # immigration
    x = np.random.binomial(1, u0)
    if x == 1:

        props = np.random.logseries(LogSeriesAlpha, 2)

        for i, prop in enumerate(props):

            COM.append(prop)
            MicYcoords.append(float(np.random.uniform(0.1*height, 0.9*height)))
            MicXcoords.append(float(np.random.uniform(0.1, 2))) # width - limit
            MicIDs.append(MicID)
            MicTimeIn.append(0)
            MicID += 1
            Q = float(np.random.uniform(100, 400))
            MicQs.append(Q)

            if prop not in microbe_color_dict:
                microbe_color_dict = get_color(prop, microbe_color_dict) # species color
            if prop not in GrowthDict:
                GrowthDict[prop] = float(np.random.uniform(0.1, 0.5)) # species growth rate
            if prop not in MaintDict:
                MaintDict[prop] = float(np.random.uniform(0.5, 1)) # species maintenance


    return [[COM, MicXcoords, MicYcoords, MicExitAge, MaintDict],
        [GrowthDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs]]



def dispersal(args1, args2):

    COM, ux, uy, MicXcoords, MicYcoords, MicExitAge, width, height = args1
    u0, MicIDs, MicID, MicTimeIn, MicQs = args2

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity


    # dispersal inside the system
    for i, val in enumerate(MicXcoords):

        X = int(round(MicXcoords[i]))
        Y = int(round(MicYcoords[i]))

        index =  int(round(X + Y * width))

        if index > len(ux) - 1: index = len(ux) - 1
        if index > len(uy) - 1: index = len(uy) - 1

        k = np.random.binomial(1, 1)
        if k == 1:
            MicXcoords[i] += ux[index]
            MicYcoords[i] += uy[index]

        y = MicYcoords[i]
        if 0 > y: MicYcoords[i] = 0
        elif y > height: MicYcoords[i] = height

        MicTimeIn[i] += 1
        if MicXcoords[i] >= width - limit:

            MicExitAge.append(MicTimeIn[i])
            MicXcoords.pop(i)
            MicYcoords.pop(i)
            MicIDs.pop(i)
            MicQs.pop(i)
            MicTimeIn.pop(i)
            COM.pop(i)

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    return [[COM, MicXcoords, MicYcoords, MicExitAge],
        [MicIDs, MicID, MicTimeIn, MicQs]]





def maintenance(args):

    COM, MicXcoords, MicYcoords, MicExitAge, microbe_color_dict, MaintDict, MicIDs, MicID, MicTimeIn, MicQs, height, width = args

    for i, val in enumerate(MicQs):

        val -= MaintDict[COM[i]]  # maint influenced by species

        if val < 5:   # starved

            MicExitAge.append(MicTimeIn[i])

            COM.pop(i)
            MicXcoords.pop(i)
            MicYcoords.pop(i)
            MicIDs.pop(i)
            MicQs.pop(i)
            MicTimeIn.pop(i)

        else: MicQs[i] = val

    return [COM, MicXcoords, MicYcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs]



def resFlow(args):

    RES, ux, uy, u0, ResXcoords, ResYcoords, width, height, ResID, ResIDs = args

    ux = np.reshape(ux, (width*height))       # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height))       # uy is the macroscopic y velocity

    for i, val in enumerate(ResXcoords):

        if len(ResXcoords) != len(ResYcoords):
            print 'different lengths:', len(ResXcoords), len(ResYcoords)
            sys.exit()

        X = int(round(ResXcoords[i]))
        Y = int(round(ResYcoords[i]))

        index =  int(round(X + Y*width))

        if index > len(ux) - 1: index = len(ux) - 1
        if index > len(uy) - 1: index = len(uy) - 1

        ResXcoords[i] += ux[index]
        ResYcoords[i] += uy[index]
        y = ResYcoords[i]

        if 0 > y: ResYcoords[i] = 0
        elif y > height: ResYcoords[i] = height

        if ResXcoords[i] >= width - limit:

            ResXcoords.pop(i)
            ResYcoords.pop(i)
            ResIDs.pop(i)
            RES.pop(i)


    # resources flow in
    x = np.random.binomial(1, 1)
    if x == 1:
        res_in = np.random.random_integers(100, 200, u0*10)

        for i, val in enumerate(res_in):

            RES.append(val)
            ResIDs.append(ResID)
            ResID += 1
            ResYcoords.append(float(np.random.uniform(0.1*height, 0.9*height)))
            ResXcoords.append(float(np.random.uniform(0.1, 20))) # width-limit

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    return [RES, ResXcoords, ResYcoords, ResID, ResIDs]



def ConsumeAndReproduce(args):

    RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict = args

    BoxesOfMicrobes = [list([]) for _ in xrange(width*height)]
    BoxesOfResources = [list([]) for _ in xrange(width*height)]

    for i, val in enumerate(MicIDs):

        roundedX = int(round(MicXcoords[i]))
        roundedY = int(round(MicYcoords[i]))
        index = int(round(roundedX + roundedY * width))

        if index > len(BoxesOfMicrobes) - 1:
            index = len(BoxesOfMicrobes) - 1

        BoxesOfMicrobes[index].append(val)

    for i, val in enumerate(ResIDs):

        roundedX = int(round(ResXcoords[i]))
        roundedY = int(round(ResYcoords[i]))
        index = int(round(roundedX + roundedY * width))

        if index > len(BoxesOfResources) - 1:
            index = len(BoxesOfResources) - 1

        BoxesOfResources[index].append(val)


    for i, MicrobeBox in enumerate(BoxesOfMicrobes):
        ResourceBox = BoxesOfResources[i]

        while MicrobeBox:
            # The resource
            if len(ResourceBox): pass
            else: break

            resID = choice(ResourceBox)

            if resID in ResIDs:  # a check
                pass
            else:
                print 'something wrong: line 565'
                print ResourceBox, resID, ResIDs
                sys.exit()

            # The food
            j = ResIDs.index(resID)
            food = RES[j]

            # The microbe
            micID = choice(MicrobeBox)
            MicrobeBox.remove(micID)
            index = MicIDs.index(micID)

            Q = MicQs[index]
            mu = GrowthDict[COM[index]]

            if food > mu * Q: # Increase microbe cell quota
                Q += mu * Q
                food -= mu * Q
            else: # Increase microbe cell quota
                Q += food
                food = 0

            MicQs[index] = Q
            BoxesOfResources[i].remove(resID)

            if food <= 0:
                ResXcoords.pop(j)
                ResYcoords.pop(j)
                RES.pop(j)
                ResIDs.pop(j)

            else: RES[j] = food

            if Q > 100:
                                                                     # reproduce
                spID = COM[index]
                X = MicXcoords[index]
                Y = MicYcoords[index]

                MicQs[index] = Q/2.0

                newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                if newX > width - limit:
                    newX = width - limit

                newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                if 0 > newY: newY = 0
                elif newY > height: newY = height

                MicQs.append(Q/2.0)
                COM.append(spID)
                MicXcoords.append(newX)
                MicYcoords.append(newY)
                MicIDs.append(MicID)
                MicTimeIn.append(0)
                MicID += 1

    return [RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords]



def MoveTracers(args):
    TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy, COM, Ns, Ss, EVs, NSs, GrowthDict, MUs, VarMUs, Maints, VarMaints, MaintDict = args

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity

    # move the inert tracer particles
    for i, val in enumerate(TracerIDs):

        # move all tracers
        X = int(round(TracerXcoords[i]))
        Y = int(round(TracerYcoords[i]))
        index =  int(round(X + Y * width))

        if index > len(ux) - 1: index = len(ux) - 1
        if index > len(uy) - 1: index = len(uy) - 1

        TracerXcoords[i] += ux[index]
        TracerYcoords[i] += uy[index]
        TracerIDs[i] += 1

        if TracerXcoords[i] >= width - limit:

            TracerExitAge.append(TracerIDs[i])
            TracerXcoords.pop(i)
            TracerYcoords.pop(i)
            TracerIDs.pop(i)

            S = len(list(set(COM)))
            N = len(COM)

            if S >= 1:
               Ss.append(S)
               Ns.append(N)

               NSs.append(float(N)/float(S))
               RAD = GetRAD(COM)
               EVs.append(e_var(RAD))

               SpList = list(set(COM))
               MuList = []
               MaintList = []

               for index, Sp in enumerate(SpList):
                   mu = GrowthDict[Sp]
                   list1 = [mu] * RAD[index]
                   MuList.extend(list1)

                   maint = MaintDict[Sp]
                   list2 = [maint] * RAD[index]
                   MaintList.extend(list2)

               AvgMu = float(np.mean(MuList))
               MUs.append(AvgMu)

               VarMu = float(np.var(MuList, ddof=1))
               VarMUs.append(VarMu)

               AvgMaint = float(np.mean(MaintList))
               Maints.append(AvgMaint)

               VarMaint = float(np.var(MaintList, ddof=1))
               VarMaints.append(VarMaint)

            elif S == 0:

               Ss.append(0)
               Ns.append(0)

            continue

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    return [TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, width, height, ux, uy, COM, Ns, Ss, EVs, NSs, GrowthDict, MUs, VarMUs, Maints, VarMaints, MaintDict]
    #return [TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, ux, uy, COM, Ns, Ss, EVs, NSs]
