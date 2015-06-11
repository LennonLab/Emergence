from __future__ import division
from random import randint, choice
import numpy as np
import sys
from scipy import stats
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics


limit = 0.5

def xcoord(low=0.2, high=8):
    return float(np.random.uniform(low, high))

def ycoord(low=1.5, high=8.5):
    return float(np.random.uniform(low, high))



def NewTracers(TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, width, height, u0, D):

    x = np.random.binomial(1, u0)

    if x == 1:
        y = coord()
        TracerYcoords.append(y)
        x = coord()
        TracerXcoords.append(x)
        TracerIDs.append(0) # used to track the age of the tracer

    if D == 3:
        z = coord()
        TracerZcoords.append(z)

    return [TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords]



def ResIn(RES, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResType, r, rmax, nr, width, height, u0, D):

    R = len(RES)
    x = np.random.binomial(1, u0)

    if x == 1:
        res_in = np.random.random_integers(1, rmax, r)
        res_types = np.random.random_integers(0, nr-1, r)

        for i, val in enumerate(res_in):

            RES.append(val)
            ResIDs.append(ResID)
            ResType.append(res_types[i])
            ResID += 1

            y = ycoord()
            ResYcoords.append(y)
            if R == 0:
                x = xcoord()
            else:
                x = xcoord()
            ResXcoords.append(x)

    return [RES, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs, ResType]



def immigration(m, COM, MicXcoords, MicYcoords, MicZcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict, nr, u0, LogSeriesAlpha, D):

    N = len(COM)
    x = np.random.binomial(1, u0)

    if x == 1:
        props = np.random.logseries(LogSeriesAlpha, m)
        #props = list(np.random.randint(1, 1000, m))

        for i, prop in enumerate(props):

            if prop <= 1000:
                COM.append(prop)

                y = ycoord()
                MicYcoords.append(y)

                if N == 0:
                    x = xcoord()
                else: x = xcoord()

                MicXcoords.append(x)

                MicIDs.append(MicID)
                MicTimeIn.append(0)
                MicID += 1
                Q = float(np.random.uniform(50, 100))
                MicQs.append(Q)

                if prop not in microbe_color_dict:
                    # species color
                    microbe_color_dict = get_color(prop, microbe_color_dict)

                    # species growth rate
                    GrowthDict[prop] = np.random.uniform(0.5, 1)

                    # species maintenance
                    MaintDict[prop] = np.random.uniform(5, 15)

                    # species active dispersal rate
                    DispParamDict[prop] = 0.0 # np.random.uniform(0.09, 0.9)

                    # species resource use efficiency
                    ResUseDict[prop] = np.random.uniform(0.1, 0.99, nr)

    return [COM, MicXcoords, MicYcoords, MicZcoords, width, height, MaintDict, GrowthDict, DispParamDict, microbe_color_dict, MicIDs, MicID, MicTimeIn, MicQs, ResUseDict]



def GetRAD(COM):
    RAD = []
    tList = list(set(COM)) # species IDs
    for i, sp in enumerate(tList):
        RAD.append(COM.count(sp)) # the abundance of each species
    return RAD, tList # the rad and the species list



def get_color(ID, color_dict): # FUNCTION TO ASSIGN COLORS TO MICROBE SPECIES

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    color_dict[ID] = color

    return color_dict


def dispersal(COM, ux, uy, MicXcoords, MicYcoords, MicZcoords, MicExitAge, width, height, u0, MicIDs, MicID, MicTimeIn, MicQs, D):

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

    return [COM, MicXcoords, MicYcoords, MicZcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs]



def maintenance(COM, MicXcoords, MicYcoords, MicZcoords, MicExitAge, microbe_color_dict, MaintDict, MicIDs, MicID, MicTimeIn, MicQs, height, width, D):

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

    return [COM, MicXcoords, MicYcoords, MicZcoords, MicExitAge, MicIDs, MicID, MicTimeIn, MicQs]




def resFlow(RES, ux, uy, u0, ResXcoords, ResYcoords, ResZcoords, width, height, ResID, ResIDs, D):

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

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    return [RES, ResXcoords, ResYcoords, ResZcoords, ResID, ResIDs]



def ConsumeAndReproduce(RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, width, height, GrowthDict, ResType, ResUseDict):

    BoxesOfMicrobes = [list([]) for _ in xrange(width*height)]
    BoxesOfResources = [list([]) for _ in xrange(width*height)]

    for i, val in enumerate(MicIDs):

        roundedX = int(round(MicXcoords[i]))
        roundedY = int(round(MicYcoords[i]))
        index = int(round(roundedX + (roundedY * width)))

        if index > len(BoxesOfMicrobes) - 1:
            index = len(BoxesOfMicrobes) - 1
        elif index < 0:
            index = 0

        BoxesOfMicrobes[index].append(val)

    for i, val in enumerate(ResIDs):

        roundedX = int(round(ResXcoords[i]))
        roundedY = int(round(ResYcoords[i]))
        index = int(round(roundedX + roundedY * width))

        if index > len(BoxesOfMicrobes) - 1:
            index = len(BoxesOfMicrobes) - 1
        elif index < 0:
            index = 0

        BoxesOfResources[index].append(val)


    for i, MicrobeBox in enumerate(BoxesOfMicrobes):
        ResourceBox = BoxesOfResources[i]

        while MicrobeBox: # The resource
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
            restype = ResType[j]
            food = RES[j]

            # The microbe
            micID = choice(MicrobeBox)
            MicrobeBox.remove(micID)
            index = MicIDs.index(micID)

            Q = MicQs[index]
            sp = COM[index]
            mu = GrowthDict[COM[index]] * ResUseDict[sp][restype]

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
                ResType.pop(j)
                ResIDs.pop(j)

            else: RES[j] = food

            if Q > 500: # reproduction

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

    BoxesOfMicrobes = [list([]) for _ in xrange(width*height)]
    BoxesOfResources = [list([]) for _ in xrange(width*height)]

    return [RES, ResIDs, ResXcoords, ResYcoords, COM, MicIDs, MicID, MicTimeIn, MicQs, MicXcoords, MicYcoords, ResType]



def MoveTracers(TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, width, height, ux, uy, T, R, RES, ResType, ResDens, ResDiv, ResRich, TracerTau, MicrobeTau, MicExitAge, COM, N, S, Mu, Maint, GrowthDict, MaintDict, Ev, ES, Nm , BP , SD , sk, D):

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity

    # move the inert tracer particles
    for i, val in enumerate(TracerIDs):

        # move all tracers
        X = int(round(TracerXcoords[i]))
        Y = int(round(TracerYcoords[i]))
        index =  int(round(X + Y * width))

        if index > len(ux) - 1: index = len(ux) - 1
        elif index < 0: index = 0

        TracerXcoords[i] += ux[index]
        TracerYcoords[i] += uy[index]
        TracerIDs[i] += 1

        if TracerXcoords[i] >= width - limit:

            TracerExitAge.append(TracerIDs[i])
            TracerXcoords.pop(i)
            TracerYcoords.pop(i)
            TracerIDs.pop(i)

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    return [TracerExitAge, TracerIDs, TracerXcoords, TracerYcoords, TracerZcoords, T, R, ResDens, ResDiv, ResRich, TracerTau, MicrobeTau, N, S, Mu, Maint, Ev, ES, Nm , BP , SD , sk]
