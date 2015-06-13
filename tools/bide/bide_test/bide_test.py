from __future__ import division
from random import randint, choice
import numpy as np
import sys

limit = 0.5


def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.pyx """

    motion = 'fluid'
    width = 10 #choice([5, 5])
    height = 10 #choice([5, 5])

    alpha = 0.99 #np.random.uniform(0.95, 0.99)

    # size of starting community
    seedcom = 0 #choice([10, 10, 100, 1000])

    # individuals immigrating per time step
    m = 8 #choice([0, 2, 4, 8])

    # resource particles flowing in per time step
    r = 100# choice([20, 40, 80, 160])

    # maximum number of resources types
    nr = 10# choice([1, 2, 4, 8, 16, 32])

    # maximum resource particle size
    rmax = 4000 #choice([500, 1000, 2000, 4000, 8000])


    return [width, height, alpha, motion, seedcom, m, r, nr, rmax]


def coord(d):
    return float(np.random.uniform(0.2*d, 0.8*d))


def GetRAD(vector):
    RAD = []
    unique = list(set(vector))
    for val in unique:
        RAD.append(vector.count(val)) # the abundance of each species

    return RAD, unique # the rad and the species list


def get_color(ID, color_dict): # FUNCTION TO ASSIGN COLORS TO SPECIES

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    color_dict[ID] = color

    return color_dict



def NewTracers(IDs, Xcoords, Ycoords, TimeIn, width, height, u0):

    x = np.random.binomial(1, u0)

    if x == 1:
        TimeIn.append(0)

        y = coord(height)
        Ycoords.append(y)

        x = coord(width)
        Xcoords.append(x)

        IDs.append(0) # used to track the age of the tracer

    return [IDs, TimeIn, Xcoords, Ycoords,]



def ResIn(Type, Vals, Xcoords, Ycoords, ID, IDs, TimeIn, numr, rmax, nr, width, height, u0):

    for r in range(numr):
        x = np.random.binomial(1, u0)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            rtype = int(np.random.random_integers(0, nr-1, 1))

            Vals.append(rval)
            IDs.append(ID)
            Type.append(rtype)
            TimeIn.append(0)
            ID += 1

            y = coord(height)
            Ycoords.append(y)

            x = coord(width)
            Xcoords.append(x)

    return [Type, Vals, Xcoords, Ycoords, IDs, ID, TimeIn]



def immigration(numin, Species, Xcoords, Ycoords, width, height, MaintDict, GrowthDict, DispParamDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict, nr, u0, alpha):

    for m in range(numin):
        x = np.random.binomial(1, u0)

        if x == 1:
            prop = int(np.random.logseries(alpha, 1))

            if prop > 1000: continue

            Species.append(prop)

            y = coord(height)
            Ycoords.append(y)

            x = coord(width)
            Xcoords.append(x)

            IDs.append(ID)
            TimeIn.append(0)
            ID += 1
            Q = float(np.random.uniform(50, 100))
            Qs.append(Q)

            if prop not in color_dict:
                # species color
                color_dict = get_color(prop, color_dict)

                # species growth rate
                GrowthDict[prop] = np.random.uniform(.5, 1)

                # species maintenance
                MaintDict[prop] = np.random.uniform(5, 15)

                # species active dispersal rate
                DispParamDict[prop] = np.random.uniform(.2, 2)

                # species resource use efficiency
                ResUseDict[prop] = np.random.uniform(0, 1, nr)

    return [Species,  Xcoords, Ycoords, MaintDict, GrowthDict, DispParamDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict]



def fluid_movement(TypeOf, List, TimeIn, ExitAge, Xcoords, Ycoords, ux, uy, width, height, u0):

    if TypeOf == 'resource' or TypeOf == 'individual':
        Type, IDs, ID, Vals = List

    else: IDs = List

    if Xcoords == []:
        if TypeOf == 'resource' or TypeOf == 'individual':
            return [Type, Xcoords, Ycoords, ExitAge, IDs, ID, TimeIn, Vals]

        elif TypeOf == 'tracer':
            return [IDs, Xcoords, Ycoords, ExitAge, TimeIn]

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity

    # dispersal inside the system
    for i, val in enumerate(Xcoords):

        X = int(round(Xcoords[i]))
        Y = int(round(Ycoords[i]))

        index =  int(round(X + Y * width))

        if index > len(ux) - 1: index = len(ux) - 1
        if index > len(uy) - 1: index = len(uy) - 1

        k = np.random.binomial(1, 1)
        if k == 1:
            Xcoords[i] += ux[index]
            Ycoords[i] += uy[index]

        y = Ycoords[i]
        if 0 > y: Ycoords[i] = 0
        elif y > height: Ycoords[i] = height

        TimeIn[i] += 1
        if Xcoords[i] >= width - limit:

            ExitAge.append(TimeIn[i])
            Xcoords.pop(i)
            Ycoords.pop(i)
            TimeIn.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    if TypeOf == 'tracer':
        return [IDs, Xcoords, Ycoords, ExitAge, TimeIn]

    elif TypeOf == 'resource' or TypeOf == 'individual':
        return [Type, Xcoords, Ycoords, ExitAge, IDs, ID, TimeIn, Vals]



def maintenance(Species, Xcoords, Ycoords, ExitAge, color_dict, MaintDict, IDs, TimeIn, Qs):

    if Species == []:
        return [Species, Xcoords, Ycoords, ExitAge, IDs, TimeIn, Qs]

    for i, val in enumerate(Qs):

        val -= MaintDict[Species[i]]  # maintanence influenced by species

        if val < 1:   # starved

            Qs.pop(i)
            ExitAge.append(TimeIn[i])
            TimeIn.pop(i)
            Species.pop(i)
            IDs.pop(i)

            Xcoords.pop(i)
            Ycoords.pop(i)

        else: Qs[i] = val

    return [Species, Xcoords, Ycoords, ExitAge, IDs, TimeIn, Qs]



def ConsumeAndReproduce(ResTypes, ResIDs, ResXcoords, ResYcoords, SpeciesIDs, IndIDs, IndID, IndTimeIn, Qs, IndXcoords, IndYcoords, ResVals, width, height, GrowthDict, ResUseDict):

    BoxesOfMicrobes = [list([]) for _ in xrange(width*height)]
    BoxesOfResources = [list([]) for _ in xrange(width*height)]

    for i, val in enumerate(IndIDs):

        roundedX = int(round(IndXcoords[i]))
        roundedY = int(round(IndYcoords[i]))
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
            restype = ResTypes[j]
            food = ResVals[j]

            # The microbe
            micID = choice(MicrobeBox)
            MicrobeBox.remove(micID)
            index = IndIDs.index(micID)

            Q = Qs[index]
            sp = SpeciesIDs[index]
            mu = GrowthDict[SpeciesIDs[index]] * ResUseDict[sp][restype]

            if food > mu * Q: # Increase microbe cell quota
                Q += mu * Q
                food -= mu * Q
            else: # Increase microbe cell quota
                Q += food
                food = 0

            Qs[index] = Q
            BoxesOfResources[i].remove(resID)

            if food <= 0:
                ResXcoords.pop(j)
                ResYcoords.pop(j)
                ResTypes.pop(j)
                ResVals.pop(j)
                ResIDs.pop(j)

            else: ResVals[j] = food

            if Q > 500: # reproduction

                spID = SpeciesIDs[index]
                X = IndXcoords[index]
                Y = IndYcoords[index]

                Qs[index] = Q/2.0

                newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                if newX > width - limit:
                    newX = width - limit

                newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                if 0 > newY: newY = 0
                elif newY > height: newY = height

                Qs.append(Q/2.0)
                SpeciesIDs.append(spID)
                IndXcoords.append(newX)
                IndYcoords.append(newY)
                IndIDs.append(IndID)
                IndTimeIn.append(0)
                IndID += 1

    BoxesOfMicrobes = [list([]) for _ in xrange(width*height)]
    BoxesOfResources = [list([]) for _ in xrange(width*height)]

    return [ResTypes, ResIDs, ResXcoords, ResYcoords, SpeciesIDs, IndIDs, IndID, IndTimeIn, Qs, IndXcoords, IndYcoords, ResVals]
