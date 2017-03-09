from __future__ import division
import numpy as np


l = 5 # x
h = 5 # y
w = 5 # z

# should be l*w*h boxes of microbes

xcoords = [0, 1, 1, 0, 4, 4] # l = 5
ycoords = [0, 1, 0, 1, 1, 4] # h = 2
zcoords = [0, 0, 1, 1, 2, 4] # w = 3

com = [6,10,20,30,7, 4]

BoM = [list([]) for _ in xrange(l*h*w)]
print len(BoM)

for i, val in enumerate(com):

    roundedX = int(round(xcoords[i]))
    roundedY = int(round(ycoords[i]))
    roundedZ = int(round(zcoords[i]))

    print roundedX, roundedY, roundedZ,' : ',

    index1 = int(round((roundedY * l * w) + (roundedX * w) + roundedZ))
    index2 = int(round((roundedZ * l * h) + (roundedY * l) + roundedX))
    index3 = int(round((roundedX * h * w) + (roundedY * w) + roundedZ))

    print index1, index2, index3
