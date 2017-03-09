from __future__ import division
import numpy as np

""" This code was initially inspired and is almost entirely based on Dan
Schroeder's http://physics.weber.edu/schroeder/fluids/, obtained in May 2014 """


def SetLattice(u0, viscosity, width, height, lefts, bottoms, barriers): # Lattice Boltzmann PARAMETERS

    four9ths = 4.0/9.0    # abbreviations for lattice-Boltzmann weight factors
    one9th   = 1.0/9.0
    one36th  = 1.0/36.0

    ############## particle densities along 9 directions
    n0 = four9ths * (np.ones((height,width)) - 1.5*u0**2)
    nN = one9th * (np.ones((height,width)) - 1.5*u0**2)
    nS = one9th * (np.ones((height,width)) - 1.5*u0**2)
    nE = one9th * (np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nW = one9th * (np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nNE = one36th * (np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nSE = one36th * (np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nNW = one36th * (np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nSW = one36th * (np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

    rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW   # macroscopic density
    ux = (nE + nNE + nSE - nW - nNW - nSW) / rho	# macroscopic x velocity
    uy = (nN + nNE + nNW - nS - nSE - nSW) / rho	# macroscopic y velocity

    barrier = np.zeros((height, width), bool)  # Initialize barriers
                                               # True wherever there's a barrier

    BarrierHeight = 0.1
    BarrierWidth = 0.1

    for i in range(barriers):

        left = lefts[i]
        bottom = bottoms[i]

        barrier[bottom*height: (bottom+BarrierHeight)*height, left*width: (left+BarrierWidth)*width] = True

    barrierN = np.roll(barrier,  1, axis=0)       # sites just north of barriers
    barrierS = np.roll(barrier, -1, axis=0)       # sites just south of barriers
    barrierE = np.roll(barrier,  1, axis=1)       # etc.
    barrierW = np.roll(barrier, -1, axis=1)
    barrierNE = np.roll(barrierN,  1, axis=1)
    barrierNW = np.roll(barrierN, -1, axis=1)
    barrierSE = np.roll(barrierS,  1, axis=1)
    barrierSW = np.roll(barrierS, -1, axis=1)

    return [n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier, rho, ux, uy, barrierN, barrierS, barrierE, barrierW, barrierNE, barrierNW, barrierSE, barrierSW]



############################################### LBM FUNCTIONS ##################
def curl(ux, uy): # function to compute curl of the macroscopic velocity field
    curlXY = np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)
    curlXY = curlXY - np.roll(ux, -1, axis=0) + np.roll(ux, 1, axis=0)

    return curlXY



def stream(args):
    nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier = args

    # particle barriers

    barrierN = np.roll(barrier,  1, axis=0)	  # sites just north of barriers
    barrierS = np.roll(barrier, -1, axis=0)	  # sites just south of barriers
    barrierE = np.roll(barrier,  1, axis=1)				  # etc.
    barrierW = np.roll(barrier, -1, axis=1)
    barrierNE = np.roll(barrierN,  1, axis=1)
    barrierNW = np.roll(barrierN, -1, axis=1)
    barrierSE = np.roll(barrierS,  1, axis=1)
    barrierSW = np.roll(barrierS, -1, axis=1)

    # function to move particles one step along their directions of motion

    nN  = np.roll(nN,   1, axis=0) # axis 0 is north-south; + direction is north
    nNE = np.roll(nNE,  1, axis=0)
    nNW = np.roll(nNW,  1, axis=0)
    nS  = np.roll(nS,  -1, axis=0)
    nSE = np.roll(nSE, -1, axis=0)
    nSW = np.roll(nSW, -1, axis=0)
    nE  = np.roll(nE,   1, axis=1)    # axis 1 is east-west; + direction is east
    nNE = np.roll(nNE,  1, axis=1)
    nSE = np.roll(nSE,  1, axis=1)
    nW  = np.roll(nW,  -1, axis=1)
    nNW = np.roll(nNW, -1, axis=1)
    nSW = np.roll(nSW, -1, axis=1)

    # Use tricky boolean arrays to handle barrier collisions (bounce-back)

    nN[barrierN] = nS[barrier]
    nS[barrierS] = nN[barrier]
    nE[barrierE] = nW[barrier]
    nW[barrierW] = nE[barrier]
    nNE[barrierNE] = nSW[barrier]
    nNW[barrierNW] = nSE[barrier]
    nSE[barrierSE] = nNW[barrier]
    nSW[barrierSW] = nNE[barrier]


    return [nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier]



def collide(viscosity, rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW, u0):
    # Collide particles within each cell to redistribute velocities

    four9ths = 4.0/9.0    # abbreviations for lattice-Boltzmann weight factors
    one9th   = 1.0/9.0
    one36th  = 1.0/36.0
    omega = 1 / (3 * viscosity + 0.5)  # relaxation parameter

    rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW
    ux = (nE + nNE + nSE - nW - nNW - nSW) / rho
    uy = (nN + nNE + nNW - nS - nSE - nSW) / rho

    ux2 = ux * ux			  # pre-compute terms used repeatedly...
    uy2 = uy * uy
    u2 = ux2 + uy2

    omu215 = 1 - 1.5*u2			                # one minus u2 times 1.5
    uxuy = ux * uy
    OxR = omega * rho

    n0 = (1-omega)*n0 + OxR * four9ths * omu215
    nN = (1-omega)*nN + OxR * one9th * (omu215 + 3*uy + 4.5*uy2)
    nS = (1-omega)*nS + OxR * one9th * (omu215 - 3*uy + 4.5*uy2)
    nE = (1-omega)*nE + OxR * one9th * (omu215 + 3*ux + 4.5*ux2)
    nW = (1-omega)*nW + OxR * one9th * (omu215 - 3*ux + 4.5*ux2)
    nNE = (1-omega)*nNE + OxR * one36th * (omu215 + 3*(ux+uy) + 4.5*(u2+2*uxuy))
    nNW = (1-omega)*nNW + OxR * one36th * (omu215 + 3*(-ux+uy) +4.5*(u2-2*uxuy))
    nSE = (1-omega)*nSE + OxR * one36th * (omu215 + 3*(ux-uy) + 4.5*(u2-2*uxuy))
    nSW = (1-omega)*nSW + OxR * one36th * (omu215 + 3*(-ux-uy) +4.5*(u2+2*uxuy))

     # Force steady rightward flow at ends (no need to set 0, N, & S components)
    nE[:,0] = one9th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nW[:,0] = one9th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nNE[:,0] = one36th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nSE[:,0] = one36th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nNW[:,0] = one36th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nSW[:,0] = one36th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

    return [rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW]
