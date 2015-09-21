#!/bin/sh

""":"
exec python $0 ${1+"$@"}
"""

# ---------------------- HelloBlt.py ----------------------
# 
# This program demonstrates more basic functionality like 
# making buttons to control showing of symbols, grids etc.
# The example also shows how to make animated graphs and 
# how to plot the graphs more or less smooth.
# 

from Tkinter import *        # The Tk package
import Pmw                   # The Python MegaWidget package
import math                  # import the sin-function


def animate():
    # This function is completely pointless, but demonstrates
    # that it's easy to update a graph "runtime".

    for t in range(31):                   # In 31 steps...
        for c in range(ncurves):          # ...on each curve
            for x in range(npoints):      # on each point...
                vector_y[c][x] = math.sin(c*x*0.5 +math.pi*t/15)

        master.after(20)            # wait 0.02 second
        master.update_idletasks()   # update screen

def symbolsOnOff():
    global symbols
    symbols = not symbols

    for curvename in g.element_show():
        if symbols:
            g.element_configure(curvename, symbol='diamond')
        else:
            g.element_configure(curvename, symbol='')
            

def smooth():
    global smoothing
    
    if smoothing == 'linear': smoothing='quadratic'
    elif smoothing == 'quadratic': smoothing='natural'
    elif smoothing == 'natural': smoothing='step'
    else: smoothing = 'linear'

    for curvename in g.element_show():
        g.element_configure(curvename, smooth=smoothing)


master = Tk()                # build Tk-environment
ncurves = 4                  # draw 4 curves
npoints = 7                  # use  8 points on each curve

smoothing='linear'
symbols  = 0

# In this example we use Pmw.Blt.Vectors. These can mostly be used like 
# a normal list, but changes will be updated in the graph automatically.
# Using Pmw.Blt.Vectors is often slower, but in this case very convenient.
vector_x = Pmw.Blt.Vector()   
vector_y = []

for y in range(ncurves):
   vector_y.append(Pmw.Blt.Vector())          # make vector for y-axis

for x in range(npoints):                      # for each point...
   vector_x.append(x*0.1)                     # make an x-value

   # fill vectors with cool graphs
   for c in range(ncurves):                   # for each curve...
      vector_y[c].append(math.sin(c*x*0.5))   # make an y-value

g = Pmw.Blt.Graph(master)                     # make a new graph area
g.pack(expand=1, fill='both')

color = ['red', '#ff9900', 'blue', '#00cc00', 'black', 'grey']

for c in range(ncurves):                      # for each curve...
   curvename = 'sin(' + str(c) +'x)'          # make a curvename
   g.line_create(curvename,                   # and create the graph
                 xdata=vector_x,              # with x data,
                 ydata=vector_y[c],           # and  y data
                 color=color[c],              # and a color
                 dashes=0,                    # and no dashed line
                 linewidth=2,                 # and 2 pixels wide
                 symbol='')                   # ...and no disks
   
g.configure(title='Hello BLT World')          # enter a title
   
# make s row of buttons
buttons = Pmw.ButtonBox(master, labelpos='n', label_text='Options')
buttons.pack(fill='both', expand=1, padx=10, pady=10)

buttons.add('Grid',       command=g.grid_toggle)
buttons.add('Symbols',    command=symbolsOnOff)
buttons.add('Smooth',     command=smooth)
buttons.add('Animate',    command=animate)
buttons.add('Quit',       command=master.quit)

master.mainloop()                              # ...and wait for input


