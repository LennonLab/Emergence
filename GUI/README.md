#HydroBIDE GUI

Code for a graphical user interface based (for the moment) on a template provided by Brendan Griffen at MIT.

http://brendangriffen.com/blog/creating-a-GUI-in-Python/

**Everything below is also from Brendan**

## Description
A basic template for making a GUI in Python using Traits. It gives you access to a matplotlib canvas object and Mayavi 3D rendering object which can be dynamically updated by adjusting the objects on the right panel. Tabbed viewing is also enabled by default.

## Requirements
You will need Enthought's EPD distribution (32-bit, **not** 64-bit as it does not include Mayavi) which can be found [here](https://www.enthought.com/repo/epd/installers/). You will need an academic license to access it for free.

## Installing EPD:
`bash epd-7.3-2-rh5-x86.sh` ('rh3' for older operating systems) or just use the `.dmg` image for OSX.

## Running GUI
`ipython` then `run main.py` or `python main.py` outright.

## Adding your own features.
There is very comprehensive documentation on using Traits [here](http://code.enthought.com/projects/traits/documentation.php) and a lot of [tutorials](http://docs.enthought.com/traitsui/tutorials/index.html). I will be giving a few introductory sessions on my [blog](http://bgriffen.scripts.mit.edu/www/) as well.
