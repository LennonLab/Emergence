from matplotlib.pyplot import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

import traitsui.api
from traitsui.api import View,UItem, Item, Group, Heading, Label, \
        HSplit, Handler, CheckListEditor, EnumEditor, TableEditor, \
        ListEditor, Tabbed, VGroup, HGroup, RangeEditor, Spring, spring

from traitsui.api import NoButtons
from traitsui.wx.editor import Editor
from traitsui.wx.basic_editor_factory import BasicEditorFactory

from traits.api import HasTraits

import numpy as np
from enable.api import ColorTrait

from matplotlib import *
import os
import sys
import wx
import platform
import socket
import random

class _MPLFigureEditor(Editor):

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # The panel lets us add additional controls.

        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)

        # matplotlib commands to create a canvas
        mpl_control = FigureCanvas(panel, -1, self.value)
        sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        toolbar = NavigationToolbar2Wx(mpl_control)
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.value.canvas.SetMinSize((10,10))

        return panel

class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor
