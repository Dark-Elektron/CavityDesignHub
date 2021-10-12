import ast
import json
import os

import matplotlib.pyplot as plt
import matplotlib
import scipy
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtWidgets import *
from utils.file_reader import FileReader
from modules.data_module.slans_data import SLANSData
from modules.data_module.abci_data import ABCIData

# from Plotter.zoom_pan import ZoomPan

matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg#, NavigationToolbar2QT as NavigationToolbar
from utils.CustomNavBar import MyCustomToolbar as NavigationToolbar
from matplotlib.figure import Figure
import utils.misc_functions as f
# from matplotlib_annotation_objects import DraggableText
import matplotlib as mpl
plt.rcParams['toolbar'] = 'toolmanager'
from matplotlib.backend_tools import ToolBase, ToolToggleBase


fr = FileReader()


class Plot(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=12, height=9, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.set_tight_layout(True)
        self.ax = self.fig.add_subplot(111)

        mpl.rcParams['savefig.facecolor'] = 'ffffff'

        super(Plot, self).__init__(self.fig)

        self.widgetUI = parent
        # self.plotUI = self.plot_widget.plotUI
        # self.MROT = self.ui.cb_Polarization_ABCI.currentIndex()

        self.GLOBAL_STATE = 0
        # plot appearance
        # self.aesthetics()
        self.plot_settings('presentation')

        # handle events
        # self.event_handler()

        # zoom and pan
        # scale = 1.1
        # zp = ZoomPan()
        # figZoom = zp.zoom_factory(self.ax, base_scale=scale)
        # # figPan = zp.pan_factory(self.axes)

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        self.toolbar = NavigationToolbar(self, None)

        # create new axis
        self.create_secondary_axes()

        # # create new widgets
        # self.create_new_widgets()

        # redefine save_image function
        self.toolbar.sizeHint()
        self.widgetUI.gl_Plot_Area.addWidget(self.toolbar)

        # handle events
        self.event_handler()

        self.signals()

    def signals(self):
        pass

    def change_scale(self, arg):
        self.ax.set_yscale(arg)
        self.ax_right.set_yscale(arg)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def plot_settings(self, mode):
        plt.style.use(['science', 'no-latex'])
        if mode == 'presentation':
            # plt.rcParams["figure.figsize"] = (11, 6)
            #
            # mpl.rcParams['lines.markersize'] = 10
            #
            # # axes
            # mpl.rcParams['axes.labelsize'] = 20
            # mpl.rcParams['axes.titlesize'] = 25
            #
            # mpl.rcParams['xtick.labelsize'] = 15
            # mpl.rcParams['ytick.labelsize'] = 15
            #

            # plot settings
            mpl.rcParams['xtick.labelsize'] = 15
            mpl.rcParams['ytick.labelsize'] = 15

            mpl.rcParams['axes.labelsize'] = 15
            mpl.rcParams['axes.titlesize'] = 15
            mpl.rcParams['legend.fontsize'] = 'large'

            mpl.rcParams['figure.figsize'] = [9.8, 6]
            mpl.rcParams['figure.dpi'] = 100

    def create_secondary_axes(self):
        # Create multiple axes
        # twin object for two different y-axis on the sample plot
        self.ax_right = self.ax.twinx()
        # turn off axis ticks and labels
        self.ax_right.axes.xaxis.set_visible(False)
        self.ax_right.axes.yaxis.set_visible(False)

        # # top axis # code incomplete
        # self.ax_top = self.ax.twinx()
        # self.ax_top.spines['top'].set_color(colors[0])
        #
        # # bottom axis # code incomplete
        # self.ax_bottom = self.ax.twinx()
        # self.ax_bottom.spines['bottom'].set_color(colors[0])

    def event_handler(self):
        self.object = None
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        cid2 = self.fig.canvas.mpl_connect('pick_event', self.onpick)

    def onclick(self, event):
        # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
        #       ('double' if event.dblclick else 'single', event.button,
        #        event.x, event.y, event.xdata, event.ydata))

        # handle identifying objects
        if event.button == 1:
            self.bt = 1
            # print("Left click")

        if event.button == 2:
            self.bt = 2
            # print("MMB click")
        if event.button == 3:
            self.bt = 3
            print("Right click")
            if not isinstance(self.object, matplotlib.text.Annotation):
                self.contextMenuEvent(event)
            # selecting objects right click
            if isinstance(self.object, matplotlib.text.Annotation):
                # print("We got ourselves a right click")
                self.draggableTextContextMenuEvent(event)

    def onpick(self, event):
        print("Pick")
        # legend toggle
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        print("legline: ", legline)
        origline = self.lined[legline]
        print("origline: ", origline)
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)

        self.fig.canvas.draw()

    def contextMenuEvent(self, event):
        contextMenu = QMenu(self)
        textAct = contextMenu.addAction("Add Text")
        arrowAct = contextMenu.addAction("Add Arrow")
        arrowAct = contextMenu.addAction("Move axis marker to minimum")
        arrowAct = contextMenu.addAction("Move axis marker to maximum")
        print("\t", QtCore.QPoint(event.x, event.y))
        print(self.mapToGlobal(QtCore.QPoint(event.x, -event.y)))

        # get figure size and subtract from event.y because event.y counts bottom up.
        action = contextMenu.exec_(self.mapToGlobal(QtCore.QPoint(event.x, self.fig.get_size_inches()[1]*self.fig.dpi-event.y)))

        # selected = None
        # item = self.scene.getItemAt(event.pos())
        #
        # if hasattr(item, 'edge'):
        #     selected = item.edge
        #
        # if selected and action == bezierAct: selected.edge_type = EDGE_TYPE_BEZIER
        # if selected and action == directAct: selected.edge_type = EDGE_TYPE_DIRECT
        # if action == deleteAct:
        #     self.close()


class ZoomPan:
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None


    def zoom_factory(self, ax, base_scale = 2.):
        def zoom(event):
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()

            xdata = event.xdata # get event x location
            ydata = event.ydata # get event y location

            if event.button == 'down':
                # deal with zoom in
                scale_factor = 1 / base_scale
            elif event.button == 'up':
                # deal with zoom out
                scale_factor = base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1
                print(event.button)

            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])

            ax.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            # ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
            ax.figure.canvas.draw()

        fig = ax.get_figure() # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom

    def pan_factory(self, ax):
        def onPress(event):
            if event.inaxes != ax: return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press

        def onRelease(event):
            self.press = None
            ax.figure.canvas.draw()

        def onMotion(event):
            if self.press is None: return
            if event.inaxes != ax: return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax.set_xlim(self.cur_xlim)
            ax.set_ylim(self.cur_ylim)

            ax.figure.canvas.draw()

        fig = ax.get_figure() # get the figure of interest

        # attach the call back
        fig.canvas.mpl_connect('button_press_event',onPress)
        fig.canvas.mpl_connect('button_release_event',onRelease)
        fig.canvas.mpl_connect('motion_notify_event',onMotion)

        #return the function
        return onMotion


