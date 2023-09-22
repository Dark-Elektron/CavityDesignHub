from ui_files.pp_plot import Ui_PPPlot
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
from analysis_modules.plot_module.plotter import Plot


class PPPlot:
    def __init__(self, parent):
        self.w_PPPlot = QWidget()

        self.ui = Ui_PPPlot()
        self.ui.setupUi(self.w_PPPlot)

        self.main_control = parent.main_control
        self.main_ui = parent.main_ui

        # Plot
        self.plt = Plot(self)
        self.ui.gl_Plot_Area.addWidget(self.plt)
        self.fig = self.plt.fig
        self.ax = self.plt.ax

        # # Create main window object
        # self.win = parent
        # self.main_control = parent
        # self.main_ui = parent.ppUI

    def signals(self):
        pass

    def initUI(self):
        pass

    def get_stylesheet(self):
        return self.w_PPPlot.styleSheet()
