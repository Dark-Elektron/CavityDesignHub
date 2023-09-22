from analysis_modules.plot_module.plotter import Plot
from ui_files.geometry_view import Ui_GeometryView
from utils.shared_functions import *


class GeometryViewControl:
    def __init__(self, parent):
        self.w_GeometryView = QWidget()

        self.ui = Ui_GeometryView()
        self.ui.setupUi(self.w_GeometryView)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        self.plot = Plot(self)
        # fix axis aspect ratio
        self.plot.ax.set_aspect('equal', adjustable='datalim')
        self.ui.gl_Plot_Area.addWidget(self.plot)
        self.ui.sp_Left_Right_Container.setStretchFactor(1, 1)

    def get_stylesheet(self):
        return self.w_GeometryView.styleSheet()
