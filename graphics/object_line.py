from graphics.graphics_line_2d import *
# from content_widgets.functional_streamline_content_widget import QFunctionalStreamlinesContentWidget
# from content_widgets.mode_content_widget import QModeContentWidget
# from content_widgets.node_content_widget import QNodeContentWidget
# from content_widgets.streamline_content_widget import QStreamlinesContentWidget
# from content_widgets.substructure_content_widget import QSubstructureContentWidget
# from content_widgets.utility_content_widget import QUtilityContentWidget
# from objects.node_socket import *
# from config import *
# import json
# from collections import OrderedDict
# from serialize import Serializable

DEBUG = False

# THIS IS THE OBJECT THAT HOLDS THE PAINTED LINE AND CAN BE ADDED TO THE SCENE

class Line:
    def __init__(self, main_win=None, line_type='line', pos=(0, 0), section='top', IC=None, OC=None, BP=None, color="#ffffff"):
        super().__init__()
        self.main_win = main_win
        # self.vtk = self.main_win.grVTK.vtk_widget
        # self._title = title

        self.scene = main_win.scene

        self.line_type = line_type
        self.section = section
        self.instantiateNodeType(IC, OC, BP, color)

        # THIS IS WHERE THE LINE IS ADDED, FIRST TO THE SCENE IN LINES_LIST AND THEN TO THE GRAPHICS SCENE
        self.scene.addLine(self)
        self.setPos(*pos)
        self.scene.grScene.addItem(self.grLine)

    def draw_cavity(self):
        print()

    @property
    def pos(self):
        return self.grLine.pos()  # returns QPointF

    def setPos(self, x, y):
        self.grLine.setPos(x, y)

    def remove(self):
        self.scene.grScene.removeItem(self.grLine)
        self.scene.removeLine(self)

    def instantiateNodeType(self, IC, OC, BP, color):
        self.grLine = QGraphicsLineDirect(self, color, section=self.section, IC=IC, OC=OC, BP=BP)
