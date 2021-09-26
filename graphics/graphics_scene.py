from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import math
from PyQt5.QtGui import QBrush, QPen


class GraphicsScene(QGraphicsScene):
    def __init__(self, scene, parent=None):
        super().__init__(parent)

        self.scene = scene

        # Settings
        self.gridSize = 20
        self.gridSquares = 5

        self._color_background = QColor("#ffffff")
        self._color_light = QColor("#2f2f2f")
        self._color_dark = QColor("#292929")
        self._color_center = QColor("#ffffff")

        self._pen_light = QPen(self._color_light)
        self._pen_light.setWidth(1)
        self._pen_dark = QPen(self._color_dark)
        self._pen_dark.setWidth(2)
        self._pen_center = QPen(self._color_center)
        self._pen_center.setWidth(2)

        self.setBackgroundBrush(self._color_background)

    def setGrScene(self, width, height):
        self.setSceneRect(-width//2, -height//2, width, height)

    def drawBackground(self, painter, rect):
        super().drawBackground(painter, rect)
        grid = False
        if grid:
            # create grid
            left = int(math.floor(rect.left()))
            right = int(math.ceil(rect.right()))
            top = int(math.floor(rect.top()))
            bottom = int(math.ceil(rect.bottom()))

            first_left = left - (left % self.gridSize)
            first_top = top - (top % self.gridSize)

            # Compute all lines to be drawn
            lines_light, lines_dark = [], []
            for x in range(first_left, right, self.gridSize):
                if (x % (self.gridSize*self.gridSquares)) != 0:
                    lines_light.append(QLine(x, top, x, bottom))
                else:
                    lines_dark.append(QLine(x, top, x, bottom))

            for y in range(first_top, bottom, self.gridSize):
                if (y % (self.gridSize*self.gridSquares)) != 0:
                    lines_light.append(QLine(left, y, right, y))
                else:
                    lines_dark.append(QLine(left, y, right, y))

            # Draw the lines
            painter.setPen(self._pen_light)
            painter.drawLines(*lines_light)

            painter.setPen(self._pen_dark)
            painter.drawLines(*lines_dark)

            #### DRAW CENTER AXIS
            lines_center = []
            # print("hor", first_left, 0, right, 0)
            # print("ver", 0, first_top, 0, bottom)
            center_line_H = QLine(left, 0, right, 0)
            center_line_V = QLine(0, top, 0, bottom)
            lines_center.append(center_line_H)
            lines_center.append(center_line_V)

            painter.setPen(self._pen_center)
            painter.drawLines(*lines_center)

