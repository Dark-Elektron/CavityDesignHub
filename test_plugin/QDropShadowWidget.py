from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor, QDragEnterEvent, QDropEvent
from PyQt5.QtWidgets import QGraphicsDropShadowEffect, QGridLayout


class QDropShadowWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):

        if parent is None:
            super().__init__()
        else:
            super().__init__(parent=parent)

        self.setFixedSize(150, 150)
        # gl = QGridLayout()
        # self.setLayout(gl)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.setGraphicsEffect(shadow)
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setStyleSheet('background-color: red;')

        # self.setAttribute(Qt.Widget())
        self.setWindowFlags(Qt.Widget)

    def dragEnterEvent(self, event):
        self.setStyleSheet('background-color: blue;')
        self.update()
        super(QDropShadowWidget, self).dragEnterEvent(event)

    def dropEvent(self, event):
        super().dropEvent(event)

    def dragLeaveEvent(self, event):
        super().dragLeaveEvent(event)

    def dragMoveEvent(self, event):
        super().dragMoveEvent(event)

    def mouseMoveEvent(self, event):
        self.setStyleSheet('background-color: yellow;')
        self.update()
        super(QDropShadowWidget, self).mouseMoveEvent(event)

    def enterEvent(self, event):
        super(QDropShadowWidget, self).enterEvent(event)

    def mousePressEvent(self, event):
        self.setStyleSheet('background-color: white;')
        self.update()
        super(QDropShadowWidget, self).mousePressEvent(event)

