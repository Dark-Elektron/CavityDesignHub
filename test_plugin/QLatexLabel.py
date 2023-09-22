from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtCore import Qt, pyqtSlot, pyqtProperty
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pyparsing import ParseFatalException


class QLatexLabel(QtWidgets.QWidget):

    def __init__(self, text, parent=None):

        self.text = text

        if parent is None:
            super().__init__()
        else:
            super().__init__(parent=parent)

        self.ll = QVBoxLayout(self)
        self.ll.setContentsMargins(0, 0, 0, 0)

        r, g, b, a = self.palette().base().color().getRgbF()
        self._figure = Figure(edgecolor=(r, g, b), facecolor='none')
        self._canvas = FigureCanvas(self._figure)
        # self._canvas.setStyleSheet("background-color:transparent;")
        self.ll.addWidget(self._canvas)
        self.draw_text()

    def get_text(self):
        return self.text

    @pyqtSlot(str)
    def set_text(self, value):
        self.text = value
        self.update()
        try:
            self.draw_text()
        except:
            pass

    Text = pyqtProperty(str, get_text, set_text)

    def draw_text(self):
        self._figure.clear()
        try:
            text = self._figure.suptitle(
                self.text,
                x=0.0,
                y=1.0,
                horizontalalignment='left',
                verticalalignment='top',
                size=11, c='dimgray'  # QtGui.QFont().pointSize() * 1.5
            )
        except:
            text = self._figure.suptitle(
                self.text.strip(r"$"),
                x=0.0,
                y=1.0,
                horizontalalignment='left',
                verticalalignment='top',
                size=11, c='dimgray'  # QtGui.QFont().pointSize() * 1.5
            )
        self._canvas.draw()

        (x0, y0), (x1, y1) = text.get_window_extent().get_points()
        w = x1 - x0
        h = y1 - y0

        self._figure.set_size_inches(w / 100, h / 100)
        self.setFixedSize(int(w), int(h))

# if __name__ == '__main__':
#     from sys import argv, exit
#
#
#     class Widget(QtWidgets.QWidget):
#         def __init__(self, parent=None, **kwargs):
#             super(QtWidgets.QWidget, self).__init__(parent)
#
#             super().__init__(parent)
#             l = QVBoxLayout(self)
#             mathText = r'$X_k = \sum_{n=0}^{N-1} x_n . e^{\frac{-i2\pi kn}{N}}$'
#             l.addWidget(MathTextLabel(mathText, self), alignment=Qt.AlignHCenter)
#
#
#     a = QtWidgets.QApplication(argv)
#     w = Widget()
#     w.show()
#     w.raise_()
#     exit(a.exec_())
