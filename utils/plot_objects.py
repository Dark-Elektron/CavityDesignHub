import sys

from PyQt5 import *
from PyQt5 import Qt
from PyQt5.QtCore import QPropertyAnimation
from PyQt5.QtWidgets import *

from ui_files.plot_widget_header import Ui_Header
from ui_files.abci_plot import Ui_ABCI_Plot
# from node_editor.node_editor_window import NodeEditorWindow
# nw = NodeEditorWindow()

class PlotObjects:
    def __init__(self):
        self.w_PlotObj = QWidget()

        self.ui = Ui_Header()
        self.ui.setupUi(self.w_PlotObj)

        self.w_ABCI = QWidget()
        self.ui_abci = Ui_ABCI_Plot()
        self.ui_abci.setupUi(self.w_ABCI)

        spacerItem = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.ui.gridLayout.addItem(spacerItem, 2, 0, 1, 1)

        self.ui.gridLayout.addWidget(self.w_ABCI)

        self.ui.pb_Drop.clicked.connect(lambda: self.animate_height(self.ui.pb_Drop, 0, 350, True))

    def animate_height(self, cb, min_height, standard, enable):
        if enable:
            # select widget

            if self.ui.cb_Type.currentText() == "abci":
                widget = self.w_ABCI
            else:
                return
            #### GET WIDTH
            height = widget.width()

            #### SET MAX WIDTH
            if cb.isChecked():
                heightCollapsed = min_height
                widget.setMinimumHeight(0)
            else:
                heightCollapsed = standard
                # widget.setMinimumWidth(standard)

            #### ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumHeight")
            self.animation.setDuration(200)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.start()

    def show(self):
        self.w_PlotObj.show()


if __name__ == '__main__':

    app = QApplication(sys.argv)
    w = PlotObjects()
    w.show()
    sys.exit(app.exec_())
