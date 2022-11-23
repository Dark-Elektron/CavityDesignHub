from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
from ui_files.edit_line2D import Ui_EAL


class EditALine2DDialog:
    def __init__(self, parent, object=None):
        self.w_eat = QWidget()

        # disable main window when popup is active
        self.w_eat.setWindowModality(Qt.ApplicationModal)

        self.ealUI = Ui_EAL()
        self.ealUI.setupUi(self.w_eat)

        # Create main window object
        self.parent = parent
        self.parent_ui = self.parent.widgetUI

        self.text = ""
        self.object = object

        self.initUI()
        self.signals()

    def initUI(self):
        if self.object:
            self.ealUI.dsb_X.setValue(self.object.get_xdata()[0])

    def signals(self):
        self.ealUI.pb_Apply.clicked.connect(lambda: self.apply())
        self.ealUI.pb_Cancel.clicked.connect(lambda: self.close())

    def apply(self):
        x0 = self.ealUI.dsb_X.value()
        if self.object is None:
            pass
        else:
            self.object.set_xdata([x0, x0])

        self.close()
        self.parent.draw()

    def close(self):
        self.w_eat.close()
        self.w_eat.setParent(None)
        del self.w_eat

    def show(self):
        self.w_eat.show()
