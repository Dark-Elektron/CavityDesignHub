from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
from ui_files.edit_annotated_text import Ui_EAT


class EditATextDialog:
    def __init__(self, parent, object=None):
        self.w_eat = QWidget()

        # disable main window when popup is active
        self.w_eat.setWindowModality(Qt.ApplicationModal)

        self.eatUI = Ui_EAT()
        self.eatUI.setupUi(self.w_eat)

        # Create main window object
        self.parent = parent
        self.parent_ui = self.parent.widgetUI

        self.text = ""
        self.annotate_object = object

        self.initUI()
        self.signals()

    def initUI(self):
        if self.annotate_object:
            self.eatUI.le_Text.setText(self.annotate_object.get_text())
            # self.eatUI.cb_Box.setCurrentText(self.annotate_object.get_bbox_patch())
            self.eatUI.sb_Font_Size.setValue(self.annotate_object.get_fontsize())

    def signals(self):
        self.eatUI.pb_Apply.clicked.connect(lambda: self.apply())
        self.eatUI.pb_Cancel.clicked.connect(lambda: self.close())

    def get_text(self):
        return self.text

    def apply(self):
        txt = self.eatUI.le_Text.text()
        box = self.eatUI.cb_Box.currentText()
        font_size = self.eatUI.sb_Font_Size.value()

        if self.annotate_object is None:
            # create new text
            if txt == "":
                pass
            else:
                self.parent.add_text(txt, box, font_size)
                self.parent.draw()
                self.close()
        else:
            if txt == "":
                pass
            else:
                self.annotate_object.set_text(txt)
                self.annotate_object.set_size(font_size)

                if box != "None":
                    bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
                    self.annotate_object.set_bbox(bbox_props)
                else:
                    self.annotate_object.set_bbox({'facecolor': 'white', 'edgecolor': 'none', 'pad': 10})

                self.parent.draw()
                self.close()

    def close(self):
        self.w_eat.close()
        self.w_eat.setParent(None)
        del self.w_eat

    def show(self):
        self.w_eat.show()
