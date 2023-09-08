from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
from ui_files.edit_annotated_text import Ui_EAT


class EditATextDialog:
    def __init__(self, parent, obj=None):
        self.w_eat = QWidget()

        # disable main window when popup is active
        self.w_eat.setWindowModality(Qt.ApplicationModal)

        self.eatUI = Ui_EAT()
        self.eatUI.setupUi(self.w_eat)

        # Create main window obj
        self.parent = parent
        self.parent_ui = self.parent.widgetUI

        self.text = ""
        self.annotate_obj = obj

        self.initUI()
        self.signals()

    def initUI(self):
        if self.annotate_obj:
            self.eatUI.le_Text.setText(self.annotate_obj.text.get_text())
            # self.eatUI.cb_Box.setCurrentText(self.annotate_obj.get_bbox_patch())
            self.eatUI.sb_Font_Size.setValue(int(self.annotate_obj.text.get_fontsize()))

    def signals(self):
        self.eatUI.pb_Apply.clicked.connect(lambda: self.apply())
        self.eatUI.pb_Cancel.clicked.connect(lambda: self.close())

    def get_text(self):
        return self.text

    def apply(self):
        txt = self.eatUI.le_Text.text()
        box = self.eatUI.cb_Box.currentText()
        font_size = self.eatUI.sb_Font_Size.value()
        rotation = self.eatUI.sb_Rotation.value()

        if self.annotate_obj is None:
            # create new text
            if txt == "":
                pass
            else:
                self.parent.add_text(text=txt, box=box, size=font_size, rotation=rotation)
                self.parent.draw()
                self.close()
        else:
            if txt == "":
                pass
            else:
                self.annotate_obj.text.set_text(txt)
                self.annotate_obj.text.set_size(font_size)
                self.annotate_obj.text.set_rotation(rotation)

                if box != "None":
                    bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
                    self.annotate_obj.text.set_bbox(bbox_props)
                else:
                    self.annotate_obj.text.set_bbox({'facecolor': 'none', 'edgecolor': 'none', 'pad': 10})

                self.parent.draw_idle()
                self.parent.flush_events()
                self.close()

    def close(self):
        self.w_eat.close()
        self.w_eat.setParent(None)
        del self.w_eat

    def show(self):
        self.w_eat.show()

#
# class EditPatchDialog:
#     def __init__(self, parent, obj=None):
#         self.w_epd = QWidget()
#
#         # disable main window when popup is active
#         self.w_epd.setWindowModality(Qt.ApplicationModal)
#
#         self.epdUI = Ui_EPD()
#         self.epdUI.setupUi(self.w_epd)
#
#         # Create main window obj
#         self.parent = parent
#         self.parent_ui = self.parent.widgetUI
#
#         self.text = ""
#         self.annotate_obj = obj
#
#         self.initUI()
#         self.signals()
#
#     def initUI(self):
#         if self.annotate_obj:
#             self.epdUI.le_Text.setText(self.annotate_obj.get_text())
#             # self.epdUI.cb_Box.setCurrentText(self.annotate_obj.get_bbox_patch())
#             self.epdUI.sb_Font_Size.setValue(int(self.annotate_obj.get_fontsize()))
#
#     def signals(self):
#         self.epdUI.pb_Apply.clicked.connect(lambda: self.apply())
#         self.epdUI.pb_Cancel.clicked.connect(lambda: self.close())
#
#     def get_text(self):
#         return self.text
#
#     def apply(self):
#         txt = self.epdUI.le_Text.text()
#         box = self.epdUI.cb_Box.currentText()
#         font_size = self.epdUI.sb_Font_Size.value()
#         rotation = self.epdUI.sb_Rotation.value()
#
#         if self.annotate_obj is None:
#             # crepde new text
#             if txt == "":
#                 pass
#             else:
#                 self.parent.add_text(text=txt, box=box, size=font_size, rotation=rotation)
#                 self.parent.draw()
#                 self.close()
#         else:
#             if txt == "":
#                 pass
#             else:
#                 self.annotate_obj.set_text(txt)
#                 self.annotate_obj.set_size(font_size)
#                 self.annotate_obj.set_rotation(rotation)
#
#                 if box != "None":
#                     bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
#                     self.annotate_obj.set_bbox(bbox_props)
#                 else:
#                     self.annotate_obj.set_bbox({'facecolor': 'white', 'edgecolor': 'none', 'pad': 10})
#
#                 self.parent.draw()
#                 self.close()
#
#     def close(self):
#         self.w_epd.close()
#         self.w_epd.setParent(None)
#         del self.w_epd
#
#     def show(self):
#         self.w_epd.show()
