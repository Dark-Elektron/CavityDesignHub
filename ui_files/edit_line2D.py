# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_files/edit_line2D.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_EAL(object):
    def setupUi(self, EAL):
        EAL.setObjectName("EAL")
        EAL.resize(258, 99)
        self.gridLayout = QtWidgets.QGridLayout(EAL)
        self.gridLayout.setObjectName("gridLayout")
        self.pb_Cancel = QtWidgets.QPushButton(EAL)
        self.pb_Cancel.setObjectName("pb_Cancel")
        self.gridLayout.addWidget(self.pb_Cancel, 1, 1, 1, 1)
        self.pb_Apply = QtWidgets.QPushButton(EAL)
        self.pb_Apply.setObjectName("pb_Apply")
        self.gridLayout.addWidget(self.pb_Apply, 1, 0, 1, 1)
        self.widget = QtWidgets.QWidget(EAL)
        self.widget.setObjectName("widget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.widget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label = QtWidgets.QLabel(self.widget)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.dsb_X = QtWidgets.QDoubleSpinBox(self.widget)
        self.dsb_X.setSingleStep(0.01)
        self.dsb_X.setProperty("value", 0.0)
        self.dsb_X.setObjectName("dsb_X")
        self.gridLayout_2.addWidget(self.dsb_X, 0, 1, 1, 1)
        self.gridLayout.addWidget(self.widget, 0, 0, 1, 2)

        self.retranslateUi(EAL)
        QtCore.QMetaObject.connectSlotsByName(EAL)

    def retranslateUi(self, EAL):
        _translate = QtCore.QCoreApplication.translate
        EAL.setWindowTitle(_translate("EAL", "Form"))
        self.pb_Cancel.setText(_translate("EAL", "Cancel"))
        self.pb_Apply.setText(_translate("EAL", "Apply"))
        self.label.setText(_translate("EAL", "Move to"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    EAL = QtWidgets.QWidget()
    ui = Ui_EAL()
    ui.setupUi(EAL)
    EAL.show()
    sys.exit(app.exec_())