# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_files/misc.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Misc(object):
    def setupUi(self, Misc):
        Misc.setObjectName("Misc")
        Misc.resize(1042, 625)
        Misc.setStyleSheet("\n"
"*{\n"
"    font: 16px \"Segoe UI\";\n"
"}")
        self.gridLayout = QtWidgets.QGridLayout(Misc)
        self.gridLayout.setObjectName("gridLayout")
        self.widget = QtWidgets.QWidget(Misc)
        self.widget.setObjectName("widget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.widget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.widget_4 = QtWidgets.QWidget(self.widget)
        self.widget_4.setObjectName("widget_4")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.widget_4)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.sa_Misc = QtWidgets.QScrollArea(self.widget_4)
        self.sa_Misc.setWidgetResizable(True)
        self.sa_Misc.setObjectName("sa_Misc")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 962, 545))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.pb_Waveguide_Calculations = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pb_Waveguide_Calculations.setMinimumSize(QtCore.QSize(200, 200))
        self.pb_Waveguide_Calculations.setMaximumSize(QtCore.QSize(200, 200))
        self.pb_Waveguide_Calculations.setStyleSheet("background-color: rgb(194, 108, 255);")
        self.pb_Waveguide_Calculations.setObjectName("pb_Waveguide_Calculations")
        self.gridLayout_6.addWidget(self.pb_Waveguide_Calculations, 0, 0, 1, 1)
        self.pb_Mode_Nomenclature = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pb_Mode_Nomenclature.setMinimumSize(QtCore.QSize(200, 200))
        self.pb_Mode_Nomenclature.setMaximumSize(QtCore.QSize(200, 200))
        self.pb_Mode_Nomenclature.setStyleSheet("background-color: rgb(255, 167, 66);")
        self.pb_Mode_Nomenclature.setObjectName("pb_Mode_Nomenclature")
        self.gridLayout_6.addWidget(self.pb_Mode_Nomenclature, 0, 1, 1, 1)
        self.pushButton_4 = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButton_4.setMinimumSize(QtCore.QSize(200, 200))
        self.pushButton_4.setMaximumSize(QtCore.QSize(200, 200))
        self.pushButton_4.setStyleSheet("background-color: rgb(185, 239, 255);")
        self.pushButton_4.setObjectName("pushButton_4")
        self.gridLayout_6.addWidget(self.pushButton_4, 0, 2, 1, 1)
        self.pushButton_5 = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButton_5.setMinimumSize(QtCore.QSize(200, 200))
        self.pushButton_5.setMaximumSize(QtCore.QSize(200, 200))
        self.pushButton_5.setStyleSheet("background-color: rgb(160, 255, 192);")
        self.pushButton_5.setObjectName("pushButton_5")
        self.gridLayout_6.addWidget(self.pushButton_5, 1, 0, 1, 1)
        self.sa_Misc.setWidget(self.scrollAreaWidgetContents)
        self.gridLayout_5.addWidget(self.sa_Misc, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.widget_4, 1, 0, 1, 1)
        self.gridLayout.addWidget(self.widget, 0, 0, 1, 1)

        self.retranslateUi(Misc)
        QtCore.QMetaObject.connectSlotsByName(Misc)

    def retranslateUi(self, Misc):
        _translate = QtCore.QCoreApplication.translate
        Misc.setWindowTitle(_translate("Misc", "Form"))
        self.pb_Waveguide_Calculations.setText(_translate("Misc", "WAVEGUIDE\n"
"CALCULATIONS"))
        self.pb_Mode_Nomenclature.setText(_translate("Misc", "MODE\n"
"NOMENCLATURE"))
        self.pushButton_4.setText(_translate("Misc", "Useful Formulas"))
        self.pushButton_5.setText(_translate("Misc", "PushButton"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Misc = QtWidgets.QWidget()
    ui = Ui_Misc()
    ui.setupUi(Misc)
    Misc.show()
    sys.exit(app.exec_())

