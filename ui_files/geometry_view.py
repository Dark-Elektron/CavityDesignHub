# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_files/geometry_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_GeometryView(object):
    def setupUi(self, GeometryView):
        GeometryView.setObjectName("GeometryView")
        GeometryView.resize(1182, 995)
        self.gridLayout = QtWidgets.QGridLayout(GeometryView)
        self.gridLayout.setObjectName("gridLayout")
        self.widget = QtWidgets.QWidget(GeometryView)
        self.widget.setObjectName("widget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.widget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.sp_Left_Right_Container = QtWidgets.QSplitter(self.widget)
        self.sp_Left_Right_Container.setOrientation(QtCore.Qt.Horizontal)
        self.sp_Left_Right_Container.setObjectName("sp_Left_Right_Container")
        self.widget_2 = QtWidgets.QWidget(self.sp_Left_Right_Container)
        self.widget_2.setBaseSize(QtCore.QSize(200, 0))
        self.widget_2.setObjectName("widget_2")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.widget_2)
        self.gridLayout_4.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.gl_InputWidgets = QtWidgets.QGridLayout()
        self.gl_InputWidgets.setObjectName("gl_InputWidgets")
        self.gridLayout_4.addLayout(self.gl_InputWidgets, 0, 0, 1, 1)
        self.w_Show_Cavity = QtWidgets.QWidget(self.sp_Left_Right_Container)
        self.w_Show_Cavity.setStyleSheet("border-radius: 15px;")
        self.w_Show_Cavity.setObjectName("w_Show_Cavity")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.w_Show_Cavity)
        self.gridLayout_6.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.gl_Plot_Area = QtWidgets.QGridLayout()
        self.gl_Plot_Area.setObjectName("gl_Plot_Area")
        self.gridLayout_6.addLayout(self.gl_Plot_Area, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.sp_Left_Right_Container, 0, 0, 1, 1)
        self.gridLayout.addWidget(self.widget, 0, 0, 1, 1)

        self.retranslateUi(GeometryView)
        QtCore.QMetaObject.connectSlotsByName(GeometryView)

    def retranslateUi(self, GeometryView):
        _translate = QtCore.QCoreApplication.translate
        GeometryView.setWindowTitle(_translate("GeometryView", "Form"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    GeometryView = QtWidgets.QWidget()
    ui = Ui_GeometryView()
    ui.setupUi(GeometryView)
    GeometryView.show()
    sys.exit(app.exec_())
