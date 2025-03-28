# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_files/plottypeselector.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_PlotTypeSelector(object):
    def setupUi(self, PlotTypeSelector):
        PlotTypeSelector.setObjectName("PlotTypeSelector")
        PlotTypeSelector.resize(350, 300)
        PlotTypeSelector.setMinimumSize(QtCore.QSize(350, 300))
        PlotTypeSelector.setMaximumSize(QtCore.QSize(350, 300))
        self.gridLayout = QtWidgets.QGridLayout(PlotTypeSelector)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.widget = QtWidgets.QWidget(PlotTypeSelector)
        self.widget.setMinimumSize(QtCore.QSize(345, 265))
        self.widget.setMaximumSize(QtCore.QSize(345, 265))
        self.widget.setObjectName("widget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.widget)
        self.gridLayout_3.setContentsMargins(0, 0, 0, -1)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.widget_2 = QtWidgets.QWidget(self.widget)
        self.widget_2.setStyleSheet("background-color: rgb(255, 170, 0);\n"
"border-radius: 10px;")
        self.widget_2.setObjectName("widget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.widget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.cb_code = QtWidgets.QComboBox(self.widget_2)
        self.cb_code.setMinimumSize(QtCore.QSize(75, 0))
        self.cb_code.setMaximumSize(QtCore.QSize(75, 16777215))
        self.cb_code.setObjectName("cb_code")
        self.cb_code.addItem("")
        self.cb_code.addItem("")
        self.cb_code.addItem("")
        self.gridLayout_2.addWidget(self.cb_code, 0, 0, 1, 1)
        self.pb_Drop = QtWidgets.QPushButton(self.widget_2)
        self.pb_Drop.setMinimumSize(QtCore.QSize(20, 25))
        self.pb_Drop.setMaximumSize(QtCore.QSize(20, 25))
        self.pb_Drop.setCheckable(True)
        self.pb_Drop.setObjectName("pb_Drop")
        self.gridLayout_2.addWidget(self.pb_Drop, 0, 3, 1, 1)
        self.switchControl = SwitchControl(self.widget_2)
        self.switchControl.setObjectName("switchControl")
        self.gridLayout_2.addWidget(self.switchControl, 0, 2, 1, 1)
        self.pb_Remove = QtWidgets.QPushButton(self.widget_2)
        self.pb_Remove.setMinimumSize(QtCore.QSize(30, 25))
        self.pb_Remove.setMaximumSize(QtCore.QSize(30, 16777215))
        self.pb_Remove.setObjectName("pb_Remove")
        self.gridLayout_2.addWidget(self.pb_Remove, 0, 4, 1, 1)
        self.gridLayout_3.addWidget(self.widget_2, 0, 0, 1, 1)
        self.sw_Pages = QtWidgets.QStackedWidget(self.widget)
        self.sw_Pages.setObjectName("sw_Pages")
        self.page = QtWidgets.QWidget()
        self.page.setObjectName("page")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.page)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_2 = QtWidgets.QLabel(self.page)
        self.label_2.setObjectName("label_2")
        self.gridLayout_4.addWidget(self.label_2, 3, 0, 1, 1)
        self.le_Folder_ABCI = QtWidgets.QLineEdit(self.page)
        self.le_Folder_ABCI.setObjectName("le_Folder_ABCI")
        self.gridLayout_4.addWidget(self.le_Folder_ABCI, 0, 1, 1, 5)
        self.label_9 = QtWidgets.QLabel(self.page)
        self.label_9.setObjectName("label_9")
        self.gridLayout_4.addWidget(self.label_9, 5, 2, 1, 1)
        self.cb_Style_ABCI = QtWidgets.QComboBox(self.page)
        self.cb_Style_ABCI.setObjectName("cb_Style_ABCI")
        self.gridLayout_4.addWidget(self.cb_Style_ABCI, 5, 5, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.page)
        self.label_5.setObjectName("label_5")
        self.gridLayout_4.addWidget(self.label_5, 1, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.page)
        self.label_4.setObjectName("label_4")
        self.gridLayout_4.addWidget(self.label_4, 2, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.page)
        self.label_7.setObjectName("label_7")
        self.gridLayout_4.addWidget(self.label_7, 4, 2, 1, 1)
        self.cb_Axis_ABCI = QtWidgets.QComboBox(self.page)
        self.cb_Axis_ABCI.setObjectName("cb_Axis_ABCI")
        self.cb_Axis_ABCI.addItem("")
        self.cb_Axis_ABCI.addItem("")
        self.gridLayout_4.addWidget(self.cb_Axis_ABCI, 5, 1, 1, 1)
        self.ccb_Id_ABCI = QCheckableComboBox(self.page)
        self.ccb_Id_ABCI.setObjectName("ccb_Id_ABCI")
        self.gridLayout_4.addWidget(self.ccb_Id_ABCI, 1, 1, 1, 5)
        self.label_8 = QtWidgets.QLabel(self.page)
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 5, 0, 1, 1)
        self.cb_Type_ABCI = QtWidgets.QComboBox(self.page)
        self.cb_Type_ABCI.setObjectName("cb_Type_ABCI")
        self.gridLayout_4.addWidget(self.cb_Type_ABCI, 5, 3, 1, 2)
        self.cb_Polarization_ABCI = QtWidgets.QComboBox(self.page)
        self.cb_Polarization_ABCI.setObjectName("cb_Polarization_ABCI")
        self.cb_Polarization_ABCI.addItem("")
        self.cb_Polarization_ABCI.addItem("")
        self.gridLayout_4.addWidget(self.cb_Polarization_ABCI, 2, 1, 1, 2)
        self.label_3 = QtWidgets.QLabel(self.page)
        self.label_3.setObjectName("label_3")
        self.gridLayout_4.addWidget(self.label_3, 0, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem, 6, 1, 1, 1)
        self.le_ScaleX_ABCI = QtWidgets.QLineEdit(self.page)
        self.le_ScaleX_ABCI.setObjectName("le_ScaleX_ABCI")
        self.gridLayout_4.addWidget(self.le_ScaleX_ABCI, 4, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.page)
        self.label_6.setObjectName("label_6")
        self.gridLayout_4.addWidget(self.label_6, 4, 0, 1, 1)
        self.le_ScaleY_ABCI = QtWidgets.QLineEdit(self.page)
        self.le_ScaleY_ABCI.setObjectName("le_ScaleY_ABCI")
        self.gridLayout_4.addWidget(self.le_ScaleY_ABCI, 4, 3, 1, 1)
        self.cb_Request_ABCI = QtWidgets.QComboBox(self.page)
        self.cb_Request_ABCI.setObjectName("cb_Request_ABCI")
        self.gridLayout_4.addWidget(self.cb_Request_ABCI, 3, 1, 1, 5)
        self.pb_Open_Folder_ABCI = QtWidgets.QPushButton(self.page)
        self.pb_Open_Folder_ABCI.setMinimumSize(QtCore.QSize(30, 25))
        self.pb_Open_Folder_ABCI.setMaximumSize(QtCore.QSize(30, 25))
        self.pb_Open_Folder_ABCI.setObjectName("pb_Open_Folder_ABCI")
        self.gridLayout_4.addWidget(self.pb_Open_Folder_ABCI, 0, 6, 1, 1)
        self.sw_Pages.addWidget(self.page)
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.page_2)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.le_ScaleY_SLANS = QtWidgets.QLineEdit(self.page_2)
        self.le_ScaleY_SLANS.setObjectName("le_ScaleY_SLANS")
        self.gridLayout_5.addWidget(self.le_ScaleY_SLANS, 4, 3, 1, 1)
        self.le_Filter_Value_SLANS = QtWidgets.QLineEdit(self.page_2)
        self.le_Filter_Value_SLANS.setObjectName("le_Filter_Value_SLANS")
        self.gridLayout_5.addWidget(self.le_Filter_Value_SLANS, 6, 4, 1, 2)
        self.ccb_Filter_Variable_SLANS = QCheckableComboBox(self.page_2)
        self.ccb_Filter_Variable_SLANS.setObjectName("ccb_Filter_Variable_SLANS")
        self.gridLayout_5.addWidget(self.ccb_Filter_Variable_SLANS, 6, 1, 1, 2)
        self.label_18 = QtWidgets.QLabel(self.page_2)
        self.label_18.setObjectName("label_18")
        self.gridLayout_5.addWidget(self.label_18, 4, 2, 1, 1)
        self.le_ScaleX_SLANS = QtWidgets.QLineEdit(self.page_2)
        self.le_ScaleX_SLANS.setObjectName("le_ScaleX_SLANS")
        self.gridLayout_5.addWidget(self.le_ScaleX_SLANS, 4, 1, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.page_2)
        self.label_21.setObjectName("label_21")
        self.gridLayout_5.addWidget(self.label_21, 6, 3, 1, 1)
        self.cb_Request_SLANS = QtWidgets.QComboBox(self.page_2)
        self.cb_Request_SLANS.setObjectName("cb_Request_SLANS")
        self.gridLayout_5.addWidget(self.cb_Request_SLANS, 3, 1, 1, 5)
        self.label_16 = QtWidgets.QLabel(self.page_2)
        self.label_16.setObjectName("label_16")
        self.gridLayout_5.addWidget(self.label_16, 5, 0, 1, 1)
        self.label_20 = QtWidgets.QLabel(self.page_2)
        self.label_20.setObjectName("label_20")
        self.gridLayout_5.addWidget(self.label_20, 5, 2, 1, 1)
        self.cb_Polarization_SLANS = QtWidgets.QComboBox(self.page_2)
        self.cb_Polarization_SLANS.setObjectName("cb_Polarization_SLANS")
        self.gridLayout_5.addWidget(self.cb_Polarization_SLANS, 2, 1, 1, 1)
        self.cb_Type_SLANS = QtWidgets.QComboBox(self.page_2)
        self.cb_Type_SLANS.setObjectName("cb_Type_SLANS")
        self.gridLayout_5.addWidget(self.cb_Type_SLANS, 5, 3, 1, 2)
        self.label_12 = QtWidgets.QLabel(self.page_2)
        self.label_12.setObjectName("label_12")
        self.gridLayout_5.addWidget(self.label_12, 3, 0, 1, 1)
        self.comboBox_7 = QtWidgets.QComboBox(self.page_2)
        self.comboBox_7.setObjectName("comboBox_7")
        self.gridLayout_5.addWidget(self.comboBox_7, 5, 5, 1, 1)
        self.cb_Axis_SLANS = QtWidgets.QComboBox(self.page_2)
        self.cb_Axis_SLANS.setObjectName("cb_Axis_SLANS")
        self.gridLayout_5.addWidget(self.cb_Axis_SLANS, 5, 1, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.page_2)
        self.label_14.setObjectName("label_14")
        self.gridLayout_5.addWidget(self.label_14, 4, 0, 1, 1)
        self.label_19 = QtWidgets.QLabel(self.page_2)
        self.label_19.setObjectName("label_19")
        self.gridLayout_5.addWidget(self.label_19, 2, 2, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.page_2)
        self.label_17.setObjectName("label_17")
        self.gridLayout_5.addWidget(self.label_17, 6, 0, 1, 1)
        self.ccb_Shape_ID_SLANS = QCheckableComboBox(self.page_2)
        self.ccb_Shape_ID_SLANS.setObjectName("ccb_Shape_ID_SLANS")
        self.gridLayout_5.addWidget(self.ccb_Shape_ID_SLANS, 2, 3, 1, 3)
        self.le_Folder_SLANS = QtWidgets.QLineEdit(self.page_2)
        self.le_Folder_SLANS.setObjectName("le_Folder_SLANS")
        self.gridLayout_5.addWidget(self.le_Folder_SLANS, 0, 1, 1, 5)
        self.label_13 = QtWidgets.QLabel(self.page_2)
        self.label_13.setObjectName("label_13")
        self.gridLayout_5.addWidget(self.label_13, 2, 0, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.page_2)
        self.label_15.setObjectName("label_15")
        self.gridLayout_5.addWidget(self.label_15, 0, 0, 1, 1)
        self.pb_Open_Folder_SLANS = QtWidgets.QPushButton(self.page_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(30)
        sizePolicy.setVerticalStretch(25)
        sizePolicy.setHeightForWidth(self.pb_Open_Folder_SLANS.sizePolicy().hasHeightForWidth())
        self.pb_Open_Folder_SLANS.setSizePolicy(sizePolicy)
        self.pb_Open_Folder_SLANS.setMinimumSize(QtCore.QSize(30, 25))
        self.pb_Open_Folder_SLANS.setObjectName("pb_Open_Folder_SLANS")
        self.gridLayout_5.addWidget(self.pb_Open_Folder_SLANS, 0, 6, 1, 1)
        self.sw_Pages.addWidget(self.page_2)
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.page_3)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.cb_Style_Other = QtWidgets.QComboBox(self.page_3)
        self.cb_Style_Other.setObjectName("cb_Style_Other")
        self.gridLayout_6.addWidget(self.cb_Style_Other, 4, 5, 1, 1)
        self.label_31 = QtWidgets.QLabel(self.page_3)
        self.label_31.setObjectName("label_31")
        self.gridLayout_6.addWidget(self.label_31, 5, 3, 1, 1)
        self.le_ScaleY_Other = QtWidgets.QLineEdit(self.page_3)
        self.le_ScaleY_Other.setObjectName("le_ScaleY_Other")
        self.gridLayout_6.addWidget(self.le_ScaleY_Other, 3, 5, 1, 1)
        self.label_26 = QtWidgets.QLabel(self.page_3)
        self.label_26.setObjectName("label_26")
        self.gridLayout_6.addWidget(self.label_26, 4, 0, 1, 1)
        self.ccb_Filter_Other = QCheckableComboBox(self.page_3)
        self.ccb_Filter_Other.setObjectName("ccb_Filter_Other")
        self.gridLayout_6.addWidget(self.ccb_Filter_Other, 5, 1, 1, 2)
        self.label_27 = QtWidgets.QLabel(self.page_3)
        self.label_27.setObjectName("label_27")
        self.gridLayout_6.addWidget(self.label_27, 5, 0, 1, 1)
        self.cb_Axis_Other = QtWidgets.QComboBox(self.page_3)
        self.cb_Axis_Other.setObjectName("cb_Axis_Other")
        self.gridLayout_6.addWidget(self.cb_Axis_Other, 4, 1, 1, 1)
        self.label_30 = QtWidgets.QLabel(self.page_3)
        self.label_30.setObjectName("label_30")
        self.gridLayout_6.addWidget(self.label_30, 4, 2, 1, 1)
        self.le_Filter_Value_Other = QtWidgets.QLineEdit(self.page_3)
        self.le_Filter_Value_Other.setObjectName("le_Filter_Value_Other")
        self.gridLayout_6.addWidget(self.le_Filter_Value_Other, 5, 4, 1, 2)
        self.ccb_Id_Other = QCheckableComboBox(self.page_3)
        self.ccb_Id_Other.setObjectName("ccb_Id_Other")
        self.gridLayout_6.addWidget(self.ccb_Id_Other, 1, 1, 1, 5)
        self.cb_Type_Other = QtWidgets.QComboBox(self.page_3)
        self.cb_Type_Other.setObjectName("cb_Type_Other")
        self.gridLayout_6.addWidget(self.cb_Type_Other, 4, 3, 1, 2)
        self.label_29 = QtWidgets.QLabel(self.page_3)
        self.label_29.setObjectName("label_29")
        self.gridLayout_6.addWidget(self.label_29, 1, 0, 1, 1)
        self.cb_X = QtWidgets.QComboBox(self.page_3)
        self.cb_X.setObjectName("cb_X")
        self.gridLayout_6.addWidget(self.cb_X, 2, 1, 1, 2)
        self.ccb_Y = QCheckableComboBox(self.page_3)
        self.ccb_Y.setObjectName("ccb_Y")
        self.gridLayout_6.addWidget(self.ccb_Y, 2, 4, 1, 2)
        self.label = QtWidgets.QLabel(self.page_3)
        self.label.setObjectName("label")
        self.gridLayout_6.addWidget(self.label, 2, 3, 1, 1)
        self.le_Filename = QtWidgets.QLineEdit(self.page_3)
        self.le_Filename.setObjectName("le_Filename")
        self.gridLayout_6.addWidget(self.le_Filename, 0, 1, 1, 5)
        self.label_22 = QtWidgets.QLabel(self.page_3)
        self.label_22.setObjectName("label_22")
        self.gridLayout_6.addWidget(self.label_22, 2, 0, 1, 1)
        self.label_25 = QtWidgets.QLabel(self.page_3)
        self.label_25.setObjectName("label_25")
        self.gridLayout_6.addWidget(self.label_25, 0, 0, 1, 1)
        self.label_24 = QtWidgets.QLabel(self.page_3)
        self.label_24.setObjectName("label_24")
        self.gridLayout_6.addWidget(self.label_24, 3, 0, 1, 1)
        self.le_ScaleX_Other = QtWidgets.QLineEdit(self.page_3)
        self.le_ScaleX_Other.setObjectName("le_ScaleX_Other")
        self.gridLayout_6.addWidget(self.le_ScaleX_Other, 3, 1, 1, 3)
        self.label_28 = QtWidgets.QLabel(self.page_3)
        self.label_28.setObjectName("label_28")
        self.gridLayout_6.addWidget(self.label_28, 3, 4, 1, 1)
        self.pb_Open_File_Other = QtWidgets.QPushButton(self.page_3)
        self.pb_Open_File_Other.setMinimumSize(QtCore.QSize(30, 25))
        self.pb_Open_File_Other.setMaximumSize(QtCore.QSize(30, 25))
        self.pb_Open_File_Other.setObjectName("pb_Open_File_Other")
        self.gridLayout_6.addWidget(self.pb_Open_File_Other, 0, 6, 1, 1)
        self.sw_Pages.addWidget(self.page_3)
        self.gridLayout_3.addWidget(self.sw_Pages, 1, 0, 1, 1)
        self.gridLayout.addWidget(self.widget, 0, 0, 1, 1)

        self.retranslateUi(PlotTypeSelector)
        self.sw_Pages.setCurrentIndex(1)
        self.pb_Remove.clicked.connect(PlotTypeSelector.close)
        self.cb_code.currentIndexChanged['int'].connect(self.sw_Pages.setCurrentIndex)
        QtCore.QMetaObject.connectSlotsByName(PlotTypeSelector)

    def retranslateUi(self, PlotTypeSelector):
        _translate = QtCore.QCoreApplication.translate
        PlotTypeSelector.setWindowTitle(_translate("PlotTypeSelector", "Form"))
        self.cb_code.setItemText(0, _translate("PlotTypeSelector", "abci"))
        self.cb_code.setItemText(1, _translate("PlotTypeSelector", "slans"))
        self.cb_code.setItemText(2, _translate("PlotTypeSelector", "other"))
        self.pb_Drop.setText(_translate("PlotTypeSelector", "v"))
        self.pb_Remove.setText(_translate("PlotTypeSelector", "-"))
        self.label_2.setText(_translate("PlotTypeSelector", "request"))
        self.label_9.setText(_translate("PlotTypeSelector", "type"))
        self.label_5.setText(_translate("PlotTypeSelector", "id"))
        self.label_4.setText(_translate("PlotTypeSelector", "pol"))
        self.label_7.setText(_translate("PlotTypeSelector", "scaley"))
        self.cb_Axis_ABCI.setItemText(0, _translate("PlotTypeSelector", "left"))
        self.cb_Axis_ABCI.setItemText(1, _translate("PlotTypeSelector", "right"))
        self.label_8.setText(_translate("PlotTypeSelector", "axis"))
        self.cb_Polarization_ABCI.setItemText(0, _translate("PlotTypeSelector", "long"))
        self.cb_Polarization_ABCI.setItemText(1, _translate("PlotTypeSelector", "trans"))
        self.label_3.setText(_translate("PlotTypeSelector", "folder"))
        self.le_ScaleX_ABCI.setText(_translate("PlotTypeSelector", "1"))
        self.label_6.setText(_translate("PlotTypeSelector", "scalex"))
        self.le_ScaleY_ABCI.setText(_translate("PlotTypeSelector", "1"))
        self.pb_Open_Folder_ABCI.setText(_translate("PlotTypeSelector", "..."))
        self.label_18.setText(_translate("PlotTypeSelector", "scaley"))
        self.label_21.setText(_translate("PlotTypeSelector", "value"))
        self.label_16.setText(_translate("PlotTypeSelector", "axis"))
        self.label_20.setText(_translate("PlotTypeSelector", "type"))
        self.label_12.setText(_translate("PlotTypeSelector", "request"))
        self.label_14.setText(_translate("PlotTypeSelector", "scalex"))
        self.label_19.setText(_translate("PlotTypeSelector", "id"))
        self.label_17.setText(_translate("PlotTypeSelector", "filter"))
        self.label_13.setText(_translate("PlotTypeSelector", "pol"))
        self.label_15.setText(_translate("PlotTypeSelector", "folder"))
        self.pb_Open_Folder_SLANS.setText(_translate("PlotTypeSelector", "..."))
        self.label_31.setText(_translate("PlotTypeSelector", "value"))
        self.label_26.setText(_translate("PlotTypeSelector", "axis"))
        self.label_27.setText(_translate("PlotTypeSelector", "filter"))
        self.label_30.setText(_translate("PlotTypeSelector", "type"))
        self.label_29.setText(_translate("PlotTypeSelector", "id"))
        self.label.setText(_translate("PlotTypeSelector", "y"))
        self.label_22.setText(_translate("PlotTypeSelector", "x"))
        self.label_25.setText(_translate("PlotTypeSelector", "file"))
        self.label_24.setText(_translate("PlotTypeSelector", "scalex"))
        self.label_28.setText(_translate("PlotTypeSelector", "scaley"))
        self.pb_Open_File_Other.setText(_translate("PlotTypeSelector", "..."))

from QCheckableComboBox import QCheckableComboBox
from QSwitchControl import SwitchControl

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    PlotTypeSelector = QtWidgets.QWidget()
    ui = Ui_PlotTypeSelector()
    ui.setupUi(PlotTypeSelector)
    PlotTypeSelector.show()
    sys.exit(app.exec_())

