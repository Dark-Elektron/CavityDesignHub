from PyQt5 import QtGui, QtWidgets, QtCore


class HouseKeeping:
    def __init__(self, win):
        self.main_control = win
        self.win = win
        self.ui = win.ui
        self.wakefieldUI = self.main_control.wakefield_widget.ui
        self.tuneUI = self.main_control.tune_widget.ui
        self.eigenmodeUI = self.main_control.eigenmode_widget.ui

        # Signals
        # self.ui.cb_Objective.currentIndexChanged.connect(lambda: self.change_image(self.ui.l_Objective, self.ui.cb_Objective))
        # self.ui.cb_Plot_Type.currentIndexChanged.connect(lambda: self.show_hide())

        # self.ui.cb_Run_Mode_Objectives.currentIndexChanged.connect(lambda: self.show_hide_(self.ui.cb_Run_Mode_Objectives, self.ui.w_No_Of_Processors_Objectives))

        # self.k_loss_polarization_init()
        # self.ui.cb_K_Loss_Polarization.currentIndexChanged.connect(lambda: self.k_loss_polarization_control())

        # self.check_init()

    def k_loss_polarization_init(self):
        self.ui.l_Z_Unit.setText('\u03A9')
        self.ui.l_0707Z_Unit.setText('\u03A9')
        self.ui.dsb_Z.valueChanged.connect(lambda: self.win.calculate_k_loss())
        self.ui.dsb_F0.valueChanged.connect(lambda: self.win.calculate_k_loss())
        self.ui.dsb_Fr.valueChanged.connect(lambda: self.win.calculate_k_loss())
        self.ui.dsb_Fl.valueChanged.connect(lambda: self.win.calculate_k_loss())

    def k_loss_polarization_control(self):
        if self.ui.cb_K_Loss_Polarization.currentText() == 'Longitudinal':
            self.ui.l_Z_Unit.setText('\u03A9')
            self.ui.l_0707Z_Unit.setText('\u03A9')
        else:
            self.ui.l_Z_Unit.setText('\u03A9/m')
            self.ui.l_0707Z_Unit.setText('\u03A9/m')

    def shape_space_check(self):
        if self.ui.cb_Shape_Space.checkState() == 2:
            self.ui.w_Mid_Cells_Manual.hide()
        else:
            self.ui.w_Mid_Cells_Manual.show()

    def show_hide(self):
        widg_dict = {
            "Heat Map": self.ui.w_Heat_Map_Settings,
            "Pairplot": self.ui.w_Pairplot_Settings,
            "Stacked Bar": self.ui.w_Stacked_Bar_Settings,
            "Scatter Plot": self.ui.w_Scatterplot_Settings
        }

        for key, widg in widg_dict.items():
            if key == self.ui.cb_Plot_Type.currentText():
                widg.show()
            else:
                widg.hide()

    def check_init(self):
        self.ui.cb_Others.setCheckState(2)
        self.ui.cb_Shape_Space.setCheckState(2)

    def change_image(self, object, widget):
        if isinstance(object, QtWidgets.QLabel):
            object.setPixmap(QtGui.QPixmap(":/general/images/{}.png".format(widget.currentText().lower())))
        elif isinstance(object, QtWidgets.QPushButton):
            object.setIcon(QtGui.QIcon(":/general/images/{}.png".format(widget.currentText().lower())))

    def show_hide_(self, wid1, wid2):
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()