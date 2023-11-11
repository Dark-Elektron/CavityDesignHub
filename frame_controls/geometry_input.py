from ui_files.geometry_input import Ui_Geometry_Input
from utils.shared_classes import *
from utils.shared_functions import *

fr = FileReader()

file_color = 'green'
DEBUG = True


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class GeometryInputControl:
    def __init__(self, parent):
        self.animation = None
        self.pause_icon = None
        self.resume_icon = None
        self.process_state = None
        self.progress_monitor_thread = None
        self.progress_list = None
        self.shape_space = None
        self.end_routine_thread = None
        self.w_GeometryInput = QWidget()

        self.ui = Ui_Geometry_Input()
        self.ui.setupUi(self.w_GeometryInput)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # Get plot object
        self.geometry_view = self.win.geometryview_widget
        self.plot = self.geometry_view.plot

        # get logger
        self.log = self.main_control.log

        self.initUI()
        self.signals()
        self.filename = None  # placeholder, made a booboo copying the end routine

        # shape space initialization
        self.loaded_shape_space = {}
        self.selected_keys = []
        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False

        # ui effects
        # self.ui_effects()

    def initUI(self):

        # init shape entry mode
        self.shape_entry_widgets_control()

        # inner cell
        self.ui.cb_Inner_Cell.setCheckState(2)
        self.ui.cb_Inner_Cell.setEnabled(False)

        # # create pause and resume1 icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # process state
        self.process_state = 'none'

        # set default boundary condition to magnetic wall at both ends
        # self.ui.cb_LBC.setCurrentIndex(2)
        # self.ui.cb_RBC.setCurrentIndex(2)

    def signals(self):
        # initial push button state
        self.ui.w_Outer_Cell_L.hide()
        self.ui.w_Outer_Cell_R.hide()
        self.ui.w_Expansion.hide()
        # self.ui.pb_Expansion.setEnabled(False)

        self.ui.le_N_Cells.editingFinished.connect(lambda: self.draw_shape_from_shape_space())
        self.ui.le_Scale.editingFinished.connect(lambda: self.draw_shape_from_shape_space())

        # load shape space
        self.ui.pb_Select_Shape_Space.clicked.connect(
            lambda: open_file(self.ui.le_Shape_Space, start_folder=str(self.main_control.projectDir / "Cavities")))
        self.ui.le_Shape_Space.textChanged.connect(
            lambda: load_shape_space(self, self.ui.le_Shape_Space, self.ui.cb_Shape_Space_Keys))

        # control shape entry mode
        self.ui.cb_Shape_Entry_Mode.currentIndexChanged.connect(lambda: self.shape_entry_widgets_control())

        # cell parameters control signals
        self.ui.cb_Outer_Cell_L.stateChanged.connect(lambda: animate_height(
            self.ui.cb_Outer_Cell_L, self.ui.w_Outer_Cell_L, 0, 160, True))
        self.ui.cb_Outer_Cell_R.stateChanged.connect(lambda: animate_height(
            self.ui.cb_Outer_Cell_R, self.ui.w_Outer_Cell_R, 0, 160, True))

        # uncomment to draw again
        self.ui.cb_Shape_Space_Keys.currentTextChanged.connect(lambda: self.draw_shape_from_shape_space())

        #
        self.ui.le_Alpha.editingFinished.connect(lambda: self.update_alpha())

        #
        self.ui.le_Req_i.editingFinished.connect(lambda: self.ui.le_Req_ol.setText(self.ui.le_Req_i.text()))
        self.ui.le_Req_i.editingFinished.connect(lambda: self.ui.le_Req_or.setText(self.ui.le_Req_i.text()))

        # self.ui.le_Scale.editingFinished.connect(lambda: validating(self.ui.le_Scale))

    def shape_entry_widgets_control(self):
        if self.ui.cb_Shape_Entry_Mode.currentIndex() == 0:
            self.ui.w_Enter_Geometry_Manual.setEnabled(False)
            self.ui.w_Select_Shape_Space.show()

            # clear cells from graphics view
            self.plot.ax.clear()
            self.plot.fig.canvas.draw()
        else:

            self.ui.w_Enter_Geometry_Manual.setEnabled(True)
            self.ui.w_Select_Shape_Space.hide()

    def draw_shape_from_shape_space(self):
        ci = 0
        # remove existing cells
        self.plot.ax.clear()
        self.plot.fig.canvas.draw()

        for key in self.loaded_shape_space.keys():
            if self.ui.cb_Shape_Space_Keys.currentText().split(', ').count(key) > 0:
                IC = self.loaded_shape_space[key]["IC"]
                OC = self.loaded_shape_space[key]["OC"]
                if 'OC_R' in self.loaded_shape_space[key].keys():
                    OC_R = self.loaded_shape_space[key]["OC_R"]
                else:
                    OC_R = OC
                BP = self.loaded_shape_space[key]["BP"]
                n_cell = int(self.ui.le_N_Cells.text())
                # boundary conditions
                bc_dict = {'Magnetic Wall En=0': 'b', 'Electric Wall Et=0': 'r',
                           'Inner Contour': 'g', 'Axis': 'k', 'Metal': 'gray'}
                bc = [bc_dict[self.win.eigenmode_widget.ui.cb_LBC.currentText()],
                      bc_dict[self.win.eigenmode_widget.ui.cb_RBC.currentText()]]

                scale = float(self.ui.le_Scale.text())

                plot_cavity_geometry(self.plot, IC, OC, OC_R, BP, n_cell, bc, scale)
                ci += 1
            if ci > 4:  # maximum of only 10 plots
                break

    def prompt(self, fid):
        path = self.main_control.projectDir / fr'SimulationData\SLANS\{fid}'
        # print(path)
        # path = os.path.join(path, fr"{}\{code}\{fid}")
        if os.path.exists(path):
            # print_("Simulation data already exists. Do you want to overwrite it?")
            msg = QMessageBox()
            msg.setWindowTitle("Folder Exist")
            msg.setText("Simulation data already exists. Do you want to overwrite it?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.YesToAll | QMessageBox.Yes | QMessageBox.No | QMessageBox.NoToAll)
            msg.setDefaultButton(QMessageBox.Yes)

            msg.buttonClicked.connect(button_clicked)
            # print_(f'msg: {msg.Yes}')

            x = msg.exec_()

            if x == msg.YesToAll:
                return 'YesToAll'
            if x == msg.Yes:
                return 'Yes'
            if x == msg.No:
                return 'No'
            if x == msg.NoToAll:
                return 'NoToAll'
        else:
            return "Does not exist"

    def update_alpha(self):
        A_i_space = text_to_list(self.ui.le_A_i.text())[0]
        B_i_space = text_to_list(self.ui.le_B_i.text())[0]
        a_i_space = text_to_list(self.ui.le_a_i.text())[0]
        b_i_space = text_to_list(self.ui.le_b_i.text())[0]
        Ri_i_space = text_to_list(self.ui.le_Ri_i.text())[0]
        L_i_space = text_to_list(self.ui.le_L_i.text())[0]
        Req_i_space = text_to_list(self.ui.le_Req_i.text())[0]

        # try:
        alpha_i_space, _ = calculate_alpha(A_i_space, B_i_space, a_i_space, b_i_space,
                                           Ri_i_space, L_i_space, Req_i_space, 0)
        self.ui.le_Alpha.setText(f"{round(alpha_i_space, 2)}")
        # except:
        #     pass

    def serialise(self, state_dict):
        """
        Serialise w_GeometryInput
        Parameters
        ----------
        state_dict: dict
            Dictionary of state of GUI widgets

        Returns
        -------

        """
        serialise(state_dict, self.w_GeometryInput, marker='Geometry')

    def deserialise(self, state_dict):
        """
        Deserialise w_GeometryInput
        Parameters
        ----------
        state_dict: dict
            Dictionary of state of GUI widgets

        Returns
        -------

        """
        deserialise(state_dict, self.w_GeometryInput, marker='Geometry')
