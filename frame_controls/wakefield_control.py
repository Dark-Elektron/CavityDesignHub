import math
import subprocess
import time
import multiprocessing as mp
from threading import Thread
import pandas as pd
from psutil import NoSuchProcess

from graphics.graphics_view import GraphicsView
from graphics.scene import Scene
from analysis_modules.data_module.abci_data import ABCIData
from analysis_modules.wakefield.ABCI.abci_geometry import ABCIGeometry
from ui_files.wakefield import Ui_Wakefield
from utils.file_reader import FileReader
from utils.shared_classes import *
from utils.shared_functions import *
import scipy.signal as sps

abci_geom = ABCIGeometry()

file_color = 'cyan'


def print_(*arg):
    print(colored(f'\t{arg}', file_color))


class WakefieldControl:
    def __init__(self, parent):
        self.operating_points = pd.DataFrame()
        self.cav_operating_points = {}
        self.end_routine_thread = None
        self.progress_monitor_thread = None
        self.progress_list = None
        self.process_state = None
        self.resume_icon = None
        self.pause_icon = None
        self.w_Wakefield = QWidget()

        self.ui = Ui_Wakefield()
        self.ui.setupUi(self.w_Wakefield)

        self.progress_bar = self.ui.pb_Progress_Bar

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui
        # ###########################

        # get logger
        self.log = self.main_control.log

        # Create Scene
        self.scene = Scene(self)

        # QGraphicsView
        self.graphicsView = GraphicsView(self, 'Wakefield')
        self.ui.vL_2D_Graphics_View.addWidget(self.graphicsView)

        # ##########################
        self.initUI()
        self.signals()
        self.exe_control()
        self.filename = None  # place holder, made a booboo copying the end routine

        # instantiate geometry
        self.abci_geom = ABCIGeometry()

        # shape space initialization
        self.shape_space = {}
        self.loaded_shape_space = {}
        self.selected_keys = []
        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False
        self.animation = None

        # ui effects
        self.ui_effects()

    def signals(self):
        # initial push button state
        self.ui.w_Outer_Cell_L.hide()
        self.ui.w_Outer_Cell_R.hide()
        self.ui.w_Expansion.hide()
        self.ui.pb_Expansion.setEnabled(False)

        # signals
        self.ui.pb_Run.clicked.connect(lambda: self.run_abci())

        # load shape space
        self.ui.pb_Select_Shape_Space.clicked.connect(
            lambda: open_file(self, self.ui.le_Shape_Space, self.ui.cb_Shape_Space_Keys,
                              start_folder=f"{self.main_control.projectDir}/Cavities"))

        self.ui.pb_Load_Machine_Parameters.clicked.connect(
            lambda: self.load_operating_points(self.ui.ccb_Operation_Points))

        # control shape entry mode
        self.ui.cb_Shape_Entry_Mode.currentIndexChanged.connect(lambda: self.shape_entry_widgets_control())

        # cancel
        self.ui.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.ui.pb_Pause_Resume.clicked.connect(
            lambda: self.pause() if self.process_state == 'running' else self.resume())

        #
        self.ui.cb_Shape_Space_Keys.currentTextChanged.connect(lambda: self.draw_shape_from_shape_space())
        self.ui.cb_Shape_Space_Keys.currentTextChanged.connect(
            lambda: self.add_row(self.ui.cb_Shape_Space_Keys.currentText().split(', ')))

        #
        self.ui.le_Req_i.editingFinished.connect(
            lambda: self.ui.le_Req_ol.setText(self.ui.le_Req_i.text()))
        self.ui.le_Req_i.editingFinished.connect(
            lambda: self.ui.le_Req_or.setText(self.ui.le_Req_i.text()))

        # self.ui.pb_Get_Eigenmode_Results.clicked.connect(lambda: self.get_eigenmode_results())

        # self.ui.le_Scale.editingFinished.connect(lambda: validating(self.ui.le_Scale))

    def initUI(self):
        self.ui.w_Save_Folder.setVisible(False)
        self.ui.w_Machine_Parameters.setEnabled(False)
        # splitter
        self.ui.sp_Left_Right_Container.setStretchFactor(1, 3)

        df = write_qtable_to_df(self.ui.tw_Operating_Points_Input)
        # init shape entry mode
        self.shape_entry_widgets_control()

        # disable expansion section for now. Feature to come later
        self.ui.cb_Expansion.setEnabled(False)

        # inner cell
        self.ui.cb_Inner_Cell.setCheckState(2)
        self.ui.cb_Inner_Cell.setEnabled(False)

        self.ui.cb_LBP.setCheckState(2)
        self.ui.cb_RBP.setCheckState(2)
        self.ui.cb_LBP.setEnabled(False)
        self.ui.cb_RBP.setEnabled(False)

        # create pause and resume icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # process state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        self.progress_bar.hide()

        # hide
        self.ui.w_UQ.hide()
        self.ui.w_Mesh.hide()
        self.ui.w_Operating_Points_QOIs.hide()

    def shape_entry_widgets_control(self):
        if self.ui.cb_Shape_Entry_Mode.currentIndex() == 0:
            self.ui.w_Enter_Geometry_Manual.setEnabled(False)
            self.ui.w_Select_Shape_Space.show()

            # clear cells from graphics view
            self.graphicsView.removeCells()
        else:
            self.ui.w_Enter_Geometry_Manual.setEnabled(True)
            self.ui.w_Select_Shape_Space.hide()

            # uncomment following lines to draw
            # clear cells from graphics view
            self.graphicsView.removeCells()

            # draw new cell
            self.graphicsView.drawCells(color=QColor(0, 0, 0, 255))

    def run_abci(self):
        # get analysis parameters
        n_cells = self.ui.sb_N_Cells.value()
        n_cells = text_to_list(self.ui.le_N_Cells.text())
        n_modules = self.ui.sb_N_Modules.value()

        WG_M = self.ui.le_LBP.text()  # half-length of beam pipe between cavities in module

        if WG_M == '':
            WG_M = ['']
        else:
            try:
                WG_M = ast.literal_eval(WG_M) * 0.001
            except ValueError:
                WG_M = ['']

        MROT = self.ui.cb_Polarization_ABCI.currentIndex()
        MT = float(
            self.ui.le_MT.text())  # number of time steps for a beam to move one cell to another default = 3
        bunch_length = float(self.ui.le_Bunch_Length.text())
        NFS = float(self.ui.le_NFS.text())  # Number of samples in FFT (max 10000)
        UBT = float(self.ui.le_Wakelength.text())
        DDZ_SIG = float(self.ui.le_DDZ_SIG.text())
        DDR_SIG = float(self.ui.le_DDR_SIG.text())
        proc_count = self.ui.sb_No_Of_Processors_ABCI.value()
        marker = self.ui.le_Marker.text()

        if self.ui.sc_Operating_Points_QOI.checkState() == 2:
            # qoi_df = write_qtable_to_df(self.ui.tw_Operating_Points_Input)
            qoi_dict = self.get_table_entries()
        else:
            # qoi_df = None
            qoi_dict = None

        # get geometric parameters
        ic(self.loaded_shape_space)
        self.shape_space = get_geometric_parameters(self, 'ABCI', text_to_list(self.ui.le_Scale.text()))
        ic(self.shape_space, self.loaded_shape_space)

        # split shape_space for different processes/ MPI share process by rank
        keys = list(self.shape_space.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

        # show progress bar
        self.progress_bar.show()
        # get the total number of simulations to be run
        if MROT == 2:
            multiplier = 2  # when longitudinal and transverse is run, two simulations are run for each cavity
        else:
            multiplier = 1
        num_sims = len(keys)
        if qoi_dict is not None:
            for k, v in qoi_dict.items():
                num_sims += len(v[0])

        num_sims = multiplier * num_sims * len(n_cells)
        self.progress_bar.setMaximum(num_sims)

        # progress list
        manager = mp.Manager()
        self.progress_list = manager.list()
        self.progress_list.append(0)

        self.processes = []
        for p in range(proc_count):
            # try:
            if p < proc_count - 1:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            processor_shape_space = {}
            for key, val in self.shape_space.items():
                if key in proc_keys_list:
                    processor_shape_space[key] = val

            service = mp.Process(target=self.run_sequential, args=(n_cells, n_modules, processor_shape_space,
                                                                   MROT, MT, NFS, UBT, bunch_length,
                                                                   DDR_SIG, DDZ_SIG,
                                                                   self.main_control.parentDir,
                                                                   self.main_control.projectDir, self.progress_list,
                                                                   WG_M, marker, qoi_dict
                                                                   ))

            service.start()
            self.processes.append(psutil.Process(service.pid))
            self.processes_id.append(service.pid)

            # except Exception as e:
            #     self.log.error(f"Exception in run_MP:: {e}")

        # display progress bar
        self.show_progress_bar = True
        self.progress_monitor_thread = ProgressMonitor(self, self.main_control.projectDir)
        self.progress_monitor_thread.sig.connect(self.update_progress_bar)
        self.progress_monitor_thread.start()

        self.log.info("Wakefield simulation started")
        # change process state to running
        self.process_state = 'running'
        self.run_pause_resume_stop_routine()

        self.end_routine_thread = EndRoutine(self, self.main_control.projectDir)
        self.end_routine_thread.start()

    def add_row(self, cavities):
        if 'All' in cavities or '' in cavities:
            pass
        else:
            # get current number of rows
            n = self.ui.tw_Operating_Points_Input.rowCount()

            # add new row
            self.ui.tw_Operating_Points_Input.setRowCount(n + 1)  # and one row in the table

            keys = list(self.cav_operating_points.keys())

            # add if len of current list greater than prevous list
            if len(cavities) < len(keys):
                for key in keys:
                    if key not in cavities:
                        # get row index
                        row = self.cav_operating_points[key]['row index']
                        self.remove_row(key, row)
            else:
                # remove otherwise
                for cav in cavities:
                    if cav not in self.cav_operating_points.keys():
                        self.create_new_row(n, self.ui.tw_Operating_Points_Input, cav)

    def create_new_row(self, row_ind, table, cavity):
        """

        Parameters
        ----------
        op_dict
        row_ind
        table

        Returns
        -------

        """
        ic(row_ind, cavity, self.ui.tw_Operating_Points_Input.rowCount())
        # Cavity
        l_cavity = QLabel()
        l_cavity.setText(f"{cavity}")
        table.setCellWidget(row_ind, 0, l_cavity)

        # Operating point
        ccb_OperatingPoint = QCheckableComboBox()
        ccb_OperatingPoint.addItem('All')
        for key, val in self.operating_points:
            ccb_OperatingPoint.addItem(f'{key}')

        table.setCellWidget(row_ind, 1, ccb_OperatingPoint)
        # signal to disable polarisation when cb_code item is 'Other'
        self.ui.ccb_Operation_Points.currentTextChanged.connect(
            lambda: self.update_cav_operating_points(ccb_OperatingPoint))

        # # Number of cells
        # sb_n_cells = QSpinBox()
        # sb_n_cells.setMaximum(9999)
        # sb_n_cells.setMinimum(1)
        # table.setCellWidget(row_ind, 2, sb_n_cells)
        #
        # # R/Q
        # dsb_R_Q = QDoubleSpinBox()
        # dsb_R_Q.setMaximum(999999.99)
        # dsb_R_Q.setMinimum(1)
        # table.setCellWidget(row_ind, 3, dsb_R_Q)

        # sigma (SR)
        le_sigma_sr = QLineEdit()
        le_sigma_sr.setReadOnly(True)
        table.setCellWidget(row_ind, 2, le_sigma_sr)
        ccb_OperatingPoint.currentTextChanged.connect(
            lambda: self.update_operating_points_widgets(le_sigma_sr, ccb_OperatingPoint, 'sigma_SR [mm]'))

        # sigma (BS)
        le_sigma_bs = QLineEdit()
        le_sigma_bs.setReadOnly(True)
        table.setCellWidget(row_ind, 3, le_sigma_bs)
        ccb_OperatingPoint.currentTextChanged.connect(
            lambda: self.update_operating_points_widgets(le_sigma_bs, ccb_OperatingPoint, 'sigma_BS [mm]'))

        # I0
        le_I0 = QLineEdit()
        le_I0.setReadOnly(True)
        table.setCellWidget(row_ind, 4, le_I0)
        ccb_OperatingPoint.currentTextChanged.connect(
            lambda: self.update_operating_points_widgets(le_I0, ccb_OperatingPoint, 'I0 [mA]'))

        # Number of bunches Nb
        le_Nb = QLineEdit()
        table.setCellWidget(row_ind, 5, le_Nb)
        ccb_OperatingPoint.currentTextChanged.connect(
            lambda: self.update_operating_points_widgets(le_Nb, ccb_OperatingPoint, 'Nb [1e11]'))

        # # Number of bunches Nb
        # dsb_freq = QDoubleSpinBox()
        # dsb_freq.setMaximum(999999.99)
        # dsb_freq.setMinimum(0.0001)
        # table.setCellWidget(row_ind, 8, dsb_freq)

        args_dict = {'row index': row_ind, 'Cavity': l_cavity, 'Operating Point': ccb_OperatingPoint,
                     'sigma (SR) [mm]': le_sigma_sr, 'sigma (BS) [mm]': le_sigma_bs,
                     'I0 [mA]': le_I0, 'Nb [1e11]': le_Nb}

        self.cav_operating_points[l_cavity.text()] = args_dict

    def remove_row(self, key, row):
        # current number of rows
        n = self.ui.tw_Operating_Points_Input.rowCount()

        # remove row from table
        self.ui.tw_Operating_Points_Input.removeRow(row)

        # remove key from dictionary
        del self.cav_operating_points[key]

        # reset number of rows
        self.ui.tw_Operating_Points_Input.setRowCount(n - 2)

        # adjust other row elements
        d = self.cav_operating_points
        tab_col = ['Cavity', 'Operating Point', 'N Cells', 'R/Q [Ohm]',
                   'sigma (SR) [mm]', 'sigma (BS) [mm]',
                   'I0 [mA]', 'Nb [1e11]', 'freq [MHz]']

        for i, (key, val) in enumerate(d.items()):
            for j, tc in enumerate(tab_col):
                self.ui.tw_Operating_Points_Input.setCellWidget(i, j, val[tc])

            # update row index
            val['row index'] = i

        self.cav_operating_points = d

    def load_operating_points(self, ccb):
        filename, _ = QFileDialog.getOpenFileName(
            None, "Open File", fr"{self.main_control.projectDir}/Cavities", "Json Files (*.json)")
        try:
            self.ui.le_Machine_Parameters_Filename.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            # populate checkboxes with key
            ccb.clear()
            ccb.addItem("All")
            for col in dd.keys():
                self.ui.ccb_Operation_Points.addItem(fr'{col}')

            self.operating_points = pd.DataFrame(dd)
            ic(self.operating_points)

        except Exception as e:
            print('Failed to open file:: ', e)

    def update_operating_points_widgets(self, le, ccb, var):
        op_points = ccb.currentText().split(', ')
        if 'All' in op_points or '' in op_points:
            pass
        else:
            txt = list(self.operating_points.loc[var, op_points])
            le.setText(str(txt))

    def update_cav_operating_points(self, ccb):
        # populate checkboxes with key
        ccb.clear()
        ccb.addItem("All")
        for k in self.ui.ccb_Operation_Points.currentText().split(', '):
            if k == 'All' or k == '':
                pass
            else:
                ccb.addItem(fr'{k}')

    def get_eigenmode_results(self):
        for key in self.cav_operating_points.keys():
            ic(key)
            if len(self.ui.le_N_Cells) == 1:
                try:
                    with open(fr'{self.main_control.projectDir}\SimulationData\SLANS\{key}\qois.json') as json_file:
                        qois_slans = json.load(json_file)
                    n_cells = qois_slans['N Cells']
                    R_Q = qois_slans['R/Q [Ohm]']
                    freq = qois_slans['freq [MHz]']
                    print("Found a corresponding SLANS file")
                    # self.cav_operating_points[key]['N Cells'].setValue(int(n_cells))
                    # self.cav_operating_points[key]['R/Q [Ohm]'].setValue(R_Q)
                    # self.cav_operating_points[key]['freq [MHz]'].setValue(freq)
                except FileNotFoundError as e:
                    # Run eigenmode analysis
                    print("Did not find a corresponding SLANS file. Please check the R/Q used to calculate the"
                          "k_FM and P_HOM", e)

    def get_table_entries(self):
        dd = {}
        for key, val in self.cav_operating_points.items():
            for scale in text_to_list(self.ui.le_Scale.text()):
                for n_cell in text_to_list(self.ui.le_N_Cells.text()):

                    if scale == 1 or scale == 0:
                        if len(text_to_list(self.ui.le_N_Cells.text())) == 1:
                            folder_name = key
                        else:
                            folder_name = f"{key}_{n_cell}"
                    else:
                        if len(text_to_list(self.ui.le_N_Cells.text())) == 1:
                            folder_name = f"{key}_scale_{scale}"
                        else:
                            folder_name = f"{key}_scale_{scale}_{n_cell}"

                    with open(fr'{self.main_control.projectDir}\SimulationData\SLANS\{folder_name}\qois.json') as json_file:
                        qois_slans = json.load(json_file)

                    op_points = val['Operating Point'].currentText().split(', ')
                    n_cells = qois_slans['N Cells']
                    R_Q = qois_slans['R/Q [Ohm]']
                    sigma_SR = ast.literal_eval(val['sigma (SR) [mm]'].text())
                    sigma_BS = ast.literal_eval(val['sigma (BS) [mm]'].text())
                    I0 = ast.literal_eval(val['I0 [mA]'].text())
                    Nb = ast.literal_eval(val['Nb [1e11]'].text())
                    # freq = val['freq [MHz]'].value()
                    freq = qois_slans['freq [MHz]']
                    dd.update({folder_name: [op_points, n_cells, R_Q, sigma_SR, sigma_BS, I0, Nb, freq]})

        ic(dd)
        return dd

    def run_pause_resume_stop_routine(self):
        if self.process_state == 'none':
            # change pause/resume icon to pause icon
            self.ui.pb_Pause_Resume.setIcon(self.pause_icon)

            # disable pause/resume and cancel buttons
            self.ui.pb_Pause_Resume.setEnabled(False)
            self.ui.pb_Cancel.setEnabled(False)

            # enable run button in case it was disabled
            self.ui.pb_Run.setEnabled(True)

        if self.process_state == "running":
            # enable run, pause/resume and cancel buttons
            self.ui.pb_Pause_Resume.setEnabled(True)
            self.ui.pb_Cancel.setEnabled(True)
            self.ui.pb_Run.setEnabled(False)

            # change pause/resume icon to pause icon
            self.ui.pb_Pause_Resume.setIcon(self.pause_icon)

        if self.process_state == 'paused':
            # disable run button
            self.ui.pb_Run.setEnabled(False)

            # change pause/resume button icon to resume icon
            self.ui.pb_Pause_Resume.setIcon(self.resume_icon)

    def pause(self):
        # self.log.info("Pausing...")
        for p in self.processes:
            p.suspend()
        self.log.info("Paused")

        self.process_state = 'paused'
        self.run_pause_resume_stop_routine()

    def resume(self):
        # self.log.info("Resuming...")
        for p in self.processes:
            p.resume()
        self.log.info("Resumed")

        self.process_state = 'running'
        self.run_pause_resume_stop_routine()

    def cancel(self):
        self.log.info("Terminating process...")
        # signal to progress bar
        self.show_progress_bar = False

        try:
            for p in self.processes:
                p.terminate()
        except NoSuchProcess:
            pass

        self.processes.clear()
        self.processes_id.clear()

        self.process_state = 'none'
        self.run_pause_resume_stop_routine()
        self.log.info("Process terminated.")

    def end_routine(self, proc_ids):
        for pid in proc_ids:
            # try:
            p = psutil.Process(pid)
            while p.is_running():
                pass
            print(fr"process {p} ended")
            # except Exception as e:
            #     print_("Exception:: ", e)

        self.cancel()

    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

        if val == 100 or not self.show_progress_bar:
            # reset progress bar
            self.progress_bar.setValue(0)
            self.progress_bar.hide()

    def prompt(self, code, fid):
        # path = os.getcwd()
        # path = os.path.join(path, fr"File\{code}\{fid}")
        path = fr'{self.main_control.projectDir}\SimulationData\ABCI\{fid}'
        if os.path.exists(path):
            print_("File already exists. Do you want to overwrite it?")
            msg = QMessageBox()
            msg.setWindowTitle("Folder Exist")
            msg.setText("File already exists. Do you want to overwrite it?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.YesToAll | QMessageBox.Yes | QMessageBox.No | QMessageBox.NoToAll)
            msg.setDefaultButton(QMessageBox.Yes)

            msg.buttonClicked.connect(button_clicked)

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
            return 'YesToAll'

    @staticmethod
    def show_hide_(wid1, wid2):
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()

    def exe_control(self):
        # Abci
        self.ui.pb_Top_Drawer.clicked.connect(lambda: self.run_abci_exe(
            fr'{self.main_control.parentDir}\exe\ABCI_exe\TopDrawer for Windows\TopDrawW.exe'))

    def ui_effects(self):
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Settings.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Inner_Cell.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Outer_Cell_L.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Outer_Cell_R.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Expansion.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Simulation_Controls.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Show_Cavity.setGraphicsEffect(shadow)
        #
        # shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        # shadow.setColor(QColor(0, 0, 0, 77))
        # self.ui.w_Load_Manual.setGraphicsEffect(shadow)
        pass

    def serialize(self, state_dict):
        # update state file
        state_dict["Wakefield_Shape_Entry_Mode"] = self.ui.cb_Shape_Entry_Mode.currentIndex()
        state_dict["Wakefield_Shape_Space"] = self.ui.le_Shape_Space.text()
        state_dict["Wakefield_Mid_Cell_CB"] = self.ui.cb_Inner_Cell.checkState()
        state_dict["Wakefield_Left_Cell_CB"] = self.ui.cb_Outer_Cell_L.checkState()
        state_dict["Wakefield_Right_Cell_CB"] = self.ui.cb_Outer_Cell_R.checkState()
        state_dict["Wakefield_Expansion_CB"] = self.ui.cb_Expansion.checkState()
        state_dict["Wakefield_LBP_CB"] = self.ui.cb_LBP.checkState()
        state_dict["Wakefield_RBP_CB"] = self.ui.cb_RBP.checkState()

        # cell parameters
        state_dict["Wakefield_A_i"] = self.ui.le_A_i.text()
        state_dict["Wakefield_B_i"] = self.ui.le_B_i.text()
        state_dict["Wakefield_a_i"] = self.ui.le_a_i.text()
        state_dict["Wakefield_b_i"] = self.ui.le_b_i.text()
        state_dict["Wakefield_Ri_i"] = self.ui.le_Ri_i.text()
        state_dict["Wakefield_L_i"] = self.ui.le_L_i.text()
        state_dict["Wakefield_Req_i"] = self.ui.le_Req_i.text()
        state_dict["Wakefield_Alpha_i"] = self.ui.le_Alpha.text()

        state_dict["Wakefield_A_ol"] = self.ui.le_A_ol.text()
        state_dict["Wakefield_B_ol"] = self.ui.le_B_ol.text()
        state_dict["Wakefield_a_ol"] = self.ui.le_a_ol.text()
        state_dict["Wakefield_b_ol"] = self.ui.le_b_ol.text()
        state_dict["Wakefield_Ri_ol"] = self.ui.le_Ri_ol.text()
        state_dict["Wakefield_L_ol"] = self.ui.le_L_ol.text()
        state_dict["Wakefield_Req_ol"] = self.ui.le_Req_ol.text()
        state_dict["Wakefield_Alpha_ol"] = self.ui.le_Alpha_ol.text()

        state_dict["Wakefield_A_or"] = self.ui.le_A_or.text()
        state_dict["Wakefield_B_or"] = self.ui.le_B_or.text()
        state_dict["Wakefield_a_or"] = self.ui.le_a_or.text()
        state_dict["Wakefield_b_or"] = self.ui.le_b_or.text()
        state_dict["Wakefield_Ri_or"] = self.ui.le_Ri_or.text()
        state_dict["Wakefield_L_or"] = self.ui.le_L_or.text()
        state_dict["Wakefield_Req_or"] = self.ui.le_Req_or.text()
        state_dict["Wakefield_Alpha_or"] = self.ui.le_Alpha_or.text()

        # settings
        state_dict["Wakefield_N_Cells"] = self.ui.sb_N_Cells.value()
        state_dict["Wakefield_N_Cells_"] = self.ui.le_N_Cells.text()
        state_dict["Wakefield_Scale"] = self.ui.le_Scale.text()
        state_dict["Wakefield_N_Modules"] = self.ui.sb_N_Modules.value()
        state_dict["Wakefield_Polarization"] = self.ui.cb_Polarization_ABCI.currentIndex()
        state_dict["Wakefield_No_Of_Processors"] = self.ui.sb_No_Of_Processors_ABCI.value()
        state_dict["Wakefield_Wakelength"] = self.ui.le_Wakelength.text()
        state_dict["Wakefield_Bunch_Length"] = self.ui.le_Bunch_Length.text()
        state_dict["Wakefield_NFS"] = self.ui.le_NFS.text()
        state_dict["Wakefield_MT"] = self.ui.le_MT.text()
        state_dict["Wakefield_DDz/SIG"] = self.ui.le_DDZ_SIG.text()
        state_dict["Wakefield_DDr/SIG"] = self.ui.le_DDR_SIG.text()

    def deserialize(self, state_dict):
        # update state file
        self.ui.cb_Shape_Entry_Mode.setCurrentIndex(state_dict["Wakefield_Shape_Entry_Mode"])
        self.ui.le_Shape_Space.setText(state_dict["Wakefield_Shape_Space"])
        self.ui.cb_Inner_Cell.setCheckState(state_dict["Wakefield_Mid_Cell_CB"])
        self.ui.cb_Outer_Cell_L.setCheckState(state_dict["Wakefield_Left_Cell_CB"])
        self.ui.cb_Outer_Cell_R.setCheckState(state_dict["Wakefield_Right_Cell_CB"])
        self.ui.cb_Expansion.setCheckState(state_dict["Wakefield_Expansion_CB"])
        self.ui.cb_LBP.setCheckState(state_dict["Wakefield_LBP_CB"])
        self.ui.cb_RBP.setCheckState(state_dict["Wakefield_RBP_CB"])

        # cell parameters
        self.ui.le_A_i.setText(state_dict["Wakefield_A_i"])
        self.ui.le_B_i.setText(state_dict["Wakefield_B_i"])
        self.ui.le_a_i.setText(state_dict["Wakefield_a_i"])
        self.ui.le_b_i.setText(state_dict["Wakefield_b_i"])
        self.ui.le_Ri_i.setText(state_dict["Wakefield_Ri_i"])
        self.ui.le_L_i.setText(state_dict["Wakefield_L_i"])
        self.ui.le_Req_i.setText(state_dict["Wakefield_Req_i"])
        self.ui.le_Alpha.setText(state_dict["Wakefield_Alpha_i"])

        self.ui.le_A_ol.setText(state_dict["Wakefield_A_ol"])
        self.ui.le_B_ol.setText(state_dict["Wakefield_B_ol"])
        self.ui.le_a_ol.setText(state_dict["Wakefield_a_ol"])
        self.ui.le_b_ol.setText(state_dict["Wakefield_b_ol"])
        self.ui.le_Ri_ol.setText(state_dict["Wakefield_Ri_ol"])
        self.ui.le_L_ol.setText(state_dict["Wakefield_L_ol"])
        self.ui.le_Req_ol.setText(state_dict["Wakefield_Req_ol"])
        self.ui.le_Alpha_ol.setText(state_dict["Wakefield_Alpha_ol"])

        self.ui.le_A_or.setText(state_dict["Wakefield_A_or"])
        self.ui.le_B_or.setText(state_dict["Wakefield_B_or"])
        self.ui.le_a_or.setText(state_dict["Wakefield_a_or"])
        self.ui.le_b_or.setText(state_dict["Wakefield_b_or"])
        self.ui.le_Ri_or.setText(state_dict["Wakefield_Ri_or"])
        self.ui.le_L_or.setText(state_dict["Wakefield_L_or"])
        self.ui.le_Req_or.setText(state_dict["Wakefield_Req_or"])
        self.ui.le_Alpha_or.setText(state_dict["Wakefield_Alpha_or"])

        # settings
        self.ui.sb_N_Cells.setValue(state_dict["Wakefield_N_Cells"])
        self.ui.le_N_Cells.setText(state_dict["Wakefield_N_Cells_"])
        self.ui.le_Scale.setText(state_dict["Wakefield_Scale"])
        self.ui.sb_N_Modules.setValue(state_dict["Wakefield_N_Modules"])
        self.ui.cb_Polarization_ABCI.setCurrentIndex(state_dict["Wakefield_Polarization"])
        self.ui.sb_No_Of_Processors_ABCI.setValue(state_dict["Wakefield_No_Of_Processors"])
        self.ui.le_Wakelength.setText(state_dict["Wakefield_Wakelength"])
        self.ui.le_Bunch_Length.setText(state_dict["Wakefield_Bunch_Length"])
        self.ui.le_NFS.setText(state_dict["Wakefield_NFS"])
        self.ui.le_MT.setText(state_dict["Wakefield_MT"])
        self.ui.le_DDZ_SIG.setText(state_dict["Wakefield_DDz/SIG"])
        self.ui.le_DDR_SIG.setText(state_dict["Wakefield_DDr/SIG"])

    @staticmethod
    def run_abci_exe(path):
        path = os.path.join(os.getcwd(), path)
        t = Thread(target=subprocess.call, args=(path,))
        t.start()

    def draw_shape_from_shape_space(self):
        colors = [[48, 162, 218, 255], [252, 79, 48, 255], [229, 174, 56, 255], [109, 144, 79, 255],
                  [139, 139, 139, 255]]
        ci = 0

        # remove existing cells
        self.graphicsView.removeCells()
        for key in self.loaded_shape_space.keys():
            if key in self.ui.cb_Shape_Space_Keys.currentText():
                IC = self.loaded_shape_space[key]["IC"]
                OC = self.loaded_shape_space[key]["OC"]
                BP = self.loaded_shape_space[key]["BP"]
                self.graphicsView.drawCells(IC, OC, BP,
                                            QColor(colors[ci][0], colors[ci][1], colors[ci][2], colors[ci][3]))

                ci += 1
            if ci > 4:  # maximum of only 10 plots
                break

    @staticmethod
    def run_sequential(n_cells, n_modules, processor_shape_space,
                       MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
                       DDR_SIG=0.1, DDZ_SIG=0.1,
                       parentDir=None, projectDir=None, progress_list=None,
                       WG_M=None, marker='', qoi_df=None):
        progress = 0
        # get length of processor
        total_no_of_shapes = len(list(processor_shape_space.keys()))
        for key, shape in processor_shape_space.items():
            # run abci code
            start_time = time.time()
            # run both polarizations if MROT == 2
            for ii in WG_M:
                if "OC_R" in list(shape.keys()):
                    OC_R = "OC_R"
                else:
                    OC_R = "OC"
                for n_cell in n_cells:
                    if len(n_cells) == 1:
                        fid = key
                    else:
                        fid = f"{key}_{n_cell}"

                    if MROT == 2:
                        for m in range(2):
                            abci_geom.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape[OC_R],
                                             fid=fid, MROT=m, MT=MT, NFS=NFS, UBT=UBT,
                                             bunch_length=bunch_length,
                                             DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                             projectDir=projectDir,
                                             WG_M=ii, marker=ii)

                    else:
                        abci_geom.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                         fid=f"{key}_{n_cell}", MROT=MROT, MT=MT, NFS=NFS, UBT=UBT,
                                         bunch_length=bunch_length,
                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir,
                                         WG_M=ii, marker=ii)

                    # update progress
                    progress_list.append(progress + 1)

                    print_(f'Cavity {fid}. Time: {time.time() - start_time}')
                    ic(qoi_df)
                    if qoi_df is not None:
                        d = {}
                        ic(qoi_df[fid])
                        #  # save qois
                    # for indx, val in qoi_df[key].items():
                        op_points, no_of_cells, R_Q, sigma_SR_list, sigma_BS_list, I0_list, Nb_list, freq = qoi_df[fid]

                        for i, op_point in enumerate(op_points):
                            WP = op_point
                            I0 = I0_list[i]
                            Nb = Nb_list[i]
                            sigma_z = [sigma_SR_list[i], sigma_BS_list[i]]
                            freq = freq
                            ic(freq)

                            bl_diff = ['SR', 'BS']

                            for j, s in enumerate(sigma_z):
                                for ii in WG_M:
                                    fid_op = f"{WP}_{bl_diff[j]}_{s}mm{ii}"
                                    if "OC_R" in list(shape.keys()):
                                        OC_R = "OC_R"
                                    else:
                                        OC_R = "OC"
                                    for m in range(2):
                                        abci_geom.cavity(no_of_cells, n_modules, shape['IC'], shape['OC'],
                                                         shape['OC_R'],
                                                         fid=fid, MROT=m, MT=MT, NFS=NFS, UBT=10 * s * 1e-3,
                                                         bunch_length=s,
                                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                                         projectDir=projectDir,
                                                         WG_M=ii, marker=ii, sub_dir=fid)

                                        dirc = fr'{projectDir}\SimulationData\ABCI\{fid}{marker}'
                                        # try:
                                        k_loss = abs(ABCIData(dirc, f'{fid_op}', 0).loss_factor['Longitudinal'])
                                        k_kick = abs(ABCIData(dirc, f'{fid_op}', 1).loss_factor['Transverse'])
                                        # except:
                                        #     k_loss = 0
                                        #     k_kick = 0

                                        d[fid_op] = get_qois_value(freq, R_Q, k_loss, k_kick, s, I0, Nb, no_of_cells)

                                        # update progress
                                        progress_list.append(progress + 1)

                        # save qoi dictionary
                        run_save_directory = fr'{projectDir}\SimulationData\ABCI\{fid}{marker}'
                        with open(fr'{run_save_directory}\qois.json', "w") as f:
                            json.dump(d, f, indent=4, separators=(',', ': '))

            print_("Done with the secondary analysis for working points")

    @staticmethod
    def uq(shape_space, objectives, solver_dict, solver_args_dict):
        for key, shape in shape_space.items():
            err = False
            result_dict_slans, result_dict_abci = {}, {}
            run_slans, run_abci = False, False
            slans_obj_list, abci_obj_list = [], []
            for o in objectives:

                if o[1] in ["Req", "freq", "Q", "E", "R/Q", "Epk/Eacc", "Bpk/Eacc"]:
                    result_dict_slans[o[1]] = {'expe': [], 'stdDev': []}
                    run_slans = True
                    slans_obj_list.append(o)

                if o[1].split(' ')[0] in ['ZL', 'ZT', 'k_loss', 'k_kick']:
                    # ic(o)
                    result_dict_abci[o[1]] = {'expe': [], 'stdDev': []}
                    run_abci = True
                    abci_obj_list.append(o)

            # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
            p_true = shape['IC'][0:5]
            # ic(p_true)
            rdim = len(p_true)  # How many variabels will be considered as random in our case 5
            degree = 1

            #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
            nodes = np.array(0)  # initialization
            weights = np.array(0)  # initialization
            flag_stroud = 1
            if flag_stroud == 1:
                nodes, weights, bpoly = quad_stroud3(rdim, degree)
                nodes = 2. * nodes - 1.
            elif flag_stroud == 2:
                nodes, weights, bpoly = quad_stroud3(rdim, degree)  # change to stroud 5 later
                nodes = 2. * nodes - 1.
            else:
                ic('flag_stroud==1 or flag_stroud==2')

            #  mean value of geometrical parameters
            p_init = np.zeros(np.shape(p_true))

            no_parm, no_sims = np.shape(nodes)
            # ic(no_sims)
            delta = 0.05  # or 0.1

            if run_abci:
                # ic("here in ANCI UQ")
                Ttab_val_f = []
                solver, solver_args = solver_dict['abci'], solver_args_dict['abci']
                n_cells = solver_args['n_cells']
                n_modules = solver_args['n_modules']
                MROT = solver_args['MROT']
                MT = solver_args['MT']
                NFS = solver_args['NFS']
                UBT = solver_args['UBT']
                bunch_length = solver_args['bunch_length']
                DDR_SIG = solver_args['DDR_SIG']
                DDZ_SIG = solver_args['DDZ_SIG']
                parentDir = solver_args['parentDir']
                projectDir = solver_args['projectDir']
                progress_list = solver_args['progress_list']
                WG_M = solver_args['WG_M']
                marker = solver_args['marker']

                proc = solver_args['proc']
                sub_dir = fr'{key}'  # the simulation runs at the quadrature points
                # are saved to the key of the mean value run
                no_error = True
                for i in range(no_sims):
                    skip = False
                    p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
                    p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
                    p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
                    p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
                    p_init[4] = p_true[4] * (1 + delta * nodes[4, i])

                    par_mid = np.append(p_init, shape['IC'][5:]).tolist()
                    par_end = par_mid

                    ok = perform_geometry_checks(par_mid, par_end)
                    if not ok:
                        no_error = False
                        break

                    fid = fr'{key}_Q{i}'

                    # check if folder exists and skip if it does
                    if os.path.exists(fr'{projectDir}\SimulationData\ABCI\{key}\{fid}'):
                        skip = True

                    if not skip:
                        #  run your model using SLANC or CST
                        # # create folders for all keys
                        solver.createFolder(fid, projectDir, subdir=sub_dir)
                        for wi in range(MROT):
                            solver.cavity(n_cells, n_modules, par_mid, par_end, par_end, fid=fid, MROT=wi,
                                          DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, beampipes=None, bunch_length=bunch_length,
                                          MT=MT, NFS=NFS, UBT=UBT,
                                          parentDir=parentDir, projectDir=projectDir, WG_M='',
                                          marker='', sub_dir=sub_dir
                                          )

                    # get objective function values
                    abci_folder = fr'{projectDir}\SimulationData\ABCI\{key}'
                    if os.path.exists(abci_folder):
                        # ic(abci_obj_list)
                        obj_result = get_wakefield_objectives_value(fid, abci_obj_list, abci_folder)
                        # ic(obj_result)

                        tab_val_f = obj_result
                        if 'error' in obj_result:
                            no_error = False
                            ic(obj_result)
                            ic("Encountered an error")
                            break
                        Ttab_val_f.append(tab_val_f)
                    else:
                        no_error = False

                if no_error:
                    v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
                    # append results to dict
                    # ic(v_expe_fobj, v_stdDev_fobj)
                    for i, o in enumerate(abci_obj_list):
                        result_dict_abci[o[1]]['expe'].append(v_expe_fobj[i])
                        result_dict_abci[o[1]]['stdDev'].append(v_stdDev_fobj[i])

                    with open(fr"{projectDir}\SimulationData\ABCI\{key}\uq.json", 'w') as file:
                        file.write(json.dumps(result_dict_abci, indent=4, separators=(',', ': ')))


def get_qois_value(f_fm, R_Q, k_loss, k_kick, sigma_z, I0, Nb, n_cell):
    c = 299792458
    w_fm = 2 * np.pi * f_fm * 1e6
    e = 1.602e-19

    k_fm = (w_fm / 4) * R_Q * np.exp(-(w_fm * sigma_z * 1e-3 / c) ** 2) * 1e-12
    k_hom = k_loss - k_fm
    p_hom = (k_hom * 1e12) * (I0 * 1e-3) * e * (Nb * 1e11)

    d = {
        "n cell": n_cell,
        "freq [MHz]": f_fm,
        "R/Q [Ohm]": R_Q,
        "k_FM [V/pC]": k_fm,
        "I0 [mA]": I0,
        "sigma_z [mm]": sigma_z,
        "Nb [1e11]": Nb,
        "|k_loss| [V/pC]": k_loss,
        "|k_kick| [V/pC/m]": k_kick,
        "P_HOM [kW]": p_hom * 1e-3
    }
    return d


def get_wakefield_objectives_value(key, obj, abci_data_dir):
    k_loss_transverse = []
    k_loss_longitudinal = []
    k_loss_M0 = []
    key_list = []

    # create list to hold Z
    Zmax_mon_list = []
    Zmax_dip_list = []
    xmax_mon_list = []
    xmax_dip_list = []
    processed_keys = []

    def calc_k_loss():
        print(f"Processing for Cavity {key}")
        abci_data_long = ABCIData(abci_data_dir, key, 0)
        abci_data_trans = ABCIData(abci_data_dir, key, 1)

        # trans
        x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
        k_loss_trans = abci_data_trans.loss_factor['Transverse']

        if math.isnan(k_loss_trans):
            print_(f"Encountered an exception: Check shape {key}")
            return [0, 0, 0]

        # long
        x, y, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
        abci_data_long.get_data('Loss Factor Spectrum Integrated up to F')

        k_M0 = abci_data_long.y_peaks[0]
        k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
        k_loss_HOM = k_loss_long - k_M0

        # append only after successful run
        k_loss_M0.append(k_M0)
        k_loss_longitudinal.append(k_loss_HOM)
        k_loss_transverse.append(k_loss_trans)

        return [k_loss_M0, k_loss_longitudinal, k_loss_transverse]

    def get_Zmax_L(mon_interval=None):
        if mon_interval is None:
            mon_interval = [0.0, 2e10]

        print(f"Processing for Cavity {key}")
        try:
            abci_data_mon = ABCIData(abci_data_dir, f"{key}", 0)

            # get longitudinal and transverse impedance plot data
            xr_mon, yr_mon, _ = abci_data_mon.get_data('Real Part of Longitudinal Impedance')
            xi_mon, yi_mon, _ = abci_data_mon.get_data('Imaginary Part of Longitudinal Impedance')

            # Zmax
            if mon_interval is None:
                mon_interval = [[0.0, 10]]

            # calculate magnitude
            ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]

            # get peaks
            peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
            xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

            for i, z_bound in enumerate(mon_interval):
                # get mask
                msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

                if len(yp_mon[msk_mon]) != 0:
                    Zmax_mon = max(yp_mon[msk_mon])

                    Zmax_mon_list[i].append(Zmax_mon)
                elif len(yp_mon) != 0:
                    Zmax_mon_list[i].append(0)
                else:
                    return ['error']

            processed_keys.append(key)
        except:
            return ['error']

        return Zmax_mon_list

    def get_Zmax_T(dip_interval=None):
        if dip_interval is None:
            dip_interval = [0.0, 2e10]

        try:
            print(f"Processing for Cavity {key}")
            abci_data_dip = ABCIData(abci_data_dir, f"{key}", 1)

            xr_dip, yr_dip, _ = abci_data_dip.get_data('Real Part of Transverse Impedance')
            xi_dip, yi_dip, _ = abci_data_dip.get_data('Imaginary Part of Transverse Impedance')

            # Zmax
            if dip_interval is None:
                dip_interval = [[0.0, 10]]

            # calculate magnitude
            ymag_dip = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_dip, yi_dip)]

            # get peaks
            peaks_dip, _ = sps.find_peaks(ymag_dip, height=0)
            xp_dip, yp_dip = np.array(xr_dip)[peaks_dip], np.array(ymag_dip)[peaks_dip]

            for i, z_bound in enumerate(dip_interval):
                # get mask
                msk_dip = [(z_bound[0] < x < z_bound[1]) for x in xp_dip]

                if len(yp_dip[msk_dip]) != 0:
                    Zmax_dip = max(yp_dip[msk_dip])

                    Zmax_dip_list[i].append(Zmax_dip)
                elif len(yp_dip) != 0:
                    Zmax_dip_list[i].append(0)
                else:
                    return ['error']

            processed_keys.append(key)
        except:
            return ['error']

        return Zmax_dip_list

    def all_(mon_interval, dip_interval):
        print(f"Processing for Cavity {key}")
        abci_data_long = ABCIData(abci_data_dir, f"{key}_", 0)
        abci_data_trans = ABCIData(abci_data_dir, f"{key}_", 1)

        # get longitudinal and transverse impedance plot data
        xr_mon, yr_mon, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
        xi_mon, yi_mon, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')

        xr_dip, yr_dip, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
        xi_dip, yi_dip, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')

        # loss factors
        # trans
        k_loss_trans = abci_data_trans.loss_factor['Transverse']

        if math.isnan(k_loss_trans):
            print_(f"Encountered an exception: Check shape {key}")
            return 0

        # long
        abci_data_long.get_data('Loss Factor Spectrum Integrated upto F')

        k_M0 = abci_data_long.y_peaks[0]
        k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
        k_loss_HOM = k_loss_long - k_M0

        # calculate magnitude
        ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]
        ymag_dip = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_dip, yi_dip)]

        # get peaks
        peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
        xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

        peaks_dip, _ = sps.find_peaks(ymag_dip, height=0)
        xp_dip, yp_dip = np.array(xr_dip)[peaks_dip], np.array(ymag_dip)[peaks_dip]

        for ii, z_bound in enumerate(mon_interval):
            # get mask
            msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

            if len(yp_mon[msk_mon]) != 0:
                Zmax_mon = max(yp_mon[msk_mon])
                xmax_mon = xp_mon[np.where(yp_mon == Zmax_mon)][0]

                Zmax_mon_list[ii].append(Zmax_mon)
                xmax_mon_list[ii].append(xmax_mon)
            elif len(yp_mon) != 0:
                Zmax_mon_list[ii].append(0.0)
                xmax_mon_list[ii].append(0.0)
            else:
                continue

        for i, z_bound in enumerate(dip_interval):
            # get mask
            msk_dip = [(z_bound[0] < x < z_bound[1]) for x in xp_dip]

            if len(yp_dip[msk_dip]) != 0:
                Zmax_dip = max(yp_dip[msk_dip])
                xmax_dip = xp_dip[np.where(yp_dip == Zmax_dip)][0]

                Zmax_dip_list[i].append(Zmax_dip)
                xmax_dip_list[i].append(xmax_dip)
            elif len(yp_dip) != 0:
                Zmax_dip_list[i].append(0.0)
                xmax_dip_list[i].append(0.0)
            else:
                continue

        # append only after successful run

        k_loss_M0.append(k_M0)
        k_loss_longitudinal.append(k_loss_HOM)
        k_loss_transverse.append(k_loss_trans)

    ZL, ZT = [], []
    freq_range_ZL, freq_range_ZT = [], []

    for i, o in enumerate(obj):

        if o[1].split(' ')[0] == 'ZL':
            freq_range_ZL.append(o[2])
        elif o[1].split(' ')[0] == 'ZT':
            freq_range_ZT.append(o[2])

        elif o[1] == "k_loss":
            pass
        elif o[1] == "k_kick":
            pass

    # ic("about to evaluate ZL", freq_range_ZL)
    if freq_range_ZL:
        for i in range(len(freq_range_ZL)):
            Zmax_mon_list.append([])
            xmax_mon_list.append([])

        ZL = get_Zmax_L(freq_range_ZL)

    if freq_range_ZT:
        for i in range(len(freq_range_ZT)):
            Zmax_dip_list.append([])
            xmax_dip_list.append([])

        ZT = get_Zmax_T(freq_range_ZT)

    ZL, ZT = np.array(ZL).T, np.array(ZT).T
    # ic(ZL, ZT)

    if ZL.size != 0 and ZT.size != 0:
        obj_result = np.hstack((ZL, ZT))
    elif ZL.size != 0:
        obj_result = ZL
    else:
        obj_result = ZT

    return list(obj_result[0])


def process_interval(interval_list):
    interval = []
    for i in range(len(interval_list) - 1):
        interval.append([interval_list[i], interval_list[i + 1]])

    return interval


def write_qtable_to_df(table):
    col_count = table.columnCount()
    row_count = table.rowCount()
    headers = [str(table.horizontalHeaderItem(i).text()) for i in range(col_count)]

    # df indexing is slow, so use lists
    df_list = []
    for row in range(row_count):
        df_list2 = []
        for col in range(col_count):
            table_item = table.item(row, col)
            df_list2.append('' if table_item is None else str(table_item.text()))
        df_list.append(df_list2)

    df = pd.DataFrame(df_list, columns=headers)

    return df
