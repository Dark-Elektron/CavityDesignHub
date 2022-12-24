import subprocess
import time
import multiprocessing as mp
from threading import Thread
from graphics.graphics_view import GraphicsView
from graphics.scene import Scene
from modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from ui_files.eigenmode import Ui_Eigenmode
from utils.file_reader import FileReader
from ui_files.geometry_input import Ui_Geometry_Input
from utils.shared_classes import *
from utils.shared_functions import *

slans_geom = SLANSGeometry()
fr = FileReader()

file_color = 'green'
DEBUG = True


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class EigenmodeControl:
    def __init__(self, parent):
        self.animation = None
        self.progress_bar = None
        self.pause_icon = None
        self.resume_icon = None
        self.process_state = None
        self.progress_monitor_thread = None
        self.progress_list = None
        self.shape_space = None
        self.end_routine_thread = None
        self.w_Eigenmode = QWidget()

        self.ui = Ui_Eigenmode()
        self.ui.setupUi(self.w_Eigenmode)

        # geometry input
        self.w_Geometry_Input = QWidget()
        self.geom_in = Ui_Geometry_Input()
        self.geom_in.setupUi(self.w_Geometry_Input)

        self.ui.gl_Enter_Geometry.addWidget(self.w_Geometry_Input)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # get logger
        self.log = self.main_control.log

        # ###########################
        # Create Scene
        self.scene = Scene(self)

        # QGraphicsView
        self.graphicsView = GraphicsView(self, 'Eigenmode')
        self.ui.vL_2D_Graphics_View.addWidget(self.graphicsView)
        # ##########################

        self.initUI()
        self.signals()
        self.exe_control()
        self.filename = None  # place holder, made a booboo copying the end routine

        # instantiate geometry
        self.slans_geom = SLANSGeometry()

        # shape space initialization
        self._shape_space = {}
        self._selected_keys = []
        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False

        # ui effects
        self.ui_effects()

    def initUI(self):

        # init shape entry mode
        self.shape_entry_widgets_control()

        # disable expansion section for now. Feature to come later
        self.geom_in.cb_Expansion.setEnabled(False)

        # inner cell
        self.geom_in.cb_Inner_Cell.setCheckState(2)
        self.geom_in.cb_Inner_Cell.setEnabled(False)

        # expand/collapse sections widgets
        if self.geom_in.cb_Expansion.checkState() == 2:
            self.geom_in.w_Expansion.setMinimumHeight(160)
        else:
            self.geom_in.w_Expansion.setMinimumHeight(0)
            self.geom_in.w_Expansion.setMaximumHeight(0)

        if self.geom_in.cb_Outer_Cell_L.checkState() == 2:
            self.geom_in.w_Outer_Cell_L.setMinimumHeight(160)
        else:
            self.geom_in.w_Outer_Cell_L.setMinimumHeight(0)
            self.geom_in.w_Outer_Cell_L.setMaximumHeight(0)

        if self.geom_in.cb_Outer_Cell_R.checkState() == 2:
            self.geom_in.w_Outer_Cell_R.setMinimumHeight(160)
        else:
            self.geom_in.w_Outer_Cell_R.setMinimumHeight(0)
            self.geom_in.w_Outer_Cell_R.setMaximumHeight(0)

        # create pause and resume icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # process state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        # set default boundary condition to magnetic wall at both ends
        self.ui.cb_LBC.setCurrentIndex(2)
        self.ui.cb_RBC.setCurrentIndex(2)

        # create progress bar object and add to widget
        self.progress_bar = QProgressBar(self.ui.w_Simulation_Controls)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.ui.gl_Simulation_Controls.addWidget(self.progress_bar, 0, 4, 1, 1)
        self.progress_bar.hide()

    def signals(self):
        # run eigenmode solver
        self.ui.pb_Run.clicked.connect(lambda: self.run_SLANS())

        # load shape space
        self.geom_in.pb_Select_Shape_Space.clicked.connect(
            lambda: self.open_file(self.geom_in.le_Shape_Space, self.geom_in.cb_Shape_Space_Keys))

        # control shape entry mode
        self.geom_in.cb_Shape_Entry_Mode.currentIndexChanged.connect(lambda: self.shape_entry_widgets_control())

        # cell parameters control signals
        self.geom_in.cb_Outer_Cell_L.stateChanged.connect(lambda: self.animate_height(
            self.geom_in.cb_Outer_Cell_L, self.geom_in.w_Outer_Cell_L, 0, 160, True))
        self.geom_in.cb_Outer_Cell_R.stateChanged.connect(lambda: self.animate_height(
            self.geom_in.cb_Outer_Cell_R, self.geom_in.w_Outer_Cell_R, 0, 160, True))
        self.geom_in.cb_Expansion.stateChanged.connect(lambda: self.animate_height(
            self.geom_in.cb_Expansion, self.geom_in.w_Expansion, 0, 160, True))

        # cancel
        self.ui.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.ui.pb_Pause_Resume.clicked.connect(lambda: self.pause() if self.process_state == 'running' else self.resume())

        # uncomment to draw again
        # self.ui.cb_Shape_Space_Keys.currentTextChanged.connect(lambda: self.draw_shape_from_shape_space())

        #
        self.geom_in.le_Alpha.editingFinished.connect(lambda: self.update_alpha())

        #
        self.geom_in.le_Req_i.editingFinished.connect(lambda: self.geom_in.le_Req_ol.setText(self.geom_in.le_Req_i.text()))
        self.geom_in.le_Req_i.editingFinished.connect(lambda: self.geom_in.le_Req_or.setText(self.geom_in.le_Req_i.text()))

    def shape_entry_widgets_control(self):
        if self.geom_in.cb_Shape_Entry_Mode.currentIndex() == 0:
            self.main_control.animate_height(self.geom_in.w_Select_Shape_Space, 0, 50, True)

            self.ui.w_Enter_Geometry_Manual.setMinimumHeight(0)
            self.ui.w_Enter_Geometry_Manual.setMaximumHeight(0)

            # clear cells from graphics view
            self.graphicsView.removeCells()
        else:
            self.main_control.animate_height(self.ui.w_Enter_Geometry_Manual, 0, 375, True)

            self.geom_in.w_Select_Shape_Space.setMinimumHeight(0)
            self.geom_in.w_Select_Shape_Space.setMaximumHeight(0)

            # uncomment following lines to draw
            # # clear cells from graphics view
            # self.graphicsView.removeCells()
            #
            # # draw new cell
            # self.graphicsView.drawCells(color=QColor(0, 0, 0, 255))

    def run_SLANS(self):
        # get analysis parameters
        n_cells = self.geom_in.sb_N_Cells.value()
        n_modules = self.geom_in.sb_N_Modules.value()
        f_shift = float(self.ui.le_Freq_Shift.text())
        n_modes = float(self.ui.le_No_Of_Modes.text())

        # boundary conditions
        lbc = self.ui.cb_LBC.currentIndex()+1
        rbc = self.ui.cb_RBC.currentIndex()+1
        bc = 10*lbc + rbc

        proc_count = self.ui.sb_No_Of_Processors_SLANS.value()

        # uq
        if self.ui.cb_UQ.isChecked():
            UQ = True
        else:
            UQ = False

        # get geometric parameters
        self.shape_space = self.get_geometric_parameters('SLANS')

        # split shape_space for different processes/ MPI share process by rank
        keys = list(self.shape_space.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

        # show progress bar
        self.progress_bar.show()

        # progress list
        manager = mp.Manager()
        self.progress_list = manager.list()
        self.progress_list.append(0)

        for p in range(proc_count):
            try:
                if p < proc_count - 1:
                    proc_keys_list = keys[p * share:p * share + share]
                else:
                    proc_keys_list = keys[p * share:]

                processor_shape_space = {}
                for key, val in self.shape_space.items():
                    if key in proc_keys_list:
                        processor_shape_space[key] = val
                # print(f'Processor {p}: {processor_shape_space}')
                print_("it's here")
                print_(n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, self.main_control.parentDir,
                       self.main_control.projectDir, self.progress_list, self.ui.le_Run_Save_Folder.text())

                service = mp.Process(target=self.run_sequential, args=(
                    n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, self.main_control.parentDir,
                    self.main_control.projectDir, self.progress_list, self.ui.le_Run_Save_Folder.text(), UQ))

                print_("it's now hersssse")
                service.start()
                print_("it's now here")

                self.processes.append(psutil.Process(service.pid))
                print_(self.processes_id)
                self.processes_id.append(service.pid)

            except Exception as e:
                self.log.error(fr"Exception in run_MP:: {e}")
                # print_("Exception in run_MP::", e)

        # display progress bar
        self.show_progress_bar = True
        self.progress_monitor_thread = ProgressMonitor(self, self.main_control.projectDir)
        self.progress_monitor_thread.sig.connect(self.update_progress_bar)
        self.progress_monitor_thread.start()

        self.log.info("Eigenmode simulation started")
        # change process state to running
        self.process_state = 'running'
        self.run_pause_resume_stop_routine()

        self.end_routine_thread = EndRoutine(self, self.main_control.projectDir)
        self.end_routine_thread.start()

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
        except:
            pass

        self.processes.clear()
        self.processes_id.clear()

        self.process_state = 'none'
        self.run_pause_resume_stop_routine()
        self.log.info("Process terminated.")

    def end_routine(self, proc_ids):
        print(proc_ids, type(proc_ids))
        for pid in proc_ids:
            try:
                p = psutil.Process(pid)
                while p.is_running():
                    pass

                print(fr"process {p} ended")
            except:
                pass

        self.cancel()

    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

        if val == 100 or not self.show_progress_bar:
            # reset progress bar
            self.progress_bar.setValue(0)
            self.progress_bar.hide()

    def get_geometric_parameters(self, code):
        self.shape_space = {}
        print_('Getting geometric parameters')
        if self.geom_in.cb_Shape_Entry_Mode.currentIndex() == 0:
            print_("Test worked")
            try:
                # self._shape_space = self.load_shape_space(shape_space_name)
                # print_(self._shape_space)

                # get selected keys
                self._selected_keys = self.geom_in.cb_Shape_Space_Keys.currentText()
                print("selected keys: ", self.geom_in.cb_Shape_Space_Keys.currentText())
                # print("Selected keys: ", self._selected_keys, type(self._selected_keys[0]))

                # check keys of shape space if results already exist
                toall = None
                for key, val in self._shape_space.items():
                    # process for only keys selected in combobox
                    if self.geom_in.cb_Shape_Space_Keys.currentText() == "" \
                            or self.geom_in.cb_Shape_Space_Keys.currentText() == "All":
                        pass
                    else:
                        if key not in self._selected_keys:
                            continue

                    if not toall:
                        ans = self.prompt(code, key)
                        if ans == 'Yes':
                            self.shape_space[key] = val

                        if ans == 'No':
                            continue

                        if ans == 'YesToAll':
                            self.shape_space[key] = val
                            toall = 'YesToAll'

                        if ans == 'NoToAll':
                            toall = 'NoToAll'

                        if ans == "Does not exist":
                            self.shape_space[key] = val
                            toall = None
                    else:
                        if toall == 'YesToAll':
                            self.shape_space[key] = val
                        else:
                            path = f'{self.main_control.projectDir}/SimulationData/{code}/Cavity{key}'
                            if os.path.exists(path):
                                continue
                            else:
                                self.shape_space[key] = val

                # print_(self.shape_space)
                return self.shape_space
            except Exception as e:
                print_(f"File not found, check path:: {e}")
        else:
            if self.geom_in.cb_Inner_Cell.checkState() == 2:
                # Middle Ellipse data
                A_i_space = self.text_to_list(self.geom_in.le_A_i.text())
                B_i_space = self.text_to_list(self.geom_in.le_B_i.text())
                a_i_space = self.text_to_list(self.geom_in.le_a_i.text())
                b_i_space = self.text_to_list(self.geom_in.le_b_i.text())
                Ri_i_space = self.text_to_list(self.geom_in.le_Ri_i.text())
                L_i_space = self.text_to_list(self.geom_in.le_L_i.text())
                Req_i_space = self.text_to_list(self.geom_in.le_Req_i.text())
                alpha_i_space = self.text_to_list(self.geom_in.le_Alpha.text())

                inner_cell_space = [A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space, alpha_i_space]
            else:
                inner_cell_space = [[0], [0], [0], [0], [0], [0], [0], [0]]

            if self.geom_in.cb_Outer_Cell_L.checkState() == 2:
                # Middle Ellipse data
                A_ol_space = self.text_to_list(self.geom_in.le_A_ol.text())
                B_ol_space = self.text_to_list(self.geom_in.le_B_ol.text())
                a_ol_space = self.text_to_list(self.geom_in.le_a_ol.text())
                b_ol_space = self.text_to_list(self.geom_in.le_b_ol.text())
                Ri_ol_space = self.text_to_list(self.geom_in.le_Ri_ol.text())
                L_ol_space = self.text_to_list(self.geom_in.le_L_ol.text())
                Req_ol_space = self.text_to_list(self.geom_in.le_Req_ol.text())
                alpha_ol_space = self.text_to_list(self.geom_in.le_Alpha_ol.text())

                outer_cell_L_space = [A_ol_space, B_ol_space, a_ol_space, b_ol_space, Ri_ol_space, L_ol_space, Req_ol_space, alpha_ol_space]
            else:
                outer_cell_L_space = inner_cell_space

            if self.geom_in.cb_Outer_Cell_R.checkState() == 2:
                # Middle Ellipse data
                A_or_space = self.text_to_list(self.geom_in.le_A_or.text())
                B_or_space = self.text_to_list(self.geom_in.le_B_or.text())
                a_or_space = self.text_to_list(self.geom_in.le_a_or.text())
                b_or_space = self.text_to_list(self.geom_in.le_b_or.text())
                Ri_or_space = self.text_to_list(self.geom_in.le_Ri_or.text())
                L_or_space = self.text_to_list(self.geom_in.le_L_or.text())
                Req_or_space = self.text_to_list(self.geom_in.le_Req_or.text())
                alpha_or_space = self.text_to_list(self.geom_in.le_Alpha_or.text())

                outer_cell_R_space = [A_or_space, B_or_space, a_or_space, b_or_space, Ri_or_space, L_or_space, Req_or_space, alpha_or_space]
            else:
                outer_cell_R_space = inner_cell_space
            count = 0
            for A_i in inner_cell_space[0]:
                for B_i in inner_cell_space[1]:
                    for a_i in inner_cell_space[2]:
                        for b_i in inner_cell_space[3]:
                            for Ri_i in inner_cell_space[4]:
                                for L_i in inner_cell_space[5]:
                                    for Req_i in inner_cell_space[6]:
                                        if outer_cell_L_space == inner_cell_space:
                                            inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                            outer_cell_L = inner_cell

                                            if self.geom_in.cb_LBP.checkState() == 2 and self.geom_in.cb_RBP.checkState() == 2:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'both', 'FREQ': None}
                                            elif self.geom_in.cb_LBP.checkState() == 2 and self.geom_in.cb_RBP.checkState() == 0:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'left', 'FREQ': None}
                                            elif self.geom_in.cb_LBP.checkState() == 0 and self.geom_in.cb_RBP.checkState() == 2:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'right', 'FREQ': None}
                                            else:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'none', 'FREQ': None}

                                            count += 1
                                        else:
                                            for A_ol in outer_cell_L_space[0]:
                                                for B_ol in outer_cell_L_space[1]:
                                                    for a_ol in outer_cell_L_space[2]:
                                                        for b_ol in outer_cell_L_space[3]:
                                                            for Ri_ol in outer_cell_L_space[4]:
                                                                for L_ol in outer_cell_L_space[5]:
                                                                    # for Req_ol in outer_cell_L_space[6]:
                                                                    if outer_cell_L_space == outer_cell_R_space:
                                                                        inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                                                        outer_cell_L = [A_ol, B_ol, a_ol, b_ol, Ri_ol, L_ol, Req_i, 0]
                                                                        outer_cell_R = outer_cell_L
                                                                        if self.geom_in.cb_LBP.checkState() == 2 and self.geom_in.cb_RBP.checkState() == 0:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                        elif self.geom_in.cb_LBP.checkState() == 0 and self.geom_in.cb_RBP.checkState() == 2:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                        elif self.geom_in.cb_LBP.checkState() == 2 and self.geom_in.cb_RBP.checkState() == 2:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'both', 'FREQ': None}
                                                                        else:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'none', 'FREQ': None}

                                                                        count += 1
                                                                    else:
                                                                        for A_or in outer_cell_R_space[0]:
                                                                            for B_or in outer_cell_R_space[1]:
                                                                                for a_or in outer_cell_R_space[2]:
                                                                                    for b_or in outer_cell_R_space[3]:
                                                                                        for Ri_or in outer_cell_R_space[4]:
                                                                                            for L_or in outer_cell_R_space[5]:
                                                                                                # for Req_or in outer_cell_R_space[6]:
                                                                                                inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                                                                                outer_cell_L = [A_ol, B_ol, a_ol, b_ol, Ri_ol, L_ol, Req_i, 0]
                                                                                                outer_cell_R = [A_or, B_or, a_or, b_or, Ri_or, L_or, Req_i, 0]
                                                                                                if self.geom_in.cb_LBP.checkState() == 2 and self.geom_in.cb_RBP.checkState() == 0:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                                                elif self.geom_in.cb_LBP.checkState() == 0 and self.geom_in.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                                                elif self.geom_in.cb_LBP.checkState() == 2 and self.geom_in.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'both', 'FREQ': None}
                                                                                                else:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'none', 'FREQ': None}

                                                                                                count += 1
            return self.shape_space

    def prompt(self, code, fid):
        path = fr'{self.main_control.projectDir}\SimulationData\SLANS\Cavity{fid}'
        # print(path)
        # path = os.path.join(path, fr"{}\{code}\Cavity{fid}")
        if os.path.exists(path):
            print_("Simulation data already exists. Do you want to overwrite it?")
            msg = QMessageBox()
            msg.setWindowTitle("Folder Exist")
            msg.setText("Simulation data already exists. Do you want to overwrite it?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.YesToAll | QMessageBox.Yes | QMessageBox.No | QMessageBox.NoToAll)
            msg.setDefaultButton(QMessageBox.Yes)

            msg.buttonClicked.connect(self.button_clicked)
            # print_(f'msg: {msg.Yes}')

            x = msg.exec_()

            print_("got here")
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

    def open_file(self, le, cb):
        # clear combobox
        self.geom_in.cb_Shape_Space_Keys.clear()
        self.geom_in.cb_Shape_Space_Keys.addItem('All')
        # self._selected_keys.clear()

        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Json Files (*.json)")
        try:
            le.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            # populate checkboxes with key
            for col in dd.keys():
                cb.addItem(fr'{col}')

            self._shape_space = dd

        except Exception as e:
            print('Failed to open file:: ', e)

    def animate_width(self, cb, widget, min_width, standard, enable, reverse=False):
        if enable:
            # GET Height
            width = widget.width()
            # SET MAX WIDTH

            if not reverse:
                if cb.checkState() != 2:
                    widthCollapsed = min_width
                    widget.setMinimumWidth(0)
                else:
                    widthCollapsed = standard
                    # widget.setMinimumWidth(standard)
            else:
                if cb.checkState() == 2:
                    widthCollapsed = min_width
                    widget.setMinimumWidth(0)
                else:
                    widthCollapsed = standard
                    # widget.setMinimumWidth(standard)

            # ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumWidth")
            self.animation.setDuration(200)
            self.animation.setStartValue(width)
            self.animation.setEndValue(widthCollapsed)
            self.animation.start()

    def animate_height(self, cb, widget, min_height, standard, enable):
        if enable:
            # GET WIDTH
            height = widget.width()
            # SET MAX WIDTH
            if cb.checkState() != 2:
                heightCollapsed = min_height
                widget.setMinimumHeight(0)
            else:
                heightCollapsed = standard
                # widget.setMinimumWidth(standard)

            # ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumHeight")
            self.animation.setDuration(200)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.start()

    def exe_control(self):
        # Slans
        self.ui.pb_Genmsh.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\genmesh2.exe"))
        self.ui.pb_Sl.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\Sl.exe"))
        self.ui.pb_Slansc.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\slansc.exe"))
        self.ui.pb_Slansm.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\slansm.exe"))
        self.ui.pb_Slanss.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\slanss.exe"))
        self.ui.pb_Slansre.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\slansre.exe"))
        self.ui.pb_MTFView.clicked.connect(lambda: self.run_slans_exe(r"SLANS_exe\Mtfview\mtfview.exe"))

    def run_slans_exe(self, path, filename=None):
        path = fr"{self.main_control.parentDir}\em_codes\{path}"
        t = Thread(target=subprocess.call, args=(path,))
        t.start()

    def draw_shape_from_shape_space(self):
        colors = [[48, 162, 218, 255], [252, 79, 48, 255], [229, 174, 56, 255],
                  [109, 144, 79, 255], [139, 139, 139, 255]]
        ci = 0

        # remove existing cells
        self.graphicsView.removeCells()
        for key in self._shape_space.keys():
            if key in self.geom_in.cb_Shape_Space_Keys.currentText():
                IC = self._shape_space[key]["IC"]
                OC = self._shape_space[key]["OC"]
                BP = self._shape_space[key]["BP"]
                self.graphicsView.drawCells(IC, OC, BP, QColor(colors[ci][0], colors[ci][1], colors[ci][2], colors[ci][3]))

                ci += 1

    def ui_effects(self):

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.ui.w_Settings.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.geom_in.w_Inner_Cell.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.geom_in.w_Outer_Cell_L.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.geom_in.w_Outer_Cell_R.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.geom_in.w_Expansion.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.ui.w_Simulation_Controls.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.ui.w_Show_Cavity.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))
        self.geom_in.w_Load_Manual.setGraphicsEffect(shadow)

    def serialize(self, state_dict):
        # update state file
        state_dict["Eigen_Shape_Entry_Mode"] = self.geom_in.cb_Shape_Entry_Mode.currentIndex()
        state_dict["Eigen_Shape Space"] = self.geom_in.le_Shape_Space.text()
        state_dict["Eigen_Mid_Cell_CB"] = self.geom_in.cb_Inner_Cell.checkState()
        state_dict["Eigen_Left_Cell_CB"] = self.geom_in.cb_Outer_Cell_L.checkState()
        state_dict["Eigen_Right_Cell_CB"] = self.geom_in.cb_Outer_Cell_R.checkState()
        state_dict["Eigen_Expansion_CB"] = self.geom_in.cb_Expansion.checkState()
        state_dict["Eigen_LBP_CB"] = self.geom_in.cb_LBP.checkState()
        state_dict["Eigen_RBP_CB"] = self.geom_in.cb_RBP.checkState()

        # cell parameters
        state_dict["Eigen_A_i"] = self.geom_in.le_A_i.text()
        state_dict["Eigen_B_i"] = self.geom_in.le_B_i.text()
        state_dict["Eigen_a_i"] = self.geom_in.le_a_i.text()
        state_dict["Eigen_b_i"] = self.geom_in.le_b_i.text()
        state_dict["Eigen_Ri_i"] = self.geom_in.le_Ri_i.text()
        state_dict["Eigen_L_i"] = self.geom_in.le_L_i.text()
        state_dict["Eigen_Req_i"] = self.geom_in.le_Req_i.text()
        state_dict["Eigen_Alpha_i"] = self.geom_in.le_Alpha.text()

        state_dict["Eigen_A_ol"] = self.geom_in.le_A_ol.text()
        state_dict["Eigen_B_ol"] = self.geom_in.le_B_ol.text()
        state_dict["Eigen_a_ol"] = self.geom_in.le_a_ol.text()
        state_dict["Eigen_b_ol"] = self.geom_in.le_b_ol.text()
        state_dict["Eigen_Ri_ol"] = self.geom_in.le_Ri_ol.text()
        state_dict["Eigen_L_ol"] = self.geom_in.le_L_ol.text()
        state_dict["Eigen_Req_ol"] = self.geom_in.le_Req_ol.text()
        state_dict["Eigen_Alpha_ol"] = self.geom_in.le_Alpha_ol.text()

        state_dict["Eigen_A_or"] = self.geom_in.le_A_or.text()
        state_dict["Eigen_B_or"] = self.geom_in.le_B_or.text()
        state_dict["Eigen_a_or"] = self.geom_in.le_a_or.text()
        state_dict["Eigen_b_or"] = self.geom_in.le_b_or.text()
        state_dict["Eigen_Ri_or"] = self.geom_in.le_Ri_or.text()
        state_dict["Eigen_L_or"] = self.geom_in.le_L_or.text()
        state_dict["Eigen_Req_or"] = self.geom_in.le_Req_or.text()
        state_dict["Eigen_Alpha_or"] = self.geom_in.le_Alpha_or.text()

        # settings
        state_dict["Eigen_N_Cells"] = self.geom_in.sb_N_Cells.value()
        state_dict["Eigen_N_Modules"] = self.geom_in.sb_N_Modules.value()
        state_dict["Eigen_Polarization"] = self.ui.cb_Polarization_SLANS.currentIndex()

        state_dict["Eigen_Freq_Shift"] = self.ui.le_Freq_Shift.text()
        state_dict["Eigen_No_Of_Modes"] = self.ui.le_No_Of_Modes.text()
        state_dict["Eigen_LBC"] = self.ui.cb_LBC.currentIndex()
        state_dict["Eigen_RBC"] = self.ui.cb_RBC.currentIndex()
        state_dict["Eigen_No_Of_Processors"] = self.ui.sb_No_Of_Processors_SLANS.value()

    def deserialize(self, state_dict):
        # update state file
        self.geom_in.cb_Shape_Entry_Mode.setCurrentIndex(state_dict["Eigen_Shape_Entry_Mode"])
        self.geom_in.le_Shape_Space.setText(state_dict["Eigen_Shape Space"])
        self.geom_in.cb_Inner_Cell.setCheckState(state_dict["Eigen_Mid_Cell_CB"])
        self.geom_in.cb_Outer_Cell_L.setCheckState(state_dict["Eigen_Left_Cell_CB"])
        self.geom_in.cb_Outer_Cell_R.setCheckState(state_dict["Eigen_Right_Cell_CB"])
        self.geom_in.cb_Expansion.setCheckState(state_dict["Eigen_Expansion_CB"])
        self.geom_in.cb_LBP.setCheckState(state_dict["Eigen_LBP_CB"])
        self.geom_in.cb_RBP.setCheckState(state_dict["Eigen_RBP_CB"])

        # cell parameters
        self.geom_in.le_A_i.setText(state_dict["Eigen_A_i"])
        self.geom_in.le_B_i.setText(state_dict["Eigen_B_i"])
        self.geom_in.le_a_i.setText(state_dict["Eigen_a_i"])
        self.geom_in.le_b_i.setText(state_dict["Eigen_b_i"])
        self.geom_in.le_Ri_i.setText(state_dict["Eigen_Ri_i"])
        self.geom_in.le_L_i.setText(state_dict["Eigen_L_i"])
        self.geom_in.le_Req_i.setText(state_dict["Eigen_Req_i"])
        self.geom_in.le_Alpha.setText(state_dict["Eigen_Alpha_i"])

        self.geom_in.le_A_ol.setText(state_dict["Eigen_A_ol"])
        self.geom_in.le_B_ol.setText(state_dict["Eigen_B_ol"])
        self.geom_in.le_a_ol.setText(state_dict["Eigen_a_ol"])
        self.geom_in.le_b_ol.setText(state_dict["Eigen_b_ol"])
        self.geom_in.le_Ri_ol.setText(state_dict["Eigen_Ri_ol"])
        self.geom_in.le_L_ol.setText(state_dict["Eigen_L_ol"])
        self.geom_in.le_Req_ol.setText(state_dict["Eigen_Req_ol"])
        self.geom_in.le_Alpha_ol.setText(state_dict["Eigen_Alpha_ol"])

        self.geom_in.le_A_or.setText(state_dict["Eigen_A_or"])
        self.geom_in.le_B_or.setText(state_dict["Eigen_B_or"])
        self.geom_in.le_a_or.setText(state_dict["Eigen_a_or"])
        self.geom_in.le_b_or.setText(state_dict["Eigen_b_or"])
        self.geom_in.le_Ri_or.setText(state_dict["Eigen_Ri_or"])
        self.geom_in.le_L_or.setText(state_dict["Eigen_L_or"])
        self.geom_in.le_Req_or.setText(state_dict["Eigen_Req_or"])
        self.geom_in.le_Alpha_or.setText(state_dict["Eigen_Alpha_or"])

        # settings
        self.geom_in.sb_N_Cells.setValue(state_dict["Eigen_N_Cells"])
        self.geom_in.sb_N_Modules.setValue(state_dict["Eigen_N_Modules"])
        self.ui.cb_Polarization_SLANS.setCurrentIndex(state_dict["Eigen_Polarization"])
        self.ui.le_Freq_Shift.setText(state_dict["Eigen_Freq_Shift"])
        self.ui.le_No_Of_Modes.setText(state_dict["Eigen_No_Of_Modes"])
        self.ui.cb_LBC.setCurrentIndex(state_dict["Eigen_LBC"])
        self.ui.cb_RBC.setCurrentIndex(state_dict["Eigen_RBC"])
        self.ui.sb_No_Of_Processors_SLANS.setValue(state_dict["Eigen_No_Of_Processors"])

    @staticmethod
    def load_shape_space(filename):
        dirc = filename

        # check if extension is included
        if dirc.split('.')[-1] != 'json':
            dirc = f'{dir}.json'

        df = fr.json_reader(dirc)
        # print_(df)

        return df.to_dict()

    @staticmethod
    def show_hide_(wid1, wid2):
        print('here')
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()

    @staticmethod
    def button_clicked(i):
        return i.text()

    @staticmethod
    def text_to_list(txt):
        if "range" in txt:
            txt = txt.replace('range', '')
            l = ast.literal_eval(txt)
            return range(l[0], l[1], l[2])
        elif 'linspace' in txt:
            l = eval(f'np.{txt}')
            print(l)
            return l
        else:
            l = ast.literal_eval(txt)
            if isinstance(l, int) or isinstance(l, float):
                return [l]
            else:
                return list(l)

    @staticmethod
    def run_sequential(n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, parentDir, projectDir,
                       progress_list, sub_dir='', UQ=False):
        progress = 0
        # get length of processor
        total_no_of_shapes = len(list(processor_shape_space.keys()))

        for key, shape in processor_shape_space.items():
            # # create folders for all keys
            slans_geom.createFolder(key, projectDir, subdir=sub_dir)

            write_cst_paramters(key, shape['IC'], shape['OC'], projectDir=projectDir, cell_type="None")

            # run SLANS code
            start_time = time.time()
            try:
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
            except KeyError:
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

            # run UQ
            if UQ:
                print_("UQ")
                uq(key, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
                   n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
                   f_shift=f_shift, bc=bc, parentDir=parentDir, projectDir=projectDir)

            print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')

            # update progress
            progress_list.append((progress+1)/total_no_of_shapes)

    def update_alpha(self):
        A_i_space = self.text_to_list(self.geom_in.le_A_i.text())[0]
        B_i_space = self.text_to_list(self.geom_in.le_B_i.text())[0]
        a_i_space = self.text_to_list(self.geom_in.le_a_i.text())[0]
        b_i_space = self.text_to_list(self.geom_in.le_b_i.text())[0]
        Ri_i_space = self.text_to_list(self.geom_in.le_Ri_i.text())[0]
        L_i_space = self.text_to_list(self.geom_in.le_L_i.text())[0]
        Req_i_space = self.text_to_list(self.geom_in.le_Req_i.text())[0]
        print(A_i_space)
        try:
            alpha_i_space = calculate_alpha(A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space, 0)
            self.geom_in.le_Alpha.setText(f"{round(alpha_i_space, 2)}")
        except:
            pass


def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, parentDir, projectDir):
    err = False
    result_dict_slans = {}
    slans_obj_list = qois
    for o in qois:
        result_dict_slans[o] = {'expe': [], 'stdDev': []}

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = shape['IC'][0:5]
    # ic(p_true)
    rdim = len(p_true)  # How many variabels will be considered as random in our case 5
    degree = 1

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
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
    ic(nodes)
    delta = 0.005  # or 0.1

    Ttab_val_f = []
    print_('3')
    sub_dir = fr'Cavity{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    for i in range(no_sims):
        skip = False
        p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
        p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
        p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
        p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
        p_init[4] = p_true[4] * (1 + delta * nodes[4, i])

        par_mid = np.append(p_init, shape['IC'][5:]).tolist()
        par_end = par_mid
        print_(par_mid)

        # perform checks on geometry
        ok = perform_geometry_checks(par_mid, par_end)
        print_("OK", ok)
        if not ok:
            err = True
            break
        fid = fr'{key}_Q{i}'

        # check if folder exists and skip if it does
        print_(fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}')
        if os.path.exists(fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}'):
            skip = True
            # ic("Skipped: ", fid, fr'{projectDir}\SimulationData\ABCI\Cavity{key}\Cavity{fid}')

        # skip analysis if folder already exists.
        if not skip:
            #  run model using SLANS or CST
            # # create folders for all keys
            slans_geom.createFolder(fid, projectDir, subdir=sub_dir)
            try:
                slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
            except KeyError:
                slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

        filename = fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}\cavity_{bc}.svl'
        if os.path.exists(filename):
            params = fr.svl_reader(filename)
            norm_length = 2*n_cells*shape['IC'][5]
            ic(n_cells, norm_length)
            qois_result = get_qoi_value(params, slans_obj_list, n_cells, norm_length)
            print_(qois_result)
            # sometimes some degenerate shapes are still generated and the solver returns zero
            # for the objective functions, such shapes are considered invalid
            for objr in qois_result:
                if objr == 0:
                    # skip key
                    err = True
                    break

            tab_val_f = qois_result

            Ttab_val_f.append(tab_val_f)
        else:
            err = True

    # # add original point
    # filename = fr'{projectDir}\SimulationData\SLANS\Cavity{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
        # Ttab_val_f.append(tab_val_f)

    print_("Error: ", err)
    ic(Ttab_val_f)
    # import matplotlib.pyplot as plt
    if not err:
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
        # ic(v_expe_fobj, v_stdDev_fobj)
        # append results to dict
        ic(v_expe_fobj, v_stdDev_fobj)
        for i, o in enumerate(slans_obj_list):
            result_dict_slans[o]['expe'].append(v_expe_fobj[i])
            result_dict_slans[o]['stdDev'].append(v_stdDev_fobj[i])

            # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
            # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

        # plt.show()

        with open(fr"{projectDir}\SimulationData\SLANS\Cavity{key}\uq.json", 'w') as file:
            file.write(json.dumps(result_dict_slans, indent=4, separators=(',', ': ')))
    else:
        print_(fr"There was a problem running UQ analysis for Cavity{key}")


def get_qoi_value(d, obj, n_cells, norm_length):
    Req = d['CAVITY RADIUS'][n_cells-1] * 10  # convert to mm
    Freq = d['FREQUENCY'][n_cells-1]
    E_stored = d['STORED ENERGY'][n_cells-1]
    Rsh = d['SHUNT IMPEDANCE'][n_cells-1]  # MOhm
    Q = d['QUALITY FACTOR'][n_cells-1]
    Epk = d['MAXIMUM ELEC. FIELD'][n_cells-1]  # MV/m
    Hpk = d['MAXIMUM MAG. FIELD'][n_cells-1]  # A/m
    # Vacc = dict['ACCELERATION'][0]
    Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells-1]  # MV/m
    Rsh_Q = d['EFFECTIVE IMPEDANCE'][n_cells-1]  # Ohm

    Vacc = np.sqrt(
        2 * Rsh_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6
    # factor of 2, remember circuit and accelerator definition
    # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
    Eacc = Vacc / (norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
    Epk_Eacc = Epk / Eacc
    Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc

    d = {
        "Req": Req,
        "freq": Freq,
        "Q": Q,
        "E": E_stored,
        "R/Q": 2 * Rsh_Q,
        "Epk/Eacc": Epk_Eacc,
        "Bpk/Eacc": Bpk_Eacc
    }

    objective = []

    # append objective functions
    for o in obj:
        if o in d.keys():
            objective.append(d[o])

    return objective
