import ast
import json
import os
import subprocess
import threading
import time
import multiprocessing as mp
from threading import Thread
from PyQt5 import QtCore
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QPropertyAnimation
from termcolor import colored
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from graphics.graphics_view import GraphicsView
from graphics.scene import Scene
from simulation_codes.ABCI.abci_geometry import ABCIGeometry
from ui_files.wakefield import Ui_Wakefield
from utils.file_reader import FileReader
import psutil

fr = FileReader()
abci_geom = ABCIGeometry()

file_color = 'red'


def print_(*arg):
    print(colored(f'\t{arg}', file_color))


class WakefieldControl:
    def __init__(self, parent):
        self.w_Wakefield = QWidget()

        self.wakefieldUI = Ui_Wakefield()
        self.wakefieldUI.setupUi(self.w_Wakefield)

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
        self.wakefieldUI.vL_2D_Graphics_View.addWidget(self.graphicsView)

        # ##########################

        self.initUI()
        self.signals()
        self.exe_control()
        # self.parentDir = self.main_control.parentDir
        # self.projectDir = self.main_control.projectDir

        # instantiate geometry
        self.abci_geom = ABCIGeometry()

        # shape space initialization
        self._shape_space = {}
        self._selected_keys = []
        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False

        # ui effects
        self.ui_effects()

    def signals(self):
        # signals
        self.wakefieldUI.pb_Run.clicked.connect(lambda: self.run_ABCI())

        # load shape space
        self.wakefieldUI.pb_Select_Shape_Space.clicked.connect(lambda: self.open_file(self.wakefieldUI.le_Shape_Space, self.cb_Shape_Space_Keys))

        # control shape entry mode
        self.wakefieldUI.cb_Shape_Entry_Mode.currentIndexChanged.connect(lambda: self.shape_entry_widgets_control())

        # cell parameters control signals
        self.wakefieldUI.cb_Outer_Cell_L.stateChanged.connect(lambda: self.animate_width(self.wakefieldUI.cb_Outer_Cell_L, self.wakefieldUI.w_Outer_Cell_L, 0, 300, True))
        self.wakefieldUI.cb_Outer_Cell_R.stateChanged.connect(lambda: self.animate_width(self.wakefieldUI.cb_Outer_Cell_R, self.wakefieldUI.w_Outer_Cell_R, 0, 300, True))
        self.wakefieldUI.cb_Expansion.stateChanged.connect(lambda: self.animate_width(self.wakefieldUI.cb_Expansion, self.wakefieldUI.w_Expansion, 0, 300, True))

        # cancel
        self.wakefieldUI.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.wakefieldUI.pb_Pause_Resume.clicked.connect(lambda: self.pause() if self.process_state == 'running' else self.resume())

        #
        self.cb_Shape_Space_Keys.currentTextChanged.connect(lambda: self.draw_shape_from_shape_space())

    def initUI(self):
        # add checkable combobox
        self.cb_Shape_Space_Keys = CheckableComboBox()
        # self.cb_Shape_Space_Keys.addItem('All')
        self.cb_Shape_Space_Keys.setMinimumWidth(75)
        self.wakefieldUI.gl_Shape_Space_Keys.addWidget(self.cb_Shape_Space_Keys)

        # init shape entry mode
        self.shape_entry_widgets_control()

        # disable expansion section for now. Feature to come later
        self.wakefieldUI.cb_Expansion.setEnabled(False)

        # inner cell
        self.wakefieldUI.cb_Inner_Cell.setCheckState(2)
        self.wakefieldUI.cb_Inner_Cell.setEnabled(False)

        # expand/collapse sections widgets
        if self.wakefieldUI.cb_Expansion.checkState() == 2:
            self.wakefieldUI.w_Expansion.setMinimumWidth(300)
        else:
            self.wakefieldUI.w_Expansion.setMinimumWidth(0)
            self.wakefieldUI.w_Expansion.setMaximumWidth(0)

        if self.wakefieldUI.cb_Outer_Cell_L.checkState() == 2:
            self.wakefieldUI.w_Outer_Cell_L.setMinimumWidth(300)
        else:
            self.wakefieldUI.w_Outer_Cell_L.setMinimumWidth(0)
            self.wakefieldUI.w_Outer_Cell_L.setMaximumWidth(0)

        if self.wakefieldUI.cb_Outer_Cell_R.checkState() == 2:
            self.wakefieldUI.w_Outer_Cell_R.setMinimumWidth(300)
        else:
            self.wakefieldUI.w_Outer_Cell_R.setMinimumWidth(0)
            self.wakefieldUI.w_Outer_Cell_R.setMaximumWidth(0)

        # wakefield analysis always uses beam pipes
        self.wakefieldUI.cb_LBP.setCheckState(2)
        self.wakefieldUI.cb_RBP.setCheckState(2)
        self.wakefieldUI.cb_LBP.setEnabled(False)
        self.wakefieldUI.cb_RBP.setEnabled(False)

        # create pause and resume icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # process state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        # create progress bar object and add to widget
        self.progress_bar = QProgressBar(self.wakefieldUI.w_Simulation_Controls)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.wakefieldUI.gl_Simulation_Controls.addWidget(self.progress_bar, 0, 4, 1, 1)
        self.progress_bar.hide()

    def shape_entry_widgets_control(self):
        if self.wakefieldUI.cb_Shape_Entry_Mode.currentIndex() == 0:
            self.main_control.animate_height(self.wakefieldUI.w_Select_Shape_Space, 0, 50, True)

            self.wakefieldUI.w_Enter_Geometry_Manual.setMinimumHeight(0)
            self.wakefieldUI.w_Enter_Geometry_Manual.setMaximumHeight(0)

            # clear cells from graphics view
            self.graphicsView.removeCells()
        else:
            self.main_control.animate_height(self.wakefieldUI.w_Enter_Geometry_Manual, 0, 375, True)

            self.wakefieldUI.w_Select_Shape_Space.setMinimumHeight(0)
            self.wakefieldUI.w_Select_Shape_Space.setMaximumHeight(0)

            # clear cells from graphics view
            self.graphicsView.removeCells()

            # draw new cell
            self.graphicsView.drawCells(color=QColor(0, 0, 0, 255))

    def run_ABCI(self):
        # get analysis parameters
        n_cells = self.wakefieldUI.sb_N_Cells.value()
        n_modules = self.wakefieldUI.sb_N_Modules.value()
        MROT = self.wakefieldUI.cb_Polarization_ABCI.currentIndex()
        MT = float(self.wakefieldUI.le_MT.text())  # number of time steps for a beam to move one cell to another default = 3
        bunch_length = float(self.wakefieldUI.le_Bunch_Length.text())
        NFS = float(self.wakefieldUI.le_NFS.text()) # Number of samples in FFT (max 10000)
        UBT = float(self.wakefieldUI.le_Wakelength.text())
        DDZ_SIG = float(self.wakefieldUI.le_DDZ_SIG.text())
        DDR_SIG = float(self.wakefieldUI.le_DDR_SIG.text())
        proc_count = self.wakefieldUI.sb_No_Of_Processors_ABCI.value()

        # get geometric parameters
        shape_space = self.get_geometric_parameters('ABCI')

        # split shape_space for different processes/ MPI share process by rank
        keys = list(shape_space.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

        # show progress bar
        self.progress_bar.show()

        # progress list
        manager = mp.Manager()
        self.progress_list = manager.list()
        self.progress_list.append(0)

        self.processes = []
        for p in range(proc_count):
            try:
                if p < proc_count - 1:
                    proc_keys_list = keys[p * share:p * share + share]
                else:
                    proc_keys_list = keys[p * share:]

                processor_shape_space = {}
                for key, val in shape_space.items():
                    if key in proc_keys_list:
                        processor_shape_space[key] = val

                service = mp.Process(target=self.run_sequential, args=(n_cells, n_modules, processor_shape_space,
                                                                       MROT, MT, NFS, UBT, bunch_length,
                                                                       DDR_SIG, DDZ_SIG,
                                                                       self.main_control.parentDir,
                                                                       self.main_control.projectDir, self.progress_list
                                                                       ))

                service.start()
                self.processes.append(psutil.Process(service.pid))
                self.processes_id.append(service.pid)

            except Exception as e:
                self.log.error(f"Exception in run_MP:: {e}")

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

    def run_pause_resume_stop_routine(self):
        if self.process_state == 'none':
            # change pause/resume icon to pause icon
            self.wakefieldUI.pb_Pause_Resume.setIcon(self.pause_icon)

            # disable pause/resume and cancel buttons
            self.wakefieldUI.pb_Pause_Resume.setEnabled(False)
            self.wakefieldUI.pb_Cancel.setEnabled(False)

            # enable run button in case it was disabled
            self.wakefieldUI.pb_Run.setEnabled(True)

        if self.process_state == "running":
            # enable run, pause/resume and cancel buttons
            self.wakefieldUI.pb_Pause_Resume.setEnabled(True)
            self.wakefieldUI.pb_Cancel.setEnabled(True)
            self.wakefieldUI.pb_Run.setEnabled(False)

            # change pause/resume icon to pause icon
            self.wakefieldUI.pb_Pause_Resume.setIcon(self.pause_icon)

        if self.process_state == 'paused':
            # disable run button
            self.wakefieldUI.pb_Run.setEnabled(False)

            # change pause/resume button icon to resume icon
            self.wakefieldUI.pb_Pause_Resume.setIcon(self.resume_icon)

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
        if self.wakefieldUI.cb_Shape_Entry_Mode.currentIndex() == 0:
            print_("Test worked")
            try:
                # self._shape_space = self.load_shape_space(shape_space_name)
                print_(self._shape_space)

                # get selected keys
                self._selected_keys = self.cb_Shape_Space_Keys.currentText()
                # print("Selected keys: ", self._selected_keys, type(self._selected_keys[0]))


                # check keys of shape space if results already exist
                toall = None
                for key, val in self._shape_space.items():
                    # process for only keys selected in combobox
                    if self.cb_Shape_Space_Keys.currentText() == "":
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
                    else:
                        if toall == 'YesToAll':
                            self.shape_space[key] = val
                        else:
                            path = f'{self.main_control.projectDir}/SimulationData/{code}/Cavity{key}'
                            if os.path.exists(path):
                                continue
                            else:
                                self.shape_space[key] = val

                print_(self.shape_space)
                return self.shape_space
            except Exception as e:
                print_(f"File not found, check path:: {e}")
        else:
            if self.wakefieldUI.cb_Inner_Cell.checkState() == 2:
                # Middle Ellipse data
                A_i_space = self.text_to_list(self.wakefieldUI.le_A_i.text())
                B_i_space = self.text_to_list(self.wakefieldUI.le_B_i.text())
                a_i_space = self.text_to_list(self.wakefieldUI.le_a_i.text())
                b_i_space = self.text_to_list(self.wakefieldUI.le_b_i.text())
                Ri_i_space = self.text_to_list(self.wakefieldUI.le_Ri_i.text())
                L_i_space = self.text_to_list(self.wakefieldUI.le_L_i.text())
                Req_i_space = self.text_to_list(self.wakefieldUI.le_Req_i.text())
                alpha_i_space = self.text_to_list(self.wakefieldUI.le_Alpha.text())

                inner_cell_space = [A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space, alpha_i_space]
            else:
                inner_cell_space = [[0], [0], [0], [0], [0], [0], [0], [0]]


            if self.wakefieldUI.cb_Outer_Cell_L.checkState() == 2:
                # Middle Ellipse data
                A_ol_space = self.text_to_list(self.wakefieldUI.le_A_ol.text())
                B_ol_space = self.text_to_list(self.wakefieldUI.le_B_ol.text())
                a_ol_space = self.text_to_list(self.wakefieldUI.le_a_ol.text())
                b_ol_space = self.text_to_list(self.wakefieldUI.le_b_ol.text())
                Ri_ol_space = self.text_to_list(self.wakefieldUI.le_Ri_ol.text())
                L_ol_space = self.text_to_list(self.wakefieldUI.le_L_ol.text())
                Req_ol_space = self.text_to_list(self.wakefieldUI.le_Req_ol.text())
                alpha_ol_space = self.text_to_list(self.wakefieldUI.le_Alpha_ol.text())

                outer_cell_L_space = [A_ol_space, B_ol_space, a_ol_space, b_ol_space, Ri_ol_space, L_ol_space, Req_ol_space, alpha_ol_space]
            else:
                outer_cell_L_space = inner_cell_space

            if self.wakefieldUI.cb_Outer_Cell_R.checkState() == 2:
                # Middle Ellipse data
                A_or_space = self.text_to_list(self.wakefieldUI.le_A_or.text())
                B_or_space = self.text_to_list(self.wakefieldUI.le_B_or.text())
                a_or_space = self.text_to_list(self.wakefieldUI.le_a_or.text())
                b_or_space = self.text_to_list(self.wakefieldUI.le_b_or.text())
                Ri_or_space = self.text_to_list(self.wakefieldUI.le_Ri_or.text())
                L_or_space = self.text_to_list(self.wakefieldUI.le_L_or.text())
                Req_or_space = self.text_to_list(self.wakefieldUI.le_Req_or.text())
                alpha_or_space = self.text_to_list(self.wakefieldUI.le_Alpha_or.text())

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

                                            if self.wakefieldUI.cb_LBP.checkState() == 2 and self.wakefieldUI.cb_RBP.checkState() == 2:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'both', 'FREQ': None}
                                            elif self.wakefieldUI.cb_LBP.checkState() == 2 and self.wakefieldUI.cb_RBP.checkState() == 0:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'left', 'FREQ': None}
                                            elif self.wakefieldUI.cb_LBP.checkState() == 0 and self.wakefieldUI.cb_RBP.checkState() == 2:
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
                                                                        if self.wakefieldUI.cb_LBP.checkState() == 2 and self.wakefieldUI.cb_RBP.checkState() == 0:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                        elif self.wakefieldUI.cb_LBP.checkState() == 0 and self.wakefieldUI.cb_RBP.checkState() == 2:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                        elif self.wakefieldUI.cb_LBP.checkState() == 2 and self.wakefieldUI.cb_RBP.checkState() == 2:
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
                                                                                                if self.wakefieldUI.cb_LBP.checkState() == 2 and self.wakefieldUI.cb_RBP.checkState() == 0:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                                                elif self.wakefieldUI.cb_LBP.checkState() == 0 and self.wakefieldUI.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                                                elif self.wakefieldUI.cb_LBP.checkState() == 2 and self.wakefieldUI.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'both', 'FREQ': None}
                                                                                                else:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'none', 'FREQ': None}

                                                                                                count += 1
            return self.shape_space

    def load_shape_space(self, filename):
        fr = FileReader()
        dir = filename

        # check if extension is included
        if dir.split('.')[-1] != 'json':
            dir = f'{dir}.json'

        df = fr.json_reader(dir)
        # print_(df)

        return df.to_dict()

    def prompt(self, code, fid):
        path = os.getcwd()
        path = os.path.join(path, fr"Data\{code}\Cavity{fid}")
        if os.path.exists(path):
            print_("Data already exists. Do you want to overwrite it?")
            msg = QMessageBox()
            msg.setWindowTitle("Folder Exist")
            msg.setText("Data already exists. Do you want to overwrite it?")
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
            return 'YesToAll'

    def show_hide_(self, wid1, wid2):
        print('here')
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()

    @staticmethod
    def button_clicked(i):
        return i.text()

    def open_file(self, le, cb):
        # clear combobox
        self.cb_Shape_Space_Keys.clear()
        self.cb_Shape_Space_Keys.addItem('All')
        self._selected_keys = []

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

    @staticmethod
    def text_to_list(txt):
        if "range" in txt:
            txt = txt.replace('range', '')
            l = ast.literal_eval(txt)
            return range(l[0], l[1], l[2])
        elif 'linspace' in txt:
            l = eval(f'np.{txt}')
            return l
        else:
            l = ast.literal_eval(txt)
            if isinstance(l, int) or isinstance(l, float):
                return [l]
            else:
                return list(l)

    def animate_width(self, cb, widget, min_width, standard, enable, reverse=False):
        if enable:
            #### GET WIDTH
            width = widget.width()
            #### SET MAX WIDTH

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

            #### ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumWidth")
            self.animation.setDuration(200)
            self.animation.setStartValue(width)
            self.animation.setEndValue(widthCollapsed)
            self.animation.start()

    def animate_height(self, cb, widget, min_height, standard, enable):
        if enable:
            #### GET WIDTH
            height = widget.width()
            #### SET MAX WIDTH
            if cb.checkState() != 2:
                heightCollapsed = min_height
                widget.setMinimumHeight(0)
            else:
                heightCollapsed = standard
                # widget.setMinimumWidth(standard)

            #### ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumHeight")
            self.animation.setDuration(200)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.start()

    def exe_control(self):
        # Abci
        self.wakefieldUI.pb_Top_Drawer.clicked.connect(lambda: self.run_abci_exe(fr'{self.main_control.parentDir}\em_codes\ABCI_exe\TopDrawer for Windows\TopDrawW.exe'))

    def ui_effects(self):
        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.wakefieldUI.w_Settings.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.wakefieldUI.scrollArea.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.wakefieldUI.w_Simulation_Controls.setGraphicsEffect(shadow)

    def serialize(self, state_dict):
        # update state file
        state_dict["Wakefield_Shape_Entry_Mode"] = self.wakefieldUI.cb_Shape_Entry_Mode.currentIndex()
        state_dict["Wakefield_Shape_Space"] = self.wakefieldUI.le_Shape_Space.text()
        state_dict["Wakefield_Mid_Cell_CB"] = self.wakefieldUI.cb_Inner_Cell.checkState()
        state_dict["Wakefield_Left_Cell_CB"] = self.wakefieldUI.cb_Outer_Cell_L.checkState()
        state_dict["Wakefield_Right_Cell_CB"] = self.wakefieldUI.cb_Outer_Cell_R.checkState()
        state_dict["Wakefield_Expansion_CB"] = self.wakefieldUI.cb_Expansion.checkState()
        state_dict["Wakefield_LBP_CB"] = self.wakefieldUI.cb_LBP.checkState()
        state_dict["Wakefield_RBP_CB"] = self.wakefieldUI.cb_RBP.checkState()

        # cell parameters
        state_dict["Wakefield_A_i"] = self.wakefieldUI.le_A_i.text()
        state_dict["Wakefield_B_i"] = self.wakefieldUI.le_B_i.text()
        state_dict["Wakefield_a_i"] = self.wakefieldUI.le_a_i.text()
        state_dict["Wakefield_b_i"] = self.wakefieldUI.le_b_i.text()
        state_dict["Wakefield_Ri_i"] = self.wakefieldUI.le_Ri_i.text()
        state_dict["Wakefield_L_i"] = self.wakefieldUI.le_L_i.text()
        state_dict["Wakefield_Req_i"] = self.wakefieldUI.le_Req_i.text()
        state_dict["Wakefield_Alpha_i"] = self.wakefieldUI.le_Alpha.text()

        state_dict["Wakefield_A_ol"] = self.wakefieldUI.le_A_ol.text()
        state_dict["Wakefield_B_ol"] = self.wakefieldUI.le_B_ol.text()
        state_dict["Wakefield_a_ol"] = self.wakefieldUI.le_a_ol.text()
        state_dict["Wakefield_b_ol"] = self.wakefieldUI.le_b_ol.text()
        state_dict["Wakefield_Ri_ol"] = self.wakefieldUI.le_Ri_ol.text()
        state_dict["Wakefield_L_ol"] = self.wakefieldUI.le_L_ol.text()
        state_dict["Wakefield_Req_ol"] = self.wakefieldUI.le_Req_ol.text()
        state_dict["Wakefield_Alpha_ol"] = self.wakefieldUI.le_Alpha_ol.text()

        state_dict["Wakefield_A_or"] = self.wakefieldUI.le_A_or.text()
        state_dict["Wakefield_B_or"] = self.wakefieldUI.le_B_or.text()
        state_dict["Wakefield_a_or"] = self.wakefieldUI.le_a_or.text()
        state_dict["Wakefield_b_or"] = self.wakefieldUI.le_b_or.text()
        state_dict["Wakefield_Ri_or"] = self.wakefieldUI.le_Ri_or.text()
        state_dict["Wakefield_L_or"] = self.wakefieldUI.le_L_or.text()
        state_dict["Wakefield_Req_or"] = self.wakefieldUI.le_Req_or.text()
        state_dict["Wakefield_Alpha_or"] = self.wakefieldUI.le_Alpha_or.text()

        # settings
        state_dict["Wakefield_N_Cells"] = self.wakefieldUI.sb_N_Cells.value()
        state_dict["Wakefield_N_Modules"] = self.wakefieldUI.sb_N_Modules.value()
        state_dict["Wakefield_Polarization"] = self.wakefieldUI.cb_Polarization_ABCI.currentIndex()
        state_dict["Wakefield_No_Of_Processors"] = self.wakefieldUI.sb_No_Of_Processors_ABCI.value()
        state_dict["Wakefield_Wakelength"] = self.wakefieldUI.le_Wakelength.text()
        state_dict["Wakefield_Bunch_Length"] = self.wakefieldUI.le_Bunch_Length.text()
        state_dict["Wakefield_NFS"] = self.wakefieldUI.le_NFS.text()
        state_dict["Wakefield_MT"] = self.wakefieldUI.le_MT.text()
        state_dict["Wakefield_DDz/SIG"] = self.wakefieldUI.le_DDZ_SIG.text()
        state_dict["Wakefield_DDr/SIG"] = self.wakefieldUI.le_DDR_SIG.text()

    def deserialize(self, state_dict):
        # update state file
        self.wakefieldUI.cb_Shape_Entry_Mode.setCurrentIndex(state_dict["Wakefield_Shape_Entry_Mode"])
        self.wakefieldUI.le_Shape_Space.setText(state_dict["Wakefield_Shape_Space"])
        self.wakefieldUI.cb_Inner_Cell.setCheckState(state_dict["Wakefield_Mid_Cell_CB"])
        self.wakefieldUI.cb_Outer_Cell_L.setCheckState(state_dict["Wakefield_Left_Cell_CB"])
        self.wakefieldUI.cb_Outer_Cell_R.setCheckState(state_dict["Wakefield_Right_Cell_CB"])
        self.wakefieldUI.cb_Expansion.setCheckState(state_dict["Wakefield_Expansion_CB"])
        self.wakefieldUI.cb_LBP.setCheckState(state_dict["Wakefield_LBP_CB"])
        self.wakefieldUI.cb_RBP.setCheckState(state_dict["Wakefield_RBP_CB"])

        # cell parameters
        self.wakefieldUI.le_A_i.setText(state_dict["Wakefield_A_i"])
        self.wakefieldUI.le_B_i.setText(state_dict["Wakefield_B_i"])
        self.wakefieldUI.le_a_i.setText(state_dict["Wakefield_a_i"])
        self.wakefieldUI.le_b_i.setText(state_dict["Wakefield_b_i"])
        self.wakefieldUI.le_Ri_i.setText(state_dict["Wakefield_Ri_i"])
        self.wakefieldUI.le_L_i.setText(state_dict["Wakefield_L_i"])
        self.wakefieldUI.le_Req_i.setText(state_dict["Wakefield_Req_i"])
        self.wakefieldUI.le_Alpha.setText(state_dict["Wakefield_Alpha_i"])

        self.wakefieldUI.le_A_ol.setText(state_dict["Wakefield_A_ol"])
        self.wakefieldUI.le_B_ol.setText(state_dict["Wakefield_B_ol"])
        self.wakefieldUI.le_a_ol.setText(state_dict["Wakefield_a_ol"])
        self.wakefieldUI.le_b_ol.setText(state_dict["Wakefield_b_ol"])
        self.wakefieldUI.le_Ri_ol.setText(state_dict["Wakefield_Ri_ol"])
        self.wakefieldUI.le_L_ol.setText(state_dict["Wakefield_L_ol"])
        self.wakefieldUI.le_Req_ol.setText(state_dict["Wakefield_Req_ol"])
        self.wakefieldUI.le_Alpha_ol.setText(state_dict["Wakefield_Alpha_ol"])
        
        self.wakefieldUI.le_A_or.setText(state_dict["Wakefield_A_or"])
        self.wakefieldUI.le_B_or.setText(state_dict["Wakefield_B_or"])
        self.wakefieldUI.le_a_or.setText(state_dict["Wakefield_a_or"])
        self.wakefieldUI.le_b_or.setText(state_dict["Wakefield_b_or"])
        self.wakefieldUI.le_Ri_or.setText(state_dict["Wakefield_Ri_or"])
        self.wakefieldUI.le_L_or.setText(state_dict["Wakefield_L_or"])
        self.wakefieldUI.le_Req_or.setText(state_dict["Wakefield_Req_or"])
        self.wakefieldUI.le_Alpha_or.setText(state_dict["Wakefield_Alpha_or"])
        
        # settings
        self.wakefieldUI.sb_N_Cells.setValue(state_dict["Wakefield_N_Cells"])
        self.wakefieldUI.sb_N_Modules.setValue(state_dict["Wakefield_N_Modules"])
        self.wakefieldUI.cb_Polarization_ABCI.setCurrentIndex(state_dict["Wakefield_Polarization"])
        self.wakefieldUI.sb_No_Of_Processors_ABCI.setValue(state_dict["Wakefield_No_Of_Processors"])
        self.wakefieldUI.le_Wakelength.setText(state_dict["Wakefield_Wakelength"])
        self.wakefieldUI.le_Bunch_Length.setText(state_dict["Wakefield_Bunch_Length"])
        self.wakefieldUI.le_NFS.setText(state_dict["Wakefield_NFS"])
        self.wakefieldUI.le_MT.setText(state_dict["Wakefield_MT"])
        self.wakefieldUI.le_DDZ_SIG.setText(state_dict["Wakefield_DDz/SIG"])
        self.wakefieldUI.le_DDR_SIG.setText(state_dict["Wakefield_DDr/SIG"])


    @staticmethod
    def run_abci_exe(path):
        path = os.path.join(os.getcwd(), path)
        t = Thread(target=subprocess.call, args=(path,))
        t.start()

    def draw_shape_from_shape_space(self):
        colors = [[48, 162, 218, 255], [252, 79, 48, 255], [229, 174, 56, 255], [109, 144, 79, 255], [139, 139, 139, 255]]
        ci = 0

        # remove existing cells
        self.graphicsView.removeCells()
        for key in self._shape_space.keys():
            if key in self.cb_Shape_Space_Keys.currentText():
                IC = self._shape_space[key]["IC"]
                OC = self._shape_space[key]["OC"]
                BP = self._shape_space[key]["BP"]
                self.graphicsView.drawCells(IC, OC, BP, QColor(colors[ci][0], colors[ci][1], colors[ci][2], colors[ci][3]))

                ci += 1


    @staticmethod
    def run_sequential(n_cells, n_modules, processor_shape_space,
                       MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
                       DDR_SIG=0.1, DDZ_SIG=0.1,
                       parentDir=None, projectDir=None, progress_list=None):
        progress = 0
        # get length of processor
        total_no_of_shapes = len(list(processor_shape_space.keys()))
        for key, shape in processor_shape_space.items():
            # run abci code
            start_time = time.time()
            try:
                abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                 fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                 DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir)
            except:
                abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                 fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                 DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir)

            print_(f'Cavity {key}. Time: {time.time() - start_time}')

            # update progress
            progress_list.append((progress+1)/total_no_of_shapes)


class ProgressMonitor(QThread):
    sig = QtCore.pyqtSignal(int)

    def __init__(self, wake_control, projectDir):
        super(QThread, self).__init__()
        self.wake_control = wake_control
        self.proc_ids = wake_control.processes_id
        self.progress_bar = wake_control.progress_bar
        self.projectDir = projectDir

    def run(self):
        self.progress_monitor(self.proc_ids)

    def progress_monitor(self, proc_ids):
        proc_count = len(proc_ids)
        while self.wake_control.show_progress_bar:
            # print("\t\t\t\t\t\t printing progress:: ", self.wake_control.progress_list)
            progress = sum(self.wake_control.progress_list)/proc_count
            self.sig.emit(round(progress*100, 10))


class EndRoutine(QThread):
    def __init__(self, control, projectDir):
        super(QThread, self).__init__()
        self.control = control
        self.proc_ids = control.processes_id
        self.projectDir = projectDir

    def run(self):
        self.end_routine()

    def end_routine(self):
        for pid in self.proc_ids:
            try:
                p = psutil.Process(pid)
                while p.is_running():
                    pass

                print(fr"process {p} ended")
            except:
                pass

        self.control.cancel()


class CheckableComboBox(QComboBox):

    # Subclass Delegate to increase item height
    class Delegate(QStyledItemDelegate):
        def sizeHint(self, option, index):
            size = super().sizeHint(option, index)
            size.setHeight(20)
            return size

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.lineEdit().setReadOnly(True)
        # Make the lineedit the same color as QPushButton
        palette = qApp.palette()
        palette.setBrush(QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # Use custom delegate
        self.setItemDelegate(CheckableComboBox.Delegate())

        # Update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # Hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)
        self.closeOnLineEditClick = False

        # Prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

    def resizeEvent(self, event):
        # Recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, object, event):

        if object == self.lineEdit():
            if event.type() == QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if object == self.view().viewport():
            if event.type() == QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                item = self.model().item(index.row())

                if item.checkState() == Qt.Checked:
                    item.setCheckState(Qt.Unchecked)

                    if item == self.model().item(0):
                        # deselect all items if item check is all
                        for i in range(1, self.model().rowCount()):
                            item = self.model().item(i)
                            item.setCheckState(Qt.Unchecked)
                else:
                    item.setCheckState(Qt.Checked)

                    if item == self.model().item(0):
                        # deselect all items if item check is all
                        for i in range(1, self.model().rowCount()):
                            item = self.model().item(i)
                            item.setCheckState(Qt.Checked)

                return True
        return False

    def showPopup(self):
        super().showPopup()
        # When the popup is displayed, a click on the lineedit should close it
        self.closeOnLineEditClick = True

    def hidePopup(self):
        super().hidePopup()
        # Used to prevent immediate reopening when clicking on the lineEdit
        self.startTimer(100)
        # Refresh the display text when closing
        self.updateText()

    def timerEvent(self, event):
        # After timeout, kill timer, and reenable click on line edit
        self.killTimer(event.timerId())
        self.closeOnLineEditClick = False

    def updateText(self):
        texts = []
        for i in range(1, self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                texts.append(self.model().item(i).text())
        text = ", ".join(texts)
        self.lineEdit().setText(text)

        # # Compute elided text (with "...")
        # metrics = QFontMetrics(self.lineEdit().font())
        # elidedText = metrics.elidedText(text, Qt.ElideRight, self.lineEdit().width())
        # self.lineEdit().setText(elidedText)

    def addItem(self, text, data=None):
        item = QStandardItem()
        item.setText(text)
        if data is None:
            item.setData(text)
        else:
            item.setData(data)
        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
        item.setData(Qt.Unchecked, Qt.CheckStateRole)
        self.model().appendRow(item)

    def addItems(self, texts, datalist=None):
        for i, text in enumerate(texts):
            try:
                data = datalist[i]
            except (TypeError, IndexError):
                data = None
            self.addItem(text, data)

    def currentData(self):
        # Return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                res.append(self.model().item(i).data())
        return res



    # def run_ABCI_sequential(self, n_cells, n_modules, shape_space,
    #                         MROT=0, MT=4, NFS=10000.0, UBT=50.0, bunch_length=20.0,
    #                         DDR_SIG=0.1, DDZ_SIG=0.1):
    #     try:
    #         start = time.time()
    #         for key, shape in shape_space.items():
    #             print_(f'key, shape: {key, shape}')
    #             # create folder
    #             self.abci_geom.createFolder(key, self.main_control.projectDir)
    #             start = time.time()
    #
    #             try:
    #                 self.abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
    #                                       fid= key, MROT= MROT, MT= MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
    #                                       DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=self.main_control.parentDir, projectDir=self.main_control.projectDir)
    #             except:
    #                 self.abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
    #                                       fid= key, MROT= MROT, MT= MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
    #                                       DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=self.main_control.parentDir, projectDir=self.main_control.projectDir)
    #
    #             end = time.time()
    #             print_(f"Done running {key}, {end-start}s")
    #
    #             print_(f'Runtime/shape: {time.time() - start}')
    #         end = time.time()
    #         runtime = end - start
    #         print_(f'Runtime: {runtime}')
    #     except Exception as e:
    #         print_("ABCI code ran into an exception starting simulation::", e)

    # def run_ABCI_parallel(self, proc_count, n_cells, n_modules, shape_space,
    #                       MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
    #                       DDR_SIG=0.1, DDZ_SIG=0.1):
    #
    #     wakefield_parameters = [MROT, MT, NFS, UBT, bunch_length]
    #     mesh_parameters = [DDR_SIG, DDZ_SIG]
    #     try:
    #         start = time.time()
    #         # run abci code
    #         command = ["mpiexec", "-np", f"{proc_count}", "python", f"{self.main_control.parentDir}/simulation_codes/ABCI/abci_mpi.py",
    #                    f"{n_cells}", f"{n_modules}", f"{shape_space}",
    #                    f"{wakefield_parameters}", f"{mesh_parameters}", f'{self.main_control.parentDir}', f'{self.main_control.projectDir}']
    #
    #         print_(command)
    #         sp = subprocess.run(command, stdout=sys.stdout, stderr=subprocess.STDOUT)
    #         # sp = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #
    #         end = time.time()
    #
    #         runtime = end - start
    #         print_(f"Runtime: {runtime}s")
    #
    #     except Exception as e:
    #         print("Exception occured in ABCI Parallel:: ", e)

    # def run_ABCI_parallel_MP(self, proc_count, n_cells, n_modules, shape_space,
    #                       MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
    #                       DDR_SIG=0.1, DDZ_SIG=0.1):
    #
    #     wakefield_parameters = [MROT, MT, NFS, UBT, bunch_length]
    #     mesh_parameters = [DDR_SIG, DDZ_SIG]
    #
    #     try:
    #         # save shape space as temporary file to be read from the parallel code
    #         shape_space_name = '_temp_shape_space.json'
    #         with open(fr'{self.main_control.projectDir}\Cavities\{shape_space_name}', "w") as outfile:
    #             json.dump(shape_space, outfile, indent=4, separators=(',', ': '))
    #
    #         # run abci code
    #         command = ["python", f"{self.main_control.parentDir}/simulation_codes/ABCI/abci_MP.py",
    #                    f"{n_cells}", f"{n_modules}", f"{shape_space_name}",
    #                    f"{wakefield_parameters}", f"{mesh_parameters}", f"{proc_count}", f'{self.main_control.parentDir}', f'{self.main_control.projectDir}']
    #         print_(command)
    #
    #         start_time = time.time()
    #         self.sp = subprocess.run(command, stdout=sys.stdout, stderr=subprocess.STDOUT)
    #         # sp = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #         print()
    #         print_(f'Total runtime: {time.time() - start_time}')
    #
    #     except Exception as e:
    #         print_(f'Exception encountered in parallel code -> {e}')

