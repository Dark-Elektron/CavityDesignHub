import ast
import json
import os
import subprocess
import time
import multiprocessing as mp
from threading import Thread
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from termcolor import colored
from simulation_codes.SLANS.slans_geom_par import SLANSGeometry
from graphics.graphics_view import GraphicsView
from graphics.scene import Scene
from simulation_codes.SLANS.slans_geometry import SLANSGeometry
from ui_files.eigenmode import Ui_Eigenmode
from utils.file_reader import FileReader
import psutil

# evol = Evolution()
slans_geom = SLANSGeometry()
fr = FileReader()

file_color = 'red'
DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))


class EigenmodeControl:
    def __init__(self, parent):
        self.w_Eigenmode = QWidget()

        self.eigenmodeUI = Ui_Eigenmode()
        self.eigenmodeUI.setupUi(self.w_Eigenmode)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # get logger
        self.log = self.main_control.log

        # ###########################
        #
        # Create Scene
        self.scene = Scene(self)

        # QGraphicsView
        self.graphicsView = GraphicsView(self, 'Eigenmode')
        self.eigenmodeUI.vL_2D_Graphics_View.addWidget(self.graphicsView)
        #
        # ##########################

        self.initUI()
        self.signals()
        self.exe_control()

        # instantiate geometry
        self.slans_geom = SLANSGeometry()

        # shape space initialization
        self._shape_space = {}
        self._selected_keys = []

    def initUI(self):
        # add checkable combobox
        self.cb_Shape_Space_Keys = CheckableComboBox()
        # self.cb_Shape_Space_Keys.addItem('All')
        self.cb_Shape_Space_Keys.setMinimumWidth(75)
        self.eigenmodeUI.gl_Shape_Space_Keys.addWidget(self.cb_Shape_Space_Keys)

        # init shape entry mode
        self.shape_entry_widgets_control()

        # disable expansion section for now. Feature to come later
        self.eigenmodeUI.cb_Expansion.setEnabled(False)

        # inner cell
        self.eigenmodeUI.cb_Inner_Cell.setCheckState(2)
        self.eigenmodeUI.cb_Inner_Cell.setEnabled(False)

        # expand/collapse sections widgets
        if self.eigenmodeUI.cb_Expansion.checkState() == 2:
            self.eigenmodeUI.w_Expansion.setMinimumWidth(300)
        else:
            self.eigenmodeUI.w_Expansion.setMinimumWidth(0)
            self.eigenmodeUI.w_Expansion.setMaximumWidth(0)

        if self.eigenmodeUI.cb_Outer_Cell_L.checkState() == 2:
            self.eigenmodeUI.w_Outer_Cell_L.setMinimumWidth(300)
        else:
            self.eigenmodeUI.w_Outer_Cell_L.setMinimumWidth(0)
            self.eigenmodeUI.w_Outer_Cell_L.setMaximumWidth(0)

        if self.eigenmodeUI.cb_Outer_Cell_R.checkState() == 2:
            self.eigenmodeUI.w_Outer_Cell_R.setMinimumWidth(300)
        else:
            self.eigenmodeUI.w_Outer_Cell_R.setMinimumWidth(0)
            self.eigenmodeUI.w_Outer_Cell_R.setMaximumWidth(0)

        # create pause and resume icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # process state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        # set default boundary condition to magnetic wall at both ends
        self.eigenmodeUI.cb_LBC.setCurrentIndex(2)
        self.eigenmodeUI.cb_RBC.setCurrentIndex(2)

    def signals(self):
        # run eigenmode solver
        self.eigenmodeUI.pb_Run.clicked.connect(lambda: self.run_SLANS())

        # load shape space
        self.eigenmodeUI.pb_Select_Shape_Space.clicked.connect(lambda: self.open_file(self.eigenmodeUI.le_Shape_Space, self.cb_Shape_Space_Keys))

        # control shape entry mode
        self.eigenmodeUI.cb_Shape_Entry_Mode.currentIndexChanged.connect(lambda: self.shape_entry_widgets_control())

        # cell parameters control signals
        self.eigenmodeUI.cb_Outer_Cell_L.stateChanged.connect(lambda: self.animate_width(self.eigenmodeUI.cb_Outer_Cell_L, self.eigenmodeUI.w_Outer_Cell_L, 0, 300, True))
        self.eigenmodeUI.cb_Outer_Cell_R.stateChanged.connect(lambda: self.animate_width(self.eigenmodeUI.cb_Outer_Cell_R, self.eigenmodeUI.w_Outer_Cell_R, 0, 300, True))
        self.eigenmodeUI.cb_Expansion.stateChanged.connect(lambda: self.animate_width(self.eigenmodeUI.cb_Expansion, self.eigenmodeUI.w_Expansion, 0, 300, True))

        # cancel
        self.eigenmodeUI.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.eigenmodeUI.pb_Pause_Resume.clicked.connect(lambda: self.pause() if self.process_state == 'running' else self.resume())

        #
        self.cb_Shape_Space_Keys.currentTextChanged.connect(lambda: self.draw_shape_from_shape_space())

    def shape_entry_widgets_control(self):
        if self.eigenmodeUI.cb_Shape_Entry_Mode.currentIndex() == 0:
            self.main_control.animate_height(self.eigenmodeUI.w_Select_Shape_Space, 0, 50, True)

            self.eigenmodeUI.w_Enter_Geometry_Manual.setMinimumHeight(0)
            self.eigenmodeUI.w_Enter_Geometry_Manual.setMaximumHeight(0)

            # clear cells from graphics view
            self.graphicsView.removeCells()
        else:
            self.main_control.animate_height(self.eigenmodeUI.w_Enter_Geometry_Manual, 0, 375, True)

            self.eigenmodeUI.w_Select_Shape_Space.setMinimumHeight(0)
            self.eigenmodeUI.w_Select_Shape_Space.setMaximumHeight(0)

            # clear cells from graphics view
            self.graphicsView.removeCells()

            # draw new cell
            self.graphicsView.drawCells(color=QColor(0, 0, 0, 255))

    def run_SLANS(self):
        # get analysis parameters
        n_cells = self.eigenmodeUI.sb_N_Cells.value()
        n_modules = self.eigenmodeUI.sb_N_Modules.value()
        f_shift = float(self.eigenmodeUI.le_Freq_Shift.text())
        n_modes = float(self.eigenmodeUI.le_No_Of_Modes.text())

        # boundary conditions
        lbc = self.eigenmodeUI.cb_LBC.currentIndex()+1
        rbc = self.eigenmodeUI.cb_RBC.currentIndex()+1
        bc = 10*lbc + rbc

        proc_count = self.eigenmodeUI.sb_No_Of_Processors_SLANS.value()
        # get geometric parameters
        shape_space = self.get_geometric_parameters('SLANS')

        # split shape_space for different processes/ MPI share process by rank
        keys = list(shape_space.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

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
                # print(f'Processor {p}: {processor_shape_space}')

                service = mp.Process(target=run_sequential, args=(
                n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, self.main_control.parentDir, self.main_control.projectDir))
                service.start()
                self.processes.append(psutil.Process(service.pid))
                # print("Done")

            except Exception as e:
                self.log.error(fr"Exception in run_MP:: {e}")
                # print_("Exception in run_MP::", e)

            self.log.info("Eigenmode simulation started")
            # change process state to running
            self.process_state = 'running'
            self.run_pause_resume_stop_routine()

    def run_pause_resume_stop_routine(self):
        if self.process_state == 'none':
            # change pause/resume icon to pause icon
            self.eigenmodeUI.pb_Pause_Resume.setIcon(self.pause_icon)

            # disable pause/resume and cancel buttons
            self.eigenmodeUI.pb_Pause_Resume.setEnabled(False)
            self.eigenmodeUI.pb_Cancel.setEnabled(False)

            # enable run button in case it was disabled
            self.eigenmodeUI.pb_Run.setEnabled(True)

        if self.process_state == "running":
            # enable run, pause/resume and cancel buttons
            self.eigenmodeUI.pb_Pause_Resume.setEnabled(True)
            self.eigenmodeUI.pb_Cancel.setEnabled(True)
            self.eigenmodeUI.pb_Run.setEnabled(False)

            # change pause/resume icon to pause icon
            self.eigenmodeUI.pb_Pause_Resume.setIcon(self.pause_icon)

        if self.process_state == 'paused':
            # disable run button
            self.eigenmodeUI.pb_Run.setEnabled(False)

            # change pause/resume button icon to resume icon
            self.eigenmodeUI.pb_Pause_Resume.setIcon(self.resume_icon)

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
        try:
            for p in self.processes:
                p.terminate()
        except:
            pass

        self.processes.clear()

        self.process_state = 'none'
        self.run_pause_resume_stop_routine()
        self.log.info("Process terminated.")

    def get_geometric_parameters(self, code):
        self.shape_space = {}
        print_('Getting geometric parameters')
        if self.eigenmodeUI.cb_Shape_Entry_Mode.currentIndex() == 0:
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
            if self.eigenmodeUI.cb_Inner_Cell.checkState() == 2:
                # Middle Ellipse data
                A_i_space = self.text_to_list(self.eigenmodeUI.le_A_i.text())
                B_i_space = self.text_to_list(self.eigenmodeUI.le_B_i.text())
                a_i_space = self.text_to_list(self.eigenmodeUI.le_a_i.text())
                b_i_space = self.text_to_list(self.eigenmodeUI.le_b_i.text())
                Ri_i_space = self.text_to_list(self.eigenmodeUI.le_Ri_i.text())
                L_i_space = self.text_to_list(self.eigenmodeUI.le_L_i.text())
                Req_i_space = self.text_to_list(self.eigenmodeUI.le_Req_i.text())
                alpha_i_space = self.text_to_list(self.eigenmodeUI.le_Alpha.text())

                inner_cell_space = [A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space, alpha_i_space]
            else:
                inner_cell_space = [[0], [0], [0], [0], [0], [0], [0], [0]]


            if self.eigenmodeUI.cb_Outer_Cell_L.checkState() == 2:
                # Middle Ellipse data
                A_ol_space = self.text_to_list(self.eigenmodeUI.le_A_ol.text())
                B_ol_space = self.text_to_list(self.eigenmodeUI.le_B_ol.text())
                a_ol_space = self.text_to_list(self.eigenmodeUI.le_a_ol.text())
                b_ol_space = self.text_to_list(self.eigenmodeUI.le_b_ol.text())
                Ri_ol_space = self.text_to_list(self.eigenmodeUI.le_Ri_ol.text())
                L_ol_space = self.text_to_list(self.eigenmodeUI.le_L_ol.text())
                Req_ol_space = self.text_to_list(self.eigenmodeUI.le_Req_ol.text())
                alpha_ol_space = self.text_to_list(self.eigenmodeUI.le_Alpha_ol.text())

                outer_cell_L_space = [A_ol_space, B_ol_space, a_ol_space, b_ol_space, Ri_ol_space, L_ol_space, Req_ol_space, alpha_ol_space]
            else:
                outer_cell_L_space = inner_cell_space

            if self.eigenmodeUI.cb_Outer_Cell_R.checkState() == 2:
                # Middle Ellipse data
                A_or_space = self.text_to_list(self.eigenmodeUI.le_A_or.text())
                B_or_space = self.text_to_list(self.eigenmodeUI.le_B_or.text())
                a_or_space = self.text_to_list(self.eigenmodeUI.le_a_or.text())
                b_or_space = self.text_to_list(self.eigenmodeUI.le_b_or.text())
                Ri_or_space = self.text_to_list(self.eigenmodeUI.le_Ri_or.text())
                L_or_space = self.text_to_list(self.eigenmodeUI.le_L_or.text())
                Req_or_space = self.text_to_list(self.eigenmodeUI.le_Req_or.text())
                alpha_or_space = self.text_to_list(self.eigenmodeUI.le_Alpha_or.text())

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

                                            if self.eigenmodeUI.cb_LBP.checkState() == 2 and self.eigenmodeUI.cb_RBP.checkState() == 2:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'both', 'FREQ': None}
                                            elif self.eigenmodeUI.cb_LBP.checkState() == 2 and self.eigenmodeUI.cb_RBP.checkState() == 0:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'left', 'FREQ': None}
                                            elif self.eigenmodeUI.cb_LBP.checkState() == 0 and self.eigenmodeUI.cb_RBP.checkState() == 2:
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
                                                                        if self.eigenmodeUI.cb_LBP.checkState() == 2 and self.eigenmodeUI.cb_RBP.checkState() == 0:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                        elif self.eigenmodeUI.cb_LBP.checkState() == 0 and self.eigenmodeUI.cb_RBP.checkState() == 2:
                                                                            self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                        elif self.eigenmodeUI.cb_LBP.checkState() == 2 and self.eigenmodeUI.cb_RBP.checkState() == 2:
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
                                                                                                if self.eigenmodeUI.cb_LBP.checkState() == 2 and self.eigenmodeUI.cb_RBP.checkState() == 0:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                                                elif self.eigenmodeUI.cb_LBP.checkState() == 0 and self.eigenmodeUI.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                                                elif self.eigenmodeUI.cb_LBP.checkState() == 2 and self.eigenmodeUI.cb_RBP.checkState() == 2:
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

    def button_clicked(self, i):
        return i.text()

    def open_file(self, le, cb):
        # clear combobox
        self.cb_Shape_Space_Keys.clear()
        self.cb_Shape_Space_Keys.addItem('All')
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

    def text_to_list(self, txt):
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
        # Slans
        self.eigenmodeUI.pb_Genmsh.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\genmesh2.exe"))
        self.eigenmodeUI.pb_Sl.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\Sl.exe"))
        self.eigenmodeUI.pb_Slansc.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slansc.exe"))
        self.eigenmodeUI.pb_Slansm.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slansm.exe"))
        self.eigenmodeUI.pb_Slanss.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slanss.exe"))
        self.eigenmodeUI.pb_Slansre.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slansre.exe"))
        self.eigenmodeUI.pb_MTFView.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\Mtfview\mtfview.exe"))

    def run_slans_exe(self, path, filename=None):
        path = fr"{self.main_control.parentDir}\em_codes\{path}"
        t = Thread(target=subprocess.call, args=(path,))
        t.start()

    def draw_shape_from_shape_space(self):
        colors = [[48,162,218, 255], [252,79,48, 255], [229,174,56, 255], [109,144,79, 255], [139,139,139, 255]]
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


def run_sequential(n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, parentDir, projectDir):
    for key, shape in processor_shape_space.items():
        try:
            # # create folders for all keys
            slans_geom.createFolder(key, projectDir)

            # run slans code
            start_time = time.time()
            try:
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir)
            except:
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir)

            print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')
        except Exception as e:
            print(f'Error in slans_mpi_mp:: run_sequential -> {e}')

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

    # def run_SLANS_sequential(self, n_cells, n_modules, shape_space, f_shift, n_modes):
    #     try:
    #         start = time.time()
    #         for key, shape in shape_space.items():
    #             print_(f'key, shape: {key, shape}')
    #
    #             # create folder
    #             self.slans_geom.createFolder(key, self.main_control.projectDir)
    #
    #             try:
    #                 self.slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'], fid=key, f_shift=f_shift, beampipes=shape['BP'], n_modes=n_modes+1,
    #                                        parentDir=self.main_control.parentDir, projectDir=self.main_control.projectDir)
    #             except:
    #                 self.slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'], fid=key, f_shift=f_shift, beampipes=shape['BP'], n_modes=n_modes+1,
    #                                        parentDir=self.main_control.parentDir, projectDir=self.main_control.projectDir)
    #
    #
    #             print_(f'Runtime/shape: {time.time() - start}')
    #         end = time.time()
    #         runtime = end - start
    #         print_(f'Runtime: {runtime}')
    #     except Exception as e:
    #         print_("SLANS code ran into an exception starting simulation::", e)
    #
    # def run_SLANS_parallel(self, proc_count, n_cells, n_modules, shape_space, f_shift, n_modes):
    #     try:
    #
    #         # create folders for all keys
    #         for key in shape_space.keys():
    #             self.slans_geom.createFolder(key, self.main_control.projectDir)
    #
    #         # run slans
    #         command = ["mpiexec", "-np", f"{proc_count}", "python", fr"{self.win.parentDir}/SLANS_code/slans_mpi.py",
    #                    f"{n_cells}", f'{n_modules}', f'{shape_space}', f'{f_shift}', f'{n_modes+1}',
    #                    f'{self.main_control.parentDir}', f'{self.main_control.projectDir}']
    #
    #         print_(command)
    #         sp = subprocess.run(command)
    #
    #         # Combining dictionaries
    #         print_('Combining dictionaries')
    #
    #         result = {}
    #         for index in range(proc_count):
    #             with open(f'{self.main_control.parentDir}\Cavities\shape_space{index}.json', "r") as infile:
    #                 result.update(json.load(infile))
    #
    #         # check if extension is included
    #         if self.filename.split('.')[-1] != 'json':
    #             self.filename = f'{self.filename}.json'
    #
    #         with open(f'{self.win.parentDir}/Cavity Population/{self.filename}', "w") as outfile:
    #             json.dump(result, outfile, indent=4, separators=(',', ': '))
    #
    #     except Exception as e:
    #         print_(f'Exception encountered in parallel code -> {e}')

    # def run_SLANS_parallel_MP(self, proc_count, n_cells, n_modules, shape_space, f_shift, n_modes):
    #     try:
    #         # save shape space as temporaty file to be read from the parallel code
    #         shape_space_name = '_temp_shape_space.json'
    #         with open(fr'{self.main_control.projectDir}\Cavities\{shape_space_name}', "w") as outfile:
    #             json.dump(shape_space, outfile, indent=4, separators=(',', ': '))
    #
    #         # run in thread ##change the program not to depend on Ri start position
    #         command = ["python", fr"{self.main_control.parentDir}\simulation_codes\SLANS\slans_mpi_MP.py",
    #                    f"{n_cells}", f'{n_modules}', fr"{shape_space_name}",
    #                    f'{f_shift}', f'{n_modes+1}', fr"{proc_count}",
    #                    fr'{self.main_control.parentDir}', fr'{self.main_control.projectDir}']
    #
    #         start_time = time.time()
    #         sp = subprocess.run(command)
    #         print()
    #         print_(f'Total runtime: {time.time() - start_time}')
    #
    #     except Exception as e:
    #         print_(f'Exception encountered in parallel code -> {e}')


    # def run_SLANS_parallel_MP_1(self, proc_count, n_cells, n_modules, shape_space, f_shift, n_modes):
    #
    #     # get dictionary from json file
    #     dirc = fr'{projectDir}\Cavities\{shape_space_name}'
    #     shape_space = fr.json_reader(dirc)
    #
    #     # split shape_space for different processes/ MPI share process by rank
    #     keys = list(shape_space.keys())
    #     shape_space_len = len(keys)
    #     share = round(shape_space_len / proc_count)
    #
    #     processes = []
    #     for p in range(proc_count):
    #         try:
    #             if p < proc_count-1:
    #                 proc_keys_list = keys[p * share:p * share + share]
    #             else:
    #                 proc_keys_list = keys[p * share:]
    #
    #             processor_shape_space = {}
    #             for key, val in shape_space.items():
    #                 if proc_keys_list[0] <= int(key) <= proc_keys_list[-1]:
    #                     processor_shape_space[key] = val
    #             # print(f'Processor {p}: {processor_shape_space}')
    #
    #             service = mp.Process(target=run_sequential, args=(n_cells, n_modules, processor_shape_space, n_modes, f_shift, parentDir, projectDir))
    #             service.start()
    #             processes.append(service)
    #             # print("Done")
    #
    #         except Exception as e:
    #             print_("Exception in run_MP::", e)
