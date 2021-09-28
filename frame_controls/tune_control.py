import ctypes
import shutil
import stat
from distutils.dir_util import copy_tree
from math import floor
import threading
from PyQt5.QtCore import QPropertyAnimation
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtWidgets import *
from scipy.optimize import fsolve
from modules.tune_module.tuners.tuner import Tuner
from simulation_codes.SLANS.slans_tune import Tune
from ui_files.run_tune import Ui_w_Tune
import ast
import json
import os
import time
import multiprocessing as mp
from termcolor import colored
from PyQt5.QtWidgets import *
import random as r
from simulation_codes.SLANS.slans_geom_par import SLANSGeometry
import numpy as np
from utils.file_reader import FileReader
import pyqtgraph as pg
from distutils import dir_util
import psutil

file_color = 'red'
DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))

fr = FileReader()
slans_geom = SLANSGeometry()

AN_DURATION = 200


class TuneControl:
    def __init__(self, parent):
        self.w_Tune = QWidget()

        self.tuneUI = Ui_w_Tune()
        self.tuneUI.setupUi(self.w_Tune)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # get logger
        self.log = self.main_control.log

        self.initUI()
        self.signals()

        self.pg_list = []
        self.pygraph_list = []
        self.tune_ended = False
        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False

    def initUI(self):
        self.tuneUI.cb_Inner_Cell.setCheckState(2)
        self.tuneUI.cb_Inner_Cell.setEnabled(False)

        # expand/collapse sections widgets
        if self.tuneUI.cb_Expansion.checkState() == 2:
            self.tuneUI.w_Expansion.setMinimumWidth(300)
        else:
            self.tuneUI.w_Expansion.setMinimumWidth(0)
            self.tuneUI.w_Expansion.setMaximumWidth(0)

        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
            self.tuneUI.w_Outer_Cell.setMinimumWidth(300)
        else:
            self.tuneUI.w_Outer_Cell.setMinimumWidth(0)
            self.tuneUI.w_Outer_Cell.setMaximumWidth(0)

        self.tuneUI.w_Tune_Settings.setMinimumHeight(0)
        self.tuneUI.w_Tune_Settings.setMaximumHeight(0)

        # disable expansion section for now. Feature to come later
        self.tuneUI.cb_Expansion.setEnabled(False)

        # tuner initial state
        self.tuner_routine()
        if self.tuneUI.cb_Tuner.currentText() == 'SLANS':
            self.tuneUI.w_Iteration_Settings.setEnabled(False)

        # initial loaded shape spaces
        self.loaded_inner_dict = {}
        self.loaded_outer_dict = {}

        # initial show/hide
        self.show_hide()

        # show/hide convergence plot
        self.animate_width(self.tuneUI.cb_Monitor_Convergence, self.tuneUI.w_PyqtGraph, 0, 375, True)

        # create pause and resume icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # # initially set pause icon
        # self.tuneUI.pb_Pause_Resume.setIcon(self.pause_icon)

        # process state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        # set default boundary condition to magnetic wall at both ends
        self.tuneUI.cb_LBC.setCurrentIndex(2)
        self.tuneUI.cb_RBC.setCurrentIndex(2)

    def signals(self):

        # signals/slots
        self.tuneUI.cb_Shape_Space_Generation_Algorithm.currentIndexChanged.connect(lambda: self.show_hide())
        self.tuneUI.pb_Run.clicked.connect(lambda: self.generate_shape_space())

        # tune variable
        self.tuneUI.cb_Tune_Option.currentTextChanged.connect(lambda: self.change_tune_option(self.tuneUI.cb_Tune_Option.currentText()))

        # cell parameters control
        self.tuneUI.cb_Outer_Cell.stateChanged.connect(lambda: self.animate_width(self.tuneUI.cb_Outer_Cell, self.tuneUI.w_Outer_Cell, 0, 300, True))
        self.tuneUI.cb_Expansion.stateChanged.connect(lambda: self.animate_width(self.tuneUI.cb_Expansion, self.tuneUI.w_Expansion, 0, 300, True))
        self.tuneUI.cb_Tuner.currentTextChanged.connect(lambda: self.tuner_routine())
        self.tuneUI.cb_Cell_Type.currentTextChanged.connect(lambda: self.tuner_routine())
        self.tuneUI.cb_Tune_Option.currentTextChanged.connect(lambda: self.tuner_routine())

        # change cavity display image
        self.tuneUI.cb_Cell_Type.currentTextChanged.connect(lambda: self.change_cavity_image())
        self.tuneUI.cb_LBP.stateChanged.connect(lambda: self.change_cavity_image())
        self.tuneUI.cb_Outer_Cell.stateChanged.connect(lambda: self.change_cavity_image())

        # collapse other settings
        self.tuneUI.pb_Tune_Settings.clicked.connect(lambda: self.main_control.animate_height(self.tuneUI.w_Tune_Settings, 0, 200, True))

        # disable secondary settings if Tuner is SLANS
        self.tuneUI.cb_Tuner.currentTextChanged.connect(lambda: self.tuneUI.w_Iteration_Settings.setEnabled(False) if self.tuneUI.cb_Tuner.currentText() == 'SLANS'
                                                        else self.tuneUI.w_Iteration_Settings.setEnabled(True))

        # control to ensure that SLANS mid cell tuner is always set to tune for Req and end cell tuner is always set to L
        self.tuneUI.cb_Tuner.currentTextChanged.connect(lambda: self.slans_tuners_control())
        self.tuneUI.cb_Cell_Type.currentTextChanged.connect(lambda: self.slans_tuners_control())

        # expand/collapse convergence plot
        self.tuneUI.cb_Monitor_Convergence.stateChanged.connect(lambda: self.animate_width(self.tuneUI.cb_Monitor_Convergence, self.tuneUI.w_PyqtGraph, 0, 375, True))

        # cancel
        self.tuneUI.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.tuneUI.pb_Pause_Resume.clicked.connect(lambda: self.pause() if self.process_state == 'running' else self.resume())

        # change variable value
        self.tuneUI.le_Tune_Variable.textChanged.connect(lambda: self.tuneUI.le_Tune_Variable_End_Cell.setText(self.tuneUI.le_Tune_Variable.text())
                                                         if self.tuneUI.cb_Tune_Option.currentIndex() == 1
                                                         else None)

    def createPyGraph(self, row, column):
        print(row, column)
        # add pyqgraph
        pg_ = pg
        pygraph = pg_.PlotWidget()
        pygraph.setBackground('w')
        print("here here")
        self.tuneUI.gl_PyqtGraph.addWidget(pygraph, row, column, 1, 1)

        self.pg_list.append(pg_)
        self.pygraph_list.append(pygraph)

    def monitor_convergence(self):
        print("in here monitor")
        self.plot_list = [None for i in range(self.tuneUI.sb_No_Of_Processors_Tune.value())]

        for p in self.pygraph_list:
            if p:
                p.clear()

        while not self.tune_ended:
            # if update_count%update_interval == 0:
            # print("In while loop")
            # read data from file
            i = 0
            for pg, pygraph in zip(self.pg_list, self.pygraph_list):
                filename = f"{self.main_control.projectDir}\SimulationData\SLANS\Cavity_process_{i}\convergence_output.json"

                if os.path.exists(filename):
                    df = fr.json_reader(filename)
                    L_list, freq_list = df["tv"].to_list(), df["freq"].to_list()

                    if self.tuneUI.cb_Monitor_Convergence.checkState() == 2:
                        pygraph.addLine(x=None, y=400.79, pen=pg.mkPen('r', width=1))

                        if not self.plot_list[i]:
                            self.plot_list[i] = pygraph.plot(L_list, freq_list, pen=None, symbol='o')
                        else:
                            print("Plotted calues", i)
                            print("\t\t", L_list)
                            print("\t\t", freq_list)
                            self.plot_list[i].setData(L_list, freq_list, pen=None, symbol='o')

                i += 1

            time.sleep(1)

        self.tune_ended = False

    def tuner_routine(self):
        if self.tuneUI.cb_Cell_Type.currentText() == 'Mid Cell':
            self.tuneUI.cb_Outer_Cell.setCheckState(0)
            self.tuneUI.cb_Outer_Cell.setEnabled(False)
            self.tuneUI.cb_LBP.setCheckState(0)
            self.tuneUI.cb_LBP.setEnabled(False)
            self.tuneUI.w_End_Cell_Tune_Extra_Variable.hide()
        elif self.tuneUI.cb_Cell_Type.currentText() == 'End Cell':
            self.tuneUI.cb_Outer_Cell.setCheckState(2)
            self.tuneUI.cb_Outer_Cell.setEnabled(False)
            self.tuneUI.cb_LBP.setCheckState(2)
            self.tuneUI.cb_LBP.setEnabled(False)

            if self.tuneUI.cb_Tune_Option.currentText() == 'L':
                self.tuneUI.w_End_Cell_Tune_Extra_Variable.show()
            else:
                self.tuneUI.w_End_Cell_Tune_Extra_Variable.hide()
        else:
            self.tuneUI.cb_LBP.setEnabled(True)
            self.tuneUI.cb_Outer_Cell.setEnabled(True)

    def generate_shape_space(self):
        # check if filename is entered or not
        if self.tuneUI.le_Generated_Shape_Space_Name.text() == '':
            self.log.error("Hey chief, seems you forgot to give a name to the shape_space file.")
            # print_("Hey chief, seems you forgot to give a name to the shape_space file.")
            return
        else:
            self.filename = self.proof_filename(self.tuneUI.le_Generated_Shape_Space_Name.text())
            # check if shape space already generated
            resume = self.continue_check()

            # check if pseudo_shape_space is to be updated
            self.existing_keys = []
            self.pseudo_shape_space = {}
            if resume == "Yes":
                # check if value set is already written. This is to enable continuation in case of break in program
                if os.path.exists(f'{self.main_control.projectDir}\Cavities\pseudo_{self.filename}'):
                    self.pseudo_shape_space = json.load(open(fr'{self.main_control.projectDir}\Cavities\pseudo_{self.filename}', 'r'))

                    self.existing_keys = list(self.pseudo_shape_space.keys())
                    # self.log(f'last saved key: {self.existing_keys}')
                    print("Existing keys", self.existing_keys)
            elif resume == "Cancel":
                return

        freq = float(self.tuneUI.le_Freq.text())
        # get variables from ui or from pseudo shape space
        A_i = self.text_to_list(self.tuneUI.le_A_i.text())
        B_i = self.text_to_list(self.tuneUI.le_B_i.text())
        a_i = self.text_to_list(self.tuneUI.le_a_i.text())
        b_i = self.text_to_list(self.tuneUI.le_b_i.text())
        Ri_i = self.text_to_list(self.tuneUI.le_Ri_i.text())

        # check tune type
        if self.tuneUI.cb_Tune_Option.currentText() == "L":
            Req_i = self.text_to_list(self.tuneUI.le_Tune_Variable.text())
            L_i = None
            inner_half_cell_parameters = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i]
        else:
            L_i = self.text_to_list(self.tuneUI.le_Tune_Variable.text())
            Req_i = None
            inner_half_cell_parameters = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i]

        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
            A_o = self.text_to_list(self.tuneUI.le_A_o.text())
            B_o = self.text_to_list(self.tuneUI.le_B_o.text())
            a_o = self.text_to_list(self.tuneUI.le_a_o.text())
            b_o = self.text_to_list(self.tuneUI.le_b_o.text())
            Ri_o =self.text_to_list(self.tuneUI.le_Ri_o.text())

            if self.tuneUI.cb_Tune_Option.currentText() == 'L':
                # update mid cell L_i
                inner_half_cell_parameters[5] = self.text_to_list(self.tuneUI.le_L_Mid_Cell.text())
                L_o = None
            else:
                L_o = self.text_to_list(self.tuneUI.le_Tune_Variable_End_Cell.text())

            Req_o = Req_i
            outer_half_cell_parameters = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o]
        else:
            outer_half_cell_parameters = inner_half_cell_parameters

        print(inner_half_cell_parameters, outer_half_cell_parameters)
        pseudo_shape_space = self.generate_pseudo_shape_space(freq, inner_half_cell_parameters, outer_half_cell_parameters)

        if pseudo_shape_space:
            self.run_tune(pseudo_shape_space, resume)

    def run_tune(self, pseudo_shape_space, resume):
        # set tuner
        self.tuner = self.tuneUI.cb_Tuner.currentText()
        proc_count = self.tuneUI.sb_No_Of_Processors_Tune.value()
        tune_variable = self.tuneUI.cb_Tune_Option.currentText().split(' ')[-1]

        # get iteration settings
        iter_method = self.tuneUI.cb_Iterative_Method.currentText()
        tolerance = self.text_to_list(self.tuneUI.le_Tolerance.text())[0]
        max_iter = self.text_to_list(self.tuneUI.sb_Max_Iteration.text())[0]
        iter_set = [iter_method, tolerance, max_iter]

        # cell type
        cell_type = self.tuneUI.cb_Cell_Type.currentText()

        # boundary conditions
        lbc = self.tuneUI.cb_LBC.currentIndex()+1
        rbc = self.tuneUI.cb_RBC.currentIndex()+1
        bc = 10*lbc + rbc
        print("1")

        # try:
        if True:
            self.start = time.time()

            # split shape_space for different processes/ MPI share process by rank
            keys = list(pseudo_shape_space.keys())

            # check if number of processors selected is greater than the number of keys in the pseudo shape space
            if proc_count > len(keys):
                proc_count = len(keys)
                print(f"Changed processor count to {proc_count}.")

            shape_space_len = len(keys)
            share = floor(shape_space_len / proc_count)

            for i in reversed(range(self.tuneUI.gl_PyqtGraph.count())):
                self.tuneUI.gl_PyqtGraph.itemAt(i).widget().setParent(None)
            print("\t1b")

            # insert graphs for convergence monitor
            # divide processors into ratio 1 to 3
            col_count = int(proc_count/3)
            row_count = proc_count - col_count

            print("row column ", row_count, col_count, proc_count)

            self.pg_list = []
            self.pygraph_list = []
            print("pglist before ", self.pg_list)
            print("pygraph before ", self.pygraph_list)

            if col_count == 0:
                for i in range(row_count):
                    self.createPyGraph(i, 0)
            else:
                for i in range(row_count):
                    for j in range(col_count):
                        self.createPyGraph(i, j)

            print("pglist after ", self.pg_list)
            print("pygraph after ", self.pygraph_list)

            # create progress bar object and add to widget
            self.progress_bar = QProgressBar(self.tuneUI.widget_4)
            self.progress_bar.setMaximum(100)
            self.progress_bar.setValue(10)
            self.tuneUI.gridLayout_10.addWidget(self.progress_bar, 0, 7, 1, 1)

            for p in range(proc_count):
                if True:
                    if p < proc_count-1:
                        proc_keys_list = keys[p * share:p * share + share]
                    else:
                        proc_keys_list = keys[p * share:]

                    print("\t1c", p, proc_keys_list)

                    self.overwriteFolder(p, self.main_control.projectDir)
                    print("\t\t1ci")
                    self.copyFiles(p, self.main_control.parentDir, self.main_control.projectDir)
                    print("\t1d")

                    processor_shape_space = {}
                    for key, val in pseudo_shape_space.items():
                        if key in proc_keys_list:
                            processor_shape_space[key] = val
                    # print(f'Processor {p}: {processor_shape_space}')
                    print("\t1e:: ", resume, p, bc, self.main_control.parentDir,
                                               self.main_control.projectDir, self.filename, self.tuner,
                                               tune_variable, iter_set, cell_type)

                    service = mp.Process(target=self.run_sequential,
                                         args=(processor_shape_space, resume, p, bc, self.main_control.parentDir,
                                               self.main_control.projectDir, self.filename, self.tuner,
                                               tune_variable, iter_set, cell_type))
                    print("\t1d")
                    service.start()
                    self.processes.append(psutil.Process(service.pid))
                    self.processes_id.append(service.pid)

                # except Exception as e:
                #     self.log.info("Exception in run_MP::", e)
                #     # write to log

            # display progress bar
            self.show_progress_bar = True

            t_progress_monitor = threading.Thread(target=self.progress_monitor, args=(self.processes_id, ))
            t_progress_monitor.start()

            # start convergence monitor
            t_monitor_convergence = threading.Thread(target=self.monitor_convergence)
            t_monitor_convergence.start()

            self.log.info("Tune started")
            # change process state to running
            self.process_state = 'running'
            self.run_pause_resume_stop_routine()

            # initiate end routine
            filename = self.filename
            projectDir = self.main_control.projectDir

            t_end_routine = threading.Thread(target=self.end_routine, args=(self.processes_id, filename, projectDir))
            t_end_routine.start()

        # except Exception as e:
        #     self.log.error(fr"TUNE CONTROL:: run_tune:: {e}")
        #     print_("TUNE CONTROL:: run_tune:: ", e)

    def run_pause_resume_stop_routine(self):
        if self.process_state == 'none':
            # change pause/resume icon to pause icon
            self.tuneUI.pb_Pause_Resume.setIcon(self.pause_icon)

            # disable pause/resume and cancel buttons
            self.tuneUI.pb_Pause_Resume.setEnabled(False)
            self.tuneUI.pb_Cancel.setEnabled(False)

            # enable run button in case it was disabled
            self.tuneUI.pb_Run.setEnabled(True)

        if self.process_state == "running":
            # enable run, pause/resume and cancel buttons
            self.tuneUI.pb_Pause_Resume.setEnabled(True)
            self.tuneUI.pb_Cancel.setEnabled(True)
            self.tuneUI.pb_Run.setEnabled(False)

            # change pause/resume icon to pause icon
            self.tuneUI.pb_Pause_Resume.setIcon(self.pause_icon)

        if self.process_state == 'paused':
            # disable run button
            self.tuneUI.pb_Run.setEnabled(False)

            # change pause/resume button icon to resume icon
            self.tuneUI.pb_Pause_Resume.setIcon(self.resume_icon)

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

        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        self.log.info("Process terminated.")

        # set exit signal
        self.tune_ended = True

    def end_routine(self, proc_ids, filename, projectDir):
        proc_count = len(proc_ids)
        for pid in proc_ids:
            try:
                p = psutil.Process(pid)
                while p.is_running():
                    pass

                print(fr"process {p} ended")
            except:
                pass

        # combine dictionaries
        try:
            self.combine_dict(proc_count, filename, projectDir)
            self.delete_process_dict(proc_count, projectDir)
        except Exception as e:
            self.log.error(f"Some error occurred -> {e}")

        self.cancel()

    def progress_monitor(self, proc_ids):
        # read progress files and update progress
        while self.show_progress_bar:
            try:
                progress = 0
                for i in range(len(proc_ids)):
                    if os.path.exists(fr'{self.main_control.projectDir}\SimulationData\SLANS\Cavity_process_{i}\progress_file.txt'):
                        with open(fr'{self.main_control.projectDir}\SimulationData\SLANS\Cavity_process_{i}\progress_file.txt', "r") as f:
                            a = f.readline()
                            progress += eval(a)
                self.progress_bar.setValue(progress*100/len(proc_ids))
            except:
                print("Error in progress update")
                pass

        try:
            # remove progress bar
            self.tuneUI.gridLayout_10.removeWidget(self.progress_bar)
            self.progress_bar = None
        except:
            pass

    def change_cell_parameters(self, d, key, par):
        try:
            if par == 'inner':
                self.tuneUI.le_A_i.setText(fr"{d[key]['IC'][0]}")
                self.tuneUI.le_B_i.setText(fr"{d[key]['IC'][1]}")
                self.tuneUI.le_a_i.setText(fr"{d[key]['IC'][2]}")
                self.tuneUI.le_b_i.setText(fr"{d[key]['IC'][3]}")
                self.tuneUI.le_Ri_i.setText(fr"{d[key]['IC'][4]}")
                # check tune variable
                if self.tuneUI.cb_Tune_Option.currentText() == 'Req':
                    self.tuneUI.le_Tune_Variable.setText(fr"{d[key]['IC'][5]}")
                else:
                    self.tuneUI.le_Tune_Variable.setText(fr"{d[key]['IC'][6]}")
            else:
                self.tuneUI.le_A_o.setText(fr"{d[key]['OC'][0]}")
                self.tuneUI.le_B_o.setText(fr"{d[key]['OC'][1]}")
                self.tuneUI.le_a_o.setText(fr"{d[key]['OC'][2]}")
                self.tuneUI.le_b_o.setText(fr"{d[key]['OC'][3]}")
                self.tuneUI.le_Ri_o.setText(fr"{d[key]['OC'][4]}")
                self.tuneUI.le_Tune_Variable_End_Cell.setText(fr"{d[key]['OC'][5]}")
        except Exception as e:
            print('Exception: ', e)

    def generate_pseudo_shape_space(self, freq, ihc, ohc):
        if self.tuneUI.cb_Shape_Space_Generation_Algorithm.currentText() == 'Monte Carlo':
            n_shapes = self.tuneUI.sb_No_Of_Shapes_Monte_Carlo.value()
            # pseudo_shape_space = {}

            # check input to avoid an infinite loop
            check = self.check_input()
            if check:
                count = 0
                loop_escape = 0
                while count < n_shapes and loop_escape < 15:
                    # mid cell
                    A_i = self.process_range(ihc[0])
                    B_i = self.process_range(ihc[1])
                    a_i = self.process_range(ihc[2])
                    b_i = self.process_range(ihc[3])
                    Ri_i =self.process_range(ihc[4])

                    if self.tuneUI.cb_Tune_Option.currentText() == 'L':
                        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
                            L_i = self.process_range(ihc[5])
                        else:
                            L_i = A_i + a_i

                        Req_i = self.process_range(ihc[6])

                    else:
                        L_i = self.process_range(ihc[5])
                        Req_i = B_i + b_i + Ri_i

                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                    print(inner_cell)

                    if self.tuneUI.cb_Outer_Cell.checkState() == 2:
                        # This also checks if the right and left bounds are equal in which case it returns a single value
                        A_o = self.process_range(ohc[0])
                        B_o = self.process_range(ohc[1])
                        a_o = self.process_range(ohc[2])
                        b_o = self.process_range(ohc[3])
                        Ri_o =self.process_range(ohc[4])

                        if self.tuneUI.cb_Tune_Option.currentText() == 'L':
                            L_o = A_o + a_o
                        else:
                            L_o = self.process_range(ohc[5])

                        Req_o = Req_i

                        # update L_i
                        # L_i = self.process_range(self.text_to_list(self.tuneUI.le_Tune_Variable_End_Cell.text()))
                        # inner_cell[5] = L_i

                        other_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, 0]
                        print(other_cell)
                    else:
                        other_cell = inner_cell

                    if other_cell[0] + other_cell[2] > other_cell[5]:
                        loop_escape += 1
                        print_("Encountered a shape with A + a > L")
                        continue

                    if other_cell[3] + other_cell[4] > other_cell[6] or inner_cell[3] + inner_cell[4] > inner_cell[6]:
                        loop_escape += 1
                        print_("Encountered a shape with B + b > Req")
                        continue

                    # print_('6e')
                    if self.tuneUI.cb_LBP.checkState() == 2:
                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                        if key not in self.existing_keys:
                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': other_cell, 'BP': 'left', 'FREQ': freq}
                    else:
                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                        if key not in self.existing_keys:
                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': other_cell, 'BP': 'none', 'FREQ': freq}
                    count += 1
                    loop_escape = 0 # reset loop escape
            else:
                print('Please check input parameters.')
        else:
            ########## UPDATE CODE TO CHECK FOR ALL THE ITEMS IN THE SHAPE SPACE FOR SHAPES WITH A + a > L
            print_("Grid shape space generate")

            check = self.check_input()

            print_('6c')
            A_o_space = ohc[0]
            B_o_space = ohc[1]
            a_o_space = ohc[2]
            b_o_space = ohc[3]
            Ri_o_space =ohc[4]

            count = 0
            print("Ao", A_o_space)
            for A_o in A_o_space:
                for B_o in B_o_space:
                    for a_o in a_o_space:
                        for b_o in b_o_space:
                            for Ri_o in Ri_o_space:

                                if self.tuneUI.cb_Tune_Option.currentText() == 'Req':
                                    L_o_space = ohc[5]
                                    for L_o in L_o_space:
                                        Req_o = B_o + b_o + Ri_o
                                        outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, 0]

                                        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
                                            # This also checks if the right and left bounds are equal in which case it returns a single value
                                            A_i_space = ihc[0]
                                            B_i_space = ihc[1]
                                            a_i_space = ihc[2]
                                            b_i_space = ihc[3]
                                            Ri_i_space =ihc[4]
                                            L_i_space = ihc[5]
                                            Req_i = Req_o

                                            for A_i in A_i_space:
                                                for B_i in B_i_space:
                                                    for a_i in a_i_space:
                                                        for b_i in b_i_space:
                                                            for Ri_i in Ri_i_space:
                                                                for L_i in L_i_space:
                                                                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                                                    if self.tuneUI.cb_LBP.checkState() == 2:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'left', 'FREQ': freq}
                                                                    else:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'none', 'FREQ': freq}
                                                                    count += 1
                                        else:
                                            inner_cell = outer_cell

                                            if self.tuneUI.cb_LBP.checkState() == 2:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'left', 'FREQ': freq}
                                            else:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'none', 'FREQ': freq}
                                            count += 1

                                else:
                                    Req_o_space = ohc[6]
                                    print('1', Req_o_space)
                                    for Req_o in Req_o_space:
                                        print('2')
                                        L_o = A_o + a_o
                                        outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, 0]
                                        print('3')

                                        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
                                            A_i_space = ihc[0]
                                            print('5', A_i_space)
                                            B_i_space = ihc[1]
                                            print('6')
                                            a_i_space = ihc[2]
                                            print('7')
                                            b_i_space = ihc[3]
                                            print('8')
                                            Ri_i_space = ihc[4]
                                            print('9', ihc)
                                            L_i_space = ihc[5]
                                            print('10')
                                            Req_i = Req_o
                                            print('11')

                                            for A_i in A_i_space:
                                                for B_i in B_i_space:
                                                    for a_i in a_i_space:
                                                        for b_i in b_i_space:
                                                            for Ri_i in Ri_i_space:
                                                                for L_i in L_i_space:
                                                                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                                                    print('4')

                                                                    if self.tuneUI.cb_LBP.checkState() == 2:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'left', 'FREQ': freq}
                                                                    else:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'none', 'FREQ': freq}
                                                                    count += 1
                                        else:
                                            inner_cell = outer_cell

                                            if self.tuneUI.cb_LBP.checkState() == 2:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'left', 'FREQ': freq}
                                            else:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell, 'BP': 'none', 'FREQ': freq}
                                            count += 1

        if check:
            print_(f"pseudo shape space check: {self.pseudo_shape_space}")

            pseudo_shape_space_name = f'{self.main_control.projectDir}/Cavities/pseudo_{self.proof_filename(self.tuneUI.le_Generated_Shape_Space_Name.text())}'
            with open(pseudo_shape_space_name, 'w') as file:
                file.write(json.dumps(self.pseudo_shape_space, indent=4, separators=(',', ': ')))

            return self.pseudo_shape_space
        else:
            return check

    def show_hide(self):
        if self.tuneUI.cb_Shape_Space_Generation_Algorithm.currentText() == 'Monte Carlo':
            self.tuneUI.w_No_Of_Shapes_Monte_Carlo.show()
        else:
            self.tuneUI.w_No_Of_Shapes_Monte_Carlo.hide()

        if self.tuneUI.cb_Tuner.currentText() == 'PyTune':
            self.tuneUI.w_BC.show()
        else:
            self.tuneUI.w_BC.hide()

    def animate_width(self, cb, widget, min_width, standard, enable):
        if enable:
            #### GET WIDTH
            width = widget.width()
            #### SET MAX WIDTH
            if cb.checkState() != 2:
                widthCollapsed = min_width
                widget.setMinimumWidth(0)
            else:
                widthCollapsed = standard
                # widget.setMinimumWidth(standard)

            #### ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumWidth")
            self.animation.setDuration(AN_DURATION)
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
            self.animation.setDuration(AN_DURATION)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.start()

    def change_cavity_image(self):
        if self.tuneUI.cb_Cell_Type.currentText() == 'Mid Cell':
            self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/mid_cell.png"))
        elif self.tuneUI.cb_Cell_Type.currentText() == 'End Cell':
            self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/end_cell.png"))
        else:
            if self.tuneUI.cb_LBP.checkState() == 0 and self.tuneUI.cb_Outer_Cell.checkState() == 0:
                self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/mid_cell.png"))
            elif self.tuneUI.cb_LBP.checkState() == 2 and self.tuneUI.cb_Outer_Cell.checkState() == 2:
                self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/end_cell.png"))
            elif self.tuneUI.cb_LBP.checkState() == 0 and self.tuneUI.cb_Outer_Cell.checkState() == 2:
                self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/one_cell.png"))
            elif self.tuneUI.cb_LBP.checkState() == 2 and self.tuneUI.cb_Outer_Cell.checkState() == 0:
                self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/end_same_cell.png"))

    def check_input(self):
        x = True
        A_i = self.text_to_list(self.tuneUI.le_A_i.text())
        a_i = self.text_to_list(self.tuneUI.le_a_i.text())
        B_i = self.text_to_list(self.tuneUI.le_B_i.text())
        b_i = self.text_to_list(self.tuneUI.le_b_i.text())
        Ri_i = self.text_to_list(self.tuneUI.le_Ri_i.text())
        tune_var = self.text_to_list(self.tuneUI.le_Tune_Variable.text())

        if self.tuneUI.cb_Tune_Option.currentText() == "Req":
            L = tune_var
            if A_i[0] + a_i[0] > L[0]:
                x = False

            if len(L) > 1:
                if L[1] < L[0]:
                    x = False
        else:
            Req = tune_var
            if b_i[0] + Ri_i[0] > Req[0]:
                x = False

            if len(Req) > 1:
                Req[0] = b_i[0] + Ri_i[0]
                if Req[1] < Req[0]:
                    x = False

        return x

    def slans_tuners_control(self):
        if self.tuneUI.cb_Tuner.currentText() == 'SLANS':
            if self.tuneUI.cb_Cell_Type.currentText() == 'Mid Cell':
                self.tuneUI.cb_Tune_Option.setCurrentIndex(0)
                self.tuneUI.l_Tune_Alternate_Variable.setText('L')
                self.tuneUI.cb_Tune_Option.setEnabled(False)
            else:
                self.tuneUI.cb_Tune_Option.setCurrentIndex(1)
                self.tuneUI.l_Tune_Alternate_Variable.setText('Req')
                self.tuneUI.cb_Tune_Option.setEnabled(False)

            # disable monitor convergence
            if self.tuneUI.cb_Monitor_Convergence.checkState() == 2:
                self.tuneUI.cb_Monitor_Convergence.setCheckState(0)
            self.tuneUI.cb_Monitor_Convergence.setEnabled(False)

            # hide boundary condition
            self.tuneUI.w_BC.hide()
        else:
            self.tuneUI.cb_Tune_Option.setEnabled(True)

            # enable monitor convergence
            self.tuneUI.cb_Monitor_Convergence.setEnabled(True)

            # show boundary conditions
            self.tuneUI.w_BC.show()

    def serialize(self, state_dict):
        # update state file
        state_dict["Frequency"] = self.tuneUI.le_Freq.text()
        state_dict["Cell_Type"] = self.tuneUI.cb_Cell_Type.currentIndex()
        state_dict["Tune_Option"] = self.tuneUI.cb_Tune_Option.currentIndex()
        state_dict["Method"] = self.tuneUI.cb_Shape_Space_Generation_Algorithm.currentIndex()

        state_dict["No_Of_Shapes_Monte_Carlo"] = self.tuneUI.sb_No_Of_Shapes_Monte_Carlo.value()

        state_dict["Tuner"] = self.tuneUI.cb_Tuner.currentIndex()
        state_dict["LBC"] = self.tuneUI.cb_LBC.currentIndex()
        state_dict["RBC"] = self.tuneUI.cb_RBC.currentIndex()

        state_dict["Inner_Cell"] = self.tuneUI.cb_Inner_Cell.checkState()
        state_dict["Outer_Cell"] = self.tuneUI.cb_Outer_Cell.checkState()
        state_dict["Expansion"] = self.tuneUI.cb_Expansion.checkState()
        state_dict["LBP"] = self.tuneUI.cb_LBP.checkState()

        # cell parameters
        state_dict["A_i"] = self.tuneUI.le_A_i.text()
        state_dict["B_i"] = self.tuneUI.le_B_i.text()
        state_dict["a_i"] = self.tuneUI.le_a_i.text()
        state_dict["b_i"] = self.tuneUI.le_b_i.text()
        state_dict["Ri_i"] = self.tuneUI.le_Ri_i.text()
        state_dict["L_i"] = self.tuneUI.le_L_Mid_Cell.text()
        state_dict["Tune_Variable"] = self.tuneUI.le_Tune_Variable.text()

        state_dict["A_o"] = self.tuneUI.le_A_o.text()
        state_dict["B_o"] = self.tuneUI.le_B_o.text()
        state_dict["a_o"] = self.tuneUI.le_a_o.text()
        state_dict["b_o"] = self.tuneUI.le_b_o.text()
        state_dict["Ri_o"] = self.tuneUI.le_Ri_o.text()
        state_dict["Tune_Variable_End_Cell"] = self.tuneUI.le_Tune_Variable_End_Cell.text()

        # settings
        state_dict["No_Of_Processors"] = self.tuneUI.sb_No_Of_Processors_Tune.value()
        state_dict["Iterative_Method"] = self.tuneUI.cb_Iterative_Method.currentIndex()
        state_dict["Tolerance"] = self.tuneUI.le_Tolerance.text()
        state_dict["Max_Iteration"] = self.tuneUI.sb_Max_Iteration.value()

    @staticmethod
    def run_sequential(pseudo_shape_space_proc, resume, p, bc, parentDir, projectDir, filename, tuner, tune_variable, iter_set, cell_type):
        tune(pseudo_shape_space_proc, bc, parentDir, projectDir, filename, resume=resume, proc=p, tuner=tuner,
                       tune_variable=tune_variable, iter_set=iter_set, cell_type=cell_type) # last_key=last_key This would have to be tested again #val2

    def deserialize(self, state_dict):
        # update state file
        self.tuneUI.le_Freq.setText(state_dict["Frequency"])
        self.tuneUI.cb_Cell_Type.setCurrentIndex(state_dict["Cell_Type"])
        self.tuneUI.cb_Tune_Option.setCurrentIndex(state_dict["Tune_Option"])
        self.tuneUI.cb_Shape_Space_Generation_Algorithm.setCurrentIndex(state_dict["Method"])

        self.tuneUI.sb_No_Of_Shapes_Monte_Carlo.setValue(state_dict["No_Of_Shapes_Monte_Carlo"])

        self.tuneUI.cb_Tuner.setCurrentIndex(state_dict["Tuner"])
        self.tuneUI.cb_LBC.setCurrentIndex(state_dict["LBC"])
        self.tuneUI.cb_RBC.setCurrentIndex(state_dict["RBC"])

        self.tuneUI.cb_Inner_Cell.setCheckState(state_dict["Inner_Cell"])
        self.tuneUI.cb_Outer_Cell.setCheckState(state_dict["Outer_Cell"])
        self.tuneUI.cb_Expansion.setCheckState(state_dict["Expansion"])
        self.tuneUI.cb_LBP.setCheckState(state_dict["LBP"])

        # cell parameters
        self.tuneUI.le_A_i.setText(state_dict["A_i"])
        self.tuneUI.le_B_i.setText(state_dict["B_i"])
        self.tuneUI.le_a_i.setText(state_dict["a_i"])
        self.tuneUI.le_b_i.setText(state_dict["b_i"])
        self.tuneUI.le_Ri_i.setText(state_dict["Ri_i"])
        self.tuneUI.le_L_Mid_Cell.setText(state_dict["L_i"])
        self.tuneUI.le_Tune_Variable.setText(state_dict["Tune_Variable"])

        self.tuneUI.le_A_o.setText(state_dict["A_o"])
        self.tuneUI.le_B_o.setText(state_dict["B_o"])
        self.tuneUI.le_a_o.setText(state_dict["a_o"])
        self.tuneUI.le_b_o.setText(state_dict["b_o"])
        self.tuneUI.le_Ri_o.setText(state_dict["Ri_o"])
        self.tuneUI.le_Tune_Variable_End_Cell.setText(state_dict["Tune_Variable_End_Cell"])

        # settings
        self.tuneUI.sb_No_Of_Processors_Tune.setValue(state_dict["No_Of_Processors"])
        self.tuneUI.cb_Iterative_Method.setCurrentIndex(state_dict["Iterative_Method"])
        self.tuneUI.le_Tolerance.setText(state_dict["Tolerance"])
        self.tuneUI.sb_Max_Iteration.setValue(state_dict["Max_Iteration"])

    @staticmethod
    def overwriteFolder(invar, projectDir):
        path = f"{projectDir}\SimulationData\SLANS\Cavity_process_{invar}"

        if os.path.exists(path):
            shutil.rmtree(path)
            dir_util._path_created = {}

        os.makedirs(path)

    @staticmethod
    def copyFiles(invar, parentDir, projectDir):
        src = fr"{parentDir}\em_codes\SLANS_exe"
        dst = fr"{projectDir}\SimulationData\SLANS\Cavity_process_{invar}\SLANS_exe"

        dir_util.copy_tree(src, dst)

    @staticmethod
    def cb_show_hide(wid1, wid2):
        if wid1.checkState() == 2:
            wid2.show()
        else:
            wid2.hide()

    @staticmethod
    def cb_toggle(wid1, wid2, wid3):
        if wid1.checkState() == 2:
            wid2.show()
            wid3.hide()
        else:
            wid2.hide()
            wid3.show()

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
                print('here', list(l))
                return list(l)

    def continue_check(self):
        path = f'{self.main_control.projectDir}/Cavities/pseudo_{self.proof_filename(self.tuneUI.le_Generated_Shape_Space_Name.text())}'
        print(path)
        if os.path.exists(path):
            msg = QMessageBox()
            msg.setWindowTitle("Resume Simulation")
            msg.setText("The pseudo shape space file already exist. Would you love to update or overwrite its content?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)

            buttonY = msg.button(QMessageBox.Yes)
            buttonY.setText('Update')
            buttonN = msg.button(QMessageBox.No)
            buttonN.setText('Overwrite')

            msg.setDefaultButton(buttonY)

            msg.buttonClicked.connect(self.button_clicked)

            x = msg.exec_()

            if x == 16384:
                return "Yes"
            elif x == 65536:
                return "No"
            elif x == 4194304:
                return "Cancel"
        else:
            return "Cancel"

    @staticmethod
    def button_clicked(i):
        return i.text()

    @staticmethod
    def load_shape_space(filename, arg=None):
        fr = FileReader()
        dir = f'{filename}'

        # check if extension is included
        if dir.split('.')[-1] != 'json':
            dir = f'{dir}.json'

        print_(dir)
        df = fr.json_reader(dir)
        print(df)

        if arg:
            return df.loc[arg].to_dict()
        else:
            return df.to_dict(orient='list').values()

    @staticmethod
    def open_file(le, cb, d):
        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Json Files (*.json)")
        try:
            le.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            # populate checkboxes with key
            for col in dd.keys():
                cb.addItem(fr'{col}')

            d.update(dd)

        except Exception as e:
            print('Failed to open file:: ', e)

    @staticmethod
    def input_control(wid_l_bound, wid_r_bound):
        if wid_r_bound.value() < wid_l_bound.value():
            wid_r_bound.setValue(wid_l_bound.value())

    @staticmethod
    def proof_filename(dir):
        # check if extension is included
        if dir.split('.')[-1] != 'json':
            dir = f'{dir}.json'

        return dir

    def change_tune_option(self, txt):
        if txt == 'Req':
            self.tuneUI.l_Tune_Alternate_Variable.setText("L")
            self.tuneUI.l_Tune_Alternate_Variable_End_Cell.setText('L')
            self.tuneUI.le_Tune_Variable_End_Cell.setEnabled(True)
        else:
            self.tuneUI.l_Tune_Alternate_Variable.setText("Req")
            self.tuneUI.l_Tune_Alternate_Variable_End_Cell.setText("Req")

            # set end cell Req equal to mid cell Req
            self.tuneUI.le_Tune_Variable_End_Cell.setText(self.tuneUI.le_Tune_Variable.text())
            self.tuneUI.le_Tune_Variable_End_Cell.setEnabled(False)

    @staticmethod
    def process_range(l):
        if isinstance(l, list):
            if len(l) > 1:
                val = (l[0] if l[0] == l[-1] else round(r.uniform(l[0], l[-1]), 2))
                return val
            else:
                return l[0]
        elif isinstance(l, int) or isinstance(l, float):
            return l
        else:
            print("Seems something is wrong with the input.")

    def change_app_status(self):
        # Application status
        movie = QMovie(':/general/icons/GIF/spinner-icon.gif')
        movie.setScaledSize(self.ui.l_Status.size())
        self.ui.l_Status.setMovie(movie)
        self.ui.l_Status_Text.setText("Busy")
        movie.start()

        t_status = Thread(target=self.app_status, args=(movie,))
        t_status.start()

    def app_status(self, movie):
        print("Second thread started")
        # keep on while t is alive
        for t in self.thread_list:
            # while t.isAlive():
            while t.is_alive():
                pass

            t.join()

        movie.stop()
        self.ui.l_Status_Text.setText("Ready")
        self.thread_list = []

    @staticmethod
    def combine_dict(proc_count, filename, projectDir):
        # Combining dictionaries
        print_('Combining dictionaries')

        result = {}
        for index in range(proc_count):
            with open(fr'{projectDir}\Cavities\shape_space{index}.json', "r") as infile:
                result.update(json.load(infile))

        # check if extension is included
        if filename.split('.')[-1] != 'json':
            filename = f'{filename}.json'

        with open(fr'{projectDir}\Cavities\{filename}', "w") as outfile:
            json.dump(result, outfile, indent=4, separators=(',', ': '))

        print_('Done combining dictionaries')

    @staticmethod
    def delete_process_dict(proc_count, projectDir):
            for index in range(proc_count):
                os.remove(fr'{projectDir}\Cavities\shape_space{index}.json')


def tune(pseudo_shape_space, bc, parentDir, projectDir, filename, resume="No",
          proc=0, tuner='SLANS', tune_variable='Req', iter_set=None, cell_type='Mid Cell'):

    # tuner
    pytune = Tuner()

    if tuner == 'SLANS':
        slans_tune = Tune(parentDir, projectDir)
    else:
        slans_tune = None

    start = time.time()
    population = {}
    dip_dist = {}
    freq = {}
    total_no_of_shapes = len(list(pseudo_shape_space.keys()))

    # check for already processed shapes
    existing_keys = []

    if resume == "Yes":
        # check if value set is already written. This is to enable continuation in case of break in program
        if os.path.exists(fr'{projectDir}/Cavities/{filename}'):
            population = json.load(open(fr'{projectDir}/Cavities/{filename}', 'r'))

            existing_keys = list(population.keys())
            print_(f'Existing keys: {existing_keys}')

    start_time = time.time()

    progress = 0
    # write to progress file
    try:
        with open(fr"{projectDir}\SimulationData\SLANS\Cavity_process_{proc}\progress_file.txt", 'w') as file:
            txt = f"{progress+1}/{total_no_of_shapes}"
            file.write(txt)
    except:
        pass

    for key, pseudo_shape in pseudo_shape_space.items():

        A_i, B_i, a_i, b_i, Ri_i, L_i, Req, _ = pseudo_shape['IC'] # Req here is none but required since the shape space consists of all variables
        A_o, B_o, a_o, b_o, Ri_o, L_o, Req, _ = pseudo_shape['OC'] # Req here is none but required since the shape space consists of all variables
        beampipes = pseudo_shape['BP']
        freq = pseudo_shape['FREQ']

        # new mid cell and end cell with initial Req guess
        inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0]
        outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0]

        # edit to check for key later
        if key not in existing_keys:
            if tuner == 'SLANS' and slans_tune:
                if tune_variable == 'Req':
                    # Tune cell to get Req
                    Req, freq, alpha, h, e = slans_tune.mid_cell_tune(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, freq, proc=proc)
                else:
                    L, freq, alpha, h, e = slans_tune.end_cell_tune(inner_cell, outer_cell, freq, proc=proc)

                    if cell_type == 'Mid Cell':
                        L_i, L_o = L, L
                    else:
                        L_o = L

                inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha]
                outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha]
            else:
                if tune_variable == 'Req':
                    print('PyTune Req')
                    Req, freq = pytune.tune("Req", inner_cell, outer_cell, freq, beampipes, bc, parentDir, projectDir, iter_set=iter_set, proc=proc)
                    # round
                    Req, freq = round(Req, 4), round(freq, 2)
                else:
                    print('PyTune L')
                    L, freq = pytune.tune("L", inner_cell, outer_cell, freq, beampipes, bc, parentDir, projectDir, iter_set=iter_set, proc=proc)

                    # round
                    L, freq = round(L, 2), round(freq, 2)

                    if cell_type == 'Mid Cell':
                        L_i, L_o = L, L
                    else:
                        L_o = L

                alpha_i = calculate_alpha(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0)
                alpha_o = calculate_alpha(A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0)

                inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha_i]
                outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha_o]

            print_(f'Done Tuning Cavity {key}')

            # write to progress file
            try:
                with open(fr"{projectDir}\SimulationData\SLANS\Cavity_process_{proc}\progress_file.txt", 'w') as file:
                    txt = f"{progress+1}/{total_no_of_shapes}"
                    file.write(txt)
            except:
                pass

            if cell_type == 'Mid Cell':
                population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'none', 'FREQ': freq}
            else:
                population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'both', 'FREQ': freq}

        # Update progressbar
        progress += 1

        print_("SAVING DICTIONARIES", f"shape_space{proc}.json")
        with open(fr"{projectDir}\Cavities\shape_space{proc}.json", 'w') as file:
            file.write(json.dumps(population, indent=4, separators=(',', ': ')))
        print_("DONE SAVING")

        print_("Time run:: ", time.time() - start_time)
        start_time = time.time()

    # print_("Extracted Data::", self.population)
    end = time.time()

    runtime = end - start
    print_(f'Proc {proc} runtime: {runtime}')


def calculate_alpha(A, B, a, b, Ri, L, Req, L_bp):

    data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
            [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1, y1, x2, y2 = fsolve(ellipse_tangent,
                            np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                            args=data)
    m = (y2 - y1) / (x2 - x1)
    alpha = 180 - np.arctan(m) * 180 / np.pi
    return alpha


def ellipse_tangent(z, *data):
    coord, dim = data
    h, k, p, q = coord
    a, b, A, B = dim
    x1, y1, x2, y2 = z

    f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
    f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
    f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
    f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

    return f1, f2, f3, f4
