import math
from datetime import datetime
import random
import shutil
import sys
import scipy.signal as sps
from scipy.stats import qmc
from math import floor
import threading
import numpy as np
import oapackage
import pandas as pd
from PyQt5 import QtCore
from PyQt5.QtCore import QPropertyAnimation, QThread, pyqtSignal
from PyQt5.QtGui import QPixmap, QIcon, QColor
from icecream import ic
from modules.data_module.abci_data import ABCIData
# from modules.data_module.slans_data import SLANSDataExtraction
from simulation_codes.ABCI.abci_geometry import ABCIGeometry
from ui_files.run_tune import Ui_w_Tune
import json
import os
import time
import multiprocessing as mp
from termcolor import colored
from PyQt5.QtWidgets import *
import random as r
from utils.file_reader import FileReader
import pyqtgraph as pg
from distutils import dir_util
import psutil
from modules.tune_module.tuners.tuner import Tuner
from utils.shared_classes import *
from utils.shared_functions import *

from simulation_codes.SLANS.slans_geom_par import SLANSGeometry
slans_geom = SLANSGeometry()  # parallel implementaion of code, consider making the two separate classes 1
from simulation_codes.SLANS.slans_geometry import SLANSGeometry
slans_geom_seq = SLANSGeometry()

file_color = 'green'
DEBUG = True


def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))


fr = FileReader()
abci_geom = ABCIGeometry()
# slans_data_extraction = SLANSDataExtraction()

tuner = Tuner()

AN_DURATION = 200


class TuneControl:
    def __init__(self, parent):
        self.pseudo_shape_space = None
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

        # ui effects
        self.ui_effects()

        # evolution algorithm
        self.opt_control = OptimizationControl(self.tuneUI)

    def initUI(self):
        self.tuneUI.cb_Inner_Cell.setCheckState(2)
        self.tuneUI.cb_Inner_Cell.setEnabled(False)

        # expand/collapse sections widgets
        if self.tuneUI.cb_Expansion.checkState() == 2:
            self.tuneUI.w_Expansion.setMinimumHeight(150)
        else:
            self.tuneUI.w_Expansion.setMinimumHeight(0)
            self.tuneUI.w_Expansion.setMaximumHeight(0)

        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
            self.tuneUI.w_Outer_Cell.setMinimumHeight(150)
        else:
            self.tuneUI.w_Outer_Cell.setMinimumHeight(0)
            self.tuneUI.w_Outer_Cell.setMaximumHeight(0)

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
        # self.animate_width(self.tuneUI.cb_Monitor_Convergence, self.tuneUI.w_PyqtGraph, 0, 500, True)

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

        # create progress bar object and add to widget
        self.progress_bar = QProgressBar(self.tuneUI.w_Simulation_Controls)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.tuneUI.gl_Simulation_Controls.addWidget(self.progress_bar, 0, 7, 1, 1)
        self.progress_bar.hide()

        # disable convergence monitor checkbox
        self.tuneUI.cb_Monitor_Convergence.setEnabled(False)

    def signals(self):

        # signals/slots
        self.tuneUI.cb_Shape_Space_Generation_Algorithm.currentIndexChanged.connect(lambda: self.show_hide())
        self.tuneUI.pb_Run.clicked.connect(lambda: self.generate_shape_space())

        # tune variable
        # self.tuneUI.cb_Tune_Option.currentTextChanged.connect(lambda: self.change_tune_option(self.tuneUI.cb_Tune_Option.currentText()))

        # cell parameters control
        self.tuneUI.cb_Outer_Cell.stateChanged.connect(
            lambda: self.animate_height(self.tuneUI.cb_Outer_Cell, self.tuneUI.w_Outer_Cell, 0, 160, True))
        self.tuneUI.cb_Expansion.stateChanged.connect(
            lambda: self.animate_height(self.tuneUI.cb_Expansion, self.tuneUI.w_Expansion, 0, 160, True))
        self.tuneUI.cb_Tuner.currentTextChanged.connect(lambda: self.tuner_routine())
        self.tuneUI.cb_Cell_Type.currentTextChanged.connect(lambda: self.tuner_routine())
        self.tuneUI.cb_Tune_Option.currentTextChanged.connect(lambda: self.tuner_routine())

        # change cavity display image
        self.tuneUI.cb_Cell_Type.currentTextChanged.connect(lambda: self.change_cavity_image())
        self.tuneUI.cb_LBP.stateChanged.connect(lambda: self.change_cavity_image())
        self.tuneUI.cb_RBP.stateChanged.connect(lambda: self.change_cavity_image())
        self.tuneUI.cb_Outer_Cell.stateChanged.connect(lambda: self.change_cavity_image())

        # collapse other settings
        self.tuneUI.pb_Tune_Settings.clicked.connect(
            lambda: self.main_control.animate_height(self.tuneUI.w_Tune_Settings, 0, 200, True))

        # disable secondary settings if Tuner is SLANS
        self.tuneUI.cb_Tuner.currentTextChanged.connect(
            lambda: self.tuneUI.w_Iteration_Settings.setEnabled(False) if self.tuneUI.cb_Tuner.currentText() == 'SLANS'
            else self.tuneUI.w_Iteration_Settings.setEnabled(True))

        # control to ensure that SLANS mid cell tuner is always set to tune for Req and end cell tuner is always set to L
        self.tuneUI.cb_Tuner.currentTextChanged.connect(lambda: self.slans_tuners_control())
        self.tuneUI.cb_Cell_Type.currentTextChanged.connect(lambda: self.slans_tuners_control())

        # expand/collapse convergence plot
        self.tuneUI.cb_Monitor_Convergence.stateChanged.connect(
            lambda: self.animate_width(self.tuneUI.cb_Monitor_Convergence, self.tuneUI.w_PyqtGraph, 0, 500, True))

        # cancel
        self.tuneUI.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.tuneUI.pb_Pause_Resume.clicked.connect(
            lambda: self.pause() if self.process_state == 'running' else self.resume())

        # change variable value
        # self.tuneUI.le_Tune_Variable.textChanged.connect(lambda: self.tuneUI.le_Tune_Variable_End_Cell.setText(self.tuneUI.le_Tune_Variable.text())
        #                                                  if self.tuneUI.cb_Tune_Option.currentIndex() == 1
        #                                                  else None)

    def tuner_routine(self):
        if self.tuneUI.cb_Cell_Type.currentText() == 'Mid Cell':
            self.tuneUI.cb_Outer_Cell.setCheckState(0)
            self.tuneUI.cb_Outer_Cell.setEnabled(False)
            self.tuneUI.cb_LBP.setCheckState(0)
            self.tuneUI.cb_LBP.setEnabled(False)
            # self.tuneUI.w_End_Cell_Tune_Extra_Variable.hide()
        elif self.tuneUI.cb_Cell_Type.currentText() == 'End Cell':
            self.tuneUI.cb_Outer_Cell.setCheckState(2)
            self.tuneUI.cb_Outer_Cell.setEnabled(False)
            self.tuneUI.cb_LBP.setCheckState(2)
            self.tuneUI.cb_LBP.setEnabled(False)

            # if self.tuneUI.cb_Tune_Option.currentText() == 'L':
            #     self.tuneUI.w_End_Cell_Tune_Extra_Variable.show()
            # else:
            #     self.tuneUI.w_End_Cell_Tune_Extra_Variable.hide()
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
                    self.pseudo_shape_space = json.load(
                        open(fr'{self.main_control.projectDir}\Cavities\pseudo_{self.filename}', 'r'))

                    self.existing_keys = list(self.pseudo_shape_space.keys())
                    # self.log(f'last saved key: {self.existing_keys}')
            elif resume == "Cancel":
                return

        self.freq = float(self.tuneUI.le_Freq.text())
        # get variables from ui or from pseudo shape space
        A_i = self.check_input(self.tuneUI.le_A_i.text())
        B_i = self.check_input(self.tuneUI.le_B_i.text())
        a_i = self.check_input(self.tuneUI.le_a_i.text())
        b_i = self.check_input(self.tuneUI.le_b_i.text())
        Ri_i = self.check_input(self.tuneUI.le_Ri_i.text())
        L_i = self.check_input(self.tuneUI.le_L_i.text())
        Req_i = self.check_input(self.tuneUI.le_Req_i.text())

        inner_half_cell_parameters = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i]
        # print(inner_half_cell_parameters)

        # ihc = self.create_pseudo_shape_space(inner_half_cell_parameters, lock_list)

        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
            A_o = self.check_input(self.tuneUI.le_A_o.text())
            B_o = self.check_input(self.tuneUI.le_B_o.text())
            a_o = self.check_input(self.tuneUI.le_a_o.text())
            b_o = self.check_input(self.tuneUI.le_b_o.text())
            Ri_o = self.check_input(self.tuneUI.le_Ri_o.text())
            L_o = self.check_input(self.tuneUI.le_L_o.text())
            Req_o = Req_i

            outer_half_cell_parameters = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o]
            # ihc = self.create_pseudo_shape_space(outer_half_cell_parameters, lock_list)
        else:
            outer_half_cell_parameters = inner_half_cell_parameters

        # print(outer_half_cell_parameters)
        lock_list = [False, False, False, False, False, False, False]

        ihc = self.create_pseudo_shape_space(inner_half_cell_parameters, lock_list, "Mid Cell")
        ohc = self.create_pseudo_shape_space(outer_half_cell_parameters, lock_list, "End Cell")

        pseudo_shape_space = self.generate_pseudo_shape_space(self.freq, ihc, ohc)

        if pseudo_shape_space:
            self.run_tune(pseudo_shape_space, resume)

    def generate_shape_space_old(self):
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
                    self.pseudo_shape_space = json.load(
                        open(fr'{self.main_control.projectDir}\Cavities\pseudo_{self.filename}', 'r'))

                    self.existing_keys = list(self.pseudo_shape_space.keys())
                    # self.log(f'last saved key: {self.existing_keys}')
            elif resume == "Cancel":
                return

        self.freq = float(self.tuneUI.le_Freq.text())
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
            Ri_o = self.text_to_list(self.tuneUI.le_Ri_o.text())

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

        pseudo_shape_space = self.generate_pseudo_shape_space(self.freq, inner_half_cell_parameters,
                                                              outer_half_cell_parameters)

        if pseudo_shape_space:
            self.run_tune(pseudo_shape_space, resume)

    def run_tune(self, pseudo_shape_space, resume):
        # set tuner
        self.tuner_option = self.tuneUI.cb_Tuner.currentText()
        proc_count = self.tuneUI.sb_No_Of_Processors_Tune.value()
        tune_variable = self.tuneUI.cb_Tune_Option.currentText().split(' ')[-1]

        # get iteration settings
        iter_method = self.tuneUI.cb_Iterative_Method.currentText()
        tolerance = self.check_input(self.tuneUI.le_Tolerance.text())[0]
        max_iter = self.check_input(self.tuneUI.sb_Max_Iteration.text())[0]
        iter_set = [iter_method, tolerance, max_iter]

        # cell type
        n_cells = self.tuneUI.sb_Tune_N_Cells.value()
        cell_type = self.tuneUI.cb_Cell_Type.currentText()

        # boundary conditions
        lbc = self.tuneUI.cb_LBC.currentIndex() + 1
        rbc = self.tuneUI.cb_RBC.currentIndex() + 1
        bc = 10 * lbc + rbc

        # save last
        save_last = self.tuneUI.cb_Save_Last.isChecked()

        # try:
        if True:
            self.start = time.time()

            # split shape_space for different processes/ MPI share process by rank
            keys = list(pseudo_shape_space.keys())

            # check if number of processors selected is greater than the number of keys in the pseudo shape space
            if proc_count > len(keys):
                proc_count = len(keys)

            shape_space_len = len(keys)
            share = floor(shape_space_len / proc_count)

            # for i in reversed(range(self.tuneUI.gl_PyqtGraph.count())):
            #     self.tuneUI.gl_PyqtGraph.itemAt(i).widget().setParent(None)
            #
            # # insert graphs for convergence monitor
            # self.pg_list = []
            # self.pygraph_list = []
            # graph_order = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (3, 0), (3, 1), (4, 0), (4, 1), (5, 0), (5, 1), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (6, 0), (6, 1), (6, 2), (7, 0), (7, 1), (7, 2), (8, 0), (8, 1), (8, 2), (0, 3), (1, 3), (2, 3), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3)]
            # # print("pglist before ", self.pg_list)
            # # print("pygraph before ", self.pygraph_list)
            #
            # for p in range(proc_count):
            #     self.createPyGraph(graph_order[p][0], graph_order[p][1])

            # print("pglist after ", self.pg_list)
            # print("pygraph after ", self.pygraph_list)

            # show progress bar
            self.progress_bar.show()
            time.sleep(1)

            # progress list
            manager = mp.Manager()
            self.progress_list = manager.list()
            self.progress_list.append(0)

            # convergence list
            self.convergence_list = manager.list()
            for i in range(proc_count):
                self.convergence_list.append([])

            for p in range(proc_count):
                if True:
                    if p < proc_count - 1:
                        proc_keys_list = keys[p * share:p * share + share]
                    else:
                        proc_keys_list = keys[p * share:]

                    self.overwriteFolder(p, self.main_control.projectDir)
                    self.copyFiles(p, self.main_control.parentDir, self.main_control.projectDir)

                    processor_shape_space = {}
                    for key, val in pseudo_shape_space.items():
                        if key in proc_keys_list:
                            processor_shape_space[key] = val

                    service = mp.Process(target=self.run_sequential,
                                         args=(processor_shape_space, resume, p, bc, self.main_control.parentDir,
                                               self.main_control.projectDir, self.filename, self.tuner_option,
                                               tune_variable, iter_set, cell_type, self.progress_list,
                                               self.convergence_list, save_last, n_cells))
                    service.start()
                    self.processes.append(psutil.Process(service.pid))
                    self.processes_id.append(service.pid)

                # except Exception as e:
                #     self.log.info("Exception in run_MP::", e)
                #     # write to log

            # display progress bar
            self.show_progress_bar = True
            self.progress_monitor_thread = ProgressMonitor(self, self.main_control.projectDir)
            self.progress_monitor_thread.sig.connect(self.update_progress_bar)
            self.progress_monitor_thread.start()

            # # start convergence monitor
            # self.convergence_monitor_thread = MonitorConvergence(self)
            # self.convergence_monitor_thread.sig.connect(self.monitor_convergence)
            # self.convergence_monitor_thread.start()

            self.log.info("Tune started")
            # change process state to running
            self.process_state = 'running'
            self.run_pause_resume_stop_routine()

            self.end_routine_thread = EndRoutine(self, self.main_control.projectDir)
            self.end_routine_thread.start()

        # except Exception as e:
        #     self.log.error(fr"TUNE CONTROL:: run_tune:: {e}")
        #     print_("TUNE CONTROL:: run_tune:: ", e)

    def createPyGraph(self, row, column):
        # add pyqgraph
        pg_ = pg
        pygraph = pg_.PlotWidget()
        pygraph.setBackground('w')
        self.tuneUI.gl_PyqtGraph.addWidget(pygraph, row, column, 1, 1)

        self.pg_list.append(pg_)
        self.pygraph_list.append(pygraph)

    def monitor_convergence(self, conv):
        l = len(self.convergence_list)
        self.plot_list = [None for i in range(l)]

        for p in self.pygraph_list:
            if p:
                p.clear()

        if self.tuneUI.cb_Monitor_Convergence.checkState() == 2:
            for i in range(l):
                self.pygraph_list[i].addLine(x=None, y=self.freq, pen=pg.mkPen('r', width=1))

                if not self.plot_list[i]:
                    self.plot_list[i] = self.pygraph_list[i].plot(self.convergence_list[i][0],
                                                                  self.convergence_list[i][1], pen=None, symbol='o')
                else:
                    self.plot_list[i].setData(self.convergence_list[i][0], self.convergence_list[i][1], pen=None,
                                              symbol='o')

    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

        if val == 100 or not self.show_progress_bar:
            # reset progress bar
            self.progress_bar.setValue(0)
            self.progress_bar.hide()

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

        # lear processes
        self.processes.clear()
        self.processes_id.clear()

        # reset application state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        self.log.info("Process terminated.")

        # set exit signal
        self.tune_ended = True
        print("Tune ended")

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
            combined_dict = self.combine_dict(proc_count, filename, projectDir)
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
                    if os.path.exists(
                            fr'{self.main_control.projectDir}\SimulationData\SLANS\Cavity_process_{i}\progress_file.txt'):
                        with open(
                                fr'{self.main_control.projectDir}\SimulationData\SLANS\Cavity_process_{i}\progress_file.txt',
                                "r") as f:
                            a = f.readline()
                            progress += eval(a)
                self.progress_bar.setValue(progress * 100 / len(proc_ids))
            except:
                print("Error in progress update")
                pass

        try:
            # reset progress bar
            self.progress_bar.setValue(0)
            self.progress_bar.hide()
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

        BP = 'none'
        if self.tuneUI.cb_LBP.checkState() == 2:
            BP = 'left'

        if self.tuneUI.cb_LBP.checkState() == 2 and self.tuneUI.cb_RBP.checkState() == 2:
            BP = 'both'

        if self.tuneUI.cb_Inner_Cell.isChecked() and self.tuneUI.cb_Outer_Cell.isChecked():
            print("IO")
            key = 0
            for indx1, inner_cell in ihc.iterrows():
                inner_cell = inner_cell.tolist()
                for indx2, other_cell in ohc.iterrows():
                    # enforce Req_inner_cell == Req_outer_cell
                    other_cell = other_cell.tolist()
                    other_cell[-2] = inner_cell[-2]

                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': other_cell, 'BP': BP, 'FREQ': freq}
                    key += 1

        elif self.tuneUI.cb_Inner_Cell.isChecked() and not self.tuneUI.cb_Outer_Cell.isChecked():
            print("II")
            key = 0
            for indx, inner_cell in ihc.iterrows():
                self.pseudo_shape_space[key] = {'IC': inner_cell.tolist(), 'OC': inner_cell.tolist(), 'BP': BP,
                                                'FREQ': freq}
                key += 1

        # remove duplicates from generated space
        self.pseudo_shape_space = self.remove_duplicate_values(self.pseudo_shape_space)

        pseudo_shape_space_name = f'{self.main_control.projectDir}/Cavities/pseudo_{self.proof_filename(self.tuneUI.le_Generated_Shape_Space_Name.text())}'
        with open(pseudo_shape_space_name, 'w') as file:
            file.write(json.dumps(self.pseudo_shape_space, indent=4, separators=(',', ': ')))

        return self.pseudo_shape_space

    def create_pseudo_shape_space(self, var_list, lock_list, cell):
        A, B, a, b, Ri, L, Req = var_list
        A_LOCKED, B_LOCKED, a_LOCKED, b_LOCKED, Ri_LOCKED, L_LOCKED, Req_LOCKED = lock_list

        if self.tuneUI.cb_Tune_Option.currentText() == "Req":
            if Req[0] == 0:
                AA, BB, aa, bb, RiRi, LL = np.meshgrid(A, B, a, b, Ri, L)

                space = []
                for j in range(len(A)):  # for some reason the first dimension of A is the dimension of B
                    for i in range(len(B)):
                        for k in range(len(a)):
                            for m in range(len(b)):
                                for n in range(len(Ri)):
                                    for o in range(len(L)):
                                        # enforce A + a <= L
                                        if AA[i][j][k][m][n][o] + aa[i][j][k][m][n][o] <= LL[i][j][k][m][n][o]:
                                            space.append(
                                                [AA[i][j][k][m][n][o], BB[i][j][k][m][n][o], aa[i][j][k][m][n][o],
                                                 bb[i][j][k][m][n][o], RiRi[i][j][k][m][n][o],
                                                 LL[i][j][k][m][n][o],
                                                 RiRi[i][j][k][m][n][o] + BB[i][j][k][m][n][o] + bb[i][j][k][m][n][
                                                     o],
                                                 0,0])
            else:
                AA, BB, aa, bb, RiRi, LL, ReqReq = np.meshgrid(A, B, a, b, Ri, L, Req)
                space = []
                for j in range(len(A)):  # for some reason the first dimension of A is the dimension of B
                    for i in range(len(B)):
                        for k in range(len(a)):
                            for m in range(len(b)):
                                for n in range(len(Ri)):
                                    for o in range(len(L)):
                                        for p in range(len(Req)):
                                            # enforce A + a < L
                                            if AA[i][j][k][m][n][o][p] + aa[i][j][k][m][n][o][p] <= \
                                                    LL[i][j][k][m][n][o][p]:
                                                space.append([AA[i][j][k][m][n][o][p], BB[i][j][k][m][n][o][p],
                                                              aa[i][j][k][m][n][o][p],
                                                              bb[i][j][k][m][n][o][p], RiRi[i][j][k][m][n][o][p],
                                                              LL[i][j][k][m][n][o][p],
                                                              ReqReq[i][j][k][m][n][o][p], 0, 0])

        elif self.tuneUI.cb_Tune_Option.currentText() == "L":
            if L[0] == 0:
                AA, BB, aa, bb, RiRi, ReqReq = np.meshgrid(A, B, a, b, Ri, Req)

                space = []
                for j in range(len(A)):  # for some reason the first dimension of A is the dimension of B
                    for i in range(len(B)):
                        for k in range(len(a)):
                            for m in range(len(b)):
                                for n in range(len(Ri)):
                                    for o in range(len(L)):
                                        # enforce Ri + B + b <= Req
                                        if RiRi[i][j][k][m][n][o] + BB[i][j][k][m][n][o] + bb[i][j][k][m][n][o] <= \
                                                ReqReq[i][j][k][m][n][o]:
                                            space.append(
                                                [AA[i][j][k][m][n][o], BB[i][j][k][m][n][o], aa[i][j][k][m][n][o],
                                                 bb[i][j][k][m][n][o], RiRi[i][j][k][m][n][o],
                                                 AA[i][j][k][m][n][o] + aa[i][j][k][m][n][o],
                                                 ReqReq[i][j][k][m][n][o], 0, 0])
            else:
                AA, BB, aa, bb, RiRi, LL, ReqReq = np.meshgrid(A, B, a, b, Ri, L, Req)

                space = []
                for j in range(len(A)):  # for some reason the first dimension of A is the dimension of B
                    for i in range(len(B)):
                        for k in range(len(a)):
                            for m in range(len(b)):
                                for n in range(len(Ri)):
                                    for o in range(len(L)):
                                        for p in range(len(Req)):
                                            # enforce Ri + B + b <= Req
                                            if RiRi[i][j][k][m][n][o][p] + BB[i][j][k][m][n][o][p] + \
                                                    bb[i][j][k][m][n][o][p] <= ReqReq[i][j][k][m][n][o][p]:
                                                space.append([AA[i][j][k][m][n][o][p], BB[i][j][k][m][n][o][p],
                                                              aa[i][j][k][m][n][o][p], bb[i][j][k][m][n][o][p],
                                                              RiRi[i][j][k][m][n][o][p],
                                                              LL[i][j][k][m][n][o][p],
                                                              ReqReq[i][j][k][m][n][o][p], 0, 0])

        else:
            print("You have to set the field value of the variable to be tuned to zero,")
            return 1

        space = np.array(space)
        # print(space)
        # create list of locked variables
        LOCKED_LIST = np.array([A_LOCKED, B_LOCKED, a_LOCKED, b_LOCKED, Ri_LOCKED, L_LOCKED])
        count = np.count_nonzero(LOCKED_LIST)
        print(count)
        if count >= 2:
            # for ll in LOCKED_LIST:
            # check if length of all locked lists are same
            dummy_list = []
            for i in LOCKED_LIST.nonzero()[0]:
                dummy_list.append(len(var_list[i]))
                lock_len = len(var_list[i])
            dummy_list = np.array(dummy_list)

            if np.all(dummy_list == dummy_list[0]):
                print("Perfect")
            else:
                print("Please make sure locked variables are of equal lengths.")
                return

            if not A_LOCKED:
                A = np.ones(lock_len) * (-1)
            if not B_LOCKED:
                B = np.ones(lock_len) * (-1)
            if not a_LOCKED:
                a = np.ones(lock_len) * (-1)
            if not b_LOCKED:
                b = np.ones(lock_len) * (-1)
            if not Ri_LOCKED:
                Ri = np.ones(lock_len) * (-1)
            if not L_LOCKED:
                L = np.ones(lock_len) * (-1)
            if not Req_LOCKED:
                L = np.ones(lock_len) * (-1)

            lock = list(zip(A, B, a, b, Ri, L))
            print(lock)

            slice = []
            for s in space:
                ll = []
                for z in lock:
                    ll = [i for (i, j) in zip(s, z) if i == j]

                    if len(ll) == count:
                        slice.append(s)

            df = pd.DataFrame(slice, columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha_i", "alpha_o"])
            # print(df)

        df = pd.DataFrame(space, columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha_i", "alpha_o"])
        # print(df)

        return df

    def text_to_list(self, s):
        # s = "range(16, 23, 10)"
        # s = "randrange(16, 23, 10)"
        # s = "[16, 23, 10]"

        if "range" in s and "rand" not in s:
            s = s.replace('range', '')
            try:
                l = ast.literal_eval(s)
                return np.linspace(l[0], l[1], l[2])
            except:
                print("Please check inputs.")
        elif "randrange" in s:
            s = s.replace('randrange', '')
            try:
                l = ast.literal_eval(s)
                ll = np.random.uniform(l[0], l[1], l[2])
                return ll
            except:
                print("Please check inputs.")
        else:
            try:
                ll = ast.literal_eval(s)
                return ll
            except:
                print("Please check inputs.")

        return 1

    def generate_pseudo_shape_space_old(self, freq, ihc, ohc):
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
                    if ihc[1] == "A":
                        B_i = A_i
                    else:
                        B_i = self.process_range(ihc[1])
                    a_i = self.process_range(ihc[2])
                    if ihc[3] == "a":
                        b_i = a_i
                    else:
                        b_i = self.process_range(ihc[3])
                    Ri_i = self.process_range(ihc[4])

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

                    if self.tuneUI.cb_Outer_Cell.checkState() == 2:
                        # This also checks if the right and left bounds are equal in which case it returns a single value
                        A_o = self.process_range(ohc[0])
                        if ohc[1] == "A":
                            B_o = A_o
                        else:
                            B_o = self.process_range(ohc[1])
                        a_o = self.process_range(ohc[2])

                        if ohc[3] == "a":
                            b_o = a_o
                        else:
                            b_o = self.process_range(ohc[3])

                        Ri_o = self.process_range(ohc[4])

                        if self.tuneUI.cb_Tune_Option.currentText() == 'L':
                            L_o = A_o + a_o
                        else:
                            L_o = self.process_range(ohc[5])

                        Req_o = Req_i

                        # update L_i
                        # L_i = self.process_range(self.text_to_list(self.tuneUI.le_Tune_Variable_End_Cell.text()))
                        # inner_cell[5] = L_i

                        other_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, 0]
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
                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': other_cell, 'BP': 'left',
                                                            'FREQ': freq}
                    else:
                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                        if key not in self.existing_keys:
                            self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': other_cell, 'BP': 'none',
                                                            'FREQ': freq}
                    count += 1
                    loop_escape = 0  # reset loop escape
            else:
                print('Please check input parameters.')
        else:
            ########## UPDATE CODE TO CHECK FOR ALL THE ITEMS IN THE SHAPE SPACE FOR SHAPES WITH A + a > L
            check = self.check_input()

            A_o_space = ohc[0]
            B_o_space = ohc[1]
            a_o_space = ohc[2]
            b_o_space = ohc[3]
            Ri_o_space = ohc[4]

            count = 0
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
                                            Ri_i_space = ihc[4]
                                            L_i_space = ihc[5]
                                            Req_i = Req_o

                                            for A_i in A_i_space:
                                                for B_i in B_i_space:
                                                    for a_i in a_i_space:
                                                        for b_i in b_i_space:
                                                            for Ri_i in Ri_i_space:
                                                                for L_i in L_i_space:
                                                                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i,
                                                                                  0]
                                                                    if self.tuneUI.cb_LBP.checkState() == 2:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {
                                                                                'IC': inner_cell, 'OC': outer_cell,
                                                                                'BP': 'left', 'FREQ': freq}
                                                                    else:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {
                                                                                'IC': inner_cell, 'OC': outer_cell,
                                                                                'BP': 'none', 'FREQ': freq}
                                                                    count += 1
                                        else:
                                            inner_cell = outer_cell

                                            if self.tuneUI.cb_LBP.checkState() == 2:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell,
                                                                                    'BP': 'left', 'FREQ': freq}
                                            else:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell,
                                                                                    'BP': 'none', 'FREQ': freq}
                                            count += 1

                                else:
                                    Req_o_space = ohc[6]
                                    for Req_o in Req_o_space:
                                        L_o = A_o + a_o
                                        outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, 0]

                                        if self.tuneUI.cb_Outer_Cell.checkState() == 2:
                                            A_i_space = ihc[0]
                                            B_i_space = ihc[1]
                                            a_i_space = ihc[2]
                                            b_i_space = ihc[3]
                                            Ri_i_space = ihc[4]
                                            L_i_space = ihc[5]
                                            Req_i = Req_o

                                            for A_i in A_i_space:
                                                for B_i in B_i_space:
                                                    for a_i in a_i_space:
                                                        for b_i in b_i_space:
                                                            for Ri_i in Ri_i_space:
                                                                for L_i in L_i_space:
                                                                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i,
                                                                                  0]

                                                                    if self.tuneUI.cb_LBP.checkState() == 2:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {
                                                                                'IC': inner_cell, 'OC': outer_cell,
                                                                                'BP': 'left', 'FREQ': freq}
                                                                    else:
                                                                        key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                                        if key not in self.existing_keys:
                                                                            self.pseudo_shape_space[key] = {
                                                                                'IC': inner_cell, 'OC': outer_cell,
                                                                                'BP': 'none', 'FREQ': freq}
                                                                    count += 1
                                        else:
                                            inner_cell = outer_cell

                                            if self.tuneUI.cb_LBP.checkState() == 2:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell,
                                                                                    'BP': 'left', 'FREQ': freq}
                                            else:
                                                key = f"{self.tuneUI.le_Marker.text()}_{count}"
                                                if key not in self.existing_keys:
                                                    self.pseudo_shape_space[key] = {'IC': inner_cell, 'OC': outer_cell,
                                                                                    'BP': 'none', 'FREQ': freq}
                                            count += 1

        # remove duplicates from generated space
        self.pseudo_shape_space = self.remove_duplicate_values(self.pseudo_shape_space)

        if check:
            pseudo_shape_space_name = f'{self.main_control.projectDir}/Cavities/pseudo_{self.proof_filename(self.tuneUI.le_Generated_Shape_Space_Name.text())}'
            with open(pseudo_shape_space_name, 'w') as file:
                file.write(json.dumps(self.pseudo_shape_space, indent=4, separators=(',', ': ')))

            return self.pseudo_shape_space
        else:
            return check

    def remove_duplicate_values(self, d):
        temp = []
        res = dict()
        for key, val in d.items():
            if val not in temp:
                temp.append(val)
                res[key] = val
        return res

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
            # GET WIDTH
            width = widget.width()
            # SET MAX WIDTH
            if cb.checkState() != 2:
                widthCollapsed = min_width
                widget.setMinimumWidth(0)
            else:
                widthCollapsed = standard
                # widget.setMinimumWidth(standard)

            # ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumWidth")
            self.animation.setDuration(AN_DURATION)
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
            elif self.tuneUI.cb_LBP.checkState() == 2 and self.tuneUI.cb_RBP.checkState() == 2:
                self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/one_cell.png"))
            elif self.tuneUI.cb_LBP.checkState() == 2 and self.tuneUI.cb_Outer_Cell.checkState() == 0:
                self.tuneUI.l_Cavity_Image.setPixmap(QPixmap(f":/imgs/images/end_same_cell.png"))

    def check_input(self, s):
        # s = "range(16, 23, 10)"
        # s = "randrange(16, 23, 10)"
        # s = "[16, 23, 10]"
        # s = 1, 2, 3
        # s = 2

        if "r" in s and "rr" not in s:
            s = s.replace('r', '')
            try:
                l = eval(s)
                return np.linspace(l[0], l[1], l[2])
            except:
                print("Please check inputs.")
        elif "rr" in s:
            s = s.replace('rr', '')
            try:
                l = eval(s)
                ll = np.random.uniform(l[0], l[1], l[2])
                return ll
            except:
                print("Please check inputs.")
        else:
            try:
                ll = eval(s)
                if isinstance(ll, int) or isinstance(ll, float):
                    ll = [ll]
                return ll
            except:
                print("Please check inputs.")

        return 1

    def check_input_old(self):
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
                # self.tuneUI.l_Tune_Alternate_Variable.setText('L')
                self.tuneUI.cb_Tune_Option.setEnabled(False)
            else:
                self.tuneUI.cb_Tune_Option.setCurrentIndex(1)
                # self.tuneUI.l_Tune_Alternate_Variable.setText('Req')
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
        state_dict["L_i"] = self.tuneUI.le_L_i.text()
        state_dict["Req_i"] = self.tuneUI.le_Req_i.text()

        state_dict["A_o"] = self.tuneUI.le_A_o.text()
        state_dict["B_o"] = self.tuneUI.le_B_o.text()
        state_dict["a_o"] = self.tuneUI.le_a_o.text()
        state_dict["b_o"] = self.tuneUI.le_b_o.text()
        state_dict["Ri_o"] = self.tuneUI.le_Ri_o.text()
        state_dict["L_o"] = self.tuneUI.le_L_o.text()

        # settings
        state_dict["No_Of_Processors"] = self.tuneUI.sb_No_Of_Processors_Tune.value()
        state_dict["Iterative_Method"] = self.tuneUI.cb_Iterative_Method.currentIndex()
        state_dict["Tolerance"] = self.tuneUI.le_Tolerance.text()
        state_dict["Max_Iteration"] = self.tuneUI.sb_Max_Iteration.value()

        ## Optimization control
        state_dict["Optimization_Algorithm"] = self.tuneUI.cb_Optimization_Algorithm.currentText()
        state_dict["Cell_Type_Optimization"] = self.tuneUI.cb_Cell_Type_Optimization.currentText()
        state_dict["UQ_Check"] = self.tuneUI.cb_UQ.checkState()

        # mid cell parameters
        state_dict["A_i_opt"] = self.tuneUI.le_A_i_opt.text()
        state_dict["B_i_opt"] = self.tuneUI.le_B_i_opt.text()
        state_dict["a_i_opt"] = self.tuneUI.le_a_i_opt.text()
        state_dict["b_i_opt"] = self.tuneUI.le_b_i_opt.text()
        state_dict["Ri_i_opt"] = self.tuneUI.le_Ri_i_opt.text()
        state_dict["L_i_opt"] = self.tuneUI.le_L_i_opt.text()
        state_dict["Req_i_opt"] = self.tuneUI.le_Req_i_opt.text()

        state_dict["Tune_Variable"] = self.tuneUI.cb_Tune_Variable.currentText()
        state_dict["Tune_Frequency"] = self.tuneUI.dsb_Tune_Frequency.value()
        state_dict["Norm_Length"] = self.tuneUI.db_Norm_Length.value()
        state_dict["N_Cells"] = self.tuneUI.sb_Norm_Length_N_Cells.value()
        state_dict["Processors_Count"] = self.tuneUI.sb_Processors_Count.value()

        state_dict["Initial_Points"] = self.tuneUI.sb_Initial_Points.value()
        state_dict["Max_Table_Size"] = self.tuneUI.sb_Max_Table_Size.value()
        state_dict["N_Generation"] = self.tuneUI.sb_N_Generation.value()

        state_dict["Init_Generation_Method"] = self.tuneUI.cb_Init_Generation_Method.currentText()
        state_dict["Sobol_Sequence_Index"] = self.tuneUI.sb_Sobol_Sequence_Index.value()
        state_dict["Optimize_By"] = self.tuneUI.cb_Optimize_By.currentText()
        state_dict["Crossover_Factor"] = self.tuneUI.sb_Crossover_Factor.value()
        state_dict["N_Elites_To_Cross"] = self.tuneUI.sb_N_Elites_To_Cross.value()
        state_dict["Mutation_Factor"] = self.tuneUI.sb_Mutation_Factor.value()
        state_dict["Chaos_Factor"] = self.tuneUI.sb_Chaos_Factor.value()

        state_dict["Populate_Objectives"] = self.tuneUI.ccb_Populate_Objectives.currentText()
        state_dict["Populate_Constraints"] = self.tuneUI.ccb_Populate_Constraints.currentText()

        state_dict["Expansion"] = self.tuneUI.cb_Expansion.checkState()
        state_dict["LBP"] = self.tuneUI.cb_LBP.checkState()

    def deserialize(self, state_dict):
        print("here at deserializin")
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
        self.tuneUI.le_L_i.setText(state_dict["L_i"])
        self.tuneUI.le_Req_i.setText(state_dict["Req_i"])

        self.tuneUI.le_A_o.setText(state_dict["A_o"])
        self.tuneUI.le_B_o.setText(state_dict["B_o"])
        self.tuneUI.le_a_o.setText(state_dict["a_o"])
        self.tuneUI.le_b_o.setText(state_dict["b_o"])
        self.tuneUI.le_Ri_o.setText(state_dict["Ri_o"])
        self.tuneUI.le_L_o.setText(state_dict["L_i"])

        # settings
        self.tuneUI.sb_No_Of_Processors_Tune.setValue(state_dict["No_Of_Processors"])
        self.tuneUI.cb_Iterative_Method.setCurrentIndex(state_dict["Iterative_Method"])
        self.tuneUI.le_Tolerance.setText(state_dict["Tolerance"])
        self.tuneUI.sb_Max_Iteration.setValue(state_dict["Max_Iteration"])

        ## Optimization control
        self.tuneUI.cb_Optimization_Algorithm.setCurrentText(state_dict["Optimization_Algorithm"])
        self.tuneUI.cb_Cell_Type_Optimization.setCurrentText(state_dict["Cell_Type_Optimization"])
        self.tuneUI.cb_UQ.setCheckState(state_dict["UQ_Check"])
        print("here at deserializin optimization")

        # mid cell parameters
        self.tuneUI.le_A_i_opt.setText(state_dict["A_i_opt"])
        self.tuneUI.le_B_i_opt.setText(state_dict["B_i_opt"])
        self.tuneUI.le_a_i_opt.setText(state_dict["a_i_opt"])
        self.tuneUI.le_b_i_opt.setText(state_dict["b_i_opt"])
        self.tuneUI.le_Ri_i_opt.setText(state_dict["Ri_i_opt"])
        self.tuneUI.le_L_i_opt.setText(state_dict["L_i_opt"])
        self.tuneUI.le_Req_i_opt.setText(state_dict["Req_i_opt"])

        self.tuneUI.cb_Tune_Variable.setCurrentText(state_dict["Tune_Variable"])
        self.tuneUI.dsb_Tune_Frequency.setValue(state_dict["Tune_Frequency"])
        self.tuneUI.db_Norm_Length.setValue(state_dict["Norm_Length"])
        self.tuneUI.sb_Norm_Length_N_Cells.setValue(state_dict["N_Cells"])
        self.tuneUI.sb_Processors_Count.setValue(state_dict["Processors_Count"])

        self.tuneUI.sb_Initial_Points.setValue(state_dict["Initial_Points"])
        self.tuneUI.sb_Max_Table_Size.setValue(state_dict["Max_Table_Size"])
        self.tuneUI.sb_N_Generation.setValue(state_dict["N_Generation"])
        self.tuneUI.cb_Init_Generation_Method.setCurrentText(state_dict["Init_Generation_Method"])
        self.tuneUI.sb_Sobol_Sequence_Index.setValue(state_dict["Sobol_Sequence_Index"])
        self.tuneUI.cb_Optimize_By.setCurrentText(state_dict["Optimize_By"])
        self.tuneUI.sb_Crossover_Factor.setValue(state_dict["Crossover_Factor"])
        self.tuneUI.sb_N_Elites_To_Cross.setValue(state_dict["N_Elites_To_Cross"])
        self.tuneUI.sb_Mutation_Factor.setValue(state_dict["Mutation_Factor"])
        self.tuneUI.sb_Chaos_Factor.setValue(state_dict["Chaos_Factor"])

        self.tuneUI.ccb_Populate_Objectives.setCurrentText(state_dict["Populate_Objectives"])
        self.tuneUI.ccb_Populate_Constraints.setCurrentText(state_dict["Populate_Constraints"])

        self.tuneUI.cb_Expansion.setCheckState(state_dict["Expansion"])
        self.tuneUI.cb_LBP.setCheckState(state_dict["LBP"])

    def ui_effects(self):

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.tuneUI.w_Tune_Input.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.tuneUI.w_Simulation_Controls.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.tuneUI.w_Inner_Cell.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.tuneUI.w_Outer_Cell.setGraphicsEffect(shadow)

        shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
        shadow.setColor(QColor(0, 0, 0, 77))

        self.tuneUI.w_Cell_Component_Control.setGraphicsEffect(shadow)

    @staticmethod
    def run_sequential(pseudo_shape_space_proc, resume, p, bc, parentDir, projectDir, filename, tuner_option,
                       tune_variable, iter_set, cell_type, progress_list, convergence_list, save_last, n_cells):

        tuner.tune(pseudo_shape_space_proc, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                   tuner_option=tuner_option, tune_variable=tune_variable, iter_set=iter_set,
                   cell_type=cell_type, progress_list=progress_list, convergence_list=convergence_list,
                   save_last=save_last, n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

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
    def text_to_list_old(txt):
        if "range" in txt:
            txt = txt.replace('range', '')
            l = ast.literal_eval(txt)
            return range(l[0], l[1], l[2])
        elif 'linspace' in txt:
            l = eval(f'np.{txt}')
            return l
        elif txt in "ABabRiLReq":
            return txt
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
                print("Done checking path. Yes selected.")
                return "Yes"
            elif x == 65536:
                print("Done checking path. No selected.")
                return "No"
            elif x == 4194304:
                print("Done checking path. Cancel selected.")
                return "Cancel"
        else:
            print("Path does not yet exist. Creating...")
            return "Yes"

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
    def proof_filename(dirc):
        # check if extension is included
        if dirc.split('.')[-1] != 'json':
            dirc = f'{dirc}.json'

        return dirc

    def change_tune_option(self, txt):
        if txt == 'Req':
            self.tuneUI.l_Tune_Alternate_Variable.setText("L")
            self.tuneUI.l_Tune_Alternate_Variable_End_Cell.setText('L')
            # self.tuneUI.le_Tune_Variable_End_Cell.setEnabled(True)
        else:
            self.tuneUI.l_Tune_Alternate_Variable.setText("Req")
            self.tuneUI.l_Tune_Alternate_Variable_End_Cell.setText("Req")

            # set end cell Req equal to mid cell Req
            # self.tuneUI.le_Tune_Variable_End_Cell.setText(self.tuneUI.le_Tune_Variable.text())
            # self.tuneUI.le_Tune_Variable_End_Cell.setEnabled(False)

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

    # @staticmethod
    def combine_dict(self, proc_count, filename, projectDir):
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

        return result

    @staticmethod
    def delete_process_dict(proc_count, projectDir):
        for index in range(proc_count):
            os.remove(fr'{projectDir}\Cavities\shape_space{index}.json')


class OptimizationControl:
    def __init__(self, tuneUI):
        self.tune_freq = None
        self.df_global = pd.DataFrame()
        self.tuneUI = tuneUI

        self.objs_dict = {}
        self.constraints_dict = {}

        self.ng_max = self.tuneUI.sb_N_Generation.value()
        self.objectives = []
        self.constraints = []
        self.weights = []
        self.processes = []
        self.df = None
        self.sbd = {}  # shape bounds dictionary
        self.poc = 5  # pareto optimal count
        self.parentDir = r"D:\Dropbox\CavityDesignHub"
        self.projectDir = fr"D:\Dropbox\CavityDesignHub\Cavity800"

        self.initUI()
        self.signals()

    def initUI(self):
        self.populate_parameters()
        self.tuneUI.w_Mid_Cell_Opt.hide()
        if self.tuneUI.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':
            self.tuneUI.w_Mid_Cell_Opt.show()

        self.tuneUI.w_Sobol_Sequence_Parameters.hide()
        if self.tuneUI.cb_Init_Generation_Method.currentText() == 'Sobol Sequence':
            self.tuneUI.w_Sobol_Sequence_Parameters.show()

    def signals(self):
        # add signal to show mid cell parameters for end-mid-cell optimization
        self.tuneUI.cb_Cell_Type_Optimization.currentTextChanged.connect(lambda: self.show_hide_mid_cell_params())

        # populate objective
        self.tuneUI.ccb_Populate_Objectives.currentTextChanged.connect(lambda: self.populate_objectives(
            self.tuneUI.ccb_Populate_Objectives.currentText().split(', '), self.tuneUI.ccb_Populate_Objectives,
            self.tuneUI.tw_Objectives, self.objs_dict))

        # add signal to update parameter influence combobox
        self.tuneUI.ccb_Populate_Objectives.currentTextChanged.connect(lambda: self.update_influence_ccb())

        # populate constraints
        self.tuneUI.ccb_Populate_Constraints.currentTextChanged.connect(lambda: self.populate_constraints(
            self.tuneUI.ccb_Populate_Constraints.currentText().split(', '), self.tuneUI.ccb_Populate_Constraints,
            self.tuneUI.tw_Constraints, self.constraints_dict))

        self.tuneUI.pb_Run_Optimization.clicked.connect(lambda: self.ea(0))

        self.tuneUI.cb_Init_Generation_Method.currentTextChanged.connect(
            lambda: self.tuneUI.w_Sobol_Sequence_Parameters.show()
            if self.tuneUI.cb_Init_Generation_Method.currentText() == 'Sobol Sequence'
            else self.tuneUI.w_Sobol_Sequence_Parameters.hide())

        # sobol sequence
        self.tuneUI.sb_Sobol_Sequence_Index.valueChanged.connect(lambda: self.tuneUI.sb_Initial_Points.setValue(int(2**self.tuneUI.sb_Sobol_Sequence_Index.value())))

    def ea(self, n):
        if n == 0:
            # update lists
            self.ng_max = self.tuneUI.sb_N_Generation.value()
            self.get_constraints()
            self.get_objectives()
            self.df = self.generate_first_men()

            # clear folder to avoid reading from previous optimization attempt
            folders = [fr"D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS",
                       fr"D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\ABCI"]

            # folders = [fr"D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS"]
            for folder in folders:
                for filename in os.listdir(folder):
                    try:
                        shutil.rmtree(fr"{folder}\{filename}")
                    except NotADirectoryError:
                        os.remove(fr"{folder}\{filename}")

        # optimize by page rank
        # remove the lowest ranking members
        df = self.df

        # compare with global dict and remove duplicates
        compared_cols = ['A', 'B', 'a', 'b', 'Ri']
        if not self.df_global.empty:
            df = df.loc[~df.set_index(compared_cols).index.isin(
                self.df_global.set_index(compared_cols).index)]  # this line of code removes duplicates

        # if self.tuneUI.cb_Optimize_By.currentText() == "Page Rankf":
        #     if len(self.df) > 15:
        #         df = self.df.iloc[:, 0:15]
        #     else:
        #         df = self.df.iloc[:, :]

        # naming convention G<generation number>_C<cavity number>_<type>
        # type refers to mutation M or crossover C

        # first_shapes dictionary
        # convert dataframe to tune dictionary style
        pseudo_shape_space = {}
        self.tune_freq = self.tuneUI.dsb_Tune_Frequency.value()
        for i, row in df.iterrows():
            # key, A, B, a, b, Ri, L, Req, alpha =
            rr = row.tolist()
            if self.tuneUI.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':

                A_i = self.check_input(self.tuneUI.le_A_i_opt.text())[0]
                B_i = self.check_input(self.tuneUI.le_B_i_opt.text())[0]
                a_i = self.check_input(self.tuneUI.le_a_i_opt.text())[0]
                b_i = self.check_input(self.tuneUI.le_b_i_opt.text())[0]
                Ri_i = self.check_input(self.tuneUI.le_Ri_i_opt.text())[0]
                L_i = self.check_input(self.tuneUI.le_L_i_opt.text())[0]
                Req_i = self.check_input(self.tuneUI.le_Req_i_opt.text())[0]
                alpha_i = self.check_input(self.tuneUI.le_Alpha_opt.text())[0]

                IC = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, alpha_i, alpha_i]
                # print(type(IC))
                # print(type(IC))
                # print(type(rr[1:]))

                pseudo_shape_space[rr[0]] = {'IC': IC, 'OC': rr[1:], 'BP': 'left', 'FREQ': self.tune_freq}
            else:
                pseudo_shape_space[rr[0]] = {'IC': rr[1:], 'OC': rr[1:], 'BP': 'none', 'FREQ': self.tune_freq}

        pseudo_shape_space = self.remove_duplicate_values(pseudo_shape_space)

        # pseudo_shape_space_name = f'{self.projectDir}/Cavities/pseudo_{self.proof_filename("EA.json")}'
        #
        # with open(pseudo_shape_space_name, 'w') as file:
        #     file.write(json.dumps(pseudo_shape_space, indent=4, separators=(',', ': ')))

        ############################
        # run tune
        n_cells = self.tuneUI.sb_Norm_Length_N_Cells.value()
        norm_length = self.tuneUI.db_Norm_Length.value()

        self.run_tune_parallel(pseudo_shape_space, n_cells)

        for o in self.objectives:
            if o[1] in ["freq", "Epk/Eacc", "Bpk/Eacc", "R/Q", "G", "Q0"]:
                # process tune results
                params = {}
                obj_result = []
                tune_result = []
                processed_keys = []
                for key in pseudo_shape_space.keys():
                    filename = fr'{self.projectDir}\SimulationData\SLANS\Cavity{key}\cavity_33.svl'
                    try:
                        params = fr.svl_reader(filename)
                        obj = self.get_objectives_value(params, self.objectives, norm_length, n_cells)
                        # print(obj, tr)
                        # read tune results
                        d_tune_res = {'Req': 0, 'L': 0, 'freq': 0}
                        with open(fr'{self.projectDir}\SimulationData\SLANS\Cavity{key}\tune_res.json', "r") as infile:
                            d_tune_res = json.load(infile)

                        obj_result.append(obj)
                        tune_result.append(list(d_tune_res.values()))
                        processed_keys.append(key)
                    except FileNotFoundError:
                        pass

                # after removing duplicates, dataframe might change size
                df = df.loc[df['key'].isin(processed_keys)]

                # replace Req with tuned Req
                df[['Req', 'L', 'alpha_i', 'alpha_o', 'freq']] = tune_result

                # append objective column to dataframe only from slans results
                obj_slans = [o[1] for o in self.objectives if 'Z' not in o[1]]
                ic(obj_result)
                ic(obj_slans)
                df[obj_slans] = obj_result
            # print("In here to evaluate wakfield")
            elif "ZL" in o[1] or "ZT" in o[1] or "k_loss" in o[1] or "k_kick" in o[1]:
                print("In here to evaluate wakfield0")
                # run wakefield analysis and return shape space
                wake_shape_space = self.run_wakefield_parallel(df)
                # print("In here to evaluate wakfield1")
                # process wakefield results
                df_wake, processed_keys = self.get_wakefield_objectives_value(wake_shape_space,
                                                                                 fr'{self.projectDir}\SimulationData\ABCI')
                ic(df_wake)
                df = df.merge(df_wake, on='key', how='inner')
                ic(df)
                # print("In here to evaluate wakfield2")
                # remove unprocessed keys from dataframe
                # df = df.loc[df['key'].isin(processed_keys)]
                #
                # # append wakefield objective column to dataframe
                # df.loc[:, np.array(obj_abci)[:, 1]] = obj_result

                # print("In here to evaluate wakfield3")
                break

        # apply UQ
        if self.tuneUI.cb_UQ.isChecked():
            solver_dict = {'slans': slans_geom_seq, 'abci': abci_geom}
            beampipes = {'Mid Cell': 'none', 'End-End Cell': 'left', 'End-Mid Cell': 'left', 'Single Cell': 'both'}
            beampipe = beampipes[self.tuneUI.cb_Cell_Type_Optimization.currentText()]

            # if self.tuneUI.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':
            #     n_cells = 5
            # else:
            #     n_cells = 5

            solver_args_dict = {'slans': {'n_cells': n_cells, 'n_modules': 1, 'f_shift': 0, 'bc': 33,
                                          'beampipes': beampipe,
                                          'norm_length': self.tuneUI.db_Norm_Length.value(),
                                          'parentDir': self.parentDir, 'projectDir': self.projectDir},
                                'abci': {'n_cells': n_cells, 'n_modules': 1,
                                         'MROT': 2, 'MT': 4, 'NFS': 10000, 'UBT': 50, 'bunch_length': 25,
                                         'DDR_SIG': 0.1, 'DDZ_SIG': 0.1,
                                         'parentDir': self.parentDir, 'projectDir': self.projectDir, 'progress_list': None,
                                         'WG_M': None, 'marker': ''}
                                }

            shape_space = self.uq_parallel(df, self.objectives, solver_dict, solver_args_dict)
            ic(shape_space)

            # get uq_parameters
            uq_result_dict = {}
            for key in shape_space.keys():
                filename_slans = fr'{self.projectDir}\SimulationData\SLANS\Cavity{key}\uq.json'
                filename_abci = fr'{self.projectDir}\SimulationData\ABCI\Cavity{key}\uq.json'
                if os.path.exists(filename_slans):# and os.path.exists(filename_abci):
                    uq_result_dict[key] = []
                    with open(filename_slans, "r") as infile:
                        uq_d = json.load(infile)
                        for o in self.objectives:
                            if o[1] in ["Req", "freq", "Q", "E", "R/Q", "GR/Q", "Epk/Eacc", "Bpk/Eacc"]:
                                uq_result_dict[key].append(uq_d[o[1]]['expe'][0])
                                uq_result_dict[key].append(uq_d[o[1]]['stdDev'][0])
                                if o[0] == 'min':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] + 6*uq_d[o[1]]['stdDev'][0])
                                elif o[0] == 'max':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] - 6*uq_d[o[1]]['stdDev'][0])
                                else:
                                    # for equal, calculate |expected_value - design_value| + 6sigma
                                    uq_result_dict[key].append(np.abs(uq_d[o[1]]['expe'][0] - o[2]) + uq_d[o[1]]['stdDev'][0])

                                # ic(uq_result_dict)

                # filename = fr'{self.projectDir}\SimulationData\ABCI\Cavity{key}\uq.json'
                if os.path.exists(filename_abci):
                    if key not in uq_result_dict:
                        uq_result_dict[key] = []

                    with open(filename_abci, "r") as infile:
                        uq_d = json.load(infile)
                        for o in self.objectives:
                            if o[1] not in ["Req", "freq", "Q", "E", "R/Q", "Epk/Eacc", "Bpk/Eacc"]:
                                uq_result_dict[key].append(uq_d[o[1]]['expe'][0])
                                uq_result_dict[key].append(uq_d[o[1]]['stdDev'][0])
                                if o[0] == 'min':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] + 6*uq_d[o[1]]['stdDev'][0])
                                elif o[0] == 'max':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] - 6*uq_d[o[1]]['stdDev'][0])

            # ic(uq_result_dict)
            uq_column_names = []
            for o in self.objectives:
                uq_column_names.append(fr'E[{o[1]}]')
                uq_column_names.append(fr'std[{o[1]}]')
                if o[0] == 'min':
                    uq_column_names.append(fr'E[{o[1]}] + 6*std[{o[1]}]')
                elif o[0] == 'max':
                    uq_column_names.append(fr'E[{o[1]}] - 6*std[{o[1]}]')
                else:
                    uq_column_names.append(fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]')

            ic(df)
            df_uq = pd.DataFrame.from_dict(uq_result_dict, orient='index')
            ic(df_uq)
            df_uq.columns = uq_column_names
            ic(uq_column_names)
            df_uq.index.name = 'key'
            df_uq.reset_index(inplace=True)
            ic(df_uq)
            df = df.merge(df_uq, on='key', how='inner')
            ic(df)
        ######################################################
        # filter by constraints
        for const in self.constraints:
            c = const.split(" ")

            if c[1] == '>':
                df = df.loc[(df[f'{c[0]}'] > float(c[2]))]
            elif c[1] == '<':
                df = df.loc[(df[f'{c[0]}'] < float(c[2]))]
            elif c[1] == '<=':
                df = df.loc[(df[f'{c[0]}'] <= float(c[2]))]
            elif c[1] == '>=':
                df = df.loc[(df[f'{c[0]}'] >= float(c[2]))]
            elif c[1] == '==':
                df = df.loc[(df[f'{c[0]}'] == float(c[2]))]

        # update with global dataframe
        if not self.df_global.empty:
            df = pd.concat([self.df_global, df], ignore_index=True)

            # # drop duplicates
            # df = df.drop_duplicates(subset=['A', 'B', 'a', 'b', 'Ri', 'L'], keep='first')

        # reset total rank
        df['total_rank'] = 0

        # rank shapes by objectives
        for i, obj in enumerate(self.objectives):
            if self.tuneUI.cb_UQ.isChecked():
                if obj[0] == "min":
                    df[f'rank_E[{obj[1]}] + 6*std[{obj[1]}]'] = df[fr'E[{obj[1]}] + 6*std[{obj[1]}]'].rank() * self.weights[i]
                    ic(df[f'rank_E[{obj[1]}] + 6*std[{obj[1]}]'])
                elif obj[0] == "max":
                    df[f'rank_E[{obj[1]}] - 6*std[{obj[1]}]'] = df[fr'E[{obj[1]}] - 6*std[{obj[1]}]'].rank(ascending=False) * self.weights[i]
                    ic(df[f'rank_E[{obj[1]}] - 6*std[{obj[1]}]'])
                elif obj[0] == "equal":
                    df[fr'rank_|E[{obj[1]}] - {obj[2]}| + std[{obj[1]}]'] = df[fr'|E[{obj[1]}] - {obj[2]}| + std[{obj[1]}]'].rank() * self.weights[i]
                    ic(df[fr'rank_|E[{obj[1]}] - {obj[2]}| + std[{obj[1]}]'])
    
                # if 'total_rank' in df.columns:
                if obj[0] == 'min':
                    df[f'total_rank'] = df[f'total_rank'] + df[f'rank_E[{obj[1]}] + 6*std[{obj[1]}]']
                elif obj[0] == 'max':
                    df[f'total_rank'] = df[f'total_rank'] + df[f'rank_E[{obj[1]}] - 6*std[{obj[1]}]']
                else:
                    df[f'total_rank'] = df[f'total_rank'] + df[fr'rank_|E[{obj[1]}] - {obj[2]}| + std[{obj[1]}]']
                # else:
                #     if obj[0] == 'min':
                #         df[f'total_rank'] = df[f'rank_E[{obj[1]}] + 6*std[{obj[1]}]']
                #     elif obj[0] == 'max':
                #         df[f'total_rank'] = df[f'rank_E[{obj[1]}] - 6*std[{obj[1]}]']
                #     else:
                #         df[f'total_rank'] = df[f'rank_std[{obj[1]}]']
            else:
                if obj[0] == "min":
                    df[f'rank_{obj[1]}'] = df[obj[1]].rank() * self.weights[i]
                elif obj[0] == "max":
                    df[f'rank_{obj[1]}'] = df[obj[1]].rank(ascending=False) * self.weights[i]
                elif obj[0] == "equal":  # define properly later
                    continue
    
                # if 'total_rank' in df.columns:
                df[f'total_rank'] = df[f'total_rank'] + df[f'rank_{obj[1]}']
                # else:
                #     df[f'total_rank'] = df[f'rank_{obj[1]}']

        # reorder
        ic(df)
        tot = df.pop(f'total_rank')
        ic(tot)
        df[f'total_rank'] = tot/sum(self.weights)  # normalize by sum of weights
        ic(df)

        # order shapes by rank
        df = df.sort_values(by=['total_rank'])
        df = df.reset_index(drop=True)

        # pareto condition
        reorder_indx = self.pareto_front(df)
        df = df.loc[reorder_indx, :]
        # reset index
        df = df.dropna().reset_index(drop=True)
        ic(df)

        # update global
        if len(df) > self.tuneUI.sb_Max_Table_Size.value():
            # self.df_global = df.loc[0:self.tuneUI.cb_Max_Table_Size.value(), :]
            self.df_global = df
        else:
            self.df_global = df
        ic(self.df_global)
        
        # check if df_global is empty
        if self.df_global.shape[0] == 0:
            ic("Unfortunately, none survived the constraints and the program has to end. I can't even say that this was a good run.")
            return

        # save dataframe
        filename = fr"{self.projectDir}\SimulationData\SLANS\Generation{n}.xlsx"
        self.recursive_save(self.df_global, filename, reorder_indx)
        # ic(self.df_global)

        # birth next generation
        # crossover
        print("Crossover")
        df_cross = self.crossover(df, n, self.tuneUI.sb_Crossover_Factor.value())  # , elites["GR/Q
        # ic(df_cross)

        # mutation
        print("Mutation")
        df_mutation = self.mutation(df, n, self.tuneUI.sb_Mutation_Factor.value())
        ic(df_mutation)

        # chaos
        print("Chaos")
        df_chaos = self.chaos(self.tuneUI.sb_Chaos_Factor.value(), n)
        # ic(df_chaos)

        # take elites from previous generation over to next generation
        df_ng = pd.concat([df_cross, df_mutation, df_chaos], ignore_index=True)

        # update dictionary
        self.df = df_ng

        n += 1
        print(n)
        print("=" * 80)
        if n < self.ng_max:
            return self.ea(n)
        else:
            return

    def uq_parallel(self, df, objectives, solver_dict, solver_args_dict):
        proc_count = self.tuneUI.sb_Processors_Count.value()
        # marker = self.wakefieldUI.le_Marker.text()

        # get geometric parameters
        df = df.loc[:, ['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"]]
        shape_space = {}

        df = df.set_index('key')
        for index, row in df.iterrows():
            rw = row.tolist()
            if self.tuneUI.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':

                A_i = self.check_input(self.tuneUI.le_A_i_opt.text())[0]
                B_i = self.check_input(self.tuneUI.le_B_i_opt.text())[0]
                a_i = self.check_input(self.tuneUI.le_a_i_opt.text())[0]
                b_i = self.check_input(self.tuneUI.le_b_i_opt.text())[0]
                Ri_i = self.check_input(self.tuneUI.le_Ri_i_opt.text())[0]
                L_i = self.check_input(self.tuneUI.le_L_i_opt.text())[0]
                Req_i = self.check_input(self.tuneUI.le_Req_i_opt.text())[0]
                alpha_i = self.check_input(self.tuneUI.le_Alpha_opt.text())[0]

                IC = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, alpha_i]

                shape_space[f'{index}'] = {'IC': IC, 'OC': rw, 'OC_R': rw}
            else:
                shape_space[f'{index}'] = {'IC': rw, 'OC': rw, 'OC_R': rw}

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

                solver_args_dict['slans']['proc'] = p
                solver_args_dict['abci']['proc'] = p
                service = mp.Process(target=self.uq,
                                     args=(processor_shape_space, objectives, solver_dict, solver_args_dict))
                service.start()
                self.processes.append(service)
            except Exception as e:
                print(f"Exception in uq_parallel analysis:: {e}")

        for p in self.processes:
            p.join()

        self.processes = []
        return shape_space

    @staticmethod
    def uq(shape_space, objectives, solver_dict, solver_args_dict):
        for key, shape in shape_space.items():
            err = False
            result_dict_slans, result_dict_abci = {}, {}
            run_slans, run_abci = False, False
            slans_obj_list, abci_obj_list = [], []
            for o in objectives:

                if o[1] in ["Req", "freq", "Q", "E", "R/Q", "GR/Q", "Epk/Eacc", "Bpk/Eacc"]:
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

            if run_slans:
                Ttab_val_f = []
                solver, solver_args = solver_dict['slans'], solver_args_dict['slans']
                bc = solver_args['bc']
                beampipes = solver_args['beampipes']
                norm_length = solver_args['norm_length']
                n_cells = solver_args['n_cells']
                proc = solver_args['proc']
                parentDir = solver_args['parentDir']
                projectDir = solver_args['projectDir']
                sub_dir = fr'Cavity{key}'  # the simulation runs at the quadrature points are saved to the key of the mean value run
                for i in range(no_sims):
                    skip = False
                    p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
                    p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
                    p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
                    p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
                    p_init[4] = p_true[4] * (1 + delta * nodes[4, i])

                    par_mid = np.append(p_init, shape['IC'][5:]).tolist()
                    par_end = par_mid

                    # perform checks on geometry
                    ok = perform_geometry_checks(par_mid, par_end)
                    if not ok:
                        err = True
                        break

                    fid = fr'{key}_Q{i}'

                    # check if folder exists and skip if it does
                    if os.path.exists(fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}'):
                        skip = True
                        # ic("Skipped: ", fid, fr'{projectDir}\SimulationData\ABCI\Cavity{key}\Cavity{fid}')

                    # skip analysis if folder already exists.
                    if not skip:
                        #  run model using SLANS or CST
                        # # create folders for all keys
                        solver.createFolder(fid, projectDir, subdir=sub_dir)

                        solver.cavity(n_cells, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, fid=fid,
                                      parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
                    filename = fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}\cavity_33.svl'

                    if os.path.exists(filename):
                        params = fr.svl_reader(filename)
                        obj_result, tune_result = get_objectives_value(params, slans_obj_list, norm_length, n_cells)

                        # sometimes some degenerate shapes are still generated and the solver returns zero
                        # for the objective functions, such shapes are considered invalid
                        for objr in obj_result:
                            if objr == 0:
                                # skip key
                                err = True
                                break

                        # ic('SLANS', obj_result)
                        tab_val_f = obj_result

                        Ttab_val_f.append(tab_val_f)
                    else:
                        err = True

                # # add original point
                # filename = fr'{projectDir}\SimulationData\SLANS\Cavity{key}\cavity_33.svl'
                # params = fr.svl_reader(filename)
                # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
                # tab_val_f = obj_result
                    # Ttab_val_f.append(tab_val_f)

                if err:
                    break

                v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
                # ic(v_expe_fobj, v_stdDev_fobj)
                # append results to dict

                for i, o in enumerate(slans_obj_list):
                    result_dict_slans[o[1]]['expe'].append(v_expe_fobj[i])
                    result_dict_slans[o[1]]['stdDev'].append(v_stdDev_fobj[i])

                with open(fr"{projectDir}\SimulationData\SLANS\Cavity{key}\uq.json", 'w') as file:
                    file.write(json.dumps(result_dict_slans, indent=4, separators=(',', ': ')))

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
                sub_dir = fr'Cavity{key}'  # the simulation runs at the quadrature points are saved to the key of the mean value run
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
                    if os.path.exists(fr'{projectDir}\SimulationData\ABCI\Cavity{key}\Cavity{fid}'):
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
                    abci_folder = fr'{projectDir}\SimulationData\ABCI\Cavity{key}'
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

                    with open(fr"{projectDir}\SimulationData\ABCI\Cavity{key}\uq.json", 'w') as file:
                        file.write(json.dumps(result_dict_abci, indent=4, separators=(',', ': ')))

    def stroud(self, p):
        # Stroud-3 method
        #
        # Input parameters:
        #  p   number of dimensions
        # Output parameters:
        #  nodes   nodes of quadrature rule in [0,1]^p (column-wise)
        #

        nodes = np.zeros((p, 2 * p))
        coeff = np.pi / p
        fac = np.sqrt(2 / 3)

        for i in range(2 * p):
            for r in range(int(np.floor(0.5 * p))):
                k = 2 * r
                nodes[k, i] = fac * np.cos((k + 1) * (i + 1) * coeff)
                nodes[k + 1, i] = fac * np.sin((k + 1) * (i + 1) * coeff)

            if 0.5 * p != np.floor(0.5 * p):
                nodes[-1, i] = ((-1) ** (i + 1)) / np.sqrt(3)

        # transform nodes from [-1,+1]^p to [0,1]^p
        nodes = 0.5 * nodes + 0.5

        return nodes

    def quad_stroud3(self, rdim, degree):
        # data for Stroud-3 quadrature in [0,1]^k
        # nodes and weights
        nodes = self.stroud(rdim)
        nodestr = 2. * nodes - 1.
        weights = (1 / (2 * rdim)) * np.ones((2 * rdim, 1))

        # evaluation of Legendre polynomials
        bpoly = np.zeros((degree + 1, rdim, 2 * rdim))
        for l in range(rdim):
            for j in range(2 * rdim):
                bpoly[0, l, j] = 1
                bpoly[1, l, j] = nodestr[l, j]
                for i in range(1, degree):
                    bpoly[i + 1, l, j] = ((2 * (i + 1) - 1) * nodestr[l, j] * bpoly[i, l, j] - i * bpoly[
                        i - 1, l, j]) / (i + 1)

        # standardisation of Legendre polynomials
        for i in range(1, degree + 1):
            bpoly[i, :, :] = bpoly[i, :, :] * np.sqrt(2 * (i + 1) - 1)

        return nodes, weights, bpoly

    def weighted_mean_obj(self, tab_var, weights):
        rows_sims_no, cols = np.shape(tab_var)
        no_weights, dummy = np.shape(weights)  # z funckji quadr_stroud wekt columnowy

        if rows_sims_no == no_weights:
            expe = np.zeros((cols, 1))
            outvar = np.zeros((cols, 1))
            for i in range(cols):
                expe[i, 0] = np.dot(tab_var[:, i], weights)
                outvar[i, 0] = np.dot(tab_var[:, i]**2, weights)
            stdDev = np.sqrt(outvar - expe**2)
        else:
            expe = 0
            stdDev = 0
            ic('Cols_sims_no != No_weights')

        return list(expe.T[0]), list(stdDev.T[0])

    def run_tune_parallel(self, pseudo_shape_space, n_cells):
        proc_count = self.tuneUI.sb_Processors_Count.value()

        # split shape_space for different processes/ MPI share process by rank
        keys = list(pseudo_shape_space.keys())

        # check if number of processors selected is greater than the number of keys in the pseudo shape space
        if proc_count > len(keys):
            proc_count = len(keys)

        shape_space_len = len(keys)
        share = floor(shape_space_len / proc_count)

        self.processes = []
        for p in range(proc_count):
            if True:
                if p < proc_count - 1:
                    proc_keys_list = keys[p * share:p * share + share]
                else:
                    proc_keys_list = keys[p * share:]

                self.overwriteFolder(p, self.projectDir)
                self.copyFiles(p, self.parentDir, self.projectDir)

                processor_shape_space = {}
                for key, val in pseudo_shape_space.items():
                    if key in proc_keys_list:
                        # check if folder alsready exists
                        if not os.path.exists(fr'{self.projectDir}\SimulationData\SLANS\Cavity{key}\cavity_33.svl'):
                            processor_shape_space[key] = val

                if 'End' in self.tuneUI.cb_Cell_Type_Optimization.currentText():
                    cell_type = 'End Cell'
                else:
                    cell_type = 'Mid Cell'

                tune_variable = self.tuneUI.cb_Tune_Variable.currentText()
                service = mp.Process(target=self.run_sequential,
                                     args=(processor_shape_space, "Yes", p, 33, r"D:\Dropbox\CavityDesignHub",
                                           r"D:\Dropbox\CavityDesignHub\Cavity800", "EA.json", 'PyTuner',
                                           tune_variable, ["Linear Interpolation", 1e-3, 10], cell_type, [],
                                           None, True, n_cells))
                service.start()
                self.processes.append(service)

        for p in self.processes:
            p.join()

        self.processes = []

    def run_wakefield_parallel(self, df):
        # get analysis parameters
        n_cells = 5

        n_modules = 1

        # change later
        WG_M = ['']  # half length of beam pipe between cavities in module

        # change all of these later
        MROT = 2  # run both longitudinal and transverse wakefield analysis
        MT = 4  # number of time steps for a beam to move one cell to another default = 3
        bunch_length = 25
        NFS = 10000  # Number of samples in FFT (max 10000)
        UBT = 50  # Wakelength in m
        DDZ_SIG = 0.1
        DDR_SIG = 0.1
        proc_count = self.tuneUI.sb_Processors_Count.value()
        # marker = self.wakefieldUI.le_Marker.text()

        # get geometric parameters
        df = df.loc[:, ['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"]]
        shape_space = {}

        df = df.set_index('key')
        for index, row in df.iterrows():
            rw = row.tolist()
            if self.tuneUI.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':

                A_i = self.check_input(self.tuneUI.le_A_i_opt.text())[0]
                B_i = self.check_input(self.tuneUI.le_B_i_opt.text())[0]
                a_i = self.check_input(self.tuneUI.le_a_i_opt.text())[0]
                b_i = self.check_input(self.tuneUI.le_b_i_opt.text())[0]
                Ri_i = self.check_input(self.tuneUI.le_Ri_i_opt.text())[0]
                L_i = self.check_input(self.tuneUI.le_L_i_opt.text())[0]
                Req_i = self.check_input(self.tuneUI.le_Req_i_opt.text())[0]
                alpha_i = self.check_input(self.tuneUI.le_Alpha_opt.text())[0]

                IC = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, alpha_i]

                shape_space[f'{index}'] = {'IC': IC, 'OC': rw, 'OC_R': rw}
            else:
                shape_space[f'{index}'] = {'IC': rw, 'OC': rw, 'OC_R': rw}

        # ic(shape_space)

        # shape_space = self.get_geometric_parameters('ABCI')

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

                service = mp.Process(target=run_sequential_wakefield,
                                     args=(n_cells, n_modules, processor_shape_space,
                                           MROT, MT, NFS, UBT, bunch_length,
                                           DDR_SIG, DDZ_SIG,
                                           self.parentDir,
                                           self.projectDir, [],
                                           '', ''
                                           ))

                service.start()
                self.processes.append(service)
            except Exception as e:
                print(f"Exception in run_MP for wakefield analysis:: {e}")

        for p in self.processes:
            p.join()

        self.processes = []
        return shape_space

    def generate_first_men(self):
        # get range from table
        tw = self.tuneUI.tw_Parameters
        for i in range(self.tuneUI.tw_Parameters.rowCount()):
            self.sbd[i] = [tw.cellWidget(i, 1).value(), tw.cellWidget(i, 2).value()]

        # populate
        initial_points = self.tuneUI.sb_Initial_Points.value()
        if self.tuneUI.cb_Init_Generation_Method.currentText() == "Uniform":
            data = {'key': [f"G0_C{i}_P" for i in range(initial_points)],
                    'A': np.linspace(self.sbd[0][0], self.sbd[0][1], initial_points),
                    'B': np.linspace(self.sbd[1][0], self.sbd[1][1], initial_points),
                    'a': np.linspace(self.sbd[2][0], self.sbd[2][1], initial_points),
                    'b': np.linspace(self.sbd[3][0], self.sbd[3][1], initial_points),
                    'Ri': np.linspace(self.sbd[4][0], self.sbd[4][1], initial_points),
                    'L': np.linspace(self.sbd[5][0], self.sbd[5][1], initial_points),
                    'Req': np.linspace(self.sbd[6][0], self.sbd[6][1] + 1, initial_points),
                    'alpha_i': np.zeros(initial_points),
                    'alpha_o': np.zeros(initial_points)}
            return pd.DataFrame.from_dict(data)
        elif self.tuneUI.cb_Init_Generation_Method.currentText() == "Random":
            data = {'key': [f"G0_C{i}_P" for i in range(initial_points)],
                    'A': random.sample(list(np.linspace(self.sbd[0][0], self.sbd[0][1], initial_points * 2)),
                                       initial_points),
                    'B': random.sample(list(np.linspace(self.sbd[1][0], self.sbd[1][1], initial_points * 2)),
                                       initial_points),
                    'a': random.sample(list(np.linspace(self.sbd[2][0], self.sbd[2][1], initial_points * 2)),
                                       initial_points),
                    'b': random.sample(list(np.linspace(self.sbd[3][0], self.sbd[3][1], initial_points * 2)),
                                       initial_points),
                    'Ri': random.sample(list(np.linspace(self.sbd[4][0], self.sbd[4][1], initial_points * 2)),
                                        initial_points),
                    'L': random.sample(list(np.linspace(self.sbd[5][0], self.sbd[5][1], initial_points * 2)),
                                       initial_points),
                    'Req': random.sample(list(np.linspace(self.sbd[6][0], self.sbd[6][1] + 1, initial_points * 2)),
                                         initial_points),
                    'alpha_i': np.zeros(initial_points),
                    'alpha_o': np.zeros(initial_points)}
            return pd.DataFrame.from_dict(data)
        else:
            ic(self.sbd.values())
            columns = ['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']
            dim = len(columns)
            index = self.tuneUI.sb_Sobol_Sequence_Index.value()
            l_bounds = np.array(list(self.sbd.values()))[:, 0]
            u_bounds = np.array(list(self.sbd.values()))[:, 1]

            const_var = []
            for i in range(dim-1, -1, -1):
                if l_bounds[i] == u_bounds[i]:
                    print("it got here")
                    const_var.append([columns[i], l_bounds[i]])
                    del columns[i]
                    l_bounds = np.delete(l_bounds, i)
                    u_bounds = np.delete(u_bounds, i)

            reduced_dim = len(columns)
            sampler = qmc.Sobol(d=reduced_dim, scramble=False)
            _ = sampler.reset()
            sample = sampler.random_base2(m=index)
            ic(qmc.discrepancy(sample))
            ic(l_bounds, u_bounds)
            sample = qmc.scale(sample, l_bounds, u_bounds)

            df = pd.DataFrame()
            df['key'] = [f"G0_C{i}_P" for i in range(initial_points)]
            df[columns] = sample

            for i in range(len(const_var)-1, -1, -1):
                df[const_var[i][0]] = np.ones(initial_points)*const_var[i][1]

            df['alpha_i'] = np.zeros(initial_points)
            df['alpha_o'] = np.zeros(initial_points)
            ic("sobol seq", df)

            return df

    def get_objectives(self):
        # objectives, weights, constraints, n, ng_max = [["min", "Epk/Eacc"], ["min", "Bpk/Eacc"], ["max", "R/Q"]], [5, 1, 1], \
        #                                               ["Bpk/Eacc < 6.5", "freq > 400.58", 'freq < 400.99'], 1, 4
        obj_vars = self.tuneUI.ccb_Populate_Objectives.currentText().split(', ')
        self.objectives = []
        self.weights = []

        for i, obj_var in enumerate(obj_vars):
            if obj_var == "ZL" or obj_var == "ZT":
                goal = self.tuneUI.tw_Objectives.cellWidget(i, 1).currentText()
                freq_ranges = self.process_interval(
                    self.text_to_list(self.tuneUI.tw_Objectives.cellWidget(i, 2).text()))
                for f in freq_ranges:
                    self.objectives.append([goal, f"{obj_var} [max({f[0]}<f<{f[1]})]", f])
                    self.weights.append(1)
            else:
                goal = self.tuneUI.tw_Objectives.cellWidget(i, 1).currentText()
                if goal == 'equal':
                    value = np.float(self.tuneUI.tw_Objectives.item(i, 2).text())
                    self.objectives.append([goal, obj_var, value])
                else:
                    self.objectives.append([goal, obj_var])
                self.weights.append(1)

        ic(self.objectives)
        ic(self.weights)

    def get_constraints(self):
        const_vars = self.tuneUI.ccb_Populate_Constraints.currentText().split(', ')
        self.constraints = []

        if not (const_vars == ['All'] or const_vars == ['']):
            for i, const_var in enumerate(const_vars):
                lower_bound = self.tuneUI.tw_Constraints.cellWidget(i, 1).value()
                upper_bound = self.tuneUI.tw_Constraints.cellWidget(i, 2).value()

                if upper_bound != lower_bound:
                    self.constraints.append(fr'{const_var} > {lower_bound}')
                    self.constraints.append(fr'{const_var} < {upper_bound}')
                else:
                    self.constraints.append(fr'{const_var} = {lower_bound}')

        # parameter bounds as constraints
        tw = self.tuneUI.tw_Parameters
        for i in range(tw.rowCount()):
            par = tw.cellWidget(i, 0).text()
            lower_bound = tw.cellWidget(i, 1).value()
            upper_bound = tw.cellWidget(i, 2).value()
            if upper_bound != lower_bound:
                self.constraints.append(fr'{par} > {lower_bound}')
                self.constraints.append(fr'{par} < {upper_bound}')
            else:
                self.constraints.append(fr'{par} = {lower_bound}')

        ic(self.constraints)

    def get_objectives_value(self, d, obj, norm_length, n_cells):
        Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
        L = d['LENGTH'][n_cells - 1] * 10  # convert to mm
        Freq = d['FREQUENCY'][n_cells - 1]
        E_stored = d['STORED ENERGY'][n_cells - 1]
        Rsh = d['SHUNT IMPEDANCE'][n_cells - 1]  # MOhm
        Q = d['QUALITY FACTOR'][n_cells - 1]
        Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
        Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
        # Vacc = dict['ACCELERATION'][n_cells - 1]
        Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells - 1]  # MV/m
        r_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm
        G = 0.00948*Q*(Freq/1300)
        GR_Q = G * 2 * r_Q

        Vacc = np.sqrt(
            2 * r_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6  # factor of 2, remember circuit and accelerator definition
        # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
        Eacc = Vacc / (
                n_cells * norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
        Epk_Eacc = Epk / Eacc
        Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc

        d = {
            "Req": Req,
            "L": L,
            "freq": Freq,
            "Q": Q,
            "E": E_stored,
            "R/Q": 2 * r_Q,
            "Epk/Eacc": Epk_Eacc,
            "Bpk/Eacc": Bpk_Eacc,
            "G": G,
            "GR/Q": GR_Q
        }

        objective = []
        # tune_result = []

        # # append freq and Req
        # if self.tuneUI.cb_Tune_Variable.currentText() == "Req":
        #     tune_result.append(Req)
        # else:
        #     tune_result.append(L)
        #
        # tune_result.append(Freq)

        # append objective functions
        for o in obj:
            if o[1] in d.keys():
                objective.append(d[o[1]])

        return objective #, tune_result

    def get_wakefield_objectives_value(self, d, abci_data_dir):
        k_loss_array_transverse = []
        k_loss_array_longitudinal = []
        k_loss_M0 = []
        key_list = []

        # create list to hold Z
        Zmax_mon_list = []
        Zmax_dip_list = []
        xmax_mon_list = []
        xmax_dip_list = []
        processed_keys_mon = []
        processed_keys_dip = []

        def calc_k_loss():
            for key, value in d.items():
                print(f"Processing for Cavity {key}")
                abci_data_long = ABCIData(abci_data_dir, key, 0)
                abci_data_trans = ABCIData(abci_data_dir, key, 1)

                # trans
                x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                k_loss_trans = abci_data_trans.loss_factor['Transverse']

                if math.isnan(k_loss_trans):
                    print_(f"Encountered an exception: Check shape {key}")
                    continue

                # long
                x, y, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                abci_data_long.get_data('Loss Factor Spectrum Integrated up to F')

                k_M0 = abci_data_long.y_peaks[0]
                k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
                k_loss_HOM = k_loss_long - k_M0

                # append only after successful run
                k_loss_M0.append(k_M0)
                k_loss_array_longitudinal.append(k_loss_HOM)
                k_loss_array_transverse.append(k_loss_trans)

            return [k_loss_M0, k_loss_array_longitudinal, k_loss_array_transverse]

        def get_Zmax_L(mon_interval=None):
            # print("2a")
            if mon_interval is None:
                mon_interval = [0.0, 2e10]
            # print("2b")

            for key, value in d.items():
                # print("2c")
                print(f"Processing for Cavity {key}")
                try:
                    abci_data_mon = ABCIData(abci_data_dir, f"Cavity{key}", 0)

                    # get longitudinal and transverse impedance plot data
                    xr_mon, yr_mon, _ = abci_data_mon.get_data('Real Part of Longitudinal Impedance')
                    xi_mon, yi_mon, _ = abci_data_mon.get_data('Imaginary Part of Longitudinal Impedance')

                    # print("2d")
                    # Zmax
                    if mon_interval is None:
                        mon_interval = [[0.0, 10]]

                    # calculate magnitude
                    ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]

                    # print("2e")
                    # get peaks
                    peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
                    xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

                    # print("2f", mon_interval)
                    for i, z_bound in enumerate(mon_interval):
                        # get mask
                        msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

                        if len(yp_mon[msk_mon]) != 0:
                            Zmax_mon = max(yp_mon[msk_mon])

                            Zmax_mon_list[i].append(Zmax_mon)
                        elif len(yp_mon) != 0:
                            Zmax_mon_list[i].append(0)
                        else:
                            ic("skipped, yp_mon = [], raise exception")
                            raise Exception()

                    processed_keys_mon.append(key)
                except:
                    ic("skipped, yp_mon = []")
                    # for i, z_bound in enumerate(mon_interval):
                    #     Zmax_mon_list[i].append(-1)

            # print("2g", Zmax_mon_list)

            return Zmax_mon_list

        def get_Zmax_T(dip_interval=None):
            if dip_interval is None:
                dip_interval = [0.0, 2e10]

            for key, value in d.items():
                try:
                    print(f"Processing for Cavity {key}")
                    abci_data_dip = ABCIData(abci_data_dir, f"Cavity{key}", 1)

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
                            ic("skipped, yp_dip = [], raise exception")
                            raise Exception()

                    processed_keys_dip.append(key)
                except:
                    ic("skipped, yp_dip = []")
                    # for i, z_bound in enumerate(dip_interval):
                    #     Zmax_dip_list[i].append(-1)

            return Zmax_dip_list

        def all(mon_interval, dip_interval):
            for key, value in d.items():
                print(f"Processing for Cavity {key}")
                abci_data_long = ABCIData(abci_data_dir, f"Cavity{key}_", 0)
                abci_data_trans = ABCIData(abci_data_dir, f"Cavity{key}_", 1)

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
                    continue

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

                for i, z_bound in enumerate(mon_interval):
                    # get mask
                    msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

                    if len(yp_mon[msk_mon]) != 0:
                        Zmax_mon = max(yp_mon[msk_mon])
                        xmax_mon = xp_mon[np.where(yp_mon == Zmax_mon)][0]

                        Zmax_mon_list[i].append(Zmax_mon)
                        xmax_mon_list[i].append(xmax_mon)
                    elif len(yp_mon) != 0:
                        Zmax_mon_list[i].append(0.0)
                        xmax_mon_list[i].append(0.0)
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
                k_loss_array_longitudinal.append(k_loss_HOM)
                k_loss_array_transverse.append(k_loss_trans)

        obj = self.tuneUI.ccb_Populate_Objectives.currentText().split(', ')
        ZL, ZT = [], []
        df_ZL, df_ZT = pd.DataFrame(), pd.DataFrame()
        # print("here here here")
        for i, o in enumerate(obj):
            # print("print ptint pting")
            if "ZL" in o:
                # print("1")
                freq_range = self.process_interval(self.text_to_list(self.tuneUI.tw_Objectives.cellWidget(i, 2).text()))
                for i in range(len(freq_range)):
                    Zmax_mon_list.append([])
                    xmax_mon_list.append([])
                    df_ZL[f"{o} [max({freq_range[i][0]}<f<{freq_range[i][1]})]"] = 0

                ZL = get_Zmax_L(freq_range)

            elif "ZT" in o:
                freq_range = self.process_interval(self.text_to_list(self.tuneUI.tw_Objectives.cellWidget(i, 2).text()))

                for i in range(len(freq_range)):
                    Zmax_dip_list.append([])
                    xmax_dip_list.append([])
                    df_ZT[f"{o} [max({freq_range[i][0]}<f<{freq_range[i][1]})]"] = 0

                ZT = get_Zmax_T(freq_range)

            elif o[1] == "k_loss":
                pass
            elif o[1] == "k_kick":
                pass

        # create dataframes from list
        ic(np.array(ZL).T)
        ic(np.array(ZT).T)
        ic(processed_keys_mon, processed_keys_dip)
        df_ZL.loc[:, :] = np.array(ZL).T
        df_ZT.loc[:, :] = np.array(ZT).T
        ic(df_ZL)
        df_ZL['key'] = processed_keys_mon
        df_ZT['key'] = processed_keys_dip

        ic(df_ZL, df_ZT)
        processed_keys = list(set(processed_keys_mon) & set(processed_keys_dip))

        # ZL, ZT = np.array(ZL).T, np.array(ZT).T

        if len(ZL) != 0 and len(ZT) != 0:
            df_wake = df_ZL.merge(df_ZT, on='key', how='inner')
            # obj_result = np.hstack((ZL, ZT))
        elif len(ZL) != 0:
            df_wake = df_ZL
            # obj_result = ZL
        else:
            df_wake = df_ZT
            # obj_result = ZT

        return df_wake, processed_keys

    def crossover(self, df, generation, f):  # , rq, grq
        ic(df)
        elites = {}
        for i, o in enumerate(self.objectives):

            if self.tuneUI.cb_UQ.isChecked():
                if o[0] == "min":
                    elites[f'E[{o[1]}] + 6*std[{o[1]}]'] = df.sort_values(f'E[{o[1]}] + 6*std[{o[1]}]')
                elif o[0] == "max":
                    elites[f'E[{o[1]}] - 6*std[{o[1]}]'] = df.sort_values(f'E[{o[1]}] - 6*std[{o[1]}]', ascending=False)
                elif o[0] == "equal":
                    elites[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] = df.sort_values(fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]')
            else:
                if o[0] == "min":
                    elites[f'{o[1]}'] = df.sort_values(f'{o[1]}')
                elif o[0] == "max":
                    elites[f'{o[1]}'] = df.sort_values(f'{o[1]}', ascending=False)
                elif o[0] == "equal":
                    continue
        ic(elites)
        obj_dict = {}
        for o in self.objectives:
            if self.tuneUI.cb_UQ.isChecked():
                if o[0] == 'min':
                    obj_dict[fr'E[{o[1]}] + 6*std[{o[1]}]'] = elites[fr'E[{o[1]}] + 6*std[{o[1]}]']
                elif o[0] == 'max':
                    obj_dict[fr'E[{o[1]}] - 6*std[{o[1]}]'] = elites[fr'E[{o[1]}] - 6*std[{o[1]}]']
                else:
                    obj_dict[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] = elites[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]']
            else:
                if o[0] != 'equal':
                    obj_dict[o[1]] = elites[o[1]]

        obj = {}
        for key, o in obj_dict.items():
            obj[key] = o.reset_index(drop=True)

        # e, b, rq = obj_list
        # e = e.reset_index(drop=True)
        # b = b.reset_index(drop=True)
        # rq = rq.reset_index(drop=True)

        # naming convention G<generation number>_C<cavity number>_<type>
        # type refers to mutation M or crossover C
        df_co = pd.DataFrame(columns=["key", 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"])

        # select only best characteristics
        A_inf = self.tuneUI.tw_Parameters.cellWidget(0, 3).currentText().split(', ')
        B_inf = self.tuneUI.tw_Parameters.cellWidget(1, 3).currentText().split(', ')
        a_inf = self.tuneUI.tw_Parameters.cellWidget(2, 3).currentText().split(', ')
        b_inf = self.tuneUI.tw_Parameters.cellWidget(3, 3).currentText().split(', ')
        Ri_inf = self.tuneUI.tw_Parameters.cellWidget(4, 3).currentText().split(', ')
        L_inf = self.tuneUI.tw_Parameters.cellWidget(5, 3).currentText().split(', ')
        Req_inf = self.tuneUI.tw_Parameters.cellWidget(6, 3).currentText().split(', ')

        ic(elites)
        inf_dict = {"A": A_inf, "B": B_inf, "a": a_inf, "b": b_inf, "Ri": Ri_inf, "L": L_inf, "Req": Req_inf}
        for key, influence in inf_dict.items():
            if influence == [''] or influence == ['All']:
                if self.tuneUI.cb_UQ.isChecked():
                    ll = []
                    for o in self.objectives:
                        if o[0] == 'min':
                            ll.append(fr'E[{o[1]}] + 6*std[{o[1]}]')
                        elif o[0] == 'max':
                            ll.append(fr'E[{o[1]}] - 6*std[{o[1]}]')
                        else:
                            ll.append(fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]')
                    inf_dict[key] = ll
                else:
                    inf_dict[key] = [o[1] for o in self.objectives if o[0] != 'equal']

        ic(elites.keys())
        n_elites_to_cross = self.tuneUI.sb_N_Elites_To_Cross.value()
        ic(n_elites_to_cross, len(elites))
        for i in range(f):
            # (<obj>[<rank>][<variable>] -> (b[c[1]][0]

            df_co.loc[i] = [f"G{generation}_C{i}_CO",
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["A"] for key in inf_dict["A"]]) / len(inf_dict["A"]),  # A
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["B"] for key in inf_dict["B"]]) / len(inf_dict["B"]),  # B
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["a"] for key in inf_dict["a"]]) / len(inf_dict["a"]),  # a
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["b"] for key in inf_dict["b"]]) / len(inf_dict["b"]),  # b
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["Ri"] for key in inf_dict["Ri"]]) / len(inf_dict["Ri"]),  # Ri
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["L"] for key in inf_dict["L"]]) / len(inf_dict["L"]),  # L
                            sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0]-1)]["Req"] for key in inf_dict["Req"]]) / len(inf_dict["Req"]),
                            0,
                            0
                            ]
        return df_co

    def mutation(self, df, n, f):
        # get range from table
        tw = self.tuneUI.tw_Parameters
        for i in range(self.tuneUI.tw_Parameters.rowCount()):
            self.sbd[i] = [tw.cellWidget(i, 1).value(), tw.cellWidget(i, 2).value()]

        # get list based on mutation length
        if df.shape[0] < f:
            ic("It is here")
            ml = np.arange(df.shape[0])
        else:
            ml = np.arange(f)

        ic(ml)
        df_ng_mut = pd.DataFrame(columns=['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"])
        if self.sbd[0][0] == self.sbd[0][1]:
            df_ng_mut.loc[:, 'A'] = df.loc[ml, "A"]
        else:
            df_ng_mut.loc[:, 'A'] = df.loc[ml, "A"] * random.uniform(0.85, 1.5)

        if self.sbd[1][0] == self.sbd[1][1]:
            df_ng_mut.loc[:, 'B'] = df.loc[ml, "B"]
        else:
            df_ng_mut.loc[:, 'B'] = df.loc[ml, "B"] * random.uniform(0.85, 1.5)

        if self.sbd[2][0] == self.sbd[2][1]:
            df_ng_mut.loc[:, 'a'] = df.loc[ml, "a"]
        else:
            df_ng_mut.loc[:, 'a'] = df.loc[ml, "a"] * random.uniform(0.85, 1.5)

        if self.sbd[3][0] == self.sbd[3][1]:
            df_ng_mut.loc[:, 'b'] = df.loc[ml, "b"]
        else:
            df_ng_mut.loc[:, 'b'] = df.loc[ml, "b"] * random.uniform(0.85, 1.5)

        if self.sbd[4][0] == self.sbd[4][1]:
            df_ng_mut.loc[:, 'Ri'] = df.loc[ml, "Ri"]
        else:
            df_ng_mut.loc[:, 'Ri'] = df.loc[ml, "Ri"] * random.uniform(0.85, 1.5)

        if self.sbd[5][0] == self.sbd[5][1]:
            df_ng_mut.loc[:, 'L'] = df.loc[ml, "L"]
        else:
            df_ng_mut.loc[:, 'L'] = df.loc[ml, "L"] * random.uniform(0.85, 1.5)

        if self.sbd[6][0] == self.sbd[6][1]:
            df_ng_mut.loc[:, 'Req'] = df.loc[ml, "Req"]
        else:
            df_ng_mut.loc[:, 'Req'] = df.loc[ml, "Req"] * random.uniform(0.85, 1.5)

        df_ng_mut.loc[:, ["alpha_i", "alpha_o"]] = df.loc[ml, ["alpha_i", "alpha_o"]]

        key1, key2 = [], []
        for i in range(len(df_ng_mut)):
            key1.append(f"G{n}_C{i}_M")

        df_ng_mut.loc[:, 'key'] = key1

        return df_ng_mut

    def chaos(self, f, n):
        data = {'key': [f"G{n}_C{i}_CH" for i in range(f)],
                'A': random.sample(list(np.linspace(self.sbd[0][0], self.sbd[0][1], f * 10)), f),
                'B': random.sample(list(np.linspace(self.sbd[1][0], self.sbd[1][1], f * 10)), f),
                'a': random.sample(list(np.linspace(self.sbd[2][0], self.sbd[2][1], f * 10)), f),
                'b': random.sample(list(np.linspace(self.sbd[3][0], self.sbd[3][1], f * 10)), f),
                'Ri': random.sample(list(np.linspace(self.sbd[4][0], self.sbd[4][1], f * 10)), f),
                'L': random.sample(list(np.linspace(self.sbd[5][0], self.sbd[5][1], f * 10)), f),
                'Req': random.sample(list(np.linspace(self.sbd[6][0], self.sbd[6][1], f * 10)), f),
                'alpha_i': np.zeros(f),
                'alpha_o': np.zeros(f)}

        df = pd.DataFrame.from_dict(data)
        return df

    def remove_duplicate_values(self, d):
        temp = []
        res = dict()
        for key, val in d.items():
            if val not in temp:
                temp.append(val)
                res[key] = val
        return res

    def proof_filename(self, dirc):
        # check if extension is included
        if dirc.split('.')[-1] != 'json':
            dirc = f'{dirc}.json'

        return dirc

    def recursive_save(self, df, filename, pareto_index):
        styler = self.color_pareto(df, self.poc)
        try:
            styler.to_excel(filename)
            # df.to_excel(filename)
        except PermissionError:
            filename = filename.split('.xlsx')[0]
            filename = fr'{filename}_1.xlsx'
            self.recursive_save(df, filename, pareto_index)

    def populate_objectives(self, f, ccb, tw, d):
        tw.setRowCount(len(f))  # and one row in the table

        for i, x in enumerate(f):
            label = QLabel(x)
            tw.setCellWidget(i, 0, label)

            cb = QComboBox()
            cb.addItems(['min', 'max', 'equal'])
            tw.setCellWidget(i, 1, cb)

            # check for impedance
            if x == "ZL":
                le = QLineEdit()
                le.setPlaceholderText("Enter interval(s) to evaluate max impedance")
                le.setText("1, 2, 5")
                tw.setCellWidget(i, 2, le)

            if x == "ZT":
                le = QLineEdit()
                le.setPlaceholderText("Enter interval(s) to evaluate max impedance")
                le.setText("0.7, 1.3, 2.25, 5")
                tw.setCellWidget(i, 2, le)

    def populate_constraints(self, f, ccb, tw, d):
        tw.setRowCount(len(f))  # and one row in the table

        for i, x in enumerate(f):
            label = QLabel(x)
            tw.setCellWidget(i, 0, label)

            dsb = QDoubleSpinBox()
            dsb.setMinimum(1)
            dsb.setMaximum(10000000000)
            tw.setCellWidget(i, 1, dsb)

            dsb2 = QDoubleSpinBox()
            dsb2.setMinimum(1)
            dsb2.setMaximum(10000000000)
            tw.setCellWidget(i, 2, dsb2)

    def populate_parameters(self):
        dd = {'A': [20, 80], 'B': [20, 80], 'a': [10, 60], 'b': [10, 60], 'Ri': [60, 85], 'L': [93.5, 93.5],
              'Req': [170, 170], }
        # dd = {'A': [10, 60], 'B': [10, 60], 'a': [5, 40], 'b': [5, 40], 'Ri': [35, 35], 'L': [57.7, 57.7],
        #       'Req': [103, 103], }
        tw = self.tuneUI.tw_Parameters

        tw.setRowCount(len(dd))  # and one row in the table
        i = 0
        self.influence_ccb = {}
        for key, val in dd.items():
            label = QLabel(key)
            tw.setCellWidget(i, 0, label)

            dsb = QDoubleSpinBox()
            dsb.setMinimum(1)
            dsb.setMaximum(10000000000)
            dsb.setValue(val[0])
            tw.setCellWidget(i, 1, dsb)

            dsb2 = QDoubleSpinBox()
            dsb2.setMinimum(1)
            dsb2.setMaximum(10000000000)
            dsb2.setValue(val[1])
            tw.setCellWidget(i, 2, dsb2)

            ccb = QCheckableComboBox()
            ccb.addItem("All")

            self.influence_ccb[key] = ccb
            tw.setCellWidget(i, 3, ccb)

            cb = QCheckableComboBox()
            cb.addItem("No")
            cb.addItem("Yes")
            tw.setCellWidget(i, 4, cb)

            i += 1

    def update_influence_ccb(self):
        for ccb in self.influence_ccb.values():
            ccb.clear()
            ccb.addItem('All')
            ccb.addItems(self.tuneUI.ccb_Populate_Objectives.currentText().split(', '))

    def pareto_front(self, df):

        # datapoints = np.array([reverse_list(x), reverse_list(y), reverse_list(z)])
        # reverse list or not based on objective goal: minimize or maximize
        # datapoints = [self.negate_list(df.loc[:, o[1]], o[0]) for o in self.objectives]

        if self.tuneUI.cb_UQ.isChecked():
            obj = []
            for o in self.objectives:
                if o[0] == 'min':
                    obj.append(fr'E[{o[1]}] + 6*std[{o[1]}]')
                elif o[0] == 'max':
                    obj.append(fr'E[{o[1]}] - 6*std[{o[1]}]')
                elif o[0] == 'equal':
                    obj.append(fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]')

            # ic(obj)
            datapoints = df.loc[:, obj]
        else:
            datapoints = df.loc[:, [o[1] for o in self.objectives]]

        # ic(datapoints)
        for o in self.objectives:
            if o[0] == 'min':
                if self.tuneUI.cb_UQ.isChecked():
                    datapoints[fr'E[{o[1]}] + 6*std[{o[1]}]'] = datapoints[fr'E[{o[1]}] + 6*std[{o[1]}]'] * (-1)
                else:
                    datapoints[o[1]] = datapoints[o[1]] * (-1)
            elif o[0] == "equal":
                if self.tuneUI.cb_UQ.isChecked():
                    datapoints[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] = datapoints[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] * (-1)
        # ic(datapoints)
        # convert datapoints to numpy array

        pareto = oapackage.ParetoDoubleLong()
        # ic(datapoints)
        for ii in range(0, datapoints.shape[0]):
            w = oapackage.doubleVector(tuple(datapoints.iloc[ii].values))
            pareto.addvalue(w, ii)
        pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs
        # ic(lst)
        self.poc = len(lst)
        # ic(self.poc)
        reorder_idx = list(lst) + [i for i in range(len(df)) if i not in lst]
        # ic(reorder_idx)
        # ic(lst)
        # optimal_datapoints = df.loc[lst, :]
        # ic(optimal_datapoints)
        # ic([optimal_datapoints.loc[i, :] for i in range(len(lst))])

        # return [optimal_datapoints[i, :] for i in range(datapoints.shape[0])]
        return reorder_idx

    def check_input(self, s):
        # s = "range(16, 23, 10)"
        # s = "randrange(16, 23, 10)"
        # s = "[16, 23, 10]"
        # s = 1, 2, 3
        # s = 2

        if "r" in s and "rr" not in s:
            s = s.replace('r', '')
            try:
                l = eval(s)
                return np.linspace(l[0], l[1], l[2])
            except:
                print("Please check inputs.")
        elif "rr" in s:
            s = s.replace('rr', '')
            try:
                l = eval(s)
                ll = np.random.uniform(l[0], l[1], l[2])
                return ll
            except:
                print("Please check inputs.")
        else:
            try:
                ll = eval(s)
                if isinstance(ll, int) or isinstance(ll, float):
                    ll = [ll]
                return ll
            except:
                print("Please check inputs.")

        return 1

    @staticmethod
    def negate_list(l, arg):
        if arg == 'max':
            return l  # to find the pareto maxima
        else:
            return [-x for x in l]  # to find the pareto minima

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
    def run_sequential(pseudo_shape_space_proc, resume, p, bc, parentDir, projectDir, filename, tuner_option,
                       tune_variable, iter_set, cell_type, progress_list, convergence_list, save_last, n_cells):

        tuner.tune(pseudo_shape_space_proc, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                   tuner_option=tuner_option, tune_variable=tune_variable, iter_set=iter_set,
                   cell_type=cell_type, progress_list=progress_list, convergence_list=convergence_list,
                   save_last=save_last, n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    @staticmethod
    def text_to_list(l):
        if l == '':
            return None
        else:
            l = ast.literal_eval(l)
            if isinstance(l, int) or isinstance(l, float):
                return [l, 2e10]
            else:
                return list(l)

    @staticmethod
    def process_interval(interval_list):
        interval = []
        for i in range(len(interval_list) - 1):
            interval.append([interval_list[i], interval_list[i + 1]])

        return interval

    # def continue_check(self):
    #     path = f'{self.main_control.projectDir}/Cavities/pseudo_{self.proof_filename(self.tuneUI.le_Generated_Shape_Space_Name.text())}'
    #     print(path)
    #     if os.path.exists(path):
    #         msg = QMessageBox()
    #         msg.setWindowTitle("Resume Simulation")
    #         msg.setText("The pseudo shape space file already exist. Would you love to update or overwrite its content?")
    #         msg.setIcon(QMessageBox.Question)
    #         msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
    #
    #         buttonY = msg.button(QMessageBox.Yes)
    #         buttonY.setText('Update')
    #         buttonN = msg.button(QMessageBox.No)
    #         buttonN.setText('Overwrite')
    #
    #         msg.setDefaultButton(buttonY)
    #
    #         msg.buttonClicked.connect(self.button_clicked)
    #
    #         x = msg.exec_()
    #
    #         if x == 16384:
    #             print("Done checking path. Yes selected.")
    #             return "Yes"
    #         elif x == 65536:
    #             print("Done checking path. No selected.")
    #             return "No"
    #         elif x == 4194304:
    #             print("Done checking path. Cancel selected.")
    #             return "Cancel"
    #     else:
    #         print("Path does not yet exist. Creating...")
    #         return "Yes"

    @staticmethod
    def color_pareto(df, no_pareto_optimal):
        def color(row):
            # if row.isnull().values.any():
            if row[0] in df['key'].tolist()[0:no_pareto_optimal]:
                return ['background-color: #6bbcd1'] * len(row)
            return [''] * len(row)

        # Save Styler Object for Later
        styler = df.style
        # Apply Styles (This can be chained or on separate lines)
        styler.apply(color, axis=1)
        # Export the styler to excel
        return styler

    def show_hide_mid_cell_params(self):
        if self.tuneUI.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':
            self.tuneUI.w_Mid_Cell_Opt.show()
        else:
            self.tuneUI.w_Mid_Cell_Opt.hide()


def run_sequential_wakefield(n_cells, n_modules, processor_shape_space,
                             MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
                             DDR_SIG=0.1, DDZ_SIG=0.1,
                             parentDir=None, projectDir=None, progress_list=None,
                             WG_M=None, marker=''):
    progress = 0
    # get length of processor
    total_no_of_shapes = len(list(processor_shape_space.keys()))
    for key, shape in processor_shape_space.items():
        skip = False
        if os.path.exists(fr'{projectDir}\SimulationData\ABCI\Cavity{key}'):
            skip = True

        start_time = time.time()
        if not skip:
            # run abci code
            # run both polarizations if MROT == 2
            try:
                if MROT == 2:
                    for m in range(2):
                        abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                         fid=key, MROT=m, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                         projectDir=projectDir,
                                         WG_M='', marker='')

                else:
                    abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                     fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                     DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir,
                                     WG_M='', marker='')
            except:
                if MROT == 2:
                    for m in range(2):
                        abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                         fid=key, MROT=m, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                         projectDir=projectDir,
                                         WG_M='', marker='')
                else:
                    abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                     fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                     DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir,
                                     WG_M='', marker='')

        print_(f'Cavity {key}. Time: {time.time() - start_time}')

            # update progress
        progress_list.append((progress + 1) / total_no_of_shapes)


# def stroud(p):
#     # Stroud-3 method
#     #
#     # Input parameters:
#     #  p   number of dimensions
#     # Output parameters:
#     #  nodes   nodes of quadrature rule in [0,1]^p (column-wise)
#     #
#
#     nodes = np.zeros((p, 2 * p))
#     coeff = np.pi / p
#     fac = np.sqrt(2 / 3)
#
#     for i in range(2 * p):
#         for r in range(int(np.floor(0.5 * p))):
#             k = 2 * r
#             nodes[k, i] = fac * np.cos((k + 1) * (i + 1) * coeff)
#             nodes[k + 1, i] = fac * np.sin((k + 1) * (i + 1) * coeff)
#
#         if 0.5 * p != np.floor(0.5 * p):
#             nodes[-1, i] = ((-1) ** (i + 1)) / np.sqrt(3)
#
#     # transform nodes from [-1,+1]^p to [0,1]^p
#     nodes = 0.5 * nodes + 0.5
#
#     return nodes
#
#
# def quad_stroud3(rdim, degree):
#     # data for Stroud-3 quadrature in [0,1]^k
#     # nodes and weights
#     nodes = stroud(rdim)
#     nodestr = 2. * nodes - 1.
#     weights = (1 / (2 * rdim)) * np.ones((2 * rdim, 1))
#
#     # evaluation of Legendre polynomials
#     bpoly = np.zeros((degree + 1, rdim, 2 * rdim))
#     for l in range(rdim):
#         for j in range(2 * rdim):
#             bpoly[0, l, j] = 1
#             bpoly[1, l, j] = nodestr[l, j]
#             for i in range(1, degree):
#                 bpoly[i + 1, l, j] = ((2 * (i + 1) - 1) * nodestr[l, j] * bpoly[i, l, j] - i * bpoly[
#                     i - 1, l, j]) / (i + 1)
#
#     # standardisation of Legendre polynomials
#     for i in range(1, degree + 1):
#         bpoly[i, :, :] = bpoly[i, :, :] * np.sqrt(2 * (i + 1) - 1)
#
#     return nodes, weights, bpoly


# def weighted_mean_obj(tab_var, weights):
#         rows_sims_no, cols = np.shape(tab_var)
#         no_weights, dummy = np.shape(weights)  # z funckji quadr_stroud wekt columnowy
#
#         if rows_sims_no == no_weights:
#             expe = np.zeros((cols, 1))
#             outvar = np.zeros((cols, 1))
#             for i in range(cols):
#                 expe[i, 0] = np.dot(tab_var[:, i], weights)
#                 outvar[i, 0] = np.dot(tab_var[:, i]**2, weights)
#             stdDev = np.sqrt(abs(outvar - expe**2))
#         else:
#             expe = 0
#             stdDev = 0
#             ic('Cols_sims_no != No_weights')
#
#         return list(expe.T[0]), list(stdDev.T[0])


def get_objectives_value(d, obj, norm_length, n_cells):
    Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
    Freq = d['FREQUENCY'][n_cells - 1]
    E_stored = d['STORED ENERGY'][n_cells - 1]
    Rsh = d['SHUNT IMPEDANCE'][n_cells - 1]  # MOhm
    Q = d['QUALITY FACTOR'][n_cells - 1]
    Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
    Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
    # Vacc = dict['ACCELERATION'][n_cells - 1]
    Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells - 1]  # MV/m
    r_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm
    G = 0.00948*Q*(Freq/1300)
    GR_Q = G * 2 * r_Q

    Vacc = np.sqrt(
        2 * r_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6  # factor of 2, remember circuit and accelerator definition
    # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
    Eacc = Vacc / (n_cells*norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
    Epk_Eacc = Epk / Eacc
    Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc

    d = {
        "Req": Req,
        "freq": Freq,
        "Q": Q,
        "E": E_stored,
        "R/Q": 2 * r_Q,
        "Epk/Eacc": Epk_Eacc,
        "Bpk/Eacc": Bpk_Eacc,
        "G": G,
        "GR/Q": GR_Q
    }

    objective = []
    # append freq and Req
    tune_result = [Req, Freq]


    # append objective functions
    for o in obj:
        if o[1] in d.keys():
            objective.append(d[o[1]])

    return objective, tune_result


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
        # print("2a")
        if mon_interval is None:
            mon_interval = [0.0, 2e10]

        print(f"Processing for Cavity {key}")
        try:
            abci_data_mon = ABCIData(abci_data_dir, f"Cavity{key}", 0)

            # get longitudinal and transverse impedance plot data
            xr_mon, yr_mon, _ = abci_data_mon.get_data('Real Part of Longitudinal Impedance')
            xi_mon, yi_mon, _ = abci_data_mon.get_data('Imaginary Part of Longitudinal Impedance')

            # print("2d")
            # Zmax
            if mon_interval is None:
                mon_interval = [[0.0, 10]]

            # calculate magnitude
            ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]

            # print("2e")
            # get peaks
            peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
            xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

            # print("2f", mon_interval)
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

        # print("2g", Zmax_mon_list)

        return Zmax_mon_list

    def get_Zmax_T(dip_interval=None):
        if dip_interval is None:
            dip_interval = [0.0, 2e10]

        try:
            print(f"Processing for Cavity {key}")
            abci_data_dip = ABCIData(abci_data_dir, f"Cavity{key}", 1)

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

    def all(mon_interval, dip_interval):
        print(f"Processing for Cavity {key}")
        abci_data_long = ABCIData(abci_data_dir, f"Cavity{key}_", 0)
        abci_data_trans = ABCIData(abci_data_dir, f"Cavity{key}_", 1)

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

        for i, z_bound in enumerate(mon_interval):
            # get mask
            msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

            if len(yp_mon[msk_mon]) != 0:
                Zmax_mon = max(yp_mon[msk_mon])
                xmax_mon = xp_mon[np.where(yp_mon == Zmax_mon)][0]

                Zmax_mon_list[i].append(Zmax_mon)
                xmax_mon_list[i].append(xmax_mon)
            elif len(yp_mon) != 0:
                Zmax_mon_list[i].append(0.0)
                xmax_mon_list[i].append(0.0)
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
    # print("here here here")
    for i, o in enumerate(obj):
        # print("print ptint pting")
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

    # ic("about to evaluate ZT", freq_range_ZT)
    if freq_range_ZT:
        for i in range(len(freq_range_ZT)):
            Zmax_dip_list.append([])
            xmax_dip_list.append([])

        ZT = get_Zmax_T(freq_range_ZT)

    # ic("Donbr rvasluating impedances")
    ZL, ZT = np.array(ZL).T, np.array(ZT).T
    # ic(ZL, ZT)

    if ZL.size != 0 and ZT.size != 0:
        obj_result = np.hstack((ZL, ZT))
    elif ZL.size != 0:
        obj_result = ZL
    else:
        obj_result = ZT

    # ic(obj_result)

    return list(obj_result[0])


def process_interval(interval_list):
        interval = []
        for i in range(len(interval_list) - 1):
            interval.append([interval_list[i], interval_list[i + 1]])

        return interval
