import os.path
import subprocess
import time
import multiprocessing as mp
from threading import Thread
from pathlib import Path

import pandas as pd
from psutil import NoSuchProcess
from scipy.stats import qmc

from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from ui_files.eigenmode import Ui_Eigenmode
from utils.shared_classes import *
from utils.shared_functions import *
from analysis_modules.eigenmode.NGSolve.eigen_ngsolve import NGSolveMEVP
from analysis_modules.eigenmode.customEig.run_field_solver import Model

slans_geom = SLANSGeometry()
ngsolve_mevp = NGSolveMEVP()
fr = FileReader()

file_color = 'green'
DEBUG = True


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class EigenmodeControl:
    def __init__(self, parent):
        self.pol = None
        self.animation = None
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

        self.progress_bar = self.ui.pb_Progress_Bar

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # Get plot object
        self.geometry_view = self.win.geometryview_widget
        self.plot = self.geometry_view.plot

        # geometry input ui
        self.geo_control = self.main_control.geometryinput_widget
        self.geo_ui = self.geo_control.ui

        # get logger
        self.log = self.main_control.log

        self.initUI()
        self.signals()
        self.exe_control()
        self.filename = None  # placeholder, made a booboo copying the end routine

        # instantiate geometry
        self.slans_geom = SLANSGeometry()

        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False

        # ui effects
        # self.ui_effects()

    def initUI(self):
        self.ui.w_EXE.setVisible(False)
        self.ui.w_UQ.setVisible(False)

        # create pause and resume1 icons to avoid creating them over and over again
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
        self.progress_bar.setValue(0)
        self.progress_bar.hide()

        self.ui.pb_Refresh.clicked.connect(lambda: self.write_to_text_edit())

        self.change_mesh_input()

    def signals(self):
        # run eigenmode solver
        self.ui.pb_Run.clicked.connect(lambda: self.run_slans())

        # cancel
        self.ui.pb_Cancel.clicked.connect(lambda: self.cancel())
        self.ui.pb_Pause_Resume.clicked.connect(
            lambda: self.pause() if self.process_state == 'running' else self.resume())
        self.ui.cb_LBC.currentTextChanged.connect(lambda: self.geo_control.draw_shape_from_shape_space())
        self.ui.cb_RBC.currentTextChanged.connect(lambda: self.geo_control.draw_shape_from_shape_space())

        self.ui.cb_Eigenproblem_Solver.currentTextChanged.connect(lambda: self.change_mesh_input())

    def run_slans(self):
        # get analysis parameters
        n_cells = text_to_list(self.geo_ui.le_N_Cells.text())
        n_modules = self.geo_ui.sb_N_Modules.value()
        f_shift = float(self.ui.le_Freq_Shift.text())
        n_modes = float(self.ui.le_No_Of_Modes.text())
        self.pol = self.ui.cb_Polarization_SLANS.currentText()

        # boundary conditions
        lbc = self.ui.cb_LBC.currentIndex() + 1
        rbc = self.ui.cb_RBC.currentIndex() + 1
        bc = 10 * lbc + rbc

        # solver
        select_solver = self.ui.cb_Eigenproblem_Solver.currentText()

        if select_solver.lower() == 'slans':
            mesh_args = [self.ui.sb_Mesh_Jxy_Cell.value(), self.ui.sb_Mesh_Jx_BP.value(), self.ui.sb_Mesh_Jy_BP.value(),
                         self.ui.sb_Max_Iteration.value()]
        else:
            mesh_args = [self.ui.sb_Max_Cells_Per_Wavelength.value(),
                         self.ui.sb_Max_Iteration.value()]

        proc_count = self.ui.sb_No_Of_Processors_SLANS.value()

        # uq
        if self.ui.cb_UQ.isChecked():
            UQ = True
        else:
            UQ = False

        # get geometric parameters
        self.shape_space = get_geometric_parameters(self.geo_control, 'SLANS',
                                                    text_to_list(self.geo_ui.le_Scale.text()))

        # split shape_space for different processes/ MPI share process by rank
        keys = list(self.shape_space.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

        # get the total number of simulations to be run
        num_sims = shape_space_len * len(n_cells)
        self.progress_bar.setMaximum(num_sims)
        self.progress_bar.show()

        # progress list
        manager = mp.Manager()
        self.progress_list = manager.list()
        self.progress_list.append(0)

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

            service = mp.Process(target=self.run_sequential, args=(
                n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, self.pol, self.main_control.parentDir,
                self.main_control.projectDir, self.progress_list, '',
                UQ, select_solver, mesh_args))

            service.start()

            self.processes.append(psutil.Process(service.pid))
            self.processes_id.append(service.pid)

            # except Exception as e:
            #     self.log.error(fr"Exception in run_MP:: {e}")
            #     # print_("Exception in run_MP::", e)

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

    def change_mesh_input(self):
        if self.ui.cb_Eigenproblem_Solver.currentText().lower() == 'slans':
            self.ui.w_SLANS_Mesh.show()
            self.ui.w_NGSolve_Mesh.hide()

            # change max iteration
            self.ui.sb_Max_Iteration.setValue(50)
        else:
            self.ui.w_SLANS_Mesh.hide()
            self.ui.w_NGSolve_Mesh.show()
            # change max iteration
            self.ui.sb_Max_Iteration.setValue(20)

    def write_to_text_edit(self):
        """
        Writes information about the eigenmode simulation to QTextEdit
        Returns
        -------

        """
        try:
            # reset QTextEdit
            self.ui.textEdit.clear()
            errorFormat = '&nbsp;&nbsp;&nbsp;&nbsp;<span style="color:red;">{}</span>'
            warningFormat = '&nbsp;&nbsp;&nbsp;&nbsp;<span style="color:orange;">{}</span>'
            validFormat = '&nbsp;&nbsp;&nbsp;&nbsp;<span style="color:green;">{}</span>'
            for cav in self.shape_space.keys():
                # print out results form svl file to QTextEdit
                lbc = self.ui.cb_LBC.currentIndex() + 1
                rbc = self.ui.cb_RBC.currentIndex() + 1

                # check if file exists
                if self.pol == 'Monopole':
                    filename = fr'{self.main_control.projectDir}/SimulationData/SLANS/' \
                               fr'{cav}_n{self.geo_ui.le_N_Cells.text()}/monopole/cavity_{lbc}{rbc}.svl'
                else:
                    filename = fr'{self.main_control.projectDir}/SimulationData/SLANS/' \
                               fr'{cav}_n{self.geo_ui.le_N_Cells.text()}/dipole/cavity_{lbc}{rbc}_2.sv2'

                if os.path.exists(filename):
                    # Open the file for reading
                    with open(filename, 'r') as file:
                        # Read all lines
                        lines = file.readlines()

                    # Filter lines containing the keyword 'ACCURACY'
                    if self.pol == 'Monopole':
                        filtered_lines = [line.strip() for line in lines if ('ACCURACY' in line or 'FREQUENCY' in line)]
                    else:
                        filtered_lines = [line.strip() for line in lines if ('accuracy' in line or 'Frequency' in line)]
                    accuracy = filtered_lines[::2]
                    frequency = filtered_lines[1::2]

                    # cavity name
                    self.ui.textEdit.append(fr'<b>{cav}</b>')
                    # Print the filtered lines
                    for acc, freq in zip(accuracy, frequency):
                        # check if accuracy is below threshold
                        a = float(acc.split(' ')[-1])
                        f = freq.split(' ')
                        if a > 1e-3:
                            self.ui.textEdit.append(errorFormat.format(
                                f'\tFreq.: {f[-2].strip()} {f[-1].strip()}, Accuracy: {a}'))
                        elif 1e-5 < a <= 1e-3:
                            self.ui.textEdit.append(warningFormat.format(
                                f'\tFreq.: {f[-2].strip()} {f[-1].strip()}, Accuracy: {a}'))
                        else:
                            self.ui.textEdit.append(validFormat.format(
                                f'\tFreq.: {f[-2].strip()} {f[-1].strip()}, Accuracy: {a}'))

        except Exception as e:
            ic('Exception in slans, run_pauseroutine: ', e)

    def run_pause_resume_stop_routine(self):
        """
        Controls the running, pausing and stopping of simulation
        Returns
        -------

        """
        if self.process_state == 'none':
            # change pause/resume icon to pause icon
            self.ui.pb_Pause_Resume.setIcon(self.pause_icon)

            # disable pause/resume and cancel buttons
            self.ui.pb_Pause_Resume.setEnabled(False)
            self.ui.pb_Cancel.setEnabled(False)

            # enable run button in case it was disabled
            self.ui.pb_Run.setEnabled(True)

            # click message buttion
            self.ui.pb_Refresh.click()

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
            # except:
            #     pass

        self.cancel()

    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

        if val == self.progress_bar.maximum() or not self.show_progress_bar:
            # reset progress bar
            self.progress_bar.setValue(0)
            self.progress_bar.hide()

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

    def exe_control(self):
        """
        Control SLANS executable files
        Returns
        -------

        """
        # Slans
        self.ui.pb_Genmsh.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\genmesh2.exe")))
        self.ui.pb_Sl.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\Sl.exe")))
        self.ui.pb_Slansc.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\slansc.exe")))
        self.ui.pb_Slansm.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\slansm.exe")))
        self.ui.pb_Slanss.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\slanss.exe")))
        self.ui.pb_Slansre.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\slansre.exe")))
        self.ui.pb_MTFView.clicked.connect(lambda: self.run_slans_exe(Path(r"SLANS_exe\Mtfview\mtfview.exe")))

    def run_slans_exe(self, path):
        """
        Runs SLANS executable files
        Parameters
        ----------
        path

        Returns
        -------

        """
        path = self.main_control.parentDir / fr"exe\{path}"
        t = Thread(target=subprocess.call, args=(path,))
        t.start()

    def serialise(self, state_dict):
        """
        Serialises eigenmode setting UI control widgets
        Parameters
        ----------
        state_dict: dict
            Python dictionary to record the state of each widget associated with eigenmode UI

        Returns
        -------

        """
        serialise(state_dict, self.w_Eigenmode, marker='eigenmode')

    def deserialise(self, state_dict):
        """
        Deserialises eigenmode setting UI control widgets
        Parameters
        ----------
        state_dict: dict
            Python dictionary to record the state of each widget associated with eigenmode UI

        Returns
        -------

        """
        deserialise(state_dict, self.w_Eigenmode, marker='eigenmode')

    @staticmethod
    def run_sequential(n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, pol, parentDir, projectDir,
                       progress_list, sub_dir='', UQ=False, select_solver='slans', mesh_args=None):
        """
        Runs a single instance of SLANS (eigenmode analysis)
        Parameters
        ----------
        n_cells: int, list
            Number of cavity cells
        n_modules: int
            Number of cavity modules
        processor_shape_space: dict
            Dictionary containing geometric dimensions of cavity geometry
        n_modes: int
            Number of eigenmodes to be calculated
        f_shift: float
            Since the eigenmode solver uses the power method, a shift can be provided
        bc: int
            Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
            bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
        pol: int {Monopole, Dipole}
            Defines whether to calculate for monopole or dipole modes
        parentDir: str | path
            Parent directory
        projectDir: str|path
            Project directory
        progress_list: list
            Global list to record the progress of each parallel simulation thread
        sub_dir: dir
            Sub directory in which to write simulation results
        UQ: bool
            Toggles between performing uncertainty quantification or not in addition to the nominal solution
        select_solver: str {slans, ngsolve}
            Select the eigenmode solver to use. Default is the SLANS solver
        mesh_args: list [Jxy, Jxy_bp, Jxy_bp_y]
            Mesh definition for logical mesh:
            Jxy -> Number of elements of logical mesh along JX and JY
            Jxy_bp -> Number of elements of logical mesh along JX in beampipe
            Jxy_bp_y -> Number of elements of logical mesh along JY in beampipe

        Returns
        -------

        """
        progress = 0
        # get length of processor
        # total_no_of_shapes = len(list(processor_shape_space.keys()))

        for key, shape in processor_shape_space.items():
            if select_solver.lower() == 'slans':
                solver = slans_geom
            else:
                solver = ngsolve_mevp

            # run SLANS code
            start_time = time.time()
            expansion = None
            expansion_r = None

            if "EXPANSION" in shape.keys():
                expansion = shape['EXPANSION']

            if 'EXPANSION_R' in shape.keys():
                expansion_r = shape['EXPANSION_R']

            # if len(n_cells) == 1:
            #     n_cells = n_cells[0]
            #     # # create folders for all keys
            #     slans_geom.createFolder(key, projectDir, subdir=sub_dir)
            #
            #     write_cst_paramters(key, shape['IC'], shape['OC'], projectDir=projectDir, cell_type="None")
            #
            # if 'OC_R' in shape.keys(): slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'],
            # shape['OC_R'], n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
            # parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, expansion=expansion,
            # expansion_r=expansion_r) else: slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'],
            # shape['OC'], n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
            # parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, expansion=expansion,
            # expansion_r=expansion_r) else:
            for n_cell in n_cells:
                # # create folders for all keys
                solver.createFolder(f"{key}_n{n_cell}", projectDir, subdir=sub_dir, pol=pol)

                if "CELL TYPE" in shape.keys():
                    if shape['CELL TYPE'] == 'flattop':
                        if 'OC_R' in shape.keys():
                            write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC_R'],
                                                projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                            solver.cavity_flattop(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                          n_modes=n_modes, fid=f"{key}_n{n_cell}", f_shift=f_shift,
                                          bc=bc, pol=pol, beampipes=shape['BP'],
                                          parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                          expansion=expansion, expansion_r=expansion_r, mesh_args=mesh_args)
                        else:
                            write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC'],
                                                projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                            solver.cavity_flattop(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                          n_modes=n_modes, fid=f"{key}_n{n_cell}", f_shift=f_shift,
                                          bc=bc, pol=pol, beampipes=shape['BP'],
                                          parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                          expansion=expansion, expansion_r=expansion_r, mesh_args=mesh_args)
                else:
                    if 'OC_R' in shape.keys():
                        write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC_R'],
                                            projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                        solver.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                      n_modes=n_modes, fid=f"{key}_n{n_cell}", f_shift=f_shift,
                                      bc=bc, pol=pol, beampipes=shape['BP'],
                                      parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                      expansion=expansion, expansion_r=expansion_r, mesh_args=mesh_args)
                    else:
                        write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC'],
                                            projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                        solver.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                      n_modes=n_modes, fid=f"{key}_n{n_cell}", f_shift=f_shift,
                                      bc=bc, pol=pol, beampipes=shape['BP'],
                                      parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                      expansion=expansion, expansion_r=expansion_r, mesh_args=mesh_args)

                # run UQ
                if UQ:
                    uq(f"{key}_n{n_cell}", shape, ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"],
                       n_cells=n_cell, n_modules=n_modules, n_modes=n_modes,
                       f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir,
                       mesh_args=mesh_args, select_solver=select_solver.lower())

                    # uq_ngsolve(f"{key}_n{n_cell}", shape, ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"],
                    #    n_cells=n_cell, n_modules=n_modules, n_modes=n_modes,
                    #    f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir,
                    #    mesh_args=mesh_args, select_solver=select_solver.lower())

                    # uq_ngsolve_parallel(f"{key}_n{n_cell}", shape, ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"],
                    #    n_cells=n_cell, n_modules=n_modules, n_modes=n_modes,
                    #    f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir,
                    #    mesh_args=mesh_args, select_solver=select_solver.lower())

                    # uq_ngsolve_parallel_multicell(f"{key}_n{n_cell}", shape, ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"],
                    #    n_cells=n_cell, n_modules=n_modules, n_modes=n_modes,
                    #    f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir,
                    #    mesh_args=mesh_args, select_solver=select_solver.lower())

                    # uq_multicell(f"{key}_n{n_cell}", shape, ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"],
                    #    n_cells=n_cell, n_modules=n_modules, n_modes=n_modes,
                    #    f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir,
                    #    mesh_args=mesh_args, select_solver=select_solver.lower())

            print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')

            # update progress
            progress_list.append(progress + 1)


def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args,
       select_solver='slans'):
    """

    Parameters
    ----------
    key: str | int
        Cavity geomery identifier
    shape: dict
        Dictionary containing geometric dimensions of cavity geometry
    qois: list
        Quantities of interest considered in uncertainty quantification
    n_cells: int
        Number of cavity cells
    n_modules: int
        Number of modules
    n_modes: int
        Number of eigenmodes to be calculated
    f_shift: float
        Since the eigenmode solver uses the power method, a shift can be provided
    bc: int
        Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
        bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    pol: int {Monopole, Dipole}
        Defines whether to calculate for monopole or dipole modes
    parentDir: str | path
        Parent directory
    projectDir: str|path
        Project directory

    Returns
    -------

    """

    if select_solver.lower() == 'slans':
        uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    else:
        uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'

    err = False
    result_dict_eigen = {}
    eigen_obj_list = qois
    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = shape['IC'][0:5]
    # p_true_el = shape['OC'][0:5]
    print(shape)

    rdim = len(p_true) #+ len(p_true_el)  # How many variabels will be considered as random in our case 5
    degree = 1

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    flag_stroud = 'load_from_file'
    if flag_stroud == 'stroud3':
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
        nodes_ = 2. * nodes_ - 1.
        # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    elif flag_stroud == 'stroud5':
        nodes_, weights_ = cn_leg_05_2(rdim)
    elif flag_stroud == 'cn_gauss':
        nodes_, weights_ = cn_gauss(rdim, 2)
    elif flag_stroud == 'load_from_file':
        nodes_ = pd.read_csv(fr'C:\Users\sosoho\DakotaProjects\Cavity\C3794_cubature5_5\sim_result_table.dat', sep='\s+').iloc[:, 2:7]
        nodes_ = nodes_.to_numpy().T
        weights_ = np.ones((nodes_.shape[1], 1))
    else:
        ic('flag_stroud==1 or flag_stroud==2')
        return 0
    # save nodes
    data_table = pd.DataFrame(nodes_.T)
    data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')

    #  mean value of geometrical parameters
    p_init = np.zeros(np.shape(p_true))
    # p_init_el = np.zeros(np.shape(p_true_el))

    no_parm, no_sims = np.shape(nodes_)
    delta = 0.01  # or 0.1

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    par_end = shape['OC']

    for i in range(no_sims):
        skip = False
        if flag_stroud == 'load_from_file':
            p_init[0] = nodes_[0, i]  # <- A
            p_init[1] = nodes_[1, i]  # <- B
            p_init[2] = nodes_[2, i]  # <- a
            p_init[3] = nodes_[3, i]  # <- b
            p_init[4] = nodes_[4, i]  # <- Ri
        else:
            #
            # p_init[0] = p_true[0] * (1 + delta * nodes_[0, i])  # <- A
            # p_init[1] = p_true[1] * (1 + delta * nodes_[1, i])  # <- B
            # p_init[2] = p_true[2] * (1 + delta * nodes_[2, i])  # <- a
            # p_init[3] = p_true[3] * (1 + delta * nodes_[3, i])  # <- b
            # p_init[4] = p_true[4] * (1 + delta * nodes_[4, i])  # <- Ri

            p_init[0] = p_true[0] + nodes_[0, i]  # <- A
            p_init[1] = p_true[1] + nodes_[1, i]  # <- B
            p_init[2] = p_true[2] + nodes_[2, i]  # <- a
            p_init[3] = p_true[3] + nodes_[3, i]  # <- b
            p_init[4] = p_true[4] + nodes_[4, i]  # <- Ri
            # p_init[5] = p_true[5] + nodes_[5, i]  # <- L
            # p_init[6] = p_true[6] + nodes_[6, i]  # <- Req
            pass

        par_mid = list(np.append(p_init, shape['IC'][5:]))

        # perform checks on geometry
        ok = perform_geometry_checks(par_mid, par_end)
        if not ok:
            err = True
            break
        fid = fr'{key}_Q{i}'

        # skip analysis if folder already exists.
        if not skip:
            if select_solver.lower() == 'slans':
                solver = slans_geom
            else:
                print(' ngsolve selected')
                solver = ngsolve_mevp
            #  run model using SLANS or CST
            # # create folders for all keys
            solver.createFolder(fid, projectDir, subdir=sub_dir)

            if "CELL TYPE" in shape.keys():
                if shape['CELL TYPE'] == 'flattop':
                    # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
                    #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                    try:
                        print(' in flattop')
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end, par_end,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
                    except KeyError:
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end, par_end,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
            else:
                try:
                    solver.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)
                except KeyError:
                    solver.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)

        filename = uq_path / f'{fid}/monopole/qois.json'
        print(filename)
        if os.path.exists(filename):
            # params = fr.svl_reader(filename)
            # norm_length = 2 * n_cells * shape['IC'][5]

            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # print_(qois_result)
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

        data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
        data_table.to_csv(uq_path / 'table.csv', index=False, sep='\t', float_format='%.32f')
        data_table.to_excel(uq_path / 'table.xlsx', index=False)
    # # add original point
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)

    # import matplotlib.pyplot as plt
    print(np.atleast_2d(Ttab_val_f), weights_)
    if not err:
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)

        # append results to dict
        for i, o in enumerate(eigen_obj_list):
            result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
            result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

            # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
            # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

        # plt.show()
        print(result_dict_eigen)
        with open(uq_path / fr"uq.json", 'w') as file:
            file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))
    else:
        print_(fr"There was a problem running UQ analysis for {key}")


def uq_multicell(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args,
       select_solver='slans'):
    """

    Parameters
    ----------
    key: str | int
        Cavity geomery identifier
    shape: dict
        Dictionary containing geometric dimensions of cavity geometry
    qois: list
        Quantities of interest considered in uncertainty quantification
    n_cells: int
        Number of cavity cells
    n_modules: int
        Number of modules
    n_modes: int
        Number of eigenmodes to be calculated
    f_shift: float
        Since the eigenmode solver uses the power method, a shift can be provided
    bc: int
        Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
        bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    pol: int {Monopole, Dipole}
        Defines whether to calculate for monopole or dipole modes
    parentDir: str | path
        Parent directory
    projectDir: str|path
        Project directory

    Returns
    -------

    """

    if select_solver.lower() == 'slans':
        uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    else:
        uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'

    err = False
    result_dict_eigen = {}
    eigen_obj_list = qois
    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    # expected input
    # l_end_cell = np.array([73.52, 131.75, 106.25, 118.7, 150, 187, 350])
    #
    # mid_cell = np.array([[[60, 60], [56, 54], [76, 56]],  # <- A
    #                      [[60, 43], [73, 65], [65, 45]],  # <- B
    #                      [[34, 50], [54, 50], [37, 40]],  # <- a
    #                      [[40, 36], [56, 57], [54, 56]],  # <- b
    #                      [[150, 151], [170, 165], [150, 155]],  # <- Ri
    #                      [[187, 187], [184, 176], [178, 170]],  # <- L
    #                      [[369.6321578127116, 345], [340, 350], [350, 360]]])
    #
    # r_end_cell = np.array([70, 120, 100, 110, 130, 176, 340])
    cav_var_list = ['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']
    midcell_var_dict = dict()
    for i1 in range(len(cav_var_list)):
        for i2 in range(n_cells):
            for i3 in range(2):
                midcell_var_dict[f'{cav_var_list[i1]}_{i2}_m{i3}'] = [i1, i2, i3]

    # create random variables
    multicell_mid_vars = create_multicell_random_variables(n_cells, np.atleast_2d(np.array(shape['IC'])[:7]).T)

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = [np.array(shape['OC'])[:7], multicell_mid_vars, np.array(shape['OC_R'])[:7]]
    print(shape)

    rdim = len(np.array(shape['OC'])[:7]) + multicell_mid_vars.size + len(np.array(shape['OC_R'])[:7])  # How many variabels will be considered as random in our case 5
    # ic(rdim, multicell_mid_vars.size)
    degree = 1

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    flag_stroud = 'stroud3'
    if flag_stroud == 'stroud3':
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
        nodes_ = 2. * nodes_ - 1.
        # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    elif flag_stroud == 'stroud5':
        nodes_, weights_ = cn_leg_05_2(rdim)
    elif flag_stroud == 'cn_gauss':
        nodes_, weights_ = cn_gauss(rdim, 2)
    else:
        ic('flag_stroud==1 or flag_stroud==2')
        return 0

    no_parm, no_sims = np.shape(nodes_)
    delta = 0.01  # or 0.1

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    # par_end = shape['OC']
    ic(nodes_)
    # save nodes
    data_table = pd.DataFrame(nodes_.T, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')

    for i in range(no_sims):
        skip = False
        # ic(nodes_[0:len(p_true[0]), i])
        # ic(nodes_[len(p_true[0]):len(p_true[0])+p_true[1].size, i].reshape(np.shape(p_true[1])))
        # ic(nodes_[len(p_true[0])+p_true[1].size:, i])
        p_init_el = p_true[0] + nodes_[0:len(p_true[0]), i]

        p_init_m = p_true[1] + nodes_[len(p_true[0]):len(p_true[0])+p_true[1].size, i].reshape(np.shape(p_true[1]))

        p_init_er = p_true[2] + nodes_[len(p_true[0])+p_true[1].size:, i]

        par_mid = p_init_m
        # ic(par_mid)
        par_end_l = p_init_el
        # ic(par_end_l)
        par_end_r = p_init_er
        # ic(par_end_r)

        # # perform checks on geometry
        # ok = perform_geometry_checks(par_mid, par_end)
        # if not ok:
        #     err = True
        #     break
        fid = fr'{key}_Q{i}'

        # skip analysis if folder already exists.
        if not skip:
            if select_solver.lower() == 'slans':
                solver = slans_geom
            else:
                print(' ngsolve selected')
                solver = ngsolve_mevp
            #  run model using SLANS or CST
            # # create folders for all keys
            solver.createFolder(fid, projectDir, subdir=sub_dir)

            if "CELL TYPE" in shape.keys():
                if shape['CELL TYPE'] == 'flattop':
                    # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
                    #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                    try:
                        print(' in flattop')
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
                    except KeyError:
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
            else:
                try:
                    solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)
                except KeyError:
                    solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)

        filename = uq_path / f'{fid}/monopole/qois.json'
        print(filename)
        if os.path.exists(filename):
            # params = fr.svl_reader(filename)
            # norm_length = 2 * n_cells * shape['IC'][5]

            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # print_(qois_result)
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

        data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
        data_table.to_csv(uq_path / 'table.csv', index=False, sep='\t', float_format='%.32f')

    # # add original point
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)

    # import matplotlib.pyplot as plt
    print(np.atleast_2d(Ttab_val_f), weights_)
    if not err:
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)

        # append results to dict
        for i, o in enumerate(eigen_obj_list):
            result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
            result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

            # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
            # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

        # plt.show()
        print(result_dict_eigen)
        with open(uq_path / fr"uq.json", 'w') as file:
            file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))
    else:
        print_(fr"There was a problem running UQ analysis for {key}")


def uq_ngsolve(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args,
               select_solver='slans'):
    """

    Parameters
    ----------
    key: str | int
        Cavity geomery identifier
    shape: dict
        Dictionary containing geometric dimensions of cavity geometry
    qois: list
        Quantities of interest considered in uncertainty quantification
    n_cells: int
        Number of cavity cells
    n_modules: int
        Number of modules
    n_modes: int
        Number of eigenmodes to be calculated
    f_shift: float
        Since the eigenmode solver uses the power method, a shift can be provided
    bc: int
        Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
        bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    pol: int {Monopole, Dipole}
        Defines whether to calculate for monopole or dipole modes
    parentDir: str | path
        Parent directory
    projectDir: str|path
        Project directory

    Returns
    -------

    """

    if select_solver.lower() == 'slans':
        uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    else:
        uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'

    err = False
    result_dict_eigen = {}
    eigen_obj_list = qois
    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    rdim = n_cells * 3  # How many variables will be considered as random in our case 5
    degree = 1

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    flag_stroud = 'stroud5'

    if flag_stroud == 'stroud3':
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
        nodes_ = 2. * nodes_ - 1.
        # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    elif flag_stroud == 'stroud5':
        nodes_, weights_ = cn_leg_05_2(rdim)
    elif flag_stroud == 'cn_gauss':
        nodes_, weights_ = cn_gauss(rdim, 2)
    else:
        ic('flag_stroud==1 or flag_stroud==2')
        return 0

    ic(nodes_)
    # save nodes
    data_table = pd.DataFrame(nodes_.T, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')

    #  mean value of geometrical parameters
    no_parm, no_sims = np.shape(nodes_)

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run

    for i in range(no_sims):
        skip = False
        # perform checks on geometry
        ok = perform_geometry_checks(shape['IC'], shape['OC'])
        if not ok:
            err = True
            break
        fid = fr'{key}_Q{i}'

        # skip analysis if folder already exists.
        if not skip:
            solver = ngsolve_mevp
            #  run model using SLANS or CST
            # # create folders for all keys
            solver.createFolder(fid, projectDir, subdir=sub_dir)

            if "CELL TYPE" in shape.keys():
                if shape['CELL TYPE'] == 'flattop':
                    # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
                    #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                    try:
                        print(' in flattop')
                        solver.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args,
                                              deformation_params=nodes_[:, i])
                    except KeyError:
                        solver.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args,
                                              deformation_params=nodes_[:, i])
            else:
                try:
                    solver.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args,
                                  deformation_params=nodes_[:, i])
                except KeyError:
                    solver.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args,
                                  deformation_params=nodes_[:, i])

        filename = uq_path / f'{fid}/monopole/qois.json'
        print(filename)
        if os.path.exists(filename):
            # params = fr.svl_reader(filename)
            # norm_length = 2 * n_cells * shape['IC'][5]

            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # print_(qois_result)
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
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)
    # save table
    data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / 'table.csv', index=False, sep='\t', float_format='%.32f')

    print(np.atleast_2d(Ttab_val_f), weights_)
    if not err:
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)

        # append results to dict
        for i, o in enumerate(eigen_obj_list):
            result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
            result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

            # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
            # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

        # plt.show()
        print(result_dict_eigen)
        with open(uq_path / fr"uq.json", 'w') as file:
            file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))
    else:
        print_(fr"There was a problem running UQ analysis for {key}")

def uq_ngsolve_parallel(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args, select_solver='slans'):
    """

    Parameters
    ----------
    key: str | int
        Cavity geomery identifier
    shape: dict
        Dictionary containing geometric dimensions of cavity geometry
    qois: list
        Quantities of interest considered in uncertainty quantification
    n_cells: int
        Number of cavity cells
    n_modules: int
        Number of modules
    n_modes: int
        Number of eigenmodes to be calculated
    f_shift: float
        Since the eigenmode solver uses the power method, a shift can be provided
    bc: int
        Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
        bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    pol: int {Monopole, Dipole}
        Defines whether to calculate for monopole or dipole modes
    parentDir: str | path
        Parent directory
    projectDir: str|path
        Project directory

    Returns
    -------

    """

    print("Starting parrallel")
    if select_solver.lower() == 'slans':
        uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    else:
        uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'

    err = False
    result_dict_eigen = {}
    eigen_obj_list = qois
    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    rdim = 18  # How many variables will be considered as random in our case 5
    degree = 2

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    flag_stroud = 'stroud3'

    if flag_stroud == 'stroud3':
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
        nodes_ = 2. * nodes_ - 1.
        # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    elif flag_stroud == 'stroud5':
        nodes_, weights_ = cn_leg_05_2(rdim)
    elif flag_stroud == 'cn_gauss':
        nodes_, weights_ = cn_gauss(rdim, 2)
    elif flag_stroud == 'lhc':
        sampler = qmc.LatinHypercube(d=rdim)
        _ = sampler.reset()
        nsamp = 2500
        sample = sampler.random(n=nsamp)
        # ic(qmc.discrepancy(sample))
        l_bounds = [-1, -1, -1, -1, -1, -1]
        u_bounds = [1, 1, 1, 1, 1, 1]
        sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

        nodes_, weights_ = sample_scaled.T, np.ones((nsamp, 1))
    else:
        ic('flag_stroud==1 or flag_stroud==2')
        return 0

    ic(nodes_)

    #  mean value of geometrical parameters
    no_parm, no_sims = np.shape(nodes_)

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    processes = []
    manager = mp.Manager()

    progress_list = manager.list()
    progress_list.append(0)
    proc_count = 25
    if proc_count > no_sims:
        proc_count = no_sims

    share = round(no_sims / proc_count)

    for p in range(proc_count):
        # try:
        end_already = False
        if p != proc_count - 1:
            if (p+1)*share < no_sims:
                proc_keys_list = np.arange(p * share, p * share + share)
            else:
                proc_keys_list = np.arange(p * share, no_sims)
                end_already = True

        if p == proc_count - 1 and not end_already:
            proc_keys_list = np.arange(p * share, no_sims)

        # ic(proc_keys_list)
        processor_nodes = nodes_[:, proc_keys_list]
        processor_weights = weights_[proc_keys_list]
        # ic(processor_nodes)
        # ic(processor_weights)

        skip = False
        # perform checks on geometry
        ok = perform_geometry_checks(shape['IC'], shape['OC'])
        if not ok:
            err = True
            break

        service = mp.Process(target=uq_piotr_sequential, args=(
            n_cells, n_modules, shape, qois, n_modes, f_shift, bc, pol, parentDir,
            projectDir, sub_dir, select_solver, mesh_args, key, uq_path,
            proc_keys_list, processor_nodes, processor_weights, p))

        service.start()

        processes.append(psutil.Process(service.pid))


def uq_piotr_sequential(n_cells, n_modules, shape, qois, n_modes, f_shift, bc, pol, parentDir, projectDir, sub_dir,
                  select_solver, mesh_args, key, uq_path, proc_keys_list, processor_nodes, processor_weights, proc_num):

    err = False
    result_dict_eigen = {}
    Ttab_val_f = []
    eigen_obj_list = qois
    delta = 0.01

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = shape['IC'][0:6]
    p_true_el = shape['OC'][0:6]
    p_true_er = shape['OC'][0:6]

    #  mean value of geometrical parameters
    p_init = np.zeros(np.shape(p_true))
    p_init_el = np.zeros(np.shape(p_true_el))
    p_init_er = np.zeros(np.shape(p_true_er))

    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    for i1 in proc_keys_list:
        skip = False
        p_init[0] = p_true[0] * (1 + delta * processor_nodes[0, i1-min(proc_keys_list)])
        p_init[1] = p_true[1] * (1 + delta * processor_nodes[1, i1-min(proc_keys_list)])
        p_init[2] = p_true[2] * (1 + delta * processor_nodes[2, i1-min(proc_keys_list)])
        p_init[3] = p_true[3] * (1 + delta * processor_nodes[3, i1-min(proc_keys_list)])
        p_init[4] = p_true[4] * (1 + delta * processor_nodes[4, i1-min(proc_keys_list)])

        p_init_el[0] = p_true_el[0] * (1 + delta * processor_nodes[5, i1-min(proc_keys_list)])
        p_init_el[1] = p_true_el[1] * (1 + delta * processor_nodes[6, i1-min(proc_keys_list)])
        p_init_el[2] = p_true_el[2] * (1 + delta * processor_nodes[7, i1-min(proc_keys_list)])
        p_init_el[3] = p_true_el[3] * (1 + delta * processor_nodes[8, i1-min(proc_keys_list)])
        p_init_el[4] = p_true_el[4] * (1 + delta * processor_nodes[9, i1-min(proc_keys_list)])

        # p_init[0] = p_true[0] + processor_nodes[0, i1-min(proc_keys_list)]  # <- A
        # p_init[1] = p_true[1] + processor_nodes[1, i1-min(proc_keys_list)]  # <- B
        # p_init[2] = p_true[2] + processor_nodes[2, i1-min(proc_keys_list)]  # <- a
        # p_init[3] = p_true[3] + processor_nodes[3, i1-min(proc_keys_list)]  # <- b
        # p_init[4] = p_true[4] + processor_nodes[4, i1-min(proc_keys_list)]  # <- Ri
        # p_init[5] = p_true[5] + processor_nodes[5, i1-min(proc_keys_list)]  # <- L
        # p_init[6] = p_true[6] + processor_nodes[6, i1-min(proc_keys_list)]  # <- Req
        #
        # p_init_el[0] = p_true_el[0] + processor_nodes[6, i1-min(proc_keys_list)]  # <- A
        # p_init_el[1] = p_true_el[1] + processor_nodes[7, i1-min(proc_keys_list)]  # <- B
        # p_init_el[2] = p_true_el[2] + processor_nodes[8, i1-min(proc_keys_list)]  # <- a
        # p_init_el[3] = p_true_el[3] + processor_nodes[9, i1-min(proc_keys_list)]  # <- b
        # p_init_el[4] = p_true_el[4] + processor_nodes[10, i1-min(proc_keys_list)]  # <- Ri
        # p_init_el[5] = p_true_el[5] + processor_nodes[11, i1-min(proc_keys_list)]  # <- L
        # p_init_el[6] = p_true_el[6] + processor_nodes[11, i1-min(proc_keys_list)]  # <- Req
        #
        # p_init_er[0] = p_true_er[0] + processor_nodes[12, i1-min(proc_keys_list)]  # <- A
        # p_init_er[1] = p_true_er[1] + processor_nodes[13, i1-min(proc_keys_list)]  # <- B
        # p_init_er[2] = p_true_er[2] + processor_nodes[14, i1-min(proc_keys_list)]  # <- a
        # p_init_er[3] = p_true_er[3] + processor_nodes[15, i1-min(proc_keys_list)]  # <- b
        # p_init_er[4] = p_true_er[4] + processor_nodes[16, i1-min(proc_keys_list)]  # <- Ri
        # p_init_er[5] = p_true_er[5] + processor_nodes[17, i1-min(proc_keys_list)]  # <- L
        # p_init_er[6] = p_true_er[6] + processor_nodes[18, i1-min(proc_keys_list)]  # <- Req

        par_mid = list(np.append(p_init, shape['IC'][6:]))
        par_end_l = list(np.append(p_init_el, shape['OC'][6:]))
        # par_end_r = list(np.append(p_init_er, shape['OC_R'][6:]))
        par_end_r = shape['OC_R']

        # # perform checks on geometry
        # ok = perform_geometry_checks(par_mid, par_end_l)
        # if not ok:
        #     err = True
        #     break
        fid = fr'{key}_Q{i1}'

        # skip analysis if folder already exists.
        if not skip:
            if select_solver.lower() == 'slans':
                solver = slans_geom
            else:
                print(' ngsolve selected')
                solver = ngsolve_mevp
            #  run model using SLANS or CST
            # # create folders for all keys
            solver.createFolder(fid, projectDir, subdir=sub_dir)

            if "CELL TYPE" in shape.keys():
                if shape['CELL TYPE'] == 'flattop':
                    # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
                    #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                    try:
                        print(' in flattop')
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
                    except KeyError:
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
            else:
                try:
                    solver.cavity(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)
                except KeyError:
                    solver.cavity(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)

        filename = uq_path / f'{fid}/monopole/qois.json'
        # print(filename)
        if os.path.exists(filename):
            # params = fr.svl_reader(filename)
            # norm_length = 2 * n_cells * shape['IC'][5]

            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # print_(qois_result)
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
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)
    # save table
    data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / fr'table_{proc_num}.csv', index=False, sep='\t', float_format='%.32f')


def uq_sequential(n_cells, n_modules, shape, qois, n_modes, f_shift, bc, pol, parentDir, projectDir, sub_dir,
                  select_solver, mesh_args, key, uq_path, proc_keys_list, processor_nodes, processor_weights, proc_num):

    err = False
    result_dict_eigen = {}
    Ttab_val_f = []
    eigen_obj_list = qois

    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    for i1 in proc_keys_list:
        skip = False
        # perform checks on geometry
        ok = perform_geometry_checks(shape['IC'], shape['OC'])
        if not ok:
            err = True
            break
        fid = fr'{key}_Q{i1}'

        # skip analysis if folder already exists.
        if not skip:
            solver = ngsolve_mevp
            #  run model using SLANS or CST
            # # create folders for all keys
            solver.createFolder(fid, projectDir, subdir=sub_dir)

            if "CELL TYPE" in shape.keys():
                if shape['CELL TYPE'] == 'flattop':
                    # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
                    #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                    try:
                        print(' in flattop')
                        solver.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args,
                                              deformation_params=processor_nodes[:, i1-min(proc_keys_list)])
                    except KeyError:
                        solver.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args,
                                              deformation_params=processor_nodes[:, i1-min(proc_keys_list)])
            else:
                try:
                    solver.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args,
                                  deformation_params=processor_nodes[:, i1-min(proc_keys_list)])
                except KeyError:
                    solver.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args,
                                  deformation_params=processor_nodes[:, i1-min(proc_keys_list)])

        filename = uq_path / f'{fid}/monopole/qois.json'
        # print(filename)
        if os.path.exists(filename):
            # params = fr.svl_reader(filename)
            # norm_length = 2 * n_cells * shape['IC'][5]

            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # print_(qois_result)
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
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)
    # save table
    data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / fr'table_{proc_num}.csv', index=False, sep='\t', float_format='%.32f')

    # print(np.atleast_2d(Ttab_val_f), processor_weights)


def uq_ngsolve_parallel_multicell(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args, select_solver='slans'):
    """

    Parameters
    ----------
    key: str | int
        Cavity geomery identifier
    shape: dict
        Dictionary containing geometric dimensions of cavity geometry
    qois: list
        Quantities of interest considered in uncertainty quantification
    n_cells: int
        Number of cavity cells
    n_modules: int
        Number of modules
    n_modes: int
        Number of eigenmodes to be calculated
    f_shift: float
        Since the eigenmode solver uses the power method, a shift can be provided
    bc: int
        Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
        bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    pol: int {Monopole, Dipole}
        Defines whether to calculate for monopole or dipole modes
    parentDir: str | path
        Parent directory
    projectDir: str|path
        Project directory

    Returns
    -------

    """

    print("Starting parallel")
    if select_solver.lower() == 'slans':
        uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    else:
        uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'

    err = False
    result_dict_eigen = {}
    eigen_obj_list = qois
    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    cav_var_list = ['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']
    midcell_var_dict = dict()
    for i1 in range(len(cav_var_list)):
        for i2 in range(n_cells):
            for i3 in range(2):
                midcell_var_dict[f'{cav_var_list[i1]}_{i2}_m{i3}'] = [i1, i2, i3]

    # create random variables
    multicell_mid_vars = create_multicell_random_variables(n_cells, np.atleast_2d(np.array(shape['IC'])[:7]).T)

    if n_cells == 1:
        # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
        p_true = [np.array(shape['OC'])[:7], np.array(shape['OC_R'])[:7]]
        rdim = len(np.array(shape['OC'])[:7]) + len(np.array(shape['OC_R'])[:7])  # How many variabels will be considered as random in our case 5
    else:
        # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
        p_true = [np.array(shape['OC'])[:7], multicell_mid_vars, np.array(shape['OC_R'])[:7]]
        rdim = len(np.array(shape['OC'])[:7]) + multicell_mid_vars.size + len(np.array(shape['OC_R'])[:7])  # How many variabels will be considered as random in our case 5
        # rdim = rdim - (n_cells*2 - 1)  # <- reduce dimension by making iris and equator radii to be equal
    # ic(rdim, multicell_mid_vars.size)

    # rdim = n_cells*3  # How many variables will be considered as random in our case 5
    degree = 1

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    flag_stroud = 'load_from_file'

    if flag_stroud == 'stroud3':
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
        nodes_ = 2. * nodes_ - 1.
        # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    elif flag_stroud == 'stroud5':
        nodes_, weights_ = cn_leg_05_2(rdim)
    elif flag_stroud == 'cn_gauss':
        nodes_, weights_ = cn_gauss(rdim, 2)
    elif flag_stroud == 'lhc':
        sampler = qmc.LatinHypercube(d=rdim)
        _ = sampler.reset()
        nsamp = 3000
        sample = sampler.random(n=nsamp)
        # ic(qmc.discrepancy(sample))
        l_bounds = -np.ones(rdim)
        u_bounds = np.ones(rdim)
        sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

        nodes_, weights_ = sample_scaled.T, np.ones((nsamp, 1))
    elif flag_stroud == 'load_from_file':
        nodes_ = pd.read_csv(fr'C:\Users\sosoho\DakotaProjects\Cavity\C3795_lhs\sim_result_table.dat', sep='\s+').iloc[:, 2:-2]
        nodes_ = nodes_.to_numpy().T
        weights_ = np.ones((nodes_.shape[1], 1))
    else:
        ic('flag_stroud==1 or flag_stroud==2')
        return 0

    ic(nodes_)

    # save nodes
    data_table = pd.DataFrame(nodes_.T)
    data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')

    #  mean value of geometrical parameters
    # ic(np.shape(nodes_))
    no_parm, no_sims = np.shape(nodes_)

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    processes = []
    manager = mp.Manager()

    progress_list = manager.list()
    progress_list.append(0)
    proc_count = 30
    if proc_count > no_sims:
        proc_count = no_sims

    share = round(no_sims / proc_count)

    for p in range(proc_count):
        # try:
        end_already = False
        if p != proc_count - 1:
            if (p+1)*share < no_sims:
                proc_keys_list = np.arange(p * share, p * share + share)
            else:
                proc_keys_list = np.arange(p * share, no_sims)
                end_already = True

        if p == proc_count - 1 and not end_already:
            proc_keys_list = np.arange(p * share, no_sims)

        # ic(proc_keys_list)
        processor_nodes = nodes_[:, proc_keys_list]
        processor_weights = weights_[proc_keys_list]
        # ic(processor_nodes)
        # ic(processor_weights)

        skip = False
        # perform checks on geometry
        ok = perform_geometry_checks(shape['IC'], shape['OC'])
        if not ok:
            err = True
            break

        service = mp.Process(target=uq_multicell_sequential, args=(
            n_cells, n_modules, shape, qois, n_modes, f_shift, bc, pol, parentDir,
            projectDir, sub_dir, select_solver, mesh_args, key, uq_path,
            proc_keys_list, processor_nodes, processor_weights, p, p_true))

        service.start()

        processes.append(psutil.Process(service.pid))


def uq_multicell_sequential(n_cells, n_modules, shape, qois, n_modes, f_shift, bc, pol, parentDir, projectDir, sub_dir,
                  select_solver, mesh_args, key, uq_path, proc_keys_list, processor_nodes, processor_weights, proc_num, p_true):
    start = time.time()
    err = False
    result_dict_eigen = {}
    Ttab_val_f = []
    eigen_obj_list = qois

    for o in qois:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    for i1 in proc_keys_list:
        skip = False
        if n_cells == 1:
            p_init_el = p_true[0] + processor_nodes[0:len(p_true[0]), i1-min(proc_keys_list)]
            p_init_er = p_true[1] + processor_nodes[len(p_true[0]):, i1-min(proc_keys_list)]
            par_mid = p_init_el
        else:
            proc_node = processor_nodes[:, i1-min(proc_keys_list)]
            # # one dimension of nodes_ is dimension of number of variables. The variables must be expanded to the unreduced dimension
            # # by filling the missing slots with radius values. Insert from end of list, index(Req) + 7 and index(Ri) + 6
            # moved_val_indx = 7
            # proc_nodes_len = len(processor_nodes)
            # for i2 in range(2 * n_cells - 1):
            #     if i2 % 2 == 0:
            #         # print("num", proc_nodes_len - moved_val_indx, "->", proc_nodes_len - moved_val_indx + 7)
            #         proc_node = np.insert(proc_node, proc_nodes_len - moved_val_indx + 7, proc_node[proc_nodes_len - moved_val_indx])
            #         # update index
            #         moved_val_indx += 7
            #     else:
            #         # print("num", proc_nodes_len - moved_val_indx, "->", proc_nodes_len - moved_val_indx + 6)
            #         proc_node = np.insert(proc_node, proc_nodes_len - moved_val_indx + 6, proc_node[proc_nodes_len - moved_val_indx])
            #         moved_val_indx += 5

            p_init_el = processor_nodes[0:len(p_true[0]), i1-min(proc_keys_list)]
            p_init_m = processor_nodes[len(p_true[0]):len(p_true[0])+p_true[1].size, i1-min(proc_keys_list)].reshape(np.shape(p_true[1])[::-1]).T
            p_init_er = processor_nodes[len(p_true[0])+p_true[1].size:, i1-min(proc_keys_list)]

            # p_init_el = p_true[0] + processor_nodes[0:len(p_true[0]), i1-min(proc_keys_list)]
            # p_init_m = p_true[1] + processor_nodes[len(p_true[0]):len(p_true[0])+p_true[1].size, i1-min(proc_keys_list)].reshape(np.shape(p_true[1]))
            # p_init_er = p_true[2] + processor_nodes[len(p_true[0])+p_true[1].size:, i1-min(proc_keys_list)]
            # ic(proc_node, proc_node.shape)

            # p_init_el = p_true[0] + proc_node[0:len(p_true[0])]
            # p_init_m = p_true[1] + proc_node[len(p_true[0]):len(p_true[0])+p_true[1].size].reshape(np.shape(p_true[1]))
            # p_init_er = p_true[2] + proc_node[len(p_true[0])+p_true[1].size:]

            par_mid = p_init_m

        par_end_l = p_init_el
        par_end_r = p_init_er

        fid = fr'{key}_Q{i1}'

        # check if folder already exist (simulation already completed)
        if os.path.exists(uq_path / f'{fid}/monopole/qois.json'):
            skip = True
            print(f'processor {proc_num} skipped ', fid)

        # skip analysis if folder already exists.
        if not skip:
            solver = ngsolve_mevp
            #  run model using SLANS or CST
            # # create folders for all keys
            solver.createFolder(fid, projectDir, subdir=sub_dir)

            if "CELL TYPE" in shape.keys():
                if shape['CELL TYPE'] == 'flattop':
                    # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
                    #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                    try:
                        print(' in flattop')
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
                    except KeyError:
                        solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                              n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                              mesh_args=mesh_args)
            else:
                try:
                    solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)
                except KeyError:
                    solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)

        filename = uq_path / f'{fid}/monopole/qois.json'
        # print(filename)
        if os.path.exists(filename):
            # params = fr.svl_reader(filename)
            # norm_length = 2 * n_cells * shape['IC'][5]

            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # print_(qois_result)
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
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)
    # save table
    # print(f'\tDone with proc {proc_num} ', time.time()-start)
    data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / fr'table_{proc_num}.csv', index=False, sep='\t', float_format='%.32f')

    # print(np.atleast_2d(Ttab_val_f), processor_weights)


def gather_uq(uq_path, no_of_processes):
    Ttab_val_f_list = []
    weights = []
    for i1 in range(no_of_processes):
        if os.path.exists(uq_path / fr'table_{i1}.csv'):
            Ttab_val_f_list.append(pd.read_csv(uq_path / fr'table_{i1}.csv', sep='\t').to_numpy())
            weights = np.vstack(pd.read_csv(uq_path / fr'weight_{i1}.csv', sep='\t').to_numpy())
        else:
            print(fr'Inspect result:: table_{i1}.csv')

    Ttab_val_f = pd.concat(Ttab_val_f_list, ignore_index=True)
    print(np.atleast_2d(Ttab_val_f), weights_)

    v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)

    # append results to dict
    for i, o in enumerate(eigen_obj_list):
        result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
        result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

        # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
        # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

    # plt.show()
    print(result_dict_eigen)
    with open(uq_path / fr"uq.json", 'w') as file:
        file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))


def get_qoi_value(d, obj):
    """
    Gets the quantities of interest from simulation results
    Parameters
    ----------
    d: dict
        Dictionary containing several figures of merits from eigenmode solver
    obj: list
        List of objective functions
    n_cells: int
        Number of cells
    norm_length: float
        Normalisation length for :math: `E_\mathrm{acc}`

    Returns
    -------

    """
    # Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
    # Freq = d['FREQUENCY'][n_cells - 1]
    # E_stored = d['STORED ENERGY'][n_cells - 1]
    # # Rsh = d['SHUNT IMPEDANCE'][n_cells-1]  # MOhm
    # Q = d['QUALITY FACTOR'][n_cells - 1]
    # Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
    # Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
    # # Vacc = dict['ACCELERATION'][0]
    # # Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells-1]  # MV/m
    # Rsh_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm
    #
    # Vacc = np.sqrt(
    #     2 * Rsh_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6
    # # factor of 2, remember circuit and accelerator definition
    # # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
    # Eacc = Vacc / (norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
    # Epk_Eacc = Epk / Eacc
    # Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc
    #
    # d = {
    #     "Req": Req,
    #     "freq": Freq,
    #     "Q": Q,
    #     "E": E_stored,
    #     "R/Q": 2 * Rsh_Q,
    #     "Epk/Eacc": Epk_Eacc,
    #     "Bpk/Eacc": Bpk_Eacc
    # }

    objective = []

    # append objective functions
    for o in obj:
        if o in d.keys():
            objective.append(d[o])

    return objective


def create_multicell_random_variables(n_cell, mid_cell):
    print(mid_cell)
    mid_cell_pair = np.hstack([mid_cell, mid_cell])
    # print(mid_cell_pair)
    if n_cell > 2:
        mid_muilticell = np.hstack([mid_cell_pair for _ in range(n_cell - 1)]).reshape(len(mid_cell), n_cell-1, 2)
        return mid_muilticell
    else:
        return np.expand_dims(mid_cell_pair, axis=1)

