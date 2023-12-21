import os.path
import subprocess
import time
import multiprocessing as mp
from threading import Thread
from pathlib import Path
from psutil import NoSuchProcess
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
                uq(key, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
                   n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
                   f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir)

            print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')

            # update progress
            progress_list.append(progress + 1)

            # else:
            #     # run own eigenmode code
            #     folder = projectDir / fr'SimulationData\NativeEig'
            #     mod = Model(folder=folder, name=f"{key}", parent_dir=parentDir)
            #
            #     try:
            #         # convert first to m.
            #         mid_cell = np.array(shape['IC'])*1e-3
            #         end_cell_left = np.array(shape['OC'])*1e-3
            #         end_cell_right = np.array(shape['OC_R'])*1e-3
            #
            #         mod.run(n_cells, mid_cell, end_cell_left, end_cell_right, beampipe=shape['BP'],
            #                 req_mode_num=int(n_modes), plot=False)
            #     except KeyError:
            #         # convert first to m.
            #         mid_cell = np.array(shape['IC'])*1e-3
            #         end_cell_left = np.array(shape['OC'])*1e-3
            #         end_cell_right = np.array(shape['OC'])*1e-3
            #
            #         mod.run(n_cells, mid_cell, end_cell_left, end_cell_right, beampipe=shape['BP'],
            #                 req_mode_num=int(n_modes), plot=False)


def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir):
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
    err = False
    result_dict_slans = {}
    slans_obj_list = qois
    for o in qois:
        result_dict_slans[o] = {'expe': [], 'stdDev': []}

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = shape['IC'][0:5]

    rdim = len(p_true)  # How many variabels will be considered as random in our case 5
    degree = 1

    #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    flag_stroud = 1
    if flag_stroud == 1:
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
        nodes_ = 2. * nodes_ - 1.
    elif flag_stroud == 2:
        nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)  # change to stroud 5 later
        nodes_ = 2. * nodes_ - 1.
    else:
        ic('flag_stroud==1 or flag_stroud==2')
        return 0

    #  mean value of geometrical parameters
    p_init = np.zeros(np.shape(p_true))

    no_parm, no_sims = np.shape(nodes_)
    delta = 0.005  # or 0.1

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    for i in range(no_sims):
        skip = False
        p_init[0] = p_true[0] * (1 + delta * nodes_[0, i])
        p_init[1] = p_true[1] * (1 + delta * nodes_[1, i])
        p_init[2] = p_true[2] * (1 + delta * nodes_[2, i])
        p_init[3] = p_true[3] * (1 + delta * nodes_[3, i])
        p_init[4] = p_true[4] * (1 + delta * nodes_[4, i])

        par_mid = np.append(p_init, shape['IC'][5:]).tolist()
        par_end = par_mid

        # perform checks on geometry
        ok = perform_geometry_checks(par_mid, par_end)
        if not ok:
            err = True
            break
        fid = fr'{key}_Q{i}'

        # check if folder exists and skip if it does
        if os.path.exists(projectDir / fr'SimulationData\SLANS\{key}\{fid}'):
            skip = True

        # skip analysis if folder already exists.
        if not skip:
            #  run model using SLANS or CST
            # # create folders for all keys
            slans_geom.createFolder(fid, projectDir, subdir=sub_dir)
            try:
                slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
            except KeyError:
                slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

        filename = projectDir / fr'SimulationData\SLANS\{key}\{fid}\cavity_{bc}.svl'
        if os.path.exists(filename):
            params = fr.svl_reader(filename)
            norm_length = 2 * n_cells * shape['IC'][5]

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
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)

    # import matplotlib.pyplot as plt
    if not err:
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)

        # append results to dict
        for i, o in enumerate(slans_obj_list):
            result_dict_slans[o]['expe'].append(v_expe_fobj[i])
            result_dict_slans[o]['stdDev'].append(v_stdDev_fobj[i])

            # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
            # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

        # plt.show()

        with open(projectDir / fr"SimulationData\SLANS\{key}\uq.json", 'w') as file:
            file.write(json.dumps(result_dict_slans, indent=4, separators=(',', ': ')))
    else:
        print_(fr"There was a problem running UQ analysis for {key}")


def get_qoi_value(d, obj, n_cells, norm_length):
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
    Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
    Freq = d['FREQUENCY'][n_cells - 1]
    E_stored = d['STORED ENERGY'][n_cells - 1]
    # Rsh = d['SHUNT IMPEDANCE'][n_cells-1]  # MOhm
    Q = d['QUALITY FACTOR'][n_cells - 1]
    Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
    Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
    # Vacc = dict['ACCELERATION'][0]
    # Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells-1]  # MV/m
    Rsh_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm

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
