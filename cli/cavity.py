import os
import shutil
import time
from distutils import dir_util
from pathlib import Path

import pandas as pd
from termcolor import colored
from tqdm import tqdm
import time
from analysis_modules.tune.tuners.tuner import Tuner
from analysis_modules.data_module.abci_data import ABCIData
import numpy as np
import multiprocessing as mp
from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from manual_run.slans.slansEigen import SLANSEigen
from analysis_modules.eigenmode.customEig.run_field_solver import Model
from analysis_modules.wakefield.ABCI.abci_geometry import ABCIGeometry
from utils.shared_functions import *
import psutil

slans_geom = SLANSEigen()
abci_geom = ABCIGeometry()
custom_eig = Model()
tuner = Tuner()

SOFTWARE_DIRECTORY = os.getcwd()  # str(Path().parents[0])

m0 = 9.1093879e-31
q0 = 1.6021773e-19
c0 = 2.99792458e8
mu0 = 4 * np.pi * 1e-7
eps0 = 8.85418782e-12


def error(*arg):
    print(colored(f'\t{arg}', 'red'))


def running(*arg):
    print(colored(f'\t{arg}', 'yellow'))


def info(*arg):
    print(colored(f'\t{arg}', 'blue'))


def done(*arg):
    print(colored(f'\t{arg}', 'green'))


class Cavity:
    """
    Command Line Interface module for running analysis.

    .. note::

       Still under development so some functions might not work properly
    """

    def __init__(self, n_cells, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none', name='cavity'):
        """
        Initialise cavity object. A cavity object is defined by the number of cells, the cell geometric parameters,
        if it has beampipes or not and the name. These properties could be changed and retrieved later using the
        corresponding ``set`` and ``get`` functions.

        Parameters
        ----------
        n_cells: int
            Number of cells
        mid_cell: list, array like
            Mid cell geometric parameters of the cavity
        end_cell_left: list, array like
            Left end cell geometric parameters of the cavity
        end_cell_right: list, array like
            Right end cell geometric parameters of the cavity
        beampipe: {'none', 'both', 'left', 'right'}
            Beampipe options
        name
        """
        self.n_modes = n_cells + 1
        self.n_modules = 1
        self.folder = None
        self.bc = 33
        self.name = name
        self.n_cells = n_cells
        self.mid_cell = mid_cell
        self.end_cell_left = end_cell_left
        self.end_cell_right = end_cell_right
        self.beampipe = beampipe
        self.no_of_modules = 1
        self.slans_qois = {}
        self.custom_eig_qois = {}
        self.abci_qois = {}
        self.slans_tune_res = {}
        self.wake_op_points = {}

        # slans results
        self.R_Q, self.k_fm, self.GR_Q, self.op_freq, self.e, self.b, \
        self.G, self.ff, self.k_cc, self.axis_field, self.surface_field = [0 for _ in range(11)]

        # abci results
        self.k_fm, self.k_loss, self.k_kick, self.phom, self.sigma, self.I0 = [0 for _ in range(6)]

        self.wall_material = None

        self.freq = None

        if not (isinstance(end_cell_left, np.ndarray) or isinstance(end_cell_left, list)):
            end_cell_left = mid_cell

        if not (isinstance(end_cell_right, np.ndarray) or isinstance(end_cell_right, list)):
            if not (isinstance(end_cell_left, np.ndarray) or isinstance(end_cell_left, list)):
                end_cell_right = mid_cell
            else:
                end_cell_right = end_cell_left

        self.end_cell_left = end_cell_left
        self.end_cell_right = end_cell_right

        self.A, self.B, self.a, self.b, self.Ri, self.L, self.Req = self.mid_cell
        self.A_el, self.B_el, self.a_el, self.b_el, self.Ri_el, self.L_el, self.Req_el = self.end_cell_left
        self.A_er, self.B_er, self.a_er, self.b_er, self.Ri_er, self.L_er, self.Req_er = self.end_cell_right

        # get geometric parameters
        self.shape_space = {
            "IC": self.mid_cell,
            "OC": self.end_cell_left,
            "OC_R": self.end_cell_right,
            "BP": beampipe
        }

    def set_name(self, name):
        """
        Set cavity name

        Parameters
        ----------
        name: str
            Name of cavity

        Returns
        -------

        """
        self.name = name

    def set_n_cells(self, n_cells):
        """
        Sets number of cells of cavity

        Parameters
        ----------
        n_cells: int
            Number of cavity cells

        Returns
        -------

        """
        self.n_cells = int(n_cells)

    def set_mid_cell(self, cell):
        """
        Set mid cell geometric parameters of cavity

        Parameters
        ----------
        cell: list, array like
            Geometric parameters of cells

        Returns
        -------

        """
        self.mid_cell = cell

    def set_end_cell_left(self, cell):
        """
        Set left end cell geometric parameters of cavity

        Parameters
        ----------
        cell: list, array like
            Geometric parameters of cells

        Returns
        -------

        """
        self.end_cell_left = cell

    def set_end_cell_right(self, cell):
        """
        Set right end cell geometric parameters of cavity

        Parameters
        ----------
        cell: list, array like
            Geometric parameters of cells

        Returns
        -------

        """
        self.end_cell_right = cell

    def set_boundary_conditions(self, bc):
        """
        Sets boundary conditions for the beampipes of cavity

        Parameters
        ----------
        bc: int
            Boundary condition of left and right cell/beampipe ends

        Returns
        -------

        """
        self.bc = int(bc)

    def set_beampipe(self, bp):
        """
        Set beampipe option of cavity

        Parameters
        ----------
        bp: str
            Beampipe option of cell

        Returns
        -------

        """
        self.beampipe = bp

    def load(self):
        """
        Load existing cavity project folder

        Parameters
        ----------

        Returns
        -------

        """
        pass

    def save(self, files_path):
        """
        Set folder to save cavity analysis results

        Parameters
        ----------
        files_path: str
            Save project directory

        Returns
        -------

        """

        if files_path is None:
            error('Please specify a folder to write the simulation results to.')
            return
        else:
            try:
                self.folder = files_path
                success = self._create_project()
                if not success:
                    error("Project could not be created. Please check the folder and try again.")
                    return
                else:
                    done("Project created successfully/already exists. You can proceed with analysis.")
            except Exception as e:
                error("Exception occurred: ", e)
                return

        if files_path is None:
            self.folder = os.getcwd()

    def load_shape_space(self, filepath):
        """
        Get cavity geometric parameters from shape space

        Parameters
        ----------
        filepath: str
            Shape space directory

        Returns
        -------

        """
        pass

    def save_shape_space(self, filepath=None):
        """
        Save current geometric parameters as shape space

        Parameters
        ----------
        filepath: str
            Directory to save shape space to. If no input is given, it is saved to the Cavities directory

        Returns
        -------

        """
        pass

    def run_tune(self, tune_variable, cell_type='Mid Cell', freq=None, solver='SLANS', proc=0, resume=False, n_cells=1):
        """
        Tune current cavity geometry

        Parameters
        ----------
        n_cells: int
            Number of cells used for tuning.
        resume: bool
            Option to resume tuning or not. Only for shape space with multiple entries.
        proc: int
            Processor number
        solver: {'SLANS', 'Native'}
            Solver to be used. Native solver is still under development. Results are not as accurate as that of SLANS.
        freq: float
            Reference frequency in MHz
        cell_type: {'mid cell', 'end-mid cell', 'mid-end cell', 'single cell'}
            Type of cell to tune
        tune_variable: {'Req', 'L'}
            Tune variable. Currently supports only the tuning of the equator radius ``Req`` and half-cell length ``L``

        Returns
        -------

        """

        for _ in tqdm([1]):
            iter_set = ['Linear Interpolation', 1e-4, 10]

            if freq is None:
                # calculate freq from mid cell length
                beta = 1
                freq = beta * c0 / (4 * self.mid_cell[5])
                info("Calculated freq from mid cell half length: ", freq)

            # create new shape space based on cell type
            if cell_type.lower() == 'mid cell':
                shape_space = {
                    f'{self.name}':
                        {
                            'IC': self.shape_space['IC'],
                            'OC': self.shape_space['IC'],
                            'OC_R': self.shape_space['IC'],
                            "BP": 'none',
                            'FREQ': freq
                        }
                }
            elif cell_type.lower() == 'end-mid cell' or cell_type.lower() == 'mid-end cell':
                shape_space = {
                    f'{self.name}':
                        {
                            'IC': self.shape_space['IC'],
                            'OC': self.shape_space['OC'],
                            'OC_R': self.shape_space['OC'],
                            "BP": 'left',
                            'FREQ': freq
                        }
                }
            elif cell_type.lower() == 'end-end cell':
                shape_space = {
                    f'{self.name}':
                        {
                            'IC': self.shape_space['IC'],
                            'OC': self.shape_space['IC'],
                            'OC_R': self.shape_space['IC'],
                            "BP": 'left',
                            'FREQ': freq
                        }
                }
            elif cell_type.lower() == 'single cell':
                shape_space = {
                    f'{self.name}':
                        {
                            'IC': self.shape_space['IC'],
                            'OC': self.shape_space['IC'],
                            'OC_R': self.shape_space['IC'],
                            "BP": 'both',
                            'FREQ': freq
                        }
                }
            else:
                shape_space = {
                    f'{self.name}':
                        {
                            'IC': self.shape_space['IC'],
                            'OC': self.shape_space['IC'],
                            'OC_R': self.shape_space['IC'],
                            "BP": 'none',
                            'FREQ': freq
                        }
                }

            if len(self.slans_tune_res.keys()) != 0:
                run_tune = input("This cavity has already been tuned. Run tune again? (y/N)")
                if solver.lower() == 'slans':
                    if run_tune.lower() == 'y':
                        # copy files required for simulation
                        self._overwriteFolder(proc, self.folder, self.name)
                        self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                        self.run_tune_slans(shape_space, resume, proc, self.bc,
                                            SOFTWARE_DIRECTORY, self.folder, self.name, tuner,
                                            tune_variable, iter_set, cell_type,
                                            progress_list=[], convergence_list=[], n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                        self.get_slans_qois()
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
            else:
                if solver.lower() == 'slans':
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_slans(shape_space, resume, proc, self.bc,
                                        SOFTWARE_DIRECTORY, fr'{self.folder}', self.name, tuner,
                                        tune_variable, iter_set, cell_type,
                                        progress_list=[], convergence_list=[], n_cells=n_cells)
                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                        self.get_slans_qois()
                    except FileNotFoundError:
                        error("OOps! Something went wrong. Could not find the tune results. Please run tune again.")

    @staticmethod
    def run_tune_slans(shape, resume, p, bc, parentDir, projectDir, filename, tuner_option,
                       tune_variable, iter_set, cell_type, progress_list, convergence_list, n_cells):
        tuner.tune(shape, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                   tuner_option=tuner_option, tune_variable=tune_variable, iter_set=iter_set,
                   cell_type=cell_type,
                   progress_list=progress_list, convergence_list=convergence_list,
                   save_last=True,
                   n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    def run_eigenmode(self, solver='SLANS', freq_shift=0, boundary_cond=None, subdir='',
                      UQ=False):
        """
        Run eigenmode analysis on cavity

        Parameters
        ----------
        solver: {'SLANS', 'Native'}
            Solver to be used. Native solver is still under development. Results are not as accurate as that of SLANS.
        freq_shift:
            Frequency shift. Eigenmode solver searches for eigenfrequencies around this value
        boundary_cond: int
            Boundary condition of left and right cell/beampipe ends
        subdir: str
            Sub directory to save results to
        UQ: bool
            Used to turn on or off uncertainty quantification

        Returns
        -------

        """

        for _ in tqdm([1]):
            if boundary_cond:
                self.bc = boundary_cond

            if solver.lower() == 'slans':
                self._run_slans(self.name, self.n_cells, self.n_modules, self.shape_space, self.n_modes, freq_shift,
                                self.bc, SOFTWARE_DIRECTORY, self.folder, sub_dir='', UQ=UQ)

                # load quantities of interest
                try:
                    self.get_slans_qois()
                except FileNotFoundError:
                    error("Could not find the tune results. Please rerun eigenmode analysis.")
            else:
                self._run_custom_eig(self.name, self.folder, self.n_cells,
                                     self.mid_cell, self.end_cell_left, self.end_cell_right,
                                     beampipe=self.beampipe, plot=False)
                # load quantities of interest
                self.get_custom_eig_qois()

    def run_wakefield(self, MROT=2, MT=10, NFS=10000, wakelength=50, bunch_length=25,
                      DDR_SIG=0.1, DDZ_SIG=0.1, WG_M=None, marker='', wp_dict=None, solver='ABCI'):
        """
        Run wakefield analysis on cavity

        Parameters
        ----------
        MROT: {0, 1}
            Polarisation 0 for longitudinal polarization and 1 for transversal polarization
        MT: int
            Number of time steps it takes for a beam to move from one mesh cell to the other
        NFS: int
            Number of frequency samples
        wakelength:
            Wakelength to be analysed
        bunch_length: float
            Length of the bunch
        DDR_SIG: float
            Mesh to bunch length ration in the r axis
        DDZ_SIG: float
            Mesh to bunch length ration in the z axis
        WG_M:
            For module simulation. Specifies the length of the beampipe between two cavities.
        marker: str
            Marker for the cavities. Adds this to the cavity name specified in a shape space json file
        wp_dict: dict
            Python dictionary containing relevant parameters for the wakefield analysis for a specific operating point
        solver: {'ABCI'}
            Only one solver is currently available

        Returns
        -------

        """

        if wp_dict is None:
            wp_dict = {}

        if len(wp_dict.keys()) != 0:
            self.wake_op_points = wp_dict

        exist = False
        # check if folders already exists
        if os.path.exists(fr'{self.folder}'):
            exist = True
            if len(wp_dict.keys()) != 0:

                for wp, vals in wp_dict.items():
                    self.sigma = vals['sigma_z (SR/BS) [mm]']

                    wp_SR = f"{wp}_SR_{self.sigma.split(r'/')[0]}mm"
                    wp_BS = f"{wp}_SR_{self.sigma.split(r'/')[1]}mm"
                    if os.path.exists(fr"{self.folder}/SimulationData/ABCI/{self.name}/{wp_SR}") and \
                            os.path.exists(fr"{self.folder}/SimulationData/ABCI/{self.name}/{wp_BS}"):
                        pass
                    else:
                        exist = False
            else:
                exist = False
        else:
            exist = False

        if exist:
            run_wake = input("Wakefield results already exist for these settings. Run wakefield again? (y/N)")
            if run_wake.lower() == 'y':
                msg = True
            else:
                msg = False

        if not exist:
            for _ in tqdm([1]):
                if solver == 'ABCI':
                    info(">> Running wakefield simulation")
                    self._run_abci(self.name, self.n_cells, self.n_modules, self.shape_space,
                                   MROT=MROT, MT=MT, NFS=NFS, UBT=wakelength, bunch_length=bunch_length,
                                   DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG,
                                   parentDir=SOFTWARE_DIRECTORY, projectDir=fr'{self.folder}', WG_M=WG_M, marker=marker,
                                   wp_dict=wp_dict, freq=self.freq, R_Q=self.R_Q)

                    try:
                        self.get_abci_qois()
                    except FileNotFoundError:
                        error("Could not find the abci wakefield results. Please rerun wakefield analysis.")
        else:
            try:
                self.get_abci_qois()
            except FileNotFoundError:
                error("Could not find the abci wakefield results. Please rerun wakefield analysis.")

    def calc_op_freq(self):
        """
        Calculates operating frequency. The operating frequency is used for tuning when a frequency is not given.
        It is advisable to always include the desired tune frequency. Example

        .. py:function:: cav.run_tune('Req', freq=1300)

        Returns
        -------

        """
        if not self.freq:
            self.freq = (c0 / 4 * self.L)

    @staticmethod
    def _run_custom_eig(name, folder, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe='both', plot=False):
        custom_eig.set_name(name)
        custom_eig.set_folder(folder)
        custom_eig.run(n_cells, mid_cell, end_cell_left, end_cell_right, beampipe='both', plot=True)

    @staticmethod
    def _run_slans(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir, sub_dir='', UQ=False):
        start_time = time.time()
        # create folders for all keys
        slans_geom.createFolder(name, projectDir, subdir=sub_dir)

        try:
            slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                              n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
        except KeyError:
            slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                              n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

        # run UQ
        if UQ:
            uq(name, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
               n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
               f_shift=f_shift, bc=bc, parentDir=parentDir, projectDir=projectDir)

        done(f'Done with Cavity {name}. Time: {time.time() - start_time}')

    @staticmethod
    def _run_abci(name, n_cells, n_modules, shape, MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
                  DDR_SIG=0.1, DDZ_SIG=0.1,
                  parentDir=None, projectDir=None,
                  WG_M=None, marker='', wp_dict=None, freq=0, R_Q=0):

        # run abci code
        if WG_M is None:
            WG_M = ['']

        start_time = time.time()
        # run both polarizations if MROT == 2
        for ii in WG_M:
            try:
                if MROT == 2:
                    for m in tqdm(range(2)):
                        abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                         fid=name, MROT=m, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                         projectDir=projectDir,
                                         WG_M=ii, marker=ii)

                else:
                    abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                     fid=name, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                     DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir,
                                     WG_M=ii, marker=ii)
            except KeyError:
                if MROT == 2:
                    for m in tqdm(range(2)):
                        abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                         fid=name, MROT=m, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                         projectDir=projectDir,
                                         WG_M=ii, marker=ii)
                else:
                    abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                     fid=name, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                     DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir,
                                     WG_M=ii, marker=ii)

        done(f'Cavity {name}. Time: {time.time() - start_time}')
        if len(wp_dict.keys()) > 0:
            try:
                if freq != 0 and R_Q != 0:
                    d = {}
                    # save qois
                    for key, vals in tqdm(wp_dict.items()):
                        WP = key
                        I0 = float(vals['I0 [mA]'])
                        Nb = float(vals['Nb [1e11]'])
                        sigma_z = [float(x) for x in vals['sigma_z (SR/BS) [mm]'].split('/')]
                        bl_diff = ['SR', 'BS']

                        info("Running wakefield analysis for given operating points.")
                        for i, s in enumerate(sigma_z):
                            for ii in WG_M:
                                fid = f"{WP}_{bl_diff[i]}_{s}mm{ii}"
                                try:
                                    for m in range(2):
                                        abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                                         fid=fid, MROT=m, MT=MT, NFS=NFS, UBT=10 * s * 1e-3,
                                                         bunch_length=s,
                                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                                         projectDir=projectDir,
                                                         WG_M=ii, marker=ii, sub_dir=f"{name}")
                                except KeyError:
                                    for m in range(2):
                                        abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                                         fid=fid, MROT=m, MT=MT, NFS=NFS, UBT=10 * s * 1e-3,
                                                         bunch_length=s,
                                                         DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                                         projectDir=projectDir,
                                                         WG_M=ii, marker=ii, sub_dir=f"{name}")

                                dirc = fr'{projectDir}\SimulationData\ABCI\{name}{marker}'
                                # try:
                                k_loss = abs(ABCIData(dirc, f'{fid}', 0).loss_factor['Longitudinal'])
                                k_kick = abs(ABCIData(dirc, f'{fid}', 1).loss_factor['Transverse'])
                                # except:
                                #     k_loss = 0
                                #     k_kick = 0

                                d[fid] = get_qois_value(freq, R_Q, k_loss, k_kick, s, I0, Nb, n_cells)

                    # save qoi dictionary
                    run_save_directory = fr'{projectDir}\SimulationData\ABCI\{name}{marker}'
                    with open(fr'{run_save_directory}\qois.json', "w") as f:
                        json.dump(d, f, indent=4, separators=(',', ': '))

                    done("Done with the secondary analysis for working points")
                else:
                    info("To run analysis for working points, eigenmode simulation has to be run first"
                         "to obtain the cavity operating frequency and R/Q")
            except KeyError:
                error('The working point entered is not valid. See below for the proper input structure.')
                show_valid_operating_point_structure()

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

    def set_wall_material(self, wm):
        self.wall_material = wm

    def get_slans_tune_res(self, tune_variable, cell_type):
        """

        Parameters
        ----------
        tune_variable: {'Req', 'L'}
            Tune variable. Currently supports only the tuning of the equator radius ``Req`` and half-cell length ``L``
        cell_type: {'mid cell', 'end-mid cell', 'mid-end cell', 'single cell'}
            Type of cell to tune

        Returns
        -------

        """
        tune_res = 'tune_res.json'
        if os.path.exists(fr"{self.folder}\SimulationData\SLANS\{self.name}\{tune_res}"):
            with open(fr"{self.folder}\SimulationData\SLANS\{self.name}\{tune_res}", 'r') as json_file:
                self.slans_tune_res = json.load(json_file)

            self.freq = self.slans_tune_res['freq']

            if tune_variable.lower() == 'req':
                self.shape_space['IC'][6] = self.slans_tune_res['Req']
                self.shape_space['OC'][6] = self.slans_tune_res['Req']
                self.shape_space['OC_R'][6] = self.slans_tune_res['Req']
                self.mid_cell[6] = self.slans_tune_res['Req']
                self.end_cell_left[6] = self.slans_tune_res['Req']
                self.end_cell_right[6] = self.slans_tune_res['Req']
            else:
                self.shape_space['IC'][5] = self.slans_tune_res['Req']
                self.shape_space['OC'][5] = self.slans_tune_res['Req']
                self.shape_space['OC_R'][5] = self.slans_tune_res['Req']
                self.mid_cell[5] = self.slans_tune_res['Req']
                self.end_cell_left[5] = self.slans_tune_res['Req']
                self.end_cell_right[5] = self.slans_tune_res['Req']

            # set alpha
            if len(self.mid_cell) == 7:
                if isinstance(self.shape_space['IC'], np.ndarray):
                    self.shape_space['IC'] = np.append(self.shape_space['IC'], self.slans_tune_res['alpha_i'])
                    self.mid_cell = np.append(self.mid_cell, self.slans_tune_res['alpha_i'])
                elif isinstance(self.shape_space['IC'], list):
                    self.shape_space['IC'].append(self.slans_tune_res['alpha_i'])
                    self.mid_cell.append(self.slans_tune_res['alpha_i'])
            elif len(self.mid_cell['IC']) == 8:
                self.shape_space['IC'][7] = self.slans_tune_res['alpha_i']
                self.mid_cell[7] = self.slans_tune_res['alpha_i']

            if cell_type.lower() == 'end-mid cell' \
                    or cell_type.lower() == 'mid-end cell' \
                    or cell_type.lower() == 'single cell':

                if len(self.end_cell_left) == 7:
                    if isinstance(self.shape_space['OC'], np.ndarray):
                        self.shape_space['OC'].append(self.slans_tune_res['alpha_o'])
                        self.end_cell_left.append(self.slans_tune_res['alpha_o'])
                    elif isinstance(self.shape_space['OC'], list):
                        self.shape_space['OC'].append(self.slans_tune_res['alpha_o'])
                        self.end_cell_left.append(self.slans_tune_res['alpha_o'])

                elif len(self.end_cell_left) == 8:
                    self.shape_space['OC'][7] = self.slans_tune_res['alpha_o']
                    self.end_cell_left[7] = self.slans_tune_res['alpha_o']

            # expand tune to be able to tune right cavity geometry also
        else:
            error("Tune results not found. Please tune the cavity")

    def get_slans_qois(self):
        """
        Get quantities of interest written by the SLANS code
        Returns
        -------

        """
        qois = 'qois.json'

        with open(fr"{self.folder}\SimulationData\SLANS\{self.name}\{qois}") as json_file:
            self.slans_qois = json.load(json_file)

        self.freq = self.slans_qois['freq [MHz]']
        self.k_cc = self.slans_qois['kcc [%]']
        self.ff = self.slans_qois['ff [%]']
        self.R_Q = self.slans_qois['R/Q [Ohm]']
        self.GR_Q = self.slans_qois['GR/Q [Ohm^2]']
        self.G = self.GR_Q / self.R_Q
        # self.Q = d_qois['Q []']
        self.e = self.slans_qois['Epk/Eacc []']
        self.b = self.slans_qois['Bpk/Eacc [mT/MV/m]']

        # get axis field
        self.axis_field = fr.txt_reader(fr"{self.folder}\SimulationData\SLANS\{self.name}\cavity_33_{self.n_cells}.af",
                                        ' ')
        # get surface field
        self.surface_field = fr.txt_reader(
            fr"{self.folder}\SimulationData\SLANS\{self.name}\cavity_33_{self.n_cells}.sf", ' ')
        # print(self.surface_field)

    def get_custom_eig_qois(self):
        """
        Get quantities of interest written by the native eigenmode analysis code
        Returns
        -------

        """
        qois = 'qois.json'

        with open(fr"{self.folder}\SimulationData\NativeEig\{self.name}\{qois}") as json_file:
            self.custom_eig_qois = json.load(json_file)

    def get_abci_qois(self, opt='SR'):
        """
        Get the quantities of interest written by the ABCI code

        Parameters
        ----------
        opt: {'SR', 'BS'}
            SR - Synchrotron radiation bunch length
            BS - Bremsstrahlung

        Returns
        -------

        """
        qois = 'qois.json'

        with open(fr"{self.folder}\SimulationData\ABCI\{self.name}\{qois}") as json_file:
            self.abci_qois = json.load(json_file)

        if len(self.wake_op_points.keys()) != 0 and self.freq != 0 and self.R_Q != 0:
            for wp, vals in self.wake_op_points.items():
                self.sigma = vals['sigma_z (SR/BS) [mm]']

                if opt == 'SR':
                    d_qois = self.abci_qois[f"{wp}_{opt}_{self.sigma.split(r'/')[0]}mm"]

                    self.k_fm = d_qois['k_FM [V/pC]']
                    self.k_loss = d_qois['|k_loss| [V/pC]']
                    self.k_kick = d_qois['|k_kick| [V/pC/m]']
                    self.phom = d_qois['P_HOM [kW]']
                    self.I0 = d_qois['I0 [mA]']
                else:
                    d_qois = self.abci_qois[f"{wp}_{opt}_{self.sigma.split(r'/')[1]}mm"]

                    self.k_fm = d_qois['k_FM [V/pC]']
                    self.k_loss = d_qois['|k_loss| [V/pC]']
                    self.k_kick = d_qois['|k_kick| [V/pC/m]']
                    self.phom = d_qois['P_HOM [kW]']
                    self.I0 = d_qois['I0 [mA]']

    def _create_project(self):
        project_name = self.name
        project_dir = self.folder

        if project_name != '':

            # check if folder already exist
            e = self._check_if_path_exists(project_dir, project_name)

            if e:
                def make_dirs_from_dict(d, current_dir=fr"{project_dir}"):
                    for key, val in d.items():
                        os.mkdir(os.path.join(current_dir, key))
                        if type(val) == dict:
                            make_dirs_from_dict(val, os.path.join(current_dir, key))

                # create project structure in folders
                project_dir_structure = {
                    f'{project_name}':
                        {
                            'Cavities': None,
                            'OperatingPoints': None,
                            'SimulationData': {
                                'SLANS': None,
                                'NativeEig': None,
                                'ABCI': None,
                                'CavitiesAnalysis': None
                            },
                            'PostprocessingData': {
                                'Plots': None,
                                'Data': None,
                                'CSTData': None
                            },
                            'Reference': None
                        }
                }
                try:
                    make_dirs_from_dict(project_dir_structure)
                    self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
                    return True
                except Exception as e:
                    self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
                    print("An exception occurred in created project: ", e)
                    return False
            else:
                self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
                return True
        else:
            print('\tPlease enter a valid project name')
            self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
            return False

    @staticmethod
    def _check_if_path_exists(directory, folder):
        path = f"{directory}/{folder}"
        if os.path.exists(path):
            x = input("Project already exists. Do you want to overwrite? y/N")

            if x.lower() == 'yes' or x.lower() == 'y':
                try:
                    directory_list = os.listdir(path)

                    if 'Cavities' in directory_list \
                            and 'PostprocessingData' in directory_list \
                            and 'SimulationData' in directory_list and len(directory_list) < 6:
                        shutil.rmtree(path)
                        return True
                    else:
                        print('\tIt seems that the folder specified is not a cavity project folder. Please check folder'
                              'again to avoid deleting important files.')
                        return False

                except Exception as e:
                    print("Exception occurred: ", e)
                    return False
            else:
                return False
        else:
            return True

    @staticmethod
    def _overwriteFolder(invar, projectDir, name):
        path = fr"{projectDir}\SimulationData\SLANS\_process_{invar}"
        if os.path.exists(path):
            shutil.rmtree(path)
            dir_util._path_created = {}

        os.makedirs(path)

    @staticmethod
    def _copyFiles(invar, parentDir, projectDir, name):
        src = fr"{parentDir}\exe\SLANS_exe"
        dst = fr"{projectDir}\SimulationData\SLANS\_process_{invar}\SLANS_exe"

        dir_util.copy_tree(src, dst)


class Cavities:
    """
    Cavities object is an object containing several Cavity objects.
    """

    def __init__(self, cavities_list=None, save_folder='None'):
        """Constructs all the necessary attributes of the Cavity object

        Parameters
        ----------
        cavities_list: list, array like
            List containing Cavity objects.

        save_folder: str
            Folder to save generated images, latex, excel, and text files.
        """

        if cavities_list is None:
            self.cavities_list = []

        self.p_qois = None
        self.fm_results = None
        self.hom_results = None
        self.save_folder = save_folder

        self.cavities_list = cavities_list

        self.returned_results = None
        self.ls = ['solid', 'dashed', 'dashdot', 'dotted',
                   'solid', 'dashed', 'dashdot', 'dotted',
                   'solid', 'dashed', 'dashdot', 'dotted']

        self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
        self.set_cavities_field()

    def add_cavity(self, cav):
        """
        Adds cavity to cavities
        Parameters
        ----------
        cav: object
            Cavity object

        Returns
        -------

        """
        self.cavities_list.append(cav)

    def set_cavities_field(self):
        """
        Sets cavities analysis field range.

        Returns
        -------

        """
        for cav in self.cavities_list:
            cav.set_Eacc(Eacc=self.E_acc)

    def compare_power(self, E_acc=None):
        if E_acc is not None:
            self.E_acc = E_acc
            self.set_cavities_field()

        self.p_qois = []
        results = []
        for i, cav in enumerate(self.cavities_list):
            # E_acc_Pin(cavity, op_field[i], ls[i], fig, ax, ax_right, ax_right2)
            results.append(self.qois(cav, cav.op_field * 1e-6, E_acc))

        self.returned_results = results

    def qois(self, cavity, op_field, E_acc):
        """

        Parameters
        ----------
        cavity: object
            Cavity object
        op_field: float
            Cavity operating field

        Returns
        -------
        Dictionary containing quantities of interest (normed optional).
        """

        ind = np.where((E_acc >= 0.99 * op_field * 1e6) & (E_acc <= 1.01 * op_field * 1e6))
        qois = {
            r"N_cav/beam": np.average(cavity.n_cav[ind]),
            r"Q0 [10^8]$": np.average(cavity.Q0[ind] * 1e-8),
            r"Rs [Ohm]$": np.average(cavity.Rs[ind]),
            r"P_stat/cav [W]": np.average(cavity.pstat[ind] / cavity.n_cav[ind]),
            r"P_dyn/cav [W]": np.average(cavity.pdyn[ind] / cavity.n_cav[ind]),
            # r"P_\mathrm{wp/cav}$ [W]": np.average(cavity.p_wp[ind]/cavity.n_cav[ind]),
            r"P_in/beam [kW]": np.average(cavity.p_in[ind]) * 1e-3,
            # r"$Q_\mathrm{0} \mathrm{[10^8]}$": np.average(cavity.Q0[ind] * 1e-8),
            # r"$Rs_\mathrm{0} \mathrm{[10^7]}$": np.average(cavity.Rs[ind])
        }
        self.p_qois.append(qois)

        qois_norm_units = {
            r"$n_\mathrm{cav/beam}$": np.average(cavity.n_cav[ind]),
            # r"$q_\mathrm{0}$": np.average(cavity.Q0[ind]),
            # r"$r_\mathrm{s}$": np.average(cavity.Rs[ind]),
            r"$p_\mathrm{stat/cav}$": np.average(cavity.pstat[ind] / cavity.n_cav[ind]),
            r"$p_\mathrm{dyn/cav}$": np.average(cavity.pdyn[ind] / cavity.n_cav[ind]),
            # r"$p_\mathrm{wp/cav}$": np.average(cavity.p_wp[ind]/cavity.n_cav[ind]),
            r"$p_\mathrm{in/beam}$": np.average(cavity.p_in[ind]),
            # r"$Q_\mathrm{0} \mathrm{[10^8]}$": np.average(cavity.Q0[ind] * 1e-8),
            # r"$Rs_\mathrm{0} \mathrm{[10^7]}$": np.average(cavity.Rs[ind])
        }

        ic(qois)
        return qois_norm_units

    def qois_fm(self):
        """
        Retrieves the fundamental mode quantities of interest

        Returns
        -------
        Dictionary containing fundamental mode quantities of interest (normed optional).
        """
        results = []
        for cav in self.cavities_list:
            results.append({
                r"$E_\mathrm{pk}/E_\mathrm{acc} [\cdot]$": cav.e,
                r"$B_\mathrm{pk}/E_\mathrm{acc} \mathrm{[mT/MV/m]}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$R/Q \mathrm{[10^2\Omega]}$": cav.R_Q * 1e-2,
                r"$G \mathrm{[10^{2}\Omega]}$": cav.G * 1e-2,
                r"$G\cdot R/Q \mathrm{[10^{5}\Omega^2]}$": cav.GR_Q * 1e-5
            })

        results_norm_units = []
        for cav in self.cavities_list:
            results_norm_units.append({
                r"$e_\mathrm{pk}$": cav.e,
                r"$b_\mathrm{pk}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$r/q$": cav.R_Q,
                r"$g$": cav.G,
                r"$g\cdot r/q $": cav.GR_Q
            })
        ic(results)
        return results_norm_units

    def qois_hom(self):
        """
        Retrieves the higher-order modes quantities of interest

        Returns
        -------
        Dictionary containing higher-order modes quantities of interest (normed optional).
        """

        results = []
        for cavity in self.cavities_list:
            results.append({
                r"$|k_\parallel| \mathrm{[V/pC]}$": cavity.k_loss,
                r"$|k_\perp| \mathrm{[V/pC/m]}$": cavity.k_kick,
                r"$P_\mathrm{HOM}/cav \mathrm{[kW]}$": cavity.phom
            })

        results_norm_units = []
        for cavity in self.cavities_list:
            results_norm_units.append({
                r"$k_\parallel$": cavity.k_loss,
                r"$k_\perp$": cavity.k_kick,
                r"$p_\mathrm{HOM}/cav$": cavity.phom
            })
        ic(results)
        return results_norm_units

    def plot_power_comparison(self, fig=None, ax_list=None):
        """
        Can be called using ``cavities.plot_power_comparison()``

        .. math::

           W^{3 \\beta}_{\delta}

        Parameters
        ----------
        fig: matplotlib figure
        ax_list: list of matplotlib axes object

        Returns
        -------

        """
        if fig is not None:
            fig = fig
            ax1, ax2, ax3 = ax_list
        else:
            # create figure
            fig = plt.figure()
            gs = fig.add_gridspec(2, 2)
            ax1 = fig.add_subplot(gs[:, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, 1])
        # def E_acc_Pin(self, cavity, E_acc, op_field, ls='-', ):

        # ax_right2._get_lines.prop_cycler = ax._get_lines.prop_cycler
        # ax_right2.spines["right"].set_position(("axes", 1.2))
        for i, cavity in enumerate(self.cavities_list):
            ax1.plot(cavity.E_acc * 1e-6, cavity.pstat / cavity.n_cav,
                     ls=self.ls[i], lw=2, c='tab:orange',
                     label=r"$P_\mathrm{static/cav}$" + fr"{cavity.name}")

            ax1.plot(cavity.E_acc * 1e-6, cavity.pdyn / cavity.n_cav,
                     ls=self.ls[i], lw=2, c='tab:blue', label=r"$P_\mathrm{dynamic/cav}$" + fr"{cavity.name}")

            # p1, = ax1.plot(cavity.E_acc * 1e-6, cavity.p_wp/cavity.n_cav,
            #                ls=self.ls[i], lw=2, c='k', label=r"$P_\mathrm{wp/beam}$" + fr"{cavity.name}")

            p2, = ax2.plot(cavity.E_acc * 1e-6, cavity.n_cav, ls=self.ls[i], lw=2, c='tab:red',
                           label=fr"{cavity.name}")

            p3, = ax3.plot(cavity.E_acc * 1e-6, cavity.p_in * 1e-3, ls=self.ls[i], lw=2, c='tab:purple',
                           label=fr"{cavity.name}")

            ax1.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax1.set_ylabel(r"$P_\mathrm{stat, dyn}$/cav [W]")
            ax2.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax2.set_ylabel(r"$N_\mathrm{cav/beam}$")
            ax3.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax3.set_ylabel(r"$P_\mathrm{in/cav}$ [kW]")

            ax1.axvline(cavity.op_field * 1e-6, ls=':', c='k')
            ax1.text(cavity.op_field * 1e-6 - 1, 0.3, f"{cavity.op_field * 1e-6} MV/m",
                     size=14, rotation=90, transform=ax1.get_xaxis_transform())
            ax2.axvline(cavity.op_field * 1e-6, ls=':', c='k')
            ax2.text(cavity.op_field * 1e-6 - 1, 0.5, f"{cavity.op_field * 1e-6} MV/m",
                     size=14, rotation=90, transform=ax2.get_xaxis_transform())
            ax3.axvline(cavity.op_field * 1e-6, ls=':', c='k')
            ax3.text(cavity.op_field * 1e-6 - 1, 0.3, f"{cavity.op_field * 1e-6} MV/m",
                     size=14, rotation=90,
                     transform=ax3.get_xaxis_transform())

            # ax.axvline(7.13, ls='--', c='k')
            # ax.axvline(10, ls='--', c='k')
            # ax.axvline(15, ls='--', c='k')
            # ax_right2.axhline(500, ls='--', c='k')
            # ax_right2.axhline(1000, ls='--', c='k')
            # ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
            # ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
            # ax_right.yaxis.set_major_locator(ticker.MultipleLocator(100))
            # ax_right2.yaxis.set_major_locator(ticker.MultipleLocator(200))

            # ax.yaxis.label.set_color(p1.get_color())
            # ax_right.yaxis.label.set_color(p2.get_color())
            # ax_right2.yaxis.label.set_color(p3.get_color())

            ax1.set_xlim(min(cavity.E_acc) * 1e-6, max(cavity.E_acc) * 1e-6)
            ax2.set_xlim(min(cavity.E_acc) * 1e-6, max(cavity.E_acc) * 1e-6)
            ax3.set_xlim(min(cavity.E_acc) * 1e-6, max(cavity.E_acc) * 1e-6)
            # # ax.set_ylim(0, 50)
            # ax_right.set_ylim(100, 400)f
            # ax_right2.set_ylim(0, 700)
            ax1.set_yscale('log')
            ax2.set_yscale('log')
            ax3.set_yscale('log')

            # tkw = dict(size=4, width=1.5)
            # ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
            # ax_right.tick_params(axis='y', colors=p2.get_color(), **tkw)
            # ax_right2.tick_params(axis='y', colors=p3.get_color(), **tkw)
            # ax.tick_params(axis='x', **tkw)

            ax1.minorticks_on()
            mplcursors.cursor(ax1)
            mplcursors.cursor(ax2)
            mplcursors.cursor(ax3)
            # ax.grid(True, which='both', axis='both')

        # dummy lines with NO entries, just to create the black style legend
        dummy_lines = []
        for b_idx, b in enumerate(self.cavities_list):
            dummy_lines.append(ax1.plot([], [], c="gray", ls=self.ls[b_idx])[0])

        lines = ax1.get_lines()
        legend1 = ax1.legend([lines[i] for i in range(3)],
                             [r"$P_\mathrm{stat}$", r"$P_\mathrm{dyn}$"], loc=3)
        legend2 = ax1.legend([dummy_lines[i] for i in range(len(self.cavities_list))],
                             [cavity.name for cavity in self.cavities_list],
                             loc=0)
        ax1.add_artist(legend1)

        # ax1.legend(ncol=len(cavities))
        ax2.legend(loc='upper left')
        ax3.legend(loc=3)

        label = [r"$\mathbf{Z^*}$", 'Z', r"$\mathbf{W^*}$", 'W']
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_power_comparison.png")

        plt.show()

    def plot_compare_bar(self):
        """
        Plots bar chart of power quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 3)
        # plot barchart
        data = np.array([list(d.values()) for d in self.returned_results])
        data_col_max = data.max(axis=0)

        x = list(self.returned_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = 0.15  # 1 / len(x)
        for i, cav in enumerate(self.cavities_list):
            print(cav.name)
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=cav.name)

        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]

        ax.axhline(1.05, c='k')
        ax.set_ylim(-0.01, 1.5 * ax.get_ylim()[-1])
        ax.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_power_comparison_bar.png")

        plt.show()

    def plot_compare_hom_bar(self):
        """
        Plot bar chart of higher-order mode's quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 3)
        # plot barchart
        self.hom_results = self.qois_hom()
        data = np.array([list(d.values()) for d in self.hom_results])
        data_col_max = data.max(axis=0)
        x = list(self.hom_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = 0.15  # 1 / (len(x)+10)
        for i, cav in enumerate(self.cavities_list):
            print(type(X), type(i), type(width), type(data[i]), data)
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=cav.name)

        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        # label = ["C3795_ttbar (5-Cell)", "FCCUROS5_ttbar (5-Cell)", "TELSA_ttbar (5-Cell)"]
        # ax.legend(label, loc="upper left")
        ax.axhline(1.05, c='k')
        ax.set_ylim(-0.01, 1.5 * ax.get_ylim()[-1])
        ax.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_hom_bar.png")

        plt.show()

    def plot_compare_fm_bar(self):
        """
        Plot bar chart of fundamental mode quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 3)
        # plot barchart
        self.fm_results = self.qois_fm()
        data = np.array([list(d.values()) for d in self.fm_results])
        data_col_max = data.max(axis=0)
        x = list(self.fm_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = 0.15  # 1 / (len(x)+10)
        for i, cav in enumerate(self.cavities_list):
            print(type(X), type(i), type(width), type(data[i]), data)
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=self.cavities_list[i].name)

        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        # label = ["C3795_ttbar (5-Cell)", "FCCUROS5_ttbar (5-Cell)", "TELSA_ttbar (5-Cell)"]
        # ax.legend(label, loc="upper left")
        ax.axhline(1.05, c='k')
        ax.set_ylim(-0.01, 1.5 * ax.get_ylim()[-1])
        ax.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_fm_bar.png")

        plt.show()

    def plot_ql_vs_pin(self):
        """
        Plot loaded quality factor versus input power for fundamental power coupler (FPC)

        Returns
        -------

        """
        label = [cav.name for cav in self.cavities_list]

        # geometry
        n_cells = [cav.n_cells for cav in self.cavities_list]  # 4
        l_cell = [cav.l_cell_mid for cav in self.cavities_list]  # m
        G = [cav.G for cav in self.cavities_list]  # 273.2 # 0.00948*Q_factor*(f/1300)**0.5
        b = [cav.b for cav in self.cavities_list]
        geometry = [n_cells, l_cell, G, b]

        # QOI
        R_Q = [cav.R_Q for cav in self.cavities_list]  # 411   # c Ohm linac definition
        f0 = [cav.op_freq for cav in self.cavities_list]
        QOI = [f0, R_Q]

        # RF
        # Vrf = [2*0.1e9, 2*0.1e9, 2*0.75e9]  #   #2*0.75e9
        # Vrf = [2*0.12e9, 2*0.12e9, 2*1e9, 2*0.44e9]  #   #2*0.75e9 update
        Vrf = [cav.v_rf for cav in self.cavities_list]  # #ttbar
        # Eacc = [20e6, 20e6, 20e6]
        Eacc = [cav.op_field for cav in self.cavities_list]  # update
        RF = [Eacc, Vrf]

        # MACHINE
        # I0 = [1390e-3, 1390e-3, 147e-3, 147e-3]  # mA
        I0 = [WP[cav.wp]['I0 [mA]'] * 1e-3 for cav in self.cavities_list]
        # rho = [10.76e3, 10.76e3, 10.76e3, 10.76e3]  # bending radius
        rho = [MACHINE['rho [m]'] for cav in self.cavities_list]  # bending radius
        E0 = [WP[cav.wp]['E [GeV]'] for cav in self.cavities_list]  # Beam energy GeV
        machine = [I0, rho, E0]

        self.ql_pin(label, geometry, RF, QOI, machine)

    def plot_cryomodule_comparison(self):
        """
        Plot cryomodule power comparison

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (9, 3)
        fig, axs = plt.subplots(1, 2)
        n_cav_per_cryomodule = np.arange(1, 11)
        for cav in self.cavities_list:
            n_cryomodules_list = []
            cryomodules_len_list = []
            for ncpc in n_cav_per_cryomodule:
                cryomodule_len = cav.n_cells * (2 * cav.l_cell_mid) * ncpc + (ncpc + 1) * 8 * cav.l_cell_mid
                cryomodules_len_list.append(cryomodule_len)

                n_cryomodules = cav.n_cav_op_field / ncpc
                n_cryomodules_list.append(n_cryomodules)

            # ic(cryomodules_len_list)
            axs[0].plot(n_cav_per_cryomodule, cryomodules_len_list, marker='o', mec='k', label=f'{cav.name}')
            axs[1].plot(n_cav_per_cryomodule, n_cryomodules_list, marker='o', mec='k', label=f'{cav.name}')
            ic(n_cav_per_cryomodule)
            ic(cryomodules_len_list, n_cryomodules_list)
            ic(n_cav_per_cryomodule)
            ic()
        axs[0].set_xlabel("$N_\mathrm{cav}$/mod.")
        axs[0].set_ylabel("$L_\mathrm{cryo}$ [m]")
        axs[1].set_xlabel("$N_\mathrm{cav}$/mod.")
        axs[1].set_ylabel("$N_\mathrm{cryo}$")
        axs[0].legend()
        axs[1].legend()
        mplcursors.cursor(axs[0])
        mplcursors.cursor(axs[1])
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_cryo.png")

        plt.show()

    def plot_cavities_contour(self, opt='mid', n_cells=1):
        """Plot geometric contour of Cavity objects

        Parameters
        ----------
        opt: {"mid", "end", "all"}
            Either plot contour for only mid cells or end cells or the entire cavity
        n_cells: int
            Option used only when opt is set to "all"

        Returns
        -------

        """
        min_x, max_x, min_y, max_y = [], [], [], []

        if opt.lower() == 'mid' or opt.lower() == 'end':
            plt.rcParams["figure.figsize"] = (4, 5)
        else:
            plt.rcParams["figure.figsize"] = (10, 4)

        for cav in self.cavities_list:
            # write contour
            self.write_contour(cav, opt)

            data = pd.read_csv(fr"{cav.slans_dir}\contour.txt", sep=r'\s+', header=None)

            plt.plot(data[1] * 1000, data[0] * 1000, lw=3., label=cav.name)
            plt.legend(loc='lower left')

            x_label = "z [mm]"
            y_label = "r [mm]"
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            min_x.append(min(data[1]))
            min_y.append(min(data[0]))
            max_x.append(max(data[1]))
            max_y.append(max(data[0]))

        if opt.lower() == 'mid' or opt.lower() == 'end':
            plt.xlim(-0.1, max(max_x) * 1e3 + 1)
            plt.ylim(-0.1, max(max_y) * 1e3 + 1)
        else:
            plt.xlim(min(min_x) * 1e3 - 1, max(max_x) * 1e3 + 1)
            plt.ylim(min(min_y) * 1e3 - 1, max(max_y) * 1e3 + 1)

        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_contour.png")

        plt.show()

    def plot_axis_fields(self):
        """
        Plot axis fields of cavities

        Returns
        -------

        """
        for cav in self.cavities_list:
            # normalize fields
            e_axis = np.abs(cav.axis_field['1'])
            e_axis_norm = e_axis / e_axis.max()

            # shift to mid
            z = cav.axis_field['0']
            z_shift = z - z.max() / 2
            plt.plot(z_shift, e_axis_norm, label=cav.name)

        plt.xlabel('$z$ [mm]')
        plt.ylabel('$|E_\mathrm{axis}|/|E_\mathrm{axis}|_\mathrm{max}$')
        plt.axhline(1.02, c='k')
        plt.ylim(-0.01, 1.5)
        plt.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_axis_fields.png")

        plt.show()

    def plot_surface_fields(self):
        """
        Plot surface fields of cavities

        Returns
        -------

        """
        for cav in self.cavities_list:
            # normalize fields
            e_surf = np.abs(cav.surface_field['0'])
            e_surf_norm = e_surf / e_surf.max()

            plt.plot(e_surf_norm, label=cav.name)

        plt.axhline(1.02, c='k')
        plt.ylim(-0.01, 1.5)
        plt.xlabel('$L_\mathrm{surf}$ [mm]')
        plt.ylabel('$|E_\mathrm{surf}|/|E_\mathrm{surf}|_\mathrm{max}$')
        plt.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_surface_fields.png")

        plt.show()

    def plot_multipac_triplot(self, folders, kind='triplot'):
        """
        Plot Multipac triplot

        Parameters
        ----------
        folders: list, array like
            List of folder to read multipacting results from
        kind

        Notes
        -----
        This will be changed later so that the multipac results will be in the same location as the SLANS and ABCI
        results

        Returns
        -------

        """

        if kind == 'triplot':
            # create figure
            fig = plt.figure()
            gs = fig.add_gridspec(3, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[2, 0])
            axs = [ax1, ax2, ax3]

        else:
            print("in here")
            fig, axs = plt.subplots(1, 1)
            axs = [axs]

        mpl.rcParams['figure.figsize'] = [6, 10]

        Eacc_list = [cav.op_field * 1e-6 for cav in self.cavities_list]
        Epk_Eacc_list = [cav.e for cav in self.cavities_list]
        labels = [cav.name for cav in self.cavities_list]
        for Eacc, Epk_Eacc, folder, label in zip(Eacc_list, Epk_Eacc_list, folders, labels):
            # load_output_data
            # files
            fnames = ["Ccounter.mat", "Acounter.mat", "Atcounter.mat", "Efcounter.mat", "param",
                      "geodata.n", "secy1", "counter_flevels.mat", "counter_initials.mat"]
            data = {}
            # files_folder = "D:\Dropbox\multipacting\MPGUI21"
            for f in fnames:
                if ".mat" in f:
                    data[f] = spio.loadmat(fr"{folder}\\{f}")
                else:
                    data[f] = pd.read_csv(fr"{folder}\\{f}", sep='\s+', header=None)

            A = data["Acounter.mat"]["A"]
            At = data["Atcounter.mat"]["At"]
            C = data["Ccounter.mat"]["C"]
            Ef = data["Efcounter.mat"]["Ef"]
            flevel = data["counter_flevels.mat"]["flevel"]
            initials = data["counter_initials.mat"]["initials"]
            secy1 = data["secy1"].to_numpy()
            Pow = flevel
            n = len(initials[:, 0]) / 2  # number of initials in the bright set
            N = int(data["param"].to_numpy()[4])  # number of impacts
            U = flevel
            Efl = flevel
            q = 1.6021773e-19
            Efq = Ef / q

            e1 = np.min(np.where(secy1[:, 1] >= 1))  # lower threshold
            e2 = np.max(np.where(secy1[:, 1] >= 1))  # upper threshold
            val, e3 = np.max(secy1[:, 1]), np.argmax(secy1[:, 1])  # maximum secondary yield

            cl = 0
            ok, ok1, ok2 = 1, 1, 1
            if ok > 0:
                if n == 0:
                    ic('Unable to plot the counters. No initial points.')
                    return

                if ok1 * ok2 == 0:
                    cl = ic('Counter functions or impact energy missing.')
                else:
                    # if ss > 0:
                    #     cl = ic(np.array(['Plotting the triplot (counter, enhanced ', 'counter and impact energy).']))

                    if kind == 'counter function' or kind == 'triplot':
                        # fig, axs = plt.subplots(3)
                        axs[0].plot(Efl / 1e6, C / n, lw=2, label=label)
                        axs[0].set_ylabel("$c_" + "{" + f"{N}" + "}/ c_0 $")
                        axs[0].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                        # axs[0].set_title(r'$\mathbf{MultiPac 2.1~~~~~Counter function~~~~}$')
                        axs[0].set_xlim(np.amin(Efl) / 1e6, np.amax(Efl) / 1e6)
                        axs[0].set_ylim(0, np.max([0.1, axs[0].get_ylim()[1]]))

                        # plot peak operating field
                        axs[0].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                        axs[0].text(np.round(Eacc * Epk_Eacc, 2) - 1.5, 0.1,
                                    f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                    size=12, rotation=90,
                                    transform=axs[0].get_xaxis_transform())

                        axs[0].minorticks_on()
                        axs[0].legend(loc='upper left')

                    if kind == 'final impact energy' or kind == 'triplot':
                        s = 0
                        if kind == 'final impact energy':
                            s = 1
                        axs[1 - s].semilogy(Efl / 1e6, Efq, lw=2, label=label)

                        # axs[1-s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e1, 0], secy1[e1, 0]], '-r')
                        e0 = sci.interp1d(secy1[0:e1 + 1, 1], secy1[0:e1 + 1, 0])(1)
                        axs[1 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [e0, e0], '-r')
                        axs[1 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e2, 0], secy1[e2, 0]], '-r')
                        axs[1 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e3, 0], secy1[e3, 0]], '--r')

                        axs[1 - s].set_ylabel("$Ef_" + "{" + f"{N}" + "}$")
                        axs[1 - s].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                        # axs[1-s].set_title('$\mathbf{Final~Impact~Energy~in~eV}$')
                        axs[1 - s].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
                        axs[1 - s].set_ylim(0, axs[1 - s].get_ylim()[1])

                        axs[1 - s].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                        axs[1 - s].text(np.round(Eacc * Epk_Eacc, 2) - 1.5, 0.1,
                                        f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                        size=12, rotation=90,
                                        transform=axs[1 - s].get_xaxis_transform())

                        axs[1 - s].minorticks_on()
                        axs[1 - s].legend(loc='upper left')
                    if kind == 'enhanced counter function' or kind == 'triplot':
                        s = 0
                        if kind == 'enhanced counter function':
                            s = 2
                        axs[2 - s].semilogy(Efl / 1e6, (A + 1) / n, lw=2, label=label)
                        axs[2 - s].set_xlabel('$V$ [MV]')
                        axs[2 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [1, 1], '-r')
                        axs[2 - s].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
                        axs[2 - s].set_ylim(np.min((A + 1) / n), axs[2 - s].get_ylim()[1])
                        axs[2 - s].set_ylabel("$e_" + "{" + f"{N}" + "}" + "/ c_0$")
                        axs[2 - s].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                        # axs[2-s].set_title('$\mathbf{Enhanced~counter~function}$')

                        axs[2 - s].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                        axs[2 - s].text(np.round(Eacc * Epk_Eacc, 2) - 1, 0.1,
                                        f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                        size=12, rotation=90,
                                        transform=axs[2 - s].get_xaxis_transform())

                        axs[2 - s].minorticks_on()
                        axs[2 - s].legend(loc='upper left')

        fig.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_{kind.replace(' ', '_')}.png")

        plt.show()

    def plot_dispersion(self):
        """
        Plot dispersion curve for the cavities

        Returns
        -------

        """
        fig, ax = plt.subplots()
        for cav in self.cavities_list:
            x = range(1, cav.n_cells + 1)
            ax.plot(x, cav.d_slans_all_results['FREQUENCY'][0:cav.n_cells], marker='o', mec='k',
                    label=f'{cav.name} (kcc={round(cav.k_cc, 2)} %)')
            ax.set_xlabel('Mode Number')
            ax.set_ylabel('Frequency [MHz]')

        plt.legend()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_dispersion.png")

        plt.show()

    def write_contour(self, cav, opt='mid', n_cells=1):
        """
        Write geometric contour for cavities

        Parameters
        ----------
        cav: Cavity object
            Cavity object
        opt: str

        n_cells: int
            Number of cavity cells

        Returns
        -------

        """

        if opt.lower() == 'mid':
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            n_cell = 1
            L_bp_l = 0.001
            L_bp_r = 0.001

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

        elif opt.lower() == 'end':
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['OC']) * 1e-3
            L_bp_l = 0.001
            L_bp_r = 1 * L_m

            n_cell = 1

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m) / 2
        else:
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.d_geom_params['OC']) * 1e-3
            try:
                A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['OC_R']) * 1e-3
            except KeyError:
                A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['OC']) * 1e-3

            L_bp_l = 4 * L_m
            L_bp_r = 4 * L_m

            n_cell = n_cells

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

        step = 2  # step in boundary points in mm
        # shift = 0
        # shift = L_m  # for end cell

        # calculate angles outside loop
        # CALCULATE x1_el, y1_el, x2_el, y2_el
        data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
                [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

        x1el, y1el, x2el, y2el = fsolve(ellipse_tangent, np.array(
            [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req_el - 0.85 * B_el]),
                                        args=data,
                                        xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # CALCULATE x1, y1, x2, y2
        data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
                [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(ellipse_tangent,
                                np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req_m - 0.85 * B_m]),
                                args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # CALCULATE x1_er, y1_er, x2_er, y2_er
        data = ([0 + L_bp_r, Ri_er + b_er, L_er + L_bp_r, Req_er - B_er],
                [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1er, y1er, x2er, y2er = fsolve(ellipse_tangent, np.array(
            [a_er + L_bp_r, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_r, Req_er - 0.85 * B_er]),
                                        args=data,
                                        xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        with open(fr'{cav.slans_dir}\contour.txt', 'w') as fil:
            # SHIFT POINT TO START POINT
            start_point = [-shift, 0]
            fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

            self.lineTo(start_point, [-shift, Ri_el], step)
            pt = [-shift, Ri_el]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            # ADD BEAM PIPE LENGTH
            self.lineTo(pt, [L_bp_l - shift, Ri_el], step)
            pt = [L_bp_l - shift, Ri_el]
            print(pt)
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            for n in range(1, n_cell + 1):
                if n == 1:
                    # DRAW ARC:
                    pts = self.arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                    pt = [-shift + x1el, y1el]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    self.lineTo(pt, [-shift + x2el, y2el], step)
                    pt = [-shift + x2el, y2el]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = self.arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt,
                                     [L_bp_l + L_el - shift, Req_el])
                    pt = [L_bp_l + L_el - shift, Req_el]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    if n_cell == 1:
                        # EQUATOR ARC TO NEXT POINT
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        ic(pt, 1)
                        pts = self.arcTo(L_el + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                         [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req_er])
                        pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                        ic(pt, 2)
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        ic(pt, 3)

                        # STRAIGHT LINE TO NEXT POINT
                        self.lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                        pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        ic(pt, 4, L_el + L_er - x1er + L_bp_l + L_bp_r - shift)

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        ic(shift)
                        pts = self.arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                         [L_bp_l + L_el + L_er - shift, y1er])
                        ic(pt, 5, L_el + L_er + L_bp_l - shift)

                        pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        ic(pt, 6)

                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_er)
                        ic(shift)
                    else:
                        print("if else")
                        # EQUATOR ARC TO NEXT POINT
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = self.arcTo(L_el + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                         [L_el + L_m - x2 + 2 * L_bp_l - shift, Req_m])
                        pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # STRAIGHT LINE TO NEXT POINT
                        self.lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                        pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = self.arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                         [L_bp_l + L_el + L_m - shift, y1])
                        pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_m)
                        ic(shift)

                elif n > 1 and n != n_cell:
                    print("elif")
                    # DRAW ARC:
                    pts = self.arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                    pt = [-shift + x1, y1]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    self.lineTo(pt, [-shift + x2, y2], step)
                    pt = [-shift + x2, y2]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt,
                                     [L_bp_l + L_m - shift, Req_m])
                    pt = [L_bp_l + L_m - shift, Req_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                     [L_m + L_m - x2 + 2 * L_bp_l - shift, Req_m])
                    pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    self.lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                     [L_bp_l + L_m + L_m - shift, y1])
                    pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                    ic(pt)
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - 2 * L_m
                else:
                    print("else")
                    # DRAW ARC:
                    pts = self.arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                    pt = [-shift + x1, y1]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    self.lineTo(pt, [-shift + x2, y2], step)
                    pt = [-shift + x2, y2]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt,
                                     [L_bp_l + L_m - shift, Req_m])
                    pt = [L_bp_l + L_m - shift, Req_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                     [L_m + L_er - x2er + 2 * L_bp_l - shift, Req_er])
                    pt = [L_m + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    self.lineTo(pt, [L_m + L_er - x1er + 2 * L_bp_l - shift, y1er], step)
                    pt = [L_m + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                     [L_bp_l + L_m + L_er - shift, y1er])
                    pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            # BEAM PIPE
            # reset shift
            print("pt before", pt)
            shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
            self.lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + L_el + L_er - shift, Ri_er], step)
            pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")
            print("pt after", pt)

            # END PATH
            self.lineTo(pt, [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0],
                        step)  # to add beam pipe to right
            pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
            # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
            # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

            # CLOSE PATH
            self.lineTo(pt, start_point, step)
            fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # plt.show()

    def ql_pin(self, labels, geometry, RF, QOI, Machine, p_data=None):
        """
        Calculate the value of input power as a function of loaded quality factor

        Parameters
        ----------
        labels: list, array like
            Descriptive labels on matplotlib plot
        geometry: list, array like
            List of grouped geometric input parameters
        RF: list, array like
            List of grouped radio-frequency (RF) properties
        QOI:
            List of quantities of interest for cavities
        Machine:
            List of grouped machine related materials
        p_data:


        Returns
        -------

        """
        # check if entries are of same length

        it = iter(geometry)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(RF)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(QOI)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(Machine)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        n_cells, l_cells, G, b = [np.array(x) for x in geometry]
        E_acc, Vrf = [np.array(x) for x in RF]

        fig, ax = plt.subplots()

        # QOI
        f0, R_Q = [np.array(x) for x in QOI]

        # Machine
        I0, rho, E0 = [np.array(x) for x in Machine]

        l_active = 2 * n_cells * l_cells
        l_cavity = l_active + 8 * l_cells

        # CALCULATED
        v_cav = E_acc * l_active

        U_loss = 88.46 * E0 ** 4 / rho * 1e-6  # GeV # energy lost per turn per beam
        v_loss = U_loss * 1e9  # V # v loss per beam

        print(v_loss, Vrf, v_loss / Vrf)
        phi = np.arccos(v_loss / Vrf)
        delta_f = -R_Q * f0 * I0 * np.sin(phi) / (2 * v_cav)  # optimal df
        QL_0_x = v_cav / (R_Q * I0 * np.cos(phi))  # optimal Q loaded

        QL_0 = np.linspace(1e4, 1e9, 1000000)

        xy_list = [(0.15, 0.13), (0.1, 0.16), (0.1, 0.19), (0.1, 0.21)]
        for i in range(len(E_acc)):
            f1_2 = f0[i] / (2 * QL_0)  # 380.6
            ic(R_Q[i], v_cav[i], Vrf[i])
            pin = v_cav[i] ** 2 / (4 * R_Q[i] * QL_0) * \
                  ((1 + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.cos(phi[i])) ** 2 +
                   ((delta_f[i] / f1_2) + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.sin(phi[i])) ** 2)

            # material/ wall power
            e_acc = np.linspace(0.5, 25, 1000) * 1e6  # MV/m

            txt = labels[i]

            if "*" in labels[i]:
                l = ax.plot(QL_0, pin * 1e-3, label=txt, lw=4,
                            ls='--')
            else:
                l = ax.plot(QL_0, pin * 1e-3, label=txt, lw=4)

            # add annotations
            ic(l_active)

            # annotext = ax.annotate(txt, xy=xy_list[i], xycoords='figure fraction', size=8, rotation=0,
            #                        c=l[0].get_color())

        if p_data:
            # plot QL with penetration
            ax_2 = ax.twinx()
            data = fr.excel_reader(p_data)
            data_ = data[list(data.keys())[0]]
            ax_2.plot(data_["QL"], data_["penetration"], lw=4)

        # plot decorations
        ax.set_xlabel(r"$Q_{L,0}$")
        ax.set_ylabel(r"$P_\mathrm{in} ~[\mathrm{kW}]$")
        ax.set_xscale('log')
        ax.set_xlim(5e3, 1e9)
        ax.set_ylim(0, 3000)
        ax.legend(loc='upper left')  #
        ax.minorticks_on()
        # ax.grid(which='both')
        fig.show()

    def run_slans(self):
        for cav in self.cavities_list:
            cav.run_slans()

    def run_abci(self):
        for cav in self.cavities_list:
            cav.run_abci()

    def run_multipacting(self):
        for cav in self.cavities_list:
            cav.run_multipacting()

    @staticmethod
    def linspace(start, stop, step=1.):
        """
        Like np.linspace but uses step instead of num
        This is inclusive to stop, so if start=1, stop=3, step=0.5
        Output is: array([1., 1.5, 2., 2.5, 3.])
        """
        if start < stop:
            ll = np.linspace(start, stop, int((stop - start) / abs(step) + 1))
            if stop not in ll:
                ll = np.append(ll, stop)

            return ll
        else:
            ll = np.linspace(stop, start, int((start - stop) / abs(step) + 1))
            if start not in ll:
                ll = np.append(ll, start)
            return ll

    @staticmethod
    def lineTo(prevPt, nextPt, step):
        if prevPt[0] == nextPt[0]:
            # vertical line
            # chwxk id nextPt is greater
            if prevPt[1] < nextPt[1]:
                py = np.linspace(prevPt[1], nextPt[1], step)
            else:
                py = np.linspace(nextPt[1], prevPt[1], step)
                py = py[::-1]
            px = np.ones(len(py)) * prevPt[0]

        elif prevPt[1] == nextPt[1]:
            # horizontal line
            if prevPt[0] < nextPt[1]:
                px = np.linspace(prevPt[0], nextPt[0], step)
            else:
                px = np.linspace(nextPt[0], prevPt[0], step)

            py = np.ones(len(px)) * prevPt[1]
        else:
            # calculate angle to get appropriate step size for x and y
            ang = np.arctan((nextPt[1] - prevPt[1]) / (nextPt[0] - prevPt[0]))
            if prevPt[0] < nextPt[0] and prevPt[1] < nextPt[1]:
                px = np.arange(prevPt[0], nextPt[0], step * np.cos(ang))
                py = np.arange(prevPt[1], nextPt[1], step * np.sin(ang))
            elif prevPt[0] > nextPt[0] and prevPt[1] < nextPt[1]:
                px = np.arange(nextPt[0], prevPt[0], step * np.cos(ang))
                px = px[::-1]
                py = np.arange(prevPt[1], nextPt[1], step * np.sin(ang))
            elif prevPt[0] < nextPt[0] and prevPt[1] > nextPt[1]:
                px = np.arange(prevPt[0], nextPt[0], step * np.cos(ang))
                py = np.arange(nextPt[1], prevPt[1], step * np.sin(ang))
                py = py[::-1]
            else:
                px = np.arange(nextPt[0], prevPt[0], step * np.cos(ang))
                px = px[::-1]
                py = np.arange(nextPt[1], prevPt[1], step * np.sin(ang))
                py = py[::-1]

        # plt.plot(px, py)

    @staticmethod
    def arcTo2(x_center, y_center, a, b, step, start_angle, end_angle):
        u = x_center  # x-position of the center
        v = y_center  # y-position of the center
        a = a  # radius on the x-axis
        b = b  # radius on the y-axis
        sa = (start_angle / 360) * 2 * np.pi  # convert angle to radians
        ea = (end_angle / 360) * 2 * np.pi  # convert angle to radians

        if ea < sa:
            # end point of curve
            x_end, y_end = u + a * np.cos(sa), v + b * np.sin(sa)

            t = np.arange(ea, sa, np.pi / 100)
            # t = np.linspace(ea, sa, 100)
            # check if end angle is included, include if not
            if sa not in t:
                t = np.append(t, sa)
            t = t[::-1]
        else:
            # end point of curve
            x_end, y_end = u + a * np.cos(ea), v + b * np.sin(ea)

            t = np.arange(sa, ea, np.pi / 100)
            # t = np.linspace(ea, sa, 100)
            if ea not in t:
                t = np.append(t, ea)

        # print("t0 ", [(u + a * np.cos(t))[0], (v + b * np.sin(t))[0]])
        # ic([u + a * np.cos(t), v + b * np.sin(t)])
        # ic()

        # plt.plot(u + a * np.cos(t), v + b * np.sin(t))

        return [x_end, y_end]

    @staticmethod
    def arcTo(x_center, y_center, a, b, step, start, end):
        u = x_center  # x-position of the center
        v = y_center  # y-position of the center
        a = a  # radius on the x-axis
        b = b  # radius on the y-axis

        t = np.arange(0, 2 * np.pi, np.pi / 100)

        x = u + a * np.cos(t)
        y = v + b * np.sin(t)
        pts = np.column_stack((x, y))
        inidx = np.all(np.logical_and(np.array(start) < pts, pts < np.array(end)), axis=1)
        inbox = pts[inidx]
        inbox = inbox[inbox[:, 0].argsort()]

        # plt.plot(inbox[:, 0], inbox[:, 1])

        return inbox

    def make_latex_summary_tables(self):
        try:
            l1 = r"\begin{table}[!htb]"
            l2 = r"\centering"
            l3 = r"\caption{Geometric parameters and QoIs of cavities.}"
            l4 = r"\resizebox{\textwidth}{!}{\begin{tabular}{l" + f"{''.join(['c' for i in self.cavities_list])}" + "}"
            l5 = r"\toprule"
            l6 = r" ".join([fr"& {cav.name} " for cav in self.cavities_list]) + r" \\"
            l7 = r"\midrule"
            l8 = r"\midrule"
            l9 = r"$A$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][0], 2)}/{round(cav.d_geom_params['OC'][0], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l10 = r"$B$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][1], 2)}/{round(cav.d_geom_params['OC'][1], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l11 = r"$a$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][2], 2)}/{round(cav.d_geom_params['OC'][2], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l12 = r"$b$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][3], 2)}/{round(cav.d_geom_params['OC'][3], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l13 = r"$R_\mathrm{i}$ " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][4], 2)}/{round(cav.d_geom_params['OC'][4], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l14 = r"$L$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][5], 2)}/{round(cav.d_geom_params['OC'][5], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l15 = r"$R_\mathrm{eq}$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][6], 2)}/{round(cav.d_geom_params['OC'][6], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l16 = r"$ \alpha [^\circ]$" + "".join(
                [fr"& {round(cav.d_geom_params['IC'][7], 2)}/{round(cav.d_geom_params['OC'][7], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            l17 = r"\midrule"
            l18 = r"\midrule"
            l19 = r"$R/Q [\Omega$] " + "".join([fr"& {round(cav.R_Q, 2)} " for cav in self.cavities_list]) + r" \\"
            l20 = r"$G [\Omega$] " + "".join([fr"& {round(cav.G, 2)} " for cav in self.cavities_list]) + r" \\"
            l21 = r"$G.R/Q [10^4\Omega^2]$ " + "".join(
                [fr"& {round(cav.GR_Q, 2)} " for cav in self.cavities_list]) + r" \\"
            l22 = r"$E_{\mathrm{pk}}/E_{\mathrm{acc}}$ " + "".join(
                [fr"& {round(cav.e, 2)} " for cav in self.cavities_list]) + r" \\"
            l23 = r"$B_{\mathrm{pk}}/E_{\mathrm{acc}} [\mathrm{\frac{mT}{MV/m}}]$ " + "".join(
                [fr"& {round(cav.b, 2)} " for cav in self.cavities_list]) + r" \\"
            l24 = r"$|k_\mathrm{FM}| \mathrm{[SR/BS]} [\mathrm{V/pC}]$ " + "".join(
                [fr"& {round(cav.k_fm, 4)} " for cav in self.cavities_list]) + r" \\"
            l25 = r"$|k_\mathrm{\parallel}| \mathrm{[SR/BS]} [\mathrm{V/pC}]$ " + "".join(
                [fr"& {round(cav.k_loss, 4)} " for cav in self.cavities_list]) + r" \\"
            l26 = r"$k_\mathrm{\perp} \mathrm{[SR/BS]} [\mathrm{V/pC/m}]$ " + "".join(
                [fr"& {round(cav.k_kick, 4)} " for cav in self.cavities_list]) + r" \\"
            l27 = r"\midrule"
            l28 = r"\midrule"

            l29 = r"$N_\mathrm{cav/beam}$ " + "".join(
                [fr"& {round(qoi[r'N_cav/beam'], 2)} " for qoi in self.p_qois]) + r" \\"
            l29a = r"$P_\mathrm{in}\mathrm{/cav} [\mathrm{kW}]$ " + "".join(
                [fr"& {round(qoi[r'P_in/beam [kW]'], 2)} " for qoi in self.p_qois]) + r" \\"
            l30 = r"$P_\mathrm{stat}\mathrm{/cav} [\mathrm{W}]$ " + "".join(
                [fr"& {round(qoi[r'P_stat/cav [W]'], 2)} " for qoi in self.p_qois]) + r" \\"
            l31 = r"$P_\mathrm{dyn}\mathrm{/cav} [\mathrm{W}]$ " + "".join(
                [fr"& {round(qoi[r'P_dyn/cav [W]'], 2)} " for qoi in self.p_qois]) + r" \\"
            l32 = r"$P_\mathrm{HOM}\mathrm{/cav} \mathrm{[SR/BS]} [\mathrm{kW}]$ " + "".join(
                [fr"& {round(cav.phom, 2)} " for cav in self.cavities_list]) + r" \\"
            l33 = r"\bottomrule"
            l34 = r"\end{tabular}}"
            l35 = r"\label{tab: selected shape}"
            l36 = r"\end{table}"

            all_lines = (l1, l2, l3, l4, l5, l6, l7, l8, l9, l10,
                         l11, l12, l13, l14, l15, l16, l17, l18, l19, l20,
                         l21, l22, l23, l24, l25, l26, l27, l28, l29, l29a, l30,
                         l31, l32, l33, l34, l35, l36)

            with open(
                    fr"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\Summaries\{self.save_folder}_latex_summary.txt",
                    'w') as f:
                for ll in all_lines:
                    f.write(ll + '\n')
        except KeyError:
            print("Either SLANS or ABCI results not available. Please use '<cav>.set_slans_qois(<folder>)' "
                  "or '<cav>.set_abci_qois(<folder>)' to fix this.")

    def make_excel_summary(self):
        try:
            data = {'Name': [cav.name for cav in self.cavities_list],
                    'Project': [cav.project for cav in self.cavities_list],
                    'Type': [cav.type for cav in self.cavities_list],
                    'CW/Pulsed': [cav.cw_pulsed for cav in self.cavities_list],
                    'Material': [cav.material for cav in self.cavities_list],
                    'N_cells': [cav.n_cells for cav in self.cavities_list],
                    'Freq [MHz]': [cav.op_freq for cav in self.cavities_list],
                    'Beta': [cav.beta for cav in self.cavities_list],
                    'T_oper [K]': [cav.op_temp for cav in self.cavities_list],
                    'I0 [mA]': [cav.I0 for cav in self.cavities_list],
                    'sigma [mm]': [cav.sigma for cav in self.cavities_list],
                    'A_i [mm]': [round(cav.d_geom_params['IC'][0], 2) for cav in self.cavities_list],
                    'B_i [mm]': [round(cav.d_geom_params['IC'][1], 2) for cav in self.cavities_list],
                    'a_i [mm]': [round(cav.d_geom_params['IC'][2], 2) for cav in self.cavities_list],
                    'b_i [mm]': [round(cav.d_geom_params['IC'][3], 2) for cav in self.cavities_list],
                    'R_i [mm]': [round(cav.d_geom_params['IC'][4], 2) for cav in self.cavities_list],
                    'L_i [mm]': [round(cav.d_geom_params['IC'][5], 2) for cav in self.cavities_list],
                    'Req [mm]': [round(cav.d_geom_params['IC'][6], 2) for cav in self.cavities_list],
                    'alpha_i [deg]': [round(cav.d_geom_params['IC'][7], 2) for cav in self.cavities_list],
                    'A_el [mm]': [round(cav.d_geom_params['OC'][0], 2) for cav in self.cavities_list],
                    'B_el [mm]': [round(cav.d_geom_params['OC'][1], 2) for cav in self.cavities_list],
                    'a_el [mm]': [round(cav.d_geom_params['OC'][2], 2) for cav in self.cavities_list],
                    'b_el [mm]': [round(cav.d_geom_params['OC'][3], 2) for cav in self.cavities_list],
                    'R_el [mm]': [round(cav.d_geom_params['OC'][4], 2) for cav in self.cavities_list],
                    'L_el [mm]': [round(cav.d_geom_params['OC'][5], 2) for cav in self.cavities_list],
                    # 'Req [mm]': [round(cav.d_geom_params['OC'][6], 2) for cav in self.cavities_list],
                    'alpha__el [deg]': [round(cav.d_geom_params['OC'][7], 2) for cav in self.cavities_list],
                    'A_er [mm]': [round(cav.d_geom_params['OC'][0], 2) for cav in self.cavities_list],
                    'B_er [mm]': [round(cav.d_geom_params['OC'][1], 2) for cav in self.cavities_list],
                    'a_er [mm]': [round(cav.d_geom_params['OC'][2], 2) for cav in self.cavities_list],
                    'b_er [mm]': [round(cav.d_geom_params['OC'][3], 2) for cav in self.cavities_list],
                    'R_er [mm]': [round(cav.d_geom_params['OC'][4], 2) for cav in self.cavities_list],
                    'L_er [mm]': [round(cav.d_geom_params['OC'][5], 2) for cav in self.cavities_list],
                    # 'Req [mm]': [round(cav.d_geom_params['OC'][6], 2) for cav in self.cavities_list],
                    'alpha_er [deg]': [round(cav.d_geom_params['OC'][7], 2) for cav in self.cavities_list],
                    'R_shunt [Ohm]': ['' for cav in self.cavities_list],
                    'R/Q [Ohm]': [cav.R_Q for cav in self.cavities_list],
                    'k_cc [%]': [cav.k_cc for cav in self.cavities_list],
                    'field flatness [%]': [cav.ff for cav in self.cavities_list],
                    'L_active [m]': [cav.l_active for cav in self.cavities_list],
                    'Epk/Eacc []': [cav.e for cav in self.cavities_list],
                    'Bpk/Eacc [mT/MV/m]': [cav.b for cav in self.cavities_list],
                    'G [Ohm]': [cav.G for cav in self.cavities_list],
                    'R/Q.G [Ohm^2]': [cav.GR_Q for cav in self.cavities_list],
                    '|k_loss| [V/pC]': [cav.k_loss for cav in self.cavities_list],
                    '|k_kick| [V/pC/m]': [cav.k_kick for cav in self.cavities_list],
                    'P_HOM/cav [kW]': [cav.phom for cav in self.cavities_list],
                    'Reference': [cav.reference for cav in self.cavities_list]
                    }

            df = pd.DataFrame.from_dict(data)
            df.to_excel(
                fr"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\Summaries\{self.save_folder}_excel_summary.xlsx",
                sheet_name='Cavities')
        except Exception as e:
            print("Either SLANS or ABCI results not available. Please use '<cav>.set_slans_qois(<folder>)' "
                  "or '<cav>.set_abci_qois(<folder>)' to fix this.")
            print(e)

    def remove_cavity(self, cav):
        """
        Removes cavity from cavity list
        Parameters
        ----------
        cav: object
            Cavity object

        Returns
        -------

        """
        self.cavities_list.remove(cav)

    def save_all_plots(self, plot_name):
        """
        Save all plots
        Parameters
        ----------
        plot_name: str
            Name of saved plot

        Returns
        -------

        """
        if self.save_folder != '':
            # check if folder exists
            if os.path.exists(fr"D:\Dropbox\Quick presentation files\{self.save_folder}"):
                save_folder = fr"D:\Dropbox\Quick presentation files\{self.save_folder}"
                plt.savefig(f"{save_folder}/{plot_name}")
            else:
                save_folder = fr"D:\Dropbox\Quick presentation files\{self.save_folder}"
                os.mkdir(save_folder)
                plt.savefig(f"{save_folder}/{plot_name}")

    def __str__(self):
        return fr"{self.cavities_list}"


class OperationPoints:
    def __init__(self, filepath):
        self.op_points = {}

        if os.path.exists(filepath):
            self.op_points = self.load_operation_point(filepath)

    def load_operation_point(self, filepath):
        with open(filepath, 'r') as f:
            op_points = json.load(f)

        self.op_points = op_points
        return op_points

    def get_default_operation_points(self):
        self.op_points = {
            "Z": {
                "freq [MHz]": 400.79,
                "E [GeV]": 45.6,
                "I0 [mA]": 1400,
                "V [GV]": 0.12,
                "Eacc [MV/m]": 5.72,
                "nu_s []": 0.0025,
                "alpha_p [1e-5]": 1.48,
                "tau_z [ms]": 424.6,
                "tau_xy [ms]": 849.2,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 56,
                "T [K]": 4.5,
                "sigma_SR [mm]": 4.32,
                "sigma_BS [mm]": 15.2
            },
            "W": {
                "freq [MHz]": 400.79,
                "E [GeV]": 80,
                "I0 [mA]": 135,
                "V [GV]": 1.0,
                "Eacc [MV/m]": 11.91,
                "nu_s []": 0.0506,
                "alpha_p [1e-5]": 1.48,
                "tau_z [ms]": 78.7,
                "tau_xy [ms]": 157.4,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 112,
                "T [K]": 4.5,
                "sigma_SR [mm]": 3.55,
                "sigma_BS [mm]": 7.02
            },
            "H": {
                "freq [MHz]": 400.79,
                "E [GeV]": 120,
                "I0 [mA]": 53.4,
                "V [GV]": 2.1,
                "Eacc [MV/m]": 11.87,
                "nu_s []": 0.036,
                "alpha_p [1e-5]": 0.73,
                "tau_z [ms]": 23.4,
                "tau_xy [ms]": 46.8,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 118,
                "T [K]": 4.5,
                "sigma_SR [mm]": 2.5,
                "sigma_BS [mm]": 4.45
            },
            "ttbar": {
                "freq [MHz]": 801.58,
                "E [GeV]": 182.5,
                "I0 [mA]": 10,
                "V [GV]": 8.8,
                "Eacc [MV/m]": 24.72,
                "nu_s []": 0.087,
                "alpha_p [1e-5]": 0.73,
                "tau_z [ms]": 6.8,
                "tau_xy [ms]": 13.6,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 140,
                "T [K]": 2,
                "sigma_SR [mm]": 1.67,
                "sigma_BS [mm]": 2.54
            },
            "MuCol RCS Stage 1": {
                "freq [MHz]": 1300,
                "E [GeV]": 250.835,
                "I0 [mA]": 20.38,
                "V [GV]": 20.87,
                "Eacc [MV/m]": 30,
                "nu_s []": 0.68,
                "alpha_p [1e-5]": 240,
                "tau_z [ms]": 78.7,
                "tau_xy [ms]": 157.4,
                "f_rev [kHz]": 50.08,
                "beta_xy [m]": 0.005,
                "N_c []": 696,
                "T [K]": 2,
                "sigma_SR [mm]": 23.1,
                "sigma_BS [mm]": 23.1
            },
            "MuCol hybrid RCS Stage 2": {
                "freq [MHz]": 1300,
                "E [GeV]": 436.15,
                "I0 [mA]": 18.78,
                "V [GV]": 11.22,
                "Eacc [MV/m]": 30,
                "nu_s []": 0.32,
                "alpha_p [1e-5]": 240,
                "tau_z [ms]": 0,
                "tau_xy [ms]": 0,
                "f_rev [kHz]": 50.08,
                "beta_xy [m]": 0.005,
                "N_c []": 374,
                "T [K]": 2,
                "sigma_SR [mm]": 23.1,
                "sigma_BS [mm]": 23.1
            },
            "MuCol hybrid RCS Stage 3": {
                "freq [MHz]": 1300,
                "E [GeV]": 750.024,
                "I0 [mA]": 9.97,
                "V [GV]": 16.07,
                "Eacc [MV/m]": 30,
                "nu_s []": 0.24,
                "alpha_p [1e-5]": 240,
                "tau_z [ms]": 0,
                "tau_xy [ms]": 0,
                "f_rev [kHz]": 28.04,
                "beta_xy [m]": 0.005,
                "N_c []": 536,
                "T [K]": 2,
                "sigma_SR [mm]": 23.1,
                "sigma_BS [mm]": 23.1
            },
            "MuCol hybrid RCS Stage 5TeV": {
                "freq [MHz]": 1300,
                "E [GeV]": 3500,
                "I0 [mA]": 3.05,
                "V [GV]": 90,
                "Eacc [MV/m]": 30,
                "nu_s []": 0.57,
                "alpha_p [1e-5]": 240,
                "tau_z [ms]": 0,
                "tau_xy [ms]": 0,
                "f_rev [kHz]": 8.57,
                "beta_xy [m]": 0.005,
                "N_c []": 3000,
                "T [K]": 2,
                "sigma_SR [mm]": 23.1,
                "sigma_BS [mm]": 23.1
            },
            "MuCol hybrid RCS Stage 5TeV LHC": {
                "freq [MHz]": 1300,
                "E [GeV]": 3500,
                "I0 [mA]": 4,
                "V [GV]": 68.75,
                "Eacc [MV/m]": 30,
                "nu_s []": 0.43,
                "alpha_p [1e-5]": 240,
                "tau_z [ms]": 0,
                "tau_xy [ms]": 0,
                "f_rev [kHz]": 11.25,
                "beta_xy [m]": 0.005,
                "N_c []": 2292,
                "T [K]": 2,
                "sigma_SR [mm]": 23.1,
                "sigma_BS [mm]": 23.1
            }
        }


def get_qois_value(f_fm, R_Q, k_loss, k_kick, sigma_z, I0, Nb, n_cell):
    c = 299792458
    w_fm = 2 * np.pi * f_fm * 1e6
    e = 1.602e-19

    k_fm = (w_fm / 4) * R_Q * np.exp(-(w_fm * sigma_z * 1e-3 / c) ** 2) * 1e-12
    k_hom = k_loss - k_fm
    p_hom = (k_hom * 1e12) * (I0 * 1e-3) * e * (Nb * 1e11)

    d = {
        "n cell": n_cell,
        # "freq [MHz]": f_fm,
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


def show_valid_operating_point_structure():
    dd = {
        '<wp1>': {
            'I0 [mA]': '<value>',
            'Nb [1e11]': '<value>',
            'sigma_z (SR/BS) [mm]': '<value>'
        },
        '<wp2>': {
            'I0 [mA]': '<value>',
            'Nb [1e11]': '<value>',
            'sigma_z (SR/BS) [mm]': '<value>'
        }
    }

    info(dd)


def get_surface_resistance(Eacc, b, m, freq, T):
    Rs_dict = {
        "Rs_NbCu_2K_400.79Mhz": 0.57 * (Eacc * 1e-6 * b) + 28.4,  # nOhm
        "Rs_NbCu_4.5K_400.79Mhz": 39.5 * np.exp(0.014 * (Eacc * 1e-6 * b)) + 27,  # nOhm
        "Rs_bulkNb_2K_400.79Mhz": (2.33 / 1000) * (Eacc * 1e-6 * b) ** 2 + 26.24,  # nOhm
        "Rs_bulkNb_4.5K_400.79Mhz": 0.0123 * (Eacc * 1e-6 * b) ** 2 + 62.53,  # nOhm

        "Rs_NbCu_2K_801.58Mhz": 1.45 * (Eacc * 1e-6 * b) + 92,  # nOhm
        "Rs_NbCu_4.5K_801.58Mhz": 50 * np.exp(0.033 * (Eacc * 1e-6 * b)) + 154,  # nOhm
        "Rs_bulkNb_2K_801.58Mhz": (16.4 + Eacc * 1e-6 * b * 0.092) * (800 / 704) ** 2,  # nOhm
        "Rs_bulkNb_4.5K_801.58Mhz": 4 * (62.7 + (Eacc * 1e-6 * b) ** 2 * 0.012)  # nOhm
    }
    if freq < 600:
        freq = 400.79

    if freq >= 600:
        freq = 801.58

    rs = Rs_dict[fr"Rs_{m}_{T}K_{freq}Mhz"]

    return rs

# if __name__ == '__main__':
#     p = str(Path(SOFTWARE_DIRECTORY).parents[0])
#     print(p)
