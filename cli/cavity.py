import os
import shutil
import time
from pathlib import Path

from termcolor import colored

from analysis_modules.data_module.abci_data import ABCIData
import numpy as np
import multiprocessing as mp
from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from analysis_modules.wakefield.ABCI.abci_geometry import ABCIGeometry
from utils.shared_functions import *
import psutil

slans_geom = SLANSGeometry()
abci_geom = ABCIGeometry()
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
    def __init__(self, n_cells, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none', name='cavity'):
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
        self.abci_qois = {}

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
        self.name = name

    def set_mid_cell(self, cell):
        self.mid_cell = cell

    def set_end_cell_left(self, cell):
        self.mid_cell = cell

    def set_end_cell_right(self, cell):
        self.end_cell_right = cell

    def set_boundary_conditions(self, bc):
        self.bc = bc

    def set_beampipe(self, bp):
        self.beampipe = bp

    def load(self):
        pass

    def save(self, files_path):

        if files_path is None:
            error('Please specify a folder to write the simulation results to.')
            return
        else:
            try:
                self.folder = files_path
                success = self.create_project()
                if not success:
                    error("Project could not be created. Please check the folder and try again.")
                    return
                else:
                    done("Project created successfully. You can proceed with analysis.")
            except Exception as e:
                error("Exception occurred: ", e)
                return

        if files_path is None:
            self.folder = os.getcwd()

    def run_eigenmode(self, solver='SLANS', freq_shift=0, boundary_cond=None, subdir='',
                      UQ=False):

        if boundary_cond:
            self.bc = boundary_cond

        if solver == 'SLANS':
            self.run_slans(self.name, self.n_cells, self.n_modules, self.shape_space, self.n_modes, freq_shift,
                           self.bc, SOFTWARE_DIRECTORY, self.folder, sub_dir='', UQ=UQ)

            # load quantities of interest
            self.get_slans_qois()
        else:
            self.run_custom_eig()

    def run_wakefield(self, MROT=2, MT=10, NFS=10000, wakelength=100, bunch_length=25,
                      DDR_SIG=0.1, DDZ_SIG=0.1, WG_M=None, marker='', qoi_df=None, solver='ABCI'):

        if solver == 'ABCI':
            info(">> Running wakefield simulation")
            self.run_abci(self.name, self.n_cells, self.n_modules, self.shape_space,
                          MROT=MROT, MT=MT, NFS=NFS, UBT=wakelength, bunch_length=bunch_length,
                          DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG,
                          parentDir=SOFTWARE_DIRECTORY, projectDir=self.folder, WG_M=WG_M, marker=marker,
                          qoi_df=qoi_df)

    def calc_op_freq(self):
        if not self.freq:
            self.freq = (c0 / 4 * self.L)

    @staticmethod
    def run_custom_eig():
        pass

    @staticmethod
    def run_slans(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir, sub_dir='', UQ=False):
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
    def run_abci(name, n_cells, n_modules, shape, MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
                 DDR_SIG=0.1, DDZ_SIG=0.1,
                 parentDir=None, projectDir=None,
                 WG_M=None, marker='', qoi_df=None):

        # run abci code
        if WG_M is None:
            WG_M = ['']

        start_time = time.time()
        # run both polarizations if MROT == 2
        for ii in WG_M:
            try:
                if MROT == 2:
                    for m in range(2):
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
                    for m in range(2):
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

        if qoi_df is not None:
            d = {}
            # # save qois
            for index, row in qoi_df.iterrows():
                WP = row['WP']
                I0 = float(row['I0 [mA]'])
                Nb = float(row['Nb [1e11]'])
                sigma_z = [float(x) for x in row['sigma_z (SR/BS) [mm]'].split('/')]
                freq = float(row['f [MHz]'])

                # try getting  R_Q from corresponding SLANS folder
                try:
                    with open(fr'{projectDir}\SimulationData\SLANS\{name}\qois.json') as json_file:
                        qois_slans = json.load(json_file)
                    R_Q = qois_slans['R/Q [Ohm]']
                    print("Found a corresponding SLANS file")
                except Exception as e:
                    R_Q = float(row['R/Q [Ohm]'])
                    print("Did not find a corresponding SLANS file. Please check the R/Q used to calculate the"
                          "k_FM and P_HOM", e)

                n_cell = int(row['n cell'])

                bl_diff = ['SR', 'BS']

                for i, s in enumerate(sigma_z):
                    for ii in WG_M:
                        fid = f"{WP}_{bl_diff[i]}_{s}mm{ii}"
                        try:
                            for m in range(2):
                                abci_geom.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                                 fid=fid, MROT=m, MT=MT, NFS=NFS, UBT=10 * s * 1e-3, bunch_length=s,
                                                 DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                                 projectDir=projectDir,
                                                 WG_M=ii, marker=ii, sub_dir=f"{name}")
                        except KeyError:
                            for m in range(2):
                                abci_geom.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                                 fid=fid, MROT=m, MT=MT, NFS=NFS, UBT=10 * s * 1e-3, bunch_length=s,
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

                        d[fid] = get_qois_value(freq, R_Q, k_loss, k_kick, s, I0, Nb, n_cell)

            # save qoi dictionary
            run_save_directory = fr'{projectDir}\SimulationData\ABCI\{name}{marker}'
            with open(fr'{run_save_directory}\qois.json', "w") as f:
                json.dump(d, f, indent=4, separators=(',', ': '))

        done("Done with the secondary analysis for working points")

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

    def get_slans_qois(self):
        qois = 'qois.json'

        with open(fr"{self.folder}\SimulationData\SLANS\{self.name}\{qois}") as json_file:
            self.slans_qois = json.load(json_file)

    def get_abci_qois(self):
        qois = 'qois.json'

        with open(fr"{self.folder}\SimulationData\ABCI\{self.name}\{qois}") as json_file:
            self.abci_qois = json.load(json_file)

    @staticmethod
    def createFolder(name, projectDir, subdir=''):
        # change save directory
        path = fr'{projectDir}\{name}'
        if subdir == '':
            pass
        else:
            new_path = fr'{projectDir}\{subdir}\{name}'
            if os.path.exists(new_path):
                path = new_path
            else:
                if not os.path.exists(fr'{projectDir}\{subdir}'):
                    os.mkdir(fr'{projectDir}\{subdir}')

                os.mkdir(new_path)
                path = fr'{projectDir}\{subdir}\{name}'

        if os.path.exists(path):
            shutil.rmtree(path)
            os.mkdir(path)
        else:
            os.mkdir(path)

    def create_project(self):
        project_name = self.name

        if project_name != '':
            project_dir = self.folder

            # check if folder already exist
            e = self.checkIfPathExist(project_dir, project_name)

            if e:
                def make_dirs_from_dict(d, current_dir=project_dir):
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
                    print("An exception occurred in created project: ", e)
                    return False
            else:
                return False
        else:
            print('\tPlease enter a valid project name')
            return False

    def checkIfPathExist(self, directory, folder):
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

# if __name__ == '__main__':
#     p = str(Path(SOFTWARE_DIRECTORY).parents[0])
#     print(p)
