import fnmatch
import os.path
import random
import shutil
import subprocess
import sys
from distutils import dir_util
from math import floor

import numpy as np
from scipy.interpolate import interp2d, LinearNDInterpolator, CloughTocher2DInterpolator
import time
import matplotlib.tri as tri
import psutil
import scipy.signal as sps
from pathlib import Path

from matplotlib import pyplot as plt, animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from IPython.core.display import Image, HTML
from IPython.core.display_functions import display
from scipy.spatial import ConvexHull, Delaunay
from scipy.special import jn_zeros, jnp_zeros
import matplotlib as mpl
import scipy.io as spio
import scipy.interpolate as sci
import mplcursors
import pandas as pd
from scipy.stats import qmc
from tqdm import trange
from tqdm.auto import tqdm
import time
import datetime
import oapackage
import multiprocessing as mp
from analysis_modules.plot_module.matplotlib_annotation_objects import DraggableText
from analysis_modules.tune.tuners.tuner import Tuner
from analysis_modules.data_module.abci_data import ABCIData
from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from analysis_modules.eigenmode.NGSolve.eigen_ngsolve import NGSolveMEVP
from analysis_modules.eigenmode.customEig.run_field_solver import Model
from analysis_modules.wakefield.ABCI.abci_geometry import ABCIGeometry
from utils.shared_functions import *

slans_geom = SLANSGeometry()
ngsolve_mevp = NGSolveMEVP()
abci_geom = ABCIGeometry()
custom_eig = Model()
tuner = Tuner()

SOFTWARE_DIRECTORY = Path(os.getcwd())  # str(Path().parents[0])
VAR_TO_INDEX_DICT = {'A': 0, 'B': 1, 'a': 2, 'b': 3, 'Ri': 4, 'L': 5, 'Req': 6}
TUNE_ACCURACY = 1e-4
DIMENSION = 'm'
DIMENSION_FACTOR = {'mm': 1, 'cm': 1e-1, 'm': 1e-3}

m0 = 9.1093879e-31
q0 = 1.6021773e-19
c0 = 2.99792458e8
mu0 = 4 * np.pi * 1e-7
eps0 = 8.85418782e-12


class Optimisation:

    def __init__(self):
        pass  # Cavities never calls to super(), fix that

    def start_optimisation(self, config):
        self.err = []
        self.pareto_history = []
        self.optimisation_config = config
        # apply optimisation settings
        self.parentDir = 'D:\Dropbox\CavityDesignHub'
        self.projectDir = config['project dir.']
        self.initial_points = config['initial points']
        self.tune_freq = config['tune freq.']
        self.ng_max = config['no. of generation']
        self.objectives_unprocessed = config['objectives']
        self.objectives, weights = self.process_objectives(config['objectives'])
        self.objective_vars = [obj[1] for obj in self.objectives]
        if 'weights' in config.keys():
            self.weights = config['weights']
            assert len(self.weights) == len(weights), \
                ("Length of delta must be equal to the length of the variables. For impedance Z entries, one less than"
                 "the length of the interval list weights are needed. Eg. for ['min', 'ZL', [1, 2, 3]], two weights are"
                 "required. ")
        else:
            self.weights = weights

        self.bounds = config['bounds']
        self.cell_type = config['cell type']
        self.constraints = self.process_constraints(config['constraints'])
        self.processes_count = config['processes']
        self.norm_length = config['normalisation length']
        self.method = config['method']
        self.mutation_factor = config['mutation factor']
        self.crossover_factor = config['crossover factor']
        self.elites_to_crossover = config['elites for crossover']
        self.chaos_factor = config['chaos factor']
        self.tune_variable = config['tune variable']
        self.uq_config = config['uq']
        if self.uq_config['delta']:
            assert len(self.uq_config['delta']) == len(self.uq_config['variables']), error("The number of deltas must "
                                                                                           "be equal to the number of "
                                                                                           "variables.")

        self.df = None

        # interpolation
        self.df_global = pd.DataFrame()
        self.objs_dict = {}
        self.constraints_dict = {}
        self.n_interp = 10000
        self.interp_error = []
        self.interp_error_avg = []
        bar = tqdm(total=self.ng_max)
        self.ea(0, bar)

    def ea(self, n, bar):
        if n == 0:
            # update lists
            self.df = self.generate_first_men(self.initial_points, 0)

            self.f2_interp = [np.zeros(self.n_interp) for _ in range(len(self.objectives))]

            folders = [Path(fr"{self.projectDir}\SimulationData\Optimisation"),
                       Path(fr"{self.projectDir}\SimulationData\Optimisation")]

            # clear folder to avoid reading from previous optimization attempt
            for folder in folders:
                if os.path.exists(folder):
                    for filename in os.listdir(folder):
                        try:
                            shutil.rmtree(folder / fr"{filename}")
                        except NotADirectoryError:
                            os.remove(folder / fr"{filename}")
                else:
                    os.mkdir(folder)

        # optimize by page rank
        # remove the lowest ranking members
        df = self.df

        # compare with global dict and remove duplicates
        compared_cols = ['A', 'B', 'a', 'b', 'Ri']
        if not self.df_global.empty:
            df = df.loc[~df.set_index(compared_cols).index.isin(
                self.df_global.set_index(compared_cols).index)]  # this line of code removes duplicates

        pseudo_shape_space = {}
        for index, row in df.iterrows():
            rw = row.tolist()

            if self.cell_type.lower() == 'mid cell' or self.cell_type.lower() == 'mid-cell' or self.cell_type.lower() == 'mid_cell':
                pseudo_shape_space[rw[0]] = {'IC': rw[1:], 'OC': rw[1:], 'OC_R': rw[1:], 'BP': 'none',
                                             'FREQ': self.tune_freq}

            elif self.cell_type.lower() == 'mid-end cell' or self.cell_type.lower() == 'mid-end-cell' or self.cell_type.lower() == 'mid_end_cell':

                assert 'mid cell' in list(self.optimisation_config.keys()), \
                    ("If cell type is set as 'mid-end cell', the mid cell geometry parameters must "
                     "be provided in the optimisation_config dictionary.")
                assert len(self.optimisation_config['mid cell']) > 6, ("Incomplete mid cell geometry parameter. "
                                                                       "At least 7 geometric parameters "
                                                                       "[A, B, a, b, Ri, L, Req] required.")
                IC = self.optimisation_config['mid cell']
                # check if mid-cell is not a degenerate geometry
                df = tangent_coords(*np.array(IC)[0:8], 0)
                assert df[-2] == 1, ("The mid-cell geometry dimensions given result in a degenerate geometry. "
                                     "Please check.")

                pseudo_shape_space[rw[0]] = {'IC': IC, 'OC': rw[1:], 'OC_R': rw[1:], 'BP': 'right',
                                             'FREQ': self.tune_freq}

            elif (self.cell_type.lower() == 'end-end cell' or self.cell_type.lower() == 'end-end-cell'
                  or self.cell_type.lower() == 'end_end_cell') or self.cell_type.lower() == 'end end cell':

                pseudo_shape_space[rw[0]] = {'IC': rw[1:], 'OC': rw[1:], 'OC_R': rw[1:], 'BP': 'right',
                                             'FREQ': self.tune_freq}

            else:
                pseudo_shape_space[rw[0]] = {'IC': rw[1:], 'OC': rw[1:], 'OC_R': rw[1:], 'BP': 'both',
                                             'FREQ': self.tune_freq}

        pseudo_shape_space = self.remove_duplicate_values(pseudo_shape_space)

        ############################
        # run tune
        n_cells = 1
        norm_length = self.norm_length

        self.run_tune_parallel(pseudo_shape_space, n_cells)
        # get successfully tuned geometries and filter initial generation dictionary
        processed_keys = []
        tune_result = []
        for key in pseudo_shape_space.keys():
            filename = self.projectDir / fr'SimulationData\Optimisation\{key}\tune_res.json'
            try:
                with open(filename, 'r') as file:
                    tune_res = json.load(file)

                # get only tune_variable, alpha_i, and alpha_o and freq but why these quantities
                freq = tune_res['FREQ']
                tune_variable_value = tune_res['IC'][VAR_TO_INDEX_DICT[self.tune_variable]]
                alpha_i = tune_res['IC'][7]
                alpha_o = tune_res['IC'][7]

                tune_result.append([tune_variable_value, alpha_i, alpha_o, freq])
                processed_keys.append(key)
            except FileNotFoundError:
                pass

        # after removing duplicates, dataframe might change size
        df = df.loc[df['key'].isin(processed_keys)]
        df.loc[:, [self.tune_variable, 'alpha_i', 'alpha_o', 'freq [MHz]']] = tune_result

        # eigen objective variables
        # for o in self.objectives:
        intersection = set(self.objective_vars).intersection(
            ["freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]", "G [Ohm]", "Q []"])
        if len(intersection) > 0:
            # process tune results
            obj_result = []
            processed_keys = []
            for key in pseudo_shape_space.keys():
                filename = self.projectDir / fr'SimulationData\Optimisation\{key}\monopole\qois.json'
                try:
                    with open(filename, 'r') as file:
                        qois = json.load(file)
                    # extract objectives from tune_res
                    obj = list(
                        {key: val for [key, val] in qois.items() if key in self.objective_vars}.values())

                    obj_result.append(obj)
                    # tune_result.append(list(qois.values()))
                    processed_keys.append(key)
                except FileNotFoundError as e:
                    pass

            # after removing duplicates, dataframe might change size
            if len(processed_keys) == 0:
                error("Unfortunately, none survived. \n"
                      "This is most likely due to all generated initial geometries being degenerate.\n"
                      "Check the variable bounds or increase the number of initial geometries to increase the"
                      "changes of survival. \n"
                      "Can't even say that this was a good run."
                      "Tune ended.")
                return

            df = df.loc[df['key'].isin(processed_keys)]

            obj_eigen = [o[1] for o in self.objectives if
                         o[1] in ["freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]", "G [Ohm]", "Q []"]]
            df[obj_eigen] = obj_result

        for o in self.objectives:
            if o[1] in ["mts monopole", 'mts dipole']:
                # process tune results
                obj_vars = self.ui.ccb_Populate_Objectives.currentText().split(', ')
                for i, obj_var in enumerate(obj_vars):
                    if obj_var == "mts monopole" or obj_var == "mts dipole":
                        goal = self.ui.tw_Objectives.cellWidget(i, 1).currentText()
                        if goal == 'equal':
                            fshift = float(self.ui.tw_Objectives.item(i, 2).text())

                obj_result = []
                tune_result = []
                processed_keys = []
                # run dipole simulation with frequency shift
                if o[1] == "mts monopole":
                    slans_shape_space = self.run_slans_parallel(df, n_cells, fshift, 'monopole')
                    for key, val in slans_shape_space.items():
                        filename = self.projectDir / fr'SimulationData\SLANS_opt\{key}\cavity_33.svl'
                        try:
                            params = fr.svl_reader(filename)
                            obj = self.get_objectives_value(params, self.objectives, norm_length, n_cells)

                            obj_result.append(obj)

                            df_slans_mts = pd.DataFrame(obj, columns=[key, o[1]])
                            df = df.merge(df_slans_mts, on='key', how='inner')
                        except FileNotFoundError:
                            pass
                else:
                    slans_shape_space = self.run_slans_parallel(df, n_cells, fshift, 'dipole')
                    for key, val in slans_shape_space.items():
                        filename = self.projectDir / fr'SimulationData\SLANS_opt\{key}_n1\cavity_33_2.sv2'
                        try:
                            params = fr.sv2_reader(filename)
                            obj = params['Frequency'][-1]
                            obj_result.append([key, obj])
                            processed_keys.append(key)
                        except FileNotFoundError:
                            pass

                    df_slans_mts = pd.DataFrame(obj_result, columns=['key', o[1]])
                    df = df.merge(df_slans_mts, on='key', how='inner')

        # wakefield objective variables
        for o in self.objectives:
            if "ZL" in o[1] or "ZT" in o[1] or "k_loss" in o[1] or "k_kick" in o[1]:
                # run wakefield analysis and return shape space
                wake_shape_space = self.run_wakefield_parallel(df)

                # process wakefield results
                df_wake, processed_keys = self.get_wakefield_objectives_value(wake_shape_space,
                                                                              self.projectDir / fr'SimulationData\ABCI')

                df = df.merge(df_wake, on='key', how='inner')
                break

        # apply UQ
        if self.uq_config['option']:
            solver_dict = {'ngsolvemevp': ngsolve_mevp, 'abci': abci_geom}
            # beampipes = {'Mid Cell': 'none', 'End-End Cell': 'left', 'End-Mid Cell': 'left', 'Single Cell': 'both'}
            beampipe = 'none'

            solver_args_dict = {'ngsolvemevp': {'n_cells': n_cells,
                                                'n_modules': 1,
                                                'f_shift': 0,
                                                'bc': 33,
                                                'beampipes': beampipe,
                                                # 'norm_length': self.norm_length
                                                },
                                'abci': {'n_cells': n_cells, 'n_modules': 1,
                                         'MROT': 2, 'MT': 4, 'NFS': 10000, 'UBT': 50, 'bunch_length': 25,
                                         'DDR_SIG': 0.1, 'DDZ_SIG': 0.1,
                                         'progress_list': None,
                                         'WG_M': None, 'marker': ''},
                                'parentDir': self.parentDir,
                                'projectDir': self.projectDir,
                                'analysis folder': 'Optimisation',
                                'cell type': self.cell_type,
                                'optimisation': True
                                }

            shape_space = self.uq_parallel(df, self.objective_vars, solver_dict, solver_args_dict, self.uq_config)

            # get uq_parameters
            uq_result_dict = {}
            for key in shape_space.keys():
                filename_slans = self.projectDir / fr'SimulationData\Optimisation\{key}\uq.json'
                filename_abci = self.projectDir / fr'SimulationData\ABCI\{key}\uq.json'
                if os.path.exists(filename_slans):  # and os.path.exists(filename_abci):
                    uq_result_dict[key] = []
                    with open(filename_slans, "r") as infile:
                        uq_d = json.load(infile)
                        for o in self.objectives:
                            if o[1] in ["Req", "freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]",
                                        "G [Ohm]", "Q []"]:
                                uq_result_dict[key].append(uq_d[o[1]]['expe'][0])
                                uq_result_dict[key].append(uq_d[o[1]]['stdDev'][0])
                                if o[0] == 'min':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] + 6 * uq_d[o[1]]['stdDev'][0])
                                elif o[0] == 'max':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] - 6 * uq_d[o[1]]['stdDev'][0])
                                else:
                                    # for equal, calculate |expected_value - design_value| + 6sigma
                                    uq_result_dict[key].append(
                                        np.abs(uq_d[o[1]]['expe'][0] - o[2]) + uq_d[o[1]]['stdDev'][0])

                if os.path.exists(filename_abci):
                    if key not in uq_result_dict:
                        uq_result_dict[key] = []

                    with open(filename_abci, "r") as infile:
                        uq_d = json.load(infile)
                        for o in self.objectives:
                            if o[1] not in ["Req", "freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]",
                                            "G [Ohm]", "Q []"]:
                                uq_result_dict[key].append(uq_d[o[1]]['expe'][0])
                                uq_result_dict[key].append(uq_d[o[1]]['stdDev'][0])
                                if o[0] == 'min':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] + 6 * uq_d[o[1]]['stdDev'][0])
                                elif o[0] == 'max':
                                    uq_result_dict[key].append(uq_d[o[1]]['expe'][0] - 6 * uq_d[o[1]]['stdDev'][0])

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

            df_uq = pd.DataFrame.from_dict(uq_result_dict, orient='index')
            print(uq_column_names)
            print(df_uq)
            df_uq.columns = uq_column_names
            df_uq.index.name = 'key'
            df_uq.reset_index(inplace=True)
            df = df.merge(df_uq, on='key', how='inner')

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

        # reset total rank
        df['total_rank'] = 0

        # rank shapes by objectives
        for i, obj in enumerate(self.objectives):
            if self.uq_config['option']:
                if obj[0] == "min":
                    df[f'rank_E[{obj[1]}] + 6*std[{obj[1]}]'] = df[fr'E[{obj[1]}] + 6*std[{obj[1]}]'].rank() * \
                                                                self.weights[i]
                elif obj[0] == "max":
                    df[f'rank_E[{obj[1]}] - 6*std[{obj[1]}]'] = df[fr'E[{obj[1]}] - 6*std[{obj[1]}]'].rank(
                        ascending=False) * self.weights[i]
                elif obj[0] == "equal":
                    df[fr'rank_|E[{obj[1]}] - {obj[2]}| + std[{obj[1]}]'] = df[
                                                                                fr'|E[{obj[1]}] - {obj[2]}| + std[{obj[1]}]'].rank() * \
                                                                            self.weights[i]

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
                elif obj[0] == "equal" and obj[1] != 'freq [MHz]':  # define properly later
                    df[f'rank_{obj[1]}'] = (df[obj[1]] - fshift).abs().rank() * self.weights[i]

                df[f'total_rank'] = df[f'total_rank'] + df[f'rank_{obj[1]}']

        # reorder
        tot = df.pop(f'total_rank')
        df[f'total_rank'] = tot / sum(self.weights)  # normalize by sum of weights

        # order shapes by rank
        df = df.sort_values(by=['total_rank'])
        df = df.reset_index(drop=True)

        # pareto condition
        reorder_indx, pareto_indx_list = self.pareto_front(df)

        # estimate convergence
        obj_error = []
        obj0 = self.objectives[0][1]
        for i, obj in enumerate(self.objectives):
            if i != 0:
                pareto_shapes = df.loc[pareto_indx_list, [obj0, obj[1]]]
                pareto_shapes_sorted = pareto_shapes.sort_values(obj0)
                print(pareto_shapes)
                f1 = np.linspace(min(pareto_shapes[obj0]), max(pareto_shapes[obj0]), self.n_interp)
                f2_interp = np.interp(f1, pareto_shapes_sorted[obj0], pareto_shapes_sorted[obj[1]])
                rel_error = np.linalg.norm(f2_interp - self.f2_interp[i]) / max(np.abs(f2_interp))
                obj_error.append(rel_error)

                self.f2_interp[i] = f2_interp

        # new error
        # stack previous and current pareto fronts
        if n == 0:
            pareto_shapes = df.loc[pareto_indx_list, self.objective_vars]
            self.pareto_history.append(pareto_shapes)
        else:
            pareto_shapes = df.loc[pareto_indx_list, self.objective_vars]
            pareto_stack = np.vstack([self.pareto_history[-1], pareto_shapes])

            # Compute Delaunay triangulation
            delaunay = Delaunay(pareto_stack)
            simplices = delaunay.simplices

            hypervolumes = self.calculate_hypervolumes(pareto_stack, simplices)
            self.err.append(sum(hypervolumes))

        if len(obj_error) != 0:
            self.interp_error.append(max(obj_error))
            self.interp_error_avg.append(np.average(self.interp_error))

        df = df.loc[reorder_indx, :]
        # reset index
        df = df.dropna().reset_index(drop=True)

        # update global
        self.df_global = df

        # check if df_global is empty
        if self.df_global.shape[0] == 0:
            error("Unfortunately, none survived the constraints and the program has to end. "
                  "Can't even say that this was a good run.")
            return

        # save dataframe
        filename = fr"{self.projectDir}\SimulationData\Optimisation\Generation{n}.xlsx"
        self.recursive_save(self.df_global, filename, reorder_indx)

        # birth next generation
        # crossover
        if len(df) > 1:
            df_cross = self.crossover(df, n, self.crossover_factor)
        else:
            df_cross = pd.DataFrame()

        # mutation
        df_mutation = self.mutation(df, n, self.mutation_factor)

        # chaos
        df_chaos = self.chaos(self.chaos_factor, n)

        # take elites from previous generation over to next generation
        df_ng = pd.concat([df_cross, df_mutation, df_chaos], ignore_index=True)

        # update dictionary
        self.df = df_ng

        n += 1
        print("=" * 80)
        if n < self.ng_max:
            bar.update(1)
            return self.ea(n, bar)
        else:
            bar.update(1)
            end = datetime.datetime.now()
            print("End time: ", end)
            plt.plot(self.interp_error, marker='P', label='max error')
            plt.plot(self.interp_error_avg, marker='X', label='avereage')
            plt.plot([x + 1 for x in range(len(self.err))], self.err, marker='o', label='convex hull vol')
            plt.yscale('log')
            plt.legend()

            plt.xlabel('Generation $n$')
            plt.ylabel(r"Pareto surface interp. error")
            plt.show()
            return

    def uq_parallel(self, df, objectives, solver_dict, solver_args_dict, uq_config):
        """

        Parameters
        ----------
        df: pd.DataFrame
            Pandas dataframe containing cavity geometry parameters
        objectives: list|ndarray
            List of objective functions
        solver_dict: dict
            Python dictionary of solver settings
        solver_args_dict: dict
            Python dictionary of solver arguments
        uq_config:
            Python dictionary of uncertainty quantification settings

        Returns
        -------

        """
        proc_count = uq_config['processes']

        # get geometric parameters
        df = df.loc[:, ['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"]]
        shape_space = {}

        df = df.set_index('key')
        for index, row in df.iterrows():
            rw = row.tolist()

            if self.cell_type.lower() == 'mid cell' or self.cell_type.lower() == 'mid-cell' or self.cell_type.lower() == 'mid_cell':
                shape_space[f'{index}'] = {'IC': rw, 'OC': rw, 'OC_R': rw}

            elif self.cell_type.lower() == 'mid-end cell' or self.cell_type.lower() == 'mid-end-cell' or self.cell_type.lower() == 'mid_end_cell':

                assert 'mid cell' in list(self.optimisation_config.keys()), \
                    ("If cell type is set as 'mid-end cell', the mid cell geometry parameters must "
                     "be provided in the optimisation_config dictionary.")
                assert len(self.optimisation_config['mid cell']) > 6, ("Incomplete mid cell geometry parameter. "
                                                                       "At least 7 geometric parameters "
                                                                       "[A, B, a, b, Ri, L, Req] required.")

                IC = self.optimisation_config['mid cell']
                # check if mid-cell is not a degenerate geometry
                df = tangent_coords(*np.array(IC)[0:8], 0)
                assert df[-2] == 1, ("The mid-cell geometry dimensions given result in a degenerate geometry. "
                                     "Please check.")
                shape_space[f'{index}'] = {'IC': IC, 'OC': rw, 'OC_R': rw}

            elif (self.cell_type.lower() == 'end-end cell' or self.cell_type.lower() == 'end-end-cell'
                  or self.cell_type.lower() == 'end_end_cell') or self.cell_type.lower() == 'end end cell':

                shape_space[f'{index}'] = {'IC': rw, 'OC': rw, 'OC_R': rw}

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

                solver_args_dict['ngsolvemevp']['proc'] = p
                solver_args_dict['abci']['proc'] = p
                service = mp.Process(target=uq,
                                     args=(processor_shape_space, objectives,
                                           solver_dict, solver_args_dict, uq_config))
                service.start()
                self.processes.append(service)
            except Exception as e:
                error(f"Exception in uq_parallel analysis:: {e}")

        for p in self.processes:
            p.join()

        self.processes = []
        return shape_space

    # @staticmethod
    # def uq_ngsolve(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args,
    #                select_solver='slans'):
    #     """
    #
    #     Parameters
    #     ----------
    #     key: str | int
    #         Cavity geomery identifier
    #     shape: dict
    #         Dictionary containing geometric dimensions of cavity geometry
    #     qois: list
    #         Quantities of interest considered in uncertainty quantification
    #     n_cells: int
    #         Number of cavity cells
    #     n_modules: int
    #         Number of modules
    #     n_modes: int
    #         Number of eigenmodes to be calculated
    #     f_shift: float
    #         Since the eigenmode solver uses the power method, a shift can be provided
    #     bc: int
    #         Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
    #         bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    #     pol: int {Monopole, Dipole}
    #         Defines whether to calculate for monopole or dipole modes
    #     parentDir: str | path
    #         Parent directory
    #     projectDir: str|path
    #         Project directory
    #
    #     Returns
    #     -------
    #
    #     """
    #
    #     if select_solver.lower() == 'slans':
    #         uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    #     else:
    #         uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'
    #
    #     err = False
    #     result_dict_eigen = {}
    #     eigen_obj_list = qois
    #     for o in qois:
    #         result_dict_eigen[o] = {'expe': [], 'stdDev': []}
    #
    #     rdim = n_cells * 3  # How many variables will be considered as random in our case 5
    #     degree = 1
    #
    #     #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    #     flag_stroud = 'stroud5'
    #
    #     if flag_stroud == 'stroud3':
    #         nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
    #         nodes_ = 2. * nodes_ - 1.
    #         # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    #     elif flag_stroud == 'stroud5':
    #         nodes_, weights_ = cn_leg_05_2(rdim)
    #     elif flag_stroud == 'cn_gauss':
    #         nodes_, weights_ = cn_gauss(rdim, 2)
    #     else:
    #         ic('flag_stroud==1 or flag_stroud==2')
    #         return 0
    #
    #     ic(nodes_)
    #     # save nodes
    #     data_table = pd.DataFrame(nodes_.T, columns=list(eigen_obj_list))
    #     data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')
    #
    #     #  mean value of geometrical parameters
    #     no_parm, no_sims = np.shape(nodes_)
    #
    #     Ttab_val_f = []
    #
    #     sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    #
    #     for i in range(no_sims):
    #         skip = False
    #         # perform checks on geometry
    #         ok = perform_geometry_checks(shape['IC'], shape['OC'])
    #         if not ok:
    #             err = True
    #             break
    #         fid = fr'{key}_Q{i}'
    #
    #         # skip analysis if folder already exists.
    #         if not skip:
    #             solver = ngsolve_mevp
    #             #  run model using SLANS or CST
    #             # # create folders for all keys
    #             solver.createFolder(fid, projectDir, subdir=sub_dir)
    #
    #             if "CELL TYPE" in shape.keys():
    #                 if shape['CELL TYPE'] == 'flattop':
    #                     # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
    #                     #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
    #                     try:
    #                         print(' in flattop')
    #                         solver.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
    #                                               n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
    #                                               beampipes=shape['BP'],
    #                                               parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
    #                                               mesh_args=mesh_args,
    #                                               deformation_params=nodes_[:, i])
    #                     except KeyError:
    #                         solver.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
    #                                               n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
    #                                               beampipes=shape['BP'],
    #                                               parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
    #                                               mesh_args=mesh_args,
    #                                               deformation_params=nodes_[:, i])
    #             else:
    #                 try:
    #                     solver.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
    #                                   n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
    #                                   parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args,
    #                                   deformation_params=nodes_[:, i])
    #                 except KeyError:
    #                     solver.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
    #                                   n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
    #                                   parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args,
    #                                   deformation_params=nodes_[:, i])
    #
    #         filename = uq_path / f'{fid}/monopole/qois.json'
    #         print(filename)
    #         if os.path.exists(filename):
    #             # params = fr.svl_reader(filename)
    #             # norm_length = 2 * n_cells * shape['IC'][5]
    #
    #             qois_result_dict = dict()
    #
    #             with open(filename) as json_file:
    #                 qois_result_dict.update(json.load(json_file))
    #
    #             qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
    #             # print_(qois_result)
    #             # sometimes some degenerate shapes are still generated and the solver returns zero
    #             # for the objective functions, such shapes are considered invalid
    #             for objr in qois_result:
    #                 if objr == 0:
    #                     # skip key
    #                     err = True
    #                     break
    #
    #             tab_val_f = qois_result
    #
    #             Ttab_val_f.append(tab_val_f)
    #         else:
    #             err = True
    #
    #     # # add original point
    #     # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    #     # params = fr.svl_reader(filename)
    #     # obj_result, tune_result = get_objectives_value(params, eigen_obj_list)
    #     # tab_val_f = obj_result
    #     # Ttab_val_f.append(tab_val_f)
    #     # save table
    #     data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    #     data_table.to_csv(uq_path / 'table.csv', index=False, sep='\t', float_format='%.32f')
    #
    #     print(np.atleast_2d(Ttab_val_f), weights_)
    #     if not err:
    #         v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)
    #
    #         # append results to dict
    #         for i, o in enumerate(eigen_obj_list):
    #             result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
    #             result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])
    #
    #             # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
    #             # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)
    #
    #         # plt.show()
    #         print(result_dict_eigen)
    #         with open(uq_path / fr"uq.json", 'w') as file:
    #             file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))
    #     else:
    #         print_(fr"There was a problem running UQ analysis for {key}")

    @staticmethod
    def calculate_hypervolumes(points, simplices):
        volumes = []
        for simplex in simplices:
            hull = ConvexHull(points[simplex])
            volumes.append(hull.volume)
        return volumes

    def plot_pareto(self, vars, which='last'):
        ################## plot 2d ##################################
        grid_results_folder = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'
        fig, axs = plt.subplot_mosaic([[0, 1, 2], [3, 4, 5]], figsize=(12, 3))

        columns_array = [['Epk/Eacc', 'Bpk/Eacc'], ['Epk/Eacc', 'R/Q'], ['Bpk/Eacc', 'R/Q']]
        columns_par_array = [['A', 'B'], ['a', 'b'], ['A', 'Ri']]
        cmap = matplotlib.colormaps['Pastel2_r']
        norm = matplotlib.colors.Normalize(vmin=0, vmax=49)
        ff = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData'
        for i, (columns, columns_par) in enumerate(zip(columns_array, columns_par_array)):
            for r in [49]:
                for lab, opt_code_result_folder, error_file in zip(['LHS', 'LHS2', 'Random'],
                                                                   [fr'{ff}\SLANS_LHS\Generation{r}.xlsx',
                                                                    fr'{ff}\SLANS_LHS2\Generation{r}.xlsx',
                                                                    fr'{ff}\SLANS_kwt_random1\Generation{r}.xlsx'],
                                                                   [fr'{ff}\SLANS_LHS2\inerp_error_and_average.txt',
                                                                    fr'{ff}\SLANS_LHS2\inerp_error_and_average.txt',
                                                                    fr'{ff}\SLANS_LHS2\inerp_error_and_average.txt']):
                    opt_code_result = pd.read_excel(opt_code_result_folder, 'Sheet1')

                    pareto_shapes = self.pareto_front(opt_code_result, columns, axs[i], show='none',
                                                      label=f"{lab}: g{r} ($n$={len(opt_code_result.index)}).",
                                                      kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'b'},
                                                      kwargs_pareto={'marker': 'o', 'mec': 'k', 'ms': 3},
                                                      # kwargs_pareto={'c': f'{matplotlib.colors.rgb2hex(cmap(norm(r)))}', 'marker': 'o',
                                                      #                'mec': 'k'}
                                                      )
                    # axs[i+3].plot(grid_results[columns[0]], grid_results[columns[1]], marker='o', ms=5, lw=0)
                    axs[i + 3].plot(opt_code_result[columns[0]], opt_code_result[columns[1]], marker='o', ms=5, lw=0)
                    # axs[i+3].plot(qmc_pareto_shapes[columns[0]], qmc_pareto_shapes[columns[1]], marker='o', c='r', label='qmc', ms=5, lw=0)
                    axs[i + 3].plot(pareto_shapes[columns[0]], pareto_shapes[columns[1]], marker='o', c='b', mec='k',
                                    label='ea', ms=5, lw=0)

                    # # load interpolation error
                    # error = pd.read_csv(error_file, header=None, sep='\s+')
                    # axs[3].plot(error[0], label=lab)
                    # axs[4].plot(error[1], label=lab)

            axs[i].set_xlabel(columns[0])
            axs[i].set_ylabel(columns[1])
            axs[i + 3].set_xlabel(columns_par[0])
            axs[i + 3].set_ylabel(columns_par[1])
        lines, labels = axs[0].get_legend_handles_labels()
        # axs[i].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=10, loc='lower left', mode='expand')
        fig.legend(*axs[0].get_legend_handles_labels(), loc="upper left", mode="expand", ncol=4)
        axs[3].legend()

        axs[3].set_yscale('log')
        axs[4].set_yscale('log')
        axs[3].set_xlabel('Interpolation error')
        # axs[3].set_ylabel()
        # plot error

        plt.tight_layout()
        plt.show()

        # ################### plot surface #########################
        # grid_results_folder = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'
        # opt_code_result_folder_lhs = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS_opt_KWT\Generation49.xlsx'
        # opt_code_result_folder_random = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS\Generation49.xlsx'
        #
        # grid_results = pd.read_excel(grid_results_folder, 'Sheet1')
        # opt_code_result_lhs = pd.read_excel(opt_code_result_folder_lhs, 'Sheet1')
        # opt_code_result_random = pd.read_excel(opt_code_result_folder_random, 'Sheet1')
        #
        # fig = plt.figure()
        # ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        # ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        # # ax3 = fig.add_subplot(1, 2, 3, projection='3d')
        # axs = [ax1, ax2]
        #
        # grid_pareto_surface = plot_pareto_surface(grid_results, ['Epk/Eacc', 'Bpk/Eacc', 'R/Q'], axs[0],
        #                                           {'cmap': 'gray', 'edgecolor': 'k'})
        # opt_pareto_surface_lhs = plot_pareto_surface(opt_code_result_lhs, ['Epk/Eacc', 'Bpk/Eacc', 'R/Q'], axs[1],
        #                                              {'cmap': 'RdBu', 'edgecolor': 'k'})
        # # opt_pareto_surface_random = plot_pareto_surface(opt_code_result_random, ['Epk/Eacc', 'Bpk/Eacc', 'R/Q'], axs[2],
        # #                                                 {'cmap': 'RdBu', 'edgecolor': 'k'})
        #
        # plt.show()

        # fig, axs = plt.subplot_mosaic([[0, 1, 2], [3, 4, 5]], figsize=(12, 5.5))
        # ani = animation.FuncAnimation(fig=fig, func=create_evolution_animation, frames=49, interval=1000)
        # ani.save(filename="D:\Dropbox\Quick presentation files/ffmpeg_example.gif", writer="pillow")
        # # plt.show()

    def pareto_front_(self, df, columns, ax=None, show='all', label='', kwargs_dataframe=None, kwargs_pareto=None):
        if kwargs_dataframe is None:
            kwargs_dataframe = {}
        if kwargs_pareto is None:
            kwargs_pareto = {}

        datapoints = df.loc[:, columns] * (-1)
        pareto = oapackage.ParetoDoubleLong()
        # ic(datapoints)
        for ii in range(0, datapoints.shape[0]):
            w = oapackage.doubleVector(tuple(datapoints.iloc[ii].values))
            pareto.addvalue(w, ii)
        # pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs
        poc = len(lst)  # number of pareto shapes
        reorder_idx = list(lst) + [i for i in range(len(df)) if
                                   i not in lst]  # reordered index putting pareto shapes first

        pareto_shapes = df.loc[lst, :]
        # sort pareto shapes in ascending x axis data
        pareto_shapes = pareto_shapes.sort_values(by=[columns[0]])

        if ax:
            if show == 'all':
                ax.scatter(df[columns[0]], df[columns[1]], **kwargs_dataframe)

            ax.plot(pareto_shapes[columns[0]], pareto_shapes[columns[1]], **kwargs_pareto, label=label)

        return pareto_shapes

    def plot_pareto_surface(self, df, columns, ax, kwargs=None):
        if kwargs is None:
            kwargs = {}
        pareto = self.pareto_front(df, columns)

        x, y, z = pareto['Epk/Eacc []'], pareto['Bpk/Eacc'], pareto['R/Q']
        xi, yi = np.meshgrid(np.linspace(min(x), max(x), 100),
                             np.linspace(min(y), max(y), 100))
        zi = griddata((x, y), z, (xi, yi), method='cubic')
        surf = ax.plot_surface(xi, yi, zi, antialiased=False, **kwargs)

        return surf

    def create_evolution_animation(self, frame):
        for ax_index in axs:
            axs[ax_index].clear()

        ################### plot 2d ##################################
        grid_results_folder = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'

        columns_array = [['Epk/Eacc', 'Bpk/Eacc'], ['Epk/Eacc', 'R/Q'], ['Bpk/Eacc', 'R/Q']]
        columns_par_array = [['A', 'B'], ['a', 'b'], ['A', 'Ri']]
        cmap = matplotlib.colormaps['Pastel2_r']
        norm = matplotlib.colors.Normalize(vmin=0, vmax=49)

        for i, (columns, columns_par) in enumerate(zip(columns_array, columns_par_array)):
            grid_results = pd.read_excel(grid_results_folder, 'Sheet1')
            qmc_pareto_shapes = pareto_front(grid_results, columns, axs[i], show='pareto',
                                             label=f'QMC \n({len(grid_results.index)} geoems.)',
                                             kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'grey'},
                                             kwargs_pareto={'c': 'k', 'marker': 'o', 'mec': 'k'})

            for r in [frame]:
                for lab, opt_code_result_folder, error_file in zip(['LHS2'],
                                                                   [
                                                                       fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS\Generation{r}.xlsx'],
                                                                   [
                                                                       fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS_LHS2\inerp_error_and_average.txt']):
                    opt_code_result = pd.read_excel(opt_code_result_folder, 'Sheet1')

                    pareto_shapes = pareto_front(opt_code_result, columns, axs[i], show='pareto',
                                                 label=f"{lab}: G{r} \n({len(opt_code_result.index)} geoms).",
                                                 kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'b'},
                                                 kwargs_pareto={'marker': 'o',
                                                                'mec': 'k'},
                                                 # kwargs_pareto={'c': f'{matplotlib.colors.rgb2hex(cmap(norm(r)))}', 'marker': 'o',
                                                 #                'mec': 'k'}
                                                 )
                    axs[i + 3].scatter(grid_results[columns_par[0]], grid_results[columns_par[1]], s=5)
                    axs[i + 3].scatter(opt_code_result[columns_par[0]], opt_code_result[columns_par[1]], s=5)
                    axs[i + 3].scatter(qmc_pareto_shapes[columns_par[0]], qmc_pareto_shapes[columns_par[1]], c='r',
                                       label='qmc', s=5)
                    axs[i + 3].scatter(pareto_shapes[columns_par[0]], pareto_shapes[columns_par[1]], c='b',
                                       edgecolor='k',
                                       label='ea', s=5)

                    # load interpolation error
                    error = pd.read_csv(error_file, header=None, sep='\s+')
                    # axs[3].plot(error[0], label=lab)
                    # axs[4].plot(error[1], label=lab)

            axs[i].set_xlabel(columns[0])
            axs[i].set_ylabel(columns[1])
            axs[i + 3].set_xlabel(columns_par[0])
            axs[i + 3].set_ylabel(columns_par[1])

        # axs[i].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=10, loc='lower left', mode='expand')
        axs[0].legend()
        axs[3].legend()

        # axs[3].set_yscale('log')
        # axs[4].set_yscale('log')
        # axs[3].set_xlabel('Interpolation error')
        # axs[3].set_ylabel()
        # plot error

        # plt.tight_layout()
        # plt.show()

    @staticmethod
    def stroud(p):
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
                outvar[i, 0] = np.dot(tab_var[:, i] ** 2, weights)
            stdDev = np.sqrt(outvar - expe ** 2)
        else:
            expe = 0
            stdDev = 0
            print('Cols_sims_no != No_weights')

        return list(expe.T[0]), list(stdDev.T[0])

    def run_tune_parallel(self, pseudo_shape_space, n_cells):
        # split shape_space for different processes/ MPI share process by rank
        keys = list(pseudo_shape_space.keys())

        # check if number of processors selected is greater than the number of keys in the pseudo shape space
        if self.processes_count > len(keys):
            self.processes_count = len(keys)

        shape_space_len = len(keys)
        share = floor(shape_space_len / self.processes_count)

        self.processes = []
        for p in range(self.processes_count):
            if True:
                if p < self.processes_count - 1:
                    proc_keys_list = keys[p * share:p * share + share]
                else:
                    proc_keys_list = keys[p * share:]

                self.overwriteFolder(p, self.projectDir)
                self.copyFiles(p, self.parentDir, self.projectDir)

                processor_shape_space = {}
                for key, val in pseudo_shape_space.items():
                    if key in proc_keys_list:
                        # check if folder exists
                        if not os.path.exists(
                                fr'{self.projectDir}\SimulationData\Optimisation\{key}\monopole\qois.json'):
                            processor_shape_space[key] = val

                tune_variable = self.tune_variable
                service = mp.Process(target=self.run_sequential,
                                     args=(processor_shape_space, "No", p, 33, self.parentDir,
                                           self.projectDir, "EA.json", 'Optimisation',
                                           tune_variable, ["Linear Interpolation", 1e-4, 10], self.cell_type, [],
                                           None, True, n_cells))
                service.start()
                self.processes.append(service)

        for p in self.processes:
            p.join()

        self.processes = []

    def run_slans_parallel(self, df, n_cells, fshift, pol):

        proc_count = self.ui.sb_Processors_Count.value()

        # get geometric parameters
        df = df.loc[:, ['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"]]
        shape_space = {}

        df = df.set_index('key')
        for index, row in df.iterrows():
            rw = row.tolist()
            if self.ui.cb_Cell_Type_Optimization.currentText() == 'End-Mid Cell':

                A_i = self.check_input(self.ui.le_A_i_opt.text())[0]
                B_i = self.check_input(self.ui.le_B_i_opt.text())[0]
                a_i = self.check_input(self.ui.le_a_i_opt.text())[0]
                b_i = self.check_input(self.ui.le_b_i_opt.text())[0]
                Ri_i = self.check_input(self.ui.le_Ri_i_opt.text())[0]
                L_i = self.check_input(self.ui.le_L_i_opt.text())[0]
                Req_i = self.check_input(self.ui.le_Req_i_opt.text())[0]
                alpha_i = self.check_input(self.ui.le_Alpha_opt.text())[0]

                IC = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, alpha_i]

                shape_space[f'{index}'] = {'IC': IC, 'OC': rw, 'OC_R': rw, 'BP': 'left'}
            else:
                shape_space[f'{index}'] = {'IC': rw, 'OC': rw, 'OC_R': rw, 'BP': 'left'}

        # split shape_space for different processes/ MPI share process by rank
        keys = list(shape_space.keys())

        # check if number of processors selected is greater than the number of keys in the pseudo shape space
        if proc_count > len(keys):
            proc_count = len(keys)

        shape_space_len = len(keys)
        share = floor(shape_space_len / proc_count)

        self.processes = []
        for pc in range(proc_count):
            if True:
                if pc < proc_count - 1:
                    proc_keys_list = keys[pc * share:pc * share + share]
                else:
                    proc_keys_list = keys[pc * share:]

                processor_shape_space = {}
                for key, val in shape_space.items():
                    if key in proc_keys_list:
                        # check if folder already exists
                        if not os.path.exists(
                                fr'{self.projectDir}\SimulationData\SLANS_opt\{key}\cavity_{fshift}_33.svl'):
                            processor_shape_space[key] = val

                service = mp.Process(target=self.run_sequential_slans,
                                     args=(n_cells, 1, processor_shape_space, n_cells, fshift, 33, pol,
                                           self.parentDir, self.projectDir, []))

                service.start()
                self.processes.append(service)

        for pc in self.processes:
            pc.join()

        self.processes = []
        return shape_space

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
        proc_count = self.processes_count

        # get geometric parameters
        df = df.loc[:, ['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"]]
        shape_space = {}

        df = df.set_index('key')
        for index, row in df.iterrows():
            rw = row.tolist()
            if self.cell_type.lower() == 'end-mid Cell':

                A_i = self.check_input(self.ui.le_A_i_opt.text())[0]
                B_i = self.check_input(self.ui.le_B_i_opt.text())[0]
                a_i = self.check_input(self.ui.le_a_i_opt.text())[0]
                b_i = self.check_input(self.ui.le_b_i_opt.text())[0]
                Ri_i = self.check_input(self.ui.le_Ri_i_opt.text())[0]
                L_i = self.check_input(self.ui.le_L_i_opt.text())[0]
                Req_i = self.check_input(self.ui.le_Req_i_opt.text())[0]
                alpha_i = self.check_input(self.ui.le_Alpha_opt.text())[0]

                IC = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, alpha_i]

                shape_space[f'{index}'] = {'IC': IC, 'OC': rw, 'OC_R': rw}
            else:
                shape_space[f'{index}'] = {'IC': rw, 'OC': rw, 'OC_R': rw}

        # print(shape_space)

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

    def generate_first_men(self, initial_points, n):

        if list(self.method.keys())[0] == "LHS":
            seed = self.method['LHS']['seed']
            if seed == '' or seed is None:
                seed = None

            columns = list(self.bounds.keys())
            dim = len(columns)
            l_bounds = np.array(list(self.bounds.values()))[:, 0]
            u_bounds = np.array(list(self.bounds.values()))[:, 1]

            const_var = []
            for i in range(dim - 1, -1, -1):
                if l_bounds[i] == u_bounds[i]:
                    const_var.append([columns[i], l_bounds[i]])
                    del columns[i]
                    l_bounds = np.delete(l_bounds, i)
                    u_bounds = np.delete(u_bounds, i)

            reduced_dim = len(columns)
            sampler = qmc.LatinHypercube(d=reduced_dim, scramble=False, seed=seed)
            _ = sampler.reset()
            sample = sampler.random(n=initial_points)
            self.discrepancy = qmc.discrepancy(sample)

            sample = qmc.scale(sample, l_bounds, u_bounds)

            df = pd.DataFrame()
            df['key'] = [f"G{n}_C{i}_P" for i in range(initial_points)]
            df[columns] = sample

            for i in range(len(const_var) - 1, -1, -1):
                df[const_var[i][0]] = np.ones(initial_points) * const_var[i][1]

            df['alpha_i'] = np.zeros(initial_points)
            df['alpha_o'] = np.zeros(initial_points)

            return df

        elif list(self.method.keys())[0] == "Sobol Sequence":
            seed = self.method['LHS']['seed']
            if seed == '' or seed is None:
                seed = None

            columns = list(self.bounds.keys())
            dim = len(columns)
            index = self.method["Sobol Sequence"]['index']
            l_bounds = np.array(list(list(self.bounds.values())))[:, 0]
            u_bounds = np.array(list(list(self.bounds.values())))[:, 1]

            const_var = []
            for i in range(dim - 1, -1, -1):
                if l_bounds[i] == u_bounds[i]:
                    const_var.append([columns[i], l_bounds[i]])
                    del columns[i]
                    l_bounds = np.delete(l_bounds, i)
                    u_bounds = np.delete(u_bounds, i)

            reduced_dim = len(columns)
            sampler = qmc.Sobol(d=reduced_dim, scramble=False, seed=seed)
            _ = sampler.reset()
            sample = sampler.random_base2(m=index)
            sample = qmc.scale(sample, l_bounds, u_bounds)

            df = pd.DataFrame()
            df['key'] = [f"G0_C{i}_P" for i in range(initial_points)]
            df[columns] = sample

            for i in range(len(const_var) - 1, -1, -1):
                df[const_var[i][0]] = np.ones(initial_points) * const_var[i][1]

            df['alpha_i'] = np.zeros(initial_points)
            df['alpha_o'] = np.zeros(initial_points)

            return df
        elif list(self.method.keys())[0] == "Random":
            data = {'key': [f"G0_C{i}_P" for i in range(initial_points)],
                    'A': random.sample(list(
                        np.linspace(list(self.bounds.values())[0][0], list(self.bounds.values())[0][1],
                                    initial_points * 2)),
                        initial_points),
                    'B': random.sample(list(
                        np.linspace(list(self.bounds.values())[1][0], list(self.bounds.values())[1][1],
                                    initial_points * 2)),
                        initial_points),
                    'a': random.sample(list(
                        np.linspace(list(self.bounds.values())[2][0], list(self.bounds.values())[2][1],
                                    initial_points * 2)),
                        initial_points),
                    'b': random.sample(list(
                        np.linspace(list(self.bounds.values())[3][0], list(self.bounds.values())[3][1],
                                    initial_points * 2)),
                        initial_points),
                    'Ri': random.sample(list(
                        np.linspace(list(self.bounds.values())[4][0], list(self.bounds.values())[4][1],
                                    initial_points * 2)),
                        initial_points),
                    'L': random.sample(list(
                        np.linspace(list(self.bounds.values())[5][0], list(self.bounds.values())[5][1],
                                    initial_points * 2)),
                        initial_points),
                    'Req': random.sample(list(
                        np.linspace(list(self.bounds.values())[6][0], list(self.bounds.values())[6][1] + 1,
                                    initial_points * 2)),
                        initial_points),
                    'alpha_i': np.zeros(initial_points),
                    'alpha_o': np.zeros(initial_points)}
            return pd.DataFrame.from_dict(data)
        elif list(self.method.keys())[0] == "Uniform":
            data = {'key': [f"G0_C{i}_P" for i in range(initial_points)],
                    'A': np.linspace(list(self.bounds.values())[0][0], self.bounds[0][1], initial_points),
                    'B': np.linspace(list(self.bounds.values())[1][0], self.bounds[1][1], initial_points),
                    'a': np.linspace(list(self.bounds.values())[2][0], self.bounds[2][1], initial_points),
                    'b': np.linspace(list(self.bounds.values())[3][0], self.bounds[3][1], initial_points),
                    'Ri': np.linspace(list(self.bounds.values())[4][0], list(self.bounds.values())[4][1],
                                      initial_points),
                    'L': np.linspace(list(self.bounds.values())[5][0], list(self.bounds.values())[5][1],
                                     initial_points),
                    'Req': np.linspace(list(self.bounds.values())[6][0], list(self.bounds.values())[6][1] + 1,
                                       initial_points),
                    'alpha_i': np.zeros(initial_points),
                    'alpha_o': np.zeros(initial_points)}
            return pd.DataFrame.from_dict(data)

    def process_constraints(self, constraints):
        processed_constraints = []

        for key, bounds in constraints.items():
            if isinstance(bounds, list):
                if len(bounds) == 2:
                    processed_constraints.append(fr'{key} > {bounds[0]}')
                    processed_constraints.append(fr'{key} < {bounds[1]}')
                else:
                    processed_constraints.append(fr'{key} > {bounds[0]}')
            else:
                processed_constraints.append(fr'{key} = {bounds}')

        return processed_constraints

    def process_objectives(self, objectives):
        # objectives, weights, constraints, n, ng_max = [["min", "Epk/Eacc"], ["min", "Bpk/Eacc"], ["max", "R/Q"]], [5, 1, 1], \
        #                                               ["Bpk/Eacc < 6.5", "freq > 400.58", 'freq < 400.99'], 1, 4
        # obj_vars = [obj[1] for obj in objectives]
        processed_objectives = []
        weights = []

        for i, obj in enumerate(objectives):
            if obj[1] == "ZL" or obj[1] == "ZT":
                goal = obj[0]
                freq_ranges = self.process_interval(obj[2])
                for f in freq_ranges:
                    processed_objectives.append([goal, f"{obj[1]} [max({f[0]}<f<{f[1]})]", f])
                    weights.append(1)
            else:
                goal = obj[0]
                if goal == 'equal':
                    processed_objectives.append(obj)
                else:
                    processed_objectives.append(obj)
                weights.append(1)

        return processed_objectives, weights

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
        G = 0.00948 * Q * (Freq / 1300)
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
        # if self.ui.cb_Tune_Variable.currentText() == "Req":
        #     tune_result.append(Req)
        # else:
        #     tune_result.append(L)
        #
        # tune_result.append(Freq)

        # append objective functions
        for o in obj:
            if o[1] in d.keys():
                objective.append(d[o[1]])

        return objective  # , tune_result

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
                abci_data_long = ABCIData(abci_data_dir, key, 0)
                abci_data_trans = ABCIData(abci_data_dir, key, 1)

                # trans
                x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                k_loss_trans = abci_data_trans.loss_factor['Transverse']

                if math.isnan(k_loss_trans):
                    print(f"Encountered an exception: Check shape {key}")
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
            if mon_interval is None:
                mon_interval = [0.0, 2e10]

            for key, value in d.items():
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
                            print("skipped, yp_mon = [], raise exception")
                            raise Exception()

                    processed_keys_mon.append(key)
                except:
                    print("skipped, yp_mon = []")
                    # for i, z_bound in enumerate(mon_interval):
                    #     Zmax_mon_list[i].append(-1)

            # print("2g", Zmax_mon_list)

            return Zmax_mon_list

        def get_Zmax_T(dip_interval=None):
            if dip_interval is None:
                dip_interval = [0.0, 2e10]

            for key, value in d.items():
                try:
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
                            print("skipped, yp_dip = [], raise exception")
                            raise Exception()

                    processed_keys_dip.append(key)
                except:
                    print("skipped, yp_dip = []")
                    # for i, z_bound in enumerate(dip_interval):
                    #     Zmax_dip_list[i].append(-1)

            return Zmax_dip_list

        def all(mon_interval, dip_interval):
            for key, value in d.items():
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

        ZL, ZT = [], []
        df_ZL, df_ZT = pd.DataFrame(), pd.DataFrame()
        for obj in self.objectives_unprocessed:
            if "ZL" in obj[1]:
                freq_range = self.process_interval(obj[2])
                for i in range(len(freq_range)):
                    Zmax_mon_list.append([])
                    xmax_mon_list.append([])
                    df_ZL[f"{obj[1]} [max({freq_range[i][0]}<f<{freq_range[i][1]})]"] = 0

                ZL = get_Zmax_L(freq_range)

            elif "ZT" in obj[1]:
                freq_range = self.process_interval(obj[2])

                for i in range(len(freq_range)):
                    Zmax_dip_list.append([])
                    xmax_dip_list.append([])
                    df_ZT[obj[1]] = 0

                ZT = get_Zmax_T(freq_range)

            elif obj[1] == "k_loss":
                pass
            elif obj[1] == "k_kick":
                pass

        # create dataframes from list
        print(processed_keys_mon, processed_keys_dip)
        df_ZL.loc[:, :] = np.array(ZL).T
        df_ZT.loc[:, :] = np.array(ZT).T
        df_ZL['key'] = processed_keys_mon
        df_ZT['key'] = processed_keys_dip

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
        elites = {}
        for i, o in enumerate(self.objectives):

            if self.uq_config['option']:
                if o[0] == "min":
                    elites[f'E[{o[1]}] + 6*std[{o[1]}]'] = df.sort_values(f'E[{o[1]}] + 6*std[{o[1]}]')
                elif o[0] == "max":
                    elites[f'E[{o[1]}] - 6*std[{o[1]}]'] = df.sort_values(f'E[{o[1]}] - 6*std[{o[1]}]', ascending=False)
                elif o[0] == "equal":
                    elites[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] = df.sort_values(
                        fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]')
            else:
                if o[0] == "min":
                    elites[f'{o[1]}'] = df.sort_values(f'{o[1]}')
                elif o[0] == "max":
                    elites[f'{o[1]}'] = df.sort_values(f'{o[1]}', ascending=False)
                elif o[0] == "equal":
                    # print(df[o[1]] - float(o[2]))
                    # print((df[o[1]] - float(o[2])).abs().sort_values())
                    elites[f'{o[1]}'] = df.sort_values(f'{o[1]}')

        obj_dict = {}
        for o in self.objectives:
            if self.uq_config['option']:
                if o[0] == 'min':
                    obj_dict[fr'E[{o[1]}] + 6*std[{o[1]}]'] = elites[fr'E[{o[1]}] + 6*std[{o[1]}]']
                elif o[0] == 'max':
                    obj_dict[fr'E[{o[1]}] - 6*std[{o[1]}]'] = elites[fr'E[{o[1]}] - 6*std[{o[1]}]']
                else:
                    obj_dict[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] = elites[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]']
            else:
                # if o[0] != 'equal':
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
        A_inf = ['All']
        B_inf = ['All']
        a_inf = ['All']
        b_inf = ['All']
        Ri_inf = ['All']
        L_inf = ['All']
        Req_inf = ['All']

        inf_dict = {"A": A_inf, "B": B_inf, "a": a_inf, "b": b_inf, "Ri": Ri_inf, "L": L_inf, "Req": Req_inf}
        for key, influence in inf_dict.items():
            if influence == [''] or influence == ['All']:
                if self.uq_config['option']:
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
                    # inf_dict[key] = [o[1] for o in self.objectives if o[0] != 'equal']
                    inf_dict[key] = self.objective_vars

        # print(elites.keys())
        # print(elites)
        # print(inf_dict)
        n_elites_to_cross = self.elites_to_crossover
        # print(n_elites_to_cross, len(elites), df.shape[0])
        # print(np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1))
        # print(np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1))
        # print(np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1))
        # print(np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1))
        # print(obj['mts dipole'])
        # print(obj['mts dipole'].loc[0])
        # print(obj['mts dipole'].loc[1])
        # print(obj['Epk/Eacc'].loc[1])
        # print(obj['mts dipole'].loc[np.random.randint(
        #     n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)])

        # print(sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["A"] for key in inf_dict["A"]]))
        # print(sum([obj[key].loc[np.random.randint(n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["A"] for key in inf_dict["A"]]))
        for i in range(f):
            # (<obj>[<rank>][<variable>] -> (b[c[1]][0]

            df_co.loc[i] = [f"G{generation}_C{i}_CO",
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["A"] for key
                                 in inf_dict["A"]]) / len(inf_dict["A"]),  # A
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["B"] for key
                                 in inf_dict["B"]]) / len(inf_dict["B"]),  # B
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["a"] for key
                                 in inf_dict["a"]]) / len(inf_dict["a"]),  # a
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["b"] for key
                                 in inf_dict["b"]]) / len(inf_dict["b"]),  # b
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["Ri"] for
                                 key in inf_dict["Ri"]]) / len(inf_dict["Ri"]),  # Ri
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["L"] for key
                                 in inf_dict["L"]]) / len(inf_dict["L"]),  # L
                            sum([obj[key].loc[np.random.randint(
                                n_elites_to_cross if n_elites_to_cross < df.shape[0] else df.shape[0] - 1)]["Req"] for
                                 key in inf_dict["Req"]]) / len(inf_dict["Req"]),
                            0,
                            0
                            ]
        return df_co

    def mutation(self, df, n, f):

        # get list based on mutation length
        if df.shape[0] < f:
            ml = np.arange(df.shape[0])
        else:
            ml = np.arange(f)

        df_ng_mut = pd.DataFrame(columns=['key', 'A', 'B', 'a', 'b', 'Ri', 'L', 'Req', "alpha_i", "alpha_o"])
        if list(self.bounds.values())[0][0] == list(self.bounds.values())[0][1]:
            df_ng_mut.loc[:, 'A'] = df.loc[ml, "A"]
        else:
            df_ng_mut.loc[:, 'A'] = df.loc[ml, "A"] * random.uniform(0.85, 1.5)

        if list(self.bounds.values())[1][0] == list(self.bounds.values())[1][1]:
            df_ng_mut.loc[:, 'B'] = df.loc[ml, "B"]
        else:
            df_ng_mut.loc[:, 'B'] = df.loc[ml, "B"] * random.uniform(0.85, 1.5)

        if list(self.bounds.values())[2][0] == list(self.bounds.values())[2][1]:
            df_ng_mut.loc[:, 'a'] = df.loc[ml, "a"]
        else:
            df_ng_mut.loc[:, 'a'] = df.loc[ml, "a"] * random.uniform(0.85, 1.5)

        if list(self.bounds.values())[3][0] == list(self.bounds.values())[3][1]:
            df_ng_mut.loc[:, 'b'] = df.loc[ml, "b"]
        else:
            df_ng_mut.loc[:, 'b'] = df.loc[ml, "b"] * random.uniform(0.85, 1.5)

        if list(self.bounds.values())[4][0] == list(self.bounds.values())[4][1]:
            df_ng_mut.loc[:, 'Ri'] = df.loc[ml, "Ri"]
        else:
            df_ng_mut.loc[:, 'Ri'] = df.loc[ml, "Ri"] * random.uniform(0.85, 1.5)

        if list(self.bounds.values())[5][0] == list(self.bounds.values())[5][1]:
            df_ng_mut.loc[:, 'L'] = df.loc[ml, "L"]
        else:
            df_ng_mut.loc[:, 'L'] = df.loc[ml, "L"] * random.uniform(0.85, 1.5)

        if list(self.bounds.values())[6][0] == list(self.bounds.values())[6][1]:
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
        df = self.generate_first_men(f, n)
        return df

    @staticmethod
    def remove_duplicate_values(d):
        temp = []
        res = dict()
        for key, val in d.items():
            if val not in temp:
                temp.append(val)
                res[key] = val
        return res

    @staticmethod
    def proof_filename(filepath):
        # check if extension is included
        if filepath.split('.')[-1] != 'json':
            filepath = f'{filepath}.json'

        return filepath

    def recursive_save(self, df, filename, pareto_index):
        styler = self.color_pareto(df, self.poc)
        try:
            styler.to_excel(filename)
            # df.to_excel(filename)
        except PermissionError:
            filename = filename.split('.xlsx')[0]
            filename = fr'{filename}_1.xlsx'
            self.recursive_save(df, filename, pareto_index)

    def pareto_front(self, df):

        # datapoints = np.array([reverse_list(x), reverse_list(y), reverse_list(z)])
        # reverse list or not based on objective goal: minimize or maximize
        # datapoints = [self.negate_list(df.loc[:, o[1]], o[0]) for o in self.objectives]

        if self.uq_config['option']:
            obj = []
            for o in self.objectives:
                if o[0] == 'min':
                    obj.append(fr'E[{o[1]}] + 6*std[{o[1]}]')
                elif o[0] == 'max':
                    obj.append(fr'E[{o[1]}] - 6*std[{o[1]}]')
                elif o[0] == 'equal':
                    obj.append(fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]')

            # print(obj)
            datapoints = df.loc[:, obj]
        else:
            datapoints = df.loc[:, self.objective_vars]

        for o in self.objectives:
            if o[0] == 'min':
                if self.uq_config['option']:
                    datapoints[fr'E[{o[1]}] + 6*std[{o[1]}]'] = datapoints[fr'E[{o[1]}] + 6*std[{o[1]}]'] * (-1)
                else:
                    datapoints[o[1]] = datapoints[o[1]] * (-1)
            elif o[0] == "equal":
                if self.uq_config['option']:
                    datapoints[fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] = datapoints[
                                                                             fr'|E[{o[1]}] - {o[2]}| + std[{o[1]}]'] * (
                                                                             -1)
                else:
                    datapoints[o[1]] = datapoints[o[1]] * (-1)
        # convert datapoints to numpy array

        pareto = oapackage.ParetoDoubleLong()
        for ii in range(0, datapoints.shape[0]):
            w = oapackage.doubleVector(tuple(datapoints.iloc[ii].values))
            pareto.addvalue(w, ii)
        pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs
        self.poc = len(lst)
        reorder_idx = list(lst) + [i for i in range(len(df)) if i not in lst]

        # return [optimal_datapoints[i, :] for i in range(datapoints.shape[0])]
        return reorder_idx, lst

    @staticmethod
    def check_input(s):
        # s = "range(16, 23, 10)"
        # s = "randrange(16, 23, 10)"
        # s = "[16, 23, 10]"
        # s = 1, 2, 3
        # s = 2

        if "r" in s and "rr" not in s:
            s = s.replace('r', '')
            # try:
            ll = eval(s)
            return np.linspace(ll[0], ll[1], ll[2])
            # except:
            #     print("Please check inputs.")
        elif "rr" in s:
            s = s.replace('rr', '')
            # try:
            ll = eval(s)
            ll = np.random.uniform(ll[0], ll[1], ll[2])
            return ll
            # except:
            #     print("Please check inputs.")
        else:
            # try:
            ll = eval(s)
            if isinstance(ll, int) or isinstance(ll, float):
                ll = [ll]
            return ll
            # except:
            #     print("Please check inputs.")

        # return 1

    @staticmethod
    def negate_list(ll, arg):
        if arg == 'max':
            return ll  # to find the pareto maxima
        else:
            return [-x for x in ll]  # to find the pareto minima

    @staticmethod
    def overwriteFolder(invar, projectDir):
        path = fr"{projectDir}\SimulationData\SLANS_opt\_process_{invar}"

        if os.path.exists(path):
            shutil.rmtree(path)
            dir_util._path_created = {}

        os.makedirs(path)

    @staticmethod
    def copyFiles(invar, parentDir, projectDir):
        src = fr"{parentDir}\exe\SLANS_exe"
        dst = fr"{projectDir}\SimulationData\SLANS_opt\_process_{invar}\SLANS_exe"

        dir_util.copy_tree(src, dst)

    @staticmethod
    def run_sequential(pseudo_shape_space_proc, resume, p, bc, parentDir, projectDir, filename, sim_folder,
                       tune_variable, iter_set, cell_type, progress_list, convergence_list, save_last, n_cells):
        tuner.tune_ngsolve(pseudo_shape_space_proc, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                           sim_folder=sim_folder, tune_variable=tune_variable, iter_set=iter_set,
                           cell_type=cell_type, progress_list=progress_list, convergence_list=convergence_list,
                           save_last=save_last,
                           n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    @staticmethod
    def run_sequential_slans(n_cell, n_modules, processor_shape_space, n_modes, f_shift, bc, pol, parentDir, projectDir,
                             progress_list, sub_dir='', uq_config=None, solver='slans', mesh=None):
        """
        Runs a single instance of SLANS (eigenmode analysis)
        Parameters
        ----------
        n_cell: int
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
        uq_config: None | dict
            Provides inputs required for uncertainty quantification. Default is None and disables uncertainty quantification.
        solver: str {slans, native}
            Select the eigenmode solver to use. Default is the SLANS solver
        mesh: list [Jxy, Jxy_bp, Jxy_bp_y]
            Mesh definition for logical mesh:
            Jxy -> Number of elements of logical mesh along JX and JY
            Jxy_bp -> Number of elements of logical mesh along JX in beampipe
            Jxy_bp_y -> Number of elements of logical mesh along JY in beampipe

        Returns
        -------

        """
        progress = 0

        for key, shape in processor_shape_space.items():
            if solver.lower() == 'slans':

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

                # # create folders for all keys
                slans_geom_seq.createFolder(f"{key}_n{n_cell}", projectDir, subdir=sub_dir, opt=True)

                if 'OC_R' in shape.keys():
                    write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC_R'],
                                        projectDir=projectDir, cell_type="None", opt=True)
                    slans_geom_seq.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                          n_modes=n_modes, fid=f"{key}_n{n_cell}", f_shift=f_shift,
                                          bc=bc, pol=pol, beampipes=shape['BP'],
                                          parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                          expansion=expansion, expansion_r=expansion_r, mesh=mesh, opt=True)
                else:
                    write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC'],
                                        projectDir=projectDir, cell_type="None", opt=True)
                    slans_geom_seq.cavity(n_cell, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                          n_modes=n_modes, fid=f"{key}_n{n_cell}", f_shift=f_shift,
                                          bc=bc, pol=pol, beampipes=shape['BP'],
                                          parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                          expansion=expansion, expansion_r=expansion_r, mesh=mesh, opt=True)

                # # run UQ
                # if uq_config:
                #     uq(key, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
                #        n_cells=n_cell, n_modules=n_modules, n_modes=n_modes,
                #        f_shift=f_shift, bc=bc, pol='Monopole', parentDir=parentDir, projectDir=projectDir)
                #
                # print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')
                #
                # # update progress
                # progress_list.append(progress + 1)

            # else:
            #     # run own eigenmode code
            #     folder = projectDir / fr'SimulationData\NativeEig'
            #     mod = Model(folder=folder, name=f"{key}", parent_dir=parentDir)
            #
            #     try:
            #         # convert first to m.
            #         mid_cell = np.array(shape['IC']) * 1e-3
            #         end_cell_left = np.array(shape['OC']) * 1e-3
            #         end_cell_right = np.array(shape['OC_R']) * 1e-3
            #
            #         mod.run(n_cells, mid_cell, end_cell_left, end_cell_right, beampipe=shape['BP'],
            #                 req_mode_num=int(n_modes), plot=False)
            #     except KeyError:
            #         # convert first to m.
            #         mid_cell = np.array(shape['IC']) * 1e-3
            #         end_cell_left = np.array(shape['OC']) * 1e-3
            #         end_cell_right = np.array(shape['OC']) * 1e-3
            #
            #         mod.run(n_cells, mid_cell, end_cell_left, end_cell_right, beampipe=shape['BP'],
            #                 req_mode_num=int(n_modes), plot=False)

    @staticmethod
    def process_interval(interval_list):
        interval = []
        for i in range(len(interval_list) - 1):
            interval.append([interval_list[i], interval_list[i + 1]])

        return interval

    def mode_trap_score(self, invar, mode_list_input):
        mt_score_list = []
        for mode in mode_list_input:
            x, y = [], []
            path = os.path.join(os.getcwd(),
                                r"SLANS_data_C{}\Cavity{}\cavity_mm_{}.af".format(mid_cell_cid, invar, mode))

            with open(path, 'r') as f:
                for ll in f.readlines():
                    ll = ll.strip()
                    x.append(float(ll.split(' ')[0]))
                    y.append(float(ll.split(' ')[1]))

            y_abs = [abs(t) for t in y]
            # plt.plot(x, y_abs)

            y_flat = []
            y_mx = max(y_abs)

            for v in y_abs:
                if v <= 0.1 * y_mx:
                    y_flat.append(0)
                else:
                    y_flat.append(v)

            # plt.plot(x, y_flat)

            peaks, _ = sps.find_peaks(y_flat)

            # calculate mt_score
            y_peaks = np.array(y_flat)[peaks]

            plt.scatter(np.array(x)[peaks], y_peaks, label="{},{}".format(invar, mode))

            if y_peaks[0] != 0:
                y_peaks_norm = [abs(t) / y_peaks[0] for t in y_peaks]

                print("\t\t\t", y_peaks_norm)

                if 1.05 * y_peaks[0] >= max(y_peaks):
                    mt_score = -(sum(y_peaks_norm) - 1) / (len(y_peaks_norm) - 1)
                else:
                    mt_score = (sum(y_peaks_norm) - 1) / (len(y_peaks_norm) - 1)

                print("\t\t\t", 1 / mt_score)
                mt_score_list.append(mt_score)
                print("MODE LIST:: ", mode_list_input)
            else:
                print("\t\t\tDividing by zero, skipped!")

        return mt_score_list

    @staticmethod
    def color_pareto(df, no_pareto_optimal):
        def color(row):
            # if row.isnull().values.any():
            if row.iloc[0] in df['key'].tolist()[0:no_pareto_optimal]:
                return ['background-color: #6bbcd1'] * len(row)
            return [''] * len(row)

        # Save Styler Object for Later
        styler = df.style
        # Apply Styles (This can be chained or on separate lines)
        styler.apply(color, axis=1)
        # Export the styler to excel
        return styler


class Cavity:
    """
    Command Line Interface module for running analysis.

    .. note::

       Still under development so some functions might not work properly
    """

    def __init__(self, n_cells, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none', name='cavity',
                 plot_label=None):
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

        self.uq_fm_results = None
        self.mesh = None
        if plot_label is None:
            self.plot_label = name
        else:
            self.plot_label = plot_label

        self.n_modes = n_cells + 1
        self.n_modules = 1
        self.folder = None
        self.bc = 33
        self.name = name
        self.n_cells = n_cells
        self.mid_cell = mid_cell[:7]
        self.end_cell_left = end_cell_left[:7]
        self.end_cell_right = end_cell_right[:7]
        self.beampipe = beampipe
        self.no_of_modules = 1
        self.eigenmode_qois = {}
        self.custom_eig_qois = {}
        self.abci_qois = {}
        self.slans_tune_res = {}
        self.wake_op_points = {}
        self.convergence_list = []
        self.eigenmode_tune_res = {}
        self.operating_points = None

        # slans results
        self.R_Q, self.k_fm, self.GR_Q, self.op_freq, self.e, self.b, \
            self.G, self.ff, self.k_cc, self.axis_field, self.surface_field = [0 for _ in range(11)]

        # abci results
        self.k_fm, self.k_loss, self.k_kick, self.phom, self.sigma, self.I0 = [{} for _ in range(6)]

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

        self.A, self.B, self.a, self.b, self.Ri, self.L, self.Req = self.mid_cell[:7]
        self.A_el, self.B_el, self.a_el, self.b_el, self.Ri_el, self.L_el, self.Req_el = self.end_cell_left[:7]
        self.A_er, self.B_er, self.a_er, self.b_er, self.Ri_er, self.L_er, self.Req_er = self.end_cell_right[:7]

        # get geometric parameters
        self.shape_space = {
            "IC": self.mid_cell[:7],
            "OC": self.end_cell_left[:7],
            "OC_R": self.end_cell_right[:7],
            "BP": beampipe
        }
        self.shape_space_multicell = {}
        self.to_multicell()  # <- get multicell representation

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

    def set_plot_label(self, plot_label):
        """
        Set cavity plot label

        Parameters
        ----------
        plot_label: str
            Cavity plot label

        Returns
        -------

        """

        if plot_label is None:
            self.plot_label = self.name
        else:
            self.plot_label = plot_label

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

    def sweep(self, sweep_config, which='eigenmode', how='independent'):
        self.sweep_results = {}
        if how == 'cross':
            pass
        else:
            for key, interval_def in tqdm(sweep_config.items()):
                # save nominal variable value form shape space
                current_var = self.shape_space['IC'][VAR_TO_INDEX_DICT[key]]
                par_vals = np.linspace(interval_def[0], interval_def[1], interval_def[2], endpoint=True)
                self.sweep_results[f'{key}'] = {}
                for val in par_vals:
                    if which == 'eigenmode':
                        # change value
                        self.shape_space['IC'][VAR_TO_INDEX_DICT[key]] = val
                        res = self.run_eigenmode()
                        if res:
                            self.sweep_results[f'{key}'][val] = copy.deepcopy(self.eigenmode_qois)

                # replace initial value
                self.shape_space['IC'][VAR_TO_INDEX_DICT[key]] = current_var

    def check_uq_config(self, uq_config):
        cell_type = uq_config['cell type']
        delta = uq_config['delta']
        method = uq_config['method']
        uq_vars = uq_config['variables']
        assert len(uq_vars) == len(delta), error('Ensure number of variables equal number of deltas')

        rdim = len(uq_vars)
        degree = 1

        if method[1].lower() == 'stroud3':
            nodes, weights, bpoly = quad_stroud3(rdim, degree)
            nodes = 2. * nodes - 1.
            # nodes, weights = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
        elif method[1].lower() == 'stroud5':
            nodes, weights = cn_leg_05_2(rdim)
        elif method[1].lower() == 'gaussian':
            nodes, weights = cn_gauss(rdim, 2)
        elif method[1].lower() == 'lhs':
            sampler = qmc.LatinHypercube(d=rdim)
            _ = sampler.reset()
            nsamp = uq_config['integration'][2]
            sample = sampler.random(n=nsamp)

            l_bounds = [-1 for _ in range(len(uq_vars))]
            u_bounds = [1 for _ in range(len(uq_vars))]
            sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

            nodes, weights = sample_scaled.T, np.ones((nsamp, 1))
        elif method[0].lower() == 'from file':
            if len(method) == 2:
                nodes = pd.read_csv(method[1], sep='\s+').iloc[:, method[1]]
            else:
                nodes = pd.read_csv(method[1], sep='\s+')

            nodes = nodes.to_numpy().T
            weights = np.ones((nodes.shape[1], 1))
        else:
            # issue warning
            warning('Integration method not recognised. Defaulting to Stroud3 quadrature rule!')
            nodes, weights, bpoly = quad_stroud3(rdim, degree)
            nodes = 2. * nodes - 1.

        perturbed_cell_node = np.array(cell_node)
        no_parm, no_sims = np.shape(nodes)
        if delta is None:
            delta = [0.05 for _ in range(len(uq_vars))]

        for i in range(no_sims):
            skip = False
            for j, uq_var in enumerate(uq_vars):
                uq_var_indx = VAR_TO_INDEX_DICT[uq_var]
                perturbed_cell_node[uq_var_indx] = cell_node[uq_var_indx] * (1 + delta[j] * nodes[j, i])

            if cell_type.lower() == 'mid cell' or cell_type.lower() == 'mid-cell' or cell_type.lower() == 'mid_cell':
                cell_node = shape['IC']
                mid = perturbed_cell_node
                left = perturbed_cell_node
                right = perturbed_cell_node
                beampipes = 'none'
            elif cell_type.lower() == 'mid-end cell' or cell_type.lower() == 'mid-end-cell' or cell_type.lower() == 'mid_end_cell':
                mid = shape['IC']
                left = shape['IC']
                right = perturbed_cell_node
                beampipes = 'right'
            elif (cell_type.lower() == 'end-end cell' or cell_type.lower() == 'end-end-cell'
                  or cell_type.lower() == 'end_end_cell') or cell_type.lower() == 'end end cell':
                mid = perturbed_cell_node
                left = perturbed_cell_node
                right = perturbed_cell_node
                beampipes = 'right'
            else:
                mid = perturbed_cell_node
                left = perturbed_cell_node
                right = perturbed_cell_node
                beampipes = 'both'

            enforce_Req_continuity(mid, left, right, cell_type)

    def study_convergence(self, step_refinement, max_steps):
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
        cell_type: {'mid cell', 'end-mid cell', 'mid-end cell', 'end-end cell', 'single cell'}
            Type of cell to tune
        tune_variable: {'Req', 'L'}
            Tune variable. Currently supports only the tuning of the equator radius ``Req`` and half-cell length ``L``

        Returns
        -------

        """

        for _ in tqdm([1]):
            iter_set = ['Linear Interpolation', TUNE_ACCURACY, 10]

            if freq is None:
                # calculate freq from mid cell length
                beta = 1
                freq = beta * c0 / (4 * self.mid_cell[5])
                info("Calculated freq from mid cell half length: ", freq)

            # create new shape space based on cell type
            # if cell_type.lower() == 'mid cell':
            shape_space = {
                f'{self.name}':
                    {
                        'IC': self.shape_space['IC'],
                        'OC': self.shape_space['OC'],
                        'OC_R': self.shape_space['OC_R'],
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
                                            progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
                else:
                    if run_tune.lower() == 'y':
                        # copy files required for simulation
                        self._overwriteFolder(proc, self.folder, self.name)
                        self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                        self.run_tune_ngsolve(shape_space, resume, proc, self.bc,
                                              SOFTWARE_DIRECTORY, self.folder, self.name,
                                              tune_variable, iter_set, cell_type,
                                              progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_ngsolve_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
            else:
                if solver.lower() == 'slans':
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_slans(shape_space, resume, proc, self.bc,
                                        SOFTWARE_DIRECTORY, self.folder, self.name, tuner,
                                        tune_variable, iter_set, cell_type,
                                        progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Oops! Something went wrong. Could not find the tune results. Please run tune again.")
                else:
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_ngsolve(shape_space, resume, proc, self.bc,
                                          SOFTWARE_DIRECTORY, self.folder, self.name,
                                          tune_variable, iter_set, cell_type,
                                          progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)
                    try:
                        self.get_ngsolve_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Oops! Something went wrong. Could not find the tune results. Please run tune again.")

    @staticmethod
    def run_tune_slans(shape, resume, p, bc, parentDir, projectDir, filename, tuner_option,
                       tune_variable, iter_set, cell_type, progress_list, convergence_list, n_cells):
        tuner.tune(shape, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                   tuner_option=tuner_option, tune_variable=tune_variable, iter_set=iter_set,
                   cell_type=cell_type,
                   progress_list=progress_list, convergence_list=convergence_list,
                   save_last=True,
                   n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    @staticmethod
    def run_tune_ngsolve(shape, resume, p, bc, parentDir, projectDir, filename,
                         tune_variable, iter_set, cell_type, progress_list, convergence_list, n_cells):
        tuner.tune_ngsolve(shape, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                           tune_variable=tune_variable, iter_set=iter_set,
                           cell_type=cell_type, sim_folder='Optimisation',
                           progress_list=progress_list, convergence_list=convergence_list,
                           save_last=True,
                           n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    def run_eigenmode(self, solver='ngsolve', freq_shift=0, boundary_cond=None, subdir='', uq_config=None):
        """
        Run eigenmode analysis on cavity

        Parameters
        ----------
        solver: {'SLANS', 'NGSolve'}
            Solver to be used. Native solver is still under development. Results are not as accurate as that of SLANS.
        freq_shift:
            Frequency shift. Eigenmode solver searches for eigenfrequencies around this value
        boundary_cond: int
            Boundary condition of left and right cell/beampipe ends
        subdir: str
            Sub directory to save results to
        uq_config: None | dict
            Provides inputs required for uncertainty quantification. Default is None and disables uncertainty quantification.

        Returns
        -------

        """

        for _ in tqdm([1], file=sys.stdout):
            if boundary_cond:
                self.bc = boundary_cond

            if solver.lower() == 'slans':
                self._run_slans(self.name, self.n_cells, self.n_modules, self.shape_space, self.n_modes, freq_shift,
                                self.bc, SOFTWARE_DIRECTORY, self.folder, sub_dir='', uq_config=uq_config)
                # # load quantities of interest
                # try:
                #     self.get_eigenmode_qois('SLANS')
                #     if uq_config:
                #         # load uq result
                #         self.get_uq_fm_results(fr"{self.folder}\SimulationData\SLANS\{self.name}\uq.json")
                #     return True
                # except FileNotFoundError:
                #     error("Could not find the eigenmode results. Please rerun eigenmode analysis.")
                #     return False
            else:
                uq_cell_complexity = 'simplecell'
                if 'cell complexity' in uq_config.keys():
                    uq_cell_complexity = uq_config['cell complexity']

                if uq_cell_complexity.lower() == 'multicell':
                    self._run_ngsolve(self.name, self.n_cells, self.n_modules, self.shape_space_multicell, self.n_modes,
                                      freq_shift, self.bc,
                                      SOFTWARE_DIRECTORY, self.folder, sub_dir='', uq_config=uq_config)
                else:
                    self._run_ngsolve(self.name, self.n_cells, self.n_modules, self.shape_space, self.n_modes,
                                      freq_shift, self.bc,
                                      SOFTWARE_DIRECTORY, self.folder, sub_dir='', uq_config=uq_config)

                # load quantities of interest
                # try:
                #     self.get_eigenmode_qois('NGSolveMEVP')
                #     if uq_config:
                #         self.get_uq_fm_results(fr"{self.folder}\SimulationData\NGSolveMEVP\{self.name}\uq.json")
                #     return True
                # except FileNotFoundError:
                #     error("Could not find eigenmode results. Please rerun eigenmode analysis.")
                #     return False

    def run_wakefield(self, MROT=2, MT=10, NFS=10000, wakelength=50, bunch_length=25,
                      DDR_SIG=0.1, DDZ_SIG=0.1, WG_M=None, marker='', operating_points=None, solver='ABCI'):
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

        if operating_points is None:
            wp_dict = {}
        exist = False

        if not exist:
            if solver == 'ABCI':
                self._run_abci(self.name, self.n_cells, self.n_modules, self.shape_space,
                               MROT=MROT, MT=MT, NFS=NFS, UBT=wakelength, bunch_length=bunch_length,
                               DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG,
                               parentDir=SOFTWARE_DIRECTORY, projectDir=self.folder, WG_M=WG_M, marker=marker,
                               operating_points=operating_points, freq=self.freq, R_Q=self.R_Q)

                try:
                    self.get_abci_data()
                    self.get_abci_qois()
                except FileNotFoundError:
                    error("Could not find the abci wakefield results. Please rerun wakefield analysis.")

        else:
            try:
                self.get_abci_data()
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
    def _run_slans(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir, sub_dir='',
                   uq_config=None):
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
        if uq_config:
            uq(name, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
               n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
               f_shift=f_shift, bc=bc, parentDir=parentDir, projectDir=projectDir)

        done(f'Done with Cavity {name}. Time: {time.time() - start_time}')

    @staticmethod
    def _run_ngsolve(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir, sub_dir='',
                     uq_config=None):
        parallel = True
        start_time = time.time()
        # create folders for all keys
        ngsolve_mevp.createFolder(name, projectDir, subdir=sub_dir)

        if 'OC_R' in shape.keys():
            OC_R = 'OC_R'
        else:
            OC_R = 'OC'

        # ngsolve_mevp.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape[OC_R],
        #                     n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
        #                     parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
        print(shape)
        if "CELL TYPE" in shape.keys():
            if shape['CELL TYPE'] == 'flattop':
                # write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC_R'],
                #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())

                ngsolve_mevp.cavity_flattop(n_cells, n_modules, shape['IC'], shape['OC'], shape[OC_R],
                                            n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc,
                                            beampipes=shape['BP'],
                                            parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

            elif shape['CELL TYPE'] == 'multicell':
                # write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC_R'],
                #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
                ngsolve_mevp.cavity_multicell(n_cells, n_modules, shape['IC'], shape['OC'], shape[OC_R],
                                              n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc,
                                              beampipes=shape['BP'],
                                              parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
        else:
            # write_cst_paramters(f"{key}_n{n_cell}", shape['IC'], shape['OC'], shape['OC_R'],
            #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())

            ngsolve_mevp.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape[OC_R],
                                n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

        # run UQ
        if uq_config:
            objectives = uq_config['objectives']
            solver_dict = {'ngsolvemevp': ngsolve_mevp}
            solver_args_dict = {'ngsolvemevp':
                                    {'n_cells': n_cells, 'n_modules': n_modules, 'f_shift': f_shift, 'bc': bc,
                                     'beampipes': shape['BP']
                                     },
                                'parentDir': parentDir,
                                'projectDir': projectDir,
                                'analysis folder': 'NGSolveMEVP',
                                'cell type': 'mid cell',
                                'optimisation': False
                                }

            uq_cell_complexity = 'simplecell'
            if 'cell complexity' in uq_config.keys():
                uq_cell_complexity = uq_config['cell complexity']

            if not parallel:
                # convert to proper shape space
                shape_space = {name: shape}
                if uq_cell_complexity == 'multicell':
                    uq_multicell(shape_space, objectives, solver_dict, solver_args_dict, uq_config)
                else:
                    uq(shape_space, objectives, solver_dict, solver_args_dict, uq_config)

            else:
                shape_space = {name: shape}
                if uq_cell_complexity == 'multicell':
                    uq_ngsolve_parallel_multicell(shape_space, objectives, solver_dict, solver_args_dict, uq_config)
                else:
                    uq_ngsolve_parallel(name, shape,
                                        objectives,
                                        n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
                                        f_shift=f_shift, bc=bc, pol='monopole', parentDir=parentDir,
                                        projectDir=projectDir,
                                        mesh_args=None,
                                        select_solver='ngsolvemevp')

        done(f'Done with Cavity {name}. Time: {time.time() - start_time}')

    # @staticmethod
    # def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir, mesh_args,
    #        select_solver='slans'):
    #     """
    #
    #     Parameters
    #     ----------
    #     key: str | int
    #         Cavity geomery identifier
    #     shape: dict
    #         Dictionary containing geometric dimensions of cavity geometry
    #     qois: list
    #         Quantities of interest considered in uncertainty quantification
    #     n_cells: int
    #         Number of cavity cells
    #     n_modules: int
    #         Number of modules
    #     n_modes: int
    #         Number of eigenmodes to be calculated
    #     f_shift: float
    #         Since the eigenmode solver uses the power method, a shift can be provided
    #     bc: int
    #         Boundary conditions {1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal}
    #         bc=33 means `Magnetic Wall En = 0` boundary condition at both ends
    #     pol: int {Monopole, Dipole}
    #         Defines whether to calculate for monopole or dipole modes
    #     parentDir: str | path
    #         Parent directory
    #     projectDir: str|path
    #         Project directory
    #
    #     Returns
    #     -------
    #
    #     """
    #
    #     if select_solver.lower() == 'slans':
    #         uq_path = projectDir / fr'SimulationData\SLANS\{key}'
    #     else:
    #         uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'
    #
    #     err = False
    #     result_dict_eigen = {}
    #     eigen_obj_list = qois
    #     for o in qois:
    #         result_dict_eigen[o] = {'expe': [], 'stdDev': []}
    #
    #     # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    #     p_true = shape['IC'][0:5]
    #     # p_true_el = shape['OC'][0:5]
    #     print(shape)
    #
    #     rdim = len(p_true)  # + len(p_true_el)  # How many variabels will be considered as random in our case 5
    #     degree = 1
    #
    #     flag_stroud = 'load_from_file'
    #     if flag_stroud == 'stroud3':
    #         nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
    #         nodes_ = 2. * nodes_ - 1.
    #         # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
    #     elif flag_stroud == 'stroud5':
    #         nodes_, weights_ = cn_leg_05_2(rdim)
    #     elif flag_stroud == 'cn_gauss':
    #         nodes_, weights_ = cn_gauss(rdim, 2)
    #     elif flag_stroud == 'load_from_file':
    #         nodes_ = pd.read_csv(fr'C:\Users\sosoho\DakotaProjects\Cavity\C3794_cubature5_5\sim_result_table.dat',
    #                              sep='\s+').iloc[:, 2:7]
    #         nodes_ = nodes_.to_numpy().T
    #         weights_ = np.ones((nodes_.shape[1], 1))
    #     else:
    #         ic('flag_stroud==1 or flag_stroud==2')
    #         return 0
    #     # save nodes
    #     data_table = pd.DataFrame(nodes_.T)
    #     data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')
    #
    #     #  mean value of geometrical parameters
    #     p_init = np.zeros(np.shape(p_true))
    #     # p_init_el = np.zeros(np.shape(p_true_el))
    #
    #     no_parm, no_sims = np.shape(nodes_)
    #     delta = 0.01  # or 0.1
    #
    #     Ttab_val_f = []
    #
    #     sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    #     par_end = shape['OC']
    #
    #     for i in range(no_sims):
    #         skip = False
    #         if flag_stroud == 'load_from_file':
    #             p_init[0] = nodes_[0, i]  # <- A
    #             p_init[1] = nodes_[1, i]  # <- B
    #             p_init[2] = nodes_[2, i]  # <- a
    #             p_init[3] = nodes_[3, i]  # <- b
    #             p_init[4] = nodes_[4, i]  # <- Ri
    #         else:
    #             #
    #             # p_init[0] = p_true[0] * (1 + delta * nodes_[0, i])  # <- A
    #             # p_init[1] = p_true[1] * (1 + delta * nodes_[1, i])  # <- B
    #             # p_init[2] = p_true[2] * (1 + delta * nodes_[2, i])  # <- a
    #             # p_init[3] = p_true[3] * (1 + delta * nodes_[3, i])  # <- b
    #             # p_init[4] = p_true[4] * (1 + delta * nodes_[4, i])  # <- Ri
    #
    #             p_init[0] = p_true[0] + nodes_[0, i]  # <- A
    #             p_init[1] = p_true[1] + nodes_[1, i]  # <- B
    #             p_init[2] = p_true[2] + nodes_[2, i]  # <- a
    #             p_init[3] = p_true[3] + nodes_[3, i]  # <- b
    #             p_init[4] = p_true[4] + nodes_[4, i]  # <- Ri
    #             # p_init[5] = p_true[5] + nodes_[5, i]  # <- L
    #             # p_init[6] = p_true[6] + nodes_[6, i]  # <- Req
    #             pass
    #
    #         par_mid = list(np.append(p_init, shape['IC'][5:]))
    #
    #         # perform checks on geometry
    #         ok = perform_geometry_checks(par_mid, par_end)
    #         if not ok:
    #             err = True
    #             break
    #         fid = fr'{key}_Q{i}'
    #
    #         # skip analysis if folder already exists.
    #         if not skip:
    #             if select_solver.lower() == 'slans':
    #                 solver = slans_geom
    #             else:
    #                 print(' ngsolve selected')
    #                 solver = ngsolve_mevp
    #             #  run model using SLANS or CST
    #             # # create folders for all keys
    #             solver.createFolder(fid, projectDir, subdir=sub_dir)
    #
    #             if "CELL TYPE" in shape.keys():
    #                 if shape['CELL TYPE'] == 'flattop':
    #                     # write_cst_paramters(fid, shape['IC'], shape['OC'], shape['OC_R'],
    #                     #                     projectDir=projectDir, cell_type="None", solver=select_solver.lower())
    #                     try:
    #                         print(' in flattop')
    #                         solver.cavity_flattop(n_cells, n_modules, par_mid, par_end, par_end,
    #                                               n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
    #                                               beampipes=shape['BP'],
    #                                               parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
    #                                               mesh_args=mesh_args)
    #                     except KeyError:
    #                         solver.cavity_flattop(n_cells, n_modules, par_mid, par_end, par_end,
    #                                               n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
    #                                               beampipes=shape['BP'],
    #                                               parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
    #                                               mesh_args=mesh_args)
    #             else:
    #                 try:
    #                     solver.cavity(n_cells, n_modules, par_mid, par_end, par_end,
    #                                   n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
    #                                   parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)
    #                 except KeyError:
    #                     solver.cavity(n_cells, n_modules, par_mid, par_end, par_end,
    #                                   n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol, beampipes=shape['BP'],
    #                                   parentDir=parentDir, projectDir=projectDir, subdir=sub_dir, mesh_args=mesh_args)
    #
    #         filename = uq_path / f'{fid}/monopole/qois.json'
    #         print(filename)
    #         if os.path.exists(filename):
    #             # params = fr.svl_reader(filename)
    #             # norm_length = 2 * n_cells * shape['IC'][5]
    #
    #             qois_result_dict = dict()
    #
    #             with open(filename) as json_file:
    #                 qois_result_dict.update(json.load(json_file))
    #
    #             qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
    #             # print_(qois_result)
    #             # sometimes some degenerate shapes are still generated and the solver returns zero
    #             # for the objective functions, such shapes are considered invalid
    #             for objr in qois_result:
    #                 if objr == 0:
    #                     # skip key
    #                     err = True
    #                     break
    #
    #             tab_val_f = qois_result
    #
    #             Ttab_val_f.append(tab_val_f)
    #         else:
    #             err = True
    #
    #         data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    #         data_table.to_csv(uq_path / 'table.csv', index=False, sep='\t', float_format='%.32f')
    #         data_table.to_excel(uq_path / 'table.xlsx', index=False)
    #     # # add original point
    #     # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    #     # params = fr.svl_reader(filename)
    #     # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    #     # tab_val_f = obj_result
    #     # Ttab_val_f.append(tab_val_f)
    #
    #     # import matplotlib.pyplot as plt
    #     print(np.atleast_2d(Ttab_val_f), weights_)
    #     if not err:
    #         v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights_)
    #
    #         # append results to dict
    #         for i, o in enumerate(eigen_obj_list):
    #             result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
    #             result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])
    #
    #             # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
    #             # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)
    #
    #         # plt.show()
    #         print(result_dict_eigen)
    #         with open(uq_path / fr"uq.json", 'w') as file:
    #             file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))
    #     else:
    #         error(fr"There was a problem running UQ analysis for {key}")

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

    @staticmethod
    def _run_abci(name, n_cells, n_modules, shape, MROT=0, MT=4.0, NFS=10000, UBT=50.0, bunch_length=20.0,
                  DDR_SIG=0.1, DDZ_SIG=0.1,
                  parentDir=None, projectDir=None,
                  WG_M=None, marker='', operating_points=None, freq=0, R_Q=0):

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
        if len(operating_points.keys()) > 0:
            try:
                if freq != 0 and R_Q != 0:
                    d = {}
                    # save qois
                    for key, vals in tqdm(operating_points.items()):
                        WP = key
                        I0 = float(vals['I0 [mA]'])
                        Nb = float(vals['Nb [1e11]'])
                        sigma_z = [float(vals["sigma_SR [mm]"]), float(vals["sigma_BS [mm]"])]
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

    # @staticmethod
    # def uq(shape_space, objectives, solver_dict, solver_args_dict):
    #     for key, shape in shape_space.items():
    #         err = False
    #         result_dict_eigen, result_dict_abci = {}, {}
    #         run_slans, run_abci = False, False
    #         eigen_obj_list, abci_obj_list = [], []
    #         for o in objectives:
    #
    #             if o[1] in ["Req", "freq", "Q", "E", "R/Q", "Epk/Eacc", "Bpk/Eacc"]:
    #                 result_dict_eigen[o[1]] = {'expe': [], 'stdDev': []}
    #                 run_slans = True
    #                 eigen_obj_list.append(o)
    #
    #             if o[1].split(' ')[0] in ['ZL', 'ZT', 'k_loss', 'k_kick']:
    #                 result_dict_abci[o[1]] = {'expe': [], 'stdDev': []}
    #                 run_abci = True
    #                 abci_obj_list.append(o)
    #
    #         # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    #         p_true = shape['IC'][0:5]
    #         # print(p_true)
    #         rdim = len(p_true)  # How many variabels will be considered as random in our case 5
    #         degree = 1
    #
    #         #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
    #         nodes = np.array(0)  # initialization
    #         weights = np.array(0)  # initialization
    #         flag_stroud = 1
    #         if flag_stroud == 1:
    #             nodes, weights, bpoly = quad_stroud3(rdim, degree)
    #             nodes = 2. * nodes - 1.
    #         elif flag_stroud == 2:
    #             nodes, weights, bpoly = quad_stroud3(rdim, degree)  # change to stroud 5 later
    #             nodes = 2. * nodes - 1.
    #         else:
    #             print('flag_stroud==1 or flag_stroud==2')
    #
    #         #  mean value of geometrical parameters
    #         p_init = np.zeros(np.shape(p_true))
    #
    #         no_parm, no_sims = np.shape(nodes)
    #         # print(no_sims)
    #         delta = 0.05  # or 0.1
    #
    #         if run_abci:
    #             # print("here in ANCI UQ")
    #             Ttab_val_f = []
    #             solver, solver_args = solver_dict['abci'], solver_args_dict['abci']
    #             n_cells = solver_args['n_cells']
    #             n_modules = solver_args['n_modules']
    #             MROT = solver_args['MROT']
    #             MT = solver_args['MT']
    #             NFS = solver_args['NFS']
    #             UBT = solver_args['UBT']
    #             bunch_length = solver_args['bunch_length']
    #             DDR_SIG = solver_args['DDR_SIG']
    #             DDZ_SIG = solver_args['DDZ_SIG']
    #             parentDir = solver_args['parentDir']
    #             projectDir = solver_args['projectDir']
    #             progress_list = solver_args['progress_list']
    #             WG_M = solver_args['WG_M']
    #             marker = solver_args['marker']
    #
    #             proc = solver_args['proc']
    #             sub_dir = fr'{key}'  # the simulation runs at the quadrature points
    #             # are saved to the key of the mean value run
    #             no_error = True
    #             for i in range(no_sims):
    #                 skip = False
    #                 p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
    #                 p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
    #                 p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
    #                 p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
    #                 p_init[4] = p_true[4] * (1 + delta * nodes[4, i])
    #
    #                 par_mid = np.append(p_init, shape['IC'][5:]).tolist()
    #                 par_end = par_mid
    #
    #                 ok = perform_geometry_checks(par_mid, par_end)
    #                 if not ok:
    #                     no_error = False
    #                     break
    #
    #                 fid = fr'{key}_Q{i}'
    #
    #                 # check if folder exists and skip if it does
    #                 if os.path.exists(fr'{projectDir}\SimulationData\ABCI\{key}\{fid}'):
    #                     skip = True
    #
    #                 if not skip:
    #                     #  run your model using SLANC or CST
    #                     # # create folders for all keys
    #                     solver.createFolder(fid, projectDir, subdir=sub_dir)
    #                     for wi in range(MROT):
    #                         solver.cavity(n_cells, n_modules, par_mid, par_end, par_end, fid=fid, MROT=wi,
    #                                       DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, beampipes=None, bunch_length=bunch_length,
    #                                       MT=MT, NFS=NFS, UBT=UBT,
    #                                       parentDir=parentDir, projectDir=projectDir, WG_M='',
    #                                       marker='', sub_dir=sub_dir
    #                                       )
    #
    #                 # get objective function values
    #                 abci_folder = fr'{projectDir}\SimulationData\ABCI\{key}'
    #                 if os.path.exists(abci_folder):
    #                     # print(abci_obj_list)
    #                     obj_result = get_wakefield_objectives_value(fid, abci_obj_list, abci_folder)
    #                     # print(obj_result)
    #
    #                     tab_val_f = obj_result
    #                     if 'error' in obj_result:
    #                         no_error = False
    #                         print(obj_result)
    #                         print("Encountered an error")
    #                         break
    #                     Ttab_val_f.append(tab_val_f)
    #                 else:
    #                     no_error = False
    #
    #             if no_error:
    #                 v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
    #                 # append results to dict
    #                 # print(v_expe_fobj, v_stdDev_fobj)
    #                 for i, o in enumerate(abci_obj_list):
    #                     result_dict_abci[o[1]]['expe'].append(v_expe_fobj[i])
    #                     result_dict_abci[o[1]]['stdDev'].append(v_stdDev_fobj[i])
    #
    #                 with open(fr"{projectDir}\SimulationData\ABCI\{key}\uq.json", 'w') as file:
    #                     file.write(json.dumps(result_dict_abci, indent=4, separators=(',', ': ')))

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
        if os.path.exists(fr"{self.folder}\SimulationData\SLANS_Opt\{self.name}\{tune_res}"):
            with open(fr"{self.folder}\SimulationData\SLANS_Opt\{self.name}\{tune_res}", 'r') as json_file:
                self.eigenmode_tune_res = json.load(json_file)

            self.freq = self.eigenmode_tune_res['freq']

            if tune_variable.lower() == 'req':
                self.shape_space['IC'][6] = self.eigenmode_tune_res['Req']
                self.shape_space['OC'][6] = self.eigenmode_tune_res['Req']
                self.shape_space['OC_R'][6] = self.eigenmode_tune_res['Req']
                self.mid_cell[6] = self.eigenmode_tune_res['Req']
                self.end_cell_left[6] = self.eigenmode_tune_res['Req']
                self.end_cell_right[6] = self.eigenmode_tune_res['Req']
            else:
                self.shape_space['IC'][5] = self.eigenmode_tune_res['L']
                self.shape_space['OC'][5] = self.eigenmode_tune_res['L']
                self.shape_space['OC_R'][5] = self.eigenmode_tune_res['L']
                self.mid_cell[5] = self.eigenmode_tune_res['L']
                self.end_cell_left[5] = self.eigenmode_tune_res['L']
                self.end_cell_right[5] = self.eigenmode_tune_res['L']

            # set alpha
            if len(self.mid_cell) == 7:
                if isinstance(self.shape_space['IC'], np.ndarray):
                    self.shape_space['IC'] = np.append(self.shape_space['IC'], self.eigenmode_tune_res['alpha_i'])
                    self.mid_cell = np.append(self.mid_cell, self.eigenmode_tune_res['alpha_i'])
                elif isinstance(self.shape_space['IC'], list):
                    self.shape_space['IC'] = np.append(self.shape_space['IC'], self.eigenmode_tune_res['alpha_i'])
                    self.mid_cell = np.append(self.mid_cell, self.eigenmode_tune_res['alpha_i'])
            elif len(self.mid_cell) == 8:
                self.shape_space['IC'][7] = self.eigenmode_tune_res['alpha_i']
                self.mid_cell[7] = self.eigenmode_tune_res['alpha_i']

            if cell_type.lower() == 'end-mid cell' \
                    or cell_type.lower() == 'mid-end cell' \
                    or cell_type.lower() == 'single cell':

                if len(self.end_cell_left) == 7:
                    if isinstance(self.shape_space['OC'], np.ndarray):
                        self.shape_space['OC'] = np.append(self.shape_space['OC'], self.eigenmode_tune_res['alpha_o'])
                        self.end_cell_left = np.append(self.end_cell_left, self.eigenmode_tune_res['alpha_o'])
                    elif isinstance(self.shape_space['OC'], list):
                        self.shape_space['OC'].append(self.eigenmode_tune_res['alpha_o'])
                        self.end_cell_left = np.append(self.end_cell_left, self.eigenmode_tune_res['alpha_o'])

                elif len(self.end_cell_left) == 8:
                    self.shape_space['OC'][7] = self.eigenmode_tune_res['alpha_o']
                    self.end_cell_left[7] = self.eigenmode_tune_res['alpha_o']

            # expand tune to be able to tune right cavity geometry also
        else:
            error("Tune results not found. Please tune the cavity")

    def get_ngsolve_tune_res(self, tune_variable, cell_type):
        """

        Parameters
        ----------
        tune_variable: {'A', 'B', 'a'. 'b', 'Ri', 'L', 'Req'}
            Tune variable.
        cell_type: {'mid cell', 'end-mid cell', 'mid-end cell', 'single cell'}
            Type of cell to tune

        Returns
        -------

        """
        tune_res = 'tune_res.json'
        print(self.name)
        if os.path.exists(fr"{self.folder}\SimulationData\Optimisation\{self.name}\{tune_res}"):
            with open(fr"{self.folder}\SimulationData\Optimisation\{self.name}\{tune_res}", 'r') as json_file:
                self.eigenmode_tune_res = json.load(json_file)

            self.freq = self.eigenmode_tune_res['FREQ']
            self.shape_space['IC'] = self.eigenmode_tune_res['IC']
            self.shape_space['OC'] = self.eigenmode_tune_res['OC']
            self.shape_space['OC_R'] = self.eigenmode_tune_res['OC_R']

            # expand tune to be able to tune right cavity geometry also
        else:
            error("Tune results not found. Please tune the cavity")

    def get_eigenmode_qois(self, solver_folder):
        """
        Get quantities of interest written by the SLANS code
        Returns
        -------

        """
        qois = 'qois.json'
        with open(fr"{self.folder}/SimulationData/{solver_folder}/{self.name}/monopole/{qois}") as json_file:
            self.eigenmode_qois = json.load(json_file)

        self.freq = self.eigenmode_qois['freq [MHz]']
        self.k_cc = self.eigenmode_qois['kcc [%]']
        self.ff = self.eigenmode_qois['ff [%]']
        self.R_Q = self.eigenmode_qois['R/Q [Ohm]']
        self.GR_Q = self.eigenmode_qois['GR/Q [Ohm^2]']
        self.G = self.GR_Q / self.R_Q
        # self.Q = d_qois['Q []']
        self.e = self.eigenmode_qois['Epk/Eacc []']
        self.b = self.eigenmode_qois['Bpk/Eacc [mT/MV/m]']

        # # get axis field
        # self.axis_field = fr.txt_reader(fr"{self.folder}\SimulationData\SLANS\{self.name}\cavity_33_{self.n_cells}.af",
        #                                 ' ')
        # # get surface field
        # self.surface_field = fr.txt_reader(
        #     fr"{self.folder}\SimulationData\SLANS\{self.name}\cavity_33_{self.n_cells}.sf", ' ')
        # # print(self.surface_field)

    def get_uq_fm_results(self, folder):
        # load uq result
        with open(folder, 'r') as json_file:
            self.uq_fm_results = json.load(json_file)

    # def get_custom_eig_qois(self):
    #     """
    #     Get quantities of interest written by the native eigenmode analysis code
    #     Returns
    #     -------
    #
    #     """
    #     qois = 'qois.json'
    #
    #     with open(fr"{self.folder}\SimulationData\NativeEig\{self.name}\{qois}") as json_file:
    #         self.custom_eig_qois = json.load(json_file)

    def get_abci_qois(self):
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

        if os.path.exists(fr"{self.folder}\SimulationData\ABCI\{self.name}\{qois}"):
            with open(fr"{self.folder}\SimulationData\ABCI\{self.name}\{qois}") as json_file:
                self.abci_qois = json.load(json_file)

        for key, val in self.abci_qois.items():
            self.k_fm[key] = val['k_FM [V/pC]']
            self.k_loss[key] = val['|k_loss| [V/pC]']
            self.k_kick[key] = val['|k_kick| [V/pC/m]']
            self.phom[key] = val['P_HOM [kW]']
            self.I0[key] = val['I0 [mA]']

    def get_abci_data(self):
        abci_data_dir = os.path.join(self.folder, "SimulationData", "ABCI")
        self.abci_data = {'Long': ABCIData(abci_data_dir, self.name, 0),
                          'Trans': ABCIData(abci_data_dir, self.name, 1)}

    def plot(self, what, ax=None, **kwargs):
        if what.lower() == 'geometry':
            if 'mid_cell' in kwargs.keys():
                new_kwargs = {key: val for key, val in kwargs.items() if key != 'mid_cell'}
                ax = plot_cavity_geometry_cli(self.mid_cell, self.mid_cell, self.mid_cell,
                                              'none', 1, scale=1, ax=ax, **new_kwargs)
            elif 'end_cell_left' in kwargs.keys():
                ax = plot_cavity_geometry_cli(self.end_cell_left, self.end_cell_left, self.end_cell_left,
                                              'left', 1, scale=1, ax=ax, **kwargs)
            elif 'end_cell_right' in kwargs.keys():
                ax = plot_cavity_geometry_cli(self.end_cell_right, self.end_cell_right, self.end_cell_right,
                                              'right', 1, scale=1, ax=ax, **kwargs)
            else:
                ax = plot_cavity_geometry_cli(self.mid_cell, self.end_cell_left, self.end_cell_right,
                                              self.beampipe, self.n_cells, scale=1, ax=ax, **kwargs)
            ax.set_xlabel('$z$ [mm]')
            ax.set_ylabel(r"$r$ [mm]")
            return ax

        if what.lower() == 'zl':
            if ax:
                x, y, _ = self.abci_data['Long'].get_data('Longitudinal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            else:
                fig, ax = plt.subplots(figsize=(12, 4))
                ax.margins(x=0)
                x, y, _ = self.abci_data['Long'].get_data('Longitudinal Impedance Magnitude')
                ax.plot(x * 1e3, y)

            ax.set_xlabel('f [MHz]')
            ax.set_ylabel(r"$Z_{\parallel} ~[\mathrm{k\Omega}]$")
            return ax
        if what.lower() == 'zt':
            if ax:
                x, y, _ = self.abci_data['Trans'].get_data('Transversal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            else:
                fig, ax = plt.subplots(figsize=(12, 4))
                ax.margins(x=0)
                x, y, _ = self.abci_data['Trans'].get_data('Transversal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            ax.set_xlabel('f [MHz]')
            ax.set_ylabel(r"$Z_{\perp} ~[\mathrm{k\Omega/m}]$")
            return ax

        if what.lower() == 'convergence':
            try:
                if ax:
                    self._plot_convergence(ax)
                else:
                    fig, ax = plt.subplot_mosaic([['conv', 'abs_err']], layout='constrained', figsize=(12, 4))
                    self._plot_convergence(ax)
                return ax
            except ValueError:
                info("Convergence data not available.")

    def plot_mesh(self, plotter='matplotlib'):
        ngsolve_mevp.plot_mesh(fr'{self.folder}\SimulationData\NGSolveMEVP\{self.name}\monopole', plotter=plotter)

    def plot_fields(self, mode=1, plotter='matplotlib'):
        ngsolve_mevp.plot_fields(fr'{self.folder}\SimulationData\NGSolveMEVP\{self.name}\monopole', mode, plotter)

    def _plot_convergence(self, ax):
        keys = list(ax.keys())
        # plot convergence
        conv_filepath = fr"{self.folder}\SimulationData\Optimisation\{self.name}\convergence.json"
        if os.path.exists(conv_filepath):
            with open(conv_filepath, 'r') as f:
                convergence_dict = json.load(f)
            if len(convergence_dict) > 0:
                x, y = convergence_dict[list(convergence_dict.keys())[0]], convergence_dict['freq [MHz]']
                ax[keys[0]].scatter(x, y, ec='k')

                # plot directions
                for i in range(len(x) - 1):
                    dx = x[i + 1] - x[i]
                    dy = y[i + 1] - y[i]
                    ax[keys[0]].quiver(x[i], y[i], dx, dy, ls='--', angles='xy',
                                       scale_units='xy', scale=1, color='red',
                                       width=0.005, units='width', headwidth=3, headlength=5,
                                       headaxislength=4)

        # plot absolute error
        abs_err_filepath = fr"{self.folder}\SimulationData\Optimisation\{self.name}\absolute_error.json"
        abs_err_dict = {}
        if os.path.exists(conv_filepath):
            with open(abs_err_filepath, 'r') as f:
                abs_err_dict = json.load(f)
        if len(abs_err_dict) > 0:
            ax[keys[1]].plot(abs_err_dict['abs_err'], marker='o', mec='k')

        ax[keys[0]].set_xlabel('Parameter')
        ax[keys[0]].set_ylabel(r"Value")
        ax[keys[1]].set_xlabel('Iteration')
        ax[keys[1]].set_ylabel(r"Absolute error")
        ax[keys[1]].set_yscale('log')

    def define_operating_points(self, op):
        self.operating_points = op

    def inspect(self, variation=0.2):
        import ipywidgets as widgets
        from ipywidgets import interact, HBox, VBox, Label

        mid_cell = self.shape_space['IC']
        A_, B_, a_, b_, Ri_, L_, Req_ = mid_cell[:7]

        # Define the function that plots the graph
        def plot_cavity_geometry(A, B, a, b, Ri, L, Req):
            mid_cell = np.array([A, B, a, b, Ri, L, Req])
            plot_cavity_geometry_cli(mid_cell, mid_cell, mid_cell, BP='none', n_cell=1, tangent_check=True, lw=1,
                                     ignore_degenerate=True)

            # Update the sum display
            sum_label.value = f'Sum of A + a: {A + a:.2f}, L: {L}, delta: {A + a - L}'

        # Create sliders for each variable
        A_slider = widgets.FloatSlider(min=(1 - variation) * A_, max=(1 + variation) * A_, step=0.1, value=A_,
                                       description='A')
        B_slider = widgets.FloatSlider(min=(1 - variation) * B_, max=(1 + variation) * B_, step=0.1, value=B_,
                                       description='B')
        a_slider = widgets.FloatSlider(min=(1 - variation) * a_, max=(1 + variation) * a_, step=0.1, value=a_,
                                       description='a')
        b_slider = widgets.FloatSlider(min=(1 - variation) * b_, max=(1 + variation) * b_, step=0.1, value=b_,
                                       description='b')
        Ri_slider = widgets.FloatSlider(min=(1 - variation) * Ri_, max=(1 + variation) * Ri_, step=0.1, value=Ri_,
                                        description='Ri')
        L_slider = widgets.FloatSlider(min=(1 - variation) * L_, max=(1 + variation) * L_, step=0.1, value=L_,
                                       description='L')
        Req_slider = widgets.FloatSlider(min=(1 - variation) * Req_, max=(1 + variation) * Req_, step=0.1, value=Req_,
                                         description='Req')
        # Create a label to display the sum of A + a
        sum_label = Label()

        # Arrange the sliders in a 3x3 layout
        ui = VBox([
            HBox([A_slider, B_slider, a_slider]),
            HBox([b_slider, Ri_slider, L_slider]),
            HBox([Req_slider]),
            sum_label  # Add the sum label to the layout
        ])

        # Create an interactive widget to update the plot
        out = widgets.interactive_output(plot_cavity_geometry,
                                         {'A': A_slider, 'B': B_slider, 'a': a_slider, 'b': b_slider,
                                          'Ri': Ri_slider, 'L': L_slider, 'Req': Req_slider})

        # Display the layout
        display(out, ui)

    def to_multicell(self):
        mid_cell = self.shape_space['IC']
        mid_cell_multi = np.array([[[a, a] for _ in range(self.n_cells - 1)] for a in mid_cell])

        self.shape_space_multicell['OC'] = self.shape_space['OC']
        self.shape_space_multicell['OC_R'] = self.shape_space['OC_R']
        self.shape_space_multicell['IC'] = mid_cell_multi
        self.shape_space_multicell['BP'] = self.shape_space['BP']
        self.shape_space_multicell['CELL TYPE'] = 'multicell'

    def _create_project(self, overwrite):
        project_name = self.name
        project_dir = self.folder

        if project_name != '':

            # check if folder already exist
            e = self._check_if_path_exists(project_dir, project_name, overwrite)

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
                                'SLANS_Opt': None,
                                'NGSolveMEVP': None,
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
                    error("An exception occurred in created project: ", e)
                    return False
            else:
                self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
                return True
        else:
            info('\tPlease enter a valid project name')
            self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
            return False

    @staticmethod
    def _check_if_path_exists(directory, folder, overwrite=False):
        path = f"{directory}/{folder}"
        if os.path.exists(path):
            if overwrite:
                x = 'y'
            else:
                x = 'n'

            if x == 'y':
                try:
                    directory_list = os.listdir(path)

                    if 'Cavities' in directory_list \
                            and 'PostprocessingData' in directory_list \
                            and 'SimulationData' in directory_list and len(directory_list) < 6:
                        shutil.rmtree(path)
                        return True
                    else:
                        info('\tIt seems that the folder specified is not a cavity project folder. Please check folder'
                             'again to avoid deleting important files.')
                        return False

                except Exception as e:
                    error("Exception occurred: ", e)
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


class Cavities(Optimisation):
    """
    Cavities object is an object containing several Cavity objects.
    """

    def __init__(self, cavities_list=None, names_list=None, save_folder='None'):
        """Constructs all the necessary attributes of the Cavity object

        Parameters
        ----------
        cavities_list: list, array like
            List containing Cavity objects.

        save_folder: str
            Folder to save generated images, latex, excel, and text files.
        """

        self.cavities_list = cavities_list
        self.cavities_dict = {}
        if cavities_list is None or cavities_list == []:
            self.cavities_list = []
        else:
            self.add_cavity(cavities_list, names_list)

        self.name = 'cavities'
        self.eigenmode_qois = {}
        self.eigenmode_tune_res = {}
        self.abci_qois = {}
        self.uq_fm_results = {}

        self.p_qois = None
        self.fm_results = None
        self.hom_results = None
        self.folder = save_folder
        self.operating_points = None
        self.operating_points_threshold = {}

        self.returned_results = None
        self.ls = ['solid', 'dashed', 'dashdot', 'dotted',
                   'solid', 'dashed', 'dashdot', 'dotted',
                   'solid', 'dashed', 'dashdot', 'dotted']

        self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
        self.set_cavities_field()

    def add_cavity(self, cavs, names=None, plot_labels=None):
        """
        Adds cavity to cavities
        Parameters
        ----------
        plot_label
        cav: Cavity, list
            Cavity object or list of cavity objects
        names: list, str
            Cavity name or list of cavity names

        Returns
        -------

        """

        if isinstance(cavs, Cavity):
            cavs.folder = self.folder
            if names:
                cavs.set_name(names)
            else:
                cavs.set_name(f'cav_{len(self.cavities_list)}')

            if isinstance(plot_labels, list):
                cavs.set_plot_label(plot_labels[0])
            else:
                cavs.set_plot_label(plot_labels)

            self.cavities_list.append(cavs)
            self.cavities_dict[cavs.name] = cavs
        else:
            if names is not None:
                assert len(cavs) == len(names), "Number of cavities does not correspond to number of names."
            else:
                names = [f'cav_{ii}' for ii in range(len(self.cavities_list))]

            if plot_labels is not None:
                assert len(cavs) == len(plot_labels), "Number of cavities does not correspond to number of labels."
            else:
                plot_labels = names

            for i1, cav in enumerate(cavs):
                cav.folder = self.folder
                cav.set_name(names[i1])
                cav.set_plot_label(plot_labels[i1])

                self.cavities_list.append(cav)
                self.cavities_dict[cav.name] = cav

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

    def save(self, files_path, overwrite=False):
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
                success = self._create_project(overwrite)
                if not success:
                    error(f"Project {files_path} could not be created. Please check the folder and try again.")
                    return
                else:
                    done(f"Project {files_path} created successfully/already exists.")
            except Exception as e:
                error("Exception occurred: ", e)
                return

        if files_path is None:
            self.folder = Path(os.getcwd())

    def save_plot_as_json(self, ax):
        ax_children = ax.get_children()

        ax_objects_dict = {
            'ax props': None,
            'fig props': {},
            'lines': {
                # 'line1': line_properties,
            },
            'axvline': {},
            'axhline': {},
            'scatter': {},
            'patches': {
                'text': {},
                'rectangle': {}
            },
        }

        # get axis properties
        axis_properties = self._get_axis_properties(ax)
        ax_objects_dict['ax props'] = axis_properties

        # get figure properties
        fig_props = self._get_figure_properties(ax.get_figure())
        ax_objects_dict['fig props'] = fig_props

        # List of known properties for Line2D
        line_properties = [
            'xdata', 'ydata', 'alpha', 'animated', 'antialiased', 'clip_on', 'clip_path',
            'color', 'dash_capstyle', 'dash_joinstyle', 'drawstyle',
            'label', 'linestyle', 'linewidth', 'marker', 'markeredgecolor',
            'markeredgewidth', 'markerfacecolor', 'markersize', 'path_effects', 'rasterized', 'sketch_params',
            'snap',
            'solid_capstyle',
            'solid_joinstyle', 'url', 'visible', 'zorder',
            # 'transform',
            # 'clip_box',
            # 'figure', 'picker', 'pickradius'
        ]

        for mpl_obj in ax_children:
            if isinstance(mpl_obj, matplotlib.lines.Line2D):
                # Extract properties into a dictionary
                line_properties_data = {prop: getattr(mpl_obj, f'get_{prop}')() for prop in line_properties}
                ax_objects_dict['lines'][fr'{id(mpl_obj)}'] = line_properties_data

        return ax_objects_dict

    def load_plot_from_json(self, filepath, ax=None):
        with open(filepath, 'r') as f:
            plot_config = json.load(f)

        if ax is not None:
            fig, ax = plt.subplot_mosaic([[0]])

        for key, value in plot_config.items():
            if key == 'lines':
                line = ax.plot([], [])
                line.set()
            if key == 'scatter':
                ax.scatter(value['x'], value['y'], value['kwargs'])
            if key == 'patches':
                if value['type'] == 'text':
                    ax.add_text(value['x'], value['y'], value['kwargs'])

    def plot_from_json(self, plot_config, ax=None):
        if ax is None:
            fig, ax = plt.subplot_mosaic([[0]])

        for key, value in plot_config.items():
            if key == 'lines':
                for line_keys, line_values in value.items():
                    line, = ax[0].plot([], [])
                    line.set(**line_values)
            if key == 'scatter':
                pass
                # ax[0].scatter(value['x'], value['y'], value['kwargs'])
            if key == 'patches':
                pass
                # if value['type'] == 'text':
                #     ax[0].add_text(value['x'], value['y'], value['kwargs'])
            if key == 'fig props':
                ax[0].get_figure().set(**value)

            if key == 'ax props':
                ax[0].set(**value)

        if isinstance(ax, dict):
            ax[0].relim()
            ax[0].get_figure().canvas.draw()
            ax[0].get_figure().canvas.flush_events()
            ax[0].get_figure().canvas.draw_idle()
        else:
            ax.relim()
            ax.get_figure().canvas.draw()
            ax.get_figure().canvas.flush_events()
            ax.get_figure().canvas.draw_idle()

        # plot legend wthout duplicates
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

    def _get_axis_properties(self, ax):
        # List of common properties for an Axis
        axis_properties = [
            'label', 'ticks', 'scale', 'margin', 'bound', 'aspect'
        ]

        def get_axis_properties(axis, prefix):
            properties = {}
            for prop in axis_properties:
                try:
                    # Use get_* methods dynamically if available
                    getter = getattr(axis, f'get_{prop}', None)
                    if callable(getter):
                        if prop == 'label':
                            properties[fr'{prefix}{prop}'] = getter().get_text()
                        else:
                            properties[fr'{prefix}{prop}'] = getter()
                    else:
                        # Fallback to axis.get_property()
                        properties[fr'{prefix}{prop}'] = axis.get_property(prop)
                except Exception as e:
                    pass
            return properties

        axis_properties = get_axis_properties(ax.xaxis, 'x')
        axis_properties.update(get_axis_properties(ax.yaxis, 'y'))

        return axis_properties

    def _get_figure_properties(self, fig):
        # List of properties to retrieve
        figure_properties = [
            'figwidth', 'figheight', 'dpi', 'tight_layout', 'constrained_layout'
        ]

        def get_figure_properties(fig):
            properties = {}
            for prop in figure_properties:
                getter = getattr(fig, f'get_{prop}', None)
                if callable(getter):
                    properties[prop] = getter()

                else:
                    properties[prop] = 'N/A'
            return properties

        # Get properties of the figure
        fig_size = []
        fig_properties = get_figure_properties(fig)

        return fig_properties

    def _create_project(self, overwrite):
        project_name = self.name
        project_dir = self.folder

        if project_name != '':

            # check if folder already exist
            e = self._check_if_path_exists(project_dir, project_name, overwrite)

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
                                'SLANS_Opt': None,
                                'NGSolveMEVP': None,
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
                # self.folder = os.path.join(project_dir, project_name)
                self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
                return True
        else:
            print('\tPlease enter a valid project name')
            self.folder = f2b_slashes(fr"{project_dir}\{project_name}")
            return False

    @staticmethod
    def _check_if_path_exists(directory, folder, overwrite):
        path = f"{directory}/{folder}"
        if os.path.exists(path):
            x = 'n'
            if overwrite:
                x = 'y'

            if x == 'y':
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

    def sweep(self, sweep_config):
        self.sweep_results = {}
        for cav in self.cavities_list:
            cav.sweep(sweep_config)
            self.sweep_results[cav.name] = cav.sweep_results

    def check_uq_config(self, uq_config):
        uq_ok = {}
        for cav in self.cavities_list:
            res = cav.check_uq_config(uq_config)
            uq_ok[cav.name] = res
        info(uq_ok)

    def run_tune(self, tune_variables, cell_types='Mid Cell', freqs=None, solver='SLANS', proc=0, resume=False,
                 n_cells=1, rerun=True):
        """
        Tune current cavity geometries

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
        freqs: float, list, ndarray
            Reference frequency or list of reference frequencies if different for each cavity in MHz
        cell_types: {'mid cell', 'end-mid cell', 'mid-end cell', 'single cell'}
            Type of cell to tune or list of type of cell to tune for the different cavities
        tune_variables: {'Req', 'L'}
            Tune variable or list of tune variables. Currently supports only the tuning of the equator radius ``Req`` and half-cell length ``L``

        Returns
        -------

        """

        for i, cav in enumerate(tqdm(self.cavities_list)):
            if isinstance(freqs, float) or isinstance(freqs, int):
                freq = freqs
            else:
                freq = freqs[i]

            if isinstance(tune_variables, str):
                tune_var = tune_variables
                cell_type = cell_types
            else:
                tune_var = tune_variables[i]
                cell_type = cell_types[i]

            if os.path.exists(os.path.join(self.folder, "SimulationData", "Optimisation", cav.name)):
                if rerun:
                    cav.run_tune(tune_var, freq=freq, cell_type=cell_type, solver=solver, n_cells=n_cells)
                else:
                    # check if tune results exist
                    if os.path.exists(
                            os.path.join(self.folder, "SimulationData", "Optimisation", cav.name, "tune_res.json")):
                        if solver == 'SLANS':
                            cav.get_slans_tune_res(tune_var, cell_type)
                        else:
                            cav.get_ngsolve_tune_res(tune_var, cell_type)

                    else:
                        cav.run_tune(tune_var, freq=freq, cell_type=cell_type, solver=solver, n_cells=n_cells)
            else:
                cav.run_tune(tune_var, freq=freq, cell_type=cell_type, solver=solver, n_cells=n_cells)

            self.eigenmode_tune_res[cav.name] = cav.eigenmode_tune_res

    def run_eigenmode(self, solver='SLANS', freq_shifts=0, boundary_conds=None, subdir='',
                      uq_config=None, rerun=True, procs=1):

        # perform all necessary checks
        assert procs > 0, error('Number of proceses must be greater than zero.')
        if uq_config:
            assert len(uq_config['delta']) == len(uq_config['variables']), error("The number of deltas must "
                                                                                 "be equal to the number of "
                                                                                 "variables.")

        # split shape_space for different processes/ MPI share process by rank
        keys = list(self.cavities_dict.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / procs)

        jobs = []
        for p in range(procs):
            # try:
            if p < procs - 1:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            processor_shape_space = {key: self.cavities_dict[key] for key in proc_keys_list}
            service = mp.Process(target=self.run_eigenmode_s, args=(processor_shape_space, self.folder,
                                                                    solver, freq_shifts, boundary_conds, subdir,
                                                                    uq_config, rerun))

            service.start()
            jobs.append(service)

        for job in jobs:
            job.join()

        self.get_eigenmode_qois(uq_config)

    def get_eigenmode_qois(self, uq_config):
        # get results
        for key, cav in self.cavities_dict.items():
            try:
                cav.get_eigenmode_qois('NGSolveMEVP')
                self.eigenmode_qois[cav.name] = cav.eigenmode_qois
                if uq_config:
                    cav.get_uq_fm_results(fr"{self.folder}\SimulationData\NGSolveMEVP\{cav.name}\uq.json")
                    self.uq_fm_results[cav.name] = cav.uq_fm_results
            except FileNotFoundError:
                error("Could not find the eigenmode results. Please rerun eigenmode analysis.")
                return False

    @staticmethod
    def run_eigenmode_s(shape_space, folder, solver='SLANS', freq_shifts=0, boundary_conds=None, subdir='',
                        uq_config=None, rerun=True):
        """
        Run eigenmode analysis on cavity

        Parameters
        ----------
        solver: {'SLANS', 'NGSolve'}
            Solver to be used. Native solver is still under development. Results are not as accurate as that of SLANS.
        freq_shifts:
            (List of) frequency shift. Eigenmode solver searches for eigenfrequencies around this value
        boundary_conds: int, list
            (List of) boundary condition of left and right cell/beampipe ends
        subdir: str
            Sub directory to save results to
        uq_config: None | dict
            Provides inputs required for uncertainty quantification. Default is None and disables uncertainty quantification.

        Returns
        -------

        """

        for i, cav in enumerate(tqdm(list(shape_space.values()))):
            if isinstance(freq_shifts, int) or isinstance(freq_shifts, float):
                freq_shift = freq_shifts
            else:
                freq_shift = freq_shifts[i]

            if isinstance(boundary_conds, str) or boundary_conds is None:
                boundary_cond = boundary_conds
            else:
                boundary_cond = boundary_conds[i]

            if solver.lower() == 'slans':
                solver_save_dir = 'SLANS'
            else:
                solver_save_dir = 'NGSolveMEVP'

            if os.path.exists(os.path.join(folder, "SimulationData", solver_save_dir, cav.name)):
                if rerun:
                    # delete old results
                    shutil.rmtree(os.path.join(folder, "SimulationData", solver_save_dir, cav.name))
                    cav.run_eigenmode(solver, freq_shift=freq_shift, boundary_cond=boundary_cond, uq_config=uq_config)
                    print(id(cav), cav.eigenmode_qois)
                else:
                    # check if eigenmode analysis results exist
                    if os.path.exists(os.path.join(folder, "SimulationData", solver_save_dir, cav.name, "monopole",
                                                   "qois.json")):
                        cav.get_eigenmode_qois(solver_save_dir)
                        cav.get_uq_fm_results(solver_save_dir)
                    else:
                        shutil.rmtree(os.path.join(folder, "SimulationData", solver_save_dir, cav.name))
                        cav.run_eigenmode(solver, freq_shift=freq_shift, boundary_cond=boundary_cond,
                                          uq_config=uq_config)
            else:
                cav.run_eigenmode(solver, freq_shift=freq_shift, boundary_cond=boundary_cond, uq_config=uq_config)

    def run_wakefield(self, MROT=2, MT=10, NFS=10000, wakelength=50, bunch_length=25,
                      DDR_SIG=0.1, DDZ_SIG=0.1, WG_M=None, marker='', operating_points=None,
                      solver='ABCI', rerun=True):
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

        for i, cav in enumerate(tqdm(self.cavities_list)):
            if os.path.exists(os.path.join(self.folder, "SimulationData", "ABCI", cav.name)):
                if rerun:
                    cav.run_wakefield(MROT, MT, NFS, wakelength, bunch_length,
                                      DDR_SIG, DDZ_SIG, WG_M, marker, operating_points, solver)
                else:
                    # check if eigenmode analysis results exist
                    cav.get_abci_data()
                    if os.path.exists(os.path.join(self.folder, "SimulationData", "ABCI", cav.name, "qois.json")):
                        cav.get_abci_qois()
                    else:
                        cav.run_wakefield(MROT, MT, NFS, wakelength, bunch_length,
                                          DDR_SIG, DDZ_SIG, WG_M, marker, operating_points, solver)
            else:
                cav.run_wakefield(MROT, MT, NFS, wakelength, bunch_length,
                                  DDR_SIG, DDZ_SIG, WG_M, marker, operating_points, solver)

            self.abci_qois[cav.name] = cav.abci_qois

    def plot(self, what, ax=None, **kwargs):
        for cav in self.cavities_list:
            if what.lower() == 'geometry':
                ax = cav.plot('geometry', ax, **kwargs)

            if what.lower() == 'zl':
                if ax:
                    ax = cav.plot('zl', ax)
                else:
                    fig, ax = plt.subplots(figsize=(12, 4))
                    ax.margins(x=0)
                    ax = cav.plot('zl', ax)

            if what.lower() == 'zt':
                if ax:
                    ax = cav.plot('zt', ax)
                else:
                    fig, ax = plt.subplots(figsize=(12, 4))
                    ax.margins(x=0)
                    ax = cav.plot('zt', ax)

            if what.lower() == 'convergence':
                ax = cav.plot('convergence', ax)
        return ax

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

        print(qois)
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
                r"$k_\mathrm{cc} [\cdot]$": cav.k_cc,
                r"$R/Q \mathrm{[\Omega]}$": cav.R_Q,
                r"$G \mathrm{[\Omega]}$": cav.G,
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

        return results

    def qois_hom(self, opt):
        """
        Retrieves the higher-order modes quantities of interest

        Returns
        -------
        Dictionary containing higher-order modes quantities of interest (normed optional).
        """

        results = []
        for cavity in self.cavities_list:
            cavity.get_abci_qois()
            results.append({
                r"$|k_\parallel| \mathrm{[V/pC]}$": cavity.k_loss[opt],
                r"$|k_\perp| \mathrm{[V/pC/m]}$": cavity.k_kick[opt],
                r"$P_\mathrm{HOM}/cav \mathrm{[kW]}$": cavity.phom[opt]
            })

        results_norm_units = []
        for cavity in self.cavities_list:
            cavity.get_abci_qois()
            results_norm_units.append({
                r"$k_\parallel$": cavity.k_loss[opt],
                r"$k_\perp$": cavity.k_kick[opt],
                r"$p_\mathrm{HOM}/cav$": cavity.phom[opt]
            })

        return results_norm_units

    def qois_all(self, opt):
        """
        Retrieves the fundamental mode quantities of interest

        Returns
        -------
        Dictionary containing fundamental mode quantities of interest (normed optional).
        """
        results = []
        for cav in self.cavities_list:
            results.append({
                r"$E_\mathrm{pk}/E_\mathrm{acc} ~[\cdot]$": cav.e,
                r"$B_\mathrm{pk}/E_\mathrm{acc} ~\mathrm{[mT/MV/m]}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$R/Q~ \mathrm{[\Omega]}$": cav.R_Q,
                r"$G ~\mathrm{[\Omega]}$": cav.G,
                r"$G\cdot R/Q ~\mathrm{[10^4\Omega^2]}$": cav.GR_Q * 1e-4,
                r"$|k_\parallel| ~\mathrm{[V/pC]}$": cav.k_loss[opt],
                r"$|k_\perp| ~\mathrm{[V/pC/m]}$": cav.k_kick[opt],
                r"$P_\mathrm{HOM}/cav~ \mathrm{[kW]}$": cav.phom[opt]
            })

        results_norm_units = []
        for cav in self.cavities_list:
            results_norm_units.append({
                r"$e_\mathrm{pk}/e_\mathrm{acc}$": cav.e,
                r"$b_\mathrm{pk}/e_\mathrm{acc}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$r/q$": cav.R_Q,
                r"$g$": cav.G,
                # r"$g\cdot r/q $": cav.GR_Q,
                # r"$|k_\mathrm{FM}|$": cav.k_fm,
                r"$|k_\parallel|$": cav.k_loss[opt],
                r"$k_\perp$": cav.k_kick[opt],
                r"$p_\mathrm{HOM}/cav$": cav.phom[opt]
            })

        return results

    def plot_uq_geometries(self):
        fig, axd = plt.subplot_mosaic([[cav.name for cav in self.cavities_list]], layout='constrained')

        for cav, ax in zip(self.cavities_list, axd.values()):
            # plot nominal
            cav.plot('geometry', ax=ax, mid_cell=True, zorder=10)
            directory = f'{self.folder}/SimulationData/NGSolveMEVP/{cav.name}'
            tag = f'{cav.name}_Q'
            uq_geom_folders = self.find_folders_with_tag(directory, tag)

            # load uq geometries
            for uq_geom_folder in uq_geom_folders:
                if os.path.exists(f'{uq_geom_folder}/monopole/geodata.n'):
                    # read geometry
                    cav_geom = pd.read_csv(f'{uq_geom_folder}/monopole/geodata.n', header=None,
                                           skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0, 2]]

                    cav_geom = cav_geom[[1, 0]]
                    ax.plot(cav_geom[1], cav_geom[0], ls='--', lw=1, c='gray')
            ax.set_title(cav.name)

    @staticmethod
    def find_folders_with_tag(directory, tag):
        matching_folders = []

        # Walk through the directory
        for root, dirs, files in os.walk(directory):
            for dir_name in dirs:
                # Check if the folder name contains the tag
                if fnmatch.fnmatch(dir_name, f'*{tag}*'):
                    matching_folders.append(os.path.join(root, dir_name))

        return matching_folders

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
        ax.margins(x=0)
        width = 0.15  # 1 / len(x)
        for i, cav in enumerate(self.cavities_list):
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

    def plot_compare_hom_bar(self, opt, ncols=3):
        """
        Plot bar chart of higher-order mode's quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (15 / 27 * 6 * len(self.cavities_list), 3)
        # plot barchart
        self.hom_results = self.qois_hom(opt)
        df = pd.DataFrame.from_dict(self.hom_results)
        fig, axd = plt.subplot_mosaic([list(df.columns)])
        # Plot each column in a separate subplot
        labels = [cav.plot_label for cav in self.cavities_list]
        for key, ax in axd.items():
            ax.bar(df.index, df[key], label=labels, color=matplotlib.colormaps['Set2'].colors[:len(df)],
                   edgecolor='k', width=1)
            ax.set_xticklabels([])
            ax.set_ylabel(key)
            h, l = ax.get_legend_handles_labels()
        if not ncols:
            ncols = min(4, len(self.cavities_list))
        fig.legend(h, l, loc='outside upper center', borderaxespad=0, ncol=ncols)

        # fig.set_tight_layout(True)

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_hom_bar.png")

        return axd

    def plot_compare_fm_bar(self, ncols=3, uq=False):
        """
        Plot bar chart of fundamental mode quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 4)

        if not uq:
            self.fm_results = self.qois_fm()

            df = pd.DataFrame.from_dict(self.fm_results)
            fig, axd = plt.subplot_mosaic([list(df.columns)], layout='constrained')

            labels = [cav.plot_label for cav in self.cavities_list]
            # Plot each column in a separate subplot
            for key, ax in axd.items():
                ax.bar(df.index, df[key], label=labels, color=matplotlib.colormaps['Set2'].colors[:len(df)],
                       edgecolor='k',
                       width=1)
                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.set_ylabel(key)
                h, l = ax.get_legend_handles_labels()
        else:
            self.fm_results = self.uq_fm_results

            # Step 1: Flatten the dictionary into a DataFrame
            rows = []
            for cav, metrics in self.uq_fm_results.items():
                for metric, values in metrics.items():
                    rows.append({
                        'cavity': cav,
                        'metric': metric,
                        'mean': values['expe'][0],
                        'std': values['stdDev'][0]
                    })

            df = pd.DataFrame(rows)

            labels = [cav.plot_label for cav in self.cavities_list]

            # Step 2: Create a Mosaic Plot
            metrics = df['metric'].unique()
            num_metrics = len(metrics)

            layout = [[metric for metric in metrics]]
            fig, axd = plt.subplot_mosaic(layout, layout='constrained')
            # fig, axd = plt.subplot_mosaic([list(df.columns)], layout='constrained')

            # Plot each metric on a separate subplot
            for metric, ax in axd.items():
                sub_df = df[df['metric'] == metric]
                ax.bar(sub_df['cavity'], sub_df['mean'], yerr=sub_df['std'], label=labels, capsize=5,
                       color=matplotlib.colormaps['Set2'].colors[:len(df)], edgecolor='k',
                       width=1)

                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.set_ylabel(metric)
                h, l = ax.get_legend_handles_labels()

        # fig.set_tight_layout(True)
        if not ncols:
            ncols = min(4, len(self.cavities_list))
        fig.legend(h, l, loc='outside upper center', borderaxespad=0, ncol=ncols)

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_all_bar.png")

        return axd

    def plot_compare_all_bar(self, opt, ncols=3):
        """
        Plot bar chart of fundamental mode quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (15, 3)
        # plot barchart
        self.all_results = self.qois_all(opt)

        df = pd.DataFrame.from_dict(self.all_results)
        fig, axd = plt.subplot_mosaic([list(df.columns)], layout='constrained')

        labels = [cav.plot_label for cav in self.cavities_list]
        # Plot each column in a separate subplot
        for key, ax in axd.items():
            ax.bar(df.index, df[key], label=labels, color=matplotlib.colormaps['Set2'].colors[:len(df)],
                   edgecolor='k', width=1)
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_ylabel(key)
            h, l = ax.get_legend_handles_labels()

        # fig.set_tight_layout(True)
        if not ncols:
            ncols = min(4, len(self.cavities_list))
        fig.legend(h, l, loc='outside upper center', borderaxespad=0, ncol=ncols)

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_all_bar.png")

        return axd

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

            axs[0].plot(n_cav_per_cryomodule, cryomodules_len_list, marker='o', mec='k', label=f'{cav.name}')
            axs[1].plot(n_cav_per_cryomodule, n_cryomodules_list, marker='o', mec='k', label=f'{cav.name}')
            print(n_cav_per_cryomodule)
            print(cryomodules_len_list, n_cryomodules_list)
            print(n_cav_per_cryomodule)
            print()
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

    def plot_cavities_contour(self, opt='mid'):
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

        fig, ax = plt.subplots()
        ax.margins(x=0)

        for i, cav in enumerate(self.cavities_list):
            # write contour
            # self.write_contour(cav, opt)

            if opt.lower() == 'mid':
                mid_cell = np.array(cav.shape_space['IC']) * 1e-3
                end_cell_left = np.array(cav.shape_space['IC']) * 1e-3
                end_cell_right = np.array(cav.shape_space['IC']) * 1e-3
            elif opt.lower() == 'end':
                mid_cell = np.array(cav.shape_space['IC']) * 1e-3
                end_cell_left = np.array(cav.shape_space['OC']) * 1e-3
                end_cell_right = np.array(cav.shape_space['OC']) * 1e-3
            else:
                mid_cell = np.array(cav.shape_space['IC']) * 1e-3
                end_cell_left = np.array(cav.shape_space['OC']) * 1e-3
                end_cell_right = np.array(cav.shape_space['OC']) * 1e-3

            scale = cav.op_freq / c0
            if cav.cell_type == 'flat top':
                writeCavityForMultipac_flat_top(fr'{cav.slans_dir}\contour.txt', 1, mid_cell, end_cell_left,
                                                end_cell_right, beampipe='none', unit=1, scale=scale, plot=True)
            else:
                writeCavityForMultipac(fr'{cav.slans_dir}\contour.txt', 1, mid_cell, end_cell_left, end_cell_right,
                                       beampipe='none', unit=1, scale=scale)

            data = pd.read_csv(fr"{cav.slans_dir}\contour.txt", sep=r'\s+', header=None, skiprows=3)

            # ax.plot(data[1] * 1000, data[0] * 1000, lw=3., label=cav.plot_label)
            ax.plot(data[1], data[0], lw=3., label=cav.plot_label, marker='x')
            ax.legend(loc='lower left')

            x_label = r"$z/\lambda$"
            y_label = r"$r/\lambda$"
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            min_x.append(min(data[1]))
            min_y.append(min(data[0]))
            max_x.append(max(data[1]))
            max_y.append(max(data[0]))

        if opt.lower() == 'mid' or opt.lower() == 'end':
            ax.set_xlim(0, max(max_x))
            ax.set_ylim(0, max(max_y))
        else:
            ax.set_xlim(min(min_x), max(max_x))
            ax.set_ylim(min(min_y), max(max_y))

        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_contour_{opt}.png")

        fig.show()

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
                    print('Unable to plot the counters. No initial points.')
                    return

                if ok1 * ok2 == 0:
                    cl = print('Counter functions or impact energy missing.')
                else:
                    # if ss > 0:
                    #     cl = print(np.array(['Plotting the triplot (counter, enhanced ', 'counter and impact energy).']))

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
        ax.margins(x=0)
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
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.shape_space['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.shape_space['IC']) * 1e-3
            A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.shape_space['IC']) * 1e-3
            n_cell = 1
            L_bp_l = 0.001
            L_bp_r = 0.001

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

        elif opt.lower() == 'end':
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.shape_space['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.shape_space['IC']) * 1e-3
            A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.shape_space['OC']) * 1e-3
            L_bp_l = 0.001
            L_bp_r = 1 * L_m

            n_cell = 1

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m) / 2
        else:
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.shape_space['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.shape_space['OC']) * 1e-3
            try:
                A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.shape_space['OC_R']) * 1e-3
            except KeyError:
                A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.shape_space['OC']) * 1e-3

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
                        print(pt, 1)
                        pts = self.arcTo(L_el + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                         [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req_er])
                        pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                        print(pt, 2)
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        print(pt, 3)

                        # STRAIGHT LINE TO NEXT POINT
                        self.lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                        pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        print(pt, 4, L_el + L_er - x1er + L_bp_l + L_bp_r - shift)

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        # print(shift)
                        pts = self.arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                         [L_bp_l + L_el + L_er - shift, y1er])
                        print(pt, 5, L_el + L_er + L_bp_l - shift)

                        pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        print(pt, 6)

                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_er)
                        # print(shift)
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
                        # print(shift)

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
                    print(pt)
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
        ax.margins(x=0)

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
            print(R_Q[i], v_cav[i], Vrf[i])
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
            print(l_active)

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

    def define_operating_points(self, ops):
        self.operating_points = ops

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
        # print([u + a * np.cos(t), v + b * np.sin(t)])
        # print()

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
            l4 = r"\resizebox{\textwidth}{!}{\begin{tabular}{|l|" + f"{''.join(['c' for i in self.cavities_list])}" + "|}"
            toprule = r"\toprule"
            header = r" ".join([fr"& {cav.name} " for cav in self.cavities_list]) + r" \\"
            hline = r"\hline"
            hhline = r"\hline \hline"
            A = r"$A$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][0], 2)}/{round(cav.d_geom_params['OC'][0], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            B = r"$B$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][1], 2)}/{round(cav.d_geom_params['OC'][1], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            a = r"$a$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][2], 2)}/{round(cav.d_geom_params['OC'][2], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            b = r"$b$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][3], 2)}/{round(cav.d_geom_params['OC'][3], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            Ri = r"$R_\mathrm{i}$ " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][4], 2)}/{round(cav.d_geom_params['OC'][4], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            L = r"$L$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][5], 2)}/{round(cav.d_geom_params['OC'][5], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            Req = r"$R_\mathrm{eq}$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][6], 2)}/{round(cav.d_geom_params['OC'][6], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            alpha = r"$ \alpha [^\circ]$" + "".join(
                [fr"& {round(cav.d_geom_params['IC'][7], 2)}/{round(cav.d_geom_params['OC'][7], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            rf_freq = r"RF Freq. [MHz] " + "".join(
                [fr"& {round(cav.op_freq * 1e-6, 2)} " for cav in self.cavities_list]) + r" \\"
            Vrf = r"$V_\mathrm{RF}$ [V] " + "".join(
                [fr"& {round(cav.v_rf, 2):.2E} " for cav in self.cavities_list]) + r" \\"
            Eacc = r"$E_\mathrm{acc}$ [MV/m] " + "".join(
                [fr"& {round(cav.op_field * 1e-6, 2)} " for cav in self.cavities_list]) + r" \\"
            R_Q = r"$R/Q [\Omega$] " + "".join([fr"& {round(cav.R_Q, 2)} " for cav in self.cavities_list]) + r" \\"
            G = r"$G$ [$\Omega$] " + "".join([fr"& {round(cav.G, 2)} " for cav in self.cavities_list]) + r" \\"
            GR_Q = r"$G.R/Q [10^4\Omega^2]$ " + "".join(
                [fr"& {round(cav.GR_Q * 1e-4, 2)} " for cav in self.cavities_list]) + r" \\"
            epk = r"$E_{\mathrm{pk}}/E_{\mathrm{acc}}$ " + "".join(
                [fr"& {round(cav.e, 2)} " for cav in self.cavities_list]) + r" \\"
            bpk = r"$B_{\mathrm{pk}}/E_{\mathrm{acc}} \left[\mathrm{\frac{mT}{MV/m}}\right]$ " + "".join(
                [fr"& {round(cav.b, 2)} " for cav in self.cavities_list]) + r" \\"

            kfm = r"$|k_\mathrm{FM}|$ [V/pC] " + "".join(
                [fr"& {'/'.join([str(round(k_fm, 4)) for k_fm in cav.k_fm])} " for cav in self.cavities_list]) + r" \\"
            kloss = r"$|k_\mathrm{\parallel}|$ [V/pC]" + "".join(
                [fr"& {'/'.join([str(round(k_loss, 4)) for k_loss in cav.k_loss])} " for cav in
                 self.cavities_list]) + r" \\"
            kkick = r"$k_\mathrm{\perp}$ [V/pC/m]" + "".join(
                [fr"& {'/'.join([str(round(k_kick, 4)) for k_kick in cav.k_kick])} " for cav in
                 self.cavities_list]) + r" \\"

            Ncav = r"$N_\mathrm{cav}$ " + "".join(
                [fr"& {int(cav.n_cav_op_field)} " for cav in self.cavities_list]) + r" \\"
            Q0 = r"$Q_\mathrm{0}$ " + "".join(
                [fr"& {int(cav.Q0_desired):2E} " for cav in self.cavities_list]) + r" \\"
            Pin = r"$P_\mathrm{in}\mathrm{/cav} [\mathrm{kW}]$ " + "".join(
                [fr"& {round(cav.p_in * 1e-3, 2)} " for cav in self.cavities_list]) + r" \\"

            Pstat = r"$P_\mathrm{stat}$/cav [W] " + "".join(
                [fr"& {round(cav.pstat, 2)} " for cav in self.cavities_list]) + r" \\"

            Pdyn = r"$P_\mathrm{dyn}$/cav [W] " + "".join(
                [fr"& {round(cav.pdyn, 2)} " for cav in self.cavities_list]) + r" \\"

            Pwp = r"$P_\mathrm{wp}$/cav [kW] " + "".join(
                [fr"& {round(cav.p_wp * 1e-3, 2)} " for cav in self.cavities_list]) + r" \\"

            Phom = r"$P_\mathrm{HOM}$/cav [kW] " + "".join(
                [fr"& {'/'.join([str(round(phom, 4)) for phom in cav.phom])} " for cav in self.cavities_list]) + r" \\"

            phom_values = np.array([[phom for phom in cav.phom] for cav in self.cavities_list])

            if len(phom_values) % 2 == 0:
                result_array = phom_values[1::2] - phom_values[::2]
                print(result_array)
                Phom_diff = r"$\Delta P_\mathrm{HOM}$/cav [W] & \multicolumn{2}{c|}{" + " & \multicolumn{2}{c|}{".join(
                    [fr"{'/'.join([str(round(val * 1e3, 2)) for val in arr])} " + '}' for arr in result_array]) + r" \\"

            PHOM = r"$P_\mathrm{HOM}$ [kW] " + "".join(
                [fr"& {'/'.join([str(round(phom * cav.n_cav_op_field, 2)) for phom in cav.phom])} " for cav in
                 self.cavities_list]) + r" \\"
            bottomrule = r"\bottomrule"
            l34 = r"\end{tabular}}"
            l35 = r"\label{tab: selected shape}"
            l36 = r"\end{table}"

            all_lines = (l1, l2, l3, l4,
                         toprule, hline,
                         header,
                         hhline,
                         A, B, a, b, Ri, L, Req, alpha,
                         hhline,
                         rf_freq, Vrf, Eacc, R_Q, G, GR_Q, epk, bpk,
                         hhline,
                         kfm, kloss, kkick, Phom,
                         hhline, Ncav, Q0, Pin, Pstat, Pdyn, Pwp, Phom, Phom_diff, PHOM,
                         hline,
                         bottomrule,
                         l34, l35, l36)

            # save plots
            fname = [cav.name for cav in self.cavities_list]
            fname = '_'.join(fname)
            print(fname)
            with open(fr"D:\Dropbox\Quick presentation files\{self.folder}\{fname}_latex_summary.txt", 'w') as f:
                for ll in all_lines:
                    f.write(ll + '\n')
        except KeyError as e:
            print("Either SLANS or ABCI results not available. Please use '<cav>.set_eigenmode_qois(<folder>)' "
                  "or '<cav>.set_abci_qois(<folder>)' to fix this. Error: ", e)

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
                    'A_i [mm]': [round(cav.shape_space['IC'][0], 2) for cav in self.cavities_list],
                    'B_i [mm]': [round(cav.shape_space['IC'][1], 2) for cav in self.cavities_list],
                    'a_i [mm]': [round(cav.shape_space['IC'][2], 2) for cav in self.cavities_list],
                    'b_i [mm]': [round(cav.shape_space['IC'][3], 2) for cav in self.cavities_list],
                    'R_i [mm]': [round(cav.shape_space['IC'][4], 2) for cav in self.cavities_list],
                    'L_i [mm]': [round(cav.shape_space['IC'][5], 2) for cav in self.cavities_list],
                    'Req [mm]': [round(cav.shape_space['IC'][6], 2) for cav in self.cavities_list],
                    'alpha_i [deg]': [round(cav.shape_space['IC'][7], 2) for cav in self.cavities_list],
                    'A_el [mm]': [round(cav.shape_space['OC'][0], 2) for cav in self.cavities_list],
                    'B_el [mm]': [round(cav.shape_space['OC'][1], 2) for cav in self.cavities_list],
                    'a_el [mm]': [round(cav.shape_space['OC'][2], 2) for cav in self.cavities_list],
                    'b_el [mm]': [round(cav.shape_space['OC'][3], 2) for cav in self.cavities_list],
                    'R_el [mm]': [round(cav.shape_space['OC'][4], 2) for cav in self.cavities_list],
                    'L_el [mm]': [round(cav.shape_space['OC'][5], 2) for cav in self.cavities_list],
                    # 'Req [mm]': [round(cav.shape_space['OC'][6], 2) for cav in self.cavities_list],
                    'alpha__el [deg]': [round(cav.shape_space['OC'][7], 2) for cav in self.cavities_list],
                    'A_er [mm]': [round(cav.shape_space['OC'][0], 2) for cav in self.cavities_list],
                    'B_er [mm]': [round(cav.shape_space['OC'][1], 2) for cav in self.cavities_list],
                    'a_er [mm]': [round(cav.shape_space['OC'][2], 2) for cav in self.cavities_list],
                    'b_er [mm]': [round(cav.shape_space['OC'][3], 2) for cav in self.cavities_list],
                    'R_er [mm]': [round(cav.shape_space['OC'][4], 2) for cav in self.cavities_list],
                    'L_er [mm]': [round(cav.shape_space['OC'][5], 2) for cav in self.cavities_list],
                    # 'Req [mm]': [round(cav.shape_space['OC'][6], 2) for cav in self.cavities_list],
                    'alpha_er [deg]': [round(cav.shape_space['OC'][7], 2) for cav in self.cavities_list],
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
                fr"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\Summaries\{self.folder}_excel_summary.xlsx",
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
        if self.folder != '':
            # check if folder exists
            if os.path.exists(fr"{self.folder}\PostProcessingData\Plots"):
                save_folder = fr"{self.folder}\PostProcessingData\Plots"
                plt.savefig(f"{save_folder}/{plot_name}")
            else:
                if not os.exists(fr"{self.folder}\PostProcessingData"):
                    os.mkdir(fr"{self.folder}\PostProcessingData")
                    os.mkdir(fr"{self.folder}\PostProcessingData\Plots")

                save_folder = fr"{self.folder}\PostProcessingData\Plots"
                os.mkdir(save_folder)
                plt.savefig(f"{save_folder}/{plot_name}")

    def calc_limits(self, which, selection):
        if self.operating_points is not None:
            # E0 = [45.6, 80, 120, 182.5]  # [GeV] Energy
            # nu_s = [0.025, 0.0506, 0.036, 0.087]  # Synchrotron oscillation tune
            # I0 = [1400, 135, 26.7, 5]  # [mA] Beam current5.4 * 2
            # alpha_c = [1.48, 1.48, 0.73, 0.73]  # [10−5] Momentum compaction factor
            # tau_z = [424.6, 78.7, 23.4, 6.8]  # [ms] Longitudinal damping time
            # tau_xy = [849.2, 157.4, 46.8, 13.6]  # [ms] Transverse damping time
            # f_rev = [3.07, 3.07, 3.07, 3.07]  # [kHz] Revolution frequency
            # beta_xy = 50
            # n_cav = [56, 112, 128, 140]  # 1_2_2_25

            # if selection is None or selection == '':
            #     if self.operating_points.empty:
            #         pass
            #     else:
            #         E0 = list(self.operating_points.loc["E [GeV]"])
            #         nu_s = list(self.operating_points.loc["nu_s []"])
            #         I0 = list(self.operating_points.loc["I0 [mA]"])
            #         alpha_c = list(self.operating_points.loc["alpha_p [1e-5]"])
            #         tau_z = list(self.operating_points.loc["tau_z [ms]"])
            #         tau_xy = list(self.operating_points.loc["tau_xy [ms]"])
            #         f_rev = list(self.operating_points.loc["f_rev [kHz]"])
            #         beta_xy = list(self.operating_points.loc["beta_xy [m]"])
            #         n_cav = list(self.operating_points.loc["N_c []"])
            # else:
            cols = selection

            E0 = list(self.operating_points.loc["E [GeV]", cols])
            nu_s = list(self.operating_points.loc["nu_s []", cols])
            I0 = list(self.operating_points.loc["I0 [mA]", cols])
            alpha_c = list(self.operating_points.loc["alpha_p [1e-5]", cols])
            tau_z = list(self.operating_points.loc["tau_z [ms]", cols])
            tau_xy = list(self.operating_points.loc["tau_xy [ms]", cols])
            f_rev = list(self.operating_points.loc["f_rev [kHz]", cols])
            beta_xy = list(self.operating_points.loc["beta_xy [m]", cols])
            n_cav = list(self.operating_points.loc["N_c []", cols])

            unit = {'MHz': 1e6,
                    'GHz': 1e9}

            Z_list, ZT_le = [], []
            if which == 'longitudinal':
                f_list = np.linspace(0, 10000, num=1000)
                Z_le = []
                try:
                    for i, n in enumerate(n_cav):
                        Z = [(2 * E0[i] * 1e9 * nu_s[i])
                             / (n * I0[i] * 1e-3 * alpha_c[i] * 1e-5 * tau_z[i] * 1e-3 * f)
                             if f > 1e-8 else 1e5 for f in f_list]

                        Z_le.append(round((2 * E0[i] * 1e9 * nu_s[i])
                                          / (n * I0[i] * 1e-3 * alpha_c[i] * 1e-5
                                             * tau_z[i] * 1e-3) * 1e-9 * 1e-3, 2))
                        Z_list.append(np.array(Z) * 1e-3)  # convert to kOhm

                except ZeroDivisionError:
                    print("ZeroDivisionError, check input")
                return f_list, Z_list, cols

            elif which == 'transversal':
                f_list = np.linspace(0, 10000, num=1000)
                try:
                    for i, n in enumerate(n_cav):
                        ZT = (2 * E0[i]) * 1e9 / (
                                n * I0[i] * 1e-3 * beta_xy[i] * tau_xy[i] * 1e-3 * f_rev[i] * 1e3)
                        ZT_le.append(round(ZT * 1e-3, 2))
                        Z_list.append(np.array(ZT) * 1e-3)  # convert to kOhm/m

                except ZeroDivisionError:
                    print("ZeroDivisionError, check input")

                return f_list, Z_list, cols
        else:
            print("Please load a valid operating point(s) file.")

    def plot_thresholds(self, which, selection, ax=None):
        self.threshold_texts_objects = []
        if ax is None:
            fig, ax = plt.subplots()
            ax.margins(x=0)

        if which.lower() in 'longitudinal':
            # calculate limits
            f_list, Z_list, labels = self.calc_limits('longitudinal', selection)

            # plot baselines
            for i, (z, lab) in enumerate(zip(Z_list, labels)):
                aa = ax.plot(f_list, z, ls='--', c='k')

            for i, (z, lab) in enumerate(zip(Z_list, labels)):
                pos = axis_data_coords_sys_transform(ax, 0.01, 0.5)
                indx = np.argmin(abs(f_list - pos[0]))
                x, y = axis_data_coords_sys_transform(ax, f_list[indx], z[indx], True)

                labe = lab.split('_')
                if len(labe) > 1:
                    txt = r"$\mathrm{" + fr"{lab.split('_')[0]}_" + "\mathrm{" + fr"{lab.split('_')[1]}" + r"}}$"
                else:
                    txt = r"$\mathrm{" + fr"{lab}" + r"}$"

                ab = add_text(ax, txt, box="Square", xy=(x, y), xycoords='axes fraction', size=12)
                self.threshold_texts_objects.append(ab)

        else:
            # calculate limits
            f_list, Z_list, labels = self.calc_limits('transversal', selection)

            # plot baselines
            for i, (z, lab) in enumerate(zip(Z_list, labels)):
                aa = ax.axhline(z, ls='--', c='k')

            for i, (z, lab) in enumerate(zip(Z_list, labels)):
                pos = axis_data_coords_sys_transform(ax, 0.01, 0.5)
                indx = np.argmin(abs(f_list - pos[0]))
                x, y = axis_data_coords_sys_transform(ax, f_list[indx], z, True)

                labe = lab.split('_')
                if len(labe) > 1:
                    txt = r"$\mathrm{" + fr"{lab.split('_')[0]}_" + "\mathrm{" + fr"{labels[i].split('_')[1]}" + r"}}$"
                else:
                    txt = r"$\mathrm{" + fr"{lab}" + r"}$"

                ab = add_text(ax, txt, box="Square", xy=(x, y), xycoords='axes fraction', size=12)
                self.threshold_texts_objects.append(ab)
        ax.autoscale(True, axis='y')

        # Adjust text annotations to avoid overlap
        # self._adjust_texts(threshold_texts_objects, ax)

    def _adjust_texts(self, texts, ax, separation=0.01, iterations=1):
        print('It is here')
        renderer = ax.figure.canvas.get_renderer()

        def get_text_bbox(text):
            bbox = text.get_window_extent(renderer=renderer)
            bbox_axes_coords = bbox.transformed(ax.transAxes.inverted())
            return bbox_axes_coords

        def calculate_adjustment(bbox1, bbox2):
            center1 = np.array([bbox1.x0 + bbox1.width / 2, bbox1.y0 + bbox1.height / 2])
            center2 = np.array([bbox2.x0 + bbox2.width / 2, bbox2.y0 + bbox2.height / 2])

            pos1 = np.array(text1.get_position())
            pos2 = np.array(text2.get_position())

            if bbox1.overlaps(bbox2):
                overlap_x = max(0, min(bbox1.x1, bbox2.x1) - max(bbox1.x0, bbox2.x0))
                overlap_y = max(0, min(bbox1.y1, bbox2.y1) - max(bbox1.y0, bbox2.y0))

                # if overlap_x > 0:
                #     midpoint_x = (center1[0] + center2[0]) / 2
                #     if bbox1.x1 > bbox2.x1:
                #         pos1[0] = midpoint_x + bbox1.width / 2 + separation/2
                #         pos2[0] = midpoint_x - bbox2.width / 2 - separation/2
                #     else:
                #         pos1[0] = midpoint_x - bbox1.width / 2 - separation/2
                #         pos2[0] = midpoint_x + bbox2.width / 2 + separation/2

                if overlap_y > 0:
                    midpoint_y = (center1[1] + center2[1]) / 2
                    if bbox1.y1 > bbox2.y1:
                        pos1[1] = midpoint_y + bbox1.height / 2 + separation / 2
                        pos2[1] = midpoint_y - bbox2.height / 2 - separation / 2
                    else:
                        pos1[1] = midpoint_y - bbox1.height / 2 - separation / 2
                        pos2[1] = midpoint_y + bbox2.height / 2 + separation / 2

            return pos1, pos2

        # Run the adjustment algorithm
        for _ in range(iterations):
            for i, text1 in enumerate(texts):
                bbox1 = get_text_bbox(text1)
                for j, text2 in enumerate(texts):
                    if i == j:
                        continue
                    bbox2 = get_text_bbox(text2)
                    print(bbox1, bbox2, text1.get_text(), text2.get_text())

                    if bbox1.overlaps(bbox2):
                        pos1, pos2 = calculate_adjustment(bbox1, bbox2)
                        text1.set_position((pos1[0], pos1[1]))
                        text2.set_position((pos2[0], pos2[1]))

        # Ensure the texts are within the plot limits
        for text in texts:
            x, y = text.get_position()
            x = np.clip(x, ax.get_xlim()[0], ax.get_xlim()[1])
            y = np.clip(y, ax.get_ylim()[0], ax.get_ylim()[1])
            text.set_position((x, y))

    def plot_cutoff(self, Ri_list, which, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.margins(x=0)

        f_list = self.calculate_beampipe_cutoff(Ri_list, which)

        for Ri in Ri_list:
            for mode_type, freq in zip(which, f_list[f'{Ri}']):
                vl = ax.axvline(freq, ls='--', c='k')  # label=f"{sc[0]} cutoff (Ri={sc[1]})",

                # get y text position from axis position. Position x is not used
                pos = axis_data_coords_sys_transform(ax, freq, 0.05, inverse=False)

                # ylim = self.ax.get_ylim()
                ab = add_text(ax, r"$f_\mathrm{c," + f"{mode_type}" + r"} (R_\mathrm{i} = "
                              + f"{Ri}" + r" ~\mathrm{mm}) $",
                              box="None", xy=(freq, pos[1]),
                              xycoords='data', size=10, rotation=90)

    @staticmethod
    def calculate_beampipe_cutoff(Ri_list, which):
        c = 299792458
        mode_list = {}
        f_list = {}

        mode_dict = {}
        for mode in which:
            mode_type = mode[:2]
            m = int(mode[2:3])
            n = int(mode[3:])
            mode_dict[mode] = {'type': mode_type, 'm': m, 'n': n}

        for Ri in Ri_list:
            f_list[f'{Ri}'] = []
            for mode, indices in mode_dict.items():
                m, n = indices['m'], indices['n']
                # get jacobian
                if mode_dict[mode]['type'].lower() == 'tm':
                    j = jn_zeros(m, n)[n - 1]
                else:
                    j = jnp_zeros(m, n)[n - 1]

                f = c / (2 * np.pi) * (j / (Ri * 1e-3))

                # append to f list
                f_list[f'{Ri}'].append(f * 1e-6)

        return f_list

    def run_optimisation(self, config):
        self.start_optimisation(config)

    @staticmethod
    def calc_cutoff(Ri, mode):
        # calculate frequency from Ri
        p_TM01, p_TE11 = 2.405, 1.841
        c = 299792458  # m/s

        if mode == 'TM01':
            freq = (c * p_TM01) / (2 * np.pi * Ri * 1e9) * 1e3
        else:
            freq = (c * p_TE11) / (2 * np.pi * Ri * 1e9) * 1e3

        return freq

    def __str__(self):
        return fr"{self.cavities_list}"

    def __getitem__(self, index):
        return self.cavities_list[index]


class Pillbox(Cavity):
    def __init__(self, n_cells, L, Req, Ri, S, L_bp, beampipe='none'):
        self.n_cells = n_cells
        self.n_modes = n_cells + 1
        self.n_modules = 1  # change later to enable module run
        self.L = L
        self.Req = Req
        self.Ri = Ri
        self.S = S
        self.L_bp = L_bp
        self.beampipe = beampipe
        self.bc = 33

        self.shape_space = {
            'IC': [L, Req, Ri, S, L_bp],
            'BP': beampipe
        }

    def plot(self, what, ax=None, **kwargs):
        if what.lower() == 'geometry':
            ax = plot_pillbox_geometry(self.n_cells, self.L, self.Req, self.Ri, self.S, self.L_bp, self.beampipe)
            ax.set_xlabel('$z$ [m]')
            ax.set_ylabel(r"$r$ [m]")
            return ax

        if what.lower() == 'zl':
            if ax:
                x, y, _ = self.abci_data['Long'].get_data('Longitudinal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            else:
                fig, ax = plt.subplots(figsize=(12, 4))
                ax.margins(x=0)
                x, y, _ = self.abci_data['Long'].get_data('Longitudinal Impedance Magnitude')
                ax.plot(x * 1e3, y)

            ax.set_xlabel('f [MHz]')
            ax.set_ylabel(r"$Z_{\parallel} ~[\mathrm{k\Omega}]$")
            return ax
        if what.lower() == 'zt':
            if ax:
                x, y, _ = self.abci_data['Trans'].get_data('Transversal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            else:
                fig, ax = plt.subplots(figsize=(12, 4))
                ax.margins(x=0)
                x, y, _ = self.abci_data['Trans'].get_data('Transversal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            ax.set_xlabel('f [MHz]')
            ax.set_ylabel(r"$Z_{\perp} ~[\mathrm{k\Omega/m}]$")
            return ax

        if what.lower() == 'convergence':
            try:
                if ax:
                    self._plot_convergence(ax)
                else:
                    fig, ax = plt.subplot_mosaic([['conv', 'abs_err']], layout='constrained', figsize=(12, 4))
                    self._plot_convergence(ax)
                return ax
            except ValueError:
                info("Convergence data not available.")

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
        cell_type: {'mid cell', 'end-mid cell', 'mid-end cell', 'end-end cell', 'single cell'}
            Type of cell to tune
        tune_variable: {'Req', 'L'}
            Tune variable. Currently supports only the tuning of the equator radius ``Req`` and half-cell length ``L``

        Returns
        -------

        """

        for _ in tqdm([1]):
            iter_set = ['Linear Interpolation', TUNE_ACCURACY, 10]

            if freq is None:
                # calculate freq from mid cell length
                beta = 1
                freq = beta * c0 / (4 * self.mid_cell[5])
                info("Calculated freq from mid cell half length: ", freq)

            # create new shape space based on cell type
            # if cell_type.lower() == 'mid cell':
            shape_space = {
                f'{self.name}':
                    {
                        'IC': self.shape_space['IC'],
                        'OC': self.shape_space['OC'],
                        'OC_R': self.shape_space['OC_R'],
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
                                            progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
                else:
                    if run_tune.lower() == 'y':
                        # copy files required for simulation
                        self._overwriteFolder(proc, self.folder, self.name)
                        self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                        self.run_tune_ngsolve(shape_space, resume, proc, self.bc,
                                              SOFTWARE_DIRECTORY, self.folder, self.name,
                                              tune_variable, iter_set, cell_type,
                                              progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_ngsolve_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
            else:
                if solver.lower() == 'slans':
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_slans(shape_space, resume, proc, self.bc,
                                        SOFTWARE_DIRECTORY, self.folder, self.name, tuner,
                                        tune_variable, iter_set, cell_type,
                                        progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Oops! Something went wrong. Could not find the tune results. Please run tune again.")
                else:
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_ngsolve(shape_space, resume, proc, self.bc,
                                          SOFTWARE_DIRECTORY, self.folder, self.name,
                                          tune_variable, iter_set, cell_type,
                                          progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)
                    try:
                        self.get_ngsolve_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Oops! Something went wrong. Could not find the tune results. Please run tune again.")

    @staticmethod
    def run_tune_ngsolve(shape, resume, p, bc, parentDir, projectDir, filename,
                         tune_variable, iter_set, cell_type, progress_list, convergence_list, n_cells):
        tuner.tune_ngsolve(shape, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                           tune_variable=tune_variable, iter_set=iter_set,
                           cell_type=cell_type, sim_folder='Optimisation',
                           progress_list=progress_list, convergence_list=convergence_list,
                           save_last=True,
                           n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    def run_eigenmode(self, solver='ngsolve', freq_shift=0, boundary_cond=None, subdir='',
                      uq_config=None):
        """
        Run eigenmode analysis on cavity

        Parameters
        ----------
        solver: {'SLANS', 'NGSolve'}
            Solver to be used. Native solver is still under development. Results are not as accurate as that of SLANS.
        freq_shift:
            Frequency shift. Eigenmode solver searches for eigenfrequencies around this value
        boundary_cond: int
            Boundary condition of left and right cell/beampipe ends
        subdir: str
            Sub directory to save results to
        uq_config: None | dict
            Provides inputs required for uncertainty quantification. Default is None and disables uncertainty quantification.

        Returns
        -------

        """

        for _ in tqdm([1], file=sys.stdout):
            if boundary_cond:
                self.bc = boundary_cond

            self._run_ngsolve(self.name, self.n_cells, self.n_modules, self.shape_space, self.n_modes, freq_shift,
                              self.bc, SOFTWARE_DIRECTORY, self.folder, sub_dir='', uq_config=uq_config)
            # load quantities of interest
            try:
                self.get_eigenmode_qois('NGSolveMEVP')
            except FileNotFoundError:
                error("Could not find eigenmode results. Please rerun eigenmode analysis.")

    def run_wakefield(self, MROT=2, MT=10, NFS=10000, wakelength=50, bunch_length=25,
                      DDR_SIG=0.1, DDZ_SIG=0.1, WG_M=None, marker='', operating_points=None, solver='ABCI'):
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

        if operating_points is None:
            wp_dict = {}
        exist = False

        if not exist:
            if solver == 'ABCI':
                self._run_abci(self.name, self.n_cells, self.n_modules, self.shape_space,
                               MROT=MROT, MT=MT, NFS=NFS, UBT=wakelength, bunch_length=bunch_length,
                               DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG,
                               parentDir=SOFTWARE_DIRECTORY, projectDir=self.folder, WG_M=WG_M, marker=marker,
                               operating_points=operating_points, freq=self.freq, R_Q=self.R_Q)

                try:
                    self.get_abci_data()
                    self.get_abci_qois()
                except FileNotFoundError:
                    error("Could not find the abci wakefield results. Please rerun wakefield analysis.")

        else:
            try:
                self.get_abci_data()
                self.get_abci_qois()
            except FileNotFoundError:
                error("Could not find the abci wakefield results. Please rerun wakefield analysis.")

    @staticmethod
    def _run_ngsolve(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir, sub_dir='',
                     uq_config=None):
        start_time = time.time()
        # create folders for all keys
        ngsolve_mevp.createFolder(name, projectDir, subdir=sub_dir)
        ngsolve_mevp.pillbox(n_cells, n_modules, shape['IC'],
                             n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                             parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
        # run UQ
        if uq_config:
            uq(name, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
               n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
               f_shift=f_shift, bc=bc, parentDir=parentDir, projectDir=projectDir)

        done(f'Done with Cavity {name}. Time: {time.time() - start_time}')

    @staticmethod
    def _run_abci(name, n_cells, n_modules, shape, MROT=0, MT=4.0, NFS=10000, UBT=50.0, bunch_length=20.0,
                  DDR_SIG=0.1, DDZ_SIG=0.1,
                  parentDir=None, projectDir=None,
                  WG_M=None, marker='', operating_points=None, freq=0, R_Q=0):

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
        if len(operating_points.keys()) > 0:
            try:
                if freq != 0 and R_Q != 0:
                    d = {}
                    # save qois
                    for key, vals in tqdm(operating_points.items()):
                        WP = key
                        I0 = float(vals['I0 [mA]'])
                        Nb = float(vals['Nb [1e11]'])
                        sigma_z = [float(vals["sigma_SR [mm]"]), float(vals["sigma_BS [mm]"])]
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


class RFGun(Cavity):
    def __init__(self):
        # self.shape_space = {
        #     'IC': [L, Req, Ri, S, L_bp],
        #     'BP': beampipe
        # }
        pass

    def plot(self, what, ax=None, **kwargs):
        if what.lower() == 'geometry':
            ax = plot_gun()
            ax.set_xlabel('$z$ [m]')
            ax.set_ylabel(r"$r$ [m]")
            return ax

        if what.lower() == 'zl':
            if ax:
                x, y, _ = self.abci_data['Long'].get_data('Longitudinal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            else:
                fig, ax = plt.subplots(figsize=(12, 4))
                ax.margins(x=0)
                x, y, _ = self.abci_data['Long'].get_data('Longitudinal Impedance Magnitude')
                ax.plot(x * 1e3, y)

            ax.set_xlabel('f [MHz]')
            ax.set_ylabel(r"$Z_{\parallel} ~[\mathrm{k\Omega}]$")
            return ax
        if what.lower() == 'zt':
            if ax:
                x, y, _ = self.abci_data['Trans'].get_data('Transversal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            else:
                fig, ax = plt.subplots(figsize=(12, 4))
                ax.margins(x=0)
                x, y, _ = self.abci_data['Trans'].get_data('Transversal Impedance Magnitude')
                ax.plot(x * 1e3, y)
            ax.set_xlabel('f [MHz]')
            ax.set_ylabel(r"$Z_{\perp} ~[\mathrm{k\Omega/m}]$")
            return ax

        if what.lower() == 'convergence':
            try:
                if ax:
                    self._plot_convergence(ax)
                else:
                    fig, ax = plt.subplot_mosaic([['conv', 'abs_err']], layout='constrained', figsize=(12, 4))
                    self._plot_convergence(ax)
                return ax
            except ValueError:
                info("Convergence data not available.")

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
        cell_type: {'mid cell', 'end-mid cell', 'mid-end cell', 'end-end cell', 'single cell'}
            Type of cell to tune
        tune_variable: {'Req', 'L'}
            Tune variable. Currently supports only the tuning of the equator radius ``Req`` and half-cell length ``L``

        Returns
        -------

        """

        for _ in tqdm([1]):
            iter_set = ['Linear Interpolation', TUNE_ACCURACY, 10]

            if freq is None:
                # calculate freq from mid cell length
                beta = 1
                freq = beta * c0 / (4 * self.mid_cell[5])
                info("Calculated freq from mid cell half length: ", freq)

            # create new shape space based on cell type
            # if cell_type.lower() == 'mid cell':
            shape_space = {
                f'{self.name}':
                    {
                        'IC': self.shape_space['IC'],
                        'OC': self.shape_space['OC'],
                        'OC_R': self.shape_space['OC_R'],
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
                                            progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
                else:
                    if run_tune.lower() == 'y':
                        # copy files required for simulation
                        self._overwriteFolder(proc, self.folder, self.name)
                        self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                        self.run_tune_ngsolve(shape_space, resume, proc, self.bc,
                                              SOFTWARE_DIRECTORY, self.folder, self.name,
                                              tune_variable, iter_set, cell_type,
                                              progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    # read tune results and update geometry
                    try:
                        self.get_ngsolve_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Could not find the tune results. Please run tune again.")
            else:
                if solver.lower() == 'slans':
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_slans(shape_space, resume, proc, self.bc,
                                        SOFTWARE_DIRECTORY, self.folder, self.name, tuner,
                                        tune_variable, iter_set, cell_type,
                                        progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)

                    try:
                        self.get_slans_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Oops! Something went wrong. Could not find the tune results. Please run tune again.")
                else:
                    # copy files required for simulation
                    self._overwriteFolder(proc, self.folder, self.name)
                    self._copyFiles(proc, SOFTWARE_DIRECTORY, self.folder, self.name)

                    self.run_tune_ngsolve(shape_space, resume, proc, self.bc,
                                          SOFTWARE_DIRECTORY, self.folder, self.name,
                                          tune_variable, iter_set, cell_type,
                                          progress_list=[], convergence_list=self.convergence_list, n_cells=n_cells)
                    try:
                        self.get_ngsolve_tune_res(tune_variable, cell_type)
                    except FileNotFoundError:
                        error("Oops! Something went wrong. Could not find the tune results. Please run tune again.")

    @staticmethod
    def run_tune_ngsolve(shape, resume, p, bc, parentDir, projectDir, filename,
                         tune_variable, iter_set, cell_type, progress_list, convergence_list, n_cells):
        tuner.tune_ngsolve(shape, bc, parentDir, projectDir, filename, resume=resume, proc=p,
                           tune_variable=tune_variable, iter_set=iter_set,
                           cell_type=cell_type, sim_folder='Optimisation',
                           progress_list=progress_list, convergence_list=convergence_list,
                           save_last=True,
                           n_cell_last_run=n_cells)  # last_key=last_key This would have to be tested again #val2

    def run_eigenmode(self, solver='ngsolve', freq_shift=0, boundary_cond=None, subdir='',
                      uq_config=None):
        """
        Run eigenmode analysis on cavity

        Parameters
        ----------
        solver: {'SLANS', 'NGSolve'}
            Solver to be used. Native solver is still under development. Results are not as accurate as that of SLANS.
        freq_shift:
            Frequency shift. Eigenmode solver searches for eigenfrequencies around this value
        boundary_cond: int
            Boundary condition of left and right cell/beampipe ends
        subdir: str
            Sub directory to save results to
        uq_config: None | dict
            Provides inputs required for uncertainty quantification. Default is None and disables uncertainty quantification.

        Returns
        -------

        """

        for _ in tqdm([1], file=sys.stdout):
            if boundary_cond:
                self.bc = boundary_cond

            self._run_ngsolve(self.name, self.n_cells, self.n_modules, self.shape_space, self.n_modes, freq_shift,
                              self.bc, SOFTWARE_DIRECTORY, self.folder, sub_dir='', uq_config=uq_config)
            # load quantities of interest
            try:
                self.get_eigenmode_qois('NGSolveMEVP')
            except FileNotFoundError:
                error("Could not find eigenmode results. Please rerun eigenmode analysis.")

    def run_wakefield(self, MROT=2, MT=10, NFS=10000, wakelength=50, bunch_length=25,
                      DDR_SIG=0.1, DDZ_SIG=0.1, WG_M=None, marker='', operating_points=None, solver='ABCI'):
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

        if operating_points is None:
            wp_dict = {}
        exist = False

        if not exist:
            if solver == 'ABCI':
                self._run_abci(self.name, self.n_cells, self.n_modules, self.shape_space,
                               MROT=MROT, MT=MT, NFS=NFS, UBT=wakelength, bunch_length=bunch_length,
                               DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG,
                               parentDir=SOFTWARE_DIRECTORY, projectDir=self.folder, WG_M=WG_M, marker=marker,
                               operating_points=operating_points, freq=self.freq, R_Q=self.R_Q)

                try:
                    self.get_abci_data()
                    self.get_abci_qois()
                except FileNotFoundError:
                    error("Could not find the abci wakefield results. Please rerun wakefield analysis.")

        else:
            try:
                self.get_abci_data()
                self.get_abci_qois()
            except FileNotFoundError:
                error("Could not find the abci wakefield results. Please rerun wakefield analysis.")

    @staticmethod
    def _run_ngsolve(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir, sub_dir='',
                     uq_config=None):
        start_time = time.time()
        # create folders for all keys
        ngsolve_mevp.createFolder(name, projectDir, subdir=sub_dir)
        ngsolve_mevp.pillbox(n_cells, n_modules, shape['IC'],
                             n_modes=n_modes, fid=f"{name}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                             parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
        # run UQ
        if uq_config:
            uq(name, shape, ["freq", "R/Q", "Epk/Eacc", "Bpk/Eacc"],
               n_cells=n_cells, n_modules=n_modules, n_modes=n_modes,
               f_shift=f_shift, bc=bc, parentDir=parentDir, projectDir=projectDir)

        done(f'Done with Cavity {name}. Time: {time.time() - start_time}')

    @staticmethod
    def _run_abci(name, n_cells, n_modules, shape, MROT=0, MT=4.0, NFS=10000, UBT=50.0, bunch_length=20.0,
                  DDR_SIG=0.1, DDZ_SIG=0.1,
                  parentDir=None, projectDir=None,
                  WG_M=None, marker='', operating_points=None, freq=0, R_Q=0):

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
        if len(operating_points.keys()) > 0:
            try:
                if freq != 0 and R_Q != 0:
                    d = {}
                    # save qois
                    for key, vals in tqdm(operating_points.items()):
                        WP = key
                        I0 = float(vals['I0 [mA]'])
                        Nb = float(vals['Nb [1e11]'])
                        sigma_z = [float(vals["sigma_SR [mm]"]), float(vals["sigma_BS [mm]"])]
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


class Dakota:
    def __init__(self, folder, name, scripts_folder=None):
        self.nodes = None
        self.sim_results = None

        assert f'/' in folder, error('Please ensure directory paths use forward slashes.')
        self.folder = folder
        if scripts_folder is None:
            self.scripts_folder = r'D:/Dropbox/CavityDesignHub/analysis_modules/uq/dakota_scripts'
        else:
            self.scripts_folder = scripts_folder
        self.name = name

    def write_input_file(self, **kwargs):
        keys = kwargs.keys()
        assert 'variables_config' in keys, error('please enter keyword "variables config"')
        assert 'interface_config' in keys, error('please enter keyword "interface config"')
        assert 'method_config' in keys, error('please enter keyword "method config"')

        variables_config = kwargs['variables_config']
        interface_config = kwargs['interface_config']
        method_config = kwargs['method_config']

        # check if folder exists, if not, create folder
        if not os.path.exists(os.path.join(self.folder, self.name)):
            try:
                os.mkdir(os.path.join(self.folder, self.name))
            except FileExistsError:
                print("Could not create folder. Make sure target location exists.")

        with open(os.path.join(self.folder, self.name, f'{self.name}.in'), 'w') as f:
            self.environment(f)
            self.method(f, **method_config)
            self.variables(f, **variables_config)
            self.interface(f, **interface_config)

    def environment(self, f):
        f.write('environment\n')
        f.write('\ttabular_data\n')
        f.write("\t\ttabular_data_file = 'sim_result_table.dat'\n")
        f.write('\tresults_output\n')
        f.write("\t\tresults_output_file = 'result_output_file.dat'\n")

        f.write('\n')

    def method(self, f, **kwargs):
        keys = kwargs.keys()
        assert 'method' in keys, error('Please enter "method" in "method config".')
        method = kwargs['method']

        f.write('method\n')
        f.write(f'\t{method}\n')
        f.write("\t\texport_expansion_file ='expansion_file.dat'\n")
        f.write("\t\tcubature_integrand = 3\n")
        f.write("\t\tsamples_on_emulator = 10000\n")
        f.write("\t\tseed = 12347\n")
        f.write("\t\tvariance_based_decomp interaction_order = 1\n")

        f.write('\n')

    def variables(self, f, **kwargs):
        """

        Parameters
        ----------
        f: File
            File
        kind: {'uniform_uncertain', 'normal_uncertain', 'beta_uncertain'}
            Type of distribution of the variables
        descriptors: list, ndarray
            Uncertain variable names
        kwargs: kwargs
            Other relevant arguments eg. {means: [], 'std_deviations': [], 'lower_bounds': [], 'upper_bounds': []}


        Returns
        -------

        """

        keys = kwargs.keys()
        assert 'kind' in keys, error('Please enter "kind"')
        assert 'lower_bounds' in keys, error('Please enter keyword "lower_bounds"')
        assert 'upper_bounds' in keys, error('Please enter keyword "upper bounds"')
        kind = kwargs['kind']
        lower_bounds = kwargs['lower_bounds']
        upper_bounds = kwargs['upper_bounds']

        assert len(upper_bounds) == len(lower_bounds), error("Length of upper and lower bounds must be equal.")

        if "descriptors" in keys:
            descriptors = kwargs['descriptors']
        else:
            info('"descriptors" not entered. Using default parameter labelling.')
            descriptors = [f'p{n}' for n in range(len(upper_bounds))]
        assert len(descriptors) == len(kwargs['upper_bounds'])

        f.write("variables\n")
        f.write(f"\t{kind} = {len(descriptors)}\n")
        f.write(
            "\tdescriptors       =   " + '\t\t\t'.join(['"' + descriptor + '"' for descriptor in descriptors]) + '\n')

        if 'means' in kwargs.keys():
            assert len(descriptors) == len(kwargs['means'])
            f.write("\tmeans      =   " + '\t\t\t'.join([str(mean) for mean in kwargs['means']]) + '\n')

        if 'std_deviations' in kwargs.keys():
            assert len(descriptors) == len(kwargs['std_deviations'])
            f.write("\tstd_deviations      =   " + '\t\t\t'.join([str(std) for std in kwargs['std_deviations']]) + '\n')

        if 'lower_bounds' in kwargs.keys():
            f.write("\tlower_bounds      =   " + '\t\t\t'.join([str(lb) for lb in kwargs['lower_bounds']]) + '\n')

        if 'upper_bounds' in kwargs.keys():
            f.write("\tupper_bounds      =   " + '\t\t\t'.join([str(ub) for ub in kwargs['upper_bounds']]) + '\n')

        f.write('\n')

    def interface(self, f, **kwargs):

        keys = kwargs.keys()
        assert 'analysis_driver' in keys, error('please enter keyword "analysis driver"')
        assert 'responses' in keys, error('Please enter "responses"')

        analysis_driver = kwargs['analysis_driver']
        responses = kwargs['responses']

        nodes_only = False
        if 'nodes_only' in keys:
            nodes_only = kwargs['nodes_only']
            responses = ['f1']

        processes = 1
        if 'processes' in keys:
            processes = kwargs['processes']

        f.write("interface\n")
        f.write("#\tcommon options\n")
        f.write("#\tfork\n")
        f.write("\tparameters_file = 'params.in'\n")
        f.write("\tresults_file    = 'results.out'\n")
        f.write(f"\tsystem asynchronous evaluation_concurrency = {processes}\n")
        f.write(f"\tanalysis_driver = '{analysis_driver} {len(responses)} {nodes_only} {self.scripts_folder}'\n")
        f.write("#\tparameters_file = 'params.in'\n")
        f.write("#\tresults_file    = 'results.out'\n")
        f.write("#\tfile_tag\n")
        f.write("#\tfile_save\n")
        f.write("#\taprepro\n")
        f.write('\n')

        self.responses(f, responses)

    def responses(self, f, responses):
        f.write("responses\n")
        f.write(f"\tresponse_functions = {len(responses)}\n")
        f.write("\tno_gradients\n")
        f.write("\tno_hessians\n")

        f.write('\n')

    def nodes_to_cst_sweep_input(self, partitions=1):
        # save parts
        row_partition = len(self.nodes.index) // partitions
        for i in range(partitions):
            if i < partitions - 1:
                df_part = self.nodes.loc[i * row_partition:(i + 1) * row_partition - 1]
            else:
                df_part = self.nodes.loc[i * row_partition:]

            df_part.to_csv(fr"{self.folder}/{self.name}/cst_sweep_files/cst_par_in_{i + 1}.txt", sep="\t", index=None)

    def run_analysis(self, write_cst=True, partitions=1):
        cwd = os.path.join(self.folder, self.name)
        dakota_in = f'{os.path.join(self.folder, self.name, f"{self.name}.in")}'
        dakota_out = f'{os.path.join(self.folder, self.name, f"{self.name}.out")}'

        subprocess.run(['dakota', '-i', dakota_in, '-o', dakota_out], cwd=cwd, shell=True)

        # read results
        filepath = fr"{self.folder}/{self.name}/sim_result_table.dat"
        self.sim_results = pd.read_csv(filepath, sep='\s+')

        # delete unnecessary columns
        self.nodes = self.sim_results.drop(self.sim_results.filter(regex='response|interface|eval_id').columns, axis=1)
        self.sim_results.to_excel(fr"{self.folder}/{self.name}/nodes.xlsx", index=False)

        if write_cst:
            # check if folder exist and clear
            if os.path.exists(os.path.join(self.folder, self.name, 'cst_sweep_files')):
                shutil.rmtree(os.path.join(self.folder, self.name, 'cst_sweep_files'))
                os.mkdir(os.path.join(self.folder, self.name, 'cst_sweep_files'))
            else:
                os.mkdir(os.path.join(self.folder, self.name, 'cst_sweep_files'))

            # post processes
            self.nodes_to_cst_sweep_input(partitions)
        else:
            if os.path.exists(os.path.join(self.folder, self.name, 'cst_sweep_files')):
                shutil.rmtree(os.path.join(self.folder, self.name, 'cst_sweep_files'))


class OperationPoints:
    def __init__(self, filepath=None):
        self.op_points = {}

        if filepath:
            if os.path.exists(filepath):
                self.op_points = self.load_operating_point(filepath)

    def load_operating_point(self, filepath):
        with open(filepath, 'r') as f:
            op_points = json.load(f)

        self.op_points = op_points
        return op_points

    def get_default_operating_points(self):
        self.op_points = pd.DataFrame({
            "Z_2023": {
                "freq [MHz]": 400.79,
                "E [GeV]": 45.6,
                "I0 [mA]": 1280,
                "V [GV]": 0.12,
                "Eacc [MV/m]": 5.72,
                "nu_s []": 0.0370,
                "alpha_p [1e-5]": 2.85,
                "tau_z [ms]": 354.91,
                "tau_xy [ms]": 709.82,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 56,
                "N_c []": 56,
                "T [K]": 4.5,
                "sigma_SR [mm]": 4.32,
                "sigma_BS [mm]": 15.2,
                "Nb [1e11]": 2.76
            },
            "W_2023": {
                "freq [MHz]": 400.79,
                "E [GeV]": 80,
                "I0 [mA]": 135,
                "V [GV]": 1.0,
                "Eacc [MV/m]": 10.61,
                "nu_s []": 0.0801,
                "alpha_p [1e-5]": 2.85,
                "tau_z [ms]": 65.99,
                "tau_xy [ms]": 131.98,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 132,
                "T [K]": 4.5,
                "sigma_SR [mm]": 3.55,
                "sigma_BS [mm]": 7.02,
                "Nb [1e11]": 2.29
            },
            "H_2023": {
                "freq [MHz]": 400.79,
                "E [GeV]": 120,
                "I0 [mA]": 53.4,
                "V [GV]": 2.1,
                "Eacc [MV/m]": 10.61,
                "nu_s []": 0.0328,
                "alpha_p [1e-5]": 0.733,
                "tau_z [ms]": 19.6,
                "tau_xy [ms]": 39.2,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 528,
                "T [K]": 4.5,
                "sigma_SR [mm]": 2.5,
                "sigma_BS [mm]": 4.45,
                "Nb [1e11]": 1.51
            },
            "ttbar_2023": {
                "freq [MHz]": 801.58,
                "E [GeV]": 182.5,
                "I0 [mA]": 10,
                "V [GV]": 9.2,
                "Eacc [MV/m]": 20.12,
                "nu_s []": 0.0826,
                "alpha_p [1e-5]": 0.733,
                "tau_z [ms]": 5.63,
                "tau_xy [ms]": 11.26,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 488,
                "T [K]": 2,
                "sigma_SR [mm]": 1.67,
                "sigma_BS [mm]": 2.54,
                "Nb [1e11]": 2.26
            },
            "Z_2022": {
                "freq [MHz]": 400.79,
                "E [GeV]": 45.6,
                "I0 [mA]": 1400,
                "V [GV]": 0.12,
                "Eacc [MV/m]": 5.72,
                "nu_s []": 0.0370,
                "alpha_p [1e-5]": 2.85,
                "tau_z [ms]": 354.91,
                "tau_xy [ms]": 709.82,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 56,
                "T [K]": 4.5,
                "sigma_SR [mm]": 4.32,
                "sigma_BS [mm]": 15.2,
                "Nb [1e11]": 2.76
            },
            "W_2022": {
                "freq [MHz]": 400.79,
                "E [GeV]": 80,
                "I0 [mA]": 135,
                "V [GV]": 1.0,
                "Eacc [MV/m]": 11.91,
                "nu_s []": 0.0801,
                "alpha_p [1e-5]": 2.85,
                "tau_z [ms]": 65.99,
                "tau_xy [ms]": 131.98,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 112,
                "T [K]": 4.5,
                "sigma_SR [mm]": 3.55,
                "sigma_BS [mm]": 7.02,
                "Nb [1e11]": 2.29
            },
            "H_2022": {
                "freq [MHz]": 400.79,
                "E [GeV]": 120,
                "I0 [mA]": 53.4,
                "V [GV]": 2.1,
                "Eacc [MV/m]": 10.61,
                "nu_s []": 0.0328,
                "alpha_p [1e-5]": 0.733,
                "tau_z [ms]": 19.6,
                "tau_xy [ms]": 39.2,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 528,
                "T [K]": 4.5,
                "sigma_SR [mm]": 2.5,
                "sigma_BS [mm]": 4.45,
                "Nb [1e11]": 1.51
            },
            "ttbar_2022": {
                "freq [MHz]": 801.58,
                "E [GeV]": 182.5,
                "I0 [mA]": 10,
                "V [GV]": 9.2,
                "Eacc [MV/m]": 20.12,
                "nu_s []": 0.0826,
                "alpha_p [1e-5]": 0.733,
                "tau_z [ms]": 5.63,
                "tau_xy [ms]": 11.26,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 488,
                "T [K]": 2,
                "sigma_SR [mm]": 1.67,
                "sigma_BS [mm]": 2.54,
                "Nb [1e11]": 2.26
            },
            "Z_booster_2022": {
                "freq [MHz]": 801.58,
                "E [GeV]": 45.6,
                "I0 [mA]": 128,
                "V [GV]": 0.14,
                "Eacc [MV/m]": 6.23,
                "nu_s []": 0.0370,
                "alpha_p [1e-5]": 2.85,
                "tau_z [ms]": 354.91,
                "tau_xy [ms]": 709.82,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 120,
                "T [K]": 4.5,
                "sigma_SR [mm]": 4.32,
                "sigma_BS [mm]": 15.2,
                "Nb [1e11]": 0.276
            },
            "Z_2018": {
                "freq [MHz]": 400.79,
                "E [GeV]": 45.6,
                "I0 [mA]": 1390,
                "V [GV]": 0.10,
                "Eacc [MV/m]": 5.72,
                "nu_s []": 0.025,
                "alpha_p [1e-5]": 1.48,
                "tau_z [ms]": 424.6,
                "tau_xy [ms]": 849.2,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 52,
                "T [K]": 4.5,
                "sigma_SR [mm]": 3.5,
                "sigma_BS [mm]": 12.1,
                "Nb [1e11]": 1.7
            },
            "W_2018": {
                "freq [MHz]": 400.79,
                "E [GeV]": 80,
                "I0 [mA]": 147,
                "V [GV]": 0.75,
                "Eacc [MV/m]": 11.91,
                "nu_s []": 0.0506,
                "alpha_p [1e-5]": 1.48,
                "tau_z [ms]": 78.7,
                "tau_xy [ms]": 157.4,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 52,
                "T [K]": 4.5,
                "sigma_SR [mm]": 3.0,
                "sigma_BS [mm]": 6.0,
                "Nb [1e11]": 1.5
            },
            "H_2018": {
                "freq [MHz]": 400.79,
                "E [GeV]": 120,
                "I0 [mA]": 29,
                "V [GV]": 2.0,
                "Eacc [MV/m]": 11.87,
                "nu_s []": 0.036,
                "alpha_p [1e-5]": 0.73,
                "tau_z [ms]": 23.4,
                "tau_xy [ms]": 46.8,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 136,
                "T [K]": 4.5,
                "sigma_SR [mm]": 3.15,
                "sigma_BS [mm]": 5.3,
                "Nb [1e11]": 1.8
            },
            "ttbar_2018": {
                "freq [MHz]": 801.58,
                "E [GeV]": 182.5,
                "I0 [mA]": 10.8,
                "V [GV]": 10.93,
                "Eacc [MV/m]": 24.72,
                "nu_s []": 0.087,
                "alpha_p [1e-5]": 0.73,
                "tau_z [ms]": 6.8,
                "tau_xy [ms]": 13.6,
                "f_rev [kHz]": 3.07,
                "beta_xy [m]": 50,
                "N_c []": 584,
                "T [K]": 2,
                "sigma_SR [mm]": 1.97,
                "sigma_BS [mm]": 2.54,
                "Nb [1e11]": 2.3
            }
        })


def uq(shape_space, objectives, solver_dict, solver_args_dict, uq_config):
    """

    Parameters
    ----------
    shape_space: dict
        Cavity geometry parameter space
    objectives: list | ndarray
        Array of objective functions
    solver_dict: dict
        Python dictionary of solver settings
    solver_args_dict: dict
        Python dictionary of solver arguments
    uq_config:
        Python dictionary of uncertainty quantification settings

    Returns
    -------

    """
    parentDir = solver_args_dict['parentDir']
    projectDir = solver_args_dict['projectDir']
    cell_type = uq_config['cell type']
    analysis_folder = solver_args_dict['analysis folder']
    opt = solver_args_dict['optimisation']
    delta = uq_config['delta']
    method = uq_config['method']
    uq_vars = uq_config['variables']
    assert len(uq_vars) == len(delta), error('Ensure number of variables equal number of deltas')

    for key, shape in shape_space.items():
        err = False
        result_dict_eigen, result_dict_abci = {}, {}
        run_eigen, run_abci = False, False
        eigen_obj_list, abci_obj_list = [], []

        for o in objectives:
            if o in ["Req", "freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]",
                     "G [Ohm]", "Q []", 'kcc [%]', "ff [%]"]:
                result_dict_eigen[o] = {'expe': [], 'stdDev': []}
                run_eigen = True
                eigen_obj_list.append(o)

            if o.split(' ')[0] in ['ZL', 'ZT', 'k_loss', 'k_kick']:
                result_dict_abci[o] = {'expe': [], 'stdDev': []}
                run_abci = True
                abci_obj_list.append(o)

        rdim = len(uq_vars)
        degree = 1

        if method[1].lower() == 'stroud3':
            nodes, weights, bpoly = quad_stroud3(rdim, degree)
            nodes = 2. * nodes - 1.
            # nodes, weights = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer. the nodes are not the same as the custom function
        elif method[1].lower() == 'stroud5':
            nodes, weights = cn_leg_05_2(rdim)
        elif method[1].lower() == 'gaussian':
            nodes, weights = cn_gauss(rdim, 2)
        elif method[1].lower() == 'lhs':
            sampler = qmc.LatinHypercube(d=rdim)
            _ = sampler.reset()
            nsamp = uq_config['integration'][2]
            sample = sampler.random(n=nsamp)

            l_bounds = [-1 for _ in range(len(uq_vars))]
            u_bounds = [1 for _ in range(len(uq_vars))]
            sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

            nodes, weights = sample_scaled.T, np.ones((nsamp, 1))
        elif method[0].lower() == 'from file':
            if len(method) == 2:
                nodes = pd.read_csv(method[1], sep='\s+').iloc[:, method[1]]
            else:
                nodes = pd.read_csv(method[1], sep='\s+')

            nodes = nodes.to_numpy().T
            weights = np.ones((nodes.shape[1], 1))
        else:
            # issue warning
            warning('Integration method not recognised. Defaulting to Stroud3 quadrature rule!')
            nodes, weights, bpoly = quad_stroud3(rdim, degree)
            nodes = 2. * nodes - 1.

        # save nodes
        data_table = pd.DataFrame(nodes.T, columns=uq_vars)
        data_table.to_csv(fr'{projectDir}\SimulationData\{analysis_folder}\{key}\nodes.csv',
                          index=False, sep='\t', float_format='%.32f')

        if cell_type.lower() == 'mid cell' or cell_type.lower() == 'mid-cell' or cell_type.lower() == 'mid_cell':
            cell_node = shape['IC']
        elif cell_type.lower() == 'mid-end cell' or cell_type.lower() == 'mid-end-cell' or cell_type.lower() == 'mid_end_cell':
            cell_node = shape['OC']
        elif (cell_type.lower() == 'end-end cell' or cell_type.lower() == 'end-end-cell'
              or cell_type.lower() == 'end_end_cell') or cell_type.lower() == 'end end cell':
            cell_node = shape['OC']
        else:
            cell_node = shape['OC']

        perturbed_cell_node = np.array(cell_node)
        no_parm, no_sims = np.shape(nodes)
        if delta is None:
            delta = [0.05 for _ in range(len(uq_vars))]

        if run_eigen:
            Ttab_val_f = []
            solver, solver_args = solver_dict['ngsolvemevp'], solver_args_dict['ngsolvemevp']
            bc = solver_args['bc']
            # beampipes = solver_args['beampipes']
            # norm_length = solver_args['norm_length']
            n_cells = solver_args['n_cells']
            # proc = solver_args['proc']
            sub_dir = fr'{key}'

            for i in range(no_sims):
                skip = False
                for j, uq_var in enumerate(uq_vars):
                    uq_var_indx = VAR_TO_INDEX_DICT[uq_var]
                    perturbed_cell_node[uq_var_indx] = cell_node[uq_var_indx] * (1 + delta[j] * nodes[j, i])

                if cell_type.lower() == 'mid cell' or cell_type.lower() == 'mid-cell' or cell_type.lower() == 'mid_cell':
                    cell_node = shape['IC']
                    mid = perturbed_cell_node
                    left = perturbed_cell_node
                    right = perturbed_cell_node
                    beampipes = 'none'
                elif cell_type.lower() == 'mid-end cell' or cell_type.lower() == 'mid-end-cell' or cell_type.lower() == 'mid_end_cell':
                    mid = shape['IC']
                    left = shape['IC']
                    right = perturbed_cell_node
                    beampipes = 'right'
                elif (cell_type.lower() == 'end-end cell' or cell_type.lower() == 'end-end-cell'
                      or cell_type.lower() == 'end_end_cell') or cell_type.lower() == 'end end cell':
                    mid = perturbed_cell_node
                    left = perturbed_cell_node
                    right = perturbed_cell_node
                    beampipes = 'right'
                else:
                    mid = perturbed_cell_node
                    left = perturbed_cell_node
                    right = perturbed_cell_node
                    beampipes = 'both'

                enforce_Req_continuity(mid, left, right, cell_type)

                # perform checks on geometry
                fid = fr'{key}_Q{i}'

                # check if folder exists and skip if it does
                if os.path.exists(fr'{projectDir}\SimulationData\{analysis_folder}\{key}\{fid}'):
                    skip = True
                    info(["Skipped: ", fid, fr'{projectDir}\SimulationData\ABCI\{key}\{fid}'])

                if not skip:
                    solver.createFolder(fid, projectDir, subdir=sub_dir, opt=opt)
                    # it does not seem to make sense to perform uq on a multi cell by repeating the same perturbation
                    # to all multi cells at once. For multicells, the uq_multicell option is more suitable as it creates
                    # independent perturbations to all cells individually
                    solver.cavity(1, 1, mid, left, right, f_shift=0, bc=bc, beampipes=beampipes,
                                  fid=fid,
                                  sim_folder=analysis_folder, parentDir=parentDir, projectDir=projectDir,
                                  subdir=sub_dir)
                filename = fr'{projectDir}\SimulationData\{analysis_folder}\{key}\{fid}\monopole\qois.json'
                if os.path.exists(filename):
                    with open(filename, 'r') as file:
                        qois = json.load(file)
                    # extract objectives from tune_res
                    obj_result = list(
                        {key: qois[key] for key in objectives if
                         key in qois.keys()}.values())

                    # # sometimes some degenerate shapes are still generated and the solver returns zero
                    # # for the objective functions, such shapes are considered invalid
                    # for objr in obj_result:
                    #     if objr == 0:
                    #         # skip key
                    #         err = True
                    #         break
                    tab_val_f = obj_result
                    Ttab_val_f.append(tab_val_f)
                else:
                    err = True
            if err:
                break

            v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(np.array(Ttab_val_f)), weights)

            for i, o in enumerate(eigen_obj_list):
                result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
                result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

            with open(fr"{projectDir}\SimulationData\{analysis_folder}\{key}\uq.json", 'w') as file:
                file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))

            data_table = pd.DataFrame(Ttab_val_f, columns=eigen_obj_list)
            data_table.to_csv(fr'{projectDir}\SimulationData\{analysis_folder}\{key}\table.csv',
                              index=False, sep='\t', float_format='%.32f')

        if run_abci:
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
            progress_list = solver_args['progress_list']
            WG_M = solver_args['WG_M']
            # marker = solver_args['marker']

            # proc = solver_args['proc']
            sub_dir = fr'{key}'
            no_error = True
            for i in range(no_sims):
                skip = False
                for j, uq_var in enumerate(uq_vars):
                    uq_var_indx = VAR_TO_INDEX_DICT[uq_var]
                    perturbed_cell_node[uq_var_indx] = cell_node[uq_var_indx] * (1 + delta[j] * nodes[j, i])

                if (cell_type.lower() == 'mid cell' or cell_type.lower() == 'mid-cell' or
                        cell_type.lower() == 'mid_cell'):
                    cell_node = shape['IC']
                    mid = perturbed_cell_node
                    left = perturbed_cell_node
                    right = perturbed_cell_node
                    beampipes = 'none'
                elif (cell_type.lower() == 'mid-end cell' or cell_type.lower() == 'mid-end-cell'
                      or cell_type.lower() == 'mid_end_cell'):
                    mid = shape['IC']
                    left = shape['IC']
                    right = perturbed_cell_node
                    beampipes = 'right'
                elif (cell_type.lower() == 'end-end cell' or cell_type.lower() == 'end-end-cell'
                      or cell_type.lower() == 'end_end_cell') or cell_type.lower() == 'end end cell':
                    mid = perturbed_cell_node
                    left = perturbed_cell_node
                    right = perturbed_cell_node
                    beampipes = 'right'
                else:
                    mid = perturbed_cell_node
                    left = perturbed_cell_node
                    right = perturbed_cell_node
                    beampipes = 'both'

                enforce_Req_continuity(mid, left, right, cell_type)

                fid = fr'{key}_Q{i}'

                # check if folder exists and skip if it does
                if os.path.exists(fr'{projectDir}\SimulationData\ABCI\{key}\{fid}'):
                    skip = True

                if not skip:
                    solver.createFolder(fid, projectDir, subdir=sub_dir)
                    for wi in range(MROT):
                        solver.cavity(n_cells, n_modules, mid, left, right, fid=fid, MROT=wi,
                                      DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, beampipes=beampipes,
                                      bunch_length=bunch_length,
                                      MT=MT, NFS=NFS, UBT=UBT,
                                      parentDir=parentDir, projectDir=projectDir, WG_M='',
                                      marker='', sub_dir=sub_dir
                                      )

                # get objective function values
                abci_folder = fr'{projectDir}\SimulationData\ABCI\{key}'
                if os.path.exists(abci_folder):
                    obj_result = get_wakefield_objectives_value(fid, abci_obj_list, abci_folder)

                    tab_val_f = obj_result
                    if 'error' in obj_result:
                        no_error = False
                        error("Encountered an error in one or more objective function result.")
                        break
                    Ttab_val_f.append(tab_val_f)
                else:
                    no_error = False

            if no_error:
                v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(np.array(Ttab_val_f)), weights)
                for i, o in enumerate(abci_obj_list):
                    result_dict_abci[o]['expe'].append(v_expe_fobj[i])
                    result_dict_abci[o]['stdDev'].append(v_stdDev_fobj[i])

                with open(fr"{projectDir}\SimulationData\ABCI\{key}\uq.json", 'w') as file:
                    file.write(json.dumps(result_dict_abci, indent=4, separators=(',', ': ')))

            data_table.to_csv(fr'{projectDir}\SimulationData\ABCI\{key}\nodes.csv',
                              index=False, sep='\t', float_format='%.32f')

            data_table = pd.DataFrame(Ttab_val_f, columns=abci_obj_list)
            data_table.to_csv(fr'{projectDir}\SimulationData\ABCI\{key}\table.csv',
                              index=False, sep='\t', float_format='%.32f')


def uq_multicell(shape_shape, objectives, solver_dict, solver_args_dict, uq_config):
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

    parentDir = solver_args_dict['parentDir']
    projectDir = solver_args_dict['projectDir']
    cell_type = uq_config['cell type']
    analysis_folder = solver_args_dict['analysis folder']
    opt = solver_args_dict['optimisation']
    delta = uq_config['delta']
    method = uq_config['method']
    uq_vars = uq_config['variables']
    assert len(uq_vars) == len(delta), error('Ensure number of variables equal number of deltas')

    uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'

    err = False
    result_dict_eigen, result_dict_abci = {}, {}
    run_eigen, run_abci = False, False
    eigen_obj_list, abci_obj_list = [], []

    for o in objectives:
        if o in ["Req", "freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]",
                 "G [Ohm]", "Q []", 'kcc [%]', "ff [%]"]:
            result_dict_eigen[o] = {'expe': [], 'stdDev': []}
            run_eigen = True
            eigen_obj_list.append(o)

        if o.split(' ')[0] in ['ZL', 'ZT', 'k_loss', 'k_kick']:
            result_dict_abci[o] = {'expe': [], 'stdDev': []}
            run_abci = True
            abci_obj_list.append(o)
    n_cells = shape_space['IC'].shape[2]
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
    multicell_mid_vars = create_multicell_random_variables(n_cells, np.atleast_2d(np.array(shape_space['IC'])[:7]).T)

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = [np.array(shape_space['OC'])[:7], multicell_mid_vars, np.array(shape_space['OC_R'])[:7]]
    print(shape_space)

    rdim = len(np.array(shape_space['OC'])[:7]) + multicell_mid_vars.size + len(np.array(shape_space['OC_R'])[:7])

    degree = 1

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

        p_init_m = p_true[1] + nodes_[len(p_true[0]):len(p_true[0]) + p_true[1].size, i].reshape(
            np.shape(p_true[1]))

        p_init_er = p_true[2] + nodes_[len(p_true[0]) + p_true[1].size:, i]

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
                    solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                          n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                          beampipes=shape['BP'],
                                          parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                          mesh_args=mesh_args)
                else:
                    solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                            n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                            beampipes=shape['BP'],
                                            parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                            mesh_args=mesh_args)
            else:
                solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                        n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                        beampipes=shape['BP'],
                                        parentDir=parentDir, projectDir=projectDir, subdir=sub_dir,
                                        mesh_args=mesh_args)

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
        error(fr"There was a problem running UQ analysis for {key}")


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
        error(fr"There was a problem running UQ analysis for {key}")


def uq_ngsolve_parallel(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, pol, parentDir, projectDir,
                        mesh_args, select_solver='slans'):
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
            if (p + 1) * share < no_sims:
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


def uq_ngsolve_parallel_multicell(shape_space, objectives, solver_dict, solver_args_dict, uq_config):
    """

    Parameters
    ----------
    shape_space: dict
        Cavity geometry parameter space
    objectives: list | ndarray
        Array of objective functions
    solver_dict: dict
        Python dictionary of solver settings
    solver_args_dict: dict
        Python dictionary of solver arguments
    uq_config:
        Python dictionary of uncertainty quantification settings

    Returns
    -------

    """
    parentDir = solver_args_dict['parentDir']
    projectDir = solver_args_dict['projectDir']
    cell_type = uq_config['cell type']
    analysis_folder = solver_args_dict['analysis folder']
    opt = solver_args_dict['optimisation']
    delta = uq_config['delta']
    method = uq_config['method']
    uq_vars = uq_config['variables']
    n_cells = solver_args_dict['ngsolvemevp']['n_cells']
    assert len(uq_vars) == len(delta), error('Ensure number of variables equal number of deltas')

    for key, shape in shape_space.items():
        # n_cells = shape['IC'].shape[1] + 1
        uq_path = projectDir / fr'SimulationData\NGSolveMEVP\{key}'
        err = False
        result_dict_eigen, result_dict_abci = {}, {}
        run_eigen, run_abci = False, False
        eigen_obj_list, abci_obj_list = [], []

        for o in objectives:
            if o in ["Req", "freq [MHz]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "R/Q [Ohm]",
                     "G [Ohm]", "Q []", 'kcc [%]', "ff [%]"]:
                result_dict_eigen[o] = {'expe': [], 'stdDev': []}
                run_eigen = True
                eigen_obj_list.append(o)

            if o.split(' ')[0] in ['ZL', 'ZT', 'k_loss', 'k_kick']:
                result_dict_abci[o] = {'expe': [], 'stdDev': []}
                run_abci = True
                abci_obj_list.append(o)

        cav_var_list = ['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']
        midcell_var_dict = dict()
        for i1 in range(len(cav_var_list)):
            for i2 in range(n_cells):
                for i3 in range(2):
                    midcell_var_dict[f'{cav_var_list[i1]}_{i2}_m{i3}'] = [i1, i2, i3]

        # create random variables
        multicell_mid_vars = shape['IC']

        if n_cells == 1:
            # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
            p_true = [np.array(shape['OC'])[:7], np.array(shape['OC_R'])[:7]]
            rdim = len(np.array(shape['OC'])[:7]) + len(
                np.array(shape['OC_R'])[:7])  # How many variabels will be considered as random in our case 5
        else:
            # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
            p_true = [np.array(shape['OC'])[:7], multicell_mid_vars, np.array(shape['OC_R'])[:7]]
            rdim = len(np.array(shape['OC'])[:7]) + multicell_mid_vars.size + len(
                np.array(shape['OC_R'])[:7])  # How many variabels will be considered as random in our case 5
            # rdim = rdim - (n_cells*2 - 1)  # <- reduce dimension by making iris and equator radii to be equal

        # ic(rdim, multicell_mid_vars.size)
        # rdim = n_cells*3  # How many variables will be considered as random in our case 5
        degree = 1

        flag = 'stroud3'  # default
        if 'method' in uq_config.keys():
            if len(uq_config['method']) > 1:
                flag = uq_config['method'][1]

        if flag == 'stroud3':
            nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
            nodes_ = 2. * nodes_ - 1.
            # nodes_, weights_ = cn_leg_03_1(rdim)  # <- for some reason unknown this gives a less accurate answer.
            # the nodes are not the same as the custom function
        elif flag == 'stroud5':
            nodes_, weights_ = cn_leg_05_2(rdim)
        elif flag == 'cn_gauss':
            nodes_, weights_ = cn_gauss(rdim, 2)
        elif flag == 'lhc':
            sampler = qmc.LatinHypercube(d=rdim)
            _ = sampler.reset()
            nsamp = 3000
            sample = sampler.random(n=nsamp)
            # ic(qmc.discrepancy(sample))
            l_bounds = -np.ones(rdim)
            u_bounds = np.ones(rdim)
            sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

            nodes_, weights_ = sample_scaled.T, np.ones((nsamp, 1))
        elif flag == 'load_from_file':
            nodes_ = pd.read_csv(fr'C:\Users\sosoho\DakotaProjects\Cavity\C3795_lhs\sim_result_table.dat',
                                 sep='\s+').iloc[:, 2:-2]
            nodes_ = nodes_.to_numpy().T
            weights_ = np.ones((nodes_.shape[1], 1))
        else:
            # defaults to Stroud3
            nodes_, weights_, bpoly_ = quad_stroud3(rdim, degree)
            nodes_ = 2. * nodes_ - 1.

        # save nodes
        data_table = pd.DataFrame(nodes_.T)
        data_table.to_csv(uq_path / 'nodes.csv', index=False, sep='\t', float_format='%.32f')

        #  mean value of geometrical parameters
        no_parm, no_sims = np.shape(nodes_)

        sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run

        proc_count = uq_config['processes']
        if proc_count > no_sims:
            proc_count = no_sims

        share = round(no_sims / proc_count)
        jobs = []
        for p in range(proc_count):
            # try:
            end_already = False
            if p != proc_count - 1:
                if (p + 1) * share < no_sims:
                    proc_keys_list = np.arange(p * share, p * share + share)
                else:
                    proc_keys_list = np.arange(p * share, no_sims)
                    end_already = True

            if p == proc_count - 1 and not end_already:
                proc_keys_list = np.arange(p * share, no_sims)

            processor_nodes = nodes_[:, proc_keys_list]
            processor_weights = weights_[proc_keys_list]

            service = mp.Process(target=uq_multicell_sequential, args=(
                n_cells, 1, shape, objectives, n_cells+1, 0, 33, 'monopole', parentDir,
                projectDir, sub_dir, key, uq_path,
                proc_keys_list, processor_nodes, processor_weights, p, p_true))

            service.start()
            jobs.append(service)

        for job in jobs:
            job.join()

        # combine results from processes
        qois_result_dict = {}
        Ttab_val_f = []
        keys = []
        for i1 in range(no_sims):
            if i1 == 0:
                df = pd.read_csv(uq_path / fr'table_{i1}.csv', sep='\t', engine='python')
            else:
                try:
                    df = pd.concat([df, pd.read_csv(uq_path / fr'table_{i1}.csv', sep='\t', engine='python')])
                except:
                    pass

        df.to_csv(uq_path / 'table.csv', index=False, sep='\t', float_format='%.32f')
        df.to_excel(uq_path / 'table.xlsx', index=False)

        # data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
        # ic(data_table)
        # data_table.index = keys
        # data_table = data_table.sort_index()
        # ic(data_table)
        # ic(keys)
        #
        # data_table.to_csv(uq_path / fr'table.csv', index=False, sep='\t')
        # data_table.to_excel(uq_path / fr'table.xlsx', index=False)

        Ttab_val_f = df.to_numpy()
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(Ttab_val_f, weights_)

        # append results to dict
        for i, o in enumerate(eigen_obj_list):
            result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
            result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

        with open(uq_path / fr'uq.json', 'w') as file:
            file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))


def uq_multicell_sequential(n_cells, n_modules, shape, qois, n_modes, f_shift, bc, pol, parentDir, projectDir, sub_dir,
                            key, uq_path, proc_keys_list, processor_nodes, processor_weights,
                            proc_num, p_true):
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
            p_init_el = p_true[0] + processor_nodes[0:len(p_true[0]), i1 - min(proc_keys_list)]
            p_init_er = p_true[1] + processor_nodes[len(p_true[0]):, i1 - min(proc_keys_list)]
            par_mid = p_init_el
        else:
            proc_node = processor_nodes[:, i1 - min(proc_keys_list)]
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

            # p_init_el = processor_nodes[0:len(p_true[0]), i1 - min(proc_keys_list)]
            # p_init_m = processor_nodes[len(p_true[0]):len(p_true[0]) + p_true[1].size,
            #            i1 - min(proc_keys_list)].reshape(np.shape(p_true[1])[::-1]).T
            # p_init_er = processor_nodes[len(p_true[0]) + p_true[1].size:, i1 - min(proc_keys_list)]

            p_init_el = p_true[0] + processor_nodes[0:len(p_true[0]), i1-min(proc_keys_list)]
            p_init_m = p_true[1] + processor_nodes[len(p_true[0]):len(p_true[0])+p_true[1].size, i1-min(proc_keys_list)].reshape(np.shape(p_true[1]))
            p_init_er = p_true[2] + processor_nodes[len(p_true[0])+p_true[1].size:, i1-min(proc_keys_list)]
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
            info(f'processor {proc_num} skipped ', fid, 'Result already exists.')

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
                    solver.cavity_flattop(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                          n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                          beampipes=shape['BP'],
                                          parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
                else:
                    solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                            n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                            beampipes=shape['BP'],
                                            parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
            else:
                solver.cavity_multicell(n_cells, n_modules, par_mid, par_end_l, par_end_r,
                                        n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, pol=pol,
                                        beampipes=shape['BP'],
                                        parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

        filename = uq_path / f'{fid}/monopole/qois.json'
        # print(filename)
        if os.path.exists(filename):
            qois_result_dict = dict()

            with open(filename) as json_file:
                qois_result_dict.update(json.load(json_file))

            qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
            # sometimes some degenerate shapes are still generated and the solver returns zero
            # for the objective functions, such shapes are considered invalid
            # for objr in qois_result:
            #     if objr == 0:
            #         # skip key
            #         err = True
            #         break

            tab_val_f = qois_result
            Ttab_val_f.append(tab_val_f)
        else:
            err = True

    # # add original point
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)
    # save table
    # print(f'\tDone with proc {proc_num} ', time.time()-start)
    data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    data_table.to_csv(uq_path / fr'table_{proc_num}.csv', index=False, sep='\t', float_format='%.32f')

    # print(np.atleast_2d(Ttab_val_f), processor_weights)


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
        if os.path.exists(fr'{projectDir}\SimulationData\ABCI\{key}'):
            skip = True

        start_time = time.time()
        if not skip:
            # run abci code
            # run both polarizations if MROT == 2
            if 'OC_R' in list(shape.keys()):
                OC_R = 'OC_R'
            else:
                OC_R = 'OC'

            if MROT == 2:
                for m in range(2):
                    abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape[OC_R],
                                     fid=key, MROT=m, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                     DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir,
                                     projectDir=projectDir,
                                     WG_M='', marker='')
            else:
                abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape[OC_R],
                                 fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                 DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir,
                                 WG_M='', marker='')

        print(f'Cavity {key}. Time: {time.time() - start_time}')

        # update progress
        progress_list.append((progress + 1) / total_no_of_shapes)


def get_objectives_value(d, obj, norm_length, n_cells):
    Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
    Freq = d['FREQUENCY'][n_cells - 1]
    E_stored = d['STORED ENERGY'][n_cells - 1]
    # Rsh = d['SHUNT IMPEDANCE'][n_cells - 1]  # MOhm
    Q = d['QUALITY FACTOR'][n_cells - 1]
    Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
    Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
    # Vacc = dict['ACCELERATION'][n_cells - 1]
    # Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells - 1]  # MV/m
    r_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm
    G = 0.00948 * Q * (Freq / 1300)
    GR_Q = G * 2 * r_Q

    Vacc = np.sqrt(
        2 * r_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6  # factor of 2, remember circuit and accelerator definition
    # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
    Eacc = Vacc / (n_cells * norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
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
    # k_loss_transverse = []
    # k_loss_longitudinal = []
    # k_loss_M0 = []
    # key_list = []

    # create list to hold Z
    Zmax_mon_list = []
    Zmax_dip_list = []
    xmax_mon_list = []
    xmax_dip_list = []
    processed_keys = []

    # def calc_k_loss():
    #     print(f"Processing for Cavity {key}")
    #     abci_data_long = ABCIData(abci_data_dir, key, 0)
    #     abci_data_trans = ABCIData(abci_data_dir, key, 1)
    #
    #     # trans
    #     x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
    #     k_loss_trans = abci_data_trans.loss_factor['Transverse']
    #
    #     if math.isnan(k_loss_trans):
    #         print_(f"Encountered an exception: Check shape {key}")
    #         return [0, 0, 0]
    #
    #     # long
    #     x, y, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
    #     abci_data_long.get_data('Loss Factor Spectrum Integrated up to F')
    #
    #     k_M0 = abci_data_long.y_peaks[0]
    #     k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
    #     k_loss_HOM = k_loss_long - k_M0
    #
    #     # append only after successful run
    #     k_loss_M0.append(k_M0)
    #     k_loss_longitudinal.append(k_loss_HOM)
    #     k_loss_transverse.append(k_loss_trans)
    #
    #     return [k_loss_M0, k_loss_longitudinal, k_loss_transverse]

    def get_Zmax_L(mon_interval=None):
        # print("2a")
        if mon_interval is None:
            mon_interval = [0.0, 2e10]

        print(f"Processing for Cavity {key}")
        try:
            abci_data_mon = ABCIData(abci_data_dir, f"{key}", 0)

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
            for n, z_bound in enumerate(mon_interval):
                # get mask
                msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

                if len(yp_mon[msk_mon]) != 0:
                    Zmax_mon = max(yp_mon[msk_mon])

                    Zmax_mon_list[n].append(Zmax_mon)
                elif len(yp_mon) != 0:
                    Zmax_mon_list[n].append(0)
                else:
                    return ['Error']

            processed_keys.append(key)
        except KeyError as e:
            return ['Exception occurred', e]

        # print("2g", Zmax_mon_list)

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

            for n, z_bound in enumerate(dip_interval):
                # get mask
                msk_dip = [(z_bound[0] < x < z_bound[1]) for x in xp_dip]

                if len(yp_dip[msk_dip]) != 0:
                    Zmax_dip = max(yp_dip[msk_dip])

                    Zmax_dip_list[n].append(Zmax_dip)
                elif len(yp_dip) != 0:
                    Zmax_dip_list[n].append(0)
                else:
                    return ['Error']

            processed_keys.append(key)
        except KeyError as e:
            return ['Exception occurred', e]

        return Zmax_dip_list

    # def eval_all(mon_interval, dip_interval):
    #     print(f"Processing for Cavity {key}")
    #     abci_data_long = ABCIData(abci_data_dir, f"{key}_", 0)
    #     abci_data_trans = ABCIData(abci_data_dir, f"{key}_", 1)
    #
    #     # get longitudinal and transverse impedance plot data
    #     xr_mon, yr_mon, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
    #     xi_mon, yi_mon, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')
    #
    #     xr_dip, yr_dip, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
    #     xi_dip, yi_dip, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')
    #
    #     # loss factors
    #     # trans
    #     k_loss_trans = abci_data_trans.loss_factor['Transverse']
    #
    #     if math.isnan(k_loss_trans):
    #         print_(f"Encountered an exception: Check shape {key}")
    #         return 0
    #
    #     # long
    #     abci_data_long.get_data('Loss Factor Spectrum Integrated upto F')
    #
    #     k_M0 = abci_data_long.y_peaks[0]
    #     k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
    #     k_loss_HOM = k_loss_long - k_M0
    #
    #     # calculate magnitude
    #     ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]
    #     ymag_dip = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_dip, yi_dip)]
    #
    #     # get peaks
    #     peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
    #     xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]
    #
    #     peaks_dip, _ = sps.find_peaks(ymag_dip, height=0)
    #     xp_dip, yp_dip = np.array(xr_dip)[peaks_dip], np.array(ymag_dip)[peaks_dip]
    #
    #     for n, z_bound in enumerate(mon_interval):
    #         # get mask
    #         msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]
    #
    #         if len(yp_mon[msk_mon]) != 0:
    #             Zmax_mon = max(yp_mon[msk_mon])
    #             xmax_mon = xp_mon[np.where(yp_mon == Zmax_mon)][0]
    #
    #             Zmax_mon_list[n].append(Zmax_mon)
    #             xmax_mon_list[n].append(xmax_mon)
    #         elif len(yp_mon) != 0:
    #             Zmax_mon_list[n].append(0.0)
    #             xmax_mon_list[n].append(0.0)
    #         else:
    #             continue
    #
    #     for n, z_bound in enumerate(dip_interval):
    #         # get mask
    #         msk_dip = [(z_bound[0] < x < z_bound[1]) for x in xp_dip]
    #
    #         if len(yp_dip[msk_dip]) != 0:
    #             Zmax_dip = max(yp_dip[msk_dip])
    #             xmax_dip = xp_dip[np.where(yp_dip == Zmax_dip)][0]
    #
    #             Zmax_dip_list[n].append(Zmax_dip)
    #             xmax_dip_list[n].append(xmax_dip)
    #         elif len(yp_dip) != 0:
    #             Zmax_dip_list[n].append(0.0)
    #             xmax_dip_list[n].append(0.0)
    #         else:
    #             continue
    #
    #     # append only after successful run
    #
    #     k_loss_M0.append(k_M0)
    #     k_loss_longitudinal.append(k_loss_HOM)
    #     k_loss_transverse.append(k_loss_trans)

    ZL, ZT = [], []
    freq_range_ZL, freq_range_ZT = [], []
    # print("here here here")
    for i, o in enumerate(obj):
        if o[1].split(' ')[0] == 'ZL':
            freq_range_ZL.append(o[2])
        elif o[1].split(' ')[0] == 'ZT':
            freq_range_ZT.append(o[2])

        elif o[1] == "k_loss":
            pass
        elif o[1] == "k_kick":
            pass

    # print("about to evaluate ZL", freq_range_ZL)
    if freq_range_ZL:
        for _ in range(len(freq_range_ZL)):
            Zmax_mon_list.append([])
            xmax_mon_list.append([])

        ZL = get_Zmax_L(freq_range_ZL)

    if freq_range_ZT:
        for _ in range(len(freq_range_ZT)):
            Zmax_dip_list.append([])
            xmax_dip_list.append([])

        ZT = get_Zmax_T(freq_range_ZT)

    ZL, ZT = np.array(ZL).T, np.array(ZT).T

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


def axis_data_coords_sys_transform(axis_obj_in, xin, yin, inverse=False):
    """ inverse = False : Axis => Data
                    = True  : Data => Axis
        """
    if axis_obj_in.get_yscale() == 'log':
        xlim = axis_obj_in.get_xlim()
        ylim = axis_obj_in.get_ylim()

        x_delta = xlim[1] - xlim[0]

        if not inverse:
            x_out = xlim[0] + xin * x_delta
            y_out = ylim[0] ** (1 - yin) * ylim[1] ** yin
        else:
            x_delta2 = xin - xlim[0]
            x_out = x_delta2 / x_delta
            y_out = np.log(yin / ylim[0]) / np.log(ylim[1] / ylim[0])

    else:
        xlim = axis_obj_in.get_xlim()
        ylim = axis_obj_in.get_ylim()

        x_delta = xlim[1] - xlim[0]
        y_delta = ylim[1] - ylim[0]

        if not inverse:
            x_out = xlim[0] + xin * x_delta
            y_out = ylim[0] + yin * y_delta
        else:
            x_delta2 = xin - xlim[0]
            y_delta2 = yin - ylim[0]
            x_out = x_delta2 / x_delta
            y_out = y_delta2 / y_delta

    return x_out, y_out


def _get_nodes_and_weights(uq_config, rdim, degree):
    method = uq_config['method']
    uq_vars = uq_config['variables']

    if method[1].lower() == 'stroud3':
        nodes, weights, bpoly = quad_stroud3(rdim, degree)
        nodes = 2. * nodes - 1.
        # nodes, weights = cn_leg_03_1(rdim)
    elif method[1].lower() == 'stroud5':
        nodes, weights = cn_leg_05_2(rdim)
    elif method[1].lower() == 'gaussian':
        nodes, weights = cn_gauss(rdim, 2)
    elif method[1].lower() == 'lhs':
        sampler = qmc.LatinHypercube(d=rdim)
        _ = sampler.reset()
        nsamp = uq_config['integration'][2]
        sample = sampler.random(n=nsamp)

        l_bounds = [-1 for _ in range(len(uq_vars))]
        u_bounds = [1 for _ in range(len(uq_vars))]
        sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

        nodes, weights = sample_scaled.T, np.ones((nsamp, 1))
    elif method[0].lower() == 'from file':
        if len(method) == 2:
            nodes = pd.read_csv(method[1], sep='\s+').iloc[:, method[1]]
        else:
            nodes = pd.read_csv(method[1], sep='\s+')

        nodes = nodes.to_numpy().T
        weights = np.ones((nodes.shape[1], 1))
    else:
        # issue warning
        warning('Integration method not recognised. Defaulting to Stroud3 quadrature rule!')
        nodes, weights, bpoly = quad_stroud3(rdim, degree)
        nodes = 2. * nodes - 1.

    return nodes, weights


def add_text(ax, text, box, xy=(0.5, 0.5), xycoords='data', xytext=None, textcoords='data',
             size=14, rotation=0, arrowprops=None):
    """

    Parameters
    ----------
    text: str
        Matplotlib annotation text
    box
    xy: tuple
        Coordinates of annotation text
    xycoords: str {data, axis}
        Coordinate system reference
    xytext
    textcoords
    size
    rotation: float
        Annotation text rotation
    arrowprops

    Returns
    -------

    """
    if text.strip("") == "":
        return

    # add text
    if xytext:
        bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
        annotext = ax.annotate(text, xy=xy, xycoords=xycoords,
                               xytext=xytext, textcoords=textcoords, bbox=bbox_props, fontsize=size,
                               rotation=rotation, arrowprops=arrowprops)
    else:
        if box == "None":
            annotext = ax.annotate(text, xy=xy, xycoords=xycoords, fontsize=size,
                                   rotation=rotation, arrowprops=arrowprops)
        else:
            bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
            annotext = ax.annotate(text, xy=xy, xycoords=xycoords, bbox=bbox_props, fontsize=size,
                                   rotation=rotation, arrowprops=arrowprops)

    ax.get_figure().canvas.draw_idle()
    ax.get_figure().canvas.flush_events()
    return annotext


def show_welcome():
    # display(Image(filename='.\cav_images\G0_C3_CO.png'))
    message = (
        "<style>p{background: -webkit-linear-gradient(#eee, #333);-webkit-background-clip: text; webkit-text-fill-color: transparent; }</style> <p><b>CAV-SIM-2D</b> loaded successfully!</p>")
    display(HTML(message))


# Call the function when the module is imported
show_welcome()
