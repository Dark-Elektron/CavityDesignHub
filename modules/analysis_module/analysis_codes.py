import ctypes
import gc
import json
import math
import os
import time
from threading import Thread
import numpy as np
from scipy.optimize import fsolve
from termcolor import colored
from simulation_codes.SLANS.slans_tune import Tune
# from objective_functions import *
from utils import misc_functions as f
from modules.tune_module.tuners.tuner import Tuner

file_color = 'green'

DEBUG = True
# DEBUG = False
def print_(*arg):
    if DEBUG: print(colored(f'\t\t{arg}', file_color))


class Analysis:
    def __init__(self, win=None):
        pass
        # self.win = win
        # if win:
        #     self.pytune = Tuner(win)
        #
        #     self.win = win
        #     # self.abci = win.abci_geom
        #     # self.slans = win.slans_geom
        #     self.ui = win.ui
        #     # self.MROT = self.ui.cb_Polarization_ABCI.currentIndex()
        #     # self.MROT = self.ui.cb_Polarization_SLANS.currentIndex()
        #
        #
        #     # initialise thread state
        #     self.thread_state = 'run'
        # else:

    def GSSEC(self, pseudo_shape_space, bc, parentDir, projectDir, filename, resume="No",
              proc=0, tuner='SLANS', tune_variable='Req', iter_set=None, cell_type='Mid Cell'):

        # tuner
        self.pytune = Tuner()

        if tuner == 'SLANS':
            self.slans_tune = Tune(parentDir, projectDir)

        start = time.time()
        # get Ri_space from cutoff frequency
        # Ri_space = self.calculate_Ri(A, B, a, b, Ri)
        self.population = {}
        self.dip_dist = {}
        self.freq = {}

        print_("Getting filename from window")
        progress = 0
        # if self.win:
        # # set progress bar
        # self.ui.w_Progress_Bar.show()
        # self.ui.progressBar.setMinimum(0)
        # print("here", len(list(pseudo_shape_space)))
        # self.ui.progressBar.setMaximum(len(list(pseudo_shape_space)))
        # self.ui.progressBar.setValue(0)

        # check for already processed shapes
        existing_keys = []

        print("Resume", resume)
        if resume == "Yes":
            # check if value set is already written. This is to enable continuation in case of break in program
            if os.path.exists(fr'{projectDir}/Cavities/{filename}'):
                self.population = json.load(open(fr'{projectDir}/Cavities/{filename}', 'r'))

                existing_keys = list(self.population.keys())
                print_(f'Existing keys: {existing_keys}')

        start_time = time.time()

        for key, pseudo_shape in pseudo_shape_space.items():
            A_i, B_i, a_i, b_i, Ri_i, L_i, Req, _ = pseudo_shape['IC'] # Req here is none but required since the shape space consists of all variables
            A_o, B_o, a_o, b_o, Ri_o, L_o, Req, _ = pseudo_shape['OC'] # Req here is none but required since the shape space consists of all variables
            beampipes = pseudo_shape['BP']
            freq = pseudo_shape['FREQ']

            # new mid cell and end cell with initial Req guess
            inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0]
            outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0]

            # Update progressbar
            progress += 1
            # self.ui.progressBar.setValue(progress)
            # print_("here 2:".format(key, last_saved_key))

            # edit to check for key later
            if key not in existing_keys:
                if tuner == 'SLANS':
                    if tune_variable == 'Req':
                        # Tune cell to get Req
                        Req, freq, alpha, h, e = self.slans_tune.mid_cell_tune(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, freq, proc=proc)
                    else:
                        L, freq, alpha, h, e = self.slans_tune.end_cell_tune(inner_cell, outer_cell, freq, proc=proc)

                        if cell_type == 'Mid Cell':
                            L_i, L_o = L, L
                        else:
                            L_o = L

                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha]
                    outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha]
                else:
                    if tune_variable == 'Req':
                        print('PyTune Req')
                        Req, freq = self.pytune.tune("Req", inner_cell, outer_cell, freq, beampipes, bc, parentDir, projectDir, iter_set=iter_set, proc=proc)
                        # round
                        Req, freq = round(Req, 4), round(freq, 2)
                    else:
                        print('PyTune L')
                        L, freq = self.pytune.tune("L", inner_cell, outer_cell, freq, beampipes, bc, parentDir, projectDir, iter_set=iter_set, proc=proc)

                        # round
                        L, freq = round(L, 2), round(freq, 2)

                        if cell_type == 'Mid Cell':
                            L_i, L_o = L, L
                        else:
                            L_o = L

                    alpha_i = self.calculate_alpha(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0)
                    alpha_o = self.calculate_alpha(A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0)

                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha_i]
                    outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha_o]

                print_(f'Done Tuning Cavity {key}')

                if cell_type == 'Mid Cell':
                    self.population[key] = {"IC": inner_cell, "OC":outer_cell, "BP": 'none', 'FREQ': freq}
                else:
                    self.population[key] = {"IC": inner_cell, "OC":outer_cell, "BP": 'both', 'FREQ': freq}

            print_("SAVING DICTIONARIES", f"shape_space{proc}.json")
            with open(fr"{projectDir}\Cavities\shape_space{proc}.json", 'w') as file:
                file.write(json.dumps(self.population, indent=4, separators=(',', ': ')))
            print_("DONE SAVING")

            print_("Time run:: ", time.time() - start_time)
            start_time = time.time()

        # print_("Extracted Data::", self.population)
        end = time.time()

        runtime = end - start
        print_(f'Proc {proc} runtime: {runtime}')

    def write_dict(self, filename, dict):
        with open(filename, 'w') as file:
            file.write(json.dumps(dict, indent=4, separators=(',', ': ')))

    def set_data_folder(self, fid=0, MROT=-1):
        self.data_dir = self.ui.le_Data_Directory.text()
        print_(self.data_dir)
        if MROT == -1:
            pass
        else:
            if self.data_dir == "":
                self.data_dir = 'Data\ABCI\Cavity{}\Cavity_MROT_{}.top'.format(fid, MROT)
            elif ":" in self.data_dir:
                self.data_dir = '{}\ABCI\Cavity{}\Cavity_MROT_{}.top'.format(self.data_dir, fid, MROT)
            else:
                self.data_dir = os.path.join(os.getcwd(), '{}\ABCI\Cavity{}\Cavity_MROT_{}.top'.format(self.data_dir, fid, MROT))

        print_("\t\t Data folder --> ", self.data_dir)

    def calculate_Ri(self, A, B, a, b, Ri):
        Ri_list = []

        # get frequency of first (TE11) dipole mode
        freq = self.get_freq(A, B, a, b, Ri) #GHz

        # calculate Ri from frequency
        c = 299792458 # m/s
        pp_TE_11 = 1.841 #p'(mn)
        Ri_0 = (c*pp_TE_11)/(2*np.pi * freq * 1e9) * 1e3 # convert from m to mm

        # increment Ri and store in a list
        Ri_list = [math.ceil(Ri_0) + 1 + x for x in range(5)]
        # freq_l = [(c*pp_TE_11)/(2*np.pi * Ri * 1e9) * 1e3 for Ri in Ri_list]
        print_(Ri_list)

        return Ri_list

    def get_freq(self, A, B, a, b, Ri):
        freq = 1000 # random value holder
        l = [A, B, a, b, Ri]
        freq_pop = json.load(open('Extracted Data\dipole_mode_freq.json', 'r'))

        fid = f.get_point_fid(l)
        print_(fid, freq_pop[fid][0])
        freq = freq_pop[fid][0]

        return freq

    def calculate_alpha(self, A, B, a, b, Ri, L, Req, L_bp):

        data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
                [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(self.ellipse_tangent,
                                np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                                args=data)
        m = (y2 - y1) / (x2 - x1)
        alpha = 180 - np.arctan(m) * 180 / np.pi
        return alpha

    def ellipse_tangent(self, z, *data):
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return f1, f2, f3, f4

    def save_to_file(self):
        with open('dictionary.txt', 'w') as f:
            for key, val in self.population.items():
                f.write("{}: {} ==> {}\n".format(key, val, self.dip_dist[key]))

    def set_thread_state(self, thread_state):
        print_("Thread state-> ", thread_state)
        self.thread_state = thread_state
