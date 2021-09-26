import ast
import json
import os
import subprocess
import sys
import time

# from mpi4py import MPI
# from SLANS.slans_geometry import SLANSGeometry
from geometry import Geometry
# import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.animation as anime
# from threading import Thread
# import scipy.signal as sps
# import multiprocessing as mp
# import matplotlib.pyplot as plt
import pandas as pd

if sys.platform.startswith("win"):
    # Don't display the Windows GPF dialog if the invoked program dies.
    # See comp.os.ms-windows.programmer.win32
    #  How to suppress crash notification dialog?, Jan 14,2004 -
    #     Raymond Chen's response [1]

    import ctypes

    SEM_NOGPFAULTERRORBOX = 0x0002  # From MSDN
    ctypes.windll.kernel32.SetErrorMode(SEM_NOGPFAULTERRORBOX)
    CREATE_NO_WINDOW = 0x08000000  # From Windows API
    subprocess_flags = CREATE_NO_WINDOW
else:
    subprocess_flags = 0


class TunePy(Geometry):
    def __init__(self, win=None):
        if win:
            super().__init__(win)
            self.win = win
            self.ui = win.ui
            self.slans_geom = win.slans_geom

        os.chdir(r'C:\\Users\\sosoho\\Dropbox\\2D_Codes\\ABCI_software\\Python_ABCI')
        # initialise variables
        self.start_time = time.time()
        self.count = 1
        self.proc_count = 12
        self.population_end_cell = {}
        # self.cid = cid
        #
        # path = os.path.join(os.getcwd(), r"ConvergencePlot\convergence_plot_{}.json".format(self.cid))
        # if os.path.exists(path):
        #     if response.strip().lower() == 'y' or response.strip() == '':
        #         self.convergence_dict = json.load(open(path, 'r'))
        #         # print("initial convergence dict: ", self.convergence_dict)
        # else:
        self.convergence_dict = {}

        # dictionary to plot to see the graph
        self.plot_dict = {}

        # save max to a list in the case of finite iterations to pick the max later
        self.mx_dict = {}

    def window_search(self, par_in, par_end, l_bound, r_bound, key):
        print("\t\tTUNE_PY_REQ::ITERATION NUMBER:: ", self.count)

        # Get variables from gui I'd have to keep track of the bin somehow
        par_mid = par_in  # self.mid_cell
        par_l_end = par_end  # self.left_end_cell
        par_r_end = par_end  # self.right_end_cell
        print("\t\tTUNE_PY_REQ::CELL PARAMS: ", par_mid, par_l_end, par_r_end)

        # par_mid = [52.0, 52.0, 28.5, 28.5, 169.476, 55.0, 93.5]
        # par_l_end = [50, 50, 15, 15, 169.476, 82, 65]
        # par_r_end = [50, 50, 15, 15, 169.476, 82, 65]

        # run mpi to get the freq_list
        shift = (r_bound - l_bound) / (2 * self.proc_count)

        # get Req_list
        self.Req_list = np.linspace(l_bound + shift, r_bound - shift, self.proc_count)

        print(os.getcwd())
        # run slans
        command = ["mpiexec", "-np", "{}".format(self.proc_count), "python",
                   "C:\\Users\sosoho\Dropbox\\2D_Codes\ABCI_software\Python_ABCI\mpi_script_Req.py",
                   "{}".format(par_mid),
                   "{}".format(par_l_end), "{}".format(par_r_end), "{}".format(l_bound), "{}".format(r_bound)]

        sp = subprocess.run(command)

        try:
            with open('freq_list.txt', 'r') as f:
                self.freq_list = ast.literal_eval(u'{}'.format(f.readline()))
        except:
            self.freq_list = [0]
            with open('freq_list.txt', 'r') as f:
                self.freq_list = ast.literal_eval(u'{}'.format(f.readline()))

        # get freq and Req
        Req, freq = self.get_Req()

        print("\t\tTUNE_PY_REQ::Lists ", self.Req_list, self.freq_list)
        print("\t\tTUNE_PY_REQ::Req, freq:: {}mm, {}%".format(Req, freq))

        # # update dictionary
        # for i, L in enumerate(self.Req_list):
        #     self.plot_dict[L] = self.freq_list[i]

        lb, rb = Req - shift, Req + shift
        print("\t\tTUNE_PY_REQ::LEFT BOUND: {}, RIGHT BOUND: {}".format(lb, rb))

        # check progress of simulation
        error = abs(801.58 - freq) / 801.58
        print("\t\t\tTUNE_PY_REQ:: ERROR -> ", error)

        if error <= 1e-5 or self.count == 5:
            print("\t\tTUNE_PY_REQ:: Done:: Req, Freq:: ", Req, freq)
            # reset count
            self.count = 1
            return Req, freq

        else:
            print("\t\t\tTUNE_PY_REQ:: Req, Freq:: ", Req, freq)
            self.count += 1

            # update par_end
            par_mid[-3] = Req
            par_l_end = par_mid
            par_r_end = par_mid
            return self.window_search(par_mid, par_l_end, lb, rb, key)

    def get_Req(self):
        # peaks, _ = sps.find_peaks(self.freq_list)
        print('\t\tgetting max')
        min_dist_list = [abs(801.58 - x) for x in self.freq_list]
        min_dist = min(min_dist_list)
        Req = 0
        freq = 0
        for k, mn_dt in enumerate(min_dist_list):
            if mn_dt == min_dist:
                Req = self.Req_list[k]
                freq = self.freq_list[k]

        print("\t\t\tTUNE_PY_REQ:: Get_Req:: ", Req, freq)
        return Req, freq

    def end_routine(self, key):
        # get max from dict
        freq_mx = max(list(self.mx_dict.keys()))
        L_mx = self.mx_dict[freq_mx]

        print("\tFinal result:: ", freq_mx, L_mx)
        self.end_time = time.time()
        print("\tTime:: ", self.end_time - self.start_time)

        # # save dict
        # self.convergence_dict[key] = self.plot_dict
        # with open(r"ConvergencePlot\convergence_plot_C{}.json".format(self.cid), 'w') as file:
        #     file.write(json.dumps(self.convergence_dict, indent=4, separators=(',', ': ')))

        # reset all
        self.count = 1
        self.start_time = 0
        self.end_time = 0
        self.plot_dict = {}
        self.mx_dict = {}

        # # plot dict
        # x = list(self.plot_dict.keys())
        # y = list(self.plot_dict.values())
        # plt.scatter(x, y)
        # plt.show()

        return L_mx, freq_mx


if __name__ == '__main__':
    g = TunePy()
    print(g.window_search(90))
