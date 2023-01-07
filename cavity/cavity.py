import time

import numpy as np
import multiprocessing as mp
from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from utils.shared_functions import *

slans_geom = SLANSGeometry()

import psutil

m0 = 9.1093879e-31
q0 = 1.6021773e-19
c0 = 2.99792458e8
mu0 = 4 * np.pi * 1e-7
eps0 = 8.85418782e-12


class Cavity:
    def __init__(self, n_cells, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none', name='cavity'):
        self.name = name
        self.n_cells = n_cells
        self.mid_cell = mid_cell
        self.end_cell_left = end_cell_left
        self.end_cell_right = end_cell_right
        self.beampipe = beampipe
        self.no_of_modules = 1

        # get geometric parameters
        self.shape_space = {
            "IC": self.mid_cell,
            "OC": self.end_cell_left,
            "OC_R": self.end_cell_right,
            "BP": beampipe
        }

        self.freq = None
        if not end_cell_left:
            end_cell_left = mid_cell

        if not end_cell_right:
            if not end_cell_left:
                end_cell_right = mid_cell
            else:
                end_cell_right = end_cell_left

        self.A, self.B, self.a, self.b, self.Ri, self.L, self.Req = mid_cell
        self.A_el, self.B_el, self.a_el, self.b_el, self.Ri_el, self.L_el, self.Req_el = end_cell_left
        self.A_er, self.B_er, self.a_er, self.b_er, self.Ri_er, self.L_er, self.Req_er = end_cell_right

    def run_eigenmode(self, solver='SLANS', freq_shift=0, boundary_conds=None, no_of_procs=1):
        # required
        # boundary conditions
        #
        if boundary_conds is None:
            boundary_conds = [3, 3]

        if solver == 'SLANS':
            self.run_slans()

    def run_wakefield(self, solver='ABCI'):
        if solver == 'ABCI':
            self.run_abci()

    def calc_op_freq(self):
        if not self.freq:
            self.freq = (c0/4*self.L)

    @staticmethod
    def run_sequential(name, n_cells, n_modules, shape, n_modes, f_shift, bc, parentDir, projectDir,
                       progress_list, sub_dir='', UQ=False):
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

        print(f'Done with Cavity {name}. Time: {time.time() - start_time}')

    def run_abci(self):
        pass

    def uq(self):
        pass
