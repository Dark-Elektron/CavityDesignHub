import json
import os
import shutil
import subprocess
from math import floor
from pathlib import Path

import pandas as pd
from scipy.signal import find_peaks
from termcolor import colored
from analysis_modules.eigenmode.SLANS.geometry_manual import Geometry
from analysis_modules.eigenmode.SLANS.slans_code import SLANS, SLANS_Multicell, SLANS_Multicell_full

from utils.file_reader import FileReader
from utils.shared_functions import *

fr = FileReader()

file_color = 'green'
DEBUG = False


def print_(*arg):
    if DEBUG:
        print(colored(f'\t\t\t\t{arg}', file_color))


class SLANSGeometry(Geometry):
    def __init__(self, win=None):
        if win:
            super().__init__(win)

            self.ui = win.ui
            self.slansUI = win.slansUI
            self.slans = None

    def cavity(self, no_of_cells=1, no_of_modules=1, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
               fid=None, bc=33, f_shift='default', beta=1, n_modes=None, beampipes='None',
               parentDir=None, projectDir=None, subdir='', expansion=None, expansion_r=None):
        """
        Write geometry file and run eigenmode analysis with SLANS

        Parameters
        ----------
        expansion
        no_of_cells: int
            Number of cells
        no_of_modules: int
            Number of modules
        mid_cells_par: list, array like
            Mid cell geometric parameters -> [A, B, a, b, Ri, L, Req, alpha]
        l_end_cell_par: list, array like
            Left end cell geometric parameters -> [A_el, B_el, a_el, b_el, Ri_el, L_el, Req, alpha_el]
        r_end_cell_par: list, array like
            Right end cell geometric parameters -> [A_er, B_er, a_er, b_er, Ri_er, L_er, Req, alpha_er]
        fid: int, str
            File id
        bc: int
            Boundary condition -> 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        f_shift: float
            Eigenvalue frequency shift
        beta: int, float
            Velocity ratio :math: `\\beta = \frac{v}{c}`
        n_modes: int
            Number of modes
        beampipes: {"left", "right", "both", "none"}
            Specify if beam pipe is on one or both ends or at no end at all
        parentDir: str
            Parent directory
        projectDir: str
            Project directory
        subdir: str
            Sub directory to save simulation results to

        Returns
        -------

        """

        # this checks whether input is from gui or from the optimisation
        if mid_cells_par is not None:
            self.set_geom_parameters(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par,
                                     beampipes, expansion=expansion, expansion_r=expansion_r)
        else:
            self.set_geom_parameters(no_of_cells)

        self.slans = SLANS(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell,
                           self.right_beam_pipe, self.Jxy_all, self.Jxy_all_bp)

        n = no_of_cells  # Number of cells
        axi_sym = 2  # 1: flat, 2: axis-symmetric
        unit = 3  # 1:m, 2:cm, 3:mm, 4:mkm
        name_index = 1
        sc = 1
        end_type = 1  # if end_type = 1 the end HALF cell is changed for tuning.
        # If end_type = 2 the WHOLE end cell is changed for tuning
        end_L = 1  # if end_L = 1 the type of end cell is type a (without iris)
        # if end_L = 2 the type of end cell is type b
        end_R = 1

        if expansion is not None:
            end_L = 2
        if expansion_r is not None:
            end_R = 2

        # # Beam pipe length
        # if end_L == 1:
        #     self.Rbp_L = self.ri_L
        #
        # if end_R == 1:
        #     self.Rbp_R = self.ri_R

        # print_(self.WG_L, self.WG_R)

        # Ellipse conjugate points x,y
        zr12_L, alpha_L = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r
        zr12_R, alpha_R = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r
        zr12_M, alpha_M = self.slans.rz_conjug('mid')  # zr12_R first column is z , second column is r

        if end_L == 2:
            zr12_BPL, alpha_BPL = self.slans.rz_conjug('expansion')  # zr12_R first column is z , second column is r

        if end_R == 2:
            zr12_BPR, alpha_BPR = self.slans.rz_conjug('expansion_r')  # zr12_R first column is z , second column is r

        # Set boundary conditions
        BC_Left = floor(bc / 10)  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        BC_Right = bc % 10  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal

        filename = f'cavity_{bc}'

        # change save directory
        if subdir == '':
            run_save_directory = projectDir / fr'SimulationData/SLANS/{fid}'
        else:
            run_save_directory = projectDir / fr'SimulationData/SLANS/{subdir}/{fid}'

        # Write SLANS Geometry
        print(self.Jxy, n, self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2))
        with open(Path(fr'{run_save_directory}/{filename}.geo'), 'w') as f:
            # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
            f.write('8 {:.0f} {:.0f} 2 {}\n'.format(
                self.Jxy * n + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) +
                (1 if self.WG_L > 0 else 0) * self.WG_mesh + (1 if self.WG_R > 0 else 0) * self.WG_mesh,
                self.Jxy, unit))

            f.write('10 0 0 0 0 0 0 0 0\n')
            f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.ri_L, self.Jy0, BC_Left))

            if end_L == 2:
                f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.Rbp_L, self.Jxy_all_bp[5] + self.Jxy_all_bp[6] +
                                                                  self.Jxy_all_bp[7], BC_Left))

            if self.WG_L > 0:
                if end_L == 2:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L - self.x_L, self.Rbp_L, self.WG_mesh))
                else:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L, self.Rbp_L, self.WG_mesh))

            # n == 1
            if n == 1:
                if self.Req_L != self.Req_R:
                    print_('The equator radius of left and right cell are not equal')

                # if exist('L_M') != 1:
                #     L_M = []

                if end_L == 2:
                    self.slans.slans_bp_L(n, zr12_BPL, self.WG_L, f)

                self.slans.slans_n1_L(n, zr12_L, self.WG_L, f)
                self.slans.slans_n1_R(n, zr12_R, self.WG_L, f)

                if end_R == 2:
                    self.slans.slans_bp_R(n, zr12_BPR, self.WG_L, f)

                if self.WG_R > 0:
                    if end_R == 2:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L + self.WG_R + self.L_L +
                                                                        self.L_R, self.Rbp_R,
                                                                        self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L + self.WG_R + self.L_L +
                                                                        self.L_R, self.Rbp_R,
                                                                        self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.WG_L + self.WG_R + self.L_L +
                                                                         self.L_R, self.ri_R,
                                                                         -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] +
                                                                           self.Jxy_all_bp[7]), BC_Right))

                f.write(
                    '1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(self.WG_L + self.WG_R + self.L_L +
                                                              self.L_R, -self.Jy0, BC_Right))

                f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(-(self.Jxy * n + self.Jxy_bp * (
                        (1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) + (
                                                                1 if self.WG_L > 0 else 0) * self.WG_mesh + (
                                                                1 if self.WG_R > 0 else 0) * self.WG_mesh)))
                f.write('0 0 0 0 0 0 0 0 0')

            # n>1
            if n > 1:
                if end_L == 2:
                    self.slans.slans_bp_L(n, zr12_BPL, self.WG_L, f)

                self.slans.slans_n1_L(n, zr12_L, self.WG_L, f)

                for i in range(1, n):
                    self.slans.slans_M(n, zr12_M, self.WG_L, f, i, end_type)

                self.slans.slans_n1_R(n, zr12_R, self.WG_L, f)

                if end_R == 2:
                    self.slans.slans_bp_R(n, zr12_BPR, self.WG_L, f)

                if self.WG_R > 0:
                    if end_R == 2:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R +
                            2 * (n - 1) * self.L_M, self.Rbp_R, self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R +
                            2 * (n - 1) * self.L_M, self.Rbp_R, self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, self.ri_R,
                        -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] + self.Jxy_all_bp[7]), BC_Right))

                f.write('1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(
                    self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, -self.Jy0, BC_Right))

                # gradual mesh decrease
                if self.WG_R > 0:
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L +
                                                                 self.L_R + 2 * (n - 1) * self.L_M,
                                                                 -((1 if self.WG_R > 0 else 0) * self.WG_mesh)))

                f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (n - 1) * self.L_M - self.L_M,
                                                             -(self.Jxy * 1)))

                for i in range(n - 1, 1, -1):
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L +
                                                                 self.L_L + 2 * (i - 1) * self.L_M - self.L_M,
                                                                 -(self.Jxy * 1)))

                f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L, -(self.Jxy * 1)))

                if self.WG_L > 0:
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(0, -((1 if self.WG_L > 0 else 0) * self.WG_mesh + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2))))

                # # direct mesh decrease
                # f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(
                #     -(self.Jxy*n+self.Jxy_bp*((1 if end_R == 2 else 0)/2+(1 if end_L == 2 else 0)/2) +
                #       (1 if self.WG_L > 0 else 0)*self.WG_mesh+(1 if self.WG_R > 0 else 0)*self.WG_mesh)))

                f.write('0 0 0 0 0 0 0 0 0')

        # Slans run
        genmesh_path = parentDir / fr'exe/SLANS_exe/genmesh2.exe'
        filepath = Path(fr'{run_save_directory}/{filename}')

        # folder for exe to write to
        cwd = run_save_directory

        # the next two lines suppress pop up windows from the slans codes
        # the slans codes, however, still disrupts windows operation, sadly. This is the case even for the slans tuner
        if os.name == 'nt':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            kwargs = {"startupinfo": startupinfo}
        else:
            kwargs = {}

        print(str(filepath))
        subprocess.call([genmesh_path, str(filepath), '-b'], cwd=cwd, **kwargs)
        path = run_save_directory

        beta, f_shift, n_modes = 1, 0, no_of_cells + 1

        self.write_dtr(path, filename, beta, f_shift, n_modes)

        slansc_path = parentDir / fr'exe/SLANS_exe/slansc'
        slansm_path = parentDir / fr'exe/SLANS_exe/slansm'
        slanss_path = parentDir / fr'exe/SLANS_exe/slanss'
        slansre_path = parentDir / fr'exe/SLANS_exe/slansre'

        # print(cwd)
        # check if corresponding file exists at before the executable is called
        if os.path.exists(projectDir / fr'SimulationData/SLANS/{fid}/{filename}.geo'):
            subprocess.call([slansc_path, str(filepath), '-b'], cwd=cwd, **kwargs)  # settings, number of modes, etc
            ic("Done with slansc")

            if os.path.exists(projectDir / fr'SimulationData/SLANS/{fid}/{filename}.gem'):
                subprocess.call([slansm_path, str(filepath), '-b'], cwd=cwd, **kwargs)
                ic("Done with slansm")

                if os.path.exists(projectDir / fr'SimulationData/SLANS/{fid}/aslans.mtx') \
                        and os.path.exists(projectDir / fr'SimulationData/SLANS/{fid}/bslans.mtx'):
                    subprocess.call([slanss_path, str(filepath), '-b'], cwd=cwd, **kwargs)
                    ic("Done with slanss")

                    if os.path.exists(projectDir / fr'SimulationData/SLANS/{fid}/{filename}.res'):
                        subprocess.call([slansre_path, str(filepath), '-b'], cwd=cwd, **kwargs)
                        ic("Done with slansre")

        # save json file
        shape = {'IC': update_alpha(mid_cells_par),
                 'OC': update_alpha(l_end_cell_par),
                 'OC_R': update_alpha(r_end_cell_par)}

        with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))
        try:
            filename = Path(fr'{run_save_directory}/cavity_33.svl')
            d = fr.svl_reader(filename)

            Req = d['CAVITY RADIUS'][no_of_cells - 1] * 10  # convert to mm
            Freq = d['FREQUENCY'][no_of_cells - 1]
            E_stored = d['STORED ENERGY'][no_of_cells - 1]
            Rsh = d['SHUNT IMPEDANCE'][no_of_cells - 1]  # MOhm
            Q = d['QUALITY FACTOR'][no_of_cells - 1]
            Epk = d['MAXIMUM ELEC. FIELD'][no_of_cells - 1]  # MV/m
            Hpk = d['MAXIMUM MAG. FIELD'][no_of_cells - 1]  # A/m
            # Vacc = dict['ACCELERATION'][no_of_cells - 1]
            Eavg = d['AVERAGE E.FIELD ON AXIS'][no_of_cells - 1]  # MV/m
            r_Q = d['EFFECTIVE IMPEDANCE'][no_of_cells - 1]  # Ohm
            G = 0.00948 * Q * np.sqrt(Freq / 1300)
            GR_Q = G * 2 * r_Q

            Vacc = np.sqrt(
                2 * r_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6  # factor of 2, circuit and accelerator definition
            # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
            norm_length = 2 * mid_cells_par[5]
            Eacc = Vacc / (
                    no_of_cells * norm_length * 1e-3)  # for 1 cell factor of 2, circuit and accelerator definition
            Epk_Eacc = Epk / Eacc
            Bpk = (Hpk * 4 * np.pi * 1e-7) * 1e3
            Bpk_Eacc = Bpk / Eacc

            # cel to cell coupling factor
            f_diff = d['FREQUENCY'][no_of_cells - 1] - d['FREQUENCY'][0]
            f_add = (d['FREQUENCY'][no_of_cells - 1] + d['FREQUENCY'][0])
            kcc = 2 * f_diff / f_add * 100

            try:
                # field flatness
                ax_field = self.get_axis_field_data(run_save_directory, no_of_cells)
                # get max in each cell
                peaks, _ = find_peaks(ax_field['y_abs'])
                E_abs_peaks = ax_field['y_abs'][peaks]
                ff = min(E_abs_peaks) / max(E_abs_peaks) * 100
            except FileNotFoundError:
                ff = 0

            d = {
                "Req [mm]": Req,
                "Normalization Length [mm]": norm_length,
                "N Cells": no_of_cells,
                "freq [MHz]": Freq,
                "Q []": Q,
                "E [MV/m]": E_stored,
                "Vacc [MV]": Vacc,
                "Eacc [MV/m]": Eacc,
                "Epk [MV/m]": Epk,
                "Hpk [A/m]": Hpk,
                "Bpk [mT]": Bpk,
                "kcc [%]": kcc,
                "ff [%]": ff,
                "Rsh [Ohm]": Rsh,
                "R/Q [Ohm]": 2 * r_Q,
                "Epk/Eacc []": Epk_Eacc,
                "Bpk/Eacc [mT/MV/m]": Bpk_Eacc,
                "G [Ohm]": G,
                "GR/Q [Ohm^2]": GR_Q
            }

            with open(Path(fr'{run_save_directory}/qois.json'), "w") as f:
                json.dump(d, f, indent=4, separators=(',', ': '))
        except FileNotFoundError as e:
            print("Simulation failed", e)

    def cavity_multicell(self, no_of_modules=1, cells_par=None, fid=None, bc=33,
                         f_shift='default', beta=1, n_modes=None, beampipes="None",
                         parentDir=None, projectDir=None, subdir=''):

        no_of_cells = len(cells_par)

        if not n_modes:
            n_modes = no_of_cells + 1  # intentional because the accuracy of the last mode is always low
        # perform geometric shapes
        # Req of adjacent cells must match # for now all Reqs should be equal
        a = np.array(cells_par)
        assert np.all(a[:, 6] == a[:, 6])

        cells_par = pd.DataFrame(cells_par, columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha"])

        left_cell = cells_par.loc[0, :].tolist()
        mid_cells = cells_par.loc[1:no_of_cells - 2, :].to_numpy().T
        right_cell = cells_par.loc[no_of_cells - 1, :].tolist()

        # this checks whether input is from gui or from the optimisation
        print(type(cells_par))
        if isinstance(cells_par, pd.DataFrame):
            self.set_geom_parameters_multicell(no_of_cells, mid_cells, left_cell, right_cell, beampipes)
        else:
            self.set_geom_parameters_multicell(no_of_cells)

        self.slans = SLANS_Multicell(self.left_beam_pipe, left_cell, mid_cells, right_cell,
                                     self.right_beam_pipe, self.Jxy_all, self.Jxy_all_bp)

        n = no_of_cells  # Number of cells
        axi_sym = 2  # 1: flat, 2: axis-symmetric
        unit = 3  # 1:m, 2:cm, 3:mm, 4:mkm
        name_index = 1
        sc = 1
        end_type = 1  # if end_type = 1 the end HALF cell is changed for tuning.
        # If end_type = 2 the WHOLE end cell is changed for tuning
        end_L = 1  # if end_L = 1 the type of end cell is type a (without iris)
        # if end_L = 2 the type of end cell is type b
        end_R = 1

        # Beam pipe length
        if end_L == 1:
            self.Rbp_L = self.ri_L

        if end_R == 1:
            self.Rbp_R = self.ri_R

        # print_(self.WG_L, self.WG_R)

        # Ellipse conjugate points x,y
        zr12_L, alpha_L = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r
        zr12_R, alpha_R = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r

        if end_L == 2:
            zr12_BPL, alpha_BPL = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r

        if end_R == 2:
            zr12_BPR, alpha_BPR = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r

        # Set boundary conditions
        BC_Left = floor(bc / 10)  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        BC_Right = bc % 10  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal

        filename = f'cavity_{bc}'

        # change save directory
        if subdir == '':
            run_save_directory = Path(fr'{parentDir}/{projectDir}/{fid}')
        else:
            run_save_directory = projectDir / fr'SimulationData/SLANS/{subdir}/{fid}'

        # Write Slans Geometry
        with open(fr'{run_save_directory}/{filename}.geo', 'w') as f:
            # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
            f.write('8 {:.0f} {:.0f} 2 {}\n'.format(
                self.Jxy * n + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) +
                (1 if self.WG_L > 0 else 0) * self.WG_mesh + (1 if self.WG_R > 0 else 0) * self.WG_mesh,
                self.Jxy, unit))

            f.write('10 0 0 0 0 0 0 0 0\n')
            f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.ri_L, self.Jy0, BC_Left))

            if end_L == 2:
                f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.Rbp_L, self.Jxy_all_bp[5] + self.Jxy_all_bp[6] +
                                                                  self.Jxy_all_bp[7], BC_Left))

            if self.WG_L > 0:
                if end_L == 2:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L - self.x_L, self.Rbp_L, self.WG_mesh))
                else:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L, self.Rbp_L, self.WG_mesh))

            if n > 1:
                if end_L == 2:
                    self.slans.slans_bp_L(n, zr12_BPL, self.WG_L, f)

                self.slans.slans_n1_L(zr12_L, self.WG_L, f)

                for i in range(1, n - 1):
                    zr12_M, alpha_M = self.slans.rz_conjug('mid',
                                                           i - 1)  # zr12_R first column is z , second column is r
                    self.slans.slans_M(n, zr12_M, self.WG_L, f, i - 1, end_type)

                self.slans.slans_n1_R(n, zr12_R, self.WG_L, f)

                if end_R == 2:
                    self.slans.slans_bp_R(n, zr12_BPR, self.WG_L, f)

                if self.WG_R > 0:
                    if end_R == 2:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R +
                            sum(self.L_M), self.Rbp_R, self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R +
                            sum(self.L_M), self.Rbp_R, self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + sum(self.L_M), self.ri_R,
                        -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] + self.Jxy_all_bp[7]), BC_Right))

                f.write('1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(
                    self.WG_L + self.WG_R + self.L_L + self.L_R + sum(self.L_M), -self.Jy0, BC_Right))

                # # gradual mesh decrease
                # if self.WG_R > 0:
                #     f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L +
                #                   self.L_R + 2 * (n - 1) * sum(self.L_M),
                #                                                  -((1 if self.WG_R > 0 else 0) * self.WG_mesh)))
                #
                # f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (n - 1) * sum(self.L_M) - sum(self.L_M),
                #                                              -(self.Jxy * 1)))
                #
                # for i in range(n - 1, 1, -1):
                #     f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L +
                #     self.L_L + 2 * (i - 1) * sum(self.L_M) - sum(self.L_M), -(self.Jxy * 1)))
                #
                # f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L, -(self.Jxy * 1)))
                #
                # if self.WG_L > 0:
                #     f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(0, -((1 if self.WG_L > 0 else 0) * self.WG_mesh)))

                # direct mesh decrease
                f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(
                    -(self.Jxy * n + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) +
                      (1 if self.WG_L > 0 else 0) * self.WG_mesh + (1 if self.WG_R > 0 else 0) * self.WG_mesh)))

                f.write('0 0 0 0 0 0 0 0 0')

        # Slans run
        genmesh_path = parentDir / fr'exe/SLANS_exe/genmesh2.exe'
        filepath = Path(fr'{run_save_directory}/{filename}')

        # folder for exe to write to
        cwd = run_save_directory

        # the next two lines suppress pop up windows from the slans codes
        # the slans codes, however, still disrupts windows operation, sadly. This is the case even for the slans tuner
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW

        subprocess.call([genmesh_path, filepath, '-b'], cwd=cwd, startupinfo=startupinfo)
        path = run_save_directory

        # if f_shift == 'default':
        #     # parameters delete later
        #     if self.ui.le_Beta.text() or self.ui.le_Freq_Shift.text() or self.ui.sb_No_Of_Modes.value():
        #         beta, f_shift, n_modes = float(self.ui.le_Beta.text()), float(self.ui.le_Freq_Shift.text()),
        #         self.ui.sb_No_Of_Modes.value()
        #         # print(beta, f_shift, n_modes)
        #     else:
        #         beta, f_shift, n_modes = 1, 0, 1

        beta, f_shift, n_modes = 1, 0, no_of_cells + 1

        self.write_dtr(path, filename, beta, f_shift, n_modes)

        slansc_path = parentDir / fr'exe/SLANS_exe/slansc'
        slansm_path = parentDir / fr'exe/SLANS_exe/slansm'
        slanss_path = parentDir / fr'exe/SLANS_exe/slanss'
        slansre_path = parentDir / fr'exe/SLANS_exe/slansre'

        # print(cwd)
        subprocess.call([slansc_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)
        subprocess.call([slansm_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)
        subprocess.call([slanss_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)
        subprocess.call([slansre_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)

        # save json file
        shape = {'IC': mid_cells_par,
                 'OC': l_end_cell_par,
                 'OC_R': r_end_cell_par}

        with open(fr"{run_save_directory}/geometric_parameters.json", 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))

        filename = fr'{run_save_directory}/cavity_33.svl'
        d = fr.svl_reader(filename)

        Req = d['CAVITY RADIUS'][no_of_cells - 1] * 10  # convert to mm
        Freq = d['FREQUENCY'][no_of_cells - 1]
        E_stored = d['STORED ENERGY'][no_of_cells - 1]
        Rsh = d['SHUNT IMPEDANCE'][no_of_cells - 1]  # MOhm
        Q = d['QUALITY FACTOR'][no_of_cells - 1]
        Epk = d['MAXIMUM ELEC. FIELD'][no_of_cells - 1]  # MV/m
        Hpk = d['MAXIMUM MAG. FIELD'][no_of_cells - 1]  # A/m
        # Vacc = dict['ACCELERATION'][no_of_cells - 1]
        Eavg = d['AVERAGE E.FIELD ON AXIS'][no_of_cells - 1]  # MV/m
        r_Q = d['EFFECTIVE IMPEDANCE'][no_of_cells - 1]  # Ohm
        G = 0.00948 * Q * np.sqrt(Freq / 1300)
        GR_Q = G * 2 * r_Q

        Vacc = np.sqrt(
            2 * r_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6  # factor of 2, circuit and accelerator definition
        # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
        norm_length = 2 * mid_cells_par[5]
        Eacc = Vacc / (
                no_of_cells * norm_length * 1e-3)  # for 1 cell factor of 2, circuit and accelerator definition
        Epk_Eacc = Epk / Eacc
        Bpk = (Hpk * 4 * np.pi * 1e-7) * 1e3
        Bpk_Eacc = Bpk / Eacc

        # cel to cell coupling factor
        f_diff = d['FREQUENCY'][no_of_cells - 1] - d['FREQUENCY'][0]
        f_add = (d['FREQUENCY'][no_of_cells - 1] + d['FREQUENCY'][0])
        kcc = 2 * f_diff / f_add * 100

        # field flatness
        ax_field = self.get_axis_field_data(run_save_directory, no_of_cells)
        # get max in each cell
        peaks, _ = find_peaks(ax_field['y_abs'])
        E_abs_peaks = ax_field['y_abs'][peaks]
        ff = min(E_abs_peaks) / max(E_abs_peaks) * 100

        d = {
            "Req [mm]": Req,
            "Normalization Length [mm]": norm_length,
            "freq [MHz]": Freq,
            "Q []": Q,
            "E [MV/m]": E_stored,
            "Vacc [MV]": Vacc,
            "Eacc [MV/m]": Eacc,
            "Epk [MV/m]": Epk,
            "Hpk [A/m]": Hpk,
            "Bpk [mT]": Bpk,
            "kcc [%]": kcc,
            "ff [%]": ff,
            "R/Q [Ohm]": 2 * r_Q,
            "Epk/Eacc []": Epk_Eacc,
            "Bpk/Eacc [mT/MV/m]": Bpk_Eacc,
            "G [Ohm]": G,
            "GR/Q [Ohm^2]": GR_Q
        }

        with open(Path(fr'{run_save_directory}/qois.json'), "w") as f:
            json.dump(d, f, indent=4, separators=(',', ': '))

    def cavity_multicell_full(self, no_of_modules=1, cells_par=None, fid=None, bc=33,
                              f_shift='default', beta=1, n_modes=None, beampipes="None",
                              parentDir=None, projectDir=None, subdir=''):

        no_of_cells = int(len(cells_par) / 2)
        no_of_half_cells = int(len(cells_par))

        if not n_modes:
            n_modes = int(no_of_half_cells / 2) + 1  # intentional because the accuracy of the last mode is always low

        # perform geometric checks
        # length of cells_par must be even
        assert no_of_half_cells % 2 == 0
        # Req of adjacent cells must match
        Req_pairs = [[2 * i, 2 * i + 1] for i in range(int(no_of_cells))]
        Ri_pairs = [[2 * i - 1, 2 * i] for i in range(1, int(no_of_cells))]

        a = np.array(cells_par)

        Req_check = np.all(
            [a[:, 6][Req_pairs][i][0] == a[:, 6][Req_pairs][i][1] for i in range(len(a[:, 6][Req_pairs]))])
        Ri_check = np.all([a[:, 4][Ri_pairs][i][0] == a[:, 4][Ri_pairs][i][1] for i in range(len(a[:, 4][Ri_pairs]))])
        # print(a)
        # print(Ri_pairs)
        assert Ri_check
        assert Req_check

        cells_par = pd.DataFrame(cells_par, columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha"])

        left_cell = cells_par.loc[0, :].tolist()
        mid_cells = cells_par.loc[1:no_of_half_cells - 2, :].to_numpy().T
        right_cell = cells_par.loc[no_of_half_cells - 1, :].tolist()

        # this checks whether input is from gui or from the optimisation
        if isinstance(cells_par, pd.DataFrame):
            self.set_geom_parameters_multicell(no_of_half_cells, mid_cells, left_cell, right_cell, beampipes)
        else:
            self.set_geom_parameters_multicell(no_of_half_cells)

        self.slans = SLANS_Multicell_full(self.left_beam_pipe, left_cell, mid_cells, right_cell,
                                          self.right_beam_pipe, self.Jxy_all, self.Jxy_all_bp)

        n = no_of_half_cells  # Number of cells
        axi_sym = 2  # 1: flat, 2: axis-symmetric
        unit = 3  # 1:m, 2:cm, 3:mm, 4:mkm
        name_index = 1
        sc = 1
        end_type = 1  # if end_type = 1 the end HALF cell is changed for tuning.
        # If end_type = 2 the WHOLE end cell is changed for tuning
        end_L = 1  # if end_L = 1 the type of end cell is type a (without iris)
        # if end_L = 2 the type of end cell is type b
        end_R = 1

        # Beam pipe length
        if end_L == 1:
            self.Rbp_L = self.ri_L

        if end_R == 1:
            self.Rbp_R = self.ri_R

        # print_(self.WG_L, self.WG_R)

        # Ellipse conjugate points x,y
        zr12_L, alpha_L = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r
        zr12_R, alpha_R = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r

        if end_L == 2:
            zr12_BPL, alpha_BPL = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r

        if end_R == 2:
            zr12_BPR, alpha_BPR = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r

        # Set boundary conditions
        BC_Left = floor(bc / 10)  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        BC_Right = bc % 10  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal

        filename = f'cavity_{bc}'

        # change save directory
        if subdir == '':
            run_save_directory = Path(fr'{parentDir}/{projectDir}/{fid}')
        else:
            print("it's here", fid)
            run_save_directory = parentDir / fr'SimulationData/SLANS/{subdir}/{fid}'

        # Write Slans Geometry
        with open(fr'{run_save_directory}/{filename}.geo', 'w') as f:
            # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
            f.write('8 {:.0f} {:.0f} 2 {}\n'.format(
                self.Jxy * n / 2 + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) +
                (1 if self.WG_L > 0 else 0) * self.WG_mesh + (1 if self.WG_R > 0 else 0) * self.WG_mesh,
                self.Jxy, unit))

            f.write('10 0 0 0 0 0 0 0 0\n')
            f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.ri_L, self.Jy0, BC_Left))

            if end_L == 2:
                f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.Rbp_L, self.Jxy_all_bp[5] + self.Jxy_all_bp[6] +
                                                                  self.Jxy_all_bp[7], BC_Left))

            if self.WG_L > 0:
                if end_L == 2:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L - self.x_L, self.Rbp_L, self.WG_mesh))
                else:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L, self.Rbp_L, self.WG_mesh))

            if n > 2:
                if end_L == 2:
                    self.slans.slans_bp_L(n, zr12_BPL, self.WG_L, f)

                self.slans.slans_n1_L(zr12_L, self.WG_L, f)

                for i in range(0, n - 2):
                    zr12_M, alpha_M = self.slans.rz_conjug('mid', i)  # zr12_R first column is z , second column is r
                    self.slans.slans_M(n, zr12_M, self.WG_L, f, i, end_type)

                self.slans.slans_n1_R(n, zr12_R, self.WG_R, f)

                if end_R == 2:
                    self.slans.slans_bp_R(n, zr12_BPR, self.WG_R, f)

                if self.WG_R > 0:
                    if end_R == 2:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R +
                            sum(self.L_M), self.Rbp_R, self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R +
                            sum(self.L_M), self.Rbp_R, self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + sum(self.L_M), self.ri_R,
                        -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] + self.Jxy_all_bp[7]), BC_Right))

                f.write('1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(
                    self.WG_L + self.WG_R + self.L_L + self.L_R + sum(self.L_M), -self.Jy0, BC_Right))

                # # gradual mesh decrease
                # if self.WG_R > 0:
                #     f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L +
                #                   self.L_R + 2 * (n - 1) * sum(self.L_M),
                #                                                  -((1 if self.WG_R > 0 else 0) * self.WG_mesh)))
                #
                # f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (n - 1) * sum(self.L_M) - sum(self.L_M),
                #                                              -(self.Jxy * 1)))
                #
                # for i in range(n - 1, 1, -1):
                #     f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L +
                #     self.L_L + 2 * (i - 1) * sum(self.L_M) - sum(self.L_M), -(self.Jxy * 1)))
                #
                # f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L, -(self.Jxy * 1)))
                #
                # if self.WG_L > 0:
                #     f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(0, -((1 if self.WG_L > 0 else 0) * self.WG_mesh)))

                # direct mesh decrease
                f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(
                    -(self.Jxy * n / 2 + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) +
                      (1 if self.WG_L > 0 else 0) * self.WG_mesh + (1 if self.WG_R > 0 else 0) * self.WG_mesh)))

                f.write('0 0 0 0 0 0 0 0 0')

        # Slans run
        genmesh_path = parentDir / fr'exe/SLANS_exe/genmesh2.exe'
        filepath = fr'{run_save_directory}/{filename}'

        # folder for exe to write to
        cwd = run_save_directory

        # the next two lines suppress pop up windows from the slans codes
        # the slans codes, however, still disrupts windows operation, sadly. This is the case even for the slans tuner
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        # print(genmesh_path, filepath, cwd)
        subprocess.call([genmesh_path, filepath, '-b'], cwd=cwd, startupinfo=startupinfo)
        path = run_save_directory

        # if f_shift == 'default':
        #     # parameters delete later
        #     if self.ui.le_Beta.text() or self.ui.le_Freq_Shift.text() or self.ui.sb_No_Of_Modes.value():
        #         beta, f_shift, n_modes = float(self.ui.le_Beta.text()), float(self.ui.le_Freq_Shift.text()),
        #         self.ui.sb_No_Of_Modes.value()
        #         # print(beta, f_shift, n_modes)
        #     else:
        #         beta, f_shift, n_modes = 1, 0, 1

        beta, f_shift, n_modes = 1, 0, no_of_cells + 1

        self.write_dtr(path, filename, beta, f_shift, n_modes)

        slansc_path = parentDir / fr'exe/SLANS_exe/slansc'
        slansm_path = parentDir / fr'exe/SLANS_exe/slansm'
        slanss_path = parentDir / fr'exe/SLANS_exe/slanss'
        slansre_path = parentDir / fr'exe/SLANS_exe/slansre'

        # print(cwd)
        subprocess.call([slansc_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)
        subprocess.call([slansm_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)
        subprocess.call([slanss_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)
        subprocess.call([slansre_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)

        # save json file
        shape = {'GEOM': cells_par.values.tolist(),
                 "BP": beampipes,
                 "FREQ": 0}

        with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))

        filename = Path(fr'{run_save_directory}/cavity_33.svl')
        d = fr.svl_reader(filename)

        Req = d['CAVITY RADIUS'][no_of_cells - 1] * 10  # convert to mm
        Freq = d['FREQUENCY'][no_of_cells - 1]
        E_stored = d['STORED ENERGY'][no_of_cells - 1]
        Rsh = d['SHUNT IMPEDANCE'][no_of_cells - 1]  # MOhm
        Q = d['QUALITY FACTOR'][no_of_cells - 1]
        Epk = d['MAXIMUM ELEC. FIELD'][no_of_cells - 1]  # MV/m
        Hpk = d['MAXIMUM MAG. FIELD'][no_of_cells - 1]  # A/m
        # Vacc = dict['ACCELERATION'][no_of_cells - 1]
        Eavg = d['AVERAGE E.FIELD ON AXIS'][no_of_cells - 1]  # MV/m
        r_Q = d['EFFECTIVE IMPEDANCE'][no_of_cells - 1]  # Ohm
        G = 0.00948 * Q * np.sqrt(Freq / 1300)
        GR_Q = G * 2 * r_Q

        Vacc = np.sqrt(
            2 * r_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6  # factor of 2, circuit and accelerator definition
        # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
        cells_par_list = cells_par.values.tolist()
        norm_length = 2 * cells_par_list[1][5]
        Eacc = Vacc / (
                no_of_cells * norm_length * 1e-3)  # for 1 cell factor of 2, circuit and accelerator definition
        Epk_Eacc = Epk / Eacc
        Bpk = (Hpk * 4 * np.pi * 1e-7) * 1e3
        Bpk_Eacc = Bpk / Eacc

        # cel to cell coupling factor
        f_diff = d['FREQUENCY'][no_of_cells - 1] - d['FREQUENCY'][0]
        f_add = (d['FREQUENCY'][no_of_cells - 1] + d['FREQUENCY'][0])
        kcc = 2 * f_diff / f_add * 100

        # field flatness
        ax_field = self.get_axis_field_data(run_save_directory, no_of_cells)
        # get max in each cell
        peaks, _ = find_peaks(ax_field['y_abs'])
        E_abs_peaks = ax_field['y_abs'][peaks]
        ff = min(E_abs_peaks) / max(E_abs_peaks) * 100

        d = {
            "Req [mm]": Req,
            "Normalization Length [mm]": norm_length,
            "freq [MHz]": Freq,
            "Q []": Q,
            "E [MV/m]": E_stored,
            "Vacc [MV]": Vacc,
            "Eacc [MV/m]": Eacc,
            "Epk [MV/m]": Epk,
            "Hpk [A/m]": Hpk,
            "Bpk [mT]": Bpk,
            "kcc [%]": kcc,
            "ff [%]": ff,
            "R/Q [Ohm]": 2 * r_Q,
            "Epk/Eacc []": Epk_Eacc,
            "Bpk/Eacc [mT/MV/m]": Bpk_Eacc,
            "G [Ohm]": G,
            "GR/Q [Ohm^2]": GR_Q
        }

        with open(Path(fr'{run_save_directory}/qois.json'), "w") as f:
            json.dump(d, f, indent=4, separators=(',', ': '))

    @staticmethod
    def write_dtr(path, filename, beta, f_shift, n_modes):
        with open(r"{}/{}.dtr".format(path, filename), 'w') as f:
            f.write(':          Date:02/04/16 \n')
            f.write('{:g} :number of iterative modes 1-10\n'.format(n_modes))
            f.write('{:g} :number of search modes\n'.format(n_modes - 1))
            f.write('9.99999997E-007 :convergence accuracy\n')
            f.write('50 :maximum number of iterations\n')
            f.write('0 :continue iterations or not 1,0\n')
            f.write(' {:g}. :initial frequency shift MHz\n'.format(f_shift))
            f.write('1 :wave type 1-E, 2-H\n')
            f.write(' 1 :struct. 1-cav,2-per.str,3-w.guid.,4-l.-hom.\n')
            f.write('0 :symmetry yes or not 1,0\n')
            f.write(' 1 :number of met.surfaces, then:sign and sigma\n')
            f.write('5  1.\n')
            f.write('0 : number of mark volumes,then:sign,EPS,MU,TGE,TGM\n')
            f.write('{:g} : beta (v/c)\n'.format(beta))

    @staticmethod
    def createFolder(fid, projectDir, subdir=''):
        # change save directory
        path = projectDir / fr'SimulationData/SLANS/{fid}'
        if subdir == '':
            pass
        else:
            new_path = projectDir / fr'SimulationData/SLANS/{subdir}/{fid}'
            if os.path.exists(new_path):
                path = new_path
            else:
                if not os.path.exists(projectDir / fr'SimulationData/SLANS/{subdir}'):
                    os.mkdir(projectDir / fr'SimulationData/SLANS/{subdir}')

                os.mkdir(new_path)
                path = projectDir / fr'SimulationData/SLANS/{subdir}/{fid}'

        print(path)
        print(type(path))
        if os.path.exists(path):
            shutil.rmtree(path)
            os.mkdir(path)
        else:
            os.mkdir(path)

    @staticmethod
    def button_clicked(i):
        return i.text()

    @staticmethod
    def get_axis_field_data(folder, mode):
        axis_field_data = {}

        x, y = [], []
        path = os.path.join(fr"{folder}/cavity_33_{mode}.af")
        with open(path, 'r') as f:
            for ll in f.readlines():
                ll = ll.strip()
                x.append(float(ll.split(' ')[0]))
                y.append(float(ll.split(' ')[1]))

        # get avg x
        # avg_x = sum(x) / len(x)
        # x_shift = [t - avg_x for t in x]

        y_abs = [abs(e) for e in y]

        # RETURN ABSOLUTE FIELD VALUE
        axis_field_data['x'], axis_field_data['y'], axis_field_data['y_abs'] = np.array(x), np.array(y), np.array(y_abs)

        return axis_field_data


if __name__ == '__main__':
    slg = SLANSGeometry()

    cell_pars = [
        [60, 60, 20, 20, 72, 97.9485, 172.041, 0],
        [50, 58, 30, 30, 80, 93.5, 172.041, 0],
        [30, 58, 30, 30, 80, 93.5, 172.041, 0],
        [50, 58, 30, 30, 70, 93.5, 172.041, 0],
        [50, 58, 30, 30, 70, 93.5, 172.041, 0],
        [50, 58, 30, 30, 70, 93.5, 172.041, 0],
        [30, 58, 30, 30, 70, 93.5, 172.041, 0],
        [50, 30, 30, 30, 70, 93.5, 172.041, 0],
        [70, 58, 30, 30, 70, 93.5, 172.041, 0],
        [60, 60, 20, 20, 72, 97.9485, 172.041, 0]]

    # cell_pars = [
    #         [62, 66, 30, 23, 72, 93.5, 172.041, 0],
    #         [62, 66, 30, 23, 72, 93.5, 172.041, 0],
    #         [62, 66, 30, 23, 72, 93.5, 172.041, 0],
    #         [62, 66, 30, 23, 72, 93.5, 172.041, 0]]

    parentDir = r"D:/Dropbox/Files/Test_multicell"
    projectDir = r"SimulationData/SLANS"
    slg.cavity_multicell_full(cells_par=cell_pars, parentDir=parentDir, projectDir=projectDir, fid='0')
