import json
import os
import shutil
import subprocess
from math import floor
import numpy as np
from utils.file_reader import FileReader
from geometry_manual import Geometry
from slans_code import SLANS
from utils.shared_functions import update_alpha

fr = FileReader()
file_color = 'yellow'
DEBUG = False


class SLANSEigen(Geometry):
    def run(self, no_of_cells, no_of_modules, mid_cells_par, l_end_cell_par, r_end_cell_par,
            name, bc=33, f_shift=0, beta=1, n_modes=None, beampipes=None, path=None, expansion=None, expansion_r=None):
        """
        :param no_of_cells: Number of cells in cavity <type: int>
        :param no_of_modules: Number of cavity analysis_modules <type: int>
        :param mid_cells_par: Mid cell parameters in the order [A, B, a, b, Ri, L, Req] <type: list>
        :param l_end_cell_par: Left end cell parameters in the order [A, B, a, b, Ri, L, Req] <type: list>
        :param r_end_cell_par: Right end cell parameters in the order [A, B, a, b, Ri, L, Req] <type: list>
        :param name: Name of folder to save analysis results. <type: str>
        :param bc: Boundary conditions. 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
                Example: bc = 23 implies Electric wall at left boundary and Magnetic Wall at right boundary. <type: int>
        :param f_shift: Frequency shift for eigenmode anlaysis. <type: int>
        :param beta: Particle velocity. <type: int>
        :param n_modes: No of modes to calculate for. <type: int>
        :param beampipes: "left", "right", "none", "both" <type: str>
        :param path: Path to save analysis results. <type: str>
        :return:
        """

        if not path:
            path = os.getcwd()

        if not n_modes:
            n_modes = no_of_cells + 1  # intentional because the accuracy of the last mode is always low

        self.name = name
        self.set_geom_parameters(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par,
                                 beampipes, expansion=expansion, expansion_r=expansion_r)

        self.slans = SLANS(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell,
                           self.right_beam_pipe, self.Jxy_all, self.Jxy_all_bp)

        n = no_of_cells  # Number of cells
        axi_sym = 2  # 1: flat, 2: axisymmetric
        self.unit = 3  # 1:m, 2:cm, 3:mm, 4:mkm
        name_index = 1
        sc = 1
        end_type = 1  # if end_type = 1 the end HALF cell is changed for tuning. If end_type = 2 the WHOLE end cell is changed for tuning
        end_L = 1  # if end_L = 1 the type of end cell is type a (without iris) if end_L = 2 the type of end cell is type b
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

        # Ellipse conjugate points x,y
        zr12_L, alpha_L = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r
        zr12_R, alpha_R = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r
        zr12_M, alpha_M = self.slans.rz_conjug('mid')  # zr12_R first column is z , second column is r

        if end_L == 2:
            zr12_BPL, alpha_BPL = self.slans.rz_conjug('expansion')  # zr12_R first column is z , second column is r
            # zr12_EL, alpha_EL = self.slans.rz_conjug('expansion')  # zr12_R first column is z , second column is r

        if end_R == 2:
            zr12_BPR, alpha_BPR = self.slans.rz_conjug('expansion_r')  # zr12_R first column is z , second column is r
            # zr12_ER, alpha_ER = self.slans.rz_conjug('expansion_r')  # zr12_R first column is z , second column is r

        # Set boundary conditions
        BC_Left = floor(bc / 10)  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        BC_Right = bc % 10  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal

        filename = f'cavity_{bc}'

        self.createFolder(name, path)

        # Write Slans Geometry
        with open(fr'{path}\{name}\{filename}.geo', 'w') as f:
            # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
            f.write('8 {:.0f} {:.0f} 2 {}\n'.format(
                self.Jxy * n + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) + (
                    1 if self.WG_L > 0 else 0) * self.WG_mesh + (1 if self.WG_R > 0 else 0) * self.WG_mesh, self.Jxy,
                self.unit))
            f.write('10 0 0 0 0 0 0 0 0\n')
            f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.ri_L, self.Jy0, BC_Left))

            if end_L == 2:
                f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.Rbp_L, self.Jxy_all_bp[5] + self.Jxy_all_bp[6] +
                                                                  self.Jxy_all_bp[7], BC_Left))

            if self.WG_L > 0:
                if end_L == 2:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L - self.x_L, self.Rbp_L, self.WG_mesh))
                else:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L - self.x_L, self.Rbp_L, self.WG_mesh))

            # # add expansion
            # if expansion is not None:
            #     if expansion[7] == 'left':
            #         self.slans.slans_n1_EL(n, zr12_E, self.WG_L, f)

            # n == 1
            if n == 1:
                if self.Req_L != self.Req_R:
                    print('The equator radius of left and right cell are not equal')

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
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L + self.WG_R + self.L_L + self.L_R,
                                                                        self.Rbp_R,
                                                                        self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L + self.WG_R + self.L_L + self.L_R,
                                                                        self.Rbp_R,
                                                                        self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.WG_L + self.WG_R + self.L_L + self.L_R,
                                                                         self.ri_R,
                                                                         -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] +
                                                                           self.Jxy_all_bp[7]), BC_Right))

                f.write(
                    '1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(self.WG_L + self.WG_R + self.L_L + self.L_R, -self.Jy0,
                                                              BC_Right))

                f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(-(self.Jxy * n + self.Jxy_bp * (
                        (1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) + (
                                                                1 if self.WG_L > 0 else 0) * self.WG_mesh + (
                                                                1 if self.WG_R > 0 else 0) * self.WG_mesh)))
                f.write('0 0 0 0 0 0 0 0 0')

            # n>1
            if n > 1:
                if end_L == 2:
                    # self.slans.slans_bp_L(n, zr12_BPL, self.WG_L, f)
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
                            self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, self.Rbp_R,
                            self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, self.Rbp_R,
                            self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, self.ri_R,
                        -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] + self.Jxy_all_bp[7]), BC_Right))

                f.write('1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(
                    self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, -self.Jy0, BC_Right))

                # gradual mesh decrease
                if self.WG_R > 0:
                    f.write(
                        '1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + self.L_R + 2 * (n - 1) * self.L_M,
                                                             -((1 if self.WG_R > 0 else 0) * self.WG_mesh)))

                f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (n - 1) * self.L_M - self.L_M,
                                                             -(self.Jxy * 1)))

                for i in range(n - 1, 1, -1):
                    f.write(
                        '1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (i - 1) * self.L_M - self.L_M,
                                                             -(self.Jxy * 1)))

                f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L, -(self.Jxy * 1)))

                if self.WG_L > 0:
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(0, -((1 if self.WG_L > 0 else 0) * self.WG_mesh
                                                                      + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2
                                                                                       + (
                                                                                           1 if end_L == 2 else 0) / 2))))

                # # direct mesh decrease
                # f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(-(self.Jxy * n + self.Jxy_bp * (
                #         (1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2) + (
                #                                                 1 if self.WG_L > 0 else 0) * self.WG_mesh + (
                #                                                 1 if self.WG_R > 0 else 0) * self.WG_mesh)))

                f.write('0 0 0 0 0 0 0 0 0')

        # Slans run
        cwd = os.getcwd()
        genmesh_path = fr'{cwd}\SLANS_exe\genmesh2.exe'
        filepath = fr'{path}\{name}\{filename}'

        # folder for exe to write to
        write_folder = fr'{path}\{name}'

        # Choose whether to display SLANS's run GUI
        show_slans_gui = False

        if show_slans_gui:
            startupinfo = None
            subprocess.call([genmesh_path, filepath, '-b'], cwd=write_folder, startupinfo=startupinfo)
        else:
            # The next two lines suppress pop up windows from the slans codes
            # the slans codes, however, still disrupts windows operation, sadly.
            # This is the case even for the slans tuner
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            subprocess.call([genmesh_path, filepath, '-b'], cwd=write_folder, startupinfo=startupinfo)

        self.write_dtr(path, name, filename, beta, f_shift, n_modes)

        slansc_path = fr'{cwd}\SLANS_exe\slansc'
        slansm_path = fr'{cwd}\SLANS_exe\slansm'
        slanss_path = fr'{cwd}\SLANS_exe\slanss'
        slansre_path = fr'{cwd}\SLANS_exe\slansre'

        subprocess.call([slansc_path, '{}'.format(filepath), '-b'], cwd=write_folder, startupinfo=startupinfo)
        subprocess.call([slansm_path, '{}'.format(filepath), '-b'], cwd=write_folder, startupinfo=startupinfo)
        subprocess.call([slanss_path, '{}'.format(filepath), '-b'], cwd=write_folder, startupinfo=startupinfo)
        subprocess.call([slansre_path, '{}'.format(filepath), '-b'], cwd=write_folder, startupinfo=startupinfo)
        # save json file
        shape = {'IC': update_alpha(mid_cells_par),
                 'OC': update_alpha(l_end_cell_par),
                 'OC_R': update_alpha(r_end_cell_par)}

        with open(fr"{write_folder}\geometric_parameters.json", 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))
        try:
            filename = fr'{write_folder}\cavity_33.svl'
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

            # # field flatness
            # ax_field = self.get_axis_field_data(write_folder, no_of_cells)
            # # get max in each cell
            # peaks, _ = find_peaks(ax_field['y_abs'])
            # E_abs_peaks = ax_field['y_abs'][peaks]
            # ff = min(E_abs_peaks) / max(E_abs_peaks) * 100
            # ff = (1 - ((max(E_abs_peaks) - min(E_abs_peaks))/np.average(E_abs_peaks))) * 100

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
                # "ff [%]": ff,
                "Rsh [MOhm]": Rsh,
                "R/Q [Ohm]": 2 * r_Q,
                "Epk/Eacc []": Epk_Eacc,
                "Bpk/Eacc [mT/MV/m]": Bpk_Eacc,
                "G [Ohm]": G,
                "GR/Q [Ohm^2]": GR_Q
            }

            with open(fr'{write_folder}\qois.json', "w") as f:
                json.dump(d, f, indent=4, separators=(',', ': '))
            print(d)
        except FileNotFoundError as e:
            print("Simulation failed", e)

    @staticmethod
    def write_dtr(path, name, filename, beta, f_shift, n_modes):
        with open(fr"{path}\{name}\{filename}.dtr", 'w') as f:
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
    def createFolder(name, path):
        path = fr"{path}\{name}"
        if os.path.exists(path):
            shutil.rmtree(path)
            os.mkdir(path)
        else:
            os.mkdir(path)

    @staticmethod
    def get_axis_field_data(folder, mode):
        axis_field_data = {}

        x, y = [], []
        path = os.path.join(fr"{folder}\cavity_33_{mode}.af")
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
    # create SLANSEigen object
    slanseigen = SLANSEigen()

    # initialize input arguments
    mid_cell_parameters = [43.99, 35.06, 12.53, 20.95, 35, 57.6524, 101.205]  # [A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m]
    left_end_cell_parameters = [52.1, 47.9, 9.9, 11.3, 37, 62.665, 101.205]  # [A_e, B_e, a_e, b_e, Ri_e, L_e, Req_e]
    right_end_cell_paramters = [50.9, 45.3, 8.4, 11.5, 39, 59.988, 101.205]

    # midNLSF_RE = [49, 35.30, 10.5, 17, 32.0, 57.7, 98.58]
    # endNLSF_RE = [50, 35, 10, 15, 32.0, 57.7, 98.58]

    at_L, c_L, c_R = 6, 2*(55 - 37), 2*(55 - 37)
    print(at_L, c_L, c_R)
    x_L = max(1.01*(c_L + at_L), 42.002244034809)
    x_R = 1.01*(c_R + at_L)
    print(x_L, x_R)
    expansion = [c_L, c_L, 6, 6, x_L, 55, 0]
    expansion_r = [c_R, c_R, 6, 6, x_R, 55, 0]

    # mid_cell_parameters = [42, 42, 12, 19, 35, 57.6524, 103.353]  # [A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m]
    # left_end_cell_parameters = [40.34, 40.34, 10, 13.5, 39, 55.716, 103.353]  # [A_e, B_e, a_e, b_e, Ri_e, L_e, Req_e]
    # right_end_cell_paramters = [42, 42, 9, 12.8, 39, 56.815, 103.353]
    #
    # # x_L = 1.01*(c_L + at_L)
    # expansion = [25, 25, 10, 10, 36, 55, 0]
    # expansion_r = [25, 25, 10, 10, 36, 55, 0]

    expansion_r = None
    # expansion = None
    beampipes = "both"  # other options:: "right", "both", "none"
    boundary_condition = 33  # other options: 12, 13, 23, 32, etc. See description in run function

    # run eigenmode analysis
    n_cells = 9
    slanseigen.run(n_cells, 1, mid_cell_parameters, left_end_cell_parameters, left_end_cell_parameters, "Cavity", 33, 0, 1, n_cells,
                   beampipes, expansion=expansion, expansion_r=expansion_r)

    # # run eigenmode analysis
    # slanseigen.run(3, 1, midNLSF_RE, endNLSF_RE, endNLSF_RE, "Cavity", 33, 0, 1, 9,
    #                beampipes, expansion=expansion, expansion_r=expansion_r)
