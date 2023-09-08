import json
import os
import subprocess
from math import floor

from scipy.signal import find_peaks

from analysis_modules.eigenmode.SLANS.geometry_manual import Geometry
from analysis_modules.eigenmode.SLANS.slans_code import SLANS
import numpy as np

from utils.file_reader import FileReader
from utils.shared_functions import update_alpha

fr = FileReader()


class SLANSGeometry(Geometry):
    def __init__(self, win=None):
        if win:
            super().__init__(win)

            self.ui = win.ui

    def cavity(self, no_of_cells=1, no_of_modules=1, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
               fid=None, bc=33, f_shift='default', beta=1, n_modes=2, proc=0, beampipes=None,
               parentDir=None, projectDir=None):

        if os.path.exists(projectDir):
            os.chdir(projectDir)
        else:
            print("Folder does not exist.")

        self.fid = fid

        self.set_geom_parameters(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)
        # print(mid_cells_par, l_end_cell_par, self.left_beam_pipe, self.right_beam_pipe)

        self.slans = SLANS(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell,
                           self.right_beam_pipe, self.Jxy_all, self.Jxy_all_bp)

        # create folder
        if fid != "_0":
            check = self.createFolder()

        n = no_of_cells  # Number of cells
        axi_sym = 2  # 1: flat, 2:axisymmetric
        self.unit = 3  # 1:m, 2:cm, 3:mm, 4:mkm
        name_index = 1
        sc = 1
        end_type = 1  # if end_type = 1 the end HALF cell is changed for tuning. If end_type = 2 the WHOLE end cell is changed for tuning
        end_L = 1  # if end_L = 1 the type of end cell is type a (without iris) if end_L = 2 the type of end cell is type b
        end_R = 1

        # Beam pipe length

        if end_L == 1:
            self.Rbp_L = self.ri_L

        if end_R == 1:
            self.Rbp_R = self.ri_R

        # Ellipse conjugate points x,y

        zr12_L, alpha_L = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r
        zr12_R, alpha_R = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r
        zr12_M, alpha_M = self.slans.rz_conjug('mid')  # zr12_R first column is z , second column is r

        if end_L == 2:
            zr12_BPL, alpha_BPL = self.slans.rz_conjug('left')  # zr12_R first column is z , second column is r

        if end_R == 2:
            zr12_BPR, alpha_BPR = self.slans.rz_conjug('right')  # zr12_R first column is z , second column is r

        # Set boundary conditions
        BC_Left = floor(bc/10) # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        BC_Right = bc % 10  # 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal

        filename = f'cavity_{bc}'

        # Write Slans Geometry
        # with open(f"{projectDir}/SimulationData/SLANS/{fid}/{filename}.geo", 'w') as f:
        with open(f"{projectDir}/SimulationData/SLANS/{fid}/{filename}.geo", 'w') as f:
            # print("it got here")
            # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
            # print(n)
            # print(self.WG_mesh)
            f.write('8 {:.0f} {:.0f} 2 {}\n'.format(
                self.Jxy * n + self.Jxy_bp * ((1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2)
                + (1 if self.WG_L > 0 else 0) * self.WG_mesh
                + (1 if self.WG_R > 0 else 0) * self.WG_mesh, self.Jxy, self.unit))
            f.write('10 0 0 0 0 0 0 0 0\n')
            f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.ri_L, self.Jy0, BC_Left))

            if end_L == 2:
                f.write('1 0 {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.Rbp_L, self.Jxy_all_bp[5]
                                                                  + self.Jxy_all_bp[6]
                                                                  + self.Jxy_all_bp[7], BC_Left))

            if self.WG_L > 0:
                if end_L == 2:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L - self.x_L, self.Rbp_L, self.WG_mesh))
                else:
                    f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L, self.Rbp_L, self.WG_mesh))

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
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L + self.WG_R
                                                                        + self.L_L + self.L_R,
                                                                        self.Rbp_R, self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(self.WG_L + self.WG_R
                                                                        + self.L_L + self.L_R,
                                                                        self.Rbp_R, self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(self.WG_L + self.WG_R
                                                                         + self.L_L + self.L_R, self.ri_R,
                                                                         -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6]
                                                                           + self.Jxy_all_bp[7]), BC_Right))

                f.write(
                    '1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(self.WG_L + self.WG_R
                                                              + self.L_L + self.L_R,
                                                              -self.Jy0, BC_Right))

                f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(-(self.Jxy * n + self.Jxy_bp * (
                            (1 if end_R == 2 else 0) / 2 + (1 if end_L == 2 else 0) / 2)
                                                            + (1 if self.WG_L > 0 else 0) * self.WG_mesh
                                                            + (1 if self.WG_R > 0 else 0) * self.WG_mesh)))
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
                            self.WG_L + self.WG_R + self.L_L + self.L_R
                            + 2 * (n - 1) * self.L_M, self.Rbp_R, self.WG_mesh))
                    else:
                        f.write('1 {:g} {:g} 0 1 {:.0f} 0 5 0\n'.format(
                            self.WG_L + self.WG_R + self.L_L
                            + self.L_R + 2 * (n - 1) * self.L_M, self.Rbp_R, self.WG_mesh))

                if end_R == 2:
                    f.write('1 {:g} {:g} 0 1 0 {:.0f} {:.0f} 0\n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, self.ri_R,
                        -(self.Jxy_all_bp[5] + self.Jxy_all_bp[6] + self.Jxy_all_bp[7]), BC_Right))

                f.write('1 {:g} 0 0 1 0 {:.0f} {:.0f} 0\n'.format(
                    self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M, -self.Jy0, BC_Right))

                # gradual mesh decrease
                if self.WG_R > 0:
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + self.L_R + 2 * (n - 1) * self.L_M,
                                                                 -((1 if self.WG_R > 0 else 0) * self.WG_mesh)))

                f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (n - 1) * self.L_M - self.L_M,
                                                             -(self.Jxy * 1)))

                for i in range(n - 1, 1, -1):
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L + self.L_L + 2 * (i - 1) * self.L_M - self.L_M,
                                                                 -(self.Jxy * 1)))

                f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(self.WG_L, -(self.Jxy * 1)))

                if self.WG_L > 0:
                    f.write('1 {:g} 0 0 1 {:.0f} 0 4 0\n'.format(0, -((1 if self.WG_L > 0 else 0) * self.WG_mesh)))

                # # direct mesh decrease
                # f.write('1 0 0 0 1 {:.0f} 0 4 0\n'.format(-(self.Jxy*n + self.Jxy_bp *
                #                                             ((1 if end_R == 2 else 0)/2 + (1 if end_L == 2 else 0)/2)
                #                                             + (1 if self.WG_L > 0 else 0)*self.WG_mesh
                #                                             +(1 if self.WG_R > 0 else 0)*self.WG_mesh)))

                f.write('0 0 0 0 0 0 0 0 0')

        # Slans run
        genmesh_path = fr'{parentDir}\exe\SLANS_exe\genmesh2.exe'
        # filepath = fr'{projectDir}\SimulationData\SLANS\{fid}\{filename}'
        filepath = fr'{projectDir}\SimulationData\SLANS\{fid}\{filename}'

        # folder for exe to write to
        # cwd = fr'{projectDir}\SimulationData\SLANS\{fid}'
        cwd = fr'{projectDir}\SimulationData\SLANS\{fid}'

        # the next two lines suppress pop up windows from the slans codes
        # the slans codes, however, still disrupts windows operation, sadly. This is the case even for the slans tuner
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW

        subprocess.call([genmesh_path, filepath, '-b'], cwd=cwd, startupinfo=startupinfo)
        # path = fr'{projectDir}\SimulationData\SLANS\Cavity{self.fid}'
        path = fr'{projectDir}\SimulationData\SLANS\{self.fid}'

        beta, f_shift, n_modes = 1, 0, no_of_cells + 1

        self.write_dtr(path, filename, beta, f_shift, n_modes)

        slansc_path = fr'{parentDir}\exe\SLANS_exe\slansc'
        slansm_path = fr'{parentDir}\exe\SLANS_exe\slansm'
        slanss_path = fr'{parentDir}\exe\SLANS_exe\slanss'
        slansre_path = fr'{parentDir}\exe\SLANS_exe\slansre'

        # check if corresponding file exists at before the executable is called
        if os.path.exists(fr'{projectDir}\SimulationData\SLANS\{fid}\{filename}.geo'):
            subprocess.call([slansc_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)  # settings, number of modes, etc

            if os.path.exists(fr'{projectDir}\SimulationData\SLANS\{fid}\{filename}.gem'):
                subprocess.call([slansm_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)

                if os.path.exists(fr'{projectDir}\SimulationData\SLANS\{fid}\aslans.mtx') \
                        and os.path.exists(fr'{projectDir}\SimulationData\SLANS\{fid}\bslans.mtx'):
                    subprocess.call([slanss_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)

                    if os.path.exists(fr'{projectDir}\SimulationData\SLANS\{fid}\{filename}.res'):
                        subprocess.call([slansre_path, '{}'.format(filepath), '-b'], cwd=cwd, startupinfo=startupinfo)

        run_save_directory = f"{projectDir}/SimulationData/SLANS/{fid}"

        # save json file
        shape = {'IC': update_alpha(mid_cells_par),
                 'OC': update_alpha(l_end_cell_par),
                 'OC_R': update_alpha(r_end_cell_par)}

        with open(fr"{run_save_directory}\geometric_parameters.json", 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))
        try:
            filename = fr'{run_save_directory}\cavity_33.svl'
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
            norm_length = 2*mid_cells_par[5]
            Eacc = Vacc / (
                        no_of_cells * norm_length * 1e-3)  # for 1 cell factor of 2, circuit and accelerator definition
            Epk_Eacc = Epk / Eacc
            Bpk = (Hpk * 4 * np.pi * 1e-7) * 1e3
            Bpk_Eacc = Bpk / Eacc

            # cel to cell coupling factor
            f_diff = d['FREQUENCY'][no_of_cells - 1] - d['FREQUENCY'][0]
            f_add = (d['FREQUENCY'][no_of_cells - 1] + d['FREQUENCY'][0])
            kcc = 2*f_diff/f_add * 100

            # field flatness•
            ax_field = self.get_axis_field_data(run_save_directory, no_of_cells)
            # get max in each cell
            peaks, _ = find_peaks(ax_field['y_abs'])
            E_abs_peaks = ax_field['y_abs'][peaks]
            # ff = min(E_abs_peaks)/max(E_abs_peaks) * 100
            ff = (1 - ((max(E_abs_peaks) - min(E_abs_peaks))/np.average(E_abs_peaks))) * 100

            d = {
                "Req [mm]": Req,
                "Normalization Length [mm]": norm_length,
                "N cells": no_of_cells,
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

            with open(fr'{run_save_directory}\qois.json', "w") as f:
                json.dump(d, f, indent=4, separators=(',', ': '))
        except (FileNotFoundError, ValueError) as e:
            print("Simulation failed", e)

    def write_dtr(self, path, filename, beta, f_shift, n_modes):
        with open("{}\{}.dtr".format(path, filename), 'w') as f:
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

    def createFolder(self):
        path = os.getcwd()
        path = os.path.join(path, f"SimulationData/SLANS/{self.fid}")
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)
            return "Yes"

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
