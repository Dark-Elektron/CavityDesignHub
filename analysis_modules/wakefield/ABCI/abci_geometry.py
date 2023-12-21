import json
import os
import subprocess
from pathlib import Path

from analysis_modules.wakefield.ABCI.geometry import Geometry
from analysis_modules.wakefield.ABCI.abci_code import ABCI, ABCI_flattop


class ABCIGeometry(Geometry):
    def __init__(self):
        super().__init__()

        # create node_editor folder on initialisation
        self.L_all = None
        self.abci = None
        path = os.getcwd()
        path = os.path.join(path, "node_editor")
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)

        self.fid = 0
        # initiate codes

    def cavity(self, no_of_cells, no_of_modules,
               mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
               fid="_0", MROT=0, beampipes=None,
               bunch_length=50, MT=3, NFS=5000, UBT=0,
               DDZ_SIG=0.1, DDR_SIG=0.1,
               parentDir='', projectDir='', WG_M=None, marker='', sub_dir='', **kwargs):

        # defaults
        RDRIVE, ISIG = 5e-3, 5
        LCRBW = 'F'
        BETA = 1
        LMATPR = 'F',
        LPRW, LPPW, LSVW, LSVWA, LSVWT, LSVWL, LSVF = 'T', 'T', 'T', 'F', 'T', 'T', 'F'
        LSAV, LCPUTM = 'F', 'F'

        # unpack kwargs
        for key, value in kwargs.items():
            if key == 'RADIAL BEAM OFFSET AT (RDRIVE)':
                RDRIVE = value
            if key == 'NUMBER OF WAKE POTENTIAL POINTS (NW)':
                NW = value
            if key == 'WAKE FOR A COUNTER-ROTATING BEAM (LCRBW)':
                LCRBW = value
            if key == 'VELOCITY OF THE BUNCH / C (BETA)':
                BETA = value
            if key == 'PRINTOUT OF CAVITY SHAPE USED (LMATPR)':
                LMATPR = value
            if key == 'PRINTOUT OF WAKE POTENTIALS (LPRW)':
                LPRW = value
            if key == 'LINE-PRINTER PLOT OF WAKE POT. (LPPW)':
                LPPW = value
            if key == 'SAVE WAKE POTENTIALS IN A FILE (LSVW)':
                LSVW = value
            if key == 'SAVE AZIMUTHAL WAKE IN A FILE (LSVWA)':
                LSVWA = value
            if key == 'SAVE TRANSVERSE WAKE IN A FILE (LSVWT)':
                LSVWT = value
            if key == 'SAVE LONGITUDINAL WAKE IN A FILE (LSVWL)':
                LSVWL = value
            if key == 'SAVE FFT RESULTS IN A FILE (LSVF)':
                LSVF = value
            if key == 'SAVE FIELDS INTO FILE (LSAV)':
                LSAV = value
            if key == 'CPUTIME MONITOR ACTIVE (LCPUTM)':
                LCPUTM = value

        # Adding parameter arguments here for testing purposes # fid, fileID
        self.fid = f'{fid}'

        # this checks whether input is from gui or from the optimisation
        if mid_cells_par is not None:
            self.set_geom_parameters(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par)
        else:
            self.set_geom_parameters(no_of_cells)

        if WG_M == '':
            WG_M = self.WG_L
        # print(WG_M)

        self.abci = ABCI(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell,
                         self.right_beam_pipe)

        # output_name = 5
        # SIG = 0.003  # One standard deviation of bunch length
        MROT = MROT  # MROT = 0 for monopole fields, MROT = 1 for dipole fields
        # UBT = 0.3  # The last longitudinal coordinate relative to the head of the beam,
        # up to which the potentials are calculated (defaults 10*Sig). The longer the better resolution of impedance
        bunch_length = bunch_length
        wake_offset = 5e-03
        beam_offset = 5e-03
        sig_var = [x * 1e-3 for x in [bunch_length]]  # bunch length converted to m

        # not needed for our parametrization
        end_type = 1  # if _type = 1 the  HALF cell is changed for tuning.
        # If _type = 2 the WHOLE  cell is changed for tuning
        end_L = 1  # if _L = 1 the type of  cell is type a (without iris) if _L = 2 the type of  cell is type b
        end_R = 1

        for i_out in range(1):  # check!!!
            module_nu = no_of_modules  # Number of cavities in module
            n = no_of_cells  # Number of cells
            SIG = sig_var[i_out]  # One standard deviation of bunch length
            mesh_DDR = min(self.ri_M / 10, DDR_SIG * SIG)  # /lines_per_sigma
            mesh_DDZ = min(self.ri_M / 10, DDZ_SIG * SIG)  # /lines_per_sigma

            #  mesh_DDR = 2.5*1e-3
            #  mesh_DDZ = 2.5*1e-3

            # UBT = 14*SIG
            UBT = UBT

            # # Beam pipe radius for Different type of transitions to beam pipe
            if end_L == 1:
                self.Rbp_L = self.ri_L

            if end_R == 1:
                self.Rbp_R = self.ri_R

            # #  Ellipse conjugate points x,y
            zr12_L, alpha_L = self.abci.rz_conjug('left')  # zr12_R first column is z , second column is r
            zr12_R, alpha_R = self.abci.rz_conjug('right')  # zr12_R first column is z , second column is r
            zr12_M, alpha_M = self.abci.rz_conjug('mid')  # zr12_R first column is z , second column is r

            if end_L == 2:
                zr12_BPL, alpha_BPL = self.abci.rz_conjug('left')  # zr12_R first column is z , second column is r

            if end_R == 2:
                zr12_BPR, alpha_BPR = self.abci.rz_conjug('right')  # zr12_R first column is z , second column is r

            # print("GUI_ABCI:: zr12_L", zr12_L)
            # print("GUI_ABCI:: zr12_R", zr12_R)
            # # print("GUI_ABCI:: zr12_BPL", zr12_BPL)
            # # print("GUI_ABCI:: zr12_BPR", zr12_BPR)
            # print("GUI_ABCI:: zr12_M", zr12_M)
            # #  Write ABCI code

            # create folder for file output set
            self.createFolder(self.fid, projectDir, sub_dir, marker)

            # # write parameters to folder
            # if self.ui.cb_Only_Mid_Cells.checkState() == 2:
            #     self.write_cst_paramters_mid(self.fid)
            # else:
            #     self.write_cst_paramters(self.fid)

            # change save directory
            if sub_dir == '':
                run_save_directory = projectDir / fr'SimulationData\ABCI\{fid}'
            else:
                run_save_directory = projectDir / fr'SimulationData\ABCI\{sub_dir}\{fid}'

            fname = fr'{run_save_directory}\Cavity_MROT_{MROT}.abc'
            # print('filename:: ', fname)

            L_all_increment = 0
            self.L_all = 0
            # print(fname)
            with open(fname, 'w') as f:
                f.write(f' &FILE LSAV = .{LSAV}., ITEST = 0, LREC = .F., LCPUTM = .{LCPUTM}. &END \n')
                f.write(' SAMPLE INPUT #1 A SIMPLE CAVITY STRUCTURE \n')
                f.write(' &BOUN  IZL = 3, IZR = 3  &END \n')
                f.write(' &MESH DDR = {}, DDZ = {} &END \n'.format(mesh_DDR, mesh_DDZ))
                f.write(' #CAVITYSHAPE \n')
                f.write('0. \n')
                f.write('0.000 0.000\n')

                if end_L == 2:
                    f.write('{} 0.000\n'.format(self.Rbp_L))
                else:
                    f.write('{} 0.000\n'.format(self.ri_L))

                if self.WG_L > 0:
                    if end_L == 2:
                        f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L))
                    else:
                        f.write('{} {} \n'.format(self.ri_L, self.WG_L))

                if n == 1:
                    for i_mode in range(1, module_nu + 1):
                        if i_mode > 0:
                            if self.WG_L > 0:
                                if end_L == 2:
                                    f.write(
                                        '{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L + (i_mode - 1) * self.L_all))
                                else:
                                    f.write('{} {} \n'.format(self.ri_L, self.WG_L + (i_mode - 1) * self.L_all))

                        if self.Req_L != self.Req_R:
                            print('Error:: The equator radius of left and right cell are not equal')

                        # if exist('L_M') != 1:
                        #     L_M = []

                        if end_L == 2:
                            self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + (i_mode - 1) * self.L_all, f)

                        # print("GUI_ABCI::It got here")
                        self.abci.abci_n1_L(n, zr12_L, self.WG_L + (i_mode - 1) * self.L_all, f)
                        self.abci.abci_n1_R(n, zr12_R, self.WG_L + (i_mode - 1) * self.L_all, f)

                        if end_R == 2:
                            self.abci.abci_bp_R(n, zr12_BPR, self.WG_L + (i_mode - 1) * self.L_all, f)

                        if self.WG_R > 0:
                            if end_R == 2:
                                f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R + self.L_L
                                                          + self.L_R + (i_mode - 1) * self.L_all))
                            else:
                                f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L
                                                          + self.L_R + (i_mode - 1) * self.L_all))

                    f.write(
                        '0 {} \n'.format(self.WG_L + self.WG_R + self.L_L + self.L_R + (module_nu - 1) * self.L_all))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                # #  n>1 multi-cell cavity

                if n > 1:
                    for i_mode in range(1, module_nu + 1):

                        # print("imode:", i_mode, i_mode)
                        # change waveguide length
                        if module_nu == 2:
                            if i_mode == 1:
                                self.WG_R = WG_M
                            if i_mode == 2:
                                self.WG_L = WG_M
                                self.WG_R = 4 * self.L_M

                        elif module_nu > 2:
                            if i_mode == 1:
                                self.WG_R = WG_M
                            elif 1 < i_mode < module_nu:
                                self.WG_L = WG_M
                            else:
                                self.WG_R = 4 * self.L_M
                        # Total length of each cavity
                        L_all_increment = self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M
                        # print(self.WG_L, self.WG_R, WG_M, self.L_all)

                        if i_mode > 1:
                            if self.WG_L > 0:
                                if end_L == 2:
                                    # f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L
                                    # + (i_mode-1)*self.L_all))
                                    f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L + self.L_all))
                                else:
                                    # f.write('{} {} \n'.format(self.ri_L, self.WG_L + (i_mode-1)*self.L_all))
                                    f.write('{} {} \n'.format(self.ri_L, self.WG_L + self.L_all))

                        if end_L == 2:
                            # self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + (i_mode-1)*self.L_all, f)
                            self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + self.L_all, f)

                        # self.abci.abci_n1_L(n, zr12_L, self.WG_L + (i_mode-1)*self.L_all, f)
                        self.abci.abci_n1_L(n, zr12_L, self.WG_L + self.L_all, f)

                        for i in range(1, n):
                            # self.abci.abci_M(n, zr12_M, self.WG_L + (i_mode-1)*self.L_all, f, i, end_type)
                            self.abci.abci_M(n, zr12_M, self.WG_L + self.L_all, f, i, end_type)

                        # self.abci.abci_n1_R(n, zr12_R, self.WG_L + (i_mode-1)*self.L_all, f)
                        self.abci.abci_n1_R(n, zr12_R, self.WG_L + self.L_all, f)

                        if end_R == 2:
                            # self.abci.abci_bp_R(n, zr12_BPR, self.WG_R + (i_mode-1)*self.L_all, f)
                            self.abci.abci_bp_R(n, zr12_BPR, self.WG_R + self.L_all, f)

                        if self.WG_R > 0:
                            if end_R == 2:
                                # f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R+ self.L_L + self.L_R
                                # + 2*(n-1)*self.L_M+(i_mode-1)*self.L_all))
                                f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R + self.L_L
                                                          + self.L_R + 2 * (n - 1) * self.L_M + self.L_all))
                            else:
                                # f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L
                                # + self.L_R+2*(n-1)*self.L_M + (i_mode-1)*self.L_all))
                                f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L
                                                          + self.L_R + 2 * (n - 1) * self.L_M + self.L_all))

                        if i_mode < no_of_modules:
                            self.L_all += L_all_increment

                    # f.write('0 {} \n'.format(self.WG_L + self.WG_R+ self.L_L
                    # + self.L_R+2*(n-1)*self.L_M+(module_nu-1)*self.L_all))
                    f.write('0 {} \n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M + self.L_all))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                f.write(f' &BEAM  SIG = {SIG}, ISIG = {ISIG}, RDRIVE = {RDRIVE}, MROT = {MROT}  &END \n')
                # f.write(' &BEAM  SIG = {}, MROT = {}, RDRIVE = {}  &END \n'.format(SIG, MROT, beam_offset))
                f.write(f' &TIME  MT = {int(MT)} &END \n')
                f.write(f' &WAKE  UBT = {int(UBT)}, LCRBW = .{LCRBW}. &END \n')  # , NFS = {NFS}
                # f.write(' &WAKE  UBT = {}, LCHIN = .F., LNAPOLY = .F., LNONAP = .F. &END \n'.format(UBT, wake_offset))
                # f.write(' &WAKE R  = {}   &END \n'.format(wake_offset))
                f.write(f' &PLOT  LCAVIN = .T., LCAVUS = .F., LPLW = .T., LFFT = .T., LSPEC = .T., '
                        f'LINTZ = .F., LPATH = .T. &END \n')
                f.write(f' &PRIN  LMATPR = .{LMATPR}., LPRW = .{LPRW}., LPPW = .{LPPW}., LSVW = .{LSVW}., '
                        f'LSVWA = .{LSVWA}., LSVWT = .{LSVWT}., LSVWL = .{LSVWL}.,  LSVF = .{LSVF}.   &END\n')
                f.write('\nSTOP\n')

            abci_path = os.getcwd()

            exe_path = os.path.join(abci_path, parentDir / fr'exe\ABCI_exe\ABCI_MP64+.exe')

            if LCPUTM == 'T':
                subprocess.call([exe_path, Path(fr'{run_save_directory}\Cavity_MROT_{MROT}.abc')])
            else:
                subprocess.call([exe_path, Path(fr'{run_save_directory}\Cavity_MROT_{MROT}.abc')],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # save json file
            shape = {'IC': list(mid_cells_par),
                     'OC': list(l_end_cell_par),
                     'OC_R': list(r_end_cell_par)}

            with open(Path(fr"{run_save_directory}\geometric_parameters.json"), 'w') as f:
                json.dump(shape, f, indent=4, separators=(',', ': '))

    def cavity_flattop(self, no_of_cells, no_of_modules,
               mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
               fid="_0", MROT=0, beampipes=None,
               bunch_length=50, MT=3, NFS=5000, UBT=0,
               DDZ_SIG=0.1, DDR_SIG=0.1,
               parentDir='', projectDir='', WG_M=None, marker='', sub_dir='', **kwargs):

        # defaults
        RDRIVE, ISIG = 5e-3, 5
        LCRBW = 'F'
        BETA = 1
        LMATPR = 'F',
        LPRW, LPPW, LSVW, LSVWA, LSVWT, LSVWL, LSVF = 'T', 'T', 'T', 'F', 'T', 'T', 'F'
        LSAV, LCPUTM = 'F', 'F'

        # unpack kwargs
        for key, value in kwargs.items():
            if key == 'RADIAL BEAM OFFSET AT (RDRIVE)':
                RDRIVE = value
            if key == 'NUMBER OF WAKE POTENTIAL POINTS (NW)':
                NW = value
            if key == 'WAKE FOR A COUNTER-ROTATING BEAM (LCRBW)':
                LCRBW = value
            if key == 'VELOCITY OF THE BUNCH / C (BETA)':
                BETA = value
            if key == 'PRINTOUT OF CAVITY SHAPE USED (LMATPR)':
                LMATPR = value
            if key == 'PRINTOUT OF WAKE POTENTIALS (LPRW)':
                LPRW = value
            if key == 'LINE-PRINTER PLOT OF WAKE POT. (LPPW)':
                LPPW = value
            if key == 'SAVE WAKE POTENTIALS IN A FILE (LSVW)':
                LSVW = value
            if key == 'SAVE AZIMUTHAL WAKE IN A FILE (LSVWA)':
                LSVWA = value
            if key == 'SAVE TRANSVERSE WAKE IN A FILE (LSVWT)':
                LSVWT = value
            if key == 'SAVE LONGITUDINAL WAKE IN A FILE (LSVWL)':
                LSVWL = value
            if key == 'SAVE FFT RESULTS IN A FILE (LSVF)':
                LSVF = value
            if key == 'SAVE FIELDS INTO FILE (LSAV)':
                LSAV = value
            if key == 'CPUTIME MONITOR ACTIVE (LCPUTM)':
                LCPUTM = value

        # Adding parameter arguments here for testing purposes # fid, fileID
        self.fid = f'{fid}'

        # this checks whether input is from gui or from the optimisation
        if mid_cells_par is not None:
            self.set_geom_parameters_flattop(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par)
        else:
            self.set_geom_parameters_flattop(no_of_cells)

        if WG_M == '':
            WG_M = self.WG_L
        # print(WG_M)

        self.abci = ABCI_flattop(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell,
                         self.right_beam_pipe)

        # output_name = 5
        # SIG = 0.003  # One standard deviation of bunch length
        MROT = MROT  # MROT = 0 for monopole fields, MROT = 1 for dipole fields
        # UBT = 0.3  # The last longitudinal coordinate relative to the head of the beam,
        # up to which the potentials are calculated (defaults 10*Sig). The longer the better resolution of impedance
        bunch_length = bunch_length
        wake_offset = 5e-03
        beam_offset = 5e-03
        sig_var = [x * 1e-3 for x in [bunch_length]]  # bunch length converted to m

        # not needed for our parametrization
        end_type = 1  # if _type = 1 the  HALF cell is changed for tuning.
        # If _type = 2 the WHOLE  cell is changed for tuning
        end_L = 1  # if _L = 1 the type of  cell is type a (without iris) if _L = 2 the type of  cell is type b
        end_R = 1

        for i_out in range(1):  # check!!!
            module_nu = no_of_modules  # Number of cavities in module
            n = no_of_cells  # Number of cells
            SIG = sig_var[i_out]  # One standard deviation of bunch length
            mesh_DDR = min(self.ri_M / 10, DDR_SIG * SIG)  # /lines_per_sigma
            mesh_DDZ = min(self.ri_M / 10, DDZ_SIG * SIG)  # /lines_per_sigma

            #  mesh_DDR = 2.5*1e-3
            #  mesh_DDZ = 2.5*1e-3

            # UBT = 14*SIG
            UBT = UBT

            # # Beam pipe radius for Different type of transitions to beam pipe
            if end_L == 1:
                self.Rbp_L = self.ri_L

            if end_R == 1:
                self.Rbp_R = self.ri_R

            # #  Ellipse conjugate points x,y
            zr12_L, alpha_L = self.abci.rz_conjug('left')  # zr12_R first column is z , second column is r
            zr12_R, alpha_R = self.abci.rz_conjug('right')  # zr12_R first column is z , second column is r
            zr12_M, alpha_M = self.abci.rz_conjug('mid')  # zr12_R first column is z , second column is r

            if end_L == 2:
                zr12_BPL, alpha_BPL = self.abci.rz_conjug('left')  # zr12_R first column is z , second column is r

            if end_R == 2:
                zr12_BPR, alpha_BPR = self.abci.rz_conjug('right')  # zr12_R first column is z , second column is r

            # print("GUI_ABCI:: zr12_L", zr12_L)
            # print("GUI_ABCI:: zr12_R", zr12_R)
            # # print("GUI_ABCI:: zr12_BPL", zr12_BPL)
            # # print("GUI_ABCI:: zr12_BPR", zr12_BPR)
            # print("GUI_ABCI:: zr12_M", zr12_M)
            # #  Write ABCI code

            # create folder for file output set
            self.createFolder(self.fid, projectDir, sub_dir, marker)

            # # write parameters to folder
            # if self.ui.cb_Only_Mid_Cells.checkState() == 2:
            #     self.write_cst_paramters_mid(self.fid)
            # else:
            #     self.write_cst_paramters(self.fid)

            # change save directory
            if sub_dir == '':
                run_save_directory = projectDir / fr'SimulationData\ABCI\{fid}'
            else:
                run_save_directory = projectDir / fr'SimulationData\ABCI\{sub_dir}\{fid}'

            fname = fr'{run_save_directory}\Cavity_MROT_{MROT}.abc'
            # print('filename:: ', fname)

            L_all_increment = 0
            self.L_all = 0
            # print(fname)
            with open(fname, 'w') as f:
                f.write(f' &FILE LSAV = .{LSAV}., ITEST = 0, LREC = .F., LCPUTM = .{LCPUTM}. &END \n')
                f.write(' SAMPLE INPUT #1 A SIMPLE CAVITY STRUCTURE \n')
                f.write(' &BOUN  IZL = 3, IZR = 3  &END \n')
                f.write(' &MESH DDR = {}, DDZ = {} &END \n'.format(mesh_DDR, mesh_DDZ))
                f.write(' #CAVITYSHAPE \n')
                f.write('0. \n')
                f.write('0.000 0.000\n')

                if end_L == 2:
                    f.write('{} 0.000\n'.format(self.Rbp_L))
                else:
                    f.write('{} 0.000\n'.format(self.ri_L))

                if self.WG_L > 0:
                    if end_L == 2:
                        f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L))
                    else:
                        f.write('{} {} \n'.format(self.ri_L, self.WG_L))

                if n == 1:
                    for i_mode in range(1, module_nu + 1):
                        if i_mode > 0:
                            if self.WG_L > 0:
                                if end_L == 2:
                                    f.write(
                                        '{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L + (i_mode - 1) * self.L_all))
                                else:
                                    f.write('{} {} \n'.format(self.ri_L, self.WG_L + (i_mode - 1) * self.L_all))

                        if self.Req_L != self.Req_R:
                            print('Error:: The equator radius of left and right cell are not equal')

                        # if exist('L_M') != 1:
                        #     L_M = []

                        if end_L == 2:
                            self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + (i_mode - 1) * self.L_all, f)

                        # print("GUI_ABCI::It got here")
                        self.abci.abci_n1_L(n, zr12_L, self.WG_L + (i_mode - 1) * self.L_all, f)

                        # add flattop
                        f.write('{} {} \n'.format(self.Req_L, self.WG_L + self.L_L + self.l_L + (i_mode - 1) * self.L_all))

                        self.abci.abci_n1_R(n, zr12_R, self.WG_L + (i_mode - 1) * self.L_all, f)

                        if end_R == 2:
                            self.abci.abci_bp_R(n, zr12_BPR, self.WG_L + (i_mode - 1) * self.L_all, f)

                        if self.WG_R > 0:
                            if end_R == 2:
                                f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R + self.L_L + self.l_L
                                                          + self.L_R + (i_mode - 1) * self.L_all))
                            else:
                                f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L + self.l_L + self.L_R + (i_mode - 1) * self.L_all))

                    f.write(
                        '0 {} \n'.format(self.WG_L + self.WG_R + self.L_L + self.l_L + self.L_R + (module_nu - 1) * self.L_all))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                # #  n>1 multi-cell cavity

                if n > 1:
                    for i_mode in range(1, module_nu + 1):

                        # print("imode:", i_mode, i_mode)
                        # change waveguide length
                        if module_nu == 2:
                            if i_mode == 1:
                                self.WG_R = WG_M
                            if i_mode == 2:
                                self.WG_L = WG_M
                                self.WG_R = 4 * self.L_M

                        elif module_nu > 2:
                            if i_mode == 1:
                                self.WG_R = WG_M
                            elif 1 < i_mode < module_nu:
                                self.WG_L = WG_M
                            else:
                                self.WG_R = 4 * self.L_M
                        # Total length of each cavity
                        L_all_increment = self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M
                        # print(self.WG_L, self.WG_R, WG_M, self.L_all)

                        if i_mode > 1:
                            if self.WG_L > 0:
                                if end_L == 2:
                                    # f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L
                                    # + (i_mode-1)*self.L_all))
                                    f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L + self.L_all))
                                else:
                                    # f.write('{} {} \n'.format(self.ri_L, self.WG_L + (i_mode-1)*self.L_all))
                                    f.write('{} {} \n'.format(self.ri_L, self.WG_L + self.L_all))

                        if end_L == 2:
                            # self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + (i_mode-1)*self.L_all, f)
                            self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + self.L_all, f)

                        # self.abci.abci_n1_L(n, zr12_L, self.WG_L + (i_mode-1)*self.L_all, f)
                        self.abci.abci_n1_L(n, zr12_L, self.WG_L + self.L_all, f)

                        # add flattop
                        f.write('{} {} \n'.format(self.Req_L, self.WG_L + self.L_L + self.l_L + (i_mode - 1) * self.L_all))

                        for i in range(1, n):
                            # self.abci.abci_M(n, zr12_M, self.WG_L + (i_mode-1)*self.L_all, f, i, end_type)
                            self.abci.abci_M(n, zr12_M, self.WG_L + self.L_all, f, i, end_type)

                            if i != n-1:
                                # add flattop
                                f.write('{} {} \n'.format(self.Req_M, self.WG_L + self.L_L + 2 * i * self.L_M + self.l_L + i*self.l_M + (i_mode - 1) * self.L_all))
                            else:
                                # add flattop
                                f.write('{} {} \n'.format(self.Req_M, self.WG_L + self.L_L + 2 * i * self.L_M + (i_mode - 1) * self.L_all + self.l_L + (i-1)*self.l_M + self.l_R))

                        # self.abci.abci_n1_R(n, zr12_R, self.WG_L + (i_mode-1)*self.L_all, f)
                        self.abci.abci_n1_R(n, zr12_R, self.WG_L + self.L_all, f)

                        if end_R == 2:
                            # self.abci.abci_bp_R(n, zr12_BPR, self.WG_R + (i_mode-1)*self.L_all, f)
                            self.abci.abci_bp_R(n, zr12_BPR, self.WG_R + self.L_all, f)

                        if self.WG_R > 0:
                            if end_R == 2:
                                # f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R+ self.L_L + self.L_R
                                # + 2*(n-1)*self.L_M+(i_mode-1)*self.L_all))
                                f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R + self.L_L + self.l_L + (n-2)*self.l_M + self.l_R
                                                          + self.L_R + 2 * (n - 1) * self.L_M + self.L_all))
                            else:
                                # f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L
                                # + self.L_R+2*(n-1)*self.L_M + (i_mode-1)*self.L_all))
                                f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L + self.l_L + (n-2)*self.l_M + self.l_R
                                                          + self.L_R + 2 * (n - 1) * self.L_M + self.L_all))

                        if i_mode < no_of_modules:
                            self.L_all += L_all_increment

                    # f.write('0 {} \n'.format(self.WG_L + self.WG_R+ self.L_L
                    # + self.L_R+2*(n-1)*self.L_M+(module_nu-1)*self.L_all))
                    f.write('0 {} \n'.format(
                        self.WG_L + self.WG_R + self.L_L + self.L_R + 2 * (n - 1) * self.L_M + self.L_all + self.l_L + (n-2)*self.l_M + self.l_R))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                f.write(f' &BEAM  SIG = {SIG}, ISIG = {ISIG}, RDRIVE = {RDRIVE}, MROT = {MROT}  &END \n')
                # f.write(' &BEAM  SIG = {}, MROT = {}, RDRIVE = {}  &END \n'.format(SIG, MROT, beam_offset))
                f.write(f' &TIME  MT = {int(MT)} &END \n')
                f.write(f' &WAKE  UBT = {int(UBT)}, LCRBW = .{LCRBW}. &END \n')  # , NFS = {NFS}
                # f.write(' &WAKE  UBT = {}, LCHIN = .F., LNAPOLY = .F., LNONAP = .F. &END \n'.format(UBT, wake_offset))
                # f.write(' &WAKE R  = {}   &END \n'.format(wake_offset))
                f.write(f' &PLOT  LCAVIN = .T., LCAVUS = .F., LPLW = .T., LFFT = .T., LSPEC = .T., '
                        f'LINTZ = .F., LPATH = .T. &END \n')
                f.write(f' &PRIN  LMATPR = .{LMATPR}., LPRW = .{LPRW}., LPPW = .{LPPW}., LSVW = .{LSVW}., '
                        f'LSVWA = .{LSVWA}., LSVWT = .{LSVWT}., LSVWL = .{LSVWL}.,  LSVF = .{LSVF}.   &END\n')
                f.write('\nSTOP\n')

            abci_path = os.getcwd()

            exe_path = os.path.join(abci_path, parentDir / fr'exe\ABCI_exe\ABCI_MP64+.exe')

            if LCPUTM == 'T':
                subprocess.call([exe_path, Path(fr'{run_save_directory}\Cavity_MROT_{MROT}.abc')])
            else:
                subprocess.call([exe_path, Path(fr'{run_save_directory}\Cavity_MROT_{MROT}.abc')],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # save json file
            shape = {'IC': list(mid_cells_par),
                     'OC': list(l_end_cell_par),
                     'OC_R': list(r_end_cell_par)}

            with open(Path(fr"{run_save_directory}\geometric_parameters.json"), 'w') as f:
                json.dump(shape, f, indent=4, separators=(',', ': '))

    @staticmethod
    def createFolder(fid, projectDir, subdir='', marker=''):
        # path = os.path.join(path, f"{projectDir}\SimulationData\ABCI\{fid}{marker}")
        # if os.path.exists(path):
        #     pass
        # else:
        #     os.mkdir(path)

        # change save directory
        if subdir == '':
            path = projectDir / fr'SimulationData\ABCI\{fid}{marker}'
            if os.path.exists(path):
                pass
            else:
                os.mkdir(path)
        else:
            new_path = projectDir / fr'SimulationData\ABCI\{subdir}{marker}\{fid}'
            if os.path.exists(projectDir / fr'SimulationData\ABCI\{subdir}{marker}'):
                if os.path.exists(new_path):
                    pass
                else:
                    os.mkdir(new_path)
            else:
                os.mkdir(projectDir / fr'SimulationData\ABCI\{subdir}{marker}')
                os.mkdir(new_path)
