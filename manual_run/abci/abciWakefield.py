import os
import shutil
import subprocess
from geometry import Geometry
from abci_code import ABCI


class ABCIWakefield(Geometry):
    def run(self, no_of_cells, no_of_modules, mid_cells_par, l_end_cell_par, r_end_cell_par,
               name, path=None, MROT=50, bunch_length=50, MT=4, NFS=10000, UBT=50, DDZ_SIG=0.1, DDR_SIG=0.1):

        """
        :param no_of_cells: Number of cells in cavity <type: int>
        :param no_of_modules: Number of cavity analysis_modules <type: int>
        :param mid_cells_par: Mid cell parameters in the order [A, B, a, b
        :param l_end_cell_par: Left end cell parameters in the order [A, B
        :param r_end_cell_par: Right end cell parameters in the order [A,
        :param name: Name of folder to save analysis results. <type: str>
        :param path: Path to save analysis results. <type: str>
        :param MROT: Polarisation. 0: Monople, 1: Dipole <type: int>
        :param bunch_length: Bunch length in mm <type: int>
        :param MT:
        :param NFS: Number of frequency samples <type: int>
        :param UBT: Wakelength in metres <type: float>
        :param DDZ_SIG: Ratio of mesh size to sigma in Z direction <type: float>
        :param DDR_SIG: Ratio of mesh size to sigma in R direction <type: float>
        :return:
        """

        if not path:
            path = os.getcwd()

        self.name = f'{name}'

        self.set_geom_parameters(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par)

        self.abci = ABCI(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell, self.right_beam_pipe)

        output_name = 5
        SIG = 0.003  # One standard deviation of bunch length
        MROT = MROT # MROT = 0 for monopole fields, MROT = 1 for dipole fields
        # UBT = 0.3  # The last longitudinal coordinate relative to the head of the beam, up to which the potentials are calculated (defaults 10*Sig). The longer the better resolution of impedance
        bunch_length = bunch_length
        sig_var = [x*1e-3 for x in [bunch_length]] # bunch length converted to m

        # not needed for our parametrization
        end_type = 1 # if _type = 1 the  HALF cell is changed for tuning. If _type = 2 the WHOLE  cell is changed for tuning
        end_L = 1 # if _L = 1 the type of  cell is type a (without iris) if _L = 2 the type of  cell is type b
        end_R = 1

        for i_out in range(1): # check!!!
            module_nu = no_of_modules # Number of cavities in module
            n = no_of_cells # Number of cells
            SIG = sig_var[i_out] # One standard deviation of bunch length
            mesh_DDR = DDR_SIG*SIG
            mesh_DDZ = DDZ_SIG*SIG

            # # Beam pipe radius for Different type of transitions to beam pipe
            if end_L == 1:
                self.Rbp_L = self.ri_L

            if end_R == 1:
                self.Rbp_R = self.ri_R

            # #  Ellipse conjugate points x,y
            zr12_L, alpha_L = self.abci.rz_conjug('left') # zr12_R first column is z , second column is r
            zr12_R, alpha_R = self.abci.rz_conjug('right') # zr12_R first column is z , second column is r
            zr12_M, alpha_M = self.abci.rz_conjug('mid') # zr12_R first column is z , second column is r

            if end_L == 2:
                zr12_BPL, alpha_BPL = self.abci.rz_conjug('left') # zr12_R first column is z , second column is r

            if end_R == 2:
                zr12_BPR, alpha_BPR = self.abci.rz_conjug('right') # zr12_R first column is z , second column is r

            # create folder for file output set
            self.createFolder(name, path)

            fname = fr'{path}\{name}\Cavity_MROT_{MROT}.abc'

            with open(fname, 'w') as f:
                f.write(' &FILE LSAV = .F., ITEST = 0, LREC = .F. &END \n')
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
                    for i_mode in range(1, module_nu+1):
                        if i_mode > 0:
                            if self.WG_L > 0:
                                if end_L == 2:
                                    f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L+(i_mode-1)*self.L_all))
                                else:
                                    f.write('{} {} \n'.format(self.ri_L, self.WG_L + (i_mode-1)*self.L_all))

                        if self.Req_L != self.Req_R:
                            print('Error:: The equator radius of left and right cell are not equal')

                        # if exist('L_M') != 1:
                        #     L_M = []

                        if end_L == 2:
                            self.abci.abci_bp_L(n, zr12_BPL, self.WG_L+(i_mode-1)*self.L_all, f)

                        # print("GUI_ABCI::It got here")
                        self.abci.abci_n1_L(n, zr12_L, self.WG_L+(i_mode-1)*self.L_all, f)
                        self.abci.abci_n1_R(n, zr12_R, self.WG_L+(i_mode-1)*self.L_all, f)

                        if end_R == 2:
                            self.abci.abci_bp_R(n, zr12_BPR, self.WG_L+(i_mode-1)*self.L_all, f)

                        if self.WG_R > 0:
                            if end_R == 2:
                                f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R+ self.L_L + self.L_R+(i_mode-1)*self.L_all))
                            else:
                                f.write('{} {} \n'.format(self.ri_R, self.WG_L+self.WG_R+ self.L_L + self.L_R + (i_mode-1)*self.L_all))

                    f.write('0 {} \n'.format(self.WG_L + self.WG_R + self.L_L + self.L_R + (module_nu-1)*self.L_all))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                # #  n>1 multicell cavity
                if n > 1:
                    for i_mode in range(1, module_nu+1):
                        if i_mode > 1:
                            if self.WG_L > 0:
                                if end_L == 2:
                                    f.write('{} {} \n'.format(self.Rbp_L, self.WG_L - self.x_L + (i_mode-1)*self.L_all))
                                else:
                                    f.write('{} {} \n'.format(self.ri_L, self.WG_L + (i_mode-1)*self.L_all))

                        if end_L == 2:
                            self.abci.abci_bp_L(n, zr12_BPL, self.WG_L + (i_mode-1)*self.L_all, f)

                        self.abci.abci_n1_L(n, zr12_L, self.WG_L + (i_mode-1)*self.L_all, f)

                        for i in range(1, n):
                            self.abci.abci_M(n, zr12_M, self.WG_L + (i_mode-1)*self.L_all, f, i, end_type)

                        self.abci.abci_n1_R(n, zr12_R, self.WG_L + (i_mode-1)*self.L_all, f)

                        if end_R == 2:
                            self.abci.abci_bp_R(n, zr12_BPR, self.WG_L + (i_mode-1)*self.L_all, f)

                        if self.WG_R > 0:
                            if end_R == 2:
                                f.write('{} {} \n'.format(self.Rbp_R, self.WG_L + self.WG_R+ self.L_L + self.L_R + 2*(n-1)*self.L_M+(i_mode-1)*self.L_all))
                            else:
                                f.write('{} {} \n'.format(self.ri_R, self.WG_L + self.WG_R + self.L_L+ self.L_R+2*(n-1)*self.L_M + (i_mode-1)*self.L_all))

                    f.write('0 {} \n'.format(self.WG_L + self.WG_R+ self.L_L + self.L_R+2*(n-1)*self.L_M+(module_nu-1)*self.L_all))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                f.write(' &BEAM  SIG = {}, MROT = {}  &END \n'.format(SIG, MROT))
                f.write(' &TIME  MT = {} &END \n'.format(MT))
                f.write(' &WAKE  UBT = {} &END \n'.format(UBT))
                f.write(' &WAKE   &END \n')
                f.write(' &PLOT  LCAVIN = .T., LCAVUS = .T., LPLW = .T., LFFT = .T., LSPEC = .T., LINTZ = .F., LPATH = .T., NFS = {} &END \n'.format(NFS))
                f.write(' &PRIN  LPRW = .T. ,LSVW = .T., LSVWL = .T.,  LSVF = .T.   &END\n')
                f.write('\nSTOP\n')

            abci_path = os.getcwd()

            exe_path = fr'{abci_path}\ABCI_exe\ABCI_MP64+.exe'
            subprocess.call([exe_path, fr'{abci_path}\{name}\Cavity_MROT_{MROT}.abc'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def createFolder(name, path):
        path = fr"{path}\{name}"
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)


if __name__ == '__main__':
    # create ABCIWakefield object
    abciwakefield = ABCIWakefield()

    # instantiate parameters
    no_of_cells = 1
    no_of_modules = 1
    mid_cells_parameters = [67.72, 67.72, 21.75, 21.75, 60, 93.5, 171.381] # [A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m]
    left_end_cell_parameters = [67.72, 67.72, 21.75, 21.75, 60, 93.5, 171.381] # [A_e, B_e, a_e, b_e, Ri_e, L_e, Req_e]
    right_end_cell_parameters = left_end_cell_parameters
    name = "Cavity"
    path = None
    MROT = 0
    bunch_length = 25
    MT = 4
    NFS = 10000 # Number of frequency samples
    UBT = 50 # wakelength
    DDZ_SIG = 0.1 # Ratio of mesh size to sigma in Z direction
    DDR_SIG = 0.1 # Ratio of mesh size to sigma in R direction

    # run wakefield analysis
    abciwakefield.run(no_of_cells, no_of_modules, mid_cells_parameters, left_end_cell_parameters, right_end_cell_parameters,
               name, path, MROT=MROT, bunch_length=bunch_length, MT=MT, NFS=NFS, UBT=UBT, DDZ_SIG=DDZ_SIG, DDR_SIG=DDR_SIG)