import os
import subprocess

from PyQt5.QtWidgets import QMessageBox

from geometry import Geometry
from ABCI_code.abci_code import ABCI


class ABCIGeometry(Geometry):
    def __init__(self, win=None):
        if win:
            super().__init__(win)
            self.win = win
            self.ui = win.ui

        # create node_editor folder on initialisation
        path = os.getcwd()
        self.path = os.path.join(path, "ABCI_Data_end")
        if os.path.exists(self.path):
            pass
        else:
            os.mkdir(self.path)

        # initiate codes

    def cavity(self, no_of_cells, no_of_modules, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None, fid="_0"):
        # Adding parameter arguments here for testing purposes # fid, fileID
        self.fid = fid.upper()

        # this checks whether input is from gui or from the optimisation
        # consider instead writing the parameters to file
        if mid_cells_par is not None:
            self.set_geom_parameters(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, code='abci')
        else:
            self.set_geom_parameters(no_of_cells)

        self.abci = ABCI(self.left_beam_pipe, self.left_end_cell, self.mid_cell, self.right_end_cell,
                         self.right_beam_pipe)

        print(self.mid_cell, self.left_end_cell, self.left_beam_pipe)

        mesh_DDR = 0.0007
        mesh_DDZ = 0.0007
        SIG = 0.003  # One standard deviation of bunch length
        MROT = 0  # MROT = 0 for monopole fields, MROT = 1 for dipole fields
        MT = 4  # number of time steps for a beam to move one cell to another default = 3
        UBT = 0.3  # The last longitudinal coordinate relative to the head of the beam, up to which the potentials are calculated (defaults 10*Sig). The longer the better resolution of impedance

        NFS = 10000 # Number of samples in FFT (max 10000)
        sig_var = [x*1e-3 for x in [20]] # bunch length

        # not needed for our parametrization
        end_type = 1 # if _type = 1 the  HALF cell is changed for tuning. If _type = 2 the WHOLE  cell is changed for tuning
        end_L = 1 # if _L = 1 the type of  cell is type a (without iris) if _L = 2 the type of  cell is type b
        end_R = 1
        for i_out in range(1): # check!!!
            module_nu = no_of_modules # Number of cavities in module
            n = no_of_cells # Number of cells
            SIG = sig_var[i_out] # One standard deviation of bunch length
            mesh_DDR = 0.1*SIG
            mesh_DDZ = 0.1*SIG

            #  mesh_DDR = 2.5*1e-3
            #  mesh_DDZ = 2.5*1e-3

            # UBT = 14*SIG
            UBT = 50

            # # Beam pipe radius for Different type of transitions to beam pipe
            if end_L == 1:
                self.Rbp_L = self.ri_L

            if end_R == 1:
                self.Rbp_R = self.ri_R
            # #  Ellipse conjugate points x,y
            print("\t\there again")
            zr12_L, alpha_L = self.abci.rz_conjug('left') # zr12_R first column is z , second column is r
            print("\t\t\tpass1")
            zr12_R, alpha_R = self.abci.rz_conjug('right') # zr12_R first column is z , second column is r
            print("\t\t\tpass2")
            zr12_M, alpha_M = self.abci.rz_conjug('mid') # zr12_R first column is z , second column is r
            print("\t\t\tpass3")
            if end_L == 2:
                zr12_BPL, alpha_BPL = self.abci.rz_conjug('left') # zr12_R first column is z , second column is r

            if end_R == 2:
                zr12_BPR, alpha_BPR = self.abci.rz_conjug('right') # zr12_R first column is z , second column is r

            # create folder for file output set
            print("creating folder")
            self.createFolder()
            print("done creating")
            # if check == "No":
            #     print("About to exit")
            #     exit()
            # else:
            #     print("continue")
            #     pass

            # # write parameters to folder
            # if self.ui.cb_Only_Mid_Cells.checkState() == 2:
            #     self.write_cst_paramters_mid(self.fid)
            # else:
            #     self.write_cst_paramters(self.fid)
            print(self.path)
            fname = ('{}\Cavity{}\Cavity_MROT_{}.abc'.format(self.path, self.fid, MROT))
            print(fname)
            # print(fname)
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

                        print("="*50)

                    f.write('0 {} \n'.format(self.WG_L + self.WG_R+ self.L_L + self.L_R+2*(n-1)*self.L_M+(module_nu-1)*self.L_all))
                    f.write('0 0 \n')
                    f.write('9999. 9999. \n')

                f.write(' &BEAM  SIG = {}, MROT = {}  &END \n'.format(SIG,MROT))
                f.write(' &TIME  MT = {} &END \n'.format(MT))
                f.write(' &WAKE  UBT = {} &END \n'.format(UBT))
                f.write(' &WAKE   &END \n')
                f.write(' &PLOT  LCAVIN = .T., LCAVUS = .T., LPLW = .T., LFFT = .T., LSPEC = .T., LINTZ = .F., LPATH = .T., NFS = {} &END \n'.format(NFS))
                f.write(' &PRIN  LPRW = .T. ,LSVW = .T., LSVWL = .T.,  LSVF = .T.   &END\n')
                f.write('\nSTOP\n')

            # [status, cmdout]  =  system(char(strcat(ABCI_path,{'\ABCI_MP64+.exe'},{' '},ABCI_files,{'\'},strcat(fname,'.abc'),{' '})))

            abci_path = os.getcwd()
            print("ABCI path:: ", abci_path)
            # print(abci_path)
            exe_path = os.path.join(abci_path, 'ABCI_exe\ABCI_MP64+.exe')
            print(exe_path)
            subprocess.call([exe_path, '{}\Cavity{}\Cavity_MROT_{}.abc'.format(self.path, self.fid, MROT)])
            print("Finished the program")


    def createFolder(self):
        path = os.getcwd()

        path = os.path.join(path, "ABCI_Data_end\Cavity{}".format(self.fid))
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)
            # return "Yes"

    def button_clicked(self, i):
        return i.text()