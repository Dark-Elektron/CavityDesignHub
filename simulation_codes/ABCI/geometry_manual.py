import os
from termcolor import colored
import numpy as np

file_color = 'magenta'

def print_(*arg):
    print(colored(f'\t\t\t\t\t{arg}', file_color))

class Geometry:
    def __init__(self, win=None):
        if win:
            self.win = win
            self.ui = win.ui

            # if self.ui.cb_Code.currentText() == 'SLANS':
            #     self.u = 1
            # elif self.ui.cb_Code.currentText() == 'ABCI':
            #     self.u = 1e-3

    def set_geom_parameters(self, n_cells, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None, beampipes=None):
        print_(n_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)
        # # ##### uncomment for gui
        # self.n = self.ui.sb_N_Cells.value()
        # self.cell_structure = self.ui.cb_Only_Mid_Cells.checkState()
        # print_("Structure:: ", self.cell_structure)
        # self.unit = self.ui.cb_Unit.currentText()
        #
        # if self.ui.cb_Code.currentText() == 'SLANS':
        #     self.u = 1
        # elif self.ui.cb_Code.currentText() == 'ABCI':
        #     self.u = 1e-3




        # # fix this later so you don't have to be changinf units upandan for abci and slans
        # if self.unit == 'mm':
        #     self.u = 1e-3
        # elif self.unit == 'cm':
        #     self.u = 1e-2
        # elif self.unit == 'm':
        #     self.u = 1
        # ######

        # A_M, B_M, a_M, b_M, Req_M, Ri_M, L_M = mid_cells_par
        # A_L, B_L, a_L, b_L, Req_L, Ri_L, L_L = l_end_cell_par
        # A_R, B_R, a_R, b_R, Req_R, Ri_R, L_R = r_end_cell_par

        if l_end_cell_par is None:
            l_end_cell_par = []

        # Cavity Parameters
        # this part is not needed for our parametrization. Initial guess parameters?

        if mid_cells_par:
            # check self.u for slans and abci
            self.u = 1
            self.cell_structure = 2
            self.n = n_cells

            # print_("Sets value")
            self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M = [i*self.u for i in mid_cells_par]
            self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L = [i*self.u for i in l_end_cell_par]
            self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R = [i*self.u for i in r_end_cell_par]

            self.WG_L, self.WG_R = 0, 0
            if beampipes.lower() == 'both':
                self.WG_L = 4 * self.L_M
                self.WG_R = 4 * self.L_M

            if beampipes == 'left':
                self.WG_L = 4*self.L_M # self.ui.dsb_Lbp_L.value()*self.u   # Length of the beam pipe connecting to the cavity

            if beampipes == 'right':
                self.WG_R = 4*self.L_M # Right Waveguide

            self.c_L = 0 * self.u
            self.Rbp_L = self.Req_L
            self.at_L, self.bt_L, self.x_L = 0, 0, 0
            self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
            self.right_beam_pipe = self.left_beam_pipe

            # for some strange reasons, I reversed the order eesh!!!
            self.mid_cell = [self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M]
            self.left_end_cell = [self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L]
            self.right_end_cell = [self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R]

            # print_(self.mid_cell)

        else:
            # Middle Ellipse data
            self.Req_M = self.ui.dsb_Req_M.value()*self.u  # Middle half cell equator radius
            self.ri_M = self.ui.dsb_Ri_M.value()*self.u  # Middle half cell iris radius
            self.L_M = self.ui.dsb_L_M.value()*self.u  # Middle half cell quarter length
            self.A_M = self.ui.dsb_Aeq_M.value()*self.u  # Middle half cell ellipse x
            self.B_M = self.ui.dsb_Beq_M.value()*self.u  # Middle half cell ellipse y
            self.a_M = self.ui.dsb_ai_M.value()*self.u  # Middle half cell cone ellipse x
            self.b_M = self.ui.dsb_bi_M.value()*self.u  # Middle half cell cone ellipse y

            self.mid_cell = [self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M]

            # print_(self.mid_cell)

            if self.cell_structure == 2:
                print_("Only mid cells")

                # First Ellipse data
                self.Req_L = self.Req_M
                self.ri_L = self.ri_M
                self.L_L = self.L_M
                self.A_L = self.A_M
                self.B_L = self.B_M
                self.a_L = self.a_M
                self.b_L = self.b_M

                self.left_end_cell = [self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L]

                # Right Ellipse data
                self.Req_R = self.Req_M
                self.ri_R = self.ri_M
                self.L_R = self.L_M
                self.A_R = self.A_M
                self.B_R = self.B_M
                self.a_R = self.a_M
                self.b_R = self.b_M

                self.right_end_cell = [self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R]

                # Beam Pipe
                self.WG_L = 0  # self.ui.dsb_Lbp_L.value()*self.u   # Length of the beam pipe connecting to the cavity
                self.WG_R = self.WG_L # Right Waveguide
                self.L_all = self.WG_L + self.WG_R + self.L_L + self.L_R + 2*(self.n - 1)*self.L_M  # Total length of each cavity

                self.Rbp_L = 0 * self.u
                self.at_L = 0 * self.u
                self.bt_L = 0 * self.u
                self.c_L = 0 * self.u
                self.x_L = 0 * (self.c_L + self.at_L)
                self.x_L = 0 * self.u  # x_L is the length of the transition from iris to beam pipe

                self.Rbp_R = 0 * self.u
                self.at_R = 0 * self.u
                self.bt_R = 0 * self.u
                self.c_R = 0 * self.u
                self.x_R = 0 * (self.c_R + self.at_R)
                self.x_R = 0 * self.u

                self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
                self.right_beam_pipe = [self.c_R, self.Rbp_R, self.at_R, self.bt_R, self.x_R]

            else:
                # First Ellipse data
                self.Req_L = self.ui.dsb_Req_L.value()*self.u  # Left half cell equator radius
                self.ri_L = self.ui.dsb_Ri_L.value()*self.u  # Left half cell iris radius
                self.L_L = self.ui.dsb_L_L.value()*self.u  # Left half cell quarter length
                self.A_L = self.ui.dsb_Aeq_L.value()*self.u  # Left half cell ellipse x
                self.B_L = self.ui.dsb_Beq_L.value()*self.u  # Left half cell ellipse y
                self.a_L = self.ui.dsb_ai_L.value()*self.u  # Left half cell cone ellipse x
                self.b_L = self.ui.dsb_bi_L.value()*self.u  # Left half cell cone ellipse y

                self.left_end_cell = [self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L]

                # Right Ellipse data
                self.Req_R = self.ui.dsb_Req_R.value()*self.u  # Right half cell equator radius
                self.ri_R = self.ui.dsb_Ri_R.value()*self.u # Right half cell iris radius
                self.L_R = self.ui.dsb_L_R.value()*self.u # Right half cell quarter length
                self.A_R = self.ui.dsb_Aeq_R.value()*self.u # Right half cell ellipse x
                self.B_R = self.ui.dsb_Beq_R.value()*self.u # Right half cell ellipse y
                self.a_R = self.ui.dsb_ai_R.value()*self.u # Right half cell cone ellipse x
                self.b_R = self.ui.dsb_bi_R.value()*self.u # Right half cell cone ellipse y

                self.right_end_cell = [self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R]

                # beam pipe
                self.WG_L = 4*self.L_M # self.ui.dsb_Lbp_L.value()*self.u   # Length of the beam pipe connecting to the cavity
                self.WG_R = self.WG_L # Right Waveguide
                self.L_all = self.WG_L + self.WG_R + self.L_L + self.L_R + 2*(self.n - 1)*self.L_M # Total length of each cavity

                self.Rbp_L = self.ri_L
                self.at_L = 0*self.u
                self.bt_L = 0*self.u
                self.c_L = 0*self.u
                self.x_L = 0*(self.c_L + self.at_L)
                # self.x_L = 230*self.u # x_L is the length of the transition from iris to beam pipe

                self.Rbp_R = self.ri_R
                self.at_R = 0*self.u
                self.bt_R = 0*self.u
                self.c_R = 0*self.u
                self.x_R = 0*(self.c_L + self.at_L)
                # self.x_R = 230*self.u

                self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
                self.right_beam_pipe = self.left_beam_pipe

        # Slans mesh parameters
        self.WG_mesh = round(self.WG_L / 5)*self.u # /5 for ende_type 1
        # print_(self.L_M, self.WG_L, self.WG_mesh)
        self.Jxy = 44  # 60 for end type 1
        self.Jx1 = round((19 / 50) * self.Jxy)  # 19/50 for end_type 1
        self.Jx2 = self.Jxy / 2 - self.Jx1

        self.Jy0 = round((19 / 50) * self.Jxy)  # 18, 12,14 for end type 1
        self.Jy1 = round((12 / 50) * self.Jxy)
        self.Jy2 = round((13 / 50) * self.Jxy)
        self.Jy3 = self.Jxy - self.Jy2 - self.Jy1 - self.Jy0
        self.Jxy_all = [0, self.Jxy, self.Jx1, self.Jx2, self.Jy0, self.Jy1, self.Jy2, self.Jy3]

        self.Jxy_bp = 30
        self.Jxy_bp_y = 35
        self.Jx1_bp = round((19 / 50) * self.Jxy_bp)
        self.Jx2_bp = self.Jxy_bp / 2 - self.Jx1_bp

        self.Jy0_bp = round((18 / 50) * self.Jxy_bp_y)
        self.Jy1_bp = round((12 / 50) * self.Jxy_bp_y)
        self.Jy2_bp = round((14 / 50) * self.Jxy_bp_y)
        self.Jy3_bp = self.Jxy_bp_y - self.Jy2_bp - self.Jy1_bp - self.Jy0_bp
        self.Jxy_all_bp = [0, self.Jxy_bp, self.Jx1_bp, self.Jx2_bp, self.Jy0_bp, self.Jy1_bp, self.Jy2_bp, self.Jy3_bp]

    def write_cst_paramters(self, fid):
        print_("Writing parameters to file")
        cwd = os.getcwd()
        if self.ui.cb_Code.currentText() == 'SLANS':
            path = os.path.join(cwd, "Data\SLANS\Cavity{}\cst_parameters.txt".format(fid))
        else:
            path = os.path.join(cwd, "Data\ABCI\Cavity{}\cst_parameters.txt".format(fid))

        print_(path)
        with open(path, 'w') as f:
            name_list = ['c_L', 'Rbp_L', 'at_L', 'bt_L', 'x_L',
                         'Req_L', 'ri_L', 'L_L', 'Aeq_L', 'Beq_L', 'ai_L', 'bi_L',
                         'Req_M', 'ri_M', 'L_M', 'Aeq_M', 'Beq_M', 'ai_M', 'bi_M',
                         'Req_R', 'ri_R', 'L_R', 'Aeq_R', 'Beq_R', 'ai_R', 'bi_R',
                         'c_R', 'Rbp_R', 'at_R', 'bt_R', 'x_R']

            value_list = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L,
                          self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L,
                          self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M,
                          self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R,
                          self.c_R, self.Rbp_R, self.at_R, self.bt_R, self.x_R]

            for i in range(len(name_list)):
                f.write("{}={}\n".format(name_list[i], value_list[i]))

        print_("Writing to file complete.")

    def write_cst_paramters_mid(self, fid):
        print_("Writing parameters to file")
        cwd = os.getcwd()
        print_(cwd)
        if self.ui.cb_Code.currentText() == 'SLANS':
            path = os.path.join(cwd, "Data\SLANS\Cavity{}\cst_parameters_mid.txt".format(fid))
        else:
            path = os.path.join(cwd, "Data\ABCI\Cavity{}\cst_parameters_mid.txt".format(fid))

        print_(path)
        with open(path, 'w') as f:
            name_list = ['Req_M', 'ri_M', 'L_M', 'Aeq_M', 'Beq_M', 'ai_M', 'bi_M']

            value_list = [self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M]

            for i in range(len(name_list)):
                f.write("{}={}\n".format(name_list[i], value_list[i]))

        print_("Writing to file complete.")