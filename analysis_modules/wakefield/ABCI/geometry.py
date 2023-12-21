
class Geometry:
    def __init__(self):

        self.WG_mesh = None
        self.WG_R = None
        self.WG_L = None
        self.left_beam_pipe, self.right_beam_pipe = None, None
        self.mid_cell, self.right_end_cell, self.left_end_cell = None, None, None
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M = \
            None, None, None, None, None, None, None
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R = \
            None, None, None, None, None, None, None
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L = \
            None, None, None, None, None, None, None

        self.Rbp_L, self.at_L, self.bt_L, self.c_L, self.x_L = None, None, None, None, None
        self.Rbp_R, self.at_R, self.bt_R, self.c_R, self.x_R = None, None, None, None, None

        self.Jxy, self.Jx1, self.Jx2 = None, None, None
        self.Jx1_bp, self.Jx2_bp = None, None
        self.Jy0, self.Jy1, self.Jy2, self.Jy3 = None, None, None, None
        self.Jy0_bp, self.Jy1_bp, self.Jy2_bp, self.Jy3_bp = None, None, None, None
        self.Jxy_al, self.Jxy_bp, self.Jxy_all, self.Jxy_bp_y = None, None, None, None

        self.Jxy_all_bp = None
        self.n = None
        self.cell_structure = None
        self.u = None

    def set_geom_parameters(self, n_cells, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None):
        # print(n_cells, mid_cells_par, l_end_cell_par, r_end_cell_par)

        if l_end_cell_par is None:
            l_end_cell_par = []

        # check self.u for slans and abci
        self.u = 1e-3
        self.cell_structure = 2
        self.n = n_cells

        # print("Sets value")
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M = \
            [i * self.u for i in mid_cells_par][0:7]
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L = \
            [i * self.u for i in l_end_cell_par][0:7]
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R = \
            [i * self.u for i in r_end_cell_par][0:7]

        self.c_L = 0 * self.u
        self.Rbp_L = self.Req_L
        self.at_L, self.bt_L, self.x_L = 0, 0, 0
        self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
        self.right_beam_pipe = self.left_beam_pipe

        # for some strange reasons, I reversed the order eesh!!!
        self.mid_cell = [self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M]
        self.left_end_cell = [self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L]
        self.right_end_cell = [self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R]
        # print(self.mid_cell)

        # beam pipe
        self.WG_L = 4 * self.L_M  # self.ui.dsb_Lbp_L.value()*self.u # Length of the beam pipe connecting to the cavity
        self.WG_R = self.WG_L  # Right Waveguide
        # self.L_all = self.WG_L + self.WG_R + self.L_L
        # + self.L_R + 2*(self.n - 1)*self.L_M # Total length of each cavity

        self.Rbp_L = self.ri_L * self.u
        self.at_L = 0 * self.u
        self.bt_L = 0 * self.u
        self.c_L = 0 * self.u
        self.x_L = 0 * (self.c_L + self.at_L)
        # self.x_L = 230*self.u # x_L is the length of the transition from iris to beam pipe

        self.Rbp_R = self.ri_R * self.u
        self.at_R = 0 * self.u
        self.bt_R = 0 * self.u
        self.c_R = 0 * self.u
        self.x_R = 0 * (self.c_L + self.at_L)
        # self.x_R = 230*self.u

        self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
        self.right_beam_pipe = self.left_beam_pipe

        # Slans mesh parameters
        self.WG_mesh = round(self.WG_L / 5) * self.u  # /5 for ende_type 1
        # print(self.L_M, self.WG_L, self.WG_mesh)
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

    def set_geom_parameters_flattop(self, n_cells, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None):
        # print(n_cells, mid_cells_par, l_end_cell_par, r_end_cell_par)

        if l_end_cell_par is None:
            l_end_cell_par = []

        # check self.u for slans and abci
        self.u = 1e-3
        self.cell_structure = 2
        self.n = n_cells

        # print("Sets value")
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M, self.l_M = \
            [i * self.u for i in mid_cells_par][0:8]
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L, self.l_L = \
            [i * self.u for i in l_end_cell_par][0:8]
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R, self.l_R = \
            [i * self.u for i in r_end_cell_par][0:8]

        self.c_L = 0 * self.u
        self.Rbp_L = self.Req_L
        self.at_L, self.bt_L, self.x_L = 0, 0, 0
        self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
        self.right_beam_pipe = self.left_beam_pipe

        # for some strange reasons, I reversed the order eesh!!!
        self.mid_cell = [self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M, self.l_M]
        self.left_end_cell = [self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L, self.l_L]
        self.right_end_cell = [self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R, self.l_R]
        # print(self.mid_cell)

        # beam pipe
        self.WG_L = 4 * self.L_M  # self.ui.dsb_Lbp_L.value()*self.u # Length of the beam pipe connecting to the cavity
        self.WG_R = self.WG_L  # Right Waveguide
        # self.L_all = self.WG_L + self.WG_R + self.L_L
        # + self.L_R + 2*(self.n - 1)*self.L_M # Total length of each cavity

        self.Rbp_L = self.ri_L * self.u
        self.at_L = 0 * self.u
        self.bt_L = 0 * self.u
        self.c_L = 0 * self.u
        self.x_L = 0 * (self.c_L + self.at_L)
        # self.x_L = 230*self.u # x_L is the length of the transition from iris to beam pipe

        self.Rbp_R = self.ri_R * self.u
        self.at_R = 0 * self.u
        self.bt_R = 0 * self.u
        self.c_R = 0 * self.u
        self.x_R = 0 * (self.c_L + self.at_L)
        # self.x_R = 230*self.u

        self.left_beam_pipe = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L]
        self.right_beam_pipe = self.left_beam_pipe

        # Slans mesh parameters
        self.WG_mesh = round(self.WG_L / 5) * self.u  # /5 for ende_type 1

    # def write_cst_paramters(self, fid):
    #     print("Writing parameters to file")
    #     cwd = os.getcwd()
    #     if self.ui.cb_Code.currentText() == 'SLANS':
    #         path = os.path.join(cwd, "node_editor\SLANS\Cavity{}\cst_parameters.txt".format(fid))
    #     else:
    #         path = os.path.join(cwd, "node_editor\ABCI\Cavity{}\cst_parameters.txt".format(fid))
    #
    #     print(path)
    #     with open(path, 'w') as f:
    #         name_list = ['c_L', 'Rbp_L', 'at_L', 'bt_L', 'x_L',
    #                      'Req_L', 'ri_L', 'L_L', 'Aeq_L', 'Beq_L', 'ai_L', 'bi_L',
    #                      'Req_M', 'ri_M', 'L_M', 'Aeq_M', 'Beq_M', 'ai_M', 'bi_M',
    #                      'Req_R', 'ri_R', 'L_R', 'Aeq_R', 'Beq_R', 'ai_R', 'bi_R',
    #                      'c_R', 'Rbp_R', 'at_R', 'bt_R', 'x_R']
    #
    #         value_list = [self.c_L, self.Rbp_L, self.at_L, self.bt_L, self.x_L,
    #                       self.Req_L, self.ri_L, self.L_L, self.A_L, self.B_L, self.a_L, self.b_L,
    #                       self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M,
    #                       self.Req_R, self.ri_R, self.L_R, self.A_R, self.B_R, self.a_R, self.b_R,
    #                       self.c_R, self.Rbp_R, self.at_R, self.bt_R, self.x_R]
    #
    #         for i in range(len(name_list)):
    #             f.write("{}={}\n".format(name_list[i], value_list[i]))
    #
    #     print("Writing to file complete.")
    #
    # def write_cst_paramters_mid(self, fid):
    #     print("Writing parameters to file")
    #     cwd = os.getcwd()
    #     print(cwd)
    #     if self.ui.cb_Code.currentText() == 'SLANS':
    #         path = os.path.join(cwd, "node_editor\SLANS\Cavity{}\cst_parameters_mid.txt".format(fid))
    #     else:
    #         path = os.path.join(cwd, "node_editor\ABCI\Cavity{}\cst_parameters_mid.txt".format(fid))
    #
    #     print(path)
    #     with open(path, 'w') as f:
    #         name_list = ['Req_M', 'ri_M', 'L_M', 'Aeq_M', 'Beq_M', 'ai_M', 'bi_M']
    #
    #         value_list = [self.Req_M, self.ri_M, self.L_M, self.A_M, self.B_M, self.a_M, self.b_M]
    #
    #         for i in range(len(name_list)):
    #             f.write("{}={}\n".format(name_list[i], value_list[i]))
    #
    #     print("Writing to file complete.")
