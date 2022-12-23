from scipy import optimize as scopt
import numpy as np


class SLANS:
    """
    This class contains functions for writing the cavity geometry in a format readable by SLANS genmesh.exe
    The type of cavity geometry it writes is regular elliptical cavity geometry with mid-cup dimensions equal
    and end cups different. Geometry type shown in the figure below:

    .. image:: ../images/Cavity_parametrised.png
      :width: 800
      :alt: Cavity parametrised model


    Note
    -----
        Considering renaming to better reflect this purpose.
    """
    def __init__(self, left_beam_pipe, left_end_cell, mid_cell, right_end_cell, right_beam_pipe, jxy_all, jxy_all_bp):
        """
        Parameters
        ----------
        left_beam_pipe : [A, B, a, b, Ri, L, Req, alpha], np.array([A, B, a, b, Ri, L, Req, alpha])
        Contains the dimensions of the left beam pipe

        left_end_cell: list
            Contains the dimensions of the left beam pipe

        mid_cell: list
        right_end_cell: list
        right_beam_pipe: list
        jxy_all:
        jxy_all_bp:
        """

        self.left_end_cell = left_end_cell
        self.mid_cell = mid_cell
        self.right_end_cell = right_end_cell
        self.left_beam_pipe = left_beam_pipe
        self.right_beam_pipe = right_beam_pipe
        self.Jxy_all = jxy_all
        self.Jxy_all_bp = jxy_all_bp

        # Mid
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M = self.mid_cell

        # Left
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L = self.left_end_cell

        # Right
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R = self.right_end_cell

        # beam pipe
        self.Rbp_L, self.at_L, self.bt_L, self.c_L, self.x_L = self.left_beam_pipe
        self.Rbp_R, self.at_R, self.bt_R, self.c_R, self.x_R = self.right_beam_pipe

    def slans_bp_L(self, n, zr12_BPL, WG_L, f):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells

        zr12_BPL: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_BPL to reflect what the variable is.
            Variable n is unused by the function. Consider removing

        Returns
        -------

        """

        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[0][0], zr12_BPL[0][1],
                                                           self.Jxy_all[2]))

        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L, self.Rbp_L - self.c_L, self.c_L,
                                                             - self.Jxy_all[5]))

        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[1][0], zr12_BPL[1][1],
                                                         - self.Jxy_all[6]))

        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L, - self.Jxy_all[7]))

        f.write(
            '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L + self.bt_L, self.bt_L,
                                                         self.Jxy_all[3]))

    def slans_n1_L(self, n, zr12_L, WG_L, f):
        """Writes the left end cup dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_L: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_R to reflect what the variable is.
            Variable n is unused by the function. Consider removing

        Returns
        -------

        """
        # print("\t\tSLANS_N1_L::It got here")

        f.write(
            '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L - zr12_L[1][0], zr12_L[1][1], self.Jxy_all[3]))
        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L, self.ri_L + self.b_L, self.b_L, self.Jxy_all[7]))
        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L - zr12_L[0][0], zr12_L[0][1], self.Jxy_all[6]))
        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L, self.Jxy_all[5]))
        f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_L, self.B_L,
                                                             self.Jxy_all[2]))

    def slans_n1_R(self, n, zr12_R, WG_L, f):
        """Writes the right end cup dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_R: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_R to reflect what the variable is.

        Returns
        -------

        """

        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_N1_R::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + zr12_R[0][0], zr12_R[0][1],
                                                               self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_R, self.B_R,
                                                                 - self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + zr12_R[1][0], zr12_R[1][1],
                                                             - self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.b_R, self.b_R,
                                                             self.Jxy_all[3]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[0][0],
                                                               zr12_R[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M,
                                                                 self.Req_R - self.B_R, self.B_R, - self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[1][0],
                                                             zr12_R[1][1], - self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R,
                                                               self.ri_R, -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R,
                                                                 self.ri_R + self.b_R, self.b_R, self.Jxy_all[3]))

    def slans_bp_R(self, n, zr12_BPR, WG_L, f):
        """Writes the right beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_BPR: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_BPR to reflect what the variable is.

        Returns
        -------

        """
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_BPR::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0],
                                                               zr12_BPR[2, 2], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.bt_R,
                                                                 self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0],
                                                             zr12_BPR[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R, self.Rbp_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R,
                                                                 self.Rbp_R - self.c_R, self.c_R, self.Jxy_all[2]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0] + 2 * (n - 1) * self.L_M, zr12_BPR[1][1],
                self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + 2 * (n - 1) * self.L_M,
                                                                 self.ri_R + self.bt_R, self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0] + 2 * (n - 1) * self.L_M, zr12_BPR[0][1],
                self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M, self.Rbp_R, self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M, self.Rbp_R, - self.c_R,
                self.Jxy_all[2]))

    def slans_M(self, n, zr12_M, WG_L, f, i, end_type):
        """Writes the mid cup dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_M: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        i: int
            Cell index

        end_type: int
            if end_type = 1 the end HALF cell is changed for tuning.
            I don't know what this means. Never had need for it.

        Note
        -----
            Consider renaming zr12_R to reflect what the variable is.

        Returns
        -------

        """
        # print("\t\tSLANS_M::It got here")
        # Left and right Cell
        # First Half cell

        if i == 1 and end_type == 1:
            zr12_LM, alpha_LM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_LM[0][0],
                                                               zr12_LM[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M,
                                                                 self.Req_L - self.B_M, self.B_M, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_LM[1][0],
                                                             zr12_LM[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M, self.ri_M,
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[3]))

        elif i == 1 and end_type != 1:
            zr12_L, alpha_L = self.rz_conjug('left')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_L + zr12_L[0][0],
                                                               zr12_L[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_L,
                                                                 self.Req_L - self.B_L, self.B_L, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_L + zr12_L[1][0],
                                                             zr12_L[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_L, self.ri_M,
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_L,
                                                                 self.ri_M + self.b_L, self.b_L, self.Jxy_all[3]))

        if i != 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_M[0][0],
                                                               zr12_M[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M,
                                                                 self.Req_M - self.B_M, self.B_M, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_M[1][0],
                                                             zr12_M[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M, self.ri_M,
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[3]))

        # Second half cells
        if i == n - 1 and end_type == 1:

            zr12_RM, alpha_RM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_RM[1][0], zr12_RM[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_RM[0][0], zr12_RM[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M, self.Req_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M,
                                                                 self.Req_R - self.B_M, self.B_M, self.Jxy_all[2]))

        elif i == n - 1 and end_type != 1:

            zr12_R, alpha_R = self.rz_conjug('right')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R - zr12_R[1][0], zr12_R[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_R, self.b_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R - zr12_R[0][0], zr12_R[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R,
                                                               self.Req_R, self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R,
                                                                 self.Req_R - self.B_R, self.B_R, self.Jxy_all[2]))

        if i != n - 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_M[1][0], zr12_M[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_M[0][0], zr12_M[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M, self.Req_M,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M,
                                                                 self.Req_M - self.B_M, self.B_M, self.Jxy_all[2]))

    def rz_conjug(self, cell):
        """Calculate shift in geoemetric values

        Parameters
        ----------
        cell: {"left", "right", "mid"}
            Specifies the cell cup type.

        Returns
        -------
        zr12: list
            Maybe it's a shift of some sort calculated from the coordinates of the tangent of the ellipses
            comprising the cavity

        alpha: float
            Angle of the line tangent to the iris and equator ellipses.

        """

        # global data
        if cell == 'left':
            data = [self.A_L, self.B_L, self.Req_L, self.ri_L, self.L_L, self.a_L, self.b_L]
        elif cell == 'right':
            data = [self.A_R, self.B_R, self.Req_R, self.ri_R, self.L_R, self.a_R, self.b_R]
        else:
            data = [self.A_M, self.B_M, self.Req_M, self.ri_M, self.L_M, self.a_M, self.b_M]

        A = data[0]  # ellipse x
        B = data[1]  # ellipse y
        Req = data[2]  # equator radius
        ri = data[3]  # iris radius
        L = data[4]  # quarter length
        a = data[5]  # cone ellipse x
        b = data[6]  # cone ellipse y

        # OLD CODE STILL KEPT HERE IN CASE MINE RUNS INTO A PROBLEM, HAHA!
        # x1 = 0
        # y1 = ri + b
        # x2 = L
        # y2 = Req - B

        # ell_1 = np.array([[x1, ri], [x1 + a, ri+b], [x1, ri+2*b], [x1 - a, ri + b]])
        # ell_2 = np.array([[L, Req - 2 * B], [L + A, Req - B], [L, Req], [L - A, Req - B]])
        # # print("ELL_1::", ell_1, np.shape(ell_1)[-1])
        # # print("ELL_2::", ell_2, np.shape(ell_2)[-1])
        #
        # k = 0
        # # initialize x0
        # x0 = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # # print(np.shape(x0))
        #
        # for i in range(len(ell_1)):
        #     for j in range(len(ell_2)):
        #         x0[k][:] = np.concatenate([ell_1[i][:], ell_2[j][:]])
        #         k = k + 1
        #
        # x = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # import time as t
        # for k in range(len(x0)):
        #     x[k][:] = scopt.fsolve(self.func, x0[k][:])
        #
        # ind_x0 = np.zeros(shape=(np.shape(x0)[0], 1))
        # for i in range(len(x)):
        #     ind_x0[i] = (x[i][3] > x[i][1]) and (x[i][0] > x1) and (L > x[i][2]) #already adjusted index
        # # print(ind_x0)
        # # print(np.shape(ind_x0))
        #
        # # print("After ind_x0")
        # xy_cross = x[np.where(ind_x0 == 1)[0]][-1] # x1 y1 x2 y2
        # print("\tXY_CROSS::", xy_cross)
        # print(np.shape(xy_cross))
        # plot_cavity(xy_cross, data)

        data = ([0, ri + b, L, Req - B], [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([a, ri + 0.85 * b, L - A, Req - 0.85 * B]),
                                args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([0.5*a, ri + 0.5 * b, L - A, Req - 0.5 * B]),
        #                         args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess
        # print("\tXY_CROSS2::", xy_cross)

        # #
        xy_L_ell = np.zeros(shape=(4, 2))

        xy_L_ell[0][:] = xy_cross[0:2]
        xy_L_ell[1][:] = xy_cross[2:4]
        xy_L_ell[2][:] = [-xy_cross[2] + 2 * L, xy_cross[3]]
        xy_L_ell[3][:] = [-xy_cross[0] + 2 * L, xy_cross[1]]

        rz_coor = xy_L_ell
        alpha = np.arctan((xy_L_ell[3][1] - xy_L_ell[2][1]) / (xy_L_ell[3][0] - xy_L_ell[2][0])) + 180
        zr12 = [[rz_coor[2][0] - L, rz_coor[2][1]], [rz_coor[3][0] - L, rz_coor[3][1]]]

        return zr12, alpha
    #
    # def func(self, x_in):
    #     A = data[0]  # ellipse x
    #     B = data[1]  # ellipse y
    #     Req = data[2]  # equator radius
    #     ri = data[3]  # iris radius
    #     L = data[4]  # quarter length
    #     a = data[5]  # cone ellipse x
    #     b = data[6]  # cone ellipse y
    #
    #     x1 = 0
    #     y1 = ri + b
    #     x2 = L
    #     y2 = Req - B
    #
    #     # F(1) = (b^2/a^2)*cot(x(1))-(B^2/A^2)*cot(x(2))
    #     # F(2) = B*sin(x(2))+y2-b*sin(x(1))-y1+(b^2/a^2)*cot(x(1))*(A*cos(x(2))+x2-(a*cos(x(1))+x1))
    #
    #     f1 = (x_in[0] - x1) ** 2 / (a ** 2) + (x_in[1] - y1) ** 2 / b ** 2 - 1  # already adjusted indices for python
    #     f2 = (x_in[2] - x2) ** 2 / (A ** 2) + (x_in[3] - y2) ** 2 / B ** 2 - 1
    #     f3 = (x_in[2] - x_in[0]) * (x_in[0] - x1) / (a ** 2) + (x_in[3] - x_in[1]) * (x_in[1] - y1) / b ** 2
    #     f4 = (x_in[2] - x_in[0]) * (x_in[2] - x2) / (A ** 2) + (x_in[3] - x_in[1]) * (x_in[3] - y2) / B ** 2
    #
    #     return [f1, f2, f3, f4]

    @staticmethod
    def ellipse_tangent(z, *data):
        """

        Parameters
        ----------
        z:

        data: array like
            First row contains the coordinates of the midpoints of the iris and equator ellipses
            Second row contains the dimensions of the ellipses


        Returns
        -------

        Notes
        -----
            Since this function is the same as that in shared function, consider removing it and importing shared
            function.

        """
        # print(data)
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return [f1, f2, f3, f4]


class SLANS_Flattop:
    """
    This class contains functions for writing the cavity geometry in a format readable by SLANS genmesh.exe
    The type of cavity geometry it writes is regular elliptical cavity geometry with mid-cup dimensions equal
    and end cups different. Geometry type shown in the figure below:

    .. image:: ../images/Cavity_parametrised.png
      :width: 800
      :alt: Cavity parametrised model


    Note
    -----
        Considering renaming to better reflect this purpose.
    """
    def __init__(self, left_beam_pipe, left_end_cell, mid_cell, right_end_cell, right_beam_pipe, jxy_all, jxy_all_bp):
        """
        Parameters
        ----------
        left_beam_pipe : [A, B, a, b, Ri, L, Req, alpha], np.array([A, B, a, b, Ri, L, Req, alpha])
        Contains the dimensions of the left beam pipe

        left_end_cell: list
            Contains the dimensions of the left beam pipe

        mid_cell: list
        right_end_cell: list
        right_beam_pipe: list
        jxy_all:
        jxy_all_bp:
        """

        self.left_end_cell = left_end_cell
        self.mid_cell = mid_cell
        self.right_end_cell = right_end_cell
        self.left_beam_pipe = left_beam_pipe
        self.right_beam_pipe = right_beam_pipe
        self.Jxy_all = jxy_all
        self.Jxy_all_bp = jxy_all_bp

        # Mid
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M = self.mid_cell

        # Left
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L = self.left_end_cell

        # Right
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R = self.right_end_cell

        # beam pipe
        self.Rbp_L, self.at_L, self.bt_L, self.c_L, self.x_L = self.left_beam_pipe
        self.Rbp_R, self.at_R, self.bt_R, self.c_R, self.x_R = self.right_beam_pipe

    def slans_bp_L(self, n, zr12_BPL, WG_L, f):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells

        zr12_BPL: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_BPL to reflect what the variable is.
            Variable n is unused by the function. Consider removing

        Returns
        -------

        """

        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[0][0], zr12_BPL[0][1],
                                                           self.Jxy_all[2]))

        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L, self.Rbp_L - self.c_L, self.c_L,
                                                             - self.Jxy_all[5]))

        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[1][0], zr12_BPL[1][1],
                                                         - self.Jxy_all[6]))

        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L, - self.Jxy_all[7]))

        f.write(
            '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L + self.bt_L, self.bt_L,
                                                         self.Jxy_all[3]))

    def slans_n1_L(self, n, zr12_L, WG_L, f):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_L: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_R to reflect what the variable is.
            Variable n is unused by the function. Consider removing

        Returns
        -------

        """
        # print("\t\tSLANS_N1_L::It got here")

        f.write(
            '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L - zr12_L[1][0], zr12_L[1][1], self.Jxy_all[3]))
        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L, self.ri_L + self.b_L, self.b_L, self.Jxy_all[7]))
        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L - zr12_L[0][0], zr12_L[0][1], self.Jxy_all[6]))
        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L, self.Jxy_all[5]))
        f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_L, self.B_L,
                                                             self.Jxy_all[2]))

    def slans_n1_R(self, n, zr12_R, WG_L, f):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_R: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_R to reflect what the variable is.

        Returns
        -------

        """

        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_N1_R::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + zr12_R[0][0], zr12_R[0][1],
                                                               self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_R, self.B_R,
                                                                 - self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + zr12_R[1][0], zr12_R[1][1],
                                                             - self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.b_R, self.b_R,
                                                             self.Jxy_all[3]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[0][0],
                                                               zr12_R[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M,
                                                                 self.Req_R - self.B_R, self.B_R, - self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[1][0],
                                                             zr12_R[1][1], - self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R,
                                                               self.ri_R, -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R,
                                                                 self.ri_R + self.b_R, self.b_R, self.Jxy_all[3]))

    def slans_bp_R(self, n, zr12_BPR, WG_L, f):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_BPR: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        Note
        -----
            Consider renaming zr12_BPR to reflect what the variable is.

        Returns
        -------

        """
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_BPR::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0],
                                                               zr12_BPR[2, 2], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.bt_R,
                                                                 self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0],
                                                             zr12_BPR[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R, self.Rbp_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R,
                                                                 self.Rbp_R - self.c_R, self.c_R, self.Jxy_all[2]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0] + 2 * (n - 1) * self.L_M, zr12_BPR[1][1],
                self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + 2 * (n - 1) * self.L_M,
                                                                 self.ri_R + self.bt_R, self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0] + 2 * (n - 1) * self.L_M, zr12_BPR[0][1],
                self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M, self.Rbp_R, self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M, self.Rbp_R, - self.c_R,
                self.Jxy_all[2]))

    def slans_M(self, n, zr12_M, WG_L, f, i, end_type):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells.

        zr12_M: float

        WG_L: float
            Length of left beam pipe

        f: file
            Geometry <filename>.geo file to be written to
        i: int
            Cell index

        end_type: int
            if end_type = 1 the end HALF cell is changed for tuning.
            I don't know what this means. Never had need for it.

        Note
        -----
            Consider renaming zr12_R to reflect what the variable is.

        Returns
        -------

        """
        # print("\t\tSLANS_M::It got here")
        # Left and right Cell
        # First Half cell

        if i == 1 and end_type == 1:
            zr12_LM, alpha_LM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_LM[0][0],
                                                               zr12_LM[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M,
                                                                 self.Req_L - self.B_M, self.B_M, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_LM[1][0],
                                                             zr12_LM[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M, self.ri_M,
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[3]))

        elif i == 1 and end_type != 1:
            zr12_L, alpha_L = self.rz_conjug('left')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_L + zr12_L[0][0],
                                                               zr12_L[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_L,
                                                                 self.Req_L - self.B_L, self.B_L, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_L + zr12_L[1][0],
                                                             zr12_L[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_L, self.ri_M,
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_L,
                                                                 self.ri_M + self.b_L, self.b_L, self.Jxy_all[3]))

        if i != 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_M[0][0],
                                                               zr12_M[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M,
                                                                 self.Req_M - self.B_M, self.B_M, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_M[1][0],
                                                             zr12_M[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M, self.ri_M,
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[3]))

        # Second half cells
        if i == n - 1 and end_type == 1:

            zr12_RM, alpha_RM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_RM[1][0], zr12_RM[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_RM[0][0], zr12_RM[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M, self.Req_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M,
                                                                 self.Req_R - self.B_M, self.B_M, self.Jxy_all[2]))

        elif i == n - 1 and end_type != 1:

            zr12_R, alpha_R = self.rz_conjug('right')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R - zr12_R[1][0], zr12_R[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_R, self.b_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R - zr12_R[0][0], zr12_R[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R,
                                                               self.Req_R, self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R,
                                                                 self.Req_R - self.B_R, self.B_R, self.Jxy_all[2]))

        if i != n - 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_M[1][0], zr12_M[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_M,
                                                                 self.ri_M + self.b_M, self.b_M, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_M[0][0], zr12_M[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M, self.Req_M,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_M,
                                                                 self.Req_M - self.B_M, self.B_M, self.Jxy_all[2]))

    def rz_conjug(self, cell):
        """

        Parameters
        ----------
        cell: {"left", "right", "mid"}
            Specifies the cell cup type.

        Returns
        -------
        zr12: list
            Maybe it's a shift of some sort calculated from the coordinates of the tangent of the ellipses
            comprising the cavity

        alpha: float
            Angle of the line tangent to the iris and equator ellipses.

        """

        # global data
        if cell == 'left':
            data = [self.A_L, self.B_L, self.Req_L, self.ri_L, self.L_L, self.a_L, self.b_L]
        elif cell == 'right':
            data = [self.A_R, self.B_R, self.Req_R, self.ri_R, self.L_R, self.a_R, self.b_R]
        else:
            data = [self.A_M, self.B_M, self.Req_M, self.ri_M, self.L_M, self.a_M, self.b_M]

        A = data[0]  # ellipse x
        B = data[1]  # ellipse y
        Req = data[2]  # equator radius
        ri = data[3]  # iris radius
        L = data[4]  # quarter length
        a = data[5]  # cone ellipse x
        b = data[6]  # cone ellipse y

        # OLD CODE STILL KEPT HERE IN CASE MINE RUNS INTO A PROBLEM, HAHA!
        # x1 = 0
        # y1 = ri + b
        # x2 = L
        # y2 = Req - B

        # ell_1 = np.array([[x1, ri], [x1 + a, ri+b], [x1, ri+2*b], [x1 - a, ri + b]])
        # ell_2 = np.array([[L, Req - 2 * B], [L + A, Req - B], [L, Req], [L - A, Req - B]])
        # # print("ELL_1::", ell_1, np.shape(ell_1)[-1])
        # # print("ELL_2::", ell_2, np.shape(ell_2)[-1])
        #
        # k = 0
        # # initialize x0
        # x0 = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # # print(np.shape(x0))
        #
        # for i in range(len(ell_1)):
        #     for j in range(len(ell_2)):
        #         x0[k][:] = np.concatenate([ell_1[i][:], ell_2[j][:]])
        #         k = k + 1
        #
        # x = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # import time as t
        # for k in range(len(x0)):
        #     x[k][:] = scopt.fsolve(self.func, x0[k][:])
        #
        # ind_x0 = np.zeros(shape=(np.shape(x0)[0], 1))
        # for i in range(len(x)):
        #     ind_x0[i] = (x[i][3] > x[i][1]) and (x[i][0] > x1) and (L > x[i][2]) #already adjusted index
        # # print(ind_x0)
        # # print(np.shape(ind_x0))
        #
        # # print("After ind_x0")
        # xy_cross = x[np.where(ind_x0 == 1)[0]][-1] # x1 y1 x2 y2
        # print("\tXY_CROSS::", xy_cross)
        # print(np.shape(xy_cross))
        # plot_cavity(xy_cross, data)

        data = ([0, ri + b, L, Req - B], [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([a, ri + 0.85 * b, L - A, Req - 0.85 * B]),
                                args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([0.5*a, ri + 0.5 * b, L - A, Req - 0.5 * B]),
        #                         args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess
        # print("\tXY_CROSS2::", xy_cross)

        # #
        xy_L_ell = np.zeros(shape=(4, 2))

        xy_L_ell[0][:] = xy_cross[0:2]
        xy_L_ell[1][:] = xy_cross[2:4]
        xy_L_ell[2][:] = [-xy_cross[2] + 2 * L, xy_cross[3]]
        xy_L_ell[3][:] = [-xy_cross[0] + 2 * L, xy_cross[1]]

        rz_coor = xy_L_ell
        alpha = np.arctan((xy_L_ell[3][1] - xy_L_ell[2][1]) / (xy_L_ell[3][0] - xy_L_ell[2][0])) + 180
        zr12 = [[rz_coor[2][0] - L, rz_coor[2][1]], [rz_coor[3][0] - L, rz_coor[3][1]]]

        return zr12, alpha
    #
    # def func(self, x_in):
    #     A = data[0]  # ellipse x
    #     B = data[1]  # ellipse y
    #     Req = data[2]  # equator radius
    #     ri = data[3]  # iris radius
    #     L = data[4]  # quarter length
    #     a = data[5]  # cone ellipse x
    #     b = data[6]  # cone ellipse y
    #
    #     x1 = 0
    #     y1 = ri + b
    #     x2 = L
    #     y2 = Req - B
    #
    #     # F(1) = (b^2/a^2)*cot(x(1))-(B^2/A^2)*cot(x(2))
    #     # F(2) = B*sin(x(2))+y2-b*sin(x(1))-y1+(b^2/a^2)*cot(x(1))*(A*cos(x(2))+x2-(a*cos(x(1))+x1))
    #
    #     f1 = (x_in[0] - x1) ** 2 / (a ** 2) + (x_in[1] - y1) ** 2 / b ** 2 - 1  # already adjusted indices for python
    #     f2 = (x_in[2] - x2) ** 2 / (A ** 2) + (x_in[3] - y2) ** 2 / B ** 2 - 1
    #     f3 = (x_in[2] - x_in[0]) * (x_in[0] - x1) / (a ** 2) + (x_in[3] - x_in[1]) * (x_in[1] - y1) / b ** 2
    #     f4 = (x_in[2] - x_in[0]) * (x_in[2] - x2) / (A ** 2) + (x_in[3] - x_in[1]) * (x_in[3] - y2) / B ** 2
    #
    #     return [f1, f2, f3, f4]

    @staticmethod
    def ellipse_tangent(z, *data):
        # print(data)
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return [f1, f2, f3, f4]


class SLANS_Multicell:
    def __init__(self, left_beam_pipe, left_end_cell, mid_cell, right_end_cell, right_beam_pipe, jxy_all,
                 jxy_all_bp):
        self.left_end_cell = left_end_cell
        self.mid_cell = mid_cell
        self.right_end_cell = right_end_cell
        self.left_beam_pipe = left_beam_pipe
        self.right_beam_pipe = right_beam_pipe
        self.Jxy_all = jxy_all
        self.Jxy_all_bp = jxy_all_bp

        # Mid
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M, _ = self.mid_cell

        # Left
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L, _ = self.left_end_cell

        # Right
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R, _ = self.right_end_cell

        # beam pipe
        self.Rbp_L, self.at_L, self.bt_L, self.c_L, self.x_L = self.left_beam_pipe
        self.Rbp_R, self.at_R, self.bt_R, self.c_R, self.x_R = self.right_beam_pipe

    def slans_bp_L(self, n, zr12_BPL, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign

        # print("\t\tSLANS_BPL::It got here")

        f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[0][0], zr12_BPL[0][1],
                                                           self.Jxy_all[2]))

        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L, self.Rbp_L - self.c_L, self.c_L,
                                                             - self.Jxy_all[5]))

        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[1][0], zr12_BPL[1][1],
                                                         - self.Jxy_all[6]))

        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L, -self.Jxy_all[7]))

        f.write(
            '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L + self.bt_L,
                                                         self.bt_L,
                                                         self.Jxy_all[3]))

    def slans_n1_L(self, zr12_L, WG_L, f):
        # print("\t\tSLANS_N1_L::It got here")

        f.write(
            '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L - zr12_L[1][0], zr12_L[1][1],
                                                       self.Jxy_all[3]))
        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L, self.ri_L + self.b_L, self.b_L, self.Jxy_all[7]))
        f.write(
            '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L - zr12_L[0][0], zr12_L[0][1], self.Jxy_all[6]))
        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L, self.Jxy_all[5]))
        f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_L, self.B_L,
                                                             self.Jxy_all[2]))

    def slans_n1_R(self, n, zr12_R, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_N1_R::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + zr12_R[0][0], zr12_R[0][1],
                                                               self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_R, self.B_R,
                                                                 - self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + zr12_R[1][0], zr12_R[1][1],
                                                             - self.Jxy_all[6]))
            f.write(
                '6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.b_R,
                                                             self.b_R,
                                                             self.Jxy_all[3]))

        if n > 1:
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M) + zr12_R[0][0],
                                                           zr12_R[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M),
                                                                 self.Req_R - self.B_R, self.B_R, -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M) + zr12_R[1][0],
                                                         zr12_R[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M) + self.L_R,
                                                               self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M) + self.L_R,
                                                             self.ri_R + self.b_R, self.b_R, self.Jxy_all[3]))

    def slans_bp_R(self, n, zr12_BPR, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_BPR::It got here")
        if n == 1:
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0],
                                                           zr12_BPR[2, 2], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.bt_R,
                                                                 self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0],
                                                             zr12_BPR[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R, self.Rbp_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R,
                                                                 self.Rbp_R - self.c_R, self.c_R, self.Jxy_all[2]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0] + 2 * sum(self.L_M), zr12_BPR[1][1],
                self.Jxy_all[3]))
            f.write(
                '7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + 2 * sum(self.L_M),
                                                             self.ri_R + self.bt_R, self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0] + 2 * sum(self.L_M), zr12_BPR[0][1],
                self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * sum(self.L_M), self.Rbp_R, self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * sum(self.L_M), self.Rbp_R, - self.c_R,
                self.Jxy_all[2]))

    def slans_M(self, n, zr12_M, WG_L, f, i, end_type):
        # print("\t\tSLANS_M::It got here")
        # Left and right Cell
        # First Half cell

        if i == 0 and end_type == 1:
            zr12_LM, alpha_LM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + zr12_LM[0][0],
                                                           zr12_LM[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]),
                                                                 self.Req_L - self.B_M[i], self.B_M[i], -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + zr12_LM[1][0],
                                                         zr12_LM[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i], self.ri_M[i],
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i],
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[3]))

        elif i == 0 and end_type != 1:
            zr12_L, alpha_L = self.rz_conjug('left')  # zr12_R first column is z , second column is r
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_L + zr12_L[0][0],
                                                           zr12_L[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_L,
                                                                 self.Req_L - self.B_L, self.B_L, -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i) * self.L_L + zr12_L[1][0],
                                                         zr12_L[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i + 1) * self.L_L, self.ri_M[i],
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i + 1) * self.L_L,
                                                                 self.ri_M[i] + self.b_L, self.b_L, self.Jxy_all[3]))

        if i != 0:
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + zr12_M[0][0],
                                                           zr12_M[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]),
                                                                 self.Req_M[i] - self.B_M[i], self.B_M[i], -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + zr12_M[1][0],
                                                         zr12_M[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i], self.ri_M[i],
                                                               - self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i],
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[3]))

        # Second half cells
        if i == n and end_type == 1:

            zr12_RM, alpha_RM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + 2 * sum(self.L_M[0:i+1]) - zr12_RM[1][0], zr12_RM[1][1],
                self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i],
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + 2 * sum(self.L_M[0:i+1]) - zr12_RM[0][0], zr12_RM[0][1],
                self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[1:i]), self.Req_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[1:i]),
                                                                 self.Req_R - self.B_M[i], self.B_M[i], self.Jxy_all[2]))

        elif i == n and end_type != 1:

            zr12_R, alpha_R = self.rz_conjug('right')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i] + self.L_R - zr12_R[1][0], zr12_R[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i],
                                                                 self.ri_M[i] + self.b_R, self.b_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i] + self.L_R - zr12_R[0][0], zr12_R[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i] + self.L_R,
                                                               self.Req_R, self.Jxy_all[5]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i] + self.L_R,
                                                             self.Req_R - self.B_R, self.B_R, self.Jxy_all[2]))

        if i != n:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + 2 * sum(self.L_M[0:i+1]) - zr12_M[1][0], zr12_M[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i]) + self.L_M[i],
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + 2 * sum(self.L_M[0:i+1]) - zr12_M[0][0], zr12_M[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i+1]), self.Req_M[i],
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * sum(self.L_M[0:i+1]),
                                                                 self.Req_M[i] - self.B_M[i], self.B_M[i], self.Jxy_all[2]))

    def rz_conjug(self, cell, i=0):
        # global data
        if cell == 'left':
            data = [self.A_L, self.B_L, self.Req_L, self.ri_L, self.L_L, self.a_L, self.b_L]
        elif cell == 'right':
            data = [self.A_R, self.B_R, self.Req_R, self.ri_R, self.L_R, self.a_R, self.b_R]
        else:
            data = [self.A_M[i], self.B_M[i], self.Req_M[i], self.ri_M[i], self.L_M[i], self.a_M[i], self.b_M[i]]

        A = data[0]  # ellipse x
        B = data[1]  # ellipse y
        Req = data[2]  # equator radius
        ri = data[3]  # iris radius
        L = data[4]  # quarter length
        a = data[5]  # cone ellipse x
        b = data[6]  # cone ellipse y

        # OLD CODE STILL KEPT HERE IN CASE MINE RUNS INTO A PROBLEM, HAHA!
        # x1 = 0
        # y1 = ri + b
        # x2 = L
        # y2 = Req - B

        # ell_1 = np.array([[x1, ri], [x1 + a, ri+b], [x1, ri+2*b], [x1 - a, ri + b]])
        # ell_2 = np.array([[L, Req - 2 * B], [L + A, Req - B], [L, Req], [L - A, Req - B]])
        # # print("ELL_1::", ell_1, np.shape(ell_1)[-1])
        # # print("ELL_2::", ell_2, np.shape(ell_2)[-1])
        #
        # k = 0
        # # initialize x0
        # x0 = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # # print(np.shape(x0))
        #
        # for i in range(len(ell_1)):
        #     for j in range(len(ell_2)):
        #         x0[k][:] = np.concatenate([ell_1[i][:], ell_2[j][:]])
        #         k = k + 1
        #
        # x = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # import time as t
        # for k in range(len(x0)):
        #     x[k][:] = scopt.fsolve(self.func, x0[k][:])
        #
        # ind_x0 = np.zeros(shape=(np.shape(x0)[0], 1))
        # for i in range(len(x)):
        #     ind_x0[i] = (x[i][3] > x[i][1]) and (x[i][0] > x1) and (L > x[i][2]) #already adjusted index
        # # print(ind_x0)
        # # print(np.shape(ind_x0))
        #
        # # print("After ind_x0")
        # xy_cross = x[np.where(ind_x0 == 1)[0]][-1] # x1 y1 x2 y2
        # print("\tXY_CROSS::", xy_cross)
        # print(np.shape(xy_cross))
        # plot_cavity(xy_cross, data)

        data = ([0, ri + b, L, Req - B], [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([a, ri + 0.85 * b, L - A, Req - 0.85 * B]),
                                args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([0.5*a, ri + 0.5 * b, L - A, Req - 0.5 * B]),
        #                         args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess
        # print("\tXY_CROSS2::", xy_cross)

        # #
        xy_L_ell = np.zeros(shape=(4, 2))

        xy_L_ell[0][:] = xy_cross[0:2]
        xy_L_ell[1][:] = xy_cross[2:4]
        xy_L_ell[2][:] = [-xy_cross[2] + 2 * L, xy_cross[3]]
        xy_L_ell[3][:] = [-xy_cross[0] + 2 * L, xy_cross[1]]

        rz_coor = xy_L_ell
        alpha = np.arctan((xy_L_ell[3][1] - xy_L_ell[2][1]) / (xy_L_ell[3][0] - xy_L_ell[2][0])) + 180
        zr12 = [[rz_coor[2][0] - L, rz_coor[2][1]], [rz_coor[3][0] - L, rz_coor[3][1]]]

        return zr12, alpha

    #
    # def func(self, x_in):
    #     A = data[0]  # ellipse x
    #     B = data[1]  # ellipse y
    #     Req = data[2]  # equator radius
    #     ri = data[3]  # iris radius
    #     L = data[4]  # quarter length
    #     a = data[5]  # cone ellipse x
    #     b = data[6]  # cone ellipse y
    #
    #     x1 = 0
    #     y1 = ri + b
    #     x2 = L
    #     y2 = Req - B
    #
    #     # F(1) = (b^2/a^2)*cot(x(1))-(B^2/A^2)*cot(x(2))
    #     # F(2) = B*sin(x(2))+y2-b*sin(x(1))-y1+(b^2/a^2)*cot(x(1))*(A*cos(x(2))+x2-(a*cos(x(1))+x1))
    #
    #     f1 = (x_in[0] - x1) ** 2 / (a ** 2) + (x_in[1] - y1) ** 2 / b ** 2 - 1  # already adjusted indices for python
    #     f2 = (x_in[2] - x2) ** 2 / (A ** 2) + (x_in[3] - y2) ** 2 / B ** 2 - 1
    #     f3 = (x_in[2] - x_in[0]) * (x_in[0] - x1) / (a ** 2) + (x_in[3] - x_in[1]) * (x_in[1] - y1) / b ** 2
    #     f4 = (x_in[2] - x_in[0]) * (x_in[2] - x2) / (A ** 2) + (x_in[3] - x_in[1]) * (x_in[3] - y2) / B ** 2
    #
    #     return [f1, f2, f3, f4]

    @staticmethod
    def ellipse_tangent(z, *data):
        # print(data)
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return [f1, f2, f3, f4]


class SLANS_Multicell_full:
    """
    This class contains functions for writing the cavity geometry in a format readable by SLANS genmesh.exe
    Considering renaming to better reflect this purpose

    The type of cavity geometry it writes is a fully parametrised elliptical cavity geometry with every mid-cup
    and end cup having separate dimensions.
    """
    def __init__(self, left_beam_pipe, left_end_cell, mid_cell, right_end_cell, right_beam_pipe, jxy_all,
                 jxy_all_bp):
        self.left_end_cell = left_end_cell
        self.mid_cell = mid_cell
        self.right_end_cell = right_end_cell
        self.left_beam_pipe = left_beam_pipe
        self.right_beam_pipe = right_beam_pipe
        self.Jxy_all = jxy_all
        self.Jxy_all_bp = jxy_all_bp

        # Mid
        self.A_M, self.B_M, self.a_M, self.b_M, self.ri_M, self.L_M, self.Req_M, _ = self.mid_cell

        # Left
        self.A_L, self.B_L, self.a_L, self.b_L, self.ri_L, self.L_L, self.Req_L, _ = self.left_end_cell

        # Right
        self.A_R, self.B_R, self.a_R, self.b_R, self.ri_R, self.L_R, self.Req_R, _ = self.right_end_cell

        # beam pipe
        self.Rbp_L, self.at_L, self.bt_L, self.c_L, self.x_L = self.left_beam_pipe
        self.Rbp_R, self.at_R, self.bt_R, self.c_R, self.x_R = self.right_beam_pipe

    def slans_bp_L(self, n, zr12_BPL, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign

        # print("\t\tSLANS_BPL::It got here")

        f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[0][0], zr12_BPL[0][1],
                                                           self.Jxy_all[2]))

        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L, self.Rbp_L - self.c_L, self.c_L,
                                                             - self.Jxy_all[5]))

        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[1][0], zr12_BPL[1][1],
                                                         - self.Jxy_all[6]))

        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L, -self.Jxy_all[7]))

        f.write(
            '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L + self.bt_L,
                                                         self.bt_L,
                                                         self.Jxy_all[3]))

    def slans_n1_L(self, zr12_L, WG_L, f):
        # print("\t\tSLANS_N1_L::It got here")

        f.write(
            '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L - zr12_L[1][0], zr12_L[1][1],
                                                       self.Jxy_all[3]))
        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L, self.ri_L + self.b_L, self.b_L, self.Jxy_all[7]))
        f.write(
            '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L - zr12_L[0][0], zr12_L[0][1], self.Jxy_all[6]))
        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L, self.Jxy_all[5]))
        f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_L, self.B_L,
                                                             self.Jxy_all[2]))

    def slans_n1_R(self, n, zr12_R, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_N1_R::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + zr12_R[0][0], zr12_R[0][1],
                                                               self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L, self.Req_L - self.B_R, self.B_R,
                                                                 - self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + zr12_R[1][0], zr12_R[1][1],
                                                             - self.Jxy_all[6]))
            f.write(
                '6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.b_R,
                                                             self.b_R,
                                                             self.Jxy_all[3]))

        if n > 1:
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M) + zr12_R[0][0],
                                                           zr12_R[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M),
                                                                 self.Req_R - self.B_R, self.B_R, - self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M) + zr12_R[1][0],
                                                         zr12_R[1][1], - self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M) + self.L_R,
                                                               self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M) + self.L_R,
                                                             self.ri_R + self.b_R, self.b_R, self.Jxy_all[3]))

    def slans_bp_R(self, n, zr12_BPR, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_BPR::It got here")
        if n == 1:
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0],
                                                           zr12_BPR[2, 2], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.bt_R,
                                                                 self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0],
                                                             zr12_BPR[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R, self.Rbp_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R,
                                                                 self.Rbp_R - self.c_R, self.c_R, self.Jxy_all[2]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0] + sum(self.L_M), zr12_BPR[1][1],
                self.Jxy_all[3]))
            f.write(
                '7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + sum(self.L_M),
                                                             self.ri_R + self.bt_R, self.bt_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0] + sum(self.L_M), zr12_BPR[0][1],
                self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + sum(self.L_M), self.Rbp_R, self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + sum(self.L_M), self.Rbp_R, - self.c_R,
                self.Jxy_all[2]))

    def slans_M(self, n, zr12_M, WG_L, f, i, end_type):
        # print("\t\tSLANS_M::It got here")
        # Left and right Cell
        # First Half cell

        if i == 0 and end_type == 1:
            # print("It's in here: i == 0 and end_type == 1: ", i, end_type)
            zr12_LM, alpha_LM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i]) + zr12_LM[0][0],
                                                           zr12_LM[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i]),
                                                                 self.Req_L - self.B_M[i], self.B_M[i], -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i]) + zr12_LM[1][0],
                                                         zr12_LM[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i+1]), self.ri_M[i],
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i+1]),
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[3]))

        elif i == 0 and end_type != 1:
            # print("It's in here: i == 0 and end_type != 1: ", i, end_type)
            zr12_L, alpha_L = self.rz_conjug('left')  # zr12_R first column is z , second column is r
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * (i - 1)) * self.L_L + zr12_L[0][0],
                                                           zr12_L[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * (i - 1)) * self.L_L,
                                                                 self.Req_L - self.B_L, self.B_L, -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * (i - 1)) * self.L_L + zr12_L[1][0],
                                                         zr12_L[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_L, self.ri_M[i],
                                                               -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + (2 * i - 1) * self.L_L,
                                                                 self.ri_M[i] + self.b_L, self.b_L, self.Jxy_all[3]))

        if i != 0 and i % 2 == 0:
            # print("It's in here: i != 0: ", i, end_type)
            f.write(
                '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i]) + zr12_M[0][0],
                                                           zr12_M[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i]),
                                                                 self.Req_M[i] - self.B_M[i], self.B_M[i], -self.Jxy_all[5]))
            f.write(
                '1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i]) + zr12_M[1][0],
                                                         zr12_M[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i+1]), self.ri_M[i],
                                                               - self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i+1]),
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[3]))

        # Second half cells, in this part i is never equal to n
        if i == n and end_type == 1:
            # print("It's in here: i == n and end_type == 1: ", i, end_type)

            zr12_RM, alpha_RM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + sum(self.L_M[0:2*(i+1)]) - zr12_RM[1][0], zr12_RM[1][1],
                self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i+1]),
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + sum(self.L_M[0:2*(i+1)]) - zr12_RM[0][0], zr12_RM[0][1],
                self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[2:2*i]), self.Req_R,
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[2:2*i]),
                                                                 self.Req_R - self.B_M[i], self.B_M[i], self.Jxy_all[2]))

        elif i == n and end_type != 1:
            # print("It's in here: i == n and end_type != 1: ", i, end_type)

            zr12_R, alpha_R = self.rz_conjug('right')  # zr12_R first column is z , second column is r
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + sum(self.L_M[0:2*i+1]) + self.L_R - zr12_R[1][0], zr12_R[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i+1]),
                                                                 self.ri_M[i] + self.b_R, self.b_R, self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + sum(self.L_M[0:2*i+1]) + self.L_R - zr12_R[0][0], zr12_R[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i+1]) + self.L_R,
                                                               self.Req_R, self.Jxy_all[5]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:2*i+1]) + self.L_R,
                                                             self.Req_R - self.B_R, self.B_R, self.Jxy_all[2]))

        if i != n-1 and i % 2 == 1:
            # print("It's in here: i != n: ", i, end_type)
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + sum(self.L_M[0:i+1]) - zr12_M[1][0], zr12_M[1][1], self.Jxy_all[3]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i]),
                                                                 self.ri_M[i] + self.b_M[i], self.b_M[i], self.Jxy_all[7]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + sum(self.L_M[0:i+1]) - zr12_M[0][0], zr12_M[0][1], self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i+1]), self.Req_M[i],
                                                               self.Jxy_all[5]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + sum(self.L_M[0:i+1]),
                                                                 self.Req_M[i] - self.B_M[i], self.B_M[i], self.Jxy_all[2]))

    def rz_conjug(self, cell, i=0):
        # global data
        if cell == 'left':
            data = [self.A_L, self.B_L, self.Req_L, self.ri_L, self.L_L, self.a_L, self.b_L]
        elif cell == 'right':
            data = [self.A_R, self.B_R, self.Req_R, self.ri_R, self.L_R, self.a_R, self.b_R]
        else:
            data = [self.A_M[i], self.B_M[i], self.Req_M[i], self.ri_M[i], self.L_M[i], self.a_M[i], self.b_M[i]]

        A = data[0]  # ellipse x
        B = data[1]  # ellipse y
        Req = data[2]  # equator radius
        ri = data[3]  # iris radius
        L = data[4]  # quarter length
        a = data[5]  # cone ellipse x
        b = data[6]  # cone ellipse y

        # OLD CODE STILL KEPT HERE IN CASE MINE RUNS INTO A PROBLEM, HAHA!
        # x1 = 0
        # y1 = ri + b
        # x2 = L
        # y2 = Req - B

        # ell_1 = np.array([[x1, ri], [x1 + a, ri+b], [x1, ri+2*b], [x1 - a, ri + b]])
        # ell_2 = np.array([[L, Req - 2 * B], [L + A, Req - B], [L, Req], [L - A, Req - B]])
        # # print("ELL_1::", ell_1, np.shape(ell_1)[-1])
        # # print("ELL_2::", ell_2, np.shape(ell_2)[-1])
        #
        # k = 0
        # # initialize x0
        # x0 = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # # print(np.shape(x0))
        #
        # for i in range(len(ell_1)):
        #     for j in range(len(ell_2)):
        #         x0[k][:] = np.concatenate([ell_1[i][:], ell_2[j][:]])
        #         k = k + 1
        #
        # x = np.zeros(shape=(np.shape(ell_1)[0] * np.shape(ell_2)[0], np.shape(ell_1)[-1] + np.shape(ell_2)[-1]))
        # import time as t
        # for k in range(len(x0)):
        #     x[k][:] = scopt.fsolve(self.func, x0[k][:])
        #
        # ind_x0 = np.zeros(shape=(np.shape(x0)[0], 1))
        # for i in range(len(x)):
        #     ind_x0[i] = (x[i][3] > x[i][1]) and (x[i][0] > x1) and (L > x[i][2]) #already adjusted index
        # # print(ind_x0)
        # # print(np.shape(ind_x0))
        #
        # # print("After ind_x0")
        # xy_cross = x[np.where(ind_x0 == 1)[0]][-1] # x1 y1 x2 y2
        # print("\tXY_CROSS::", xy_cross)
        # print(np.shape(xy_cross))
        # plot_cavity(xy_cross, data)

        data = ([0, ri + b, L, Req - B], [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([a, ri + 0.85 * b, L - A, Req - 0.85 * B]),
                                args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # xy_cross = scopt.fsolve(self.ellipse_tangent, np.array([0.5*a, ri + 0.5 * b, L - A, Req - 0.5 * B]),
        #                         args=data)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess
        # print("\tXY_CROSS2::", xy_cross)

        # #
        xy_L_ell = np.zeros(shape=(4, 2))

        xy_L_ell[0][:] = xy_cross[0:2]
        xy_L_ell[1][:] = xy_cross[2:4]
        xy_L_ell[2][:] = [-xy_cross[2] + 2 * L, xy_cross[3]]
        xy_L_ell[3][:] = [-xy_cross[0] + 2 * L, xy_cross[1]]

        rz_coor = xy_L_ell
        alpha = np.arctan((xy_L_ell[3][1] - xy_L_ell[2][1]) / (xy_L_ell[3][0] - xy_L_ell[2][0])) + 180
        zr12 = [[rz_coor[2][0] - L, rz_coor[2][1]], [rz_coor[3][0] - L, rz_coor[3][1]]]

        return zr12, alpha

    #
    # def func(self, x_in):
    #     A = data[0]  # ellipse x
    #     B = data[1]  # ellipse y
    #     Req = data[2]  # equator radius
    #     ri = data[3]  # iris radius
    #     L = data[4]  # quarter length
    #     a = data[5]  # cone ellipse x
    #     b = data[6]  # cone ellipse y
    #
    #     x1 = 0
    #     y1 = ri + b
    #     x2 = L
    #     y2 = Req - B
    #
    #     # F(1) = (b^2/a^2)*cot(x(1))-(B^2/A^2)*cot(x(2))
    #     # F(2) = B*sin(x(2))+y2-b*sin(x(1))-y1+(b^2/a^2)*cot(x(1))*(A*cos(x(2))+x2-(a*cos(x(1))+x1))
    #
    #     f1 = (x_in[0] - x1) ** 2 / (a ** 2) + (x_in[1] - y1) ** 2 / b ** 2 - 1  # already adjusted indices for python
    #     f2 = (x_in[2] - x2) ** 2 / (A ** 2) + (x_in[3] - y2) ** 2 / B ** 2 - 1
    #     f3 = (x_in[2] - x_in[0]) * (x_in[0] - x1) / (a ** 2) + (x_in[3] - x_in[1]) * (x_in[1] - y1) / b ** 2
    #     f4 = (x_in[2] - x_in[0]) * (x_in[2] - x2) / (A ** 2) + (x_in[3] - x_in[1]) * (x_in[3] - y2) / B ** 2
    #
    #     return [f1, f2, f3, f4]

    @staticmethod
    def ellipse_tangent(z, *data):
        # print(data)
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return [f1, f2, f3, f4]
