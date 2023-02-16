import scipy as scp
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from scipy import optimize as scopt
import numpy as np
import sympy as sym

from utils.shared_functions import tangent_coords


class SLANS:
    def __init__(self, left_beam_pipe, left_end_cell, mid_cell, right_end_cell, right_beam_pipe, jxy_all, jxy_all_bp):
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
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign

        # print("\t\tSLANS_BPL::It got here")

        f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[0][0], zr12_BPL[0][1], self.Jxy_all_bp[2]))
        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L, self.Rbp_L - self.c_L, self.c_L, -self.Jxy_all_bp[5]))
        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + zr12_BPL[1][0], zr12_BPL[1][1], -self.Jxy_all_bp[6]))
        f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L, -self.Jxy_all_bp[7]))
        f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L - self.x_L + self.x_L, self.ri_L + self.bt_L, self.bt_L, self.Jxy_all_bp[3]))

    def slans_n1_L(self, n, zr12_L, WG_L, f):
        # print("\t\tSLANS_N1_L::It got here")

        f.write(
            '6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L - zr12_L[1][0], zr12_L[1][1], self.Jxy_all[3]))
        f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L, self.ri_L + self.b_L, self.b_L, self.Jxy_all[7]))
        f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L - zr12_L[0][0], zr12_L[0][1], self.Jxy_all[6]))
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
                                                                 -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + zr12_R[1][0], zr12_R[1][1],
                                                             -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R, -self.Jxy_all[7]))
            f.write(
                '7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.b_R, self.b_R,
                                                             self.Jxy_all[3]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[0][0],
                                                               zr12_R[0][1], self.Jxy_all[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M,
                                                                 self.Req_R - self.B_R, self.B_R, -self.Jxy_all[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[1][0],
                                                             zr12_R[1][1], -self.Jxy_all[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R,
                                                               self.ri_R, -self.Jxy_all[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R,
                                                                 self.ri_R + self.b_R, self.b_R, self.Jxy_all[3]))

    def slans_bp_R(self, n, zr12_BPR, WG_L, f):
        # N1 Z R Alfa Mesh_thick Jx Jy BC_sign Vol_sign
        # print("\t\tSLANS_BPR::It got here")
        if n == 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0],
                                                               zr12_BPR[1][1], self.Jxy_all_bp[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R, self.ri_R + self.bt_R,
                                                                 self.bt_R, self.Jxy_all_bp[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0],
                                                             zr12_BPR[0][1], self.Jxy_all_bp[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R, self.Rbp_R,
                                                               self.Jxy_all_bp[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(WG_L + self.L_L + self.L_R + self.x_R,
                                                                 self.Rbp_R - self.c_R, self.c_R, self.Jxy_all_bp[3]))

        if n > 1:
            f.write('6 {:g} {:g} 0.5 1 {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0] + 2 * (n - 1) * self.L_M, zr12_BPR[1][1],
                self.Jxy_all_bp[2]))
            f.write('7 {:g} {:g} 90 {:g} 0 {:.0f} 5 0 \n'.format(WG_L + self.L_L + self.L_R + 2 * (n - 1) * self.L_M,
                                                                 self.ri_R + self.bt_R, self.bt_R, self.Jxy_all_bp[5]))
            f.write('1 {:g} {:g} 0 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0] + 2 * (n - 1) * self.L_M, zr12_BPR[0][1],
                self.Jxy_all_bp[6]))
            f.write('6 {:g} {:g} 0.5 1 0 {:.0f} 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M, self.Rbp_R, self.Jxy_all_bp[7]))
            f.write('7 {:g} {:g} 90 {:g} {:.0f} 0 5 0 \n'.format(
                WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M, self.Rbp_R - self.c_R, self.c_R,
                self.Jxy_all_bp[3]))

    def slans_M(self, n, zr12_M, WG_L, f, i, end_type):
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
        elif cell == 'mid':
            data = [self.A_M, self.B_M, self.Req_M, self.ri_M, self.L_M, self.a_M, self.b_M]
        elif cell == 'expansion':
            data = [self.c_L, self.c_L, self.Rbp_L, self.ri_L, self.x_L, self.at_L, self.bt_L]
        elif cell == 'expansion_r':
            data = [self.c_R, self.c_R, self.Rbp_R, self.ri_R, self.x_R, self.at_R, self.bt_R]

        A = data[0]  # ellipse x
        B = data[1]  # ellipse y
        Req = data[2]  # equator radius
        ri = data[3]  # iris radius
        L = data[4]  # quarter length
        a = data[5]  # cone ellipse x
        b = data[6]  # cone ellipse y

        df = tangent_coords(A, B, a, b, ri, L, Req, 0)
        x1, y1, x2, y2 = df[0]
        alpha = 180 - np.arctan2(y2 - y1, (x2 - x1)) * 180 / np.pi

        xy_cross = np.array([x1, y1, x2, y2])
        xy_L_ell = np.zeros(shape=(4, 2))

        xy_L_ell[0][:] = xy_cross[0:2]
        xy_L_ell[1][:] = xy_cross[2:4]
        xy_L_ell[2][:] = [-xy_cross[2] + 2 * L, xy_cross[3]]
        xy_L_ell[3][:] = [-xy_cross[0] + 2 * L, xy_cross[1]]

        rz_coor = xy_L_ell
        # alpha = 180 - np.arctan2(y2 - y1, (x2 - x1)) * 180 / np.pi
        zr12 = [[rz_coor[2][0] - L, rz_coor[2][1]], [rz_coor[3][0] - L, rz_coor[3][1]]]

        print(zr12, alpha)
        return zr12, alpha


def func(z, *data):
    """
    Calculates the coordinates of the tangent line that connects two ellipses

    .. _ellipse tangent:

    .. figure:: ../images/ellipse_tangent.png
       :alt: ellipse tangent
       :align: center
       :width: 200px

    Parameters
    ----------
    z: list, array like
        Contains list of tangent points coordinate's variables ``[x1, y1, x2, y2]``.
        See :numref:`ellipse tangent`
    data: list, array like
        Contains midpoint coordinates of the two ellipses and the dimensions of the ellipses
        data = ``[coords, dim]``; ``coords`` = ``[h, k, p, q]``, ``dim`` = ``[a, b, A, B]``


    Returns
    -------
    list of four non-linear functions

    Note
    -----
    The four returned non-linear functions are

    .. math::

       f_1 = \\frac{A^2b^2(x_1 - h)(y_2-q)}{a^2B^2(x_2-p)(y_1-k)} - 1

       f_2 = \\frac{(x_1 - h)^2}{a^2} + \\frac{(y_1-k)^2}{b^2} - 1

       f_3 = \\frac{(x_2 - p)^2}{A^2} + \\frac{(y_2-q)^2}{B^2} - 1

       f_4 = \\frac{-b^2(x_1-x_2)(x_1-h)}{a^2(y_1-y_2)(y_1-k)} - 1
    """

    coord, dim = data
    h, k, p, q = coord
    a, b, A, B = dim
    x1, y1, x2, y2 = z

    f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
    f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
    f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
    f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

    return f1, f2, f3, f4


def jac(z, *data):
    """
    Computes the Jacobian of the non-linear system of ellipse tangent equations

    Parameters
    ----------
    z: list, array like
        Contains list of tangent points coordinate's variables ``[x1, y1, x2, y2]``.
        See :numref:`ellipse tangent`
    data: list, array like
        Contains midpoint coordinates of the two ellipses and the dimensions of the ellipses
        data = ``[coords, dim]``; ``coords`` = ``[h, k, p, q]``, ``dim`` = ``[a, b, A, B]``

    Returns
    -------
    J: array like
        Array of the Jacobian

    """
    coord, dim = data
    h, k, p, q = coord
    a, b, A, B = dim
    x1, y1, x2, y2 = z

    # f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
    # f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
    # f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
    # f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

    df1_dx1 = A ** 2 * b ** 2 * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k))
    df1_dy1 = - A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k) ** 2)
    df1_dx2 = - A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) ** 2 * (y1 - k))
    df1_dy2 = A ** 2 * b ** 2 * (x1 - h) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k))

    df2_dx1 = 2 * (x1 - h) / a ** 2
    df2_dy1 = 2 * (y1 - k) / b ** 2
    df2_dx2 = 0
    df2_dy2 = 0

    df3_dx1 = 0
    df3_dy1 = 0
    df3_dx2 = 2 * (x2 - p) / A ** 2
    df3_dy2 = 2 * (y2 - q) / B ** 2

    df4_dx1 = -b ** 2 * ((x1 - x2) + (x1 - h)) / (a ** 2 * (y1 - y2) * (y1 - k))
    df4_dy1 = -b ** 2 * (x1 - x2) * (x1 - h) * ((y1 - y2) + (y1 - k)) / (a ** 2 * ((y1 - y2) * (y1 - k)) ** 2)
    df4_dx2 = b ** 2 * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k))
    df4_dy2 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) ** 2 * (y1 - k))

    J = [[df1_dx1, df1_dy1, df1_dx2, df1_dy2],
         [df2_dx1, df2_dy1, df2_dx2, df2_dy2],
         [df3_dx1, df3_dy1, df3_dx2, df3_dy2],
         [df4_dx1, df4_dy1, df4_dx2, df4_dy2]]

    return J
