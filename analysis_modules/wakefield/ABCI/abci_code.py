from scipy import optimize as scopt
from utils.shared_functions import *


class ABCI:
    """
    This class contains functions for writing the cavity geometry in a format readable by ABCI
    The type of cavity geometry it writes is regular elliptical cavity geometry with mid-cup dimensions equal
    and end cups different. Geometry type shown in the figure below:

    .. image:: ../images/Cavity_parametrised.png
      :width: 800
      :alt: Cavity parametrised model


    Note
    -----
        Considering renaming to better reflect this purpose.
    """
    def __init__(self, left_beam_pipe, left_end_cell, mid_cell, right_end_cell, right_beam_pipe):
        """

        Parameters
        ----------
        left_beam_pipe: list, array like

        left_end_cell: list, array like

        mid_cell: list, array like

        right_end_cell: list, array like

        right_beam_pipe: list, array like

        """
        self.left_end_cell = left_end_cell
        self.mid_cell = mid_cell
        self.right_end_cell = right_end_cell
        self.left_beam_pipe = left_beam_pipe
        self.right_beam_pipe = right_beam_pipe

        # First Ellipse data
        self.Req_L = self.left_end_cell[0]
        self.ri_L = self.left_end_cell[1]
        self.L_L = self.left_end_cell[2]
        self.A_L = self.left_end_cell[3]
        self.B_L = self.left_end_cell[4]
        self.a_L = self.left_end_cell[5]
        self.b_L = self.left_end_cell[6]

        # Right
        self.Req_R = self.right_end_cell[0]
        self.ri_R = self.right_end_cell[1]
        self.L_R = self.right_end_cell[2]
        self.A_R = self.right_end_cell[3]
        self.B_R = self.right_end_cell[4]
        self.a_R = self.right_end_cell[5]
        self.b_R = self.right_end_cell[6]

        # Middle
        self.Req_M = self.mid_cell[0]
        self.ri_M = self.mid_cell[1]
        self.L_M = self.mid_cell[2]
        self.A_M = self.mid_cell[3]
        self.B_M = self.mid_cell[4]
        self.a_M = self.mid_cell[5]
        self.b_M = self.mid_cell[6]

        # beam pipe
        self.Rbp_L = self.left_beam_pipe[0]
        self.at_L = self.left_beam_pipe[1]
        self.bt_L = self.left_beam_pipe[2]
        self.c_L = self.left_beam_pipe[3]
        self.x_L = self.left_beam_pipe[4]

        self.Rbp_R = self.right_beam_pipe[0]
        self.at_R = self.right_beam_pipe[1]
        self.bt_R = self.right_beam_pipe[2]
        self.c_R = self.right_beam_pipe[3]
        self.x_R = self.right_beam_pipe[4]

    def abci_bp_L(self, n, zr12_BPL, WG_L, f):
        """Writes the left beam pipe dimensions to a geometry file

        Parameters
        ----------
        n: int
            Number of cavity cells

        zr12_BPL: list, array like

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
        # print("\t\tABCI_BPL::It got here")
        f.write('-3., 0.000\n')
        f.write('{} {}\n'.format(self.Rbp_L - self.c_L, WG_L - self.x_L))
        f.write('{} {}\n'.format(zr12_BPL[0][1], WG_L - self.x_L + zr12_BPL[0][0]))
        f.write('{} {}\n'.format(zr12_BPL[1][1], WG_L - self.x_L + zr12_BPL[1][0]))
        f.write('-3., 0.000\n')
        f.write('{} {}\n'.format(self.ri_L + self.bt_L, WG_L - self.x_L + self.x_L))
        f.write('{} {}\n'.format(self.ri_L, WG_L - self.x_L + self.x_L))

    def abci_n1_L(self, n, zr12_L, WG_L, f):
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
        # print("\t\tABCI_N1_L::It got here")
        f.write('-3., 0.000\n')
        f.write('{} {}\n'.format(self.ri_L + self.b_L, WG_L))
        f.write('{} {}\n'.format(zr12_L[1][1], WG_L + self.L_L - zr12_L[1][0]))
        f.write('{} {}\n'.format(zr12_L[0][1], WG_L + self.L_L - zr12_L[0][0]))
        f.write('-3., 0.000\n')
        f.write('{} {}\n'.format(self.Req_L - self.B_L, WG_L + self.L_L))
        f.write('{} {}\n'.format(self.Req_L, WG_L + self.L_L))

    def abci_n1_R(self, n, zr12_R, WG_L, f):
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
        # print("\t\tABCI_N1_R::It got here")
        if n == 1:
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_L - self.B_R, WG_L + self.L_L))
            f.write('{} {}\n'.format(zr12_R[0][1], WG_L + self.L_L + zr12_R[0][0]))
            f.write('{} {}\n'.format(zr12_R[1][1], WG_L + self.L_L + zr12_R[1][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_R + self.b_R, WG_L + self.L_L + self.L_R))
            f.write('{} {}\n'.format(self.ri_R, WG_L + self.L_L + self.L_R))

        if n > 1:
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_R - self.B_R, WG_L + self.L_L + 2 * (n - 1) * self.L_M))
            f.write('{} {}\n'.format(zr12_R[0][1], WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[0][0]))
            f.write('{} {}\n'.format(zr12_R[1][1], WG_L + self.L_L + 2 * (n - 1) * self.L_M + zr12_R[1][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_R + self.b_R, WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R))
            f.write('{} {}\n'.format(self.ri_R, WG_L + self.L_L + 2 * (n - 1) * self.L_M + self.L_R))

    def abci_bp_R(self, n, zr12_BPR, WG_L, f):
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
        # print("\t\tABCI_BPR::It got here")
        if n == 1:
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_R + self.bt_R, WG_L + self.L_L + self.L_R))
            f.write('{} {}\n'.format(zr12_BPR[1][1], WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0]))
            f.write('{} {}\n'.format(zr12_BPR[0][1], WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Rbp_R - self.c_R, WG_L + self.L_L + self.L_R + self.x_R))
            f.write('{} {}\n'.format(self.Rbp_R, WG_L + self.L_L + self.L_R + self.x_R))

        if n > 1:
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_R + self.bt_R, WG_L + self.L_L + self.L_R + 2 * (n - 1) * self.L_M))
            f.write('{} {}\n'.format(zr12_BPR[1][1],
                                     WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[1][0] + 2 * (n - 1) * self.L_M))
            f.write('{} {}\n'.format(zr12_BPR[0][1],
                                     WG_L + self.L_L + self.L_R + self.x_R - zr12_BPR[0][0] + 2 * (n - 1) * self.L_M))
            f.write('-3., 0.000\n')
            f.write(
                '{} {}\n'.format(self.Rbp_R - self.c_R, WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M))
            f.write('{} {}\n'.format(self.Rbp_R, WG_L + self.L_L + self.L_R + self.x_R + 2 * (n - 1) * self.L_M))

    def abci_M(self, n, zr12_M, WG_L, f, i, end_type):
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
        # print("\t\tABCI_M::It got here")
        # Left and right Cell
        # First Half cell
        if i == 1 and end_type == 1:

            zr12_LM, alpha_LM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r

            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_L - self.B_M, WG_L + self.L_L + (2 * i - 2) * self.L_M))
            f.write('{} {}\n'.format(zr12_LM[0][1], WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_LM[0][0]))
            f.write('{} {}\n'.format(zr12_LM[1][1], WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_LM[1][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_M + self.b_M, WG_L + self.L_L + (2 * i - 1) * self.L_M))
            f.write('{} {}\n'.format(self.ri_M, WG_L + self.L_L + (2 * i - 1) * self.L_M))

        elif i == 1 and end_type != 1:

            zr12_L, alpha_L = self.rz_conjug('left')  # zr12_R first column is z , second column is r

            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_L - self.B_L, WG_L + self.L_L + (2 * i - 2) * self.L_M))
            f.write('{} {}\n'.format(zr12_L[0][1], WG_L + self.L_L + (2 * i - 2) * self.L_L + zr12_L[0][0]))
            f.write('{} {}\n'.format(zr12_L[1][1], WG_L + self.L_L + (2 * i - 2) * self.L_L + zr12_L[1][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_M + self.b_L, WG_L + self.L_L + (2 * i - 1) * self.L_L))
            f.write('{} {}\n'.format(self.ri_M, WG_L + self.L_L + (2 * i - 1) * self.L_L))

        if i != 1:
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_M - self.B_M, WG_L + self.L_L + (2 * i - 2) * self.L_M))
            f.write('{} {}\n'.format(zr12_M[0][1], WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_M[0][0]))
            f.write('{} {}\n'.format(zr12_M[1][1], WG_L + self.L_L + (2 * i - 2) * self.L_M + zr12_M[1][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_M + self.b_M, WG_L + self.L_L + (2 * i - 1) * self.L_M))
            f.write('{} {}\n'.format(self.ri_M, WG_L + self.L_L + (2 * i - 1) * self.L_M))

        # Second half cells
        if i == n - 1 and end_type == 1:
            zr12_RM, alpha_RM = self.rz_conjug('mid')  # zr12_R first column is z , second column is r

            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_M + self.b_M, WG_L + self.L_L + (2 * i - 1) * self.L_M))
            f.write(
                '{} {}\n'.format(zr12_RM[1][1], WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_RM[1][0]))
            f.write(
                '{} {}\n'.format(zr12_RM[0][1], WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_RM[0][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_R - self.B_M, WG_L + self.L_L + (2 * i) * self.L_M))
            f.write('{} {}\n'.format(self.Req_R, WG_L + self.L_L + (2 * i) * self.L_M))

        elif i == n - 1 and end_type != 1:

            zr12_R, alpha_R = self.rz_conjug('right')  # zr12_R first column is z , second column is r

            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_M + self.b_R, WG_L + self.L_L + (2 * i - 1) * self.L_M))
            f.write('{} {}\n'.format(zr12_R[1][1], WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R - zr12_R[1][0]))
            f.write('{} {}\n'.format(zr12_R[0][1], WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R - zr12_R[0][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_R - self.B_R, WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R))
            f.write('{} {}\n'.format(self.Req_R, WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_R))

        if i != n - 1:
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.ri_M + self.b_M, WG_L + self.L_L + (2 * i - 1) * self.L_M))
            f.write('{} {}\n'.format(zr12_M[1][1], WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_M[1][0]))
            f.write('{} {}\n'.format(zr12_M[0][1], WG_L + self.L_L + (2 * i - 1) * self.L_M + self.L_M - zr12_M[0][0]))
            f.write('-3., 0.000\n')
            f.write('{} {}\n'.format(self.Req_M - self.B_M, WG_L + self.L_L + (2 * i) * self.L_M))
            f.write('{} {}\n'.format(self.Req_M, WG_L + self.L_L + (2 * i) * self.L_M))

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

        # ic(zr12, alpha, df[-2])
        return zr12, alpha
