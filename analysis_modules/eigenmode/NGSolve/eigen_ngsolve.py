import itertools
import json
import os.path
import pickle
import shutil
import time
from sys import exit

from matplotlib import tri
from scipy.interpolate import CloughTocher2DInterpolator

from utils.shared_functions import *
from pathlib import Path
from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.signal import find_peaks

mu0 = 4 * pi * 1e-7
eps0 = 8.85418782e-12
c0 = 299792458


class NGSolveMEVP:
    def __init__(self):
        """
        MEVP for accelerating cavities with NGSolve

        ..math:



        """
        pass

    @staticmethod
    def ell_ips_grads(x_center, y_center, a, b, step, start, end, plot=True):
        h = x_center  # x-position of the center
        k = y_center  # y-position of the center
        a = a  # radius on the x-axis
        b = b  # radius on the y-axis

        t = np.arange(0, 2 * np.pi, np.pi / 100)

        x = h + a * np.cos(t)
        y = k + b * np.sin(t)
        pts = np.column_stack((x, y))
        inidx = np.all(np.logical_and(np.array(start) < pts, np.array(end) > pts), axis=1)
        inbox = pts[inidx]
        inbox = inbox[inbox[:, 0].argsort()]

        ips, gipts, Nps = [], {}, 5
        for n in range(1, Nps):
            # get inbox midpoint
            ip = inbox[int(n * len(inbox) / Nps)]  # <- spline interpolation point
            gipt = b ** 2 / (2 * (ip[1] - k)) - (b ** 2 / a ** 2) * (
                    (ip[0] - h) / (ip[1] - k))  # <- gradient at interpolation point

            ips.append(tuple(ip))
            gipts[n] = tuple([1, gipt])

            # x = np.array([ip[0] * 0.98, ip[0] * 1.02])
            # y = gipt * x + ip[1] - gipt * ip[0]
        #         plt.plot(x, y)
        #         plt.scatter(ip[0], ip[1], s=10, c='r')

        #     # get inbox midpoint
        #     midpoint = inbox[int(len(inbox)/2)]
        #     grad_midpoint = b**2/(2*(midpoint[1]-k)) - (b**2/a**2)*((midpoint[0]-h)/(midpoint[1]-k))

        #     # get midpoint neighbour 1
        #     n1 = inbox[int(len(inbox)/4)]
        #     grad_n1 = b**2/(2*(n1[1]-k)) - (b**2/a**2)*((n1[0]-h)/(n1[1]-k))

        #     # get midpoint
        #     n2 = inbox[int(3*len(inbox)/4)]
        #     grad_n2 = b**2/(2*(n2[1]-k)) - (b**2/a**2)*((n2[0]-h)/(n2[1]-k))

        #     x = np.array([midpoint[0]*0.98, midpoint[0]*1.02])
        #     y = grad_midpoint*x + midpoint[1] - grad_midpoint*midpoint[0]

        #     x1 = np.array([n1[0]*0.98, n1[0]*1.02])
        #     y1 = grad_n1*x1 + n1[1] - grad_n1*n1[0]

        #     x2 = np.array([n2[0]*0.98, n2[0]*1.02])
        #     y2 = grad_n2*x2 + n2[1] - grad_n2*n2[0]

        #     plt.plot(inbox[:, 0], inbox[:, 1])
        #     plt.plot(x, y)
        #     plt.plot(x1, y1)
        #     plt.plot(x2, y2)

        #     # plot box
        #     xs, ys = start
        #     xe, ye = end
        #     plt.plot([xs, xe], [ys, ys], 'b')
        #     plt.plot([xs, xs], [ys, ye], 'b')
        #     plt.plot([xs, xe], [ye, ye], 'b')
        #     plt.plot([xe, xe], [ys, ye], 'b')

        nlast = len(ips) + 1
        return ips, gipts, nlast

    def write_geometry(self, folder, n_cells, mid_cell, end_cell_left=None, end_cell_right=None,
                       beampipe='none', plot=False, cell_type='normal'):
        """
        Define geometry

        Parameters
        ----------
        folder
        n_cells: int
            Number of cavity cells
        mid_cell: array like
            Mid cell cell geometric parameters
        end_cell_left: array like
            Left end cell cell geometric parameters
        end_cell_right: array like
            Right end cell cell geometric parameters
        beampipe: {"none", "both", "left", "right"}
            Which sides to include beampipes
        plot: bool
            Show geometry after definition or not

        Returns
        -------

        """

        n_cells = n_cells
        mid_cell = mid_cell

        if end_cell_left is None:
            end_cell_left = mid_cell
        else:
            end_cell_left = end_cell_left

        if end_cell_right is None:
            if self.n_cells > 1:
                if end_cell_left is None:
                    self.end_cell_right = mid_cell
                else:
                    self.end_cell_right = end_cell_left
            else:
                self.end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_right

        # check if folder exists
        if not os.path.exists(folder):
            try:
                os.mkdir(folder)
            except FileNotFoundError:
                ic("There was a problem creating the directory for the simulation files. Please check folder.")
                exit()

        file_path = fr"{folder}\geodata.n"
        if cell_type == 'normal':
            writeCavityForMultipac(file_path, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe, plot=plot)
        else:
            writeCavityForMultipac_flat_top(file_path, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe,
                                            plot=plot)

    def write_geometry_multicell(self, folder, n_cells, mid_cell, end_cell_left=None, end_cell_right=None,
                                 beampipe='none', plot=False, cell_type='normal'):
        """
        Define geometry

        Parameters
        ----------
        folder
        n_cells: int
            Number of cavity cells
        mid_cell: array like
            Mid cell cell geometric parameters
        end_cell_left: array like
            Left end cell cell geometric parameters
        end_cell_right: array like
            Right end cell cell geometric parameters
        beampipe: {"none", "both", "left", "right"}
            Which sides to include beampipes
        plot: bool
            Show geometry after definition or not

        Returns
        -------

        """

        n_cells = n_cells
        mid_cell = mid_cell

        if end_cell_left is None:
            end_cell_left = mid_cell
        else:
            end_cell_left = end_cell_left

        if end_cell_right is None:
            if self.n_cells > 1:
                if end_cell_left is None:
                    self.end_cell_right = mid_cell
                else:
                    self.end_cell_right = end_cell_left
            else:
                self.end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_right

        beampipe = beampipe

        # check if folder exists
        if not os.path.exists(folder):
            try:
                os.mkdir(folder)
            except FileNotFoundError:
                ic("There was a problem creating the directory for the simulation files. Please check folder.")
                exit()

        file_path = fr"{folder}\geodata.n"
        if cell_type == 'normal':
            writeCavityForMultipac_multicell(file_path, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe,
                                             plot=plot)
        else:
            writeCavityForMultipac_flat_top(file_path, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe,
                                            plot=plot)

    def cavgeom_ngsolve(self, n_cell, mid_cell, end_cell_left, end_cell_right, beampipe='both', draw=False):
        """
        Write cavity geometry to be used by Multipac for multipacting analysis
        Parameters
        ----------
        n_cell: int
            Number of cavity cells
        mid_cell: list, ndarray
            Array of cavity middle cells' geometric parameters
        end_cell_left: list, ndarray
            Array of cavity left end cell's geometric parameters
        end_cell_right: list, ndarray
            Array of cavity left end cell's geometric parameters
        beampipe: str {"left", "right", "both", "none"}
            Specify if beam pipe is on one or both ends or at no end at all
        draw: bool
            If True, the cavity geometry is plotted for viewing

        Returns
        -------

        """

        if end_cell_left is None:
            end_cell_left = mid_cell

        if end_cell_right is None:
            if end_cell_left is None:
                end_cell_right = mid_cell
            else:
                end_cell_right = end_cell_left

        A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(mid_cell[:7]) * 1e-3
        A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(end_cell_left[:7]) * 1e-3
        A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(end_cell_right[:7]) * 1e-3

        step = 0.1  # step in boundary points in mm

        if beampipe.lower() == 'both':
            L_bp_l = 4 * L_m
            L_bp_r = 4 * L_m
        elif beampipe.lower() == 'none':
            L_bp_l = 0.000  # 4 * L_m  #
            L_bp_r = 0.000  # 4 * L_m  #
        elif beampipe.lower() == 'left':
            L_bp_l = 4 * L_m
            L_bp_r = 0.000
        elif beampipe.lower() == 'right':
            L_bp_l = 0.000
            L_bp_r = 4 * L_m
        else:
            L_bp_l = 0.000  # 4 * L_m  #
            L_bp_r = 0.000  # 4 * L_m  #

        # calculate shift
        shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

        # calculate angles outside loop
        # CALCULATE x1_el, y1_el, x2_el, y2_el

        df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
        x1el, y1el, x2el, y2el = df[0]

        # CALCULATE x1, y1, x2, y2
        df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
        x1, y1, x2, y2 = df[0]

        # CALCULATE x1_er, y1_er, x2_er, y2_er
        df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
        x1er, y1er, x2er, y2er = df[0]

        # start workplane
        wp = WorkPlane()

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        wp.MoveTo(-shift, 0)
        wp.LineTo(*[-shift, Ri_el])  # <-
        wp.LineTo(*[L_bp_l - shift, Ri_el])  # <- add left beampipe
        pt = [L_bp_l - shift, Ri_el]

        for n in range(1, n_cell + 1):
            if n == 1:
                ips, gipts, nlast = self.ell_ips_grads(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt,
                                                       [-shift + x1el, y1el])
                wp.Spline([*ips, (-shift + x1el, y1el)],
                          tangents={0: (1, 0), nlast: (x2el - x1el, y2el - y1el), **gipts})  # <- draw iris arc

                wp.LineTo(*[-shift + x2el, y2el])  # <- draw tangent line

                pt = [-shift + x2el, y2el]
                ips, gipts, nlast = self.ell_ips_grads(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt,
                                                       [L_bp_l + L_el - shift, Req])
                wp.Spline([*ips, (L_bp_l + L_el - shift, Req)],
                          tangents={0: (x2el - x1el, y2el - y1el), nlast: (1, 0), **gipts})  # <- draw left eq. arc
                pt = [L_bp_l + L_el - shift, Req]

                if n_cell == 1:
                    ips, gipts, nlast = self.ell_ips_grads(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step,
                                                           [pt[0], Req - B_er],
                                                           [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req])
                    wp.Spline([*ips, (L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er)],
                              tangents={0: (1, 0), nlast: (x2er - x1er, y1er - y2er), **gipts})  # <- draw left eq. arc

                    wp.LineTo(*[L_el + L_er - x1er + + L_bp_l + L_bp_r - shift, y1er])  # <- draw tangent line

                    pt = [L_el + L_er - x1er + + L_bp_l + L_bp_r - shift, y1er]
                    ips, gipts, nlast = self.ell_ips_grads(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step,
                                                           [pt[0], Ri_er],
                                                           [L_bp_l + L_el + L_er - shift, y1er])
                    wp.Spline([*ips, (L_bp_l + L_el + L_er - shift, Ri_er)],
                              tangents={0: (x2er - x1er, y1er - y2er), nlast: (1, 0), **gipts})  # <- draw left eq. arc
                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    ips, gipts, nlast = self.ell_ips_grads(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step,
                                                           [pt[0], Req - B_m],
                                                           [L_el + L_m - x2 + 2 * L_bp_l - shift, Req])
                    wp.Spline([*ips, (L_el + L_m - x2 + 2 * L_bp_l - shift, y2)],
                              tangents={0: (1, 0), nlast: (x2 - x1, y1 - y2), **gipts})  # <- draw left eq. arc

                    wp.LineTo(*[L_el + L_m - x1 + 2 * L_bp_l - shift, y1])  # <- draw tangent line

                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    ips, gipts, nlast = self.ell_ips_grads(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step,
                                                           [pt[0], Ri_m],
                                                           [L_bp_l + L_el + L_m - shift, y1])
                    wp.Spline([*ips, (L_bp_l + L_el + L_m - shift, Ri_m)],
                              tangents={0: (x2 - x1, y1 - y2), nlast: (1, 0), **gipts})  # <- draw left eq. arc
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]

                    # calculate new shift
                    shift = shift - (L_el + L_m)

            elif n > 1 and n != n_cell:
                ips, gipts, nlast = self.ell_ips_grads(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt,
                                                       [-shift + x1, y1])
                wp.Spline([*ips, (-shift + x1, y1)],
                          tangents={0: (1, 0), nlast: (x2 - x1, y2 - y1), **gipts})  # <- draw iris arc

                wp.LineTo(*[-shift + x2, y2])  # <- draw tangent line

                pt = [-shift + x2, y2]
                ips, gipts, nlast = self.ell_ips_grads(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt,
                                                       [L_bp_l + L_m - shift, Req])
                wp.Spline([*ips, (L_bp_l + L_m - shift, Req)],
                          tangents={0: (x2 - x1, y2 - y1), nlast: (1, 0), **gipts})  # <- draw left eq. arc

                pt = [L_bp_l + L_m - shift, Req]
                ips, gipts, nlast = self.ell_ips_grads(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step,
                                                       [pt[0], Req - B_m],
                                                       [L_m + L_m - x2 + 2 * L_bp_l - shift, Req])
                wp.Spline([*ips, (L_m + L_m - x2 + 2 * L_bp_l - shift, y2)],
                          tangents={0: (1, 0), nlast: (x2 - x1, y1 - y2), **gipts})  # <- draw left eq. arc

                wp.LineTo(*[L_m + L_m - x1 + 2 * L_bp_l - shift, y1])  # <- draw tangent line

                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                ips, gipts, nlast = self.ell_ips_grads(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step,
                                                       [pt[0], Ri_m],
                                                       [L_bp_l + L_m + L_m - shift, y1])
                wp.Spline([*ips, (L_bp_l + L_m + L_m - shift, Ri_m)],
                          tangents={0: (x2 - x1, y1 - y2), nlast: (1, 0), **gipts})  # <- draw left eq. arc
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]

                # calculate new shift
                shift = shift - 2 * L_m
            else:
                ips, gipts, nlast = self.ell_ips_grads(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt,
                                                       [-shift + x1, y1])
                wp.Spline([*ips, (-shift + x1, y1)],
                          tangents={0: (1, 0), nlast: (x2 - x1, y2 - y1), **gipts})  # <- draw iris arc

                wp.LineTo(*[-shift + x2, y2])  # <- draw tangent line

                pt = [-shift + x2, y2]
                ips, gipts, nlast = self.ell_ips_grads(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt,
                                                       [L_bp_l + L_m - shift, Req])
                wp.Spline([*ips, (L_bp_l + L_m - shift, Req)],
                          tangents={0: (x2 - x1, y2 - y1), nlast: (1, 0), **gipts})  # <- draw left eq. arc

                pt = [L_bp_l + L_m - shift, Req]
                ips, gipts, nlast = self.ell_ips_grads(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step,
                                                       [pt[0], Req - B_er],
                                                       [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req])
                wp.Spline([*ips, (L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er)],
                          tangents={0: (1, 0), nlast: (x2er - x1er, y1er - y2er), **gipts})  # <- draw left eq. arc

                wp.LineTo(*[L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er])  # <- draw tangent line

                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                ips, gipts, nlast = self.ell_ips_grads(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step,
                                                       [pt[0], Ri_er],
                                                       [L_bp_l + L_m + L_er - shift, y1er])
                wp.Spline([*ips, (L_bp_l + L_m + L_er - shift, Ri_er)],
                          tangents={0: (x2er - x1er, y1er - y2er), nlast: (1, 0), **gipts})  # <- draw left eq. arc
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]

        # BEAM PIPE
        # reset shift
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
        wp.LineTo(*[2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er])  # <- right beampipe
        wp.LineTo(*[2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0])

        wp.Close().Reverse()
        face = wp.Face()
        # name the boundaries
        face.edges.Max(X).name = "r"
        face.edges.Max(X).col = (1, 0, 0)
        face.edges.Min(X).name = "l"
        face.edges.Min(X).col = (1, 0, 0)
        face.edges.Min(Y).name = "b"
        face.edges.Min(Y).col = (1, 0, 0)
        # Draw(face)

        return face

    def cavity(self, no_of_cells=1, no_of_modules=1, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
               fid=None, bc=33, pol='Monopole', f_shift='default', beta=1, n_modes=None, beampipes='None',
               sim_folder='NGSolveMEVP', parentDir=None, projectDir=None, subdir='',
               expansion=None, expansion_r=None, mesh_args=None, opt=False, deformation_params=None):
        """
        Write geometry file and run eigenmode analysis with NGSolveMEVP

        Parameters
        ----------
        pol
        expansion
        no_of_cells: int
            Number of cells
        no_of_modules: int
            Number of modules
        mid_cells_par: list, array like
            Mid cell geometric parameters -> [A, B, a, b, Ri, L, Req, alpha]
        l_end_cell_par: list, array like
            Left end cell geometric parameters -> [A_el, B_el, a_el, b_el, Ri_el, L_el, Req, alpha_el]
        r_end_cell_par: list, array like
            Right end cell geometric parameters -> [A_er, B_er, a_er, b_er, Ri_er, L_er, Req, alpha_er]
        fid: int, str
            File id
        bc: int
            Boundary condition -> 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        f_shift: float
            Eigenvalue frequency shift
        beta: int, float
            Velocity ratio :math: `\\beta = \frac{v}{c}`
        n_modes: int
            Number of modes
        beampipes: {"left", "right", "both", "none"}
            Specify if beam pipe is on one or both ends or at no end at all
        parentDir: str
            Parent directory
        projectDir: str
            Project directory
        subdir: str
            Sub directory to save simulation results to
        mesh: list [Jxy, Jxy_bp, Jxy_bp_y]
            Mesh definition for logical mesh:
            Jxy -> Number of elements of logical mesh along JX and JY
            Jxy_bp -> Number of elements of logical mesh along JX in beampipe
            Jxy_bp_y -> Number of elements of logical mesh along JY in beampipe

        Returns
        -------

        """

        if pol != 'Monopole':
            pol_subdir = 'dipole'
        else:
            pol_subdir = 'monopole'

        # change save directory
        if mesh_args is None:
            mesh_args = [20, 20]
        if opt:
            # consider making better. This was just an adhoc fix
            run_save_directory = os.path.join(projectDir, fr'SimulationData/Optimisation/{fid}')
        else:
            # change save directory
            if subdir == '':
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{fid}/{pol_subdir}'
            else:
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{subdir}/{fid}/{pol_subdir}'

        # write
        self.write_geometry(run_save_directory, no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes,
                            plot=False)

        if os.path.exists(f'{run_save_directory}\geodata.n'):
            # read geometry
            cav_geom = pd.read_csv(f'{run_save_directory}\geodata.n',
                                   header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0, 2]]

            if deformation_params is not None:
                cav_geom = self.gaussian_deform(no_of_cells, cav_geom.drop_duplicates(subset=[0, 1]).to_numpy(),
                                                deformation_params)
                # save deformed cavity profile
                cav_geom.to_csv(f'{run_save_directory}\geodata_deformed.n', sep='\t')

            cav_geom = cav_geom[[1, 0]]
            # plt.plot(cav_geom[1], cav_geom[0], ls='--', lw=4)
            # plt.show()

            pnts = list(cav_geom.itertuples(index=False, name=None))
            wp = WorkPlane()
            wp.MoveTo(*pnts[0])
            for p in pnts[1:]:
                wp.LineTo(*p)
            wp.Close().Reverse()
            face = wp.Face()

            # # get face from ngsolve geom <- coming later
            # face = self.cavgeom_ngsolve(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)

            # name the boundaries
            face.edges.Max(X).name = "r"
            face.edges.Max(X).col = (1, 0, 0)
            face.edges.Min(X).name = "l"
            face.edges.Min(X).col = (1, 0, 0)
            face.edges.Min(Y).name = "b"
            face.edges.Min(Y).col = (1, 0, 0)

            geo = OCCGeometry(face, dim=2)

            # mesh
            A_m, B_m, a_m, b_m, Ri_m, L, Req = np.array(mid_cells_par[:7])
            maxh = L / mesh_args[0] * 1e-3

            # try to generate mesh
            try:
                ngmesh = geo.GenerateMesh(maxh=maxh)
            except Exception as e:
                error('Unable to generate mesh:: ', e)
                plt.plot(cav_geom[1], cav_geom[0])
                plt.show()
                exit()
            mesh = Mesh(ngmesh)
            mesh.Curve(3)

            # save mesh
            self.save_mesh(run_save_directory, mesh)

            # define finite element space
            fes = HCurl(mesh, order=1, dirichlet='default')

            u, v = fes.TnT()

            a = BilinearForm(y * curl(u) * curl(v) * dx)
            m = BilinearForm(y * u * v * dx)

            apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
            pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

            with TaskManager():
                a.Assemble()
                m.Assemble()
                apre.Assemble()
                # freedof_matrix = a.mat.CreateSmoother(fes.FreeDofs())

                # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
                gradmat, fesh1 = fes.CreateGradient()

                gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
                math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
                math1[0, 0] += 1  # fix the 1-dim kernel
                invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

                # build the Poisson projector with operator Algebra:
                proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

                projpre = proj @ pre.mat
                evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=no_of_cells + 1, maxit=mesh_args[1],
                                              printrates=False)

            # print out eigenvalues
            # print("Eigenvalues")
            freq_fes = []
            evals[0] = 1  # <- replace nan with zero
            for i, lam in enumerate(evals):
                freq_fes.append(c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6)
                # print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

            # plot results
            gfu_E = []
            gfu_H = []
            for i in range(len(evecs)):
                w = 2 * pi * freq_fes[i] * 1e6
                gfu = GridFunction(fes)
                gfu.vec.data = evecs[i]

                gfu_E.append(gfu)
                gfu_H.append(1j / (mu0 * w) * curl(gfu))

            # save fields
            self.save_fields(run_save_directory, gfu_E, gfu_H)

            # # alternative eigenvalue solver, but careful, mode numbering may change
            # u = GridFunction(fes, multidim=30, name='resonances')
            # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
            #                         list(u.vecs), shift=300)
            # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

            # save json file
            shape = {'IC': update_alpha(mid_cells_par),
                     'OC': update_alpha(l_end_cell_par),
                     'OC_R': update_alpha(r_end_cell_par)}
            # print(run_save_directory)
            with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
                json.dump(shape, f, indent=4, separators=(',', ': '))

            # qois = self.evaluate_qois(face, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
            qois = self.evaluate_qois(cav_geom, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
            # ic(qois)

            with open(fr'{run_save_directory}\qois.json', "w") as f:
                json.dump(qois, f, indent=4, separators=(',', ': '))
            return True
        else:
            error('Could not run eigenmode analysis due to error in geometry.')
            return False

    def cavity_multicell(self, no_of_cells=1, no_of_modules=1, mid_cells_par=None, l_end_cell_par=None,
                         r_end_cell_par=None,
                         fid=None, bc=33, pol='monopole', f_shift='default', beta=1, n_modes=None, beampipes='None',
                         sim_folder='NGSolveMEVP', parentDir=None, projectDir=None, subdir='',
                         expansion=None, expansion_r=None, mesh_args=None, opt=False,
                         deformation_params=None):
        """
        Write geometry file and run eigenmode analysis with NGSolveMEVP

        Parameters
        ----------
        pol
        expansion
        no_of_cells: int
            Number of cells
        no_of_modules: int
            Number of modules
        mid_cells_par: list, ndarray
            Mid cell geometric parameters -> [A, B, a, b, Ri, L, Req, alpha]
        l_end_cell_par: list, ndarray
            Left end cell geometric parameters -> [A_el, B_el, a_el, b_el, Ri_el, L_el, Req, alpha_el]
        r_end_cell_par: list, ndarray
            Right end cell geometric parameters -> [A_er, B_er, a_er, b_er, Ri_er, L_er, Req, alpha_er]
        fid: int, str
            File id
        bc: int
            Boundary condition -> 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        f_shift: float
            Eigenvalue frequency shift
        beta: int, float
            Velocity ratio :math: `\\beta = \frac{v}{c}`
        n_modes: int
            Number of modes
        beampipes: {"left", "right", "both", "none"}
            Specify if beam pipe is on one or both ends or at no end at all
        parentDir: str
            Parent directory
        projectDir: str
            Project directory
        subdir: str
            Sub directory to save simulation results to
        mesh: list [Jxy, Jxy_bp, Jxy_bp_y]
            Mesh definition for logical mesh:
            Jxy -> Number of elements of logical mesh along JX and JY
            Jxy_bp -> Number of elements of logical mesh along JX in beampipe
            Jxy_bp_y -> Number of elements of logical mesh along JY in beampipe

        Returns
        -------

        """

        start = time.time()
        if pol.lower() != 'monopole':
            pol_subdir = 'dipole'
        else:
            pol_subdir = 'monopole'

        # change save directory
        if mesh_args is None:
            mesh_args = [20, 20]
        if opt:  # consider making better. This was just an adhoc fix
            run_save_directory = projectDir / fr'SimulationData/Optimisation/{fid}'
        else:
            # change save directory
            if subdir == '':
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{fid}/{pol_subdir}'
            else:
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{subdir}/{fid}/{pol_subdir}'

        # write
        # st1 = time.time()
        self.write_geometry_multicell(run_save_directory, no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par,
                                      beampipes, plot=False)
        # ic('Time to write geom: ', time.time()-st1)

        # read geometry
        # ic(run_save_directory)
        cav_geom = pd.read_csv(f'{run_save_directory}\geodata.n',
                               header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0, 2]]

        if deformation_params is not None:
            cav_geom = self.gaussian_deform(no_of_cells, cav_geom.drop_duplicates(subset=[0, 1]).to_numpy(),
                                            deformation_params)
            # save deformed cavity profile
            cav_geom.to_csv(f'{run_save_directory}\geodata_deformed.n', sep='\t')

        cav_geom = cav_geom[[1, 0]]
        # plt.plot(cav_geom[0], cav_geom[1])
        # plt.show()

        pnts = list(cav_geom.itertuples(index=False, name=None))
        wp = WorkPlane()
        wp.MoveTo(*pnts[0])
        for p in pnts[1:]:
            try:
                wp.LineTo(*p)
            except:
                pass
        wp.Close().Reverse()
        face = wp.Face()

        # # get face from ngsolve geom <- coming later
        # face = self.cavgeom_ngsolve(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)

        # name the boundaries
        face.edges.Max(X).name = "r"
        face.edges.Max(X).col = (1, 0, 0)
        face.edges.Min(X).name = "l"
        face.edges.Min(X).col = (1, 0, 0)
        face.edges.Min(Y).name = "b"
        face.edges.Min(Y).col = (1, 0, 0)

        geo = OCCGeometry(face, dim=2)

        # mesh
        A_m, B_m, a_m, b_m, Ri_m, L, Req = np.array(l_end_cell_par[:7])
        maxh = L / mesh_args[0] * 1e-3
        # ic(maxh)
        # st1 = time.time()
        ngmesh = geo.GenerateMesh(maxh=maxh)
        mesh = Mesh(ngmesh)
        mesh.Curve(3)

        # save mesh
        self.save_mesh(run_save_directory, mesh)
        # ic('Time to write geom: ', time.time()-st1)

        # define finite element space
        fes = HCurl(mesh, order=1, dirichlet='default')

        u, v = fes.TnT()

        a = BilinearForm(y * curl(u) * curl(v) * dx)  #.Assemble()
        m = BilinearForm(y * u * v * dx)  #.Assemble()

        apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
        pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

        with TaskManager():
            a.Assemble()
            m.Assemble()
            apre.Assemble()
            # freedof_matrix = a.mat.CreateSmoother(fes.FreeDofs())

            # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
            gradmat, fesh1 = fes.CreateGradient()

            gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
            math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
            math1[0, 0] += 1  # fix the 1-dim kernel
            invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

            # build the Poisson projector with operator Algebra:
            proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

            projpre = proj @ pre.mat
            evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=no_of_cells + 1, maxit=mesh_args[1],
                                          printrates=False)

        # ic('Time to complete sim: ', time.time()-start)
        # print out eigenvalues
        # print("Eigenvalues")
        # freq_fes = []
        evals[0] = 1  # <- replace nan with zero
        # for i, lam in enumerate(evals):
        freq_fes = c0 * np.sqrt(evals) / (2 * np.pi) * 1e-6
        # print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

        # plot results
        gfu_E = []
        gfu_H = []
        for i in range(len(evecs)):
            w = 2 * pi * freq_fes[i] * 1e6
            gfu = GridFunction(fes)
            gfu.vec.data = evecs[i]

            gfu_E.append(gfu)
            gfu_H.append(1j / (mu0 * w) * curl(gfu))

        # save fields
        self.save_fields(run_save_directory, gfu_E, gfu_H)

        # # alternative eigenvalue solver, but careful, mode numbering may change
        # u = GridFunction(fes, multidim=30, name='resonances')
        # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
        #                         list(u.vecs), shift=300)
        # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

        # save json file
        shape = {'IC': mid_cells_par.tolist(),
                 'OC': l_end_cell_par.tolist(),
                 'OC_R': r_end_cell_par.tolist()}
        # print(run_save_directory)
        with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))

        # qois = self.evaluate_qois(face, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
        qois = self.evaluate_qois(cav_geom, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
        # ic(qois)

        with open(fr'{run_save_directory}\qois.json', "w") as f:
            json.dump(qois, f, indent=4, separators=(',', ': '))

    def cavity_flattop(self, no_of_cells=1, no_of_modules=1,
                       mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
                       fid=None, bc=33, pol='Monopole', f_shift='default', beta=1, n_modes=None, beampipes='None',
                       sim_folder='NGSolveMEVP', parentDir=None, projectDir=None, subdir='',
                       expansion=None, expansion_r=None, mesh_args=None, opt=False,
                       deformation_params=None):
        """
        Write geometry file and run eigenmode analysis with NGSolveMEVP

        Parameters
        ----------
        pol
        expansion
        no_of_cells: int
            Number of cells
        no_of_modules: int
            Number of modules
        mid_cells_par: list, array like
            Mid cell geometric parameters -> [A, B, a, b, Ri, L, Req, alpha]
        l_end_cell_par: list, array like
            Left end cell geometric parameters -> [A_el, B_el, a_el, b_el, Ri_el, L_el, Req, alpha_el]
        r_end_cell_par: list, array like
            Right end cell geometric parameters -> [A_er, B_er, a_er, b_er, Ri_er, L_er, Req, alpha_er]
        fid: int, str
            File id
        bc: int
            Boundary condition -> 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        f_shift: float
            Eigenvalue frequency shift
        beta: int, float
            Velocity ratio :math: `\\beta = \frac{v}{c}`
        n_modes: int
            Number of modes
        beampipes: {"left", "right", "both", "none"}
            Specify if beam pipe is on one or both ends or at no end at all
        parentDir: str
            Parent directory
        projectDir: str
            Project directory
        subdir: str
            Sub directory to save simulation results to
        mesh: list [Jxy, Jxy_bp, Jxy_bp_y]
            Mesh definition for logical mesh:
            Jxy -> Number of elements of logical mesh along JX and JY
            Jxy_bp -> Number of elements of logical mesh along JX in beampipe
            Jxy_bp_y -> Number of elements of logical mesh along JY in beampipe

        Returns
        -------

        """

        if pol != 'Monopole':
            pol_subdir = 'dipole'
        else:
            pol_subdir = 'monopole'

        # change save directory
        if mesh_args is None:
            mesh_args = [20, 20]
        if opt:  # consider making better. This was just an adhoc fix
            run_save_directory = projectDir / fr'SimulationData/Optimisation/{fid}'
        else:
            # change save directory
            if subdir == '':
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{fid}/{pol_subdir}'
            else:
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{subdir}/{fid}/{pol_subdir}'
        # write
        self.write_geometry(run_save_directory, no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes,
                            cell_type='flattop', plot=False)
        # print('done writing geometry')

        # read geometry
        cav_geom = pd.read_csv(f'{run_save_directory}\geodata.n',
                               header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0, 2]]

        if deformation_params is not None:
            cav_geom = self.gaussian_deform(no_of_cells, cav_geom.drop_duplicates(subset=[0, 1]).to_numpy(),
                                            deformation_params)
            # save deformed cavity profile
            cav_geom.to_csv(f'{run_save_directory}\geodata_deformed.n', sep='\t')

        cav_geom = cav_geom[[1, 0]]
        # plt.plot(cav_geom[0], cav_geom[1])
        # plt.show()

        pnts = list(cav_geom.itertuples(index=False, name=None))
        wp = WorkPlane()
        wp.MoveTo(*pnts[0])
        for p in pnts[1:]:
            wp.LineTo(*p)
        wp.Close().Reverse()
        face = wp.Face()

        # # get face from ngsolve geom <- coming later
        # face = self.cavgeom_ngsolve(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)

        # name the boundaries
        face.edges.Max(X).name = "r"
        face.edges.Max(X).col = (1, 0, 0)
        face.edges.Min(X).name = "l"
        face.edges.Min(X).col = (1, 0, 0)
        face.edges.Min(Y).name = "b"
        face.edges.Min(Y).col = (1, 0, 0)

        geo = OCCGeometry(face, dim=2)

        # mesh
        A_m, B_m, a_m, b_m, Ri_m, L, Req = np.array(mid_cells_par[:7])
        maxh = L / mesh_args[0] * 1e-3
        ngmesh = geo.GenerateMesh(maxh=maxh)
        mesh = Mesh(ngmesh)
        mesh.Curve(3)

        # save mesh
        self.save_mesh(run_save_directory, mesh)

        # define finite element space
        fes = HCurl(mesh, order=1, dirichlet='default')

        u, v = fes.TnT()

        a = BilinearForm(y * curl(u) * curl(v) * dx).Assemble()
        m = BilinearForm(y * u * v * dx).Assemble()

        apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
        pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

        with TaskManager():
            a.Assemble()
            m.Assemble()
            apre.Assemble()
            freedof_matrix = a.mat.CreateSmoother(fes.FreeDofs())

            # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
            gradmat, fesh1 = fes.CreateGradient()

            gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
            math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
            math1[0, 0] += 1  # fix the 1-dim kernel
            invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

            # build the Poisson projector with operator Algebra:
            proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

            projpre = proj @ pre.mat
            evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=no_of_cells + 1, maxit=mesh_args[1],
                                          printrates=False)

        # print out eigenvalues
        freq_fes = []
        evals[0] = 1  # <- replace nan with zero
        for i, lam in enumerate(evals):
            freq_fes.append(c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6)
            # print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

        # plot results
        gfu_E = []
        gfu_H = []
        for i in range(len(evecs)):
            w = 2 * pi * freq_fes[i] * 1e6
            gfu = GridFunction(fes)
            gfu.vec.data = evecs[i]

            gfu_E.append(gfu)
            gfu_H.append(1j / (mu0 * w) * curl(gfu))

        # save fields
        self.save_fields(run_save_directory, gfu_E, gfu_H)

        # # alternative eigenvalue solver, but careful, mode numbering may change
        # u = GridFunction(fes, multidim=30, name='resonances')
        # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
        #                         list(u.vecs), shift=300)
        # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

        # save json file
        shape = {'IC': update_alpha(mid_cells_par, cell_type='flattop'),
                 'OC': update_alpha(l_end_cell_par, cell_type='flattop'),
                 'OC_R': update_alpha(r_end_cell_par, cell_type='flattop')}

        with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))

        # qois = self.evaluate_qois(face, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
        qois = self.evaluate_qois(cav_geom, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
        # ic(qois)

        with open(fr'{run_save_directory}\qois.json', "w") as f:
            json.dump(qois, f, indent=4, separators=(',', ': '))

    def geometry_from_file(self, filepath):
        pass

    def geometry_from_array(self, ):
        pass

    def pillbox(self, no_of_cells=1, no_of_modules=1, cell_par=None,
                fid=None, bc=33, pol='Monopole', f_shift='default', beta=1, n_modes=None, beampipes='None',
                sim_folder='NGSolveMEVP', parentDir=None, projectDir=None, subdir='',
                expansion=None, expansion_r=None, mesh_args=None, opt=False, deformation_params=None):

        if pol != 'Monopole':
            pol_subdir = 'dipole'
        else:
            pol_subdir = 'monopole'

        # change save directory
        if mesh_args is None:
            mesh_args = [20, 20]
        if opt:
            # consider making better. This was just an adhoc fix
            run_save_directory = os.path.join(projectDir, fr'SimulationData/Optimisation/{fid}')
        else:
            # change save directory
            if subdir == '':
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{fid}/{pol_subdir}'
            else:
                run_save_directory = projectDir / fr'SimulationData/{sim_folder}/{subdir}/{fid}/{pol_subdir}'

        # write geometry
        file_path = f'{run_save_directory}\geodata.n'
        write_pillbox_geometry(file_path, no_of_cells, cell_par, beampipes)

        if os.path.exists(file_path):
            # read geometry
            cav_geom = pd.read_csv(f'{run_save_directory}\geodata.n',
                                   header=None, sep='\s+', engine='python')[[1, 0, 2]]

            if deformation_params is not None:
                cav_geom = self.gaussian_deform(no_of_cells, cav_geom.drop_duplicates(subset=[0, 1]).to_numpy(),
                                                deformation_params)
                # save deformed cavity profile
                cav_geom.to_csv(f'{run_save_directory}\geodata_deformed.n', sep='\t')

            cav_geom = cav_geom[[1, 0]]
            # plt.plot(cav_geom[1], cav_geom[0], ls='--', lw=4)
            # plt.show()

            pnts = list(cav_geom.itertuples(index=False, name=None))
            wp = WorkPlane()
            wp.MoveTo(*pnts[0])
            for p in pnts[1:]:
                wp.LineTo(*p)
            wp.Close().Reverse()
            face = wp.Face()

            # name the boundaries
            face.edges.Max(X).name = "r"
            face.edges.Max(X).col = (1, 0, 0)
            face.edges.Min(X).name = "l"
            face.edges.Min(X).col = (1, 0, 0)
            face.edges.Min(Y).name = "b"
            face.edges.Min(Y).col = (1, 0, 0)

            geo = OCCGeometry(face, dim=2)

            # mesh
            L, Req = cell_par[0], cell_par[1]
            maxh = L / mesh_args[0] * 1e-3
            print(maxh)

            # try to generate mesh
            try:
                ngmesh = geo.GenerateMesh(maxh=maxh)
            except Exception as e:
                error('Unable to generate mesh:: ', e)
                plt.plot(cav_geom[1], cav_geom[0])
                plt.show()
                exit()
            mesh = Mesh(ngmesh)
            mesh.Curve(3)

            # save mesh
            self.save_mesh(run_save_directory, mesh)

            # define finite element space
            fes = HCurl(mesh, order=1, dirichlet='default')

            u, v = fes.TnT()

            a = BilinearForm(y * curl(u) * curl(v) * dx)
            m = BilinearForm(y * u * v * dx)

            apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
            pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

            with TaskManager():
                a.Assemble()
                m.Assemble()
                apre.Assemble()

                # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
                gradmat, fesh1 = fes.CreateGradient()

                gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
                math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
                math1[0, 0] += 1  # fix the 1-dim kernel
                invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

                # build the Poisson projector with operator Algebra:
                proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

                projpre = proj @ pre.mat
                evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=no_of_cells + 1, maxit=mesh_args[1],
                                              printrates=False)

            freq_fes = []
            evals[0] = 1  # <- replace nan with zero
            for i, lam in enumerate(evals):
                freq_fes.append(c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6)
                # print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

            # plot results
            gfu_E = []
            gfu_H = []
            for i in range(len(evecs)):
                w = 2 * pi * freq_fes[i] * 1e6
                gfu = GridFunction(fes)
                gfu.vec.data = evecs[i]

                gfu_E.append(gfu)
                gfu_H.append(1j / (mu0 * w) * curl(gfu))

            # save fields
            self.save_fields(run_save_directory, gfu_E, gfu_H)

            # # alternative eigenvalue solver, but careful, mode numbering may change
            # u = GridFunction(fes, multidim=30, name='resonances')
            # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
            #                         list(u.vecs), shift=300)
            # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

            # save json file
            shape = {'IC': cell_par}
            # print(run_save_directory)
            with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
                json.dump(shape, f, indent=4, separators=(',', ': '))

            qois = self.evaluate_qois(cav_geom, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)

            with open(fr'{run_save_directory}\qois.json', "w") as f:
                json.dump(qois, f, indent=4, separators=(',', ': '))
            return True
        else:
            error('Could not run eigenmode analysis due to error in geometry.')
            return False

    def vhf_gun(self):
        # variables
        # y1, T2, R2, L3, R4, L5, R6, R8, T9, R9, T10, R10, L11, R12, L13, R14, G, L_bp
        pass

    @staticmethod
    def eigen3d(geometry_dir):
        # define geometry
        cav_geom = pd.read_csv('D:\Dropbox\multipacting\MPGUI21\geodata.n',
                               header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0]]

        pnts = list(cav_geom.itertuples(index=False, name=None))
        wp = WorkPlane()
        wp.MoveTo(*pnts[0])
        for p in pnts[1:]:
            wp.LineTo(*p)
        wp.Close().Reverse()
        face = wp.Face()

        # name the boundaries
        face.edges.Max(X).name = "r"
        face.edges.Max(X).col = (1, 0, 0)
        face.edges.Min(X).name = "l"
        face.edges.Min(X).col = (1, 0, 0)
        face.edges.Min(Y).name = "b"
        face.edges.Min(Y).col = (1, 0, 0)

        vol = face.Revolve(Axis((0, 0, 0), X), 360)
        Draw(vol)
        geo = OCCGeometry(vol, dim=3)

        # mesh
        ngmesh = geo.GenerateMesh(maxh=0.01)
        mesh = Mesh(ngmesh)

        # define finite element space
        fes = HCurl(mesh, order=1, dirichlet='default')

        u, v = fes.TnT()

        a = BilinearForm(y * curl(u) * curl(v) * dx).Assemble()
        m = BilinearForm(y * u * v * dx).Assemble()

        apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
        pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

        with TaskManager():
            a.Assemble()
            m.Assemble()
            apre.Assemble()
            freedof_matrix = a.mat.CreateSmoother(fes.FreeDofs())

            # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
            gradmat, fesh1 = fes.CreateGradient()

            gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
            math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
            math1[0, 0] += 1  # fix the 1-dim kernel
            invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

            # build the Poisson projector with operator Algebra:
            proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

            projpre = proj @ pre.mat
            evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=30, maxit=20, printrates=False);

        # print out eigenvalues
        # print("Eigenvalues")
        freq_fes = []
        for i, lam in enumerate(evals):
            freq_fes.append(c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6)
            # print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

        # plot results
        gfu = GridFunction(fes, multidim=len(evecs))
        for i in range(len(evecs)):
            gfu.vecs[i].data = evecs[i]

        # # alternative eigenvalue solver
        # u = GridFunction(fes, multidim=30, name='resonances')
        # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
        #                         list(u.vecs), shift=300)
        # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

    @staticmethod
    def evaluate_qois(cav_geom, n, Req, L, gfu_E, gfu_H, mesh, freq_fes, beta=1):
        w = 2 * pi * freq_fes[n] * 1e6

        # calculate Vacc and Eacc
        Vacc = abs(Integrate(gfu_E[n][0] * exp(1j * w / (beta * c0) * x), mesh, definedon=mesh.Boundaries('b')))
        Eacc = Vacc / (L * 1e-3 * 2 * n)
        # print(f"Vacc: {Vacc}")
        # print(f"Eacc: {Eacc}")

        # calculate U and R/Q
        U = 2 * pi * 0.5 * eps0 * Integrate(y * InnerProduct(gfu_E[n], gfu_E[n]), mesh)
        # Uh = 2 * pi * 0.5 * mu0 * Integrate(y * InnerProduct(gfu_H[n], gfu_H[n]), mesh)
        RoQ = Vacc ** 2 / (w * U)
        # print(f"U: {U}")
        # print(f"Uh: {Uh}")
        # print("R/Q: ", RoQ)

        # OLD GEOMETRY INPUT TYPE
        # calculate peak surface fields
        xpnts_surf = cav_geom[(cav_geom[0] > 0) & (cav_geom[1] > min(cav_geom[1])) & (cav_geom[1] < max(cav_geom[1]))]
        Esurf = [Norm(gfu_E[n])(mesh(xi, yi)) for indx, (xi, yi) in xpnts_surf.iterrows()]
        Epk = (max(Esurf))
        # print("Epk/Eacc", Epk / Eacc)

        # calculate peak surface fields
        xpnts = cav_geom[(cav_geom[0] > 0) & (cav_geom[1] > min(cav_geom[1])) & (cav_geom[1] < max(cav_geom[1]))]
        Hsurf = [Norm(gfu_H[n])(mesh(xi, yi)) for indx, (xi, yi) in xpnts.iterrows()]
        Hpk = (max(Hsurf))
        # print("Bpk/Eacc", mu0 * Hpk * 1e9 / Eacc)

        # # calculate peak surface fields
        # pec_boundary = mesh.Boundaries("default")
        # bel = [xx.vertices for xx in pec_boundary.Elements()]
        # bel_unique = list(set(itertools.chain(*bel)))
        # xpnts_surf = sorted([mesh.vertices[xy.nr].point for xy in bel_unique])
        # Esurf = [Norm(gfu_E[n])(mesh(xi, yi)) for xi, yi in xpnts_surf]
        # Epk = (max(Esurf))
        # # print("Epk/Eacc", Epk / Eacc)
        #
        # Hsurf = [Norm(gfu_H[n])(mesh(xi, yi)) for xi, yi in xpnts_surf]
        # Hpk = (max(Hsurf))
        # # print("Bpk/Eacc", mu0 * Hpk * 1e9 / Eacc)

        # calculate surface power loss
        sigma_cond = 5.96e7  # <- conduction of copper
        Rs = np.sqrt(mu0 * w / (2 * sigma_cond))  # Surface resistance
        Ploss = 2 * pi * 0.5 * Rs * Integrate(y * InnerProduct(CF(Hsurf), CF(Hsurf)), mesh,
                                              definedon=mesh.Boundaries('default'))
        # print(f'Ploss: {Ploss}')

        # calculate cell to cell coupling factor
        f_diff = freq_fes[n] - freq_fes[1]
        f_add = (freq_fes[n] + freq_fes[1])
        kcc = 2 * f_diff / f_add * 100
        # print(f"kcc: {kcc}")

        # calculate Q
        Q = w * U / Ploss
        # print(f"Q: {Q}")

        # OLD
        # Get axis field
        xpnts_ax = np.linspace(min(cav_geom[1].tolist()), max(cav_geom[1].tolist()), 100)
        Eax = np.array([Norm(gfu_E[n])(mesh(xi, 0.0)) for xi in xpnts_ax])

        # # Get axis field
        # xmin = face.vertices.Min(X)
        # xmax = face.vertices.Max(X)
        # xpnts_ax = np.linspace(xmin.p[0], xmax.p[0], 100)
        # Eax = np.array([Norm(gfu_E[n])(mesh(xi, 0.0)) for xi in xpnts_ax])

        # calculate field flatness
        peaks, _ = find_peaks(Eax)
        E_abs_peaks = Eax[peaks]
        # ff = min(E_abs_peaks)/max(E_abs_peaks) * 100
        ff = (1 - ((max(E_abs_peaks) - min(E_abs_peaks)) / np.average(E_abs_peaks))) * 100
        # print(f"ff: {ff}")

        # calculate G
        G = Q * Rs
        # print(f"G: {G}")

        qois = {
            "Req [mm]": Req,
            "Normalization Length [mm]": 2 * L,
            "N Cells": n,
            "freq [MHz]": freq_fes[n],
            "Q []": Q,
            "Vacc [MV]": Vacc * 1e-6,
            "Eacc [MV/m]": Eacc * 1e-6,
            "Epk [MV/m]": Epk * 1e-6,
            "Hpk [A/m]": Hpk,
            "Bpk [mT]": mu0 * Hpk * 1e3,
            "kcc [%]": kcc,
            "ff [%]": ff,
            "Rsh [MOhm]": RoQ * Q * 1e-6,
            "R/Q [Ohm]": RoQ,
            "Epk/Eacc []": Epk / Eacc,
            "Bpk/Eacc [mT/MV/m]": mu0 * Hpk * 1e9 / Eacc,
            "G [Ohm]": G,
            "GR/Q [Ohm^2]": G * RoQ
        }

        return qois

    @staticmethod
    def save_fields(project_folder, gfu_E, gfu_H):
        # save mode fields
        with open(f"{project_folder}/gfu_EH.pkl", "wb") as f:
            pickle.dump([gfu_E, gfu_H], f)

    @staticmethod
    def save_mesh(project_folder, mesh):
        # save mesh
        folder = os.path.dirname(project_folder)
        with open(f"{folder}/mesh.pkl", "wb") as f:
            pickle.dump(mesh, f)

    @staticmethod
    def load_fields(folder, mode):
        with open(f'{folder}/gfu_EH.pkl', "rb") as f:
            [gfu_E, gfu_H] = pickle.load(f)
        return gfu_E, gfu_H

    @staticmethod
    def load_mesh(folder):
        folder = os.path.dirname(folder)
        # pickle mesh and fields
        with open(f'{folder}/mesh.pkl', 'rb') as f:
            mesh = pickle.load(f)

        return mesh

    def plot_fields(self, folder, mode=1, plotter='matplotlib'):
        gfu_E, gfu_H = self.load_fields(folder, mode)

        if plotter == 'matplotlib':
            mesh = self.load_mesh(folder)
            mesh_points = mesh.vertices
            E_values_mesh = []
            mesh_points_grid = []
            for pp in mesh_points:
                mesh_points_grid.append(pp.point)
                E_values_mesh.append(Norm(gfu_E[mode])(mesh(*pp.point)))

            mesh_points_grid = np.array(mesh_points_grid)
            xxt, yyt = mesh_points_grid[:, 0], mesh_points_grid[:, 1]
            triang = tri.Triangulation(xxt, yyt)

            # plot only triangles with sidelength smaller some max_radius
            max_radius = 0.02
            triangles = triang.triangles
            # Mask unwanted triangles.
            xtri = xxt[triangles] - np.roll(xxt[triangles], 1, axis=1)
            ytri = yyt[triangles] - np.roll(yyt[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri ** 2 + ytri ** 2), axis=1)
            triang.set_mask(np.max(xtri, axis=1) > max_radius)

            plt.tricontourf(triang, E_values_mesh, cmap='jet')
            plt.gca().set_aspect('equal', 'box')
            plt.show()
        else:
            Draw(gfu_H)

    def plot_mesh(self, folder, plotter='matplotlib'):
        mesh = self.load_mesh(folder)

        if plotter == 'matplotlib':

            mesh_points = mesh.vertices
            mesh_points_grid = []
            for pp in mesh_points:
                mesh_points_grid.append(pp.point)

            mesh_points_grid = np.array(mesh_points_grid)
            xxt, yyt = mesh_points_grid[:, 0], mesh_points_grid[:, 1]

            triang = tri.Triangulation(xxt, yyt)

            # plot only triangles with sidelength smaller some max_radius
            max_radius = 0.02
            triangles = triang.triangles
            # Mask unwanted triangles.
            xtri = xxt[triangles] - np.roll(xxt[triangles], 1, axis=1)
            ytri = yyt[triangles] - np.roll(yyt[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri ** 2 + ytri ** 2), axis=1)
            triang.set_mask(maxi > max_radius)

            plt.triplot(triang, lw=1, c='k', zorder=4)
            plt.scatter(mesh_points_grid[:, 0], mesh_points_grid[:, 1])
            plt.gca().set_aspect('equal', 'box')
            plt.show()
        else:
            Draw(mesh)

    @staticmethod
    def createFolder(fid, projectDir, subdir='', filename=None, opt=False, pol='monopole'):
        if opt:
            ngsolvemevp_folder = 'Optimisation'
        else:
            ngsolvemevp_folder = 'NGSolveMEVP'
        # change save directory
        path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}'
        if pol.lower() == 'monopole':
            pol_subdir = 'monopole'
        else:
            pol_subdir = 'dipole'

        if subdir == '':
            new_path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}/{pol_subdir}'
            if os.path.exists(new_path):
                path = new_path
            else:
                if not os.path.exists(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}'):
                    os.mkdir(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}')

                os.mkdir(new_path)
                path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}/{pol_subdir}'
        else:
            new_path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}/{fid}/{pol_subdir}'
            if os.path.exists(new_path):
                path = new_path
            else:
                if not os.path.exists(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}'):
                    os.mkdir(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}')
                    os.mkdir(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}/{fid}')

                if not os.path.exists(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}/{fid}'):
                    os.mkdir(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}/{fid}')

                os.mkdir(new_path)
                path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{subdir}/{fid}/{pol_subdir}'

        if os.path.exists(path):
            try:
                shutil.rmtree(path)
                os.mkdir(path)
            except PermissionError:
                pass
        else:
            os.mkdir(path)

    def gauss(self, n=11, sigma=1.0, shift=0.0):
        r = np.linspace(-int(n / 2) + 0.5, int(n / 2) - 0.5, n)
        g = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(r - shift) ** 2 / (2 * sigma ** 2))
        return g / max(g)

    def gaussian_deform(self, n_cells, surface, deformation_params):

        # plt.plot(surface[:, 0], surface[:, 1], label='undeformed')#, marker='o', mfc='none'

        deform_matrix = np.identity(len(surface))
        deform_vector = np.zeros((len(surface), 1))
        for n_cell in range(n_cells):
            # get index of ncell surface nodes
            n_cell_surf_indx = np.where(surface[:, 2] == n_cell + 1)[0]

            mn, mx = np.min(n_cell_surf_indx), np.max(n_cell_surf_indx)

            # disp = deformation_params[n_cell*3] * 0.02
            disp = deformation_params[n_cell * 3]
            # scale sigma to between 10 and 20% of length of cavity cell
            sigma = (0.1 + (0.2 - 0.1) * (1 + deformation_params[n_cell * 3 + 1]) / 2) * len(n_cell_surf_indx)
            shift = deformation_params[n_cell * 3 + 2] * (mn + mx) / 2
            # n_cell_deform_matrix = np.identity(len(surface)) + disp * np.diag(self.gauss(len(surface), sigma, shift=shift))
            # n_cell_deform_surface = n_cell_deform_matrix @ surface
            # plt.plot(n_cell_deform_surface[:, 0], n_cell_deform_surface[:, 1], label=n_cell)#, marker='o', mfc='none'
            # n_cell_deform_vector = np.atleast_2d(disp*1e-3*self.gauss(len(surface), sigma, shift=shift)).T
            # n_cell_deform_surface = n_cell_deform_vector + surface
            # plt.plot(n_cell_deform_surface[:, 0], n_cell_deform_surface[:, 1], label=n_cell)#, marker='o', mfc='none'

            # deform_matrix += disp * np.diag(self.gauss(len(surface), sigma, shift=shift))
            deform_vector += np.atleast_2d(disp * 2e-3 * self.gauss(len(surface), sigma, shift=shift)).T

        # deform surface
        # surface_def = deform_matrix @ surface
        surface_def = deform_vector + surface
        surface_def[:, 2] = surface[:, 2]

        # enforce end plane to be equal by checking that the first two and last two points are on same plane
        surface_def[1, 0] = surface_def[0, 0]
        surface_def[-1, 0] = surface_def[-2, 0]
        surface_def[0, 1] = 0
        surface_def[-1, 1] = 0

        # plt.plot(surface_def[:, 0], surface_def[:, 1], label=f'gaussian: {disp:.4f}, {sigma:.4f}, {shift:.4f}')  # , marker='o', mfc='none'

        # plt.legend()
        # # plt.gca().set_aspect('equal')
        # plt.show()
        return pd.DataFrame(surface_def, columns=[1, 0, 2])


if __name__ == '__main__':
    eig = NGSolveMEVP()

    mid_cell = [73.52, 131.75, 106.25, 118.7, 150, 187, 369.6321578127116]
    l_end_cell = [73.52, 131.75, 106.25, 118.7, 150, 187, 369.6321578127116]
    r_end_cell = [73.52, 131.75, 106.25, 118.7, 150, 187, 369.6321578127116]
    eig.cavity(2, 1, mid_cell, l_end_cell, r_end_cell)
