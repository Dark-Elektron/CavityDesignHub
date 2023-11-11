import itertools
import json
import shutil

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
        pass

    @staticmethod
    def ell_mdpt_grad(x_center, y_center, a, b, step, start, end, plot=True):
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

        # get inbox midpoint
        midpoint = inbox[int(len(inbox) / 2)]
        grad_midpoint = b ** 2 / (2 * (midpoint[1] - k)) - (b ** 2 / a ** 2) * ((midpoint[0] - h) / (midpoint[1] - k))
        #     print(midpoint, grad_midpoint)

        # plt.plot(inbox[:, 0], inbox[:, 1])
        # x = np.array([midpoint[0] * 0.9, midpoint[0] * 1.1])
        # y = grad_midpoint * x + midpoint[1] - grad_midpoint * midpoint[0]
        # plt.plot(x, y)
        # plot box
        # x1, y1 = start
        # x2, y2 = end
        # plt.plot([x1, x2], [y1, y1], 'b')
        # plt.plot([x1, x1], [y1, y2], 'b')
        # plt.plot([x1, x2], [y2, y2], 'b')
        # plt.plot([x2, x2], [y1, y2], 'b')
        # plt.scatter(midpoint[0], midpoint[1], s=10, c='r')

        return tuple(midpoint), grad_midpoint
    
    def write_geometry(self, folder, n_cells, mid_cell, end_cell_left=None, end_cell_right=None,
                       beampipe='none', plot=False):
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
        write_cavity_for_custom_eig_solver(file_path, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe, plot)

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

        A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(mid_cell[:7])*1e-3
        A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(end_cell_left[:7])*1e-3
        A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(end_cell_right[:7])*1e-3

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
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                wp.Spline([mdpt, (-shift + x1el, y1el)],
                          tangents={0: (1, 0), 1: (1, gd_mdpt), 2: (x2el - x1el, y2el - y1el)})  # <- draw iris arc

                wp.LineTo(*[-shift + x2el, y2el])  # <- draw tangent line

                pt = [-shift + x2el, y2el]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt,
                                      [L_bp_l + L_el - shift, Req])
                wp.Spline([mdpt, (L_bp_l + L_el - shift, Req)],
                          tangents={0: (x2el - x1el, y2el - y1el), 1: (1, gd_mdpt), 2: (1, 0)})  # <- draw left eq. arc
                pt = [L_bp_l + L_el - shift, Req]

                if n_cell == 1:
                    mdpt, gd_mdpt = self.ell_mdpt_grad(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                          [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req])
                    wp.Spline([mdpt, (L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er)],
                              tangents={0: (1, 0), 1: (1, gd_mdpt),
                                        2: (x2er - x1er, y1er - y2er)})  # <- draw left eq. arc

                    wp.LineTo(*[L_el + L_er - x1er + + L_bp_l + L_bp_r - shift, y1er])  # <- draw tangent line

                    pt = [L_el + L_er - x1er + + L_bp_l + L_bp_r - shift, y1er]
                    mdpt, gd_mdpt = self.ell_mdpt_grad(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                          [L_bp_l + L_el + L_er - shift, y1er])
                    wp.Spline([mdpt, (L_bp_l + L_el + L_er - shift, Ri_er)],
                              tangents={0: (x2er - x1er, y1er - y2er), 1: (1, gd_mdpt),
                                        2: (1, 0)})  # <- draw left eq. arc
                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    mdpt, gd_mdpt = self.ell_mdpt_grad(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                          [L_el + L_m - x2 + 2 * L_bp_l - shift, Req])
                    wp.Spline([mdpt, (L_el + L_m - x2 + 2 * L_bp_l - shift, y2)],
                              tangents={0: (1, 0), 1: (1, gd_mdpt), 2: (x2 - x1, y1 - y2)})  # <- draw left eq. arc

                    wp.LineTo(*[L_el + L_m - x1 + 2 * L_bp_l - shift, y1])  # <- draw tangent line

                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    mdpt, gd_mdpt = self.ell_mdpt_grad(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                          [L_bp_l + L_el + L_m - shift, y1])
                    wp.Spline([mdpt, (L_bp_l + L_el + L_m - shift, Ri_m)],
                              tangents={0: (x2 - x1, y1 - y2), 1: (1, gd_mdpt), 2: (1, 0)})  # <- draw left eq. arc
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]

                    # calculate new shift
                    shift = shift - (L_el + L_m)

            elif n > 1 and n != n_cell:
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                wp.Spline([mdpt, (-shift + x1, y1)],
                          tangents={0: (1, 0), 1: (1, gd_mdpt), 2: (x2 - x1, y2 - y1)})  # <- draw iris arc

                wp.LineTo(*[-shift + x2, y2])  # <- draw tangent line

                pt = [-shift + x2, y2]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req])
                wp.Spline([mdpt, (L_bp_l + L_m - shift, Req)],
                          tangents={0: (x2 - x1, y2 - y1), 1: (1, gd_mdpt), 2: (1, 0)})  # <- draw left eq. arc

                pt = [L_bp_l + L_m - shift, Req]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                      [L_m + L_m - x2 + 2 * L_bp_l - shift, Req])
                wp.Spline([mdpt, (L_m + L_m - x2 + 2 * L_bp_l - shift, y2)],
                          tangents={0: (1, 0), 1: (1, gd_mdpt), 2: (x2 - x1, y1 - y2)})  # <- draw left eq. arc

                wp.LineTo(*[L_m + L_m - x1 + 2 * L_bp_l - shift, y1])  # <- draw tangent line

                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                      [L_bp_l + L_m + L_m - shift, y1])
                wp.Spline([mdpt, (L_bp_l + L_m + L_m - shift, Ri_m)],
                          tangents={0: (x2 - x1, y1 - y2), 1: (1, gd_mdpt), 2: (1, 0)})  # <- draw left eq. arc
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]

                # calculate new shift
                shift = shift - 2 * L_m
            else:
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                wp.Spline([mdpt, (-shift + x1, y1)],
                          tangents={0: (1, 0), 1: (1, gd_mdpt), 2: (x2 - x1, y2 - y1)})  # <- draw iris arc

                wp.LineTo(*[-shift + x2, y2])  # <- draw tangent line

                pt = [-shift + x2, y2]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req])
                wp.Spline([mdpt, (L_bp_l + L_m - shift, Req)],
                          tangents={0: (x2 - x1, y2 - y1), 1: (1, gd_mdpt), 2: (1, 0)})  # <- draw left eq. arc

                pt = [L_bp_l + L_m - shift, Req]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                      [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req])
                wp.Spline([mdpt, (L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er)],
                          tangents={0: (1, 0), 1: (1, gd_mdpt), 2: (x2er - x1er, y1er - y2er)})  # <- draw left eq. arc

                wp.LineTo(*[L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er])  # <- draw tangent line

                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                mdpt, gd_mdpt = self.ell_mdpt_grad(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                      [L_bp_l + L_m + L_er - shift, y1er])
                wp.Spline([mdpt, (L_bp_l + L_m + L_er - shift, Ri_er)],
                          tangents={0: (x2er - x1er, y1er - y2er), 1: (1, gd_mdpt), 2: (1, 0)})  # <- draw left eq. arc
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
        #     Draw(face)

        return face

    def cavity(self, no_of_cells=1, no_of_modules=1, mid_cells_par=None, l_end_cell_par=None, r_end_cell_par=None,
               fid=None, bc=33, pol='Monopole', f_shift='default', beta=1, n_modes=None, beampipes='None',
               parentDir=None, projectDir=None, subdir='', expansion=None, expansion_r=None, mesh=None, opt=False):
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

        # change save directory
        if opt:  # consider making better. This was just an adhoc fix
            run_save_directory = projectDir / fr'SimulationData/NGSolveMEVP_opt/{fid}'
        else:
            # change save directory
            if subdir == '':
                run_save_directory = projectDir / fr'SimulationData/NGSolveMEVP/{fid}'
            else:
                run_save_directory = projectDir / fr'SimulationData/NGSolveMEVP/{fid}/{subdir}'
        # # write
        # self.write_geometry(run_save_directory, no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)

        # # define geometry
        # cav_geom = pd.read_csv(f'{run_save_directory}\geodata.n',
        #                        header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0]]
        #
        # pnts = list(cav_geom.itertuples(index=False, name=None))
        # wp = WorkPlane()
        # wp.MoveTo(*pnts[0])
        # for p in pnts[1:]:
        #     wp.LineTo(*p)
        # wp.Close().Reverse()
        # face = wp.Face()

        # # get face from ngsolve geom
        face = self.cavgeom_ngsolve(no_of_cells, mid_cells_par, l_end_cell_par, r_end_cell_par, beampipes)

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
        maxh = L/mesh[0]*1e-3
        # print(maxh)
        ngmesh = geo.GenerateMesh(maxh=maxh)
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
            evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=no_of_cells + 3, maxit=mesh[1], printrates=False)

        # print out eigenvalues
        # print("Eigenvalues")
        freq_fes = []
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

        # # alternative eigenvalue solver
        # u = GridFunction(fes, multidim=30, name='resonances')
        # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
        #                         list(u.vecs), shift=300)
        # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

        # save json file
        shape = {'IC': update_alpha(mid_cells_par),
                 'OC': update_alpha(l_end_cell_par),
                 'OC_R': update_alpha(r_end_cell_par)}

        with open(Path(fr"{run_save_directory}/geometric_parameters.json"), 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))

        qois = self.evaluate_qois(face, no_of_cells, Req, L, gfu_E, gfu_H, mesh, freq_fes)
        ic(qois)

        with open(fr'{run_save_directory}\qois.json', "w") as f:
            json.dump(qois, f, indent=4, separators=(',', ': '))

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
    def evaluate_qois(face, n, Req, L, gfu_E, gfu_H, mesh, freq_fes, beta=1):
        w = 2 * pi * freq_fes[n] * 1e6

        # calculate Vacc and Eacc
        Vacc = abs(Integrate(gfu_E[n][0] * exp(1j * w / (beta * c0) * x), mesh, definedon=mesh.Boundaries('b')))
        Eacc = Vacc / (L*1e-3 * 2 * n)
        # print(f"Vacc: {Vacc}")
        # print(f"Eacc: {Eacc}")

        # calculate U and R/Q
        U = 2 * pi * 0.5 * eps0 * Integrate(y * InnerProduct(gfu_E[n], gfu_E[n]), mesh)
        Uh = 2 * pi * 0.5 * mu0 * Integrate(y * InnerProduct(gfu_H[n], gfu_H[n]), mesh)
        RoQ = Vacc ** 2 / (w * U)
        # print(f"U: {U}")
        # print(f"Uh: {Uh}")
        # print("R/Q: ", RoQ)

        # OLD GEOMETRY INPUT TYPE
        # # calculate peak surface fields
        # xpnts_surf = cav_geom[(cav_geom[0] > 0) & (cav_geom[1] > min(cav_geom[1])) & (cav_geom[1] < max(cav_geom[1]))]
        # Esurf = [Norm(gfu_E[n])(mesh(xi, yi)) for indx, (xi, yi) in xpnts_surf.iterrows()]
        # Epk = (max(Esurf))
        # # print("Epk/Eacc", Epk / Eacc)
        #
        # # calculate peak surface fields
        # xpnts = cav_geom[(cav_geom[0] > 0) & (cav_geom[1] > min(cav_geom[1])) & (cav_geom[1] < max(cav_geom[1]))]
        # Hsurf = [Norm(gfu_H[n])(mesh(xi, yi)) for indx, (xi, yi) in xpnts.iterrows()]
        # Hpk = (max(Hsurf))
        # # print("Bpk/Eacc", mu0 * Hpk * 1e9 / Eacc)

        # calculate peak surface fields
        pec_boundary = mesh.Boundaries("default")
        bel = [xx.vertices for xx in pec_boundary.Elements()]
        bel_unique = list(set(itertools.chain(*bel)))
        xpnts_surf = sorted([mesh.vertices[xy.nr].point for xy in bel_unique])
        Esurf = [Norm(gfu_E[n])(mesh(xi, yi)) for xi, yi in xpnts_surf]
        Epk = (max(Esurf))
        # print("Epk/Eacc", Epk / Eacc)

        Hsurf = [Norm(gfu_H[n])(mesh(xi, yi)) for xi, yi in xpnts_surf]
        Hpk = (max(Hsurf))
        # print("Bpk/Eacc", mu0 * Hpk * 1e9 / Eacc)

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
        # # Get axis field
        # xpnts_ax = np.linspace(min(cav_geom[1].tolist()), max(cav_geom[1].tolist()), 100)
        # Eax = np.array([Norm(gfu_E[n])(mesh(xi, 0.0)) for xi in xpnts_ax])

        # Get axis field
        xmin = face.vertices.Min(X)
        xmax = face.vertices.Max(X)
        xpnts_ax = np.linspace(xmin.p[0], xmax.p[0], 100)
        Eax = np.array([Norm(gfu_E[n])(mesh(xi, 0.0)) for xi in xpnts_ax])

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
            "Normalization Length [mm]": 2*L,
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
    def createFolder(fid, projectDir, subdir='', filename=None, opt=False, pol='monopole'):
        if opt:
            ngsolvemevp_folder = 'NGSolveMEVP_opt'
        else:
            ngsolvemevp_folder = 'NGSolveMEVP'
        # change save directory
        path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}'
        if pol.lower() == 'monopole':
            subdir = 'monopole'
        else:
            subdir = 'dipole'

        if subdir == '':
            pass
        else:
            new_path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}/{subdir}'
            if os.path.exists(new_path):
                path = new_path
            else:
                if not os.path.exists(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}'):
                    os.mkdir(projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}')

                os.mkdir(new_path)
                path = projectDir / fr'SimulationData/{ngsolvemevp_folder}/{fid}/{subdir}'

        if os.path.exists(path):
            try:
                shutil.rmtree(path)
                os.mkdir(path)
            except PermissionError:
                pass
        else:
            os.mkdir(path)


if __name__ == '__main__':
    eig = NGSolveMEVP()

    mid_cell = [73.52, 131.75, 106.25, 118.7, 150, 187, 369.6321578127116]
    l_end_cell = [73.52, 131.75, 106.25, 118.7, 150, 187, 369.6321578127116]
    r_end_cell = [73.52, 131.75, 106.25, 118.7, 150, 187, 369.6321578127116]
    eig.cavity(2, 1, mid_cell, l_end_cell, r_end_cell)
