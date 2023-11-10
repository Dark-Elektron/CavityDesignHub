import shutil
import subprocess
from sys import exit
from matplotlib import colors
from scipy.signal import find_peaks
from utils.shared_functions import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.io as spio
import pandas as pd
from icecream import ic

m0 = 9.1093879e-31
q0 = 1.6021773e-19
mu0 = 4 * np.pi * 1e-7
eps0 = 8.85418782e-12
c0 = 2.99792458e8
eta0 = 376.7303134111465


class Model:
    def __init__(self, folder=None, name=None, parent_dir=None):
        self.name = name
        self.Epk, self.Hpk, self.Eacc, self.Vacc, self.E_z_axis = None, None, None, None, None
        self.k_loss, self.epk, self.bpk = None, None, None
        self.ff, self.kcc, self.Q0, self.U, self.Pds = None, None, None, None, None
        self.Rs, self.R_Q, self.G, self.GR_Q = None, None, None, None
        self.eigen_frequencies = None
        self.beampipe = None
        self.mid_cell = None
        self.end_cell_left = None
        self.end_cell_right = None
        self.n_cells = None
        self.eig_freq = None
        self.fields = None
        self.job = None
        self.mesh = None
        self.fieldparam = None
        self.param = None
        self.folder = fr"{folder}\{name}"
        self.geodata = None
        self.parent_dir = parent_dir

    def set_name(self, name):
        self.name = name

    def set_folder(self, folder):
        self.folder = folder

    def create_inputs(self):
        """
        Create required inputs

        Returns
        -------

        """
        # remove old files
        files_list = ["geodata.n", "initials", "flevel", "y00", "param", "fieldparam", "counter_flevels.mat",
                      "counter_initials.mat",
                      "gene_initials.mat", "gene_temp.mat", "initials_temp.mat", "distance_flevel.mat",
                      "counter_initialsl.mat", "gene_initialsl.mat", "initials_templ.mat", "distance_flevell.mat",
                      "counter_initialsr.mat", "gene_initialsr.mat", "initials_tempr.mat", "distance_flevelr.mat",
                      "dx.mat", "alphas.mat"]
        # input_files_list = ["initials", "param", "fieldparam", "counter_flevels.mat",
        #                     "counter_initials.mat",
        #                     "gene_initials.mat"]
        for file in files_list:
            if os.path.exists(fr"{self.folder}\{file}"):
                os.remove(fr"{self.folder}\{file}")

        # write intset. sample points and set for integration
        if not os.path.exists(fr"{self.folder}\intset.mat"):
            PP, ww = self.inttri(6)
            intset = {"PP": PP,
                      "ww": ww}
            spio.savemat(f"{self.folder}/intset.mat", intset, format='4')

        # check if geometry has been written
        if self.geodata is None:
            ic("The geometry file geodata.n seems to be missing")
            return

        ic('Old input data deleted.')

        # Inputs for eigenmode solver
        gtype = 1
        epsr = 1  # relative epsilon
        strech = 0
        d1 = 1  # Grid constant
        d2 = 1  # Mesh constant
        R = -1 + 0 * 1j
        v0 = 2  # initial velocity
        N = 20  # Number of impacts
        Z = 0
        freq = 0

        # initialize fieldparam
        self.fieldparam = [gtype, freq, epsr, d1 / 1000, R.real, R.imag, strech, Z]
        df = pd.DataFrame(self.fieldparam)
        df.to_csv(fr'{self.folder}\fieldparam', index=False, header=False, float_format='%.7E')

        # parameters for the MP analysis
        param = np.zeros((7, 1))
        V0 = 1  # intensity of the EM field
        # gamma = c0 / freq / (2 * np.pi)  # distance/phase -coefficient
        v00 = c0 * np.sqrt(1 - 1 / ((v0 * q0 / (m0 * c0 ** 2) + 1) ** 2))  # initial velocity (relativistic)
        # v00   = sqrt(2*v0*q/m)                 # (classical)
        ctype = 1  # compute counter functions
        tol = 1e-3  # tolerance for the ODE solver
        emin = 1  # not required for eigenmode analysis
        emax = 10  # not required for eigenmode analysis
        # self.param = [freq, V0, gamma, v00, N, ctype, tol, emin, emax]
        # save -ascii param param
        # df = pd.DataFrame(self.param)
        # df.to_csv(fr'{self.folder}\param', index=False, header=False, float_format='%.7E')

        # grid constant for creating a new grid
        # geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        self.geodata[0, 0] = d2 / 1000
        df = pd.DataFrame(self.geodata)
        df.to_csv(fr'{self.folder}\geodata.n', sep=r' ', index=False, header=False, float_format='%.7E')

    def define_geometry(self, n_cells, mid_cell, end_cell_left=None, end_cell_right=None,
                        beampipe='none', plot=False):
        """
        Define geometry

        Parameters
        ----------
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

        self.n_cells = n_cells
        self.mid_cell = mid_cell

        if end_cell_left is None:
            self.end_cell_left = mid_cell
        else:
            self.end_cell_left = end_cell_left

        if end_cell_right is None:
            if self.n_cells > 1:
                if end_cell_left is None:
                    self.end_cell_right = mid_cell
                else:
                    self.end_cell_right = end_cell_left
            else:
                self.end_cell_right = mid_cell
        else:
            self.end_cell_right = end_cell_right

        self.beampipe = beampipe

        # check if folder exists
        if not os.path.exists(self.folder):
            try:
                os.mkdir(self.folder)
            except FileNotFoundError:
                ic("There was a problem creating the directory for the simulation files. Please check folder.")
                exit()

        file_path = fr"{self.folder}\geodata.n"
        write_cavity_for_custom_eig_solver(file_path, n_cells, mid_cell, end_cell_left, end_cell_right, beampipe, plot)

        self.geodata = pd.read_csv(file_path, sep='\s+', header=None).to_numpy()

        self.create_inputs()

    def generate_mesh(self, gridcons=0.005, epsr=1, plot=False):
        """
        Generate mesh

        Parameters
        ----------
        gridcons: float
            Maximum element width
        epsr:
            Relative permittivity
        plot:
            Show mesh after generation or not

        Returns
        -------

        """

        # Then the edges, the third column secifies the type of an edge as follows:
        #   1 : nxE = 0 on boundary
        #   3 : nxH = 0 on boundary
        #   0 : for stretching (open boundary condition?)

        if self.geodata is None:
            ic("Geometry data geodata.n missing. No geometry to mesh")
            return

        ic("Meshing geometry")

        n = len(self.geodata[:, 1]) - 1
        gr = self.geodata[3:n, 0]
        gz = self.geodata[3:n, 1]
        gn = self.geodata[3:n, 2]
        n = len(gr)

        PEC = np.where(gn == 1)[0]  # electric walls
        WIN = np.where(gn == 2)[0]  # dielectric walls
        PMC = np.where(gn == 3)[0]  # magnetic walls

        if len(WIN) > 0:
            ic('Dielectric boundary found in a cavity.')
            return

        # specify the job for the eigenvalue solver
        if len(PMC) > 0:
            ic('A cavity with magnetic walls.')
            self.job = 1
        elif len(PEC) > 0 & len(PMC) == 0:
            ic('A cavity with electric walls.')
            self.job = 0
        else:
            ic('Invalid surface boundary conditions. Please check the geometry file.')
            return

        job = {"job": self.job}
        spio.savemat(f'{self.folder}/job.mat', job, format='4')

        nodes = np.array([gz, gr]).T
        edges = np.array([np.arange(1, n + 1), np.append(np.arange(2, n + 1), 1), gn]).T
        ne = len(edges[:, 0])

        # Patches: first the edges are listed, last three values give relative epsilon, relative mu and patch type,
        # 0 = can stretch, 1-n = regulars
        patches = np.append(np.arange(1, ne + 1), [epsr, 1, 1])
        esize = gridcons

        # save model.mat
        model = {'nodes': nodes,
                 "edges": edges,
                 "patches": patches,
                 "esize": esize
                 }
        spio.savemat(f"{self.folder}/model.mat", model, format='4')

        # start the mesh generator
        if self.parent_dir is not None:
            subprocess.call(fr"{self.parent_dir}\exe\own_exe\2dgen_bin.exe", cwd=self.folder,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # subprocess.call(fr"{self.parent_dir}\exe\own_exe\2dgen_bin.exe", cwd=self.folder)
        else:
            print(os.path.exists(fr"..\..\..\exe\own_exe\2dgen_bin.exe"), os.path.exists(self.folder))
            subprocess.call(fr"..\..\..\exe\own_exe\2dgen_bin.exe", cwd=self.folder,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # subprocess.call(f"{self.folder}/2dgen_bin.exe", cwd=self.folder)

        self.mesh = spio.loadmat(fr"{self.folder}\mesh.mat")
        ic("Meshing completed.")
        if plot:
            self.plot_mesh()

    def plot_mesh(self):
        """
        Plot generated mesh

        Returns
        -------

        """
        fig, ax = plt.subplots()
        if not os.path.exists(fr"{self.folder}\mesh.mat"):
            ic(['The mesh does not exists. There might have been an error in generating the mesh.'])
            return

        if not os.path.exists(fr"{self.folder}\fieldparam"):
            ic('Parameter file fieldparam does not exist.')
            return

        # gtype = self.fieldparam[0]
        # if gtype == 1:
        #     ic('\tPlotting the mesh.')
        # else:
        #     ic('\tPlotting the mesh. Blue area for streching.')

        # plots 2D mesh in mesh.mat, which includes coord and etopol -arrays
        coord = self.mesh['coord']
        etopol = self.mesh['etopol']
        alue = self.mesh['alue']
        tyyppi = self.mesh['tyyppi']
        boundary = self.mesh['boundary']
        edges = self.mesh['edges']

        x_min = min(coord[:, 0])
        x_max = max(coord[:, 0])
        y_min = min(coord[:, 1])
        y_max = max(coord[:, 1])
        ax.plot([x_min, x_max, x_max, x_min], [y_min, y_min, y_max, y_max], 'k.')

        koko, pois = np.shape(etopol)
        ala = 0.0
        X = np.zeros((4, koko))
        Y = np.zeros((4, koko))
        C = np.zeros((0, koko))
        for i in range(koko):
            n1 = int(etopol[i, 0]) - 1
            x1 = coord[n1, 0]
            y1 = coord[n1, 1]
            n2 = int(etopol[i, 1]) - 1
            x2 = coord[n2, 0]
            y2 = coord[n2, 1]
            n3 = int(etopol[i, 1]) - 1
            x3 = coord[n3, 0]
            y3 = coord[n3, 1]
            X[:, i] = np.array([x1, x2, x3, x1]).T
            Y[:, i] = np.array([y1, y2, y3, y1]).T

            osa = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
            ala = ala + osa
            ax.plot([x1, x2, x3, x1], [y1, y2, y3, y1], 'k')

        mask = np.where(alue.T[0] == 0)[0]
        ax.fill(X[:, mask], Y[:, mask], 'b')
        mask = np.where(tyyppi.T[0] == 1)[0]
        ax.fill(X[:, mask], Y[:, mask], 'r')
        mask = np.where(tyyppi.T[0] == 2)[0]
        ax.fill(X[:, mask], Y[:, mask], 'g')

        mask = np.where(edges[:, 2] > 0)[0]  # no bouncing
        for i in range(len(mask)):
            i1 = int(edges[mask[i], 0]) - 1
            i2 = int(edges[mask[i], 1]) - 1

            ax.plot([coord[i1, 0], coord[i2, 0]], [coord[i1, 1], coord[i2, 1]], 'r')

        mask = np.where(boundary == 3)
        for i in range(len(mask)):
            ax.plot(coord[mask[i], 0], coord[(mask[i] - 1, 1)], 'b*')

        mask = np.where(boundary == 0)[0]
        ax.plot(coord[mask, 0], coord[mask, 1], 'w*')
        ala = ala / 2
        ax.set_xlabel('z axis [m]')
        ax.set_ylabel('r axis [m]')

        plt.show()

    def run_field_solver(self, n_modes=None, freq=None, req_mode_num=None, show_plots=True, search=True):
        # Function program cavity_field(s)
        # -----------------------------------------------------------------------
        # Runs the field solver for computing the EM fields in a cavity with
        # magnetic ends.
        # INPUT  s : 1 - display messages, 0 - display only a few messages
        #
        # -----------------------------------------------------------------------
        # CALLS TO : print.m, make_model_cavity.m, plot_mesh.m, eigen.m,
        #            calculate_fields.m, plot_FEM_fields.m, clear_window.m
        # 10/04/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
        # ------------------------------------------------------------------------

        if self.geodata is None:
            return

        if self.mesh is None:
            ic("No mesh found. Generating mesh with default settings.")
            self.generate_mesh(plot=show_plots)

        # load geodata.n
        # geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()

        # if self.epsr > 1:
        #     ic('Relative permittivity > 1?')
        ic(n_modes)
        if freq is None:
            # calculate freq from mid cell length
            beta = 1
            freq = beta*c0/(4*self.mid_cell[5])
            ic("calculated freq from length: ", freq)

        eigen_freq = self.eigen(n_modes, freq, req_mode_num, search)

        if search:
            if freq != 0:
                err = (abs(eigen_freq - freq) / freq) * 100
                if err > 1:
                    ic('Warning: Error in eigen frequency more than 1#.')
            else:
                ic("No frequency to compare with.")

        # compute and plot the fields
        self.calculate_fields()

        if show_plots:
            ic("Plotting fields")
            self.plot_FEM_fields('contour')

    def eigen(self, n_modes, freq, req_mode_num, search):
        """
        Perform eigenmode simulation

        Parameters
        ----------
        n_modes: int
            Number of modes
        freq: float
            Frequency to find eigenvalues around
        req_mode_num: int
            Required mode number
        search: null
            Redundant for now

        Returns
        -------

        """
        if n_modes is None:
            n_modes = self.n_cells + 1

        if req_mode_num:
            if req_mode_num > n_modes:
                req_mode_num = self.n_cells
        else:
            req_mode_num = self.n_cells

        # maara = 10  Means number of modes in original code
        offset = {"offset": 0}
        spio.savemat(f"{self.folder}/offset.mat", offset, format='4')

        # !eigenC_bin
        cwd = fr'{self.folder}'
        if self.parent_dir is not None:
            eigenCpath = fr'{self.parent_dir}\exe\own_exe\eigenC_bin.exe'
        else:
            eigenCpath = fr'..\..\..\exe\own_exe\eigenC_bin.exe'

        if os.path.exists(eigenCpath):
            subprocess.call(eigenCpath, cwd=cwd,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # subprocess.call(eigenCpath, cwd=cwd)
        ic("\tDone running eigenC_bin.exe")

        cinfo = spio.loadmat(fr"{self.folder}\cinfo.mat")['cinfo'][0]
        if cinfo != 0:
            ic('\t\tEigenvalue solver failed.')

        # load o_eigen
        # I think eigen2_bin only writes out the matrices necessary for the computation of the eigenvalues and
        # does not actually solve the eigenvalue problem
        o_eigen = spio.loadmat(fr"{self.folder}\o_eigen.mat")
        n = int(o_eigen['n'])
        ia = o_eigen['ia'].T[0] - 1  # matlab indices to python
        ja = o_eigen['ja'].T[0] - 1
        Aarvot = o_eigen['Aarvot'].T[0]
        Barvot = o_eigen['Barvot'].T[0]
        nf = o_eigen['nf']
        index = o_eigen['index']
        siirto = o_eigen['siirto'][0]

        AA = sps.csr_matrix((Aarvot, (ia, ja)), shape=(n, n))
        BB = sps.csr_matrix((Barvot, (ia, ja)), shape=(n, n))

        k0 = 2 * np.pi * freq * np.sqrt(mu0 * eps0)
        d2, u2 = spsl.eigs(AA, M=BB, k=n_modes, sigma=k0**2)
        # imaginary component of eigenvectors are zero
        u2 = u2.real
        d2 = np.absolute(d2)
        k = np.sqrt(d2)
        if search:
            k2 = k
            eigen_frequencies = np.array(k / (2 * np.pi * np.sqrt(mu0 * eps0)))
            # sort array
            sort_ind = eigen_frequencies.argsort()
            self.eigen_frequencies = eigen_frequencies[sort_ind]
            k2 = k2[sort_ind]
            u2 = u2[:, sort_ind]

            if req_mode_num:
                ind = int(req_mode_num) - 1
                eigen_freq = self.eigen_frequencies[ind]
                k = k2[ind]
                u = u2[:, ind]
            else:
                # find eigenvalue with min error.
                val = np.min(abs(k0 - k2))
                ind = np.argmin(abs(k0 - k2))
                k = k2[ind]
                u = u2[:, ind]
                # new frequency
                eigen_freq = k / (2 * np.pi * np.sqrt(mu0 * eps0))
        else:
            self.eigen_frequencies = np.sort(k / (2 * np.pi * np.sqrt(mu0 * eps0)))
            eigen_freq = self.eigen_frequencies[int(req_mode_num)-1]
            u = u2[:, int(req_mode_num)-1]

        # param = pd.read_csv(fr"{self.folder}\param", sep='\s+', header=None).to_numpy().T[0]
        # self.param[0] = eigen_freq
        # df = pd.DataFrame(self.param)
        # df.to_csv(fr'{self.folder}\param', index=False, header=False, float_format='%.7E')

        # fieldparam = pd.read_csv(fr"{self.folder}\fieldparam", sep='\s+', header=None).to_numpy().T[0]
        self.eig_freq = eigen_freq
        self.fieldparam[1] = eigen_freq
        df = pd.DataFrame(self.fieldparam)
        df.to_csv(fr'{self.folder}\fieldparam', index=False, header=False, float_format='%.7E')

        kama0 = {"index": index,
                 "u": u,
                 "k": k
                 }

        spio.savemat(f"{self.folder}/kama0.mat", kama0, format='4')
        ic("\t\tDone with eigen")
        return eigen_freq

    def calculate_fields(self):
        """
        Calculate the fields from the eigenvectors

        Returns
        -------

        """
        ic('\t\tComputing the fields.')
        if self.geodata is None:
            ic("Geometry file geodata.n does not exist.")
            return

        if not os.path.exists(fr"{self.folder}\mesh.mat"):
            ic('\t\tThe mesh does not exists. Choose Mesh Generator in menu Run.')
            return

        if not os.path.exists(fr"{self.folder}\kama0.mat"):
            ic('\t\tThe fields do not exist. Choose Field Solver in menu Run.')
            return

        # geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        gridcons = self.geodata[0, 0]
        n = len(self.geodata[:, 0])
        gr = self.geodata[3:n, 0]
        gz = self.geodata[3:n, 1]

        # load mesh
        # mesh = spio.loadmat(fr"{self.folder}\mesh.mat")
        coord = self.mesh['coord']
        etopol = self.mesh['etopol']
        alue = self.mesh['alue']
        tyyppi = self.mesh['tyyppi']
        boundary = self.mesh['boundary']
        edges = self.mesh['edges']

        z1 = min(gz)
        z2 = max(gz)
        r1 = min(gr)
        r2 = max(gr)

        zmaara = max(10, min(200, len(np.arange(z1, z2, gridcons))))
        if np.mod(zmaara, 2) == 0:
            zmaara = zmaara + 1
        rmaara = max(10, min(200, len(np.arange(r1, r2, gridcons))))

        I1 = np.arange(z1, z2 + (z2 - z1) / (zmaara - 1),
                       (z2 - z1) / (zmaara - 1))  # (z2-z1)/(zmaara - 1) included so that end pint is included
        I2 = np.arange(r1, r2 + (z2 - z1) / (zmaara - 1),
                       (r2 - r1) / (rmaara - 1))  # (z2-z1)/(zmaara - 1) included so that end pint is included

        zz, rr = np.meshgrid(I1, I2)
        m = len(I1)
        n = len(I2)
        z = zz.T.flatten()
        r = rr.T.flatten()

        alue = 0
        ax = np.where(r == 0)
        r[ax] = r[ax] + 1e-10
        ind = np.arange(0, n) * np.ones((m, n))
        rind = ind.flatten()
        ind2 = (np.ones((n, m)) * np.arange(0, m)).T
        zind = ind2.flatten()

        zr = {'z': z, 'r': r, 'alue': alue, 'rind': rind, 'zind': zind}
        spio.savemat(fr'{self.folder}\zr.mat', zr, format='4')

        # compute the fields at generated points
        # !copy /y kama0.mat kama.mat
        shutil.copyfile(fr'{self.folder}\kama0.mat', fr'{self.folder}\kama.mat')

        cwd = fr'{self.folder}'
        if self.parent_dir is not None:
            multipacPath = fr'{self.parent_dir}\exe\own_exe\Multipac.exe'
        else:
            multipacPath = fr'..\..\..\exe\own_exe\Multipac.exe'
        if os.path.exists(multipacPath):
            subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd)

        Er = spio.loadmat(fr"{self.folder}\Er.mat")['Er']
        Ez = spio.loadmat(fr"{self.folder}\Ez.mat")['Ez']
        H = spio.loadmat(fr"{self.folder}\H.mat")['H']

        self.fields = {'I1': I1, 'I2': I2, 'rr': rr, 'zz': zz, 'r': r, 'z': z, 'Er': Er, 'Ez': Ez, 'H': H}
        spio.savemat(fr'{self.folder}\fields.mat', self.fields)  # this was not saved as Mat object version 4. Ask why?

        self.calculate_QoI()

    def plot_FEM_fields(self, ptype):
        ok2 = 0
        if os.path.exists(fr'{self.folder}/geodata.n'):
            files = ['fields.mat', 'Er.mat', 'Ez.mat', 'H.mat', 'fieldparam']
            for f in files:
                if os.path.exists(fr'{self.folder}/{f}'):
                    continue
                else:
                    ok2 = 1

        if ok2 == 0:
            if ptype == 'contour':
                ic('Plotting the fields. A pcolor plot.')
            elif ptype == 'arrow':
                ic('Plotting the fields. An arrow plot.')

            # geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
            self.plot_cavity_fields(ptype)

    def plot_cavity_fields(self, ptype):
        """
        Plots the computed fields on the input geometry

        Parameters
        ----------
        ptype: {'contour', 'arrow'}
            Plot type

        Returns
        -------

        """
        fig, axs = plt.subplots(2, 1)

        # geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        n = len(self.geodata[:, 0])
        gr = self.geodata[3:n, 0]
        gz = self.geodata[3:n, 1]

        # load the field values
        # load fields Er Ez H I1 I2 zz rr z r
        # fields = spio.loadmat(fr"{self.folder}\fields.mat")
        Er = self.fields['Er']
        Ez = self.fields['Ez']
        H = self.fields['H']
        I1 = self.fields['I1']
        I2 = self.fields['I2']
        zz = self.fields['zz']
        rr = self.fields['rr']
        z = self.fields['z']
        r = self.fields['r']

        # plot the fields
        if ptype == 'contour':  # a pcolor plot
            # ic(I2.shape, I1.shape, Er.shape)
            Er = np.reshape(Er, (max(I1.shape), max(I2.shape)))
            Ez = np.reshape(Ez, (max(I1.shape), max(I2.shape)))
            EE = np.sqrt(abs(Er) ** 2 + abs(Ez) ** 2).T

            pcE = axs[0].pcolor(zz, rr, EE, cmap='RdBu')
            # pcE = axs[0].contourf(zz, rr, EE, cmap='RdBu', levels=40)
            # color = 2 * np.log(np.hypot(Ez, Er))
            # pcE = axs[0].streamplot(zz, rr, Ez, Er, linewidth=1, cmap=plt.cm.inferno,
            #               density=2, arrowstyle='->', arrowsize=1.5)
            plt.colorbar(pcE, ax=axs[0])

            axs[0].plot(gz, gr, '-b')

            axs[0].set_xlabel('z axis [m]')
            axs[0].set_ylabel('r axis [m]')
            axs[0].set_title(r'Electric field |$E$| [V/m]')

            HH = np.reshape(H, (max(I1.shape), max(I2.shape))).T
            pcH = axs[1].pcolor(zz, rr, 4e-7 * np.pi * HH, cmap='RdBu')
            plt.colorbar(pcH, ax=axs[1])

            axs[1].plot(gz, gr, '-b')
            axs[1].set_title(r'Magnetic field |$B_\phi$| [T]')
            axs[1].set_xlabel('z axis [m]')
            axs[1].set_ylabel('r axis [m]')
            plt.show()
        else:
            EE = np.sqrt(abs(Er) ** 2 + abs(Ez) ** 2)
            fig, ax = plt.subplots()
            # arrow(F,A,gridcons,'f','s','k')
            ic(EE.min(), np.min(EE))
            skip = slice(None, None, 8)
            ax.quiver(z[skip], r[skip], Ez[skip], Er[skip], EE[skip], cmap='coolwarm',
                      norm=colors.LogNorm(vmin=np.min(EE) + 0.1, vmax=np.max(EE)),
                      minlength=3)

            ax.plot(gz, gr, '-b')
            ax.set_xlabel('z axis [m]')
            ax.set_ylabel('r axis [m]')
            ax.set_title(r'Electric field |$E$| [V/m]')
            plt.show()

    def calculate_QoI(self):
        ic("\t\tCalculating QoIs")

        # geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        n = len(self.geodata[:, 0])
        gr = np.array(self.geodata[3:n, 0])
        gz = np.array(self.geodata[3:n, 1])

        # load the field values
        # fields = spio.loadmat(fr"{self.folder}\fields.mat")
        Er = self.fields['Er']
        Ez = self.fields['Ez']
        H = self.fields['H']
        I1 = self.fields['I1']
        I2 = self.fields['I2']
        zz = self.fields['zz']
        rr = self.fields['rr']
        z = self.fields['z']  # flattened zz array
        r = self.fields['r']  # flattened rr array

        # plot the fields
        Er = np.reshape(Er, (max(I1.shape), max(I2.shape)))
        Ez = np.reshape(Ez, (max(I1.shape), max(I2.shape)))

        # check electric fields
        EE = np.sqrt(abs(Er) ** 2 + abs(Ez) ** 2).T
        HH = np.reshape(abs(H), (max(I1.shape), max(I2.shape))).T

        # Calculate accelerating field
        self.E_z_axis = Ez.T[rr == 0]  # & zz>=-L_half_cell & zz <= L_half_cell
        self.E_z_axis_abs = abs(self.E_z_axis)
        z_slice = zz[0, :]
        # z_active = z_slice[(z_slice >= -L_half_cell) & (z_slice <= L_half_cell)]

        self.Vacc = self.calculate_vacc(z_slice, self.E_z_axis)
        self.Eacc = self.Vacc / (2 * self.n_cells * self.mid_cell[5])

        z_e_surf, r_e_surf, E_surf = self.get_surface_field(EE, zz, rr)
        z_h_surf, r_h_surf, H_surf = self.get_surface_field(HH, zz, rr)

        self.Epk = np.max(E_surf)
        self.Hpk = np.max(H_surf)

        self.epk = self.calculate_epk(self.Epk, self.Eacc)
        self.bpk = self.calculate_bpk(self.Hpk, self.Eacc)

        # plt.plot(np.array(z_e), E_surf/np.max(E_surf))
        # plt.plot(np.array(z_h), H_surf/np.max(H_surf))
        # plt.show()

        # calculate energy
        self.U = self.calculate_U(EE, HH, zz, rr)
        # calculate R/Q
        self.R_Q = self.calculate_rq(self.Vacc, self.U)
        # calculate loss factor
        self.k_loss = self.calculate_loss_factor(self.R_Q)
        # calculate surface resistance using copper conductivity
        self.Rs = self.calculate_surface_resistance_copper()
        self.Pds = self.calculate_disspated_power(H_surf, z_h_surf, r_h_surf)
        self.Q0 = self.calculate_q()
        self.G = self.calculate_g()

        # ic(self.eig_freq, self.epk, self.bpk, self.R_Q, self.k_loss,
        #    self.Rs, self.Pds, self.Q0, self.G)

        self.save_qois()
        # write_cst_paramters(self.name, self.mid_cell, self.end_cell_left, self.end_cell_right,
        #                     projectDir=self.folder, cell_type="None")

    @staticmethod
    def calculate_U(EE, HH, zz, rr):
        """
        Calculates the energy :math:`U_n` of eigenmode :math:`n` from either electric or magnetic field.

        .. math::
           U_n = \\frac{1}{2} \int_{z_\min}^{z_\max} \int_{r_\min}^{r_\max} \int_{0}^{2\pi} D_n \cdot E_n \mathrm{d}V
           = \\frac{1}{2}\int_{z_\min}^{z_\max} \int_{r_\min}^{r_\max} \int_{0}^{2\pi} B_n \cdot H_n \mathrm{d}V

           D_n \cdot E_n = \epsilon_0 E^2; ~B_n \cdot H_n = \mu_0 H^2

        Parameters
        ----------
        EE: 2D array
            Array of electric field magnitude in 2D space
        HH:
            Array of magnetic field magnitude in 2D space
        zz: 2D array
            Grid points z coordinate
        rr: 2D array
            Grid points r coordinate

        Returns
        -------
        U: float
            Energy over the shape volume

        """

        # calculate volume which is area under contour curve
        z = zz[0, :]
        r = rr[:, 0]
        # ic(z, r)

        U = np.pi * mu0 * np.trapz(np.trapz(HH ** 2 * rr, r, axis=0), z)
        UE = np.pi * eps0 * np.trapz(np.trapz(EE ** 2 * rr, r, axis=0), z)
        ic(U, UE)

        return U

    @staticmethod
    def calculate_lorentz_force(EE, HH, zz, rr):
        """
        Calculates the pressure  :math:`p` induced on the cavity walls due to the Lorentz force of
        eigenmode :math:`n` from the electric or magnetic field.

        .. math::
           p = \\frac{1}{4} \mu_0H^2 - \epsilon_0 E^2

        Parameters
        ----------
        EE: 2D array
            Array of electric field magnitude in 2D space
        HH:
            Array of magnetic field magnitude in 2D space
        zz: 2D array
            Grid points z coordinate
        rr: 2D array
            Grid points r coordinate

        Returns
        -------
        U: float
            Energy over the shape volume

        """

        # calculate volume which is area under contour curve
        z = zz[0, :]
        r = rr[:, 0]
        # ic(z, r)

        U = np.pi * mu0 * np.trapz(np.trapz(HH ** 2 * rr, r, axis=0), z)
        UE = np.pi * eps0 * np.trapz(np.trapz(EE ** 2 * rr, r, axis=0), z)
        ic(U, UE)

        return U

    @staticmethod
    def frequency_shift(EE, HH, zz, rr, U):
        pass

    def calculate_surface_resistance_copper(self):
        """
        Calculates the surface resistance for Copper with conductivity

        of :math:`\sigma = 5.96 \\times 10^7 \mathrm{S/m}` at eigen frequency

        .. math::
           R_\mathrm{s} = \\sqrt{\\frac{\mu \omega}{2 \sigma}}

        Returns
        -------
        Rs: float
            Surface resistance of copper at operating frequency.
        """
        sigma_Cu = 5.96e7
        Rs = np.sqrt((mu0 * 2 * np.pi * self.eig_freq) / (2 * sigma_Cu))

        return Rs

    def calculate_disspated_power(self, H_surf, z_surf, r_surf):
        """
        Calculates the dissipated power from the cavity walls

        .. math::
           P_\mathrm{ds} = \\frac{1}{2} R_\mathrm{s} \int_{z_\min}^{z_\max} \int_{0}^{2\pi} |H_\\theta|^2 \mathrm{d}s

        Parameters
        ----------
        H_surf
        z_surf
        r_surf

        Returns
        -------

        """
        dr_surf = np.diff(r_surf)
        dz_surf = np.diff(z_surf)

        # Pds = self.Rs * np.pi * np.trapz(abs(H_surf) ** 2 * r_surf, z_surf)
        H2r_surf = H_surf**2*r_surf
        Pds = self.Rs * np.pi * np.sum((H2r_surf[1:] + H2r_surf[:-1])/2 * np.sqrt(dr_surf**2 + dz_surf**2))

        dA_l = r_surf[:-1] * np.sqrt(dr_surf**2 + dz_surf**2)
        dA_r = r_surf[1:] * np.sqrt(dr_surf**2 + dz_surf**2)
        dA = (r_surf[1:] + r_surf[:-1])/2 * np.sqrt(dr_surf**2 + dz_surf**2)
        A_l = 2*np.pi*np.sum(dA_l)
        A_r = 2*np.pi*np.sum(dA_r)
        A = 2*np.pi*np.sum(dA)
        ic(A_l, A_r, A)

        return Pds

    def calculate_q(self):
        return (2 * np.pi * self.eig_freq * self.U) / self.Pds

    def calculate_g(self):
        """
        Calculates the geometry or form factor for the cavity

        .. math::
           G = Q_0 \cdot R_\mathrm{s} = \\frac{\omega U_n}{P_\mathrm{ds}}

        Parameters
        ----------

        Returns
        -------

        """

        return self.Q0 * self.Rs

    def get_surface_field(self, field_array, zz, rr):
        """
        Get the surface field values. Works by looping through the 2D field array and returning the first and last
        non-zero elements in each row

        Parameters
        ----------
        rr: 2D array
            r coordinates
        zz: 2D array
            z coordinates
        field_array: array
            Field values on grid

        Returns
        -------

        """

        # copy rr to avoid making changes to it
        rr_copy = np.copy(rr)
        rr_copy[field_array == 0] = 0  # this line alters rr hence the copy
        r_surf = np.max(rr_copy, axis=0)
        indx_r_surf = np.argmax(rr_copy, axis=0)
        z_surf = zz[0, :]
        surf_field = field_array[indx_r_surf, np.arange(len(z_surf))]

        return z_surf, r_surf, surf_field

    def find_resonance(self, freq):
        maara = 10  # number of eigenvalues
        raja = 1e-3  # error tolerance

        EKENTTA = 0
        HKENTTA = 1

        mu0 = 4 * np.pi * 1e-7
        e0 = 8.85418782e-12
        k0 = 2 * np.pi * freq * np.sqrt(mu0 * e0)  # "correct" eigenvalue

        ic('----------- E-walls ------------')
        ic(' Eigenvalue   error     shift')
        job = {'job': EKENTTA}
        ic("It's here here1")

        # save job job -v4
        spio.savemat(f"{self.folder}/job.mat", job, format='4')
        ic("It's here here2")
        ind, k1, u1 = self.search(k0, maara, 0, raja)
        ic("It's here here3")

        k = k1
        u = u1
        # load o_eigen
        o_eigen = spio.loadmat(f'{self.folder}/o_eigen.mat')
        index = o_eigen['index'][0]

        # save kama1 index -v4
        # save kama1 u -v4 -append
        # save kama1 k -v4 -append
        kama1 = {'index': index, 'u': u, 'k': k}
        spio.savemat(f"{self.folder}/kama1.mat", kama1, format='4')

        ic('----------- H-walls ------------')
        ic(' Eigenvalue   error     shift')
        job = {'job': HKENTTA}

        # save job job -v4
        spio.savemat(f"{self.folder}/job.mat", job, format='4')

        ind, k2, u2 = np.where(k0, maara, 0, raja)

        k = k2
        u = u2
        # load o_eigen
        o_eigen = spio.loadmat(f'{self.folder}/o_eigen.mat')
        index = o_eigen['index'][0]

        # save kama2 index -v4
        # save kama2 u -v4 -append
        # save kama2 k -v4 -append
        kama2 = {'index': index, 'u': u, 'k': k}
        spio.savemat(f"{self.folder}/kama1.mat", kama2, format='4')

        return k1, k2, u1, u2

    def calculate_vacc(self, z, E_z_axis):
        """
        Calculates the accelerating voltage of the accelerating mode

        .. math::
           V_\mathrm{cav} = \\vert \int_{z_\min}^{z_\max} E_z(r=0, z) \mathrm{e}^{jk_0 z}\\vert \mathrm{d}z

        Parameters
        ----------
        E_z_axis: array
            :math:`z` component of electric field on axis :math:`r=0`
        z: array
            :math:`z` coordinate points along axis :math:`r=0`
        Returns
        -------
        Vacc: float
            Accelerating voltage along axis :math:`r=0`

        """
        # calculate Vacc
        E_axis = E_z_axis * np.exp(1j * (2 * np.pi * self.eig_freq / c0) * z)
        Vacc = np.trapz(E_axis, z)

        return np.abs(Vacc)

    def calculate_rq(self, Vacc, U):
        """
        Calculates the :math:`R/Q_\parallel` of the accelerating mode.

        .. math::
           R/Q = \\frac{|V_{\parallel, n(0, z)}|^2}{\omega_n U_n}

        Parameters
        ----------
        Vacc: float
            Accelerating voltage.
        U: float
            Stored energy

        Returns
        -------

        """

        return Vacc ** 2 / (2 * np.pi * self.eig_freq * U)

    @staticmethod
    def calculate_epk(Epk, Eacc):
        """
        Calculates the peak surface electric field to accelerating field ratio

        Parameters
        ----------
        Epk: float
            Peak electric field
        Eacc: float
            Accelerating electric field

        Returns
        -------
        Epk/Eacc: float

        """
        return Epk / Eacc

    @staticmethod
    def calculate_bpk(Hpk, Eacc):
        """
        Calculates the peak surface magnetic field to accelerating field ratio

        Parameters
        ----------
        Hpk: float
            Peak wall surface magnetic field
        Eacc: float
            Accelerating electric field

        Returns
        -------
        Bpk/Eacc: float

        """
        return mu0 * Hpk * 1e3 / (Eacc * 1e-6)

    def calculate_loss_factor(self, R_Q):
        """
        Calculates fundamental mode loss factor

        .. math::
          k_{\parallel, n}  = \\frac{|V_{\parallel, n}(0, 0)|^2}{4 U_n} = \\frac{\omega_n}{4} \\frac{R}{Q}_n

        Parameters
        ----------
        R_Q: float
            R/Q of cavity geometry

        Returns
        -------
        k_loss: float
            Mode loss factor
        """
        return (2 * np.pi * self.eig_freq) / 4 * R_Q * 1e-12  # [V/pC]

    def calculate_volume(self, EE, zz, rr):
        """
        Calculates the volume of the volume of revolution of the cavity profile

        Parameters
        ----------
        EE: array like
            Electric field array
        zz: array like
            z coordinates on grid
        rr: array like
            r coordinates on grid

        Returns
        -------
        V: float
            Volume of revolution of cavity profile
        """

        z = zz[0, :]
        r = rr[:, 0]

        rr_copy = np.copy(rr)
        rr_copy[EE == 0] = 0

        V = 2*np.pi*np.trapz(np.trapz(rr_copy, r, axis=0), z)

        return V

    def save_qois(self):
        """
        Save the quantities of interest

        Returns
        -------

        """

        shape = {'IC': update_alpha(self.mid_cell*1e3),
                 'OC': update_alpha(self.end_cell_left*1e3),
                 'OC_R': update_alpha(self.end_cell_right*1e3)}

        with open(fr"{self.folder}\geometric_parameters.json", 'w') as f:
            json.dump(shape, f, indent=4, separators=(',', ': '))

        Req = self.mid_cell[6]
        self.G = self.Q0 * self.Rs
        self.GR_Q = self.G * self.R_Q

        # cel to cell coupling factor
        ic(self.eigen_frequencies)
        f_diff = self.eigen_frequencies[self.n_cells - 1] - self.eigen_frequencies[0]
        f_add = self.eigen_frequencies[self.n_cells - 1] + self.eigen_frequencies[0]
        self.kcc = 2 * f_diff / f_add * 100

        # field flatness
        # get max in each cell
        peaks, _ = find_peaks(abs(self.E_z_axis))
        plt.plot(self.E_z_axis_abs)
        plt.show()
        E_abs_peaks = self.E_z_axis_abs[peaks]
        # self.ff = min(E_abs_peaks) / max(E_abs_peaks) * 100
        ff = (1 - ((max(E_abs_peaks) - min(E_abs_peaks))/np.average(E_abs_peaks))) * 100

        d = {
            "Req [mm]": Req*1e3,
            "Normalization Length [mm]": self.mid_cell[5]*1e3,
            "freq [MHz]": self.eig_freq*1e-6,
            "Q []": self.Q0,
            "E [MV/m]": self.U,
            "Vacc [MV]": self.Vacc,
            "Eacc [MV/m]": self.Eacc,
            "Epk [MV/m]": self.Epk*1e-6,
            "Hpk [A/m]": self.Hpk,
            "Bpk [mT]": mu0 * self.Hpk*1e3,
            "kcc [%]": self.kcc,
            "ff [%]": self.ff,
            "Rsh [MOhm]": self.R_Q*self.Q0,
            "R/Q [Ohm]": self.R_Q,
            "Epk/Eacc []": self.epk,
            "Bpk/Eacc [mT/MV/m]": self.bpk,
            "G [Ohm]": self.G,
            "GR/Q [Ohm^2]": self.GR_Q
        }
        ic(d)
        with open(fr'{self.folder}\qois.json', "w") as f:
            json.dump(d, f, indent=4, separators=(',', ': '))

    def save_fields(self, wall, job):
        """
        Save mode fields

        Parameters
        ----------
        wall
        job

        Returns
        -------

        """
        fieldfile1 = pd.read_csv(fr"{self.folder}\fieldfile1.txt", sep='\s+',
                                 header=None).to_numpy()
        file = fieldfile1

        # compute the peak electric field on the boundary or the rf power
        # load job
        if wall == 0 and job == 1:
            E0 = self.peak_cavity_field()
        else:
            E0 = self.peak_coupler_field(wall)

        # normalize the fields
        my0 = 4e-7 * np.pi
        er = np.where(file[:, 0] == 0)
        ez = np.where(file[:, 0] == 1)
        bp = np.where(file[:, 0] == 2)
        file[er, 5:6] = file[er, 5:6] / E0[0]
        file[ez, 5:6] = file[ez, 5:6] / E0[0]
        file[bp, 5:6] = my0 * file[bp, 5:6] / E0[0]
        fieldfile1 = file

        if wall == 0:
            # save -ascii fieldfile1.n fieldfile1
            df = pd.DataFrame(fieldfile1)
            df.to_csv(fr'{self.folder}\fieldfile1.n', index=False, header=False, float_format='%.7E')
        elif wall == 1:
            # save -ascii fieldfileE.n fieldfile1
            df = pd.DataFrame(fieldfile1)
            df.to_csv(fr'{self.folder}\fieldfileE.n', index=False, header=False, float_format='%.7E')
        elif wall == 2:
            # save -ascii fieldfileH.n fieldfile1
            df = pd.DataFrame(fieldfile1)
            df.to_csv(fr'{self.folder}\fieldfileH.n', index=False, header=False, float_format='%.7E')

    def peak_cavity_field(self):

        # load mesh and field solution
        # load mesh
        # load kama0
        mesh = spio.loadmat(fr"{self.folder}\mesh.mat")
        kama0 = spio.loadmat(fr"{self.folder}\kama0.mat")

        # compute the field on the boundary
        # load geodata.n
        geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+',
                              header=None).to_numpy()

        n = len(geodata[:, 0])
        ind = np.where(geodata[3:n, 2] == 1)
        ind = ind[1:len(ind)]
        r = geodata[3:n, 0]  # r = r(ind)-5e-4
        r = r[ind] - 1.75e-4  # move points inside
        z = geodata[3:n, 1]
        z = z[ind]

        rind = np.arange(0, len(r))
        zind = np.arange(0, len(z))

        alue = 0

        # save zr z -v4
        # save zr r -append
        # save zr alue -append
        # save zr rind -append
        # save zr zind -append

        zr = {'z': z, 'r': r, 'alue': alue, 'rind': rind, 'zind': zind}
        spio.savemat(fr'{self.folder}\zr.mat', zr, format='4')

        # compute the field at generated points
        # !copy /y kama0.mat kama.mat
        # !Multipac fields
        cwd = fr'{self.folder}'
        multipacPath = fr'{self.folder}\Multipac.exe'
        if os.path.exists(multipacPath):
            subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        Er = spio.loadmat(fr"{self.folder}\Er.mat")['Er']
        Ez = spio.loadmat(fr"{self.folder}\Ez.mat")['Ez']

        ee = np.sqrt(abs(Er) ** 2 + abs(Ez) ** 2)
        E0 = max(ee)

        return E0

    def search(self, korig, maara, ind, raja):
        offset1 = -0.9
        offset2 = 0.0
        offset3 = 0.9

        laskuri = 1
        k1, _, _ = self.eee(offset1, maara, 0)
        k, u, _ = self.eee(offset2, maara, 0)
        k3, _, _ = self.eee(offset3, maara, 0)

        if ind == 0:
            ero2 = np.min(abs(k - korig))
            ind = np.argmin(abs(k - korig))
            if abs(k[ind]) < 1e-6:
                ind = ind + 1

        ero1 = korig - k1[ind]
        s1 = np.sign(ero1)
        ero2 = korig - k[ind]
        s2 = np.sign(ero2)
        ero3 = korig - k3[ind]
        s3 = np.sign(ero3)
        # ic([s1 s2 s3])

        if (s1 == 1) and (s2 == 1) and (s3 == 1):
            ic('?????????????????????????????')
            ic([k1[ind], k[ind], k3[ind]])
            ic("Can't compress enough")

        if (s1 == -1) and (s2 == -1) and (s3 == -1):
            jatka = 1
            while (jatka > 0) & (jatka < 4):
                ic('problem, have not stretched enough, trying to fix')
                offset1 = offset1 - 0.9

                k3 = k
                offset3 = offset2
                s3 = s2
                k = k1
                offset2 = offset1
                s2 = s1
                k1 = self.eee(offset1, maara, 0)
                ero1 = korig - k1[ind]
                s1 = np.sign(ero3)
                if s3 == 1:
                    jatka = 0
                    ic('succeeded')
                else:
                    jatka = jatka + 1

        if jatka != 0:
            ic('FAILED')

        while abs(ero2) > raja:
            if s1 != s2:
                s3 = s2
                offset3 = offset2
                offset2 = 0.5 * (offset1 + offset3)
                s3 = s2
            elif s2 != s3:
                s1 = s2
                offset1 = offset2
                offset2 = 0.5 * (offset2 + offset3)
                s1 = s2

            k, u, siirto = self.eee(offset2, maara, 0)
            ero2 = korig - k[ind]
            s2 = np.sign(ero2)
            #   #  ic([k(ind), ero2, offset2, siirto])
            ic([k[ind], ero2, siirto])
        offset1 = -0.9
        offset2 = 0
        offset3 = 0.9

        laskuri = 1
        k1, _, _ = self.eee(offset1, maara, 0)
        k, u, siirto = self.eee(offset2, maara, 0)
        k3, _, _ = self.eee(offset3, maara, 0)

        if ind == 0:
            ero2 = np.min(abs(k - korig))
            ind = np.argmin(abs(k - korig))
            if abs(k[ind]) < 1e-6:
                ind = ind + 1

        # ero1 = korig-k1(ind) s1 = np.sign(ero1)
        # ero2 = korig-k(ind) s2 = np.sign(ero2)
        # ero3 = korig-k3(ind) s3 = np.sign(ero3)
        # #ic([s1 s2 s3])

        if (s1 == 1) & (s2 == 1) & (s3 == 1):
            ic('?????????????????????????????')
            ic([k1[ind], k[ind], k3[ind]])
            ic('Can''t compress enough')

        if (s1 == -1) & (s2 == -1) & (s3 == -1):
            jatka = 1
            while (jatka > 0) & (jatka < 4):
                ic('problem, have not stretched enough, trying to fix')
                offset1 = offset1 - 0.9

                k3 = k
                offset3 = offset2
                s3 = s2
                k = k1
                offset2 = offset1
                s2 = s1
                k1 = self.eee(offset1, maara, 0)
                ero1 = korig - k1[ind]
                s1 = np.sign(ero3)
                if s3 == 1:
                    jatka = 0
                    ic('succeeded')
                else:
                    jatka = jatka + 1

            if jatka != 0:
                ic('FAILED')

        while abs(ero2) > raja:
            if s1 != s2:
                s3 = s2
                offset3 = offset2
                offset2 = 0.5 * (offset1 + offset3)
                s3 = s2
            elif s2 != s3:
                s1 = s2
                offset1 = offset2
                offset2 = 0.5 * (offset2 + offset3)
                s1 = s2

            k, u, siirto = self.eee(offset2, maara, 0)
            ero2 = korig - k[ind]
            s2 = np.sign(ero2)
            #  ic([k(ind), ero2, offset2, siirto])
            ic([k[ind], ero2, siirto])

        k = k[ind]
        u = u[:, ind]
        # # -----------------------------------------------------------------------

        return ind, k, u

    def eee(self, offset, maara, show):
        offset = {'offset': offset}
        # save offset offset -v4
        spio.savemat(f"{self.folder}/offset.mat", offset, format='4')
        ic(offset)
        # options.disp = 0

        # !eigenC_bin

        ic('self.eee: started eigenmode analysis')
        eigenCpath = fr'{self.folder}\eigenC_bin.exe'
        if os.path.exists(eigenCpath):
            subprocess.call(eigenCpath, cwd=self.folder,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        ic('self.eee: finised eigenmode analysis')

        # load o_eigen
        # load cinfo
        o_eigen = spio.loadmat(f'{self.folder}/o_eigen.mat')
        cinfo = spio.loadmat(f'{self.folder}/cinfo.mat')['cinfo'][0][0]
        ic(cinfo)

        n = int(o_eigen['n'])
        ia = o_eigen['ia'].T[0] - 1  # mnatlab indices to python
        ja = o_eigen['ja'].T[0] - 1
        Aarvot = o_eigen['Aarvot'].T[0]
        Barvot = o_eigen['Barvot'].T[0]
        nf = o_eigen['nf'][0]
        index = o_eigen['index'][0]
        siirto = o_eigen['siirto'][0]

        if cinfo < 0:
            if cinfo == -1:
                ic('ran out of memory')
            elif cinfo == -2:
                ic('unknown type element found in mesh')
            elif (cinfo == -3):
                ic('no compressible area found in mesh')
            elif cinfo == -4:
                ic('could not compress mesh that much')
            elif cinfo == -5:
                ic('could not add element in matrix')
            elif cinfo == -6:
                ic('nodes in mesh were in wrong order')
            elif cinfo == -7:
                ic('unknown job')
            else:
                ic('unknown error')

        AA = sps.csr_matrix((Aarvot, (ia, ja)), shape=(n, n))
        BB = sps.csr_matrix((Barvot, (ia, ja)), shape=(n, n))
        # ic([offset, siirto])

        d2, u = spsl.eigs(AA, M=BB, k=maara, sigma=0)
        d2 = np.absolute(d2)

        k = np.sqrt(d2)

        if show > 0:
            ic(k[1:show])

        # u = u2
        # for i = 1:n
        #   u(i,:) = u2(nf(i),:)
        # end
        # # --------------------------------------------------------------------
        return k, u, siirto

    def normalize_u(self, wall=0):
        ic("Inside normalize_u")
        job = 0
        if wall == 0:
            E0 = self.peak_cavity_field()
            # normalize u
            kama0 = spio.loadmat(fr"{self.folder}\kama0.mat")
            index = kama0['index']
            u = kama0['u']
            k = kama0['k']
            u = u / E0

            kama_n = {'index': index, 'u': u, 'k': k}
            spio.savemat(fr'{self.folder}\kama_n.mat', kama_n, format='4')

    def run(self, n_cells, mid_cell, end_cell_left=None, end_cell_right=None, beampipe=None,
            req_mode_num=None, plot=False, gridcons=0.005):
        self.define_geometry(n_cells, mid_cell, end_cell_left, end_cell_right, beampipe=beampipe, plot=plot)
        self.generate_mesh(plot=True, gridcons=gridcons)
        self.run_field_solver(req_mode_num=req_mode_num, show_plots=True)
        ic("Analysis ended.")

    @staticmethod
    def inttri(p):
        # function [X,W] = inttri(p)
        #
        # The program gives the sample points X as a (2,n)-matrix and
        # the weights W as an (n,1)-matrix for an integration over a
        # planar triangle with the vertices (0,0), (1,0) and (0,1).
        # The input argument p indicates the maximal degree of the poly-
        # nomials of two variables which will be integrated exactly
        # correctly, provided p<20. By a linear transformation, X and
        # W are obtained from the paper: D.A. Dunavant: High degree
        # efficient symmetrical gaussian quadrature rules for the
        # triangle.- International Journal for Num. Methods in Eng.,
        # vol. 21, 1129-1148(1985).
        # ----------------------------------------------------------
        # CALLS TO none
        # 20 JUNI 1993:  Jukka Sarvas, Rolf Nevanlinna Institute
        # ----------------------------------------------------------

        if p == 11:
            p = 12
        elif (p == 15) | (p == 16):
            p = 17
        elif (p == 18) | (p > 19):
            p = 19

        table = {
            '1': [3.33333333333333e-01, 3.33333333333333e-01, 5.00000000000000e-01],
            '2': [1.66666666666666e-01, 6.66666666666667e-01, 1.66666666666667e-01, 6.66666666666667e-01,
                  1.66666666666667e-01, 1.66666666666666e-01, 1.66666666666667e-01, 1.66666666666667e-01,
                  1.66666666666667e-01],
            '3': [3.33333333333333e-0, 2.00000000000000e-01, 6.00000000000000e-01, 2.00000000000000e-01,
                  3.33333333333333e-01, 6.00000000000000e-01, 2.00000000000000e-01, 2.00000000000000e-01,
                  -2.81250000000000e-01, 2.60416666666667e-01, 2.60416666666667e-01, 2.60416666666667e-01],
            '4': [4.45948490915965e-01, 9.15762135097699e-02, 1.08103018168070e-01, 8.16847572980459e-01,
                  4.45948490915965e-01, 9.15762135097710e-02, 1.08103018168070e-01, 8.16847572980459e-01,
                  4.45948490915965e-01, 9.15762135097710e-02, 4.45948490915965e-01, 9.15762135097699e-02,
                  1.11690794839006e-01, 5.49758718276610e-02, 1.11690794839006e-01, 5.49758718276610e-02,
                  1.11690794839006e-01, 5.49758718276610e-02],
            '5': [3.33333333333333e-01, 4.70142064105115e-01, 1.01286507323457e-01, 5.97158717897700e-02,
                  7.97426985353087e-01, 4.70142064105115e-01, 1.01286507323456e-01, 3.33333333333333e-01,
                  5.97158717897700e-02, 7.97426985353087e-01, 4.70142064105115e-01, 1.01286507323456e-01,
                  4.70142064105115e-01, 1.01286507323457e-01, 1.12500000000000e-01, 6.61970763942530e-02,
                  6.29695902724135e-02, 6.61970763942530e-02, 6.29695902724135e-02, 6.61970763942530e-02,
                  6.29695902724135e-02],
            '6': [2.49286745170911e-01, 6.30890144915021e-02, 5.01426509658179e-01, 8.73821971016996e-01,
                  2.49286745170910e-01, 6.30890144915020e-02, 6.36502499121399e-01, 3.10352451033784e-01,
                  5.31450498448170e-02, 5.31450498448170e-02, 3.10352451033784e-01, 6.36502499121399e-01,
                  5.01426509658179e-01, 8.73821971016996e-01, 2.49286745170910e-01, 6.30890144915020e-02,
                  2.49286745170911e-01, 6.30890144915021e-02, 5.31450498448170e-02, 5.31450498448170e-02,
                  3.10352451033784e-01, 6.36502499121399e-01, 6.36502499121399e-01, 3.10352451033784e-01,
                  5.83931378631895e-02, 2.54224531851035e-02, 5.83931378631895e-02, 2.54224531851035e-02,
                  5.83931378631895e-02, 2.54224531851035e-02, 4.14255378091870e-02, 4.14255378091870e-02,
                  4.14255378091870e-02, 4.14255378091870e-02, 4.14255378091870e-02, 4.14255378091870e-02],
            '7': [3.33333333333333e-01, 2.60345966079040e-01, 6.51301029022160e-02, 4.79308067841920e-01,
                  8.69739794195568e-01, 2.60345966079040e-01, 6.51301029022160e-02, 6.38444188569810e-01,
                  3.12865496004874e-01, 4.86903154253160e-02, 4.86903154253160e-02, 3.12865496004874e-01,
                  6.38444188569810e-01, 3.33333333333333e-01, 4.79308067841920e-01, 8.69739794195568e-01,
                  2.60345966079040e-01, 6.51301029022160e-02, 2.60345966079040e-01, 6.51301029022160e-02,
                  4.86903154253160e-02, 4.86903154253160e-02, 3.12865496004874e-01, 6.38444188569810e-01,
                  6.38444188569810e-01, 3.12865496004874e-01, -7.47850222338410e-02, 8.78076287166040e-02,
                  2.66736178044190e-02, 8.78076287166040e-02, 2.66736178044190e-02, 8.78076287166040e-02,
                  2.66736178044190e-02, 3.85568804451285e-02, 3.85568804451285e-02, 3.85568804451285e-02,
                  3.85568804451285e-02, 3.85568804451285e-02, 3.85568804451285e-02],
            '8': [3.33333333333333e-01, 4.59292588292723e-01, 1.70569307751760e-01, 5.05472283170311e-02,
                  8.14148234145540e-02, 6.58861384496480e-01, 8.98905543365938e-01, 4.59292588292723e-01,
                  1.70569307751760e-01, 5.05472283170310e-02, 7.28492392955404e-01, 2.63112829634638e-01,
                  8.39477740995800e-03, 8.39477740995800e-03, 2.63112829634638e-01, 7.28492392955404e-01,
                  3.33333333333333e-01, 8.14148234145540e-02, 6.58861384496480e-01, 8.98905543365938e-01,
                  4.59292588292723e-01, 1.70569307751760e-01, 5.05472283170310e-02, 4.59292588292723e-01,
                  1.70569307751760e-01, 5.05472283170311e-02, 8.39477740995800e-03, 8.39477740995800e-03,
                  2.63112829634638e-01, 7.28492392955404e-01, 7.28492392955404e-01, 2.63112829634638e-01,
                  7.21578038388935e-02, 4.75458171336425e-02, 5.16086852673590e-02, 1.62292488115990e-02,
                  4.75458171336425e-02, 5.16086852673590e-02, 1.62292488115990e-02, 4.75458171336425e-02,
                  5.16086852673590e-02, 1.62292488115990e-02, 1.36151570872175e-02, 1.36151570872175e-02,
                  1.36151570872175e-02, 1.36151570872175e-02, 1.36151570872175e-02, 1.36151570872175e-02],
            '9': [3.33333333333333e-01, 4.89682519198737e-01, 4.37089591492936e-01, 1.88203535619032e-01,
                  4.47295133944519e-02, 2.06349616025250e-02, 1.25820817014127e-01, 6.23592928761935e-01,
                  9.10540973211095e-01, 4.89682519198738e-01, 4.37089591492937e-01, 1.88203535619033e-01,
                  4.47295133944530e-02, 7.41198598784498e-01, 2.21962989160766e-01, 3.68384120547360e-02,
                  3.68384120547360e-02, 2.21962989160766e-01, 7.41198598784498e-01, 3.33333333333333e-01,
                  2.06349616025250e-02, 1.25820817014127e-01, 6.23592928761935e-01, 9.10540973211095e-01,
                  4.89682519198738e-01, 4.37089591492937e-01, 1.88203535619033e-01, 4.47295133944530e-02,
                  4.89682519198737e-01, 4.37089591492936e-01, 1.88203535619032e-01, 4.47295133944519e-02,
                  3.68384120547360e-02, 3.68384120547360e-02, 2.21962989160766e-01, 7.41198598784498e-01,
                  7.41198598784498e-01, 2.21962989160766e-01, 4.85678981413995e-02, 1.56673501135695e-02,
                  3.89137705023870e-02, 3.98238694636050e-02, 1.27888378293490e-02, 1.56673501135695e-02,
                  3.89137705023870e-02, 3.98238694636050e-02, 1.27888378293490e-02, 1.56673501135695e-02,
                  3.89137705023870e-02, 3.98238694636050e-02, 1.27888378293490e-02, 2.16417696886445e-02,
                  2.16417696886445e-02, 2.16417696886445e-02, 2.16417696886445e-02, 2.16417696886445e-02,
                  2.16417696886445e-02],
            '10': [3.33333333333333e-01, 4.85577633383658e-01, 1.09481575485037e-01, 2.88447332326850e-02,
                   7.81036849029926e-01, 4.85577633383657e-01, 1.09481575485037e-01, 5.50352941820999e-01,
                   7.28323904597411e-01, 9.23655933587501e-01, 3.07939838764121e-01, 2.46672560639903e-01,
                   6.68032510122000e-02, 1.41707219414880e-01, 2.50035347626860e-02, 9.54081540029900e-03,
                   1.41707219414880e-01, 2.50035347626860e-02, 9.54081540029900e-03, 3.07939838764121e-01,
                   2.46672560639903e-01, 6.68032510122000e-02, 5.50352941820999e-01, 7.28323904597411e-01,
                   9.23655933587501e-01, 3.33333333333333e-01, 2.88447332326850e-02, 7.81036849029926e-01,
                   4.85577633383657e-01, 1.09481575485037e-01, 4.85577633383658e-01, 1.09481575485037e-01,
                   1.41707219414880e-01, 2.50035347626860e-02, 9.54081540029900e-03, 1.41707219414880e-01,
                   2.50035347626860e-02, 9.54081540029900e-03, 3.07939838764121e-01, 2.46672560639903e-01,
                   6.68032510122000e-02, 5.50352941820999e-01, 7.28323904597411e-01, 9.23655933587501e-01,
                   5.50352941820999e-01, 7.28323904597411e-01, 9.23655933587501e-01, 3.07939838764121e-01,
                   2.46672560639903e-01, 6.68032510122000e-02, 4.54089951913770e-02, 1.83629788782335e-02,
                   2.26605297177640e-02, 1.83629788782335e-02, 2.26605297177640e-02, 1.83629788782335e-02,
                   2.26605297177640e-02, 3.63789584227100e-02, 1.41636212655285e-02, 4.71083348186650e-03,
                   3.63789584227100e-02, 1.41636212655285e-02, 4.71083348186650e-03, 3.63789584227100e-02,
                   1.41636212655285e-02, 4.71083348186650e-03, 3.63789584227100e-02, 1.41636212655285e-02,
                   4.71083348186650e-03, 3.63789584227100e-02, 1.41636212655285e-02, 4.71083348186650e-03,
                   3.63789584227100e-02, 1.41636212655285e-02, 4.71083348186650e-03],
            '12': [4.88217389773805e-01, 4.39724392294461e-01, 2.71210385012116e-01, 1.27576145541586e-01,
                   2.13173504532110e-02, 2.35652204523900e-02, 1.20551215411079e-01, 4.57579229975768e-01,
                   7.44847708916828e-01, 9.57365299093579e-01, 4.88217389773805e-01, 4.39724392294460e-01,
                   2.71210385012116e-01, 1.27576145541586e-01, 2.13173504532100e-02, 6.08943235779788e-01,
                   6.95836086787803e-01, 8.58014033544073e-01, 2.75713269685514e-01, 2.81325580989940e-01,
                   1.16251915907597e-01, 1.15343494534698e-01, 2.28383322222570e-02, 2.57340505483300e-02,
                   1.15343494534698e-01, 2.28383322222570e-02, 2.57340505483300e-02, 2.75713269685514e-01,
                   2.81325580989940e-01, 1.16251915907597e-01, 6.08943235779788e-01, 6.95836086787803e-01,
                   8.58014033544073e-01, 2.35652204523900e-02, 1.20551215411079e-01, 4.57579229975768e-01,
                   7.44847708916828e-01, 9.57365299093579e-01, 4.88217389773805e-01, 4.39724392294460e-01,
                   2.71210385012116e-01, 1.27576145541586e-01, 2.13173504532100e-02, 4.88217389773805e-01,
                   4.39724392294461e-01, 2.71210385012116e-01, 1.27576145541586e-01, 2.13173504532110e-02,
                   1.15343494534698e-01, 2.28383322222570e-02, 2.57340505483300e-02, 1.15343494534698e-01,
                   2.28383322222570e-02, 2.57340505483300e-02, 2.75713269685514e-01, 2.81325580989940e-01,
                   1.16251915907597e-01, 6.08943235779788e-01, 6.95836086787803e-01, 8.58014033544073e-01,
                   6.08943235779788e-01, 6.95836086787803e-01, 8.58014033544073e-01, 2.75713269685514e-01,
                   2.81325580989940e-01, 1.16251915907597e-01, 1.28655332202275e-02, 2.18462722690190e-02,
                   3.14291121089425e-02, 1.73980564653545e-02, 3.08313052577950e-03, 1.28655332202275e-02,
                   2.18462722690190e-02, 3.14291121089425e-02, 1.73980564653545e-02, 3.08313052577950e-03,
                   1.28655332202275e-02, 2.18462722690190e-02, 3.14291121089425e-02, 1.73980564653545e-02,
                   3.08313052577950e-03, 2.01857788831905e-02, 1.11783866011515e-02, 8.65811555432950e-03,
                   2.01857788831905e-02, 1.11783866011515e-02, 8.65811555432950e-03, 2.01857788831905e-02,
                   1.11783866011515e-02, 8.65811555432950e-03, 2.01857788831905e-02, 1.11783866011515e-02,
                   8.65811555432950e-03, 2.01857788831905e-02, 1.11783866011515e-02, 8.65811555432950e-03,
                   2.01857788831905e-02, 1.11783866011515e-02, 8.65811555432950e-03],
            '13': [3.33333333333333e-01, 4.95048184939704e-01, 4.68716635109574e-01, 4.14521336801276e-01,
                   2.29399572042832e-01, 1.14424495196330e-01, 2.48113913634590e-02, 9.90363012059100e-03,
                   6.25667297808520e-02, 1.70957326397447e-01, 5.41200855914337e-01, 7.71151009607340e-01,
                   9.50377217273082e-01, 4.95048184939705e-01, 4.68716635109574e-01, 4.14521336801277e-01,
                   2.29399572042831e-01, 1.14424495196330e-01, 2.48113913634590e-02, 6.36351174561660e-01,
                   6.90169159986905e-01, 8.51409537834241e-01, 2.68794997058761e-01, 2.91730066734288e-01,
                   1.26357385491669e-01, 9.48538283795790e-02, 1.81007732788070e-02, 2.22330766740900e-02,
                   9.48538283795790e-02, 1.81007732788070e-02, 2.22330766740900e-02, 2.68794997058761e-01,
                   2.91730066734288e-01, 1.26357385491669e-01, 6.36351174561660e-01, 6.90169159986905e-01,
                   8.51409537834241e-01, 3.33333333333333e-01, 9.90363012059100e-03, 6.25667297808520e-02,
                   1.70957326397447e-01, 5.41200855914337e-01, 7.71151009607340e-01, 9.50377217273082e-01,
                   4.95048184939705e-01, 4.68716635109574e-01, 4.14521336801277e-01, 2.29399572042831e-01,
                   1.14424495196330e-01, 2.48113913634590e-02, 4.95048184939704e-01, 4.68716635109574e-01,
                   4.14521336801276e-01, 2.29399572042832e-01, 1.14424495196330e-01, 2.48113913634590e-02,
                   9.48538283795790e-02, 1.81007732788070e-02, 2.22330766740900e-02, 9.48538283795790e-02,
                   1.81007732788070e-02, 2.22330766740900e-02, 2.68794997058761e-01, 2.91730066734288e-01,
                   1.26357385491669e-01, 6.36351174561660e-01, 6.90169159986905e-01, 8.51409537834241e-01,
                   6.36351174561660e-01, 6.90169159986905e-01, 8.51409537834241e-01, 2.68794997058761e-01,
                   2.91730066734288e-01, 1.26357385491669e-01, 2.62604617004010e-02, 5.64007260466500e-03,
                   1.57117591812270e-02, 2.35362512520970e-02, 2.36817932681775e-02, 1.55837645228970e-02,
                   3.98788573253700e-03, 5.64007260466500e-03, 1.57117591812270e-02, 2.35362512520970e-02,
                   2.36817932681775e-02, 1.55837645228970e-02, 3.98788573253700e-03, 5.64007260466500e-03,
                   1.57117591812270e-02, 2.35362512520970e-02, 2.36817932681775e-02, 1.55837645228970e-02,
                   3.98788573253700e-03, 1.84242013643660e-02, 8.70073165191100e-03, 7.76089341952250e-03,
                   1.84242013643660e-02, 8.70073165191100e-03, 7.76089341952250e-03, 1.84242013643660e-02,
                   8.70073165191100e-03, 7.76089341952250e-03, 1.84242013643660e-02, 8.70073165191100e-03,
                   7.76089341952250e-03, 1.84242013643660e-02, 8.70073165191100e-03, 7.76089341952250e-03,
                   1.84242013643660e-02, 8.70073165191100e-03, 7.76089341952250e-03],
            "14": [4.8896391036217808e-001, 4.1764471934045408e-001, 2.7347752830883808e-001, 1.7720553241254400e-001,
                   6.1799883090871952e-002, 1.9390961248700980e-002, 2.2072179275642996e-002, 1.6471056131909198e-001,
                   4.5304494338232296e-001, 6.4558893517491304e-001, 8.7640023381825504e-001, 9.6121807750259808e-001,
                   4.8896391036217896e-001, 4.1764471934045392e-001, 2.7347752830883896e-001, 1.7720553241254300e-001,
                   6.1799883090873000e-002, 1.9390961248701000e-002, 7.7060855477499600e-001, 5.7022229084668320e-001,
                   6.8698016780808800e-001, 8.7975717137017088e-001, 1.7226668782135598e-001, 3.3686145979634492e-001,
                   2.9837288213625800e-001, 1.1897449769695700e-001, 5.7124757403647992e-002, 9.2916249356971984e-002,
                   1.4646950055654000e-002, 1.2683309328720000e-003, 5.7124757403647992e-002, 9.2916249356971984e-002,
                   1.4646950055654000e-002, 1.2683309328720000e-003, 1.7226668782135598e-001, 3.3686145979634492e-001,
                   2.9837288213625800e-001, 1.1897449769695700e-001, 7.7060855477499600e-001, 5.7022229084668320e-001,
                   6.8698016780808800e-001, 8.7975717137017088e-001, 2.2072179275642996e-002, 1.6471056131909198e-001,
                   4.5304494338232296e-001, 6.4558893517491304e-001, 8.7640023381825504e-001, 9.6121807750259808e-001,
                   4.8896391036217896e-001, 4.1764471934045392e-001, 2.7347752830883896e-001, 1.7720553241254300e-001,
                   6.1799883090873000e-002, 1.9390961248701000e-002, 4.8896391036217808e-001, 4.1764471934045408e-001,
                   2.7347752830883808e-001, 1.7720553241254400e-001, 6.1799883090871952e-002, 1.9390961248700980e-002,
                   5.7124757403647992e-002, 9.2916249356971984e-002, 1.4646950055654000e-002, 1.2683309328720000e-003,
                   5.7124757403647992e-002, 9.2916249356971984e-002, 1.4646950055654000e-002, 1.2683309328720000e-003,
                   1.7226668782135598e-001, 3.3686145979634492e-001, 2.9837288213625800e-001, 1.1897449769695700e-001,
                   7.7060855477499600e-001, 5.7022229084668320e-001, 6.8698016780808800e-001, 8.7975717137017088e-001,
                   7.7060855477499600e-001, 5.7022229084668320e-001, 6.8698016780808800e-001, 8.7975717137017088e-001,
                   1.7226668782135598e-001, 3.3686145979634492e-001, 2.9837288213625800e-001, 1.1897449769695700e-001,
                   1.0941790684714502e-002, 1.6394176772062500e-002, 2.5887052253645996e-002, 2.1081294368496500e-002,
                   7.2168498348885008e-003, 2.4617018011999996e-003, 1.0941790684714502e-002, 1.6394176772062500e-002,
                   2.5887052253645996e-002, 2.1081294368496500e-002, 7.2168498348885008e-003, 2.4617018011999996e-003,
                   1.0941790684714502e-002, 1.6394176772062500e-002, 2.5887052253645996e-002, 2.1081294368496500e-002,
                   7.2168498348885008e-003, 2.4617018011999996e-003, 1.2332876606282000e-002, 1.9285755393530500e-002,
                   7.2181540567669984e-003, 2.5051144192504996e-003, 1.2332876606282000e-002, 1.9285755393530500e-002,
                   7.2181540567669984e-003, 2.5051144192504996e-003, 1.2332876606282000e-002, 1.9285755393530500e-002,
                   7.2181540567669984e-003, 2.5051144192504996e-003, 1.2332876606282000e-002, 1.9285755393530500e-002,
                   7.2181540567669984e-003, 2.5051144192504996e-003, 1.2332876606282000e-002, 1.9285755393530500e-002,
                   7.2181540567669984e-003, 2.5051144192504996e-003, 1.2332876606282000e-002, 1.9285755393530500e-002,
                   7.2181540567669984e-003, 2.5051144192504996e-003],
            '17': [3.3333333333333332e-001, 4.9717054055677400e-001, 4.8217632262462408e-001, 4.5023996902078096e-001,
                   4.0026623937739688e-001, 2.5214126797095200e-001, 1.6204700465846200e-001, 7.5875882260746080e-002,
                   1.5654726967821994e-002, 5.6589188864519992e-003, 3.5647354750751004e-002, 9.9520061958437008e-002,
                   1.9946752124520604e-001, 4.9571746405809496e-001, 6.7590599068307696e-001, 8.4824823547850784e-001,
                   9.6869054606435600e-001, 4.9717054055677400e-001, 4.8217632262462496e-001, 4.5023996902078208e-001,
                   4.0026623937739704e-001, 2.5214126797095300e-001, 1.6204700465846100e-001, 7.5875882260746000e-002,
                   1.5654726967822000e-002, 6.5549320380942296e-001, 5.7233759053202008e-001, 6.2600119028622704e-001,
                   7.9642721497407104e-001, 7.5235100593773008e-001, 9.0462550409560800e-001, 3.3431986736365804e-001,
                   2.9222153779694400e-001, 3.1957488542319000e-001, 1.9070422419229200e-001, 1.8048321164874596e-001,
                   8.0711313679564016e-002, 1.0186928826919000e-002, 1.3544087167103600e-001, 5.4423924290583000e-002,
                   1.2868560833637000e-002, 6.7165782413524000e-002, 1.4663182224828000e-002, 1.0186928826919000e-002,
                   1.3544087167103600e-001, 5.4423924290583000e-002, 1.2868560833637000e-002, 6.7165782413524000e-002,
                   1.4663182224828000e-002, 3.3431986736365804e-001, 2.9222153779694400e-001, 3.1957488542319000e-001,
                   1.9070422419229200e-001, 1.8048321164874596e-001, 8.0711313679564016e-002, 6.5549320380942296e-001,
                   5.7233759053202008e-001, 6.2600119028622704e-001, 7.9642721497407104e-001, 7.5235100593773008e-001,
                   9.0462550409560800e-001, 3.3333333333333332e-001, 5.6589188864519992e-003, 3.5647354750751004e-002,
                   9.9520061958437008e-002, 1.9946752124520604e-001, 4.9571746405809496e-001, 6.7590599068307696e-001,
                   8.4824823547850784e-001, 9.6869054606435600e-001, 4.9717054055677400e-001, 4.8217632262462496e-001,
                   4.5023996902078208e-001, 4.0026623937739704e-001, 2.5214126797095300e-001, 1.6204700465846100e-001,
                   7.5875882260746000e-002, 1.5654726967822000e-002, 4.9717054055677400e-001, 4.8217632262462408e-001,
                   4.5023996902078096e-001, 4.0026623937739688e-001, 2.5214126797095200e-001, 1.6204700465846200e-001,
                   7.5875882260746080e-002, 1.5654726967821994e-002, 1.0186928826919000e-002, 1.3544087167103600e-001,
                   5.4423924290583000e-002, 1.2868560833637000e-002, 6.7165782413524000e-002, 1.4663182224828000e-002,
                   1.0186928826919000e-002, 1.3544087167103600e-001, 5.4423924290583000e-002, 1.2868560833637000e-002,
                   6.7165782413524000e-002, 1.4663182224828000e-002, 3.3431986736365804e-001, 2.9222153779694400e-001,
                   3.1957488542319000e-001, 1.9070422419229200e-001, 1.8048321164874596e-001, 8.0711313679564016e-002,
                   6.5549320380942296e-001, 5.7233759053202008e-001, 6.2600119028622704e-001, 7.9642721497407104e-001,
                   7.5235100593773008e-001, 9.0462550409560800e-001, 6.5549320380942296e-001, 5.7233759053202008e-001,
                   6.2600119028622704e-001, 7.9642721497407104e-001, 7.5235100593773008e-001, 9.0462550409560800e-001,
                   3.3431986736365804e-001, 2.9222153779694400e-001, 3.1957488542319000e-001, 1.9070422419229200e-001,
                   1.8048321164874596e-001, 8.0711313679564016e-002, 1.6718599645401496e-002, 2.5467077202534996e-003,
                   7.3354322638190000e-003, 1.2175439176836000e-002, 1.5553775434484498e-002, 1.5628555609310000e-002,
                   1.2407827169832500e-002, 7.0280365352784992e-003, 1.5973380868895002e-003, 2.5467077202534996e-003,
                   7.3354322638190000e-003, 1.2175439176836000e-002, 1.5553775434484498e-002, 1.5628555609310000e-002,
                   1.2407827169832500e-002, 7.0280365352784992e-003, 1.5973380868895002e-003, 2.5467077202534996e-003,
                   7.3354322638190000e-003, 1.2175439176836000e-002, 1.5553775434484498e-002, 1.5628555609310000e-002,
                   1.2407827169832500e-002, 7.0280365352784992e-003, 1.5973380868895002e-003, 4.0598276594965000e-003,
                   1.3402871141581498e-002, 9.2299966054110000e-003, 4.2384342671640000e-003, 9.1463983850125008e-003,
                   3.3328160020824996e-003, 4.0598276594965000e-003, 1.3402871141581498e-002, 9.2299966054110000e-003,
                   4.2384342671640000e-003, 9.1463983850125008e-003, 3.3328160020824996e-003, 4.0598276594965000e-003,
                   1.3402871141581498e-002, 9.2299966054110000e-003, 4.2384342671640000e-003, 9.1463983850125008e-003,
                   3.3328160020824996e-003, 4.0598276594965000e-003, 1.3402871141581498e-002, 9.2299966054110000e-003,
                   4.2384342671640000e-003, 9.1463983850125008e-003, 3.3328160020824996e-003, 4.0598276594965000e-003,
                   1.3402871141581498e-002, 9.2299966054110000e-003, 4.2384342671640000e-003, 9.1463983850125008e-003,
                   3.3328160020824996e-003, 4.0598276594965000e-003, 1.3402871141581498e-002, 9.2299966054110000e-003,
                   4.2384342671640000e-003, 9.1463983850125008e-003, 3.3328160020824996e-003],
            '19': [3.3333333333333332e-001, 4.8960998707300704e-001, 4.5453689269789208e-001, 4.0141668064943096e-001,
                   2.5555165440309700e-001, 1.7707794215212902e-001, 1.1006105322795214e-001, 5.5528624251839048e-002,
                   1.2621863777228030e-002, 2.0780025853987000e-002, 9.0926214604214992e-002, 1.9716663870113800e-001,
                   4.8889669119380496e-001, 6.4584411569574096e-001, 7.7987789354409584e-001, 8.8894275149632096e-001,
                   9.7475627244554304e-001, 4.8960998707300600e-001, 4.5453689269789296e-001, 4.0141668064943104e-001,
                   2.5555165440309800e-001, 1.7707794215213002e-001, 1.1006105322795200e-001, 5.5528624251839992e-002,
                   1.2621863777228998e-002, 6.0063379479464496e-001, 5.5760326158878384e-001, 7.2098702581736496e-001,
                   5.9452706895587104e-001, 8.3933147368083808e-001, 7.0108797892617304e-001, 8.2293132406985696e-001,
                   9.2434425262078384e-001, 3.9575478735694296e-001, 3.0792998388043608e-001, 2.6456694840652000e-001,
                   3.5853935220595100e-001, 1.5780740596859500e-001, 7.5050596975911008e-002, 1.4242160111338298e-001,
                   6.5494628082938000e-002, 3.6114178484120000e-003, 1.3446675453078000e-001, 1.4446025776115000e-002,
                   4.6933578838178000e-002, 2.8611203505670000e-003, 2.2386142409791600e-001, 3.4647074816760000e-002,
                   1.0161119296278000e-002, 3.6114178484120000e-003, 1.3446675453078000e-001, 1.4446025776115000e-002,
                   4.6933578838178000e-002, 2.8611203505670000e-003, 2.2386142409791600e-001, 3.4647074816760000e-002,
                   1.0161119296278000e-002, 3.9575478735694296e-001, 3.0792998388043608e-001, 2.6456694840652000e-001,
                   3.5853935220595100e-001, 1.5780740596859500e-001, 7.5050596975911008e-002, 1.4242160111338298e-001,
                   6.5494628082938000e-002, 6.0063379479464496e-001, 5.5760326158878384e-001, 7.2098702581736496e-001,
                   5.9452706895587104e-001, 8.3933147368083808e-001, 7.0108797892617304e-001, 8.2293132406985696e-001,
                   9.2434425262078384e-001, 3.3333333333333332e-001, 2.0780025853987000e-002, 9.0926214604214992e-002,
                   1.9716663870113800e-001, 4.8889669119380496e-001, 6.4584411569574096e-001, 7.7987789354409584e-001,
                   8.8894275149632096e-001, 9.7475627244554304e-001, 4.8960998707300600e-001, 4.5453689269789296e-001,
                   4.0141668064943104e-001, 2.5555165440309800e-001, 1.7707794215213002e-001, 1.1006105322795200e-001,
                   5.5528624251839992e-002, 1.2621863777228998e-002, 4.8960998707300704e-001, 4.5453689269789208e-001,
                   4.0141668064943096e-001, 2.5555165440309700e-001, 1.7707794215212902e-001, 1.1006105322795214e-001,
                   5.5528624251839048e-002, 1.2621863777228030e-002, 3.6114178484120000e-003, 1.3446675453078000e-001,
                   1.4446025776115000e-002, 4.6933578838178000e-002, 2.8611203505670000e-003, 2.2386142409791600e-001,
                   3.4647074816760000e-002, 1.0161119296278000e-002, 3.6114178484120000e-003, 1.3446675453078000e-001,
                   1.4446025776115000e-002, 4.6933578838178000e-002, 2.8611203505670000e-003, 2.2386142409791600e-001,
                   3.4647074816760000e-002, 1.0161119296278000e-002, 3.9575478735694296e-001, 3.0792998388043608e-001,
                   2.6456694840652000e-001, 3.5853935220595100e-001, 1.5780740596859500e-001, 7.5050596975911008e-002,
                   1.4242160111338298e-001, 6.5494628082938000e-002, 6.0063379479464496e-001, 5.5760326158878384e-001,
                   7.2098702581736496e-001, 5.9452706895587104e-001, 8.3933147368083808e-001, 7.0108797892617304e-001,
                   8.2293132406985696e-001, 9.2434425262078384e-001, 6.0063379479464496e-001, 5.5760326158878384e-001,
                   7.2098702581736496e-001, 5.9452706895587104e-001, 8.3933147368083808e-001, 7.0108797892617304e-001,
                   8.2293132406985696e-001, 9.2434425262078384e-001, 3.9575478735694296e-001, 3.0792998388043608e-001,
                   2.6456694840652000e-001, 3.5853935220595100e-001, 1.5780740596859500e-001, 7.5050596975911008e-002,
                   1.4242160111338298e-001, 6.5494628082938000e-002, 1.6453165694459500e-002, 5.1653659456360000e-003,
                   1.1193623631508002e-002, 1.5133062934734000e-002, 1.5245483901099000e-002, 1.2079606370820500e-002,
                   8.0254017934004992e-003, 4.0422901308920000e-003, 1.0396810137425000e-003, 5.1653659456360000e-003,
                   1.1193623631508002e-002, 1.5133062934734000e-002, 1.5245483901099000e-002, 1.2079606370820500e-002,
                   8.0254017934004992e-003, 4.0422901308920000e-003, 1.0396810137425000e-003, 5.1653659456360000e-003,
                   1.1193623631508002e-002, 1.5133062934734000e-002, 1.5245483901099000e-002, 1.2079606370820500e-002,
                   8.0254017934004992e-003, 4.0422901308920000e-003, 1.0396810137425000e-003, 1.9424384524905000e-003,
                   1.2787080306011000e-002, 4.4404517866690008e-003, 8.0622733808655008e-003, 1.2459709087455000e-003,
                   9.1214200594755008e-003, 5.1292818680995000e-003, 1.8999644276510000e-003, 1.9424384524905000e-003,
                   1.2787080306011000e-002, 4.4404517866690008e-003, 8.0622733808655008e-003, 1.2459709087455000e-003,
                   9.1214200594755008e-003, 5.1292818680995000e-003, 1.8999644276510000e-003, 1.9424384524905000e-003,
                   1.2787080306011000e-002, 4.4404517866690008e-003, 8.0622733808655008e-003, 1.2459709087455000e-003,
                   9.1214200594755008e-003, 5.1292818680995000e-003, 1.8999644276510000e-003, 1.9424384524905000e-003,
                   1.2787080306011000e-002, 4.4404517866690008e-003, 8.0622733808655008e-003, 1.2459709087455000e-003,
                   9.1214200594755008e-003, 5.1292818680995000e-003, 1.8999644276510000e-003, 1.9424384524905000e-003,
                   1.2787080306011000e-002, 4.4404517866690008e-003, 8.0622733808655008e-003, 1.2459709087455000e-003,
                   9.1214200594755008e-003, 5.1292818680995000e-003, 1.8999644276510000e-003, 1.9424384524905000e-003,
                   1.2787080306011000e-002, 4.4404517866690008e-003, 8.0622733808655008e-003, 1.2459709087455000e-003,
                   9.1214200594755008e-003, 5.1292818680995000e-003, 1.8999644276510000e-003]
        }

        A = np.array(table[f'{p}'])
        n = int(A.shape[0] / 3)
        X = [A[0:n].T, A[n + 0:2 * n].T]
        W = np.array([A[2 * n:3 * n]]).T

        return X, W


if __name__ == '__main__':
    folder = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\NativeEig'
    mod = Model(folder, 'TESLA')
    n_cells = 5
    midC3795 = np.array([62.22, 66.13, 30.22, 23.11, 72, 93.5, 171.2]) * 1e-3
    endC3795 = np.array([62.58, 57.54, 17.21, 12.0, 80.0, 93.795, 171.2]) * 1e-3

    # midTESLA = np.array([42, 42, 12, 19, 35, 57.7, 103.3]) * 1e-3
    # midTESLA_800 = np.array([68.12, 68.12, 19.46, 30.81, 56.76, 93.5, 167.646]) * 1e-3
    # midFCCUROS5 = np.array([67.72, 57.45, 21.75, 35.95, 60, 93.5, 166.591]) * 1e-3
    #
    # midDegen = np.array([50.052, 36.5, 7.6, 10.0, 30.0, 57.7, 98.58, 0]) * 1e-3
    # endDegen = np.array([50.052, 36.5, 7.6, 10.0, 30.0, 57.7, 98.58, 0]) * 1e-3
    #
    # midNLSF_RE = np.array([49, 35.30, 10.5, 17, 32.0, 57.7, 98.58, 0]) * 1e-3
    # endNLSF_RE = np.array([50, 35, 10, 15, 32.0, 57.7, 98.58, 0]) * 1e-3
    #
    # mid_cell_parameters = [43.99, 35.06, 12.53, 20.95, 35, 57.6524, 101.205]  # [A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m]
    # left_end_cell_parameters = [50.9, 45.3, 8.4, 11.5, 39, 59.988, 101.205]  # [A_e, B_e, a_e, b_e, Ri_e, L_e, Req_e]
    # right_end_cell_paramters = [52.1, 47.9, 9.9, 11.3, 37, 62.665, 101.205]

    mod.run(n_cells, midC3795, endC3795, endC3795, beampipe='both', gridcons=0.005, plot=True)
