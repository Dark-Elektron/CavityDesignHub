import os
import shutil
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.io as spio
import pandas as pd
from icecream import ic

m0 = 9.1093879e-31
q0 = 1.6021773e-19
c0 = 2.99792458e8


class FS:
    # parameters for the MP analysis

    def __init__(self, folder):
        self.folder = folder
        pass

    def create_inputs(self):
        files_list = ["initials", "flevel", "y00", "param", "fieldparam", "counter_flevels.mat", "counter_initials.mat",
                      "gene_initials.mat", "gene_temp.mat", "initials_temp.mat", "distance_flevel.mat",
                      "counter_initialsl.mat", "gene_initialsl.mat", "initials_templ.mat", "distance_flevell.mat",
                      "counter_initialsr.mat", "gene_initialsr.mat", "initials_tempr.mat", "distance_flevelr.mat",
                      "dx.mat", "alphas.mat"]
        for file in files_list:
            if os.path.exists(fr"{self.folder}\{file}"):
                os.remove(fr"{self.folder}\{file}")

        ic('Old input data deleted.')

        # Inputs for eigenmode solver
        freq = 801.58e6
        gtype = 1
        # save the inputs 
        # save_values

        my0    = 4*np.pi*1e-7 
        eps0   = 8.85418782e-12 
        lambda_ = 1 / (freq*np.sqrt(eps0*my0))
        eta0   = 376.7303134111465
        epsr = 1  # relative epsilon
        strech = 0
        d1 = 1  # Grid constant
        d2 = 1  # Mesh constant
        R = -1 + 0*1j
        v0 = 2  # initial velocity
        N = 20  # Number of impacts
        Z = 0
        
        fieldparam = [gtype, freq, epsr, d1/1000, R.real, R.imag, strech, Z]
        # save -ascii fieldparam fieldparam
        df = pd.DataFrame(fieldparam)
        df.to_csv(fr'{self.folder}\fieldparam', index=False, header=False, float_format='%.7E')
        
        # field levels
        # flevel = np.arange(flmin, flmax, flstep).T*1e3
        
        # # save counter_flevels flevel
        # df = pd.DataFrame(flevels)
        # df.to_csv(fr'{self.folder}\flevels', index=False, header=False, float_format='%.7E')
        
        # parameters for the MP analysis
        m = 9.1093879e-31
        q = 1.6021773e-19
        c = 2.99792458e8
        param = np.zeros((7, 1))
        V0 = 1                               # intensity of the EM field
        # V0    = 2  error_message('POISTA TÄMÄ - TESTI (save_values.m)')
        gamma = c/freq/(2*np.pi)                   # distance/phase -coefficient
        v00 = c*np.sqrt(1 - 1/((v0 * q/(m*c**2)+1)**2))  # initial velocity (relativistic)
        # v00   = sqrt(2*v0*q/m)                 # (classical)
        ctype = 1                               # compute counter functions
        tol = 1e-3                            # tolerance for the ODE solver
        emin = 1  # not required for eigenmode analysis
        emax = 10  # not required for eigenmode analysis
        param = [freq, V0, gamma, v00, N, ctype, tol, emin, emax]
        # save -ascii param param
        df = pd.DataFrame(param)
        df.to_csv(fr'{self.folder}\param', index=False, header=False, float_format='%.7E')
        
        # grid constant for creating a new grid
        # load geodata.n
        geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()

        geodata[0, 0] = d2/1000
        # save -ascii geodata.n geodata
        df = pd.DataFrame(geodata)
        df.to_csv(fr'{self.folder}\geodata.n', sep=r' ', index=False, header=False, float_format='%.7E')
        
        # # parameters for the initial point generator
        # dx = dx/1000                         # dimensions in mm
        # dc = dx/10                           # distance from the nearest corner
        # alpha = 0                               # initial velocity angle
        # dt = dphi/freq/360                   # time step
        # initials = [-dx,dc,alpha,dt,zmin,zmax]
        # # save gene_initials initials
        # df = pd.DataFrame(initials)
        # df.to_csv(fr'{self.folder}\initials', index=False, header=False, float_format='%.7E')
        # -----------------------------------------------------------------------

    def mesh_cavity(self, geodata, gridcons, epsr, s):
        ic("inside mesh cavity")
        n = len(geodata[:, 1]) - 1
        gr = geodata[3:n, 0]
        gz = geodata[3:n, 1]
        gn = geodata[3:n, 2]
        n = len(gr)

        PEC = np.where(gn == 1)[0]       # electric walls
        WIN = np.where(gn == 2)[0]       # dielectric walls
        PMC = np.where(gn == 3)[0]       # magnetic walls

        if len(WIN) > 0:
            ic('Dielectric boundary found in a cavity.')
            return

        # specify the job for the eigenvalue solver
        if len(PMC) > 0:
            if s > 0:
                ic('A cavity with magnetic walls.')
                job = 1
            elif len(PEC) > 0 & len(PMC) == 0:
                if s > 0:
                    ic('A cavity with electric walls.')
                    job = 0
            else:
                ic('Invalid geometry.')
                return
        job = {"job": 1}
        # save job job
        spio.savemat(f'{self.folder}/job.mat', job, format='4')

        # First the nodes of the geometry
        nodes = np.array([gz, gr]).T

        # Then the edges, the third column secifies the type of an edge as follows:
        #   1 : nxE = 0 on boundary
        #   3 : nxH = 0 on boundary
        #   0 : for streching

        edges = np.array([np.arange(1, n+1), np.append(np.arange(2, n+1), 1), gn]).T
        # ic(edges)
        ne = len(edges[:, 0])

        # And finally the pacthes, first the edges are listed, last three values
        # give relative epsilon, relative mu and patch type, 0 = can strecth,
        # 1-n = regulars
        patches = np.append(np.arange(1, ne+1), [epsr, 1, 1])
        # ic(patches)

        esize = gridcons
        esize = 0.005

        # save model.mat
        model = {'nodes': nodes,
                 "edges": edges,
                 "patches": patches,
                 "esize": esize
                 }

        spio.savemat(f"{self.folder}/model.mat", model, format='4')

        # start the mesh generator
        ic("It's here")
        subprocess.call(f"{self.folder}/2dgen_bin.exe", cwd=self.folder,
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        ic("done done")
        print("Done done")

    def plot_mesh(self, s):
        # Function program plot_mesh(s)
        # --------------------------------------------------------------------
        # Plots the mesh.
        #
        # --------------------------------------------------------------------
        # CALLS TO : print.m, clear_window.m
        # xx/yy/00 : Seppo Järvempää
        # 04/04/00 : Pasi Ylä-Oijala - m-file
        # --------------------------------------------------------------------
        s = 0
        # if nargin < 1:
        #   s = 1

        cl = 0
        ok1 = os.path.exists(fr"{self.folder}\mesh.mat")
        ok2 = os.path.exists(fr"{self.folder}\fieldparam")
        print(ok1, ok2)
        if not ok1:
            print(['The mesh does not exists. Choose Mesh Generator in menu Run.'])
        elif not ok2:
            print('Parameter file fieldparam does not exist.')
        else:
            gtype = self.gtype
            if gtype == 1:
                if s > 0:
                    print('Plotting the mesh.')
                #      print('
                else:
                    if s > 0:
                        print('Plotting the mesh. Blue area for streching.')
                    #      print('                  ')

            # plots 2-d mesh in mesh.mat, which includes coord and etopol -arrays
            # load mesh
            mesh = spio.loadmat(fr"{self.folder}\mesh.mat")
            coord = mesh['coord']
            etopol = mesh['etopol']
            alue = mesh['alue']
            tyyppi = mesh['tyyppi']
            boundary = mesh['boundary']
            edges = mesh['edges']

            xmin = min(coord[:, 0])
            xmax = max(coord[:, 0])
            ymin = min(coord[:, 1])
            ymax = max(coord[:, 1])
            # self.ax.plot([xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax], 'k.')

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

                osa = (x2 -x1) * (y3 -y1) -(x3 -x1) *(y2 -y1)
                ala = ala + osa
                self.ax.plot([x1, x2, x3, x1], [y1, y2, y3, y1], 'k')

            I = np.where(alue.T[0] == 0)[0]
            self.ax.fill(X[:, I], Y[:, I], 'b')
            I = np.where(tyyppi.T[0] == 1)[0]
            self.ax.fill(X[:, I], Y[:, I], 'r')
            I = np.where(tyyppi.T[0] == 2)[0]
            self.ax.fill(X[:, I], Y[:, I], 'g')

            I = np.where(edges[:, 2] > 0)[0]  # no bouncing
            for i in range(len(I)):
                i1 = int(edges[I[i], 0]) - 1
                i2 = int(edges[I[i], 1]) - 1
                # ic(i1, np.shape(i1))
                self.ax.plot([coord[i1, 0], coord[i2, 0]], [coord[i1, 1], coord[i2, 1]], 'r')

            I = np.where(boundary == 3)
            for i in range(len(I)):
                self.ax.plot(coord[I[i], 0], coord[(I[i ] -1, 1)], 'b*')

            I = np.where(boundary == 0)[0]
            self.ax.plot(coord[I, 0], coord[I, 1], 'w*')
            ala = ala / 2
            self.ax.set_title('[MultiPac 2.0                    Mesh                  date ]')
            self.ax.set_xlabel('z axis [m]')
            self.ax.set_ylabel('r axis [m]')
            # self.plt.ax.set_axis('image')
            # self.ax.set_colormap('jet')
            self.fig.canvas.draw_idle()

    def run_field_solver(self):
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

        self.create_inputs()

        # load geodata.n
        geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        gridcons1 = geodata[0, 0]

        # load fieldparam
        fieldparam = pd.read_csv(fr"{self.folder}\fieldparam", sep='\s+', header=None).to_numpy()
        self.gtype = fieldparam[0][0]
        self.freq = fieldparam[1][0]
        self.epsr = fieldparam[2][0]
        self.gridcons = fieldparam[3][0]

        s = 0

        self.fig, self.ax = plt.subplots()

        # gtype = 'Cavity'
        # self.gtype = gtype
        # freq = 801.58e6
        # releps = 1
        # gridcons = 5e-3
        # self.epsr = releps

        if self.epsr > 1:
            print('Relative permittivity > 1?')

        # generate the mesh
        self.mesh_cavity(geodata, gridcons1, self.epsr, s)

        # plot the mesh
        self.plot_mesh(0)
        ic("Here now done with plotting mesh")

        # find resonance solution
        ic('Computing the eigen values.')
        job = 1
        freqn = self.eigen(job, self.freq)
        ic(freqn)
        # k1, k2 = self.find_resonance(freq)

        if s > 0:
            ic('Eigen values calculated.')
            ic('                        ')

        err = (abs(freqn - self.freq) / self.freq) * 100
        # ic(err)
        if err > 1:
            ic('Warning: Error in eigen frequency more than 1#.')

        # compute and plot the fields
        self.calculate_fields(0, gridcons1, 0)

        ic("Now to plot the fields")
        self.plot_FEM_fields(0, gridcons1, s)
        # ---------------------------------------------------------------------

    def eigen(self, jobl, freq):
        maara = 10
        job = {"job": 1}
        offset = {"offset": 0}
        spio.savemat(f"{self.folder}/job.mat", job, format='4')
        spio.savemat(f"{self.folder}/offset.mat", offset, format='4')
        # options.disp = 0

        # !eigenC_bin
        cwd = fr'{self.folder}'
        eigenCpath = fr'{self.folder}\eigenC_bin.exe'
        if os.path.exists(eigenCpath):
            subprocess.call(eigenCpath, cwd=cwd,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        ic("Done with eigenmode analysis")

        cinfo = spio.loadmat(fr"{self.folder}\cinfo.mat")['cinfo'][0]
        if cinfo != 0:
            ic(cinfo)
            ic('Eigenvalue solver failed.')

        # load o_eigen
        o_eigen = spio.loadmat(fr"{self.folder}\o_eigen.mat")
        n = int(o_eigen['n'])
        ia = o_eigen['ia'].T[0] - 1  # mnatlab indices to python
        ja = o_eigen['ja'].T[0] - 1
        Aarvot = o_eigen['Aarvot'].T[0]
        Barvot = o_eigen['Barvot'].T[0]
        nf = o_eigen['nf']
        index = o_eigen['index']
        siirto = o_eigen['siirto'][0]
        # ic(index.shape, Aarvot.shape, Barvot.shape)

        # ic(ia, ja, Aarvot, n)
        # ic(np.shape(ia), np.shape(ja), np.shape(Aarvot))
        AA = sps.csr_matrix((Aarvot, (ia, ja)), shape=(n, n))
        BB = sps.csr_matrix((Barvot, (ia, ja)), shape=(n, n))

        ic("Here here here here")
        d2, u2 = spsl.eigs(AA, M=BB, k=maara, sigma=0)
        ic("It is now where", u2)
        u2 = u2.real

        # ic(u2.shape)
        d2 = np.absolute(d2)

        mu0 = 4 * np.pi * 1e-7
        e0 = 8.85418782e-12
        k0 = 2 * np.pi * freq * np.sqrt(mu0 * e0)

        k = np.sqrt(d2)
        eigenvalues = k[0:maara]
        # ic(eigenvalues)
        k2 = k

        val = np.min(abs(k0 - k2))
        ind = np.argmin(abs(k0 - k2))
        k = k2[ind]
        u = u2[:, ind]
        # ic(k)

        ic("It is now where2")
        # new frequency
        freqn = k / (2 * np.pi * np.sqrt(mu0 * e0))
        ic("It is now where3")
        # ic(freqn)
        param = pd.read_csv(fr"{self.folder}\param", sep='\s+', header=None).to_numpy().T[0]
        ic("It is now where4")
        # ic(param)
        # load param
        param[0] = freqn
        # save -ascii param param
        ic("It is now where3")
        df = pd.DataFrame(param)
        df.to_csv(fr'{self.folder}\param', index=False, header=False, float_format='%.7E')

        ic("It is now where")
        fieldparam = pd.read_csv(fr"{self.folder}\fieldparam", sep='\s+', header=None).to_numpy().T[0]
        # load fieldparam
        fieldparam[1] = freqn
        # ic(freqn)
        # save -ascii fieldparam fieldparam
        df = pd.DataFrame(fieldparam)
        df.to_csv(fr'{self.folder}\fieldparam', index=False, header=False, float_format='%.7E')

        # save kama0 index -v4
        # save kama0 u -append
        # save kama0 k -append
        kama0 = {"index": index,
                 "u": u,
                 "k": k
                 }

        spio.savemat(f"{self.folder}/kama0.mat", kama0, format='4')
        ic("here now and here")
        # --------------------------
        ic(freqn)
        return freqn

    def calculate_fields(self, s, gridcons, ss=1):

        ok3 = 1
        if os.path.exists(fr"{self.folder}\fieldparam"):
            ok3 = 0

        if ok3 == 0:
            fieldparam = pd.read_csv(fr"{self.folder}\fieldparam", sep='\s+',
                                     header=None).to_numpy().T[0]
            gtype = fieldparam[0]
            s = 0

        # if not os.path.exists('s'):
        #     s = 0

        ok1 = os.path.exists(fr"{self.folder}\mesh.mat")
        if not ok1:
            ic(['The mesh does not exists. Choose Mesh Generator in menu Run.'])

        ok2 = 1
        if s == 0:
            if os.path.exists(fr"{self.folder}\kama0.mat"):
                ok2 = 0
            if ok2 == 0:
                ic(['The fields do not exist. Choose Field Solver in menu Run.'])

            elif s == 1:
                ok21 = os.path.exists(fr"{self.folder}\kama1.mat")
                ok22 = os.path.exists(fr"{self.folder}\kama2.mat")
                ok2 = ok21 * ok22
                if ok2 == 0:
                    ic(['The fields do not exist. Choose Field Solver in menu Run.'])

        if ok1 == 0 and ok2 == 0 and ok3 == 0:
            ok = os.path.exists(fr"{self.folder}\geodata.n")
            if ok:
                ic('Computing the fields.')
                #  error_message('                                  ')

                # define the grid constant
                # if gtype <= 2:
                geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
                gridcons = geodata[0, 0]
                # else:
                #     geodata = pd.read_csv(fr"{self.folder}\geodatal.n", sep='\s+',
                #                           header=None).to_numpy()
                #     gridcons = geodata[0, 0]

        # load mesh
        # ic(gtype)
        mesh = spio.loadmat(fr"{self.folder}\mesh.mat")
        coord = mesh['coord']
        etopol = mesh['etopol']
        alue = mesh['alue']
        tyyppi = mesh['tyyppi']
        boundary = mesh['boundary']
        edges = mesh['edges']

        # if gtype <= 2:
        geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+',
                              header=None).to_numpy()
        n = len(geodata[:, 0])
        gr = geodata[3:n, 0]
        gz = geodata[3:n, 1]
        # else:
        #     geodatal = pd.read_csv(fr"{self.folder}\geodatal.n", sep='\s+',
        #                            header=None).to_numpy()
        #     nl = len(geodatal[:, 0]) - 1
        #     grl = geodatal[3:nl, 0]
        #     gzl = geodatal[3:nl, 1]
        #     geodataw = pd.read_csv(fr"{self.folder}\geodataw.n", sep='\s+',
        #                            header=None).to_numpy()
        #     nw = len(geodataw[:, 0]) - 1
        #     grw = geodataw[3:nw, 0]
        #     gzw = geodataw[3:nw, 0]
        #     geodatar = pd.read_csv(fr"{self.folder}\geodatar.n", sep='\s+',
        #                            header=None).to_numpy()
        #     nr = len(geodatar[:, 0]) - 1
        #     grr = geodatar[3:nr, 0]
        #     gzr = geodatar[3:nr, 1]
        #     gz = np.sort([gzl, gzw, gzr])
        #     gr = np.sort([grl, grw, grr])

        # generate field points
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

        #     save zr r -v4 -append
        #     save zr alue -v4 -append
        #     save zr rind -v4 -append
        #     save zr zind -v4 -append
        zr = {'z': z, 'r': r, 'alue': alue, 'rind': rind, 'zind': zind}
        spio.savemat(fr'{self.folder}\zr.mat', zr, format='4')

        # compute the fields at generated points
        if s == 0:
            ic("Calculating fields")
            # !copy /y kama0.mat kama.mat
            shutil.copyfile(fr'{self.folder}\kama0.mat', fr'{self.folder}\kama.mat')

            # !Multipac fields
            cwd = fr'{self.folder}'
            multipacPath = fr'{self.folder}\Multipac.exe'
            if os.path.exists(multipacPath):
                subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd,
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            Er = spio.loadmat(fr"{self.folder}\Er.mat")['Er']
            Ez = spio.loadmat(fr"{self.folder}\Ez.mat")['Ez']
            H = spio.loadmat(fr"{self.folder}\H.mat")['H']

            fields = {'I1': I1, 'I2': I2, 'rr': rr, 'zz': zz, 'r': r, 'z': z, 'Er': Er, 'Ez': Ez, 'H': H}

            spio.savemat(fr'{self.folder}\fields.mat', fields)  # this was not saved as Mat object version 4. Ask why?

            self.save_fields(0, 1)
            ic("It's now hwere before normalize")
            self.calculate_QoI()
            # self.normalize_u(0)
        # else:
        #     job = 0
        #     # save job job -v4
        #     # !copy /y kama1.mat kama.mat
        #     # !Multipac fields
        #     cwd = fr'{self.folder}'
        #     multipacPath = fr'{self.folder}\Multipac.exe'
        #     if os.path.exists(multipacPath):
        #         subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd)
        #
        #     Er = spio.loadmat(fr"{self.folder}\Er.mat")['Er']
        #     Ez = spio.loadmat(fr"{self.folder}\Ez.mat")['Ez']
        #     H = spio.loadmat(fr"{self.folder}\H.mat")['H']
        #     ic(Er.shape)
        #
        #     fields = {'I1': I1, 'I2': I2, 'rr': rr, 'zz': zz, 'r': r, 'z': z, 'Er': Er, 'Ez': Ez, 'H': H}
        #     spio.savemat(fr'{self.folder}\fields', fields, format='4')
        #     self.save_fields(1, job)
        #     self.normalize_u(1)

        # save zr z -v4
        # save zr r -v4 -append
        # save zr alue -v4 -append
        # save zr rind -v4 -append
        # save zr zind -v4 -append
        #
        # job = 1
        # save job job -v4
        # !copy /y kama2.mat kama.mat
        # !Multipac fields
        # load Er
        # load Ez
        # load H
        # save fields2 I1 I2 rr zz r z Er Ez H
        # save_fields(2)
        # normalize_u(2)
        # end
        #
        # if ss > 0
        # cl = error_message('Fields are computed and saved.')
        # cl = error_message(['To plot the fields, choose Plot FEM Fields in '...
        # 'Menu Fields.'])
        # cl = error_message('                               ')
        #
        # end
        # end
        # end
        #
        # if cl == 1
        # clear_window
        # end

    def plot_FEM_fields(self, ptype, gridcons, s=0):
        ok1 = 1
        if os.path.exists(fr'{self.folder}/geodata.n'):
            ok1 = 0
        ok2 = 0
        if ok1 == 0:
            files = ['fields.mat', 'Er.mat', 'Ez.mat', 'H.mat', 'fieldparam']
            for f in files:
                if os.path.exists(fr'{self.folder}/{f}'):
                    continue
                else:
                    ok2 = 1

        # ic(ok1, ok2)
        if ok1 == 0 and ok2 == 0:
            if s == 0:
                if ptype == 0:
                    ic('Plotting the fields. A pcolor plot.')
                elif ptype == 1:
                    ic('Plotting the fields. An arrow plot.')

                geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
                gridcons = geodata[0, 0]
                self.plot_cavity_fields(ptype, gridcons, s)

            if s > 0:
                ic('                                                  ')

    def plot_cavity_fields(self, s, gridcons, ss):
        # ----------------------------------------------------------------------
        # Plots the fields computed by the field solver. For a cavity.
        #
        # INPUT  s : 0 - pcolor plots
        #            1 - arrow plots
        #  gridcons : grid constant, [mm]
        # ----------------------------------------------------------------------
        # CALLS TO : arrow.m, error_message.m, clear_window.m
        # 17/04/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
        # ----------------------------------------------------------------------

        geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        n = len(geodata[:, 0])
        gr = geodata[3:n, 0]
        gz = geodata[3:n, 1]

        # load the field values
        # load fields Er Ez H I1 I2 zz rr z r
        fields = spio.loadmat(fr"{self.folder}\fields.mat")
        Er = fields['Er']
        Ez = fields['Ez']
        H = fields['H']
        I1 = fields['I1']
        I2 = fields['I2']
        zz = fields['zz']
        rr = fields['rr']
        z = fields['z']
        r = fields['r']
        # ic(zz, type(zz), np.shape(zz))

        # plot the fields
        if s == 0:                                   # a pcolor plot
            # ic(I2.shape, I1.shape, Er.shape)
            Er = np.reshape(Er, (max(I1.shape), max(I2.shape)))
            Ez = np.reshape(Ez, (max(I1.shape), max(I2.shape)))
            ic(Ez)
            EE = np.sqrt(abs(Er)**2 + abs(Ez)**2).T

            # Calculate accelerating field
            L_half_cell = 0.0935  # length of cavity half cell
            E_axis = EE[rr == 0]  # & zz>=-L_half_cell & zz <= L_half_cell
            E_axis = np.sqrt(Ez**2).T[rr == 0]  # & zz>=-L_half_cell & zz <= L_half_cell
            plt.plot(E_axis)
            plt.show()
            z_slice = zz[0, :]
            z_active = z_slice[(z_slice >= -L_half_cell) & (z_slice <= L_half_cell)]

            Vacc = np.trapz(E_axis * np.exp(1j * 2 * np.pi * 801.58e6 * z_slice / c0), z_slice)
            Eacc = abs(Vacc)/(2*5*L_half_cell)
            Epk_Eacc = np.max(EE)/Eacc

            Hmag = np.sqrt(H**2).T
            Hpk_Eacc = np.max(Hmag)*1e-3/(Eacc*1e-6)

            print(np.max(EE), Eacc, Epk_Eacc)
            print(np.max(Hmag), Eacc, Hpk_Eacc)

            fig, axs = plt.subplots(2, 1)
            pcE = axs[0].pcolor(zz, rr, EE, cmap='RdBu')
            # pcE = axs[0].contourf(zz, rr, EE, cmap='RdBu', levels=40)
            # color = 2 * np.log(np.hypot(Ez, Er))
            # pcE = axs[0].streamplot(zz, rr, Ez, Er, linewidth=1, cmap=plt.cm.inferno,
            #               density=2, arrowstyle='->', arrowsize=1.5)
            plt.colorbar(pcE, ax=axs[0])

            axs[0].plot(gz, gr, '-b')

            axs[0].set_xlabel('z axis [m]')
            axs[0].set_ylabel('r axis [m]')
            axs[0].set_title('MultiPac 2.0          Electric field   abs (E)   [V/m]      ')

            HH = np.reshape(H, (max(I1.shape), max(I2.shape))).T
            pcH = axs[1].pcolor(zz, rr, 4e-7*np.pi*HH, cmap='RdBu')
            plt.colorbar(pcH, ax=axs[1])

            axs[1].plot(gz, gr, '-b')
            axs[1].set_title('Magnetic field     B_\phi  [TESLA]')
            axs[1].set_xlabel('z axis [m]')
            axs[1].set_ylabel('r axis [m]')
            plt.show()

            print(Ez)
            print(Ez.shape)
            print(zz)
            print(zz.shape)
            plt.spy(zz)
            plt.show()
            plt.spy(Ez.T)
            plt.show()
            plt.spy(Er.T)
            plt.show()

        else:                                   # an arrow plot
            F = [Ez.T, Er.T]
            A = [z.T, r.T]

            fig, ax = plt.subplots()
            # arrow(F,A,gridcons,'f','s','k')
            ax.quiver(z, r, Ez, Er)
            ax.plot(gz, gr, '-b')
            ax[1].set_xlabel('z axis [m]')
            ax[1].set_ylabel('r axis [m]')
            ax[1].set_title('MultiPac 2.0            Electric field        [V/m]         ')

        if ss > 0:
            ic('Electric field plotted. To see the magnetic field choose Pcolor.')

    def calculate_QoI(self):

        geodata = pd.read_csv(fr"{self.folder}\geodata.n", sep='\s+', header=None).to_numpy()
        n = len(geodata[:, 0])
        gr = geodata[3:n, 0]
        gz = geodata[3:n, 1]

        # load the field values
        # load fields Er Ez H I1 I2 zz rr z r
        fields = spio.loadmat(fr"{self.folder}\fields.mat")
        Er = fields['Er']
        Ez = fields['Ez']
        H = fields['H']
        I1 = fields['I1']
        I2 = fields['I2']
        zz = fields['zz']
        rr = fields['rr']
        z = fields['z']
        r = fields['r']
        ic(zz, rr, zz.shape, rr.shape)
        ic(Er, Er.shape)

        # plot the fields
        Er = np.reshape(Er, (max(I1.shape), max(I2.shape)))
        ic(Er)
        Ez = np.reshape(Ez, (max(I1.shape), max(I2.shape)))
        EE = np.sqrt(abs(Er) ** 2 + abs(Ez) ** 2).T
        U = 3

        # Calculate accelerating field
        L_half_cell = 0.187  # length of cavity half cell
        E_axis = EE[rr == 0]  # & zz>=-L_half_cell & zz <= L_half_cell
        z_slice = zz[0, :]
        z_active = z_slice[(z_slice >= -L_half_cell) & (z_slice <= L_half_cell)]
        Vacc = np.trapz(z_slice, E_axis)
        Eacc = Vacc / (2 * L_half_cell)
        Epk_Eacc = np.max(EE) / Eacc

    def find_resonance(self, freq):
        maara = 10           # number of eigenvalues
        raja = 1e-3           # error tolerance

        EKENTTA = 0
        HKENTTA = 1

        mu0 = 4*np.pi*1e-7
        e0 = 8.85418782e-12
        k0 = 2*np.pi*freq*np.sqrt(mu0*e0)            # "correct" eigenvalue

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

    def save_fields(self, wall, job):
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

        ic(E0)

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
                offset1 = offset1-0.9

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
                offset2 = 0.5*(offset1 + offset3)
                s3 = s2
            elif s2 != s3:
                s1 = s2
                offset1 = offset2
                offset2 = 0.5*(offset2 + offset3)
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
            ero2 = np.min(abs(k-korig))
            ind = np.argmin(abs(k-korig))
            if abs(k[ind]) < 1e-6:
                ind = ind+1

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
                offset1 = offset1-0.9

                k3 = k
                offset3 = offset2
                s3 = s2
                k = k1
                offset2 = offset1
                s2 = s1
                k1 = self.eee(offset1, maara, 0)
                ero1 = korig-k1[ind]
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
                offset2 = 0.5*(offset1+offset3)
                s3 = s2
            elif s2 != s3:
                s1 = s2
                offset1 = offset2
                offset2 = 0.5*(offset2+offset3)
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


if __name__ == '__main__':
    folder = fr'D:\Dropbox\CavityDesignHub\em_codes\customEig\run_files'
    fs = FS(folder)
    fs.run_field_solver()
    plt.show()