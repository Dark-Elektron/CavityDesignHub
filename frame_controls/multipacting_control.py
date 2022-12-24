import ast
import json
import os
import subprocess
import time
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.io as spio
import pandas as pd
from threading import Thread
from PyQt5 import QtCore
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from icecream import ic
from scipy.optimize import fsolve
from termcolor import colored
from graphics.scene import Scene
from modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from ui_files.multipacting import Ui_Multipacting
from utils.file_reader import FileReader
import psutil
from modules.plot_module.plotter import Plot

# evol = Evolution()
slans_geom = SLANSGeometry()
fr = FileReader()

file_color = 'red'
DEBUG = True


def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))


class MultipactingControl:
    def __init__(self, parent):
        print("Check 1: multipacting_control.py")
        self.w_Multipacting = QWidget()

        self.ui = Ui_Multipacting()
        self.ui.setupUi(self.w_Multipacting)
        print("Check 2: multipacting_control.py")

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui
        print("Check 3: multipacting_control.py")

        # get logger
        self.log = self.main_control.log

        # ###########################
        # Create Scene
        self.scene = Scene(self)

        # ##########################
        # Add 2d graph
        self.plt = Plot(self.ui)
        self.fig = self.plt.fig
        self.ax = self.plt.ax
        self.ui.gl_Plot_Area.addWidget(self.plt)
        print("Check 4: multipacting_control.py")

        # add 3D
        # self.plotter = QtInteractor()
        # self.ui.gl_Display_3D.addWidget(self.plotter.interactor)
        # self.plotter.show_grid()
        print("Check 5: multipacting_control.py")

        self.initUI()
        self.signals()
        # self.exe_control()

        # instantiate geometry
        self.slans_geom = SLANSGeometry()

        # shape space initialization
        self._shape_space = {}
        self._selected_keys = []
        self.processes = []
        self.processes_id = []
        self.show_progress_bar = False

        # ui effects
        # self.ui_effects()

    def initUI(self):

        # init shape entry mode
        self.shape_entry_widgets_control()

        # disable expansion section for now. Feature to come later
        self.ui.cb_Expansion.setEnabled(False)

        # inner cell
        self.ui.cb_Inner_Cell.setCheckState(2)
        self.ui.cb_Inner_Cell.setEnabled(False)

        # expand/collapse sections widgets
        if self.ui.cb_Expansion.checkState() == 2:
            self.ui.w_Expansion.setMinimumHeight(160)
        else:
            self.ui.w_Expansion.setMinimumHeight(0)
            self.ui.w_Expansion.setMaximumHeight(0)

        if self.ui.cb_Outer_Cell_L.checkState() == 2:
            self.ui.w_Outer_Cell_L.setMinimumHeight(160)
        else:
            self.ui.w_Outer_Cell_L.setMinimumHeight(0)
            self.ui.w_Outer_Cell_L.setMaximumHeight(0)

        if self.ui.cb_Outer_Cell_R.checkState() == 2:
            self.ui.w_Outer_Cell_R.setMinimumHeight(160)
        else:
            self.ui.w_Outer_Cell_R.setMinimumHeight(0)
            self.ui.w_Outer_Cell_R.setMaximumHeight(0)

        # create pause and resume icons to avoid creating them over and over again
        self.pause_icon = QIcon()
        self.pause_icon.addPixmap(QPixmap(f":/icons/icons/PNG/pause.png"), QIcon.Normal, QIcon.Off)
        self.resume_icon = QIcon()
        self.resume_icon.addPixmap(QPixmap(f":/icons/icons/PNG/resume.png"), QIcon.Normal, QIcon.Off)

        # process state
        self.process_state = 'none'
        self.run_pause_resume_stop_routine()

        # create progress bar object and add to widget
        # self.progress_bar = QProgressBar(self.ui.w_Simulation_Controls)
        # self.progress_bar.setMaximum(100)
        # self.progress_bar.setValue(0)
        # self.ui.gl_Simulation_Controls.addWidget(self.progress_bar, 0, 4, 1, 1)
        # self.progress_bar.hide()

    def signals(self):
        # run multipacting solver
        self.ui.pb_Run.clicked.connect(lambda: self.run_multipac())

        # load shape space
        self.ui.pb_Select_Shape_Space.clicked.connect(
            lambda: self.open_file(self.ui.le_Shape_Space, self.ui.cb_Shape_Space_Keys))

        # control shape entry mode
        self.ui.cb_Shape_Entry_Mode.currentIndexChanged.connect(lambda: self.shape_entry_widgets_control())

        # cell parameters control signals
        self.ui.cb_Outer_Cell_L.stateChanged.connect(
            lambda: self.animate_height(self.ui.cb_Outer_Cell_L, self.ui.w_Outer_Cell_L, 0, 160,
                                        True))
        self.ui.cb_Outer_Cell_R.stateChanged.connect(
            lambda: self.animate_height(self.ui.cb_Outer_Cell_R, self.ui.w_Outer_Cell_R, 0, 160,
                                        True))
        self.ui.cb_Expansion.stateChanged.connect(
            lambda: self.animate_height(self.ui.cb_Expansion, self.ui.w_Expansion, 0, 160,
                                        True))

        # cancel
        # self.ui.pb_Cancel.clicked.connect(lambda: self.cancel())
        # self.ui.pb_Pause_Resume.clicked.connect(lambda: self.pause() if self.process_state == 'running' else self.resume())

        # uncomment to draw again
        # self.ui.cb_Shape_Space_Keys.currentTextChanged.connect(lambda: self.draw_shape_from_shape_space())

        #
        self.ui.le_Alpha.editingFinished.connect(lambda: self.update_alpha())

        #
        self.ui.le_Req_i.editingFinished.connect(
            lambda: self.ui.le_Req_ol.setText(self.ui.le_Req_i.text()))
        self.ui.le_Req_i.editingFinished.connect(
            lambda: self.ui.le_Req_or.setText(self.ui.le_Req_i.text()))

    def shape_entry_widgets_control(self):
        if self.ui.cb_Shape_Entry_Mode.currentIndex() == 0:
            self.main_control.animate_height(self.ui.w_Select_Shape_Space, 0, 50, True)

            self.ui.w_Enter_Geometry_Manual.setMinimumHeight(0)
            self.ui.w_Enter_Geometry_Manual.setMaximumHeight(0)

            # clear cells from graphics view
            # self.graphicsView.removeCells()
        else:
            self.main_control.animate_height(self.ui.w_Enter_Geometry_Manual, 0, 375, True)

            self.ui.w_Select_Shape_Space.setMinimumHeight(0)
            self.ui.w_Select_Shape_Space.setMaximumHeight(0)

            # uncomment following lines to draw
            # # clear cells from graphics view
            # self.graphicsView.removeCells()
            #
            # # draw new cell
            # self.graphicsView.drawCells(color=QColor(0, 0, 0, 255))

    def get_parameters(self):
        # Program save_parameters.m
        # -------------------------------------------------------------------------
        # Save the input arguments to file mpgui_inputs.mat.
        #
        # -------------------------------------------------------------------------
        # CALLS TO : check_geodata.m
        # 15/12/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
        # -------------------------------------------------------------------------
        self.gtype = self.ui.cb_Geometry_Type.currentText()
        unit_multiplier = {'Hz': 1, 'kHz': 1e3, 'MHz': 1e6, 'GHz': 1e9}
        unit = self.ui.cb_Frequency_Unit.currentText()
        self.freq = self.ui.dsb_Frequency.value() * unit_multiplier[unit]
        self.epsr = self.ui.dsb_Relative_Epsilon.value()
        self.d1 = self.ui.dsb_Grid_Constant.value()
        self.Rre = self.ui.dsb_Reflection_Coeff_Re.value()
        self.Rim = self.ui.dsb_Reflection_Coeff_Im.value()
        self.R = self.Rre + 1j * self.Rim
        self.d2 = self.ui.dsb_Mesh_Constant.value()  # for multipacting
        self.N = self.ui.sb_Number_Of_Impacts.value()
        self.v0 = self.ui.dsb_Initial_Velocity.value()
        self.emin = self.ui.dsb_Impact_Energy_Min.value()
        self.emax = self.ui.dsb_Impact_Energy_Max.value()
        self.dphi = self.ui.dsb_Initial_Sites_Phase_Step.value()
        self.zmin = self.ui.dsb_Initial_Sites_Min.value()
        self.dx = self.ui.dsb_Initial_Sites_Step.value()
        self.zmax = self.ui.dsb_Initial_Sites_Max.value()
        self.flmin = self.ui.dsb_Field_Levels_Min.value()
        self.flstep = self.ui.dsb_Field_Levels_Step.value()
        self.flmax = self.ui.dsb_Field_Levels_Max.value()

        # save the inputs parameters to the file mpgui_inputs.mat
        # save mpgui_inputs gtype freq epsr d1 R d2 N v0 emin emax dphi zmin dx ...
        #      zmax flmin flstep flmax
        # input_parameters = {'Geometry Type': self.gtype, 'Frequency': freq, 'Relative Epsilon': epsr, 'Grid Constant': d1,
        #                     'Reflection Coefficient': R, 'Mesh Constant': d2, 'Number of Impacts': N,
        #                     'Initial Velocity': v0, 'Min Impact Energy': emin, 'Max Impact Energy': emax,
        #                     'Phase Step': dphi, 'min Z': zmin, 'Space Step': dx, 'max z': zmax,
        #                     'Min Field Levels': flmin, 'Max Field Levels': flstep, 'Field Levels Step': flmax
        #                     }
        # save to json

        # with open('inputs.json', 'w') as f:
        #     json.dump(input_parameters, f)

    def run_multipac(self):
        self.get_parameters()
        self.run_field_solver()
        self.run_mpanalysis()

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

        # load geodata.n

        # load fieldparam
        # self.gtype = fieldparam(1)
        # self.freq = fieldparam(2)
        # self.epsr = fieldparam(3)
        # gridcons = fieldparam(4)
        s = 0
        # if nargin < 1:
        #     s = 1

        gtype = self.gtype
        freq = self.freq
        releps = self.epsr
        gridcons = self.d1

        cl = 0
        if releps > 1:
            print('Relative permittivity > 1?')

        geodata = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\geodata.n", sep='\s+', header=None).to_numpy()
        gridcons1 = geodata[0, 0]
          
        # generate the mesh
        job = self.mesh_cavity(geodata, gridcons1, self.epsr, s)
        
        # plot the mesh
        self.plot_mesh(0)
        ic("Here now done with plotting mesh")
        
        # find resonance solution
        cl = print('Computing the eigen values.')
        freqn = self.eigen(job, freq)
        if s > 0:
            print('Eigen values calculated.')
            print('                        ')

        err = abs(freqn - freq)/freq*100
        if err > 1:
          print('Warning: Error in eigen frequency more than 1#.')

        # compute and plot the fields
        self.calculate_fields(0, gridcons1, 0)

        self.plot_FEM_fields(0, gridcons1, s)
        # ---------------------------------------------------------------------

    def mesh_cavity(self, geodata, gridcons, epsr, s):
        n = len(geodata[:,1]) - 1
        gr = geodata[3:n, 0]
        gz = geodata[3:n, 1]
        gn = geodata[3:n, 2]
        n = len(gr)

        PEC = np.where(gn == 1)       # electric walls
        WIN = np.where(gn == 2)       # dielectric walls
        PMC = np.where(gn == 3)       # magnetic walls

        if len(WIN) > 0:
            print('Dielectric boundary found in a cavity.')
            return

        # specify the job for the eigenvalue solver
        if len(PMC) > 0:
            if s > 0:
                print('A cavity with magnetic walls.')
                job = 1
            elif len(PEC) > 0 & len(PMC) == 0:
                if s > 0:
                    cl = print('A cavity with electric walls.')
                    job = 0
            else:
                print('Invalid geometry.')
                return

        # save job job

        # First the nodes of the geometry
        nodes = [gz, gr]

        # Then the edges, the third column secifies the type of an edge as follows:
        #   1 : nxE = 0 on boundary
        #   3 : nxH = 0 on boundary
        #   0 : for streching

        edges = np.array([list(range(n)), list(range(1, n)).append(1), gn]).T
        print(edges)
        ne = len(edges[:, 0])

        # And finally the pacthes, first the edges are listed, last three values
        # give relative epsilon, relative mu and patch type, 0 = can strecth,
        # 1-n = regulars
        patches = np.append(np.array(0, ne), [epsr, 1, 1])
    
        esize = gridcons

        # start the mesh generator
        # !2dgen_bin
        return job

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
        ok1 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\mesh.mat")
        ok2 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\fieldparam")
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
            mesh = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\mesh.mat")
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

                osa = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
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
                self.ax.plot(coord[I[i], 0], coord[(I[i]-1, 1)], 'b*')

            I = np.where(boundary == 0)[0]
            self.ax.plot(coord[I, 0], coord[I, 1], 'w*')
            ala = ala/2
            self.ax.set_title('[MultiPac 2.0                    Mesh                  date ]')
            self.ax.set_xlabel('z axis [m]')
            self.ax.set_ylabel('r axis [m]')
            # self.plt.ax.set_axis('image')
            # self.ax.set_colormap('jet')
            self.fig.canvas.draw_idle()

    def eigen(self, job, freq):
        maara = 10
        job = 1
        offset = 0
        offset = 0
        # options.disp = 0

        cwd = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes'
        eigenCpath = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes\eigenC_bin.exe'
        if os.path.exists(eigenCpath):
            subprocess.call([eigenCpath, '-b'], cwd=cwd)
        # !eigenC_bin

        cinfo = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\cinfo.mat")
        if cinfo != 0:
            print('Eigenvalue solver failed.')

        o_eigen = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\o_eigen.mat")
        n = int(o_eigen['n'])
        ia = o_eigen['ia'].T[0] - 1  # mnatlab indices to python
        ja = o_eigen['ja'].T[0] - 1
        Aarvot = o_eigen['Aarvot'].T[0]
        Barvot = o_eigen['Barvot'].T[0]
        nf = o_eigen['nf']
        index = o_eigen['index']
        siirto = o_eigen['siirto']
        # load o_eigen
        ic(ia, ja, Aarvot, n)
        ic(np.shape(ia), np.shape(ja), np.shape(Aarvot))
        AA = sps.csr_matrix((Aarvot, (ia, ja)), shape=(n, n))
        BB = sps.csr_matrix((Barvot, (ia, ja)), shape=(n, n))

        d2, u2 = spsl.eigs(AA, M=BB, k=maara, sigma=0)
        d2 = np.absolute(d2)
        #
        mu0 = 4*np.pi*1e-7
        e0 = 8.85418782e-12
        k0 = 2*np.pi*freq*np.sqrt(mu0*e0)

        k = np.sqrt(d2)
        eigenvalues = k[0:maara]
        ic(eigenvalues)
        k2 = k

        val = np.min(abs(k0-k2))
        ind = np.argmin(abs(k0-k2))
        k = k2[ind]
        u = u2[:, ind]
        ic(k)

        # new frequency
        freqn = k/(2*np.pi*np.sqrt(mu0*e0))
        ic(freqn)
        param = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\param", sep='\s+', header=None).to_numpy().T[0]
        ic(param)
        # load param
        param[0] = freqn
        # save -ascii param param

        fieldparam = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\fieldparam", sep='\s+', header=None).to_numpy().T[0]
        # load fieldparam
        fieldparam[1] = freqn
        ic(freqn)
        # save -ascii fieldparam fieldparam

        # save kama0 index -v4
        # save kama0 u -append
        # save kama0 k -append
        # save kama0 u -v4 -append
        # save kama0 k -v4 -append
        # --------------------------
        return freqn
    
    def calculate_fields(self, a, gridcons1, b):
        nargin = 0
        if nargin < 3:
            ss = 1

        ok3 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\fieldparam")
        if ok3:
            fieldparam = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\fieldparam", sep='\s+', header=None).to_numpy().T[0]
            ic(fieldparam)
            gtype = fieldparam[0]
            ic(gtype)
            if nargin < 1:
                if gtype == 1:
                    s = 0
                else:
                    s = 1

        if not os.path.exists('s'):
            s = 0
        
        cl = 0
        ok1 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\mesh.mat")
        if not ok1:
            ic(['The mesh does not exists. Choose Mesh Generator in menu Run.'])

        if s == 0:
            ok2 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\kama0.mat")
            if ok2 == 0:
                ic(['The fields do not exist. Choose Field Solver in menu Run.'])

            elif s == 1:
                ok21 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\kama1.mat")
                ok22 = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\kama2.mat")
                ok2 = ok21*ok22
                if ok2 == 0:
                    ic(['The fields do not exist. Choose Field Solver in menu Run.'])
        
        if ok1*ok2*ok3 > 0:
            ok = os.path.exists(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodata.n")
            if ok:
                ic('Computing the fields.')
                #  error_message('                                  ')
        
                # define the grid constant
                if nargin < 2:
                    if gtype <= 2:
                        geodata = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodata.n", sep='\s+', header=None).to_numpy()
                        gridcons = geodata[0, 0]
                else:
                    geodata = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodatal.n", sep='\s+', header=None).to_numpy()
                    gridcons = geodata[0, 0]
        
        # load mesh
        ic(gtype)
        mesh = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\mesh.mat")
        coord = mesh['coord']
        etopol = mesh['etopol']
        alue = mesh['alue']
        tyyppi = mesh['tyyppi']
        boundary = mesh['boundary']
        edges = mesh['edges']

        if gtype <= 2:
            geodata = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodata.n", sep='\s+', header=None).to_numpy()
            n = len(geodata[:, 0])
            gr = geodata[3:n, 0]
            gz = geodata[3:n, 1]
        else:
            geodatal = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodatal.n", sep='\s+', header=None).to_numpy()
            nl = len(geodatal[:, 0]) - 1
            grl = geodatal[3:nl, 0]
            gzl = geodatal[3:nl, 1]
            geodataw = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodataw.n", sep='\s+', header=None).to_numpy()
            nw = len(geodataw[:, 0]) - 1
            grw = geodataw[3:nw, 0]
            gzw = geodataw[3:nw, 0]
            geodatar = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodatar.n", sep='\s+', header=None).to_numpy()
            nr = len(geodatar[:, 0]) - 1
            grr = geodatar[3:nr, 0]
            gzr = geodatar[3:nr, 1]
            gz = np.sort([gzl, gzw, gzr])
            gr = np.sort([grl, grw, grr])
        
        # generate field points
        z1 = min(gz)
        z2 = max(gz)
        r1 = min(gr)
        r2 = max(gr)
        
        zmaara = max(10, min(200, len(np.arange(z1, z2, gridcons))))
        if np.mod(zmaara, 2) == 0:
            zmaara = zmaara + 1

        rmaara = max(10, min(200, len(np.arange(r1, r2, gridcons))))
        I1 = np.arange(z1, z2 + (z2-z1)/(zmaara - 1), (z2-z1)/(zmaara - 1))  # (z2-z1)/(zmaara - 1) included so that end pint is included
        I2 = np.arange(r1, r2 + (z2-z1)/(zmaara - 1), (r2-r1)/(rmaara - 1))  # (z2-z1)/(zmaara - 1) included so that end pint is included
        
        zz, rr = np.meshgrid(I1, I2)
        ic(rr)
        ic(zz)
        m = len(I1)
        n = len(I2)
        z = zz.T.flatten()
        r = rr.T.flatten()

        alue = 0

        ax = np.where(r == 0)
        ic(ax)
        r[ax] = r[ax] + 1e-10
        ic(r, len(r))
        ind = np.arange(0, n)*np.ones((m, n))
        ic(ind)
        rind = ind.flatten()
        ic(rind)
        ind2 = (np.ones((n, m))*np.arange(0, m)).T
        ic(ind2)
        zind = ind2.flatten()
        ic(zind)
        
        zr = {'z': z, 'r': r, 'alue': alue, 'rind': rind, 'zind': zind}
        spio.savemat(zr)

        #     save zr r -v4 -append
        #     save zr alue -v4 -append
        #     save zr rind -v4 -append
        #     save zr zind -v4 -append
        
        # compute the fields at generated points
        if s == 0:
            # !copy /y kama0.mat kama.mat
            # !Multipac fields

            cwd = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes'
            multipacPath = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Multipac.exe'
            if os.path.exists(multipacPath):
                subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd)

            Er = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Er.mat")
            Ez = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Ez.mat")
            H = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\H.mat")

            fields = {'I1': I1, 'I2': I2, 'rr': rr, 'zz': zz, 'r': r, 'z': z, 'Er': Er, 'Ez': Ez, 'H': H}
            spio.savemat('fields', fields)
            self.save_fields(0, 1)
            self.normalize_u(0)
        else:
            job = 0
            # save job job -v4
            # !copy /y kama1.mat kama.mat
            # !Multipac fields
            cwd = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes'
            multipacPath = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Multipac.exe'
            if os.path.exists(multipacPath):
                subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd)

            Er = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Er.mat")
            Ez = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Ez.mat")
            H = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\H.mat")

            fields = {'I1': I1, 'I2': I2, 'rr': rr, 'zz': zz, 'r': r, 'z': z, 'Er': Er, 'Ez': Ez, 'H': H}
            spio.savemat('fields', fields)
            self.save_fields(1, job)
            self.normalize_u(1)

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

    def save_fields(self, wall, job):

        fieldfile1 = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\fieldfile1.txt", sep='\s+',
                               header=None).to_numpy()
        file = fieldfile1
        
        # compute the peak electric field on the boundary or the rf power
        # load job
        if wall == 0 and job == 1:
            E0 = self.peak_cavity_field()
        else:
            E0 = self.peak_coupler_field(wall)
        
        # normalize the fields
        my0 = 4e-7*np.pi
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
            spio.savemat('fieldfile1.n', fieldfile1)
        elif wall == 1:
            # save -ascii fieldfileE.n fieldfile1
            spio.savemat('fieldfileE.n', fieldfile1)
        elif wall == 2:
            # save -ascii fieldfileH.n fieldfile1
            spio.savemat('fieldfileH.n', fieldfile1)

    def normalize_u(self, wall):
        pass

    def peak_cavity_field(self):
        
        # load mesh and field solution
        # load mesh
        # load kama0
        mesh = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\mesh.mat")
        kama0 = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\kama0.mat")
        
        # compute the field on the boundary
        # load geodata.n
        geodata = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\geodata.n", sep='\s+',
                                 header=None).to_numpy()
        
        n = len(geodata[:, 0])
        ind = np.where(geodata[3:n, 2] == 1)
        ind = ind[2:len(ind)]
        r = geodata[4:n, 1]  # r = r(ind)-5e-4
        r = r[ind]-1.75e-4                   # move points inside
        z = geodata[3:n, 1]
        z = z[ind]
        
        rind = np.arange(1, len(r))
        zind = np.arange(1, len(z))
        
        alue = 0
        
        # save zr z -v4
        # save zr r -append
        # save zr alue -append
        # save zr rind -append
        # save zr zind -append

        zr = {'z': z, 'r': r, 'alue': alue, 'rind': rind, 'zind': zind}
        spio.savemat(zr)
        # save zr r -v4 -append
        # save zr alue -v4 -append
        # save zr rind -v4 -append
        # save zr zind -v4 -append
        
        # compute the field at generated points
        # !copy /y kama0.mat kama.mat
        # !Multipac fields

        cwd = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes'
        multipacPath = fr'D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Multipac.exe'
        if os.path.exists(multipacPath):
            subprocess.call([multipacPath, 'fields', '-b'], cwd=cwd)
        
        Er = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Er.mat")['Er']
        Ez = spio.loadmat(fr"D:\Dropbox\CavityDesignHub\utils\multipacting_codes\Ez.mat")['Ez']

        ee = np.sqrt(abs(Er)**2+abs(Ez)**2)
        E0 = max(ee)

        return ee

    def peak_coupler_field(self, wall):
        pass

    def plot_FEM_fields(self, a, gridcons1, b):
        pass

    def run_mpanalysis(self):
        pass

    def run_pause_resume_stop_routine(self):
        pass
        # if self.process_state == 'none':
        #     # change pause/resume icon to pause icon
        #     self.ui.pb_Pause_Resume.setIcon(self.pause_icon)
        #
        #     # disable pause/resume and cancel buttons
        #     self.ui.pb_Pause_Resume.setEnabled(False)
        #     self.ui.pb_Cancel.setEnabled(False)
        #
        #     # enable run button in case it was disabled
        #     self.ui.pb_Run.setEnabled(True)
        #
        # if self.process_state == "running":
        #     # enable run, pause/resume and cancel buttons
        #     self.ui.pb_Pause_Resume.setEnabled(True)
        #     self.ui.pb_Cancel.setEnabled(True)
        #     self.ui.pb_Run.setEnabled(False)
        #
        #     # change pause/resume icon to pause icon
        #     self.ui.pb_Pause_Resume.setIcon(self.pause_icon)
        #
        # if self.process_state == 'paused':
        #     # disable run button
        #     self.ui.pb_Run.setEnabled(False)
        #
        #     # change pause/resume button icon to resume icon
        #     self.ui.pb_Pause_Resume.setIcon(self.resume_icon)

    def pause(self):
        # self.log.info("Pausing...")
        for p in self.processes:
            p.suspend()
        self.log.info("Paused")

        self.process_state = 'paused'
        self.run_pause_resume_stop_routine()

    def resume(self):
        # self.log.info("Resuming...")
        for p in self.processes:
            p.resume()
        self.log.info("Resumed")

        self.process_state = 'running'
        self.run_pause_resume_stop_routine()

    def cancel(self):
        self.log.info("Terminating process...")
        # signal to progress bar
        self.show_progress_bar = False

        try:
            for p in self.processes:
                p.terminate()
        except:
            pass

        self.processes.clear()
        self.processes_id.clear()

        self.process_state = 'none'
        self.run_pause_resume_stop_routine()
        self.log.info("Process terminated.")

    def end_routine(self, proc_ids):
        print(proc_ids, type(proc_ids))
        for pid in proc_ids:
            try:
                p = psutil.Process(pid)
                while p.is_running():
                    pass

                print(fr"process {p} ended")
            except:
                pass

        self.cancel()

    # def update_progress_bar(self, val):
    #     self.progress_bar.setValue(val)
    #
    #     if val == 100 or not self.show_progress_bar:
    #         # reset progress bar
    #         self.progress_bar.setValue(0)
    #         self.progress_bar.hide()

    def get_geometric_parameters(self, code):
        self.shape_space = {}
        print_('Getting geometric parameters')
        if self.ui.cb_Shape_Entry_Mode.currentIndex() == 0:
            print_("Test worked")
            try:
                # self._shape_space = self.load_shape_space(shape_space_name)
                # print_(self._shape_space)

                # get selected keys
                self._selected_keys = self.ui.cb_Shape_Space_Keys.currentText()
                print("selected keys: ", self.ui.cb_Shape_Space_Keys.currentText())
                # print("Selected keys: ", self._selected_keys, type(self._selected_keys[0]))

                # check keys of shape space if results already exist
                toall = None
                for key, val in self._shape_space.items():
                    # process for only keys selected in combobox
                    if self.ui.cb_Shape_Space_Keys.currentText() == "" or self.ui.cb_Shape_Space_Keys.currentText() == "All":
                        pass
                    else:
                        if key not in self._selected_keys:
                            continue

                    if not toall:
                        ans = self.prompt(code, key)
                        if ans == 'Yes':
                            self.shape_space[key] = val

                        if ans == 'No':
                            continue

                        if ans == 'YesToAll':
                            self.shape_space[key] = val
                            toall = 'YesToAll'

                        if ans == 'NoToAll':
                            toall = 'NoToAll'

                        if ans == "Does not exist":
                            self.shape_space[key] = val
                            toall = None
                    else:
                        if toall == 'YesToAll':
                            self.shape_space[key] = val
                        else:
                            path = f'{self.main_control.projectDir}/SimulationData/{code}/Cavity{key}'
                            if os.path.exists(path):
                                continue
                            else:
                                self.shape_space[key] = val

                # print_(self.shape_space)
                return self.shape_space
            except Exception as e:
                print_(f"File not found, check path:: {e}")
        else:
            if self.ui.cb_Inner_Cell.checkState() == 2:
                # Middle Ellipse data
                A_i_space = self.text_to_list(self.ui.le_A_i.text())
                B_i_space = self.text_to_list(self.ui.le_B_i.text())
                a_i_space = self.text_to_list(self.ui.le_a_i.text())
                b_i_space = self.text_to_list(self.ui.le_b_i.text())
                Ri_i_space = self.text_to_list(self.ui.le_Ri_i.text())
                L_i_space = self.text_to_list(self.ui.le_L_i.text())
                Req_i_space = self.text_to_list(self.ui.le_Req_i.text())
                alpha_i_space = self.text_to_list(self.ui.le_Alpha.text())

                inner_cell_space = [A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space,
                                    alpha_i_space]
            else:
                inner_cell_space = [[0], [0], [0], [0], [0], [0], [0], [0]]

            if self.ui.cb_Outer_Cell_L.checkState() == 2:
                # Middle Ellipse data
                A_ol_space = self.text_to_list(self.ui.le_A_ol.text())
                B_ol_space = self.text_to_list(self.ui.le_B_ol.text())
                a_ol_space = self.text_to_list(self.ui.le_a_ol.text())
                b_ol_space = self.text_to_list(self.ui.le_b_ol.text())
                Ri_ol_space = self.text_to_list(self.ui.le_Ri_ol.text())
                L_ol_space = self.text_to_list(self.ui.le_L_ol.text())
                Req_ol_space = self.text_to_list(self.ui.le_Req_ol.text())
                alpha_ol_space = self.text_to_list(self.ui.le_Alpha_ol.text())

                outer_cell_L_space = [A_ol_space, B_ol_space, a_ol_space, b_ol_space, Ri_ol_space, L_ol_space,
                                      Req_ol_space, alpha_ol_space]
            else:
                outer_cell_L_space = inner_cell_space

            if self.ui.cb_Outer_Cell_R.checkState() == 2:
                # Middle Ellipse data
                A_or_space = self.text_to_list(self.ui.le_A_or.text())
                B_or_space = self.text_to_list(self.ui.le_B_or.text())
                a_or_space = self.text_to_list(self.ui.le_a_or.text())
                b_or_space = self.text_to_list(self.ui.le_b_or.text())
                Ri_or_space = self.text_to_list(self.ui.le_Ri_or.text())
                L_or_space = self.text_to_list(self.ui.le_L_or.text())
                Req_or_space = self.text_to_list(self.ui.le_Req_or.text())
                alpha_or_space = self.text_to_list(self.ui.le_Alpha_or.text())

                outer_cell_R_space = [A_or_space, B_or_space, a_or_space, b_or_space, Ri_or_space, L_or_space,
                                      Req_or_space, alpha_or_space]
            else:
                outer_cell_R_space = inner_cell_space
            count = 0
            for A_i in inner_cell_space[0]:
                for B_i in inner_cell_space[1]:
                    for a_i in inner_cell_space[2]:
                        for b_i in inner_cell_space[3]:
                            for Ri_i in inner_cell_space[4]:
                                for L_i in inner_cell_space[5]:
                                    for Req_i in inner_cell_space[6]:
                                        if outer_cell_L_space == inner_cell_space:
                                            inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                            outer_cell_L = inner_cell

                                            if self.ui.cb_LBP.checkState() == 2 and self.ui.cb_RBP.checkState() == 2:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                           'OC_R': outer_cell_L, 'BP': 'both',
                                                                           'FREQ': None}
                                            elif self.ui.cb_LBP.checkState() == 2 and self.ui.cb_RBP.checkState() == 0:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                           'OC_R': outer_cell_L, 'BP': 'left',
                                                                           'FREQ': None}
                                            elif self.ui.cb_LBP.checkState() == 0 and self.ui.cb_RBP.checkState() == 2:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                           'OC_R': outer_cell_L, 'BP': 'right',
                                                                           'FREQ': None}
                                            else:
                                                self.shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                           'OC_R': outer_cell_L, 'BP': 'none',
                                                                           'FREQ': None}

                                            count += 1
                                        else:
                                            for A_ol in outer_cell_L_space[0]:
                                                for B_ol in outer_cell_L_space[1]:
                                                    for a_ol in outer_cell_L_space[2]:
                                                        for b_ol in outer_cell_L_space[3]:
                                                            for Ri_ol in outer_cell_L_space[4]:
                                                                for L_ol in outer_cell_L_space[5]:
                                                                    # for Req_ol in outer_cell_L_space[6]:
                                                                    if outer_cell_L_space == outer_cell_R_space:
                                                                        inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i,
                                                                                      Req_i, 0]
                                                                        outer_cell_L = [A_ol, B_ol, a_ol, b_ol, Ri_ol,
                                                                                        L_ol, Req_i, 0]
                                                                        outer_cell_R = outer_cell_L
                                                                        if self.ui.cb_LBP.checkState() == 2 and self.ui.cb_RBP.checkState() == 0:
                                                                            self.shape_space[count] = {'IC': inner_cell,
                                                                                                       'OC': outer_cell_L,
                                                                                                       'OC_R': outer_cell_R,
                                                                                                       'BP': 'left',
                                                                                                       'FREQ': None}
                                                                        elif self.ui.cb_LBP.checkState() == 0 and self.ui.cb_RBP.checkState() == 2:
                                                                            self.shape_space[count] = {'IC': inner_cell,
                                                                                                       'OC': outer_cell_L,
                                                                                                       'OC_R': outer_cell_R,
                                                                                                       'BP': 'right',
                                                                                                       'FREQ': None}
                                                                        elif self.ui.cb_LBP.checkState() == 2 and self.ui.cb_RBP.checkState() == 2:
                                                                            self.shape_space[count] = {'IC': inner_cell,
                                                                                                       'OC': outer_cell_L,
                                                                                                       'OC_R': outer_cell_R,
                                                                                                       'BP': 'both',
                                                                                                       'FREQ': None}
                                                                        else:
                                                                            self.shape_space[count] = {'IC': inner_cell,
                                                                                                       'OC': outer_cell_L,
                                                                                                       'OC_R': outer_cell_R,
                                                                                                       'BP': 'none',
                                                                                                       'FREQ': None}

                                                                        count += 1
                                                                    else:
                                                                        for A_or in outer_cell_R_space[0]:
                                                                            for B_or in outer_cell_R_space[1]:
                                                                                for a_or in outer_cell_R_space[2]:
                                                                                    for b_or in outer_cell_R_space[3]:
                                                                                        for Ri_or in outer_cell_R_space[
                                                                                            4]:
                                                                                            for L_or in \
                                                                                            outer_cell_R_space[5]:
                                                                                                # for Req_or in outer_cell_R_space[6]:
                                                                                                inner_cell = [A_i, B_i,
                                                                                                              a_i, b_i,
                                                                                                              Ri_i, L_i,
                                                                                                              Req_i, 0]
                                                                                                outer_cell_L = [A_ol,
                                                                                                                B_ol,
                                                                                                                a_ol,
                                                                                                                b_ol,
                                                                                                                Ri_ol,
                                                                                                                L_ol,
                                                                                                                Req_i,
                                                                                                                0]
                                                                                                outer_cell_R = [A_or,
                                                                                                                B_or,
                                                                                                                a_or,
                                                                                                                b_or,
                                                                                                                Ri_or,
                                                                                                                L_or,
                                                                                                                Req_i,
                                                                                                                0]
                                                                                                if self.ui.cb_LBP.checkState() == 2 and self.ui.cb_RBP.checkState() == 0:
                                                                                                    self.shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'left',
                                                                                                        'FREQ': None}
                                                                                                elif self.ui.cb_LBP.checkState() == 0 and self.ui.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'right',
                                                                                                        'FREQ': None}
                                                                                                elif self.ui.cb_LBP.checkState() == 2 and self.ui.cb_RBP.checkState() == 2:
                                                                                                    self.shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'both',
                                                                                                        'FREQ': None}
                                                                                                else:
                                                                                                    self.shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'none',
                                                                                                        'FREQ': None}

                                                                                                count += 1
            return self.shape_space

    def prompt(self, code, fid):
        path = fr'{self.main_control.projectDir}\SimulationData\SLANS\Cavity{fid}'
        # print(path)
        # path = os.path.join(path, fr"{}\{code}\Cavity{fid}")
        if os.path.exists(path):
            print_("Simulation data already exists. Do you want to overwrite it?")
            msg = QMessageBox()
            msg.setWindowTitle("Folder Exist")
            msg.setText("Simulation data already exists. Do you want to overwrite it?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.YesToAll | QMessageBox.Yes | QMessageBox.No | QMessageBox.NoToAll)
            msg.setDefaultButton(QMessageBox.Yes)

            msg.buttonClicked.connect(self.button_clicked)
            # print_(f'msg: {msg.Yes}')

            x = msg.exec_()

            print_("got here")
            if x == msg.YesToAll:
                return 'YesToAll'
            if x == msg.Yes:
                return 'Yes'
            if x == msg.No:
                return 'No'
            if x == msg.NoToAll:
                return 'NoToAll'
        else:
            return "Does not exist"

    def open_file(self, le, cb):
        # clear combobox
        self.ui.cb_Shape_Space_Keys.clear()
        self.ui.cb_Shape_Space_Keys.addItem('All')
        # self._selected_keys.clear()

        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Json Files (*.json)")
        try:
            le.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            # populate checkboxes with key
            for col in dd.keys():
                cb.addItem(fr'{col}')

            self._shape_space = dd

        except Exception as e:
            print('Failed to open file:: ', e)

    def animate_width(self, cb, widget, min_width, standard, enable, reverse=False):
        if enable:
            #### GET Height
            width = widget.width()
            #### SET MAX WIDTH

            if not reverse:
                if cb.checkState() != 2:
                    widthCollapsed = min_width
                    widget.setMinimumWidth(0)
                else:
                    widthCollapsed = standard
                    # widget.setMinimumWidth(standard)
            else:
                if cb.checkState() == 2:
                    widthCollapsed = min_width
                    widget.setMinimumWidth(0)
                else:
                    widthCollapsed = standard
                    # widget.setMinimumWidth(standard)

            #### ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumWidth")
            self.animation.setDuration(200)
            self.animation.setStartValue(width)
            self.animation.setEndValue(widthCollapsed)
            self.animation.start()

    def animate_height(self, cb, widget, min_height, standard, enable):
        if enable:
            #### GET WIDTH
            height = widget.width()
            #### SET MAX WIDTH
            if cb.checkState() != 2:
                heightCollapsed = min_height
                widget.setMinimumHeight(0)
            else:
                heightCollapsed = standard
                # widget.setMinimumWidth(standard)

            #### ANIMATION
            self.animation = QPropertyAnimation(widget, b"maximumHeight")
            self.animation.setDuration(200)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.start()

    # def exe_control(self):
    #     # Slans
    #     self.ui.pb_Genmsh.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\genmesh2.exe"))
    #     self.ui.pb_Sl.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\Sl.exe"))
    #     self.ui.pb_Slansc.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slansc.exe"))
    #     self.ui.pb_Slansm.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slansm.exe"))
    #     self.ui.pb_Slanss.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slanss.exe"))
    #     self.ui.pb_Slansre.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\slansre.exe"))
    #     self.ui.pb_MTFView.clicked.connect(lambda: self.run_slans_exe("SLANS_exe\Mtfview\mtfview.exe"))

    def run_slans_exe(self, path, filename=None):
        path = fr"{self.main_control.parentDir}\em_codes\{path}"
        t = Thread(target=subprocess.call, args=(path,))
        t.start()

    def draw_shape_from_shape_space(self):
        colors = [[48, 162, 218, 255], [252, 79, 48, 255], [229, 174, 56, 255], [109, 144, 79, 255],
                  [139, 139, 139, 255]]
        ci = 0

        # # remove existing cells
        # self.graphicsView.removeCells()
        # for key in self._shape_space.keys():
        #     if key in self.ui.cb_Shape_Space_Keys.currentText():
        #         IC = self._shape_space[key]["IC"]
        #         OC = self._shape_space[key]["OC"]
        #         BP = self._shape_space[key]["BP"]
        #         self.graphicsView.drawCells(IC, OC, BP, QColor(colors[ci][0], colors[ci][1], colors[ci][2], colors[ci][3]))
        #
        #         ci += 1

    # def ui_effects(self):
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Settings.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Inner_Cell.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Outer_Cell_L.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Outer_Cell_R.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Expansion.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Simulation_Controls.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Show_Cavity.setGraphicsEffect(shadow)
    #
    #     shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=5, yOffset=5)
    #     shadow.setColor(QColor(0, 0, 0, 77))
    #     self.ui.w_Load_Manual.setGraphicsEffect(shadow)

    def serialize(self, state_dict):
        # update state file
        state_dict["Eigen_Shape_Entry_Mode"] = self.ui.cb_Shape_Entry_Mode.currentIndex()
        state_dict["Eigen_Shape Space"] = self.ui.le_Shape_Space.text()
        state_dict["Eigen_Mid_Cell_CB"] = self.ui.cb_Inner_Cell.checkState()
        state_dict["Eigen_Left_Cell_CB"] = self.ui.cb_Outer_Cell_L.checkState()
        state_dict["Eigen_Right_Cell_CB"] = self.ui.cb_Outer_Cell_R.checkState()
        state_dict["Eigen_Expansion_CB"] = self.ui.cb_Expansion.checkState()
        state_dict["Eigen_LBP_CB"] = self.ui.cb_LBP.checkState()
        state_dict["Eigen_RBP_CB"] = self.ui.cb_RBP.checkState()

        # cell parameters
        state_dict["Eigen_A_i"] = self.ui.le_A_i.text()
        state_dict["Eigen_B_i"] = self.ui.le_B_i.text()
        state_dict["Eigen_a_i"] = self.ui.le_a_i.text()
        state_dict["Eigen_b_i"] = self.ui.le_b_i.text()
        state_dict["Eigen_Ri_i"] = self.ui.le_Ri_i.text()
        state_dict["Eigen_L_i"] = self.ui.le_L_i.text()
        state_dict["Eigen_Req_i"] = self.ui.le_Req_i.text()
        state_dict["Eigen_Alpha_i"] = self.ui.le_Alpha.text()

        state_dict["Eigen_A_ol"] = self.ui.le_A_ol.text()
        state_dict["Eigen_B_ol"] = self.ui.le_B_ol.text()
        state_dict["Eigen_a_ol"] = self.ui.le_a_ol.text()
        state_dict["Eigen_b_ol"] = self.ui.le_b_ol.text()
        state_dict["Eigen_Ri_ol"] = self.ui.le_Ri_ol.text()
        state_dict["Eigen_L_ol"] = self.ui.le_L_ol.text()
        state_dict["Eigen_Req_ol"] = self.ui.le_Req_ol.text()
        state_dict["Eigen_Alpha_ol"] = self.ui.le_Alpha_ol.text()

        state_dict["Eigen_A_or"] = self.ui.le_A_or.text()
        state_dict["Eigen_B_or"] = self.ui.le_B_or.text()
        state_dict["Eigen_a_or"] = self.ui.le_a_or.text()
        state_dict["Eigen_b_or"] = self.ui.le_b_or.text()
        state_dict["Eigen_Ri_or"] = self.ui.le_Ri_or.text()
        state_dict["Eigen_L_or"] = self.ui.le_L_or.text()
        state_dict["Eigen_Req_or"] = self.ui.le_Req_or.text()
        state_dict["Eigen_Alpha_or"] = self.ui.le_Alpha_or.text()

        # settings
        # state_dict["Eigen_N_Cells"] = self.ui.sb_N_Cells.value()
        # state_dict["Eigen_N_Modules"] = self.ui.sb_N_Modules.value()
        # state_dict["Eigen_Polarization"] = self.ui.cb_Polarization_SLANS.currentIndex()
        #
        # state_dict["Eigen_Freq_Shift"] = self.ui.le_Freq_Shift.text()
        # state_dict["Eigen_No_Of_Modes"] = self.ui.le_No_Of_Modes.text()
        # state_dict["Eigen_LBC"] = self.ui.cb_LBC.currentIndex()
        # state_dict["Eigen_RBC"] = self.ui.cb_RBC.currentIndex()
        # state_dict["Eigen_No_Of_Processors"] = self.ui.sb_No_Of_Processors_SLANS.value()

    def deserialize(self, state_dict):
        # update state file
        self.ui.cb_Shape_Entry_Mode.setCurrentIndex(state_dict["Eigen_Shape_Entry_Mode"])
        self.ui.le_Shape_Space.setText(state_dict["Eigen_Shape Space"])
        self.ui.cb_Inner_Cell.setCheckState(state_dict["Eigen_Mid_Cell_CB"])
        self.ui.cb_Outer_Cell_L.setCheckState(state_dict["Eigen_Left_Cell_CB"])
        self.ui.cb_Outer_Cell_R.setCheckState(state_dict["Eigen_Right_Cell_CB"])
        self.ui.cb_Expansion.setCheckState(state_dict["Eigen_Expansion_CB"])
        self.ui.cb_LBP.setCheckState(state_dict["Eigen_LBP_CB"])
        self.ui.cb_RBP.setCheckState(state_dict["Eigen_RBP_CB"])

        # cell parameters
        self.ui.le_A_i.setText(state_dict["Eigen_A_i"])
        self.ui.le_B_i.setText(state_dict["Eigen_B_i"])
        self.ui.le_a_i.setText(state_dict["Eigen_a_i"])
        self.ui.le_b_i.setText(state_dict["Eigen_b_i"])
        self.ui.le_Ri_i.setText(state_dict["Eigen_Ri_i"])
        self.ui.le_L_i.setText(state_dict["Eigen_L_i"])
        self.ui.le_Req_i.setText(state_dict["Eigen_Req_i"])
        self.ui.le_Alpha.setText(state_dict["Eigen_Alpha_i"])

        self.ui.le_A_ol.setText(state_dict["Eigen_A_ol"])
        self.ui.le_B_ol.setText(state_dict["Eigen_B_ol"])
        self.ui.le_a_ol.setText(state_dict["Eigen_a_ol"])
        self.ui.le_b_ol.setText(state_dict["Eigen_b_ol"])
        self.ui.le_Ri_ol.setText(state_dict["Eigen_Ri_ol"])
        self.ui.le_L_ol.setText(state_dict["Eigen_L_ol"])
        self.ui.le_Req_ol.setText(state_dict["Eigen_Req_ol"])
        self.ui.le_Alpha_ol.setText(state_dict["Eigen_Alpha_ol"])

        self.ui.le_A_or.setText(state_dict["Eigen_A_or"])
        self.ui.le_B_or.setText(state_dict["Eigen_B_or"])
        self.ui.le_a_or.setText(state_dict["Eigen_a_or"])
        self.ui.le_b_or.setText(state_dict["Eigen_b_or"])
        self.ui.le_Ri_or.setText(state_dict["Eigen_Ri_or"])
        self.ui.le_L_or.setText(state_dict["Eigen_L_or"])
        self.ui.le_Req_or.setText(state_dict["Eigen_Req_or"])
        self.ui.le_Alpha_or.setText(state_dict["Eigen_Alpha_or"])

        # # settings
        # self.ui.sb_N_Cells.setValue(state_dict["Eigen_N_Cells"])
        # self.ui.sb_N_Modules.setValue(state_dict["Eigen_N_Modules"])
        # self.ui.cb_Polarization_SLANS.setCurrentIndex(state_dict["Eigen_Polarization"])
        # self.ui.le_Freq_Shift.setText(state_dict["Eigen_Freq_Shift"])
        # self.ui.le_No_Of_Modes.setText(state_dict["Eigen_No_Of_Modes"])
        # self.ui.cb_LBC.setCurrentIndex(state_dict["Eigen_LBC"])
        # self.ui.cb_RBC.setCurrentIndex(state_dict["Eigen_RBC"])
        # self.ui.sb_No_Of_Processors_SLANS.setValue(state_dict["Eigen_No_Of_Processors"])

    @staticmethod
    def load_shape_space(filename):
        fr = FileReader()
        dir = filename

        # check if extension is included
        if dir.split('.')[-1] != 'json':
            dir = f'{dir}.json'

        df = fr.json_reader(dir)
        # print_(df)

        return df.to_dict()

    @staticmethod
    def show_hide_(wid1, wid2):
        print('here')
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()

    @staticmethod
    def button_clicked(i):
        return i.text()

    @staticmethod
    def text_to_list(txt):
        if "range" in txt:
            txt = txt.replace('range', '')
            l = ast.literal_eval(txt)
            return range(l[0], l[1], l[2])
        elif 'linspace' in txt:
            l = eval(f'np.{txt}')
            print(l)
            return l
        else:
            l = ast.literal_eval(txt)
            if isinstance(l, int) or isinstance(l, float):
                return [l]
            else:
                return list(l)

    # @staticmethod
    def run_sequential(self, n_cells, n_modules, processor_shape_space, n_modes, f_shift, bc, parentDir, projectDir,
                       progress_list, sub_dir=''):
        progress = 0
        # get length of processor
        print("Now it is heredddddddddd")
        total_no_of_shapes = len(list(processor_shape_space.keys()))
        print("Now it is heredddddddddd")

        for key, shape in processor_shape_space.items():
            # # create folders for all keys
            slans_geom.createFolder(key, projectDir, subdir=sub_dir)

            print("Now it is here")
            self.write_cst_paramters(key, shape['IC'], shape['OC'], projectDir)

            print("Now it is heresdfsdf")
            # run slans code
            start_time = time.time()
            try:
                print("Now it is here")
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
            except:
                print("Now it is here")
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

            print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')

            # update progress
            progress_list.append((progress + 1) / total_no_of_shapes)

    def update_alpha(self):
        A_i_space = self.text_to_list(self.ui.le_A_i.text())[0]
        B_i_space = self.text_to_list(self.ui.le_B_i.text())[0]
        a_i_space = self.text_to_list(self.ui.le_a_i.text())[0]
        b_i_space = self.text_to_list(self.ui.le_b_i.text())[0]
        Ri_i_space = self.text_to_list(self.ui.le_Ri_i.text())[0]
        L_i_space = self.text_to_list(self.ui.le_L_i.text())[0]
        Req_i_space = self.text_to_list(self.ui.le_Req_i.text())[0]
        print(A_i_space)
        try:
            alpha_i_space = self._calculate_alpha(A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space,
                                                  Req_i_space, 0)
            self.ui.le_Alpha.setText(f"{round(alpha_i_space, 2)}")
        except:
            pass

    def _calculate_alpha(self, A, B, a, b, Ri, L, Req, L_bp):

        data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
                [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(self._ellipse_tangent,
                                np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                                args=data)
        m = (y2 - y1) / (x2 - x1)
        alpha = 180 - np.arctan(m) * 180 / np.pi

        return alpha

    def _ellipse_tangent(self, z, *data):
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return f1, f2, f3, f4

    def write_cst_paramters(self, key, ic, oc, projectDir):
        # print("Writing parameters to file")
        path = fr'{projectDir}/SimulationData/SLANS/Cavity{key}/{key}.txt'

        # print(path)
        with open(path, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', 'Aeq_e', 'Beq_e', 'ai_e', 'bi_e', 'Ri_e',
                         'L_e', 'Req', 'alpha_e', 'key']

            value_list = [ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], ic[6], ic[7],
                          oc[0], oc[1], oc[2], oc[3], oc[4], oc[5], oc[6], oc[7], key]

            for i in range(len(name_list)):
                f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')


class ProgressMonitor(QThread):
    sig = QtCore.pyqtSignal(int)

    def __init__(self, eig_control, projectDir):
        super(QThread, self).__init__()
        self.eig_control = eig_control
        self.proc_ids = eig_control.processes_id
        self.progress_bar = eig_control.progress_bar
        self.projectDir = projectDir

    def run(self):
        ic(self.proc_ids)
        self.progress_monitor(self.proc_ids)

    def progress_monitor(self, proc_ids):
        proc_count = len(proc_ids)
        while self.eig_control.show_progress_bar:
            progress = sum(self.eig_control.progress_list) / proc_count
            self.sig.emit(round(progress * 100, 10))


class EndRoutine(QThread):
    def __init__(self, control, projectDir):
        super(QThread, self).__init__()
        self.control = control
        self.proc_ids = control.processes_id
        self.projectDir = projectDir

    def run(self):
        self.end_routine()

    def end_routine(self):
        for pid in self.proc_ids:
            try:
                p = psutil.Process(pid)
                while p.is_running():
                    pass

                print(fr"process {p} ended")
            except:
                pass

        self.control.cancel()
