import ast
import json
import os

# # convert matlab in folders to python
#
# #"python matlab2python/matlab2python.py SSC_single_model.m -o SSC_single_model.py"
# import subprocess
#
# dirr = fr'{os.getcwd()}\S_for_Sosoho'
# print(dirr)
#
# if os.path.exists(dirr):
#     print("yeah")
#
# for path, subdirs, files in os.walk(dirr):
#     for name in files:
#         if name.split(".")[-1] == 'm':
#             file_dirname = os.path.join(path, name).split('.')[0]
#             print(file_dirname)
#             command = fr"python matlab2python/matlab2python.py {file_dirname}.m -o {file_dirname}.py"
#             p = subprocess.run(command)
import shutil
import subprocess

import sys
import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg
from icecream import ic
import scipy as sp
# import pyvista as pv
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from utils.file_reader import FileReader

fr = FileReader()

# print(sys.path)

# filename = "D:\SSC_for_Sosoho\Modes_vtk\\1_1_x.vtk"
# mesh = pv.read(filename)
# cpos = mesh.plot()
# import numba as nb
import pandas as pd

# fig = plt.figure(figsize=(12, 6))
# ax = fig.add_subplot(projection='3d')

# @nb.njit(parallel=True)
# def isin(a, b):
#     out = np.empty(a.shape[0], dtype=nb.boolean)
#     b = set(b)
#     for i in nb.prange(a.shape[0]):
#         if a[i] in b:
#             out[i] = False
#         else:
#             out[i] = True
#     return out
#
#
# def remove_column_elements_csr_matrix(mat, indx):
#     # t1 = time.time()
#     mask = isin(mat.indices, indx)
#     # t2 = time.time()
#     # ic(t2-t1)
#     mat.data = mat.data[mask]
#     # t3 = time.time()
#     # ic(t3-t2)
#     mat.indices = mat.indices[mask]
#     # t4 = time.time()
#     # ic(t4-t3)
#     parts = [np.sum(mask[mat.indptr[i]:mat.indptr[i + 1]]) for i in range(mat.indptr.shape[0] - 1)]
#     # t5 = time.time()
#     # ic(t5-t4)
#     mat.indptr[1:] = np.cumsum(parts)
#     # t6 = time.time()
#     # ic(t6-t5)
#
#     return mat
#
#
# def ispm_iteration(A, s, num_simulations: int, q_k=None):
#     # shift A
#     A = A - s * np.identity(np.shape(A)[0])
#
#     if q_k is None:
#         q_k = np.random.rand(A.shape[1])
#
#     for _ in range(num_simulations):
#         # calculate z by solving Az = q_k
#         z = sps.linalg.spsolve(A, q_k)
#
#         # calculate the norm of z
#         z_norm = np.linalg.norm(z)
#
#         # re normalize the vector
#         q_k = z / z_norm
#
#     # to calculate the corresponding eigenvalue
#     # solve LLT y = q_k for y using forward and backward substitutions
#     y = sps.linalg.spsolve(A, q_k)
#
#     # solve q_k^T y
#     v = q_k.T @ y
#
#     # get eigenvalue by inverting v
#     l = 1 / v + s
#
#     return l, q_k
#
#
# def arnoldi_iteration(A, b, n: int):
#     """Computes a basis of the (n + 1)-Krylov subspace of A: the space
#     spanned by {b, Ab, ..., A^n b}.
#
#     Arguments
#       A: m × m array
#       b: initial vector (length m)
#       n: dimension of Krylov subspace, must be >= 1
#
#     Returns
#       Q: m x (n + 1) array, the columns are an orthonormal basis of the
#         Krylov subspace.
#       h: (n + 1) x n array, A on basis Q. It is upper Hessenberg.
#     """
#     eps = 1e-12
#     h = np.zeros((n + 1, n), dtype=np.complex_)
#     Q = np.zeros((A.shape[0], n + 1), dtype=np.complex_)
#     # Normalize the input vector
#     Q[:, 0] = b / np.linalg.norm(b, 2)  # Use it as the first Krylov vector
#     for k in range(1, n + 1):
#         v = np.dot(A, Q[:, k - 1])  # Generate a new candidate vector
#         for j in range(k):  # Subtract the projections on previous vectors
#             h[j, k - 1] = np.dot(Q[:, j].T, v)
#             v = v - h[j, k - 1] * Q[:, j]
#         h[k, k - 1] = np.linalg.norm(v, 2)
#         if h[k, k - 1] > eps:  # Add the produced vector to the list, unless
#             Q[:, k] = v / h[k, k - 1]
#         else:  # If that happens, stop iterating.
#             return Q, h
#     return Q, h, k


# [-1.1279152 +0.j         -0.04781313+0.j          5.58786416+2.42104632j 5.58786416-2.42104632j]

# A = np.array([[1, 3, 4, -4, -7], [0, 3, 5, -2, 0], [1, 2, 1, 0, -6], [0, 3.0, 2, 5, -4], [5, -5, 3, 2, 0]], dtype=np.complex_)
# A = np.random.rand(400, 400)
# # np.savetxt('test_matrix.txt', A)
# A = np.loadtxt('test_matrix.txt')
#
# As = sps.csr_matrix(A)
# n = 20
# d1, v1 = np.linalg.eig(A)
# d, v = sps.linalg.eigs(As, k=n, sigma=1)
# # l, q_k = ispm_iteration(As, s=4-2j, num_simulations=500)
# ic(np.sort_complex(d1))
# # print(d)
# # print(l)
#
# Q, h, k = arnoldi_iteration(A, np.random.rand(A.shape[1]), n=n)
# print()
# ic(np.sort_complex(np.linalg.eigvals(h[0:n,0:n])))

# import matplotlib.pyplot as plt
# data1 = pd.read_excel(r"D:\CST Studio\Multipacting\E_field_equator_offset_Efield.xlsx", "E3794")
# data2 = pd.read_excel(r"D:\CST Studio\Multipacting\E_field_equator_offset_Efield.xlsx", "E3794_HC")
# data3 = pd.read_excel(r"D:\CST Studio\Multipacting\E_field_equator_offset_Efield.xlsx", "E3794_HC_W_PORTS")
#
# plt.scatter(data1["R"], data1["E3794"], marker="+", label="E3794")
# plt.scatter(data2["R"], data2["E3794_HC"], marker="+", label="E3794_HC")
# plt.scatter(data3["R"], data3["E3794_HC_W_PORTS"], marker="+", label="E3794_HC_W_PORTS")
# plt.xlim([-50, 50])
# plt.ylim([1.7e6, 1.725e6])
# plt.axvline(0)
# plt.legend()
# plt.show()

# import scipy.io as spio
# import os
#
#
# print(os.name)
# matfile = "D:\Dropbox\multipacting\MPGUI21\H.mat"
# matdata = spio.loadmat(matfile)
# spio.savemat(r"D:\Dropbox\multipacting\MPGUI21\multipac\texst.mat", matdata)

# (b'MATLAB 5.0 MAT-file Platform: nt, Created on: Tue Apr 26 20:05:47 2022', 0, 256, b'IM') <class 'numpy.ndarray'> <class 'numpy.ndarray'>
# b'MATLAB 5.0 MAT-file Platform: nt, Created on: Tue Apr 26 20:05:47 2022\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01IM'
#
#
# def create_pseudo_shape_space(var_list, lock_list):
#     A, B, a, b, Ri, L, Req = var_list
#     A_LOCKED, B_LOCKED, a_LOCKED, b_LOCKED, Ri_LOCKED, L_LOCKED, Req_LOCKED = lock_list
#
#     if Req[0] == 0:
#         AA, BB, aa, bb, RiRi, LL = np.meshgrid(A, B, a, b, Ri, L)
#
#         space = []
#         for j in range(len(A)):  # for some reason the first dimension of A is the dimension of B
#             for i in range(len(B)):
#                 for k in range(len(a)):
#                     for m in range(len(b)):
#                         for n in range(len(Ri)):
#                             for o in range(len(L)):
#                                 space.append([AA[i][j][k][m][n][o], BB[i][j][k][m][n][o], aa[i][j][k][m][n][o],
#                                               bb[i][j][k][m][n][o], RiRi[i][j][k][m][n][o], LL[i][j][k][m][n][o],
#                                               RiRi[i][j][k][m][n][o] + BB[i][j][k][m][n][o] + bb[i][j][k][m][n][o], 0])
#     else:
#         AA, BB, aa, bb, RiRi, ReqReq = np.meshgrid(A, B, a, b, Ri, Req)
#
#         space = []
#         for j in range(len(A)):  # for some reason the first dimension of A is the dimension of B
#             for i in range(len(B)):
#                 for k in range(len(a)):
#                     for m in range(len(b)):
#                         for n in range(len(Ri)):
#                             for o in range(len(L)):
#                                 space.append([AA[i][j][k][m][n][o], BB[i][j][k][m][n][o], aa[i][j][k][m][n][o],
#                                               bb[i][j][k][m][n][o], RiRi[i][j][k][m][n][o],
#                                               AA[i][j][k][m][n][o] + aa[i][j][k][m][n][o],
#                                               ReqReq[i][j][k][m][n][o], 0])
#
#     space = np.array(space)
#     # print(space)
#     # create list of locked variables
#     LOCKED_LIST = np.array([A_LOCKED, B_LOCKED, a_LOCKED, b_LOCKED, Ri_LOCKED, L_LOCKED])
#     count = np.count_nonzero(LOCKED_LIST)
#     print(count)
#     if count >= 2:
#         # for ll in LOCKED_LIST:
#         # check if length of all locked lists are same
#         dummy_list = []
#         for i in LOCKED_LIST.nonzero()[0]:
#             dummy_list.append(len(var_list[i]))
#             lock_len = len(var_list[i])
#         dummy_list = np.array(dummy_list)
#
#         if np.all(dummy_list == dummy_list[0]):
#             print("Perfect")
#         else:
#             print("Please make sure locked variables are of equal lengths.")
#             return
#
#         if not A_LOCKED:
#             A = np.ones(lock_len)*(-1)
#         if not B_LOCKED:
#             B = np.ones(lock_len)*(-1)
#         if not a_LOCKED:
#             a = np.ones(lock_len)*(-1)
#         if not b_LOCKED:
#             b = np.ones(lock_len)*(-1)
#         if not Ri_LOCKED:
#             Ri = np.ones(lock_len)*(-1)
#         if not L_LOCKED:
#             L = np.ones(lock_len)*(-1)
#         if not Req_LOCKED:
#             L = np.ones(lock_len)*(-1)
#
#         lock = list(zip(A, B, a, b, Ri, L))
#         print(lock)
#
#         slice = []
#         for s in space:
#             ll = []
#             for z in lock:
#                 ll = [i for (i, j) in zip(s, z) if i == j]
#
#                 if len(ll) == count:
#                     slice.append(s)
#
#         df = pd.DataFrame(slice, columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha"])
#         print(df)
#
#     df = pd.DataFrame(space, columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha"])
#     print(df)
#
#     return df
#
#
# def check_input(s):
#     # s = "range(16, 23, 10)"
#     # s = "randrange(16, 23, 10)"
#     # s = "[16, 23, 10]"
#     # s = 1, 2, 3
#     # s = 2
#
#     if "range" in s and "rand" not in s:
#         s = s.replace('range', '')
#         try:
#             l = ast.literal_eval(s)
#             return np.linspace(l[0], l[1], l[2])
#         except:
#             print("Please check inputs.")
#     elif "randrange" in s:
#         s = s.replace('randrange', '')
#         try:
#             l = ast.literal_eval(s)
#             ll = np.random.uniform(l[0], l[1], l[2])
#             return ll
#         except:
#             print("Please check inputs.")
#     else:
#         try:
#             ll = ast.literal_eval(s)
#             if isinstance(ll, int):
#                 ll = [ll]
#             return ll
#         except:
#             print("Please check inputs.")
#
#     return 1
#
#
# iBP = 'left'
# BP = 'none'
# type = "IO"
# freq = 400.79
# pseudo_shape_space = {}
#
# A = [1]
# B = [4, 6]
# a = [7]
# b = [9]
# Ri = [2]
# L = [0]
# Req = [10]
#
#
# A = check_input('range(3, 17, 4)')
# A = np.around(A, decimals=2)
# a = check_input('1')
# a = np.around(a, decimals=2)
# ic(a)
#
# var_list = [A, B, a, b, Ri, L, Req]
# lock_list = [False, False, True, False, False, False, False]
# ihc = create_pseudo_shape_space(var_list, lock_list)
# # ihc = pd.DataFrame([[1, 2, 3, 4, 5, 6, 7, 7]], columns=["A", "B", "a", "b", "Ri", "L", "Req", "alpha"])
# ohc = create_pseudo_shape_space(var_list, lock_list)
#
#
# if ihc is not None and ohc is not None:
#     if type == "IO":
#         print("IO")
#         key = 0
#         for indx1, inner_cell in ihc.iterrows():
#             inner_cell = inner_cell.tolist()
#             for indx2, other_cell in ohc.iterrows():
#                 # enforce Req_inner_cell == Req_outer_cell
#                 other_cell = other_cell.tolist()
#                 other_cell[-2] = inner_cell[-2]
#
#                 pseudo_shape_space[key] = {'IC': inner_cell, 'OC': other_cell, 'BP': BP, 'FREQ': freq}
#                 key += 1
#
#     else:
#         print("II")
#         key = 0
#         for indx, inner_cell in ihc.iterrows():
#             pseudo_shape_space[key] = {'IC': inner_cell.tolist(), 'OC': inner_cell.tolist(), 'BP': BP, 'FREQ': freq}
#             key += 1
#
#     with open("test.json", 'w') as file:
#         file.write(json.dumps(pseudo_shape_space, indent=4, separators=(',', ': ')))
# else:
#     print("Please check inputs.")
#
# print(eval('np.linspace(73.52-2.5, 73.52+2.5, 6)'))

# alpha = np.pi/2
#
# A = np.array([1, 0, 0])
# B = np.array([-1, 0, 0])
# a = np.array([0, 1, 0])
# b = np.array([0, -1, 0])
# Ri = np.array([0, 0, 1])
# L = np.array([0, 0, -1])
# Req = np.array([np.cos(alpha), np.cos(alpha), np.cos(alpha)])
#
# ll = [A, B, a, b, Ri, L, Req]
# for l in ll:
#     ax.plot([0, 1*l[0]], [0, 1*l[1]], [0, 1*l[2]])
# # plt.show()
# # df1 = fr.json_reader(r"D:\Dropbox\CEMCodesHub\SampleProject_s\Cavities\C3794_1Cell.json")
# # df2 = fr.json_reader(r"D:\Dropbox\CEMCodesHub\C538\Cavities\C538_1Cell.json")
# df2 = pd.ExcelFile(r"D:\Dropbox\CEMCodesHub\C538\PostprocessingData\Data\C538_1Cell_ABCI_SLANS.xlsx")
# df2 = df2.parse("Sheet1")
# df_list = [df2]
# for df in df_list:
#     # print(df)
#     df = df/df.max()
#     # print(data)
#
#     x = []
#     y = []
#     z = []
#     c = []  # color
#     fmx = [1092, 865, 901, 1128, 1227, 867, 1155, 1200, 1083, 1104]
#     for ind, val in df.iterrows():
#         # print(val, type(val["A"]))
#         # p = val["A"]*A + val["B"]*B + val["ai"]*a + val["bi"]*b + val["Ri"]*Ri + val["L"]*L + val["Req"]*Req
#         p = val["A"]*A + val["B"]*B + val["a2"]*a + val["b3"]*b + val["Ri"]*Ri + val["L"]*L + val["Req"]*Req
#         x.append(p[0])
#         y.append(p[1])
#         z.append(p[2])
#         c.append(val["fmax[max(0.7<f<0.77)]"])
#
#         if ind in fmx:
#             ax.scatter(p[0], p[1], p[2], c='k')
#         if ind == 1093 or ind == 1094 or ind == 1095 or ind == 1096:
#             ax.scatter(p[0], p[1], p[2], c='r')
#
#     p = ax.scatter(x, y, z, c=c, alpha=.1, cmap='jet')
#
# plt.colorbar(p)
# plt.show()
# Zt = 109.98
# # Yt = 1/Zt
# # f = 519.1e6
# # S21 = 0.22
# # c = 3e8
# #
# # Y11 = (2*Yt*(np.exp(-1j*2*2*np.pi*f/c) - S21))/S21
# # Y11_mag = abs(Y11)
# # print(Y11_mag)

# import numpy as np
# from itertools import product
#
#
# def squiggle_xy(a, b, c, d, i=np.arange(0.0, 2*np.pi, 0.05)):
#     return np.sin(i*a)*np.cos(i*b), np.sin(i*c)*np.cos(i*d)
#
#
# fig11 = plt.figure(figsize=(8, 8), constrained_layout=False)
#
# # gridspec inside gridspec
# outer_grid = fig11.add_gridspec(4, 4, wspace=0.0, hspace=0.0)
#
# for i in range(16):
#     inner_grid = outer_grid[i].subgridspec(3, 3, wspace=0.0, hspace=0.0)
#     a, b = int(i/4)+1, i % 4+1
#     for j, (c, d) in enumerate(product(range(1, 4), repeat=2)):
#         ax = fig11.add_subplot(inner_grid[j])
#         ax.plot(*squiggle_xy(a, b, c, d))
#         ax.set_xticks([])
#         ax.set_yticks([])
#         fig11.add_subplot(ax)

# all_axes = fig11.get_axes()

# # show only the outside spines
# for ax in all_axes:
#     for sp in ax.spines.values():
#         sp.set_visible(False)
#     if ax.is_first_row():
#         ax.spines['top'].set_visible(True)
#     if ax.is_last_row():
#         ax.spines['bottom'].set_visible(True)
#     if ax.is_first_col():
#         ax.spines['left'].set_visible(True)
#     if ax.is_last_col():
#         ax.spines['right'].set_visible(True)

# plt.show()

# folder = fr"D:\Dropbox\CEMCodesHub\C800MHz\SimulationData\SLANS"
# folders = os.listdir(folder)
# for d in folders:
#     # delete SLANS_EXE folder
#     if os.path.exists(fr'{folder}/{d}/SLANS_exe'):
#         shutil.rmtree(fr'{folder}/{d}/SLANS_exe')
#         print(fr"Removed from {d}")
# from PyQt5 import  QtWidgets
# import os
# import numpy as np
# from numpy import cos
# from mayavi.mlab import contour3d
#
# os.environ['ETS_TOOLKIT'] = 'qt5'
# from pyface.qt import QtGui, QtCore
# from traits.api import HasTraits, Instance, on_trait_change
# from traitsui.api import View, Item
# from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
#
# ## create Mayavi Widget and show
#
#
# class Visualization(HasTraits):
#     scene = Instance(MlabSceneModel, ())
#
#     @on_trait_change('scene.activated')
#     def update_plot(self):
#     ## PLot to Show
#         x, y, z = np.ogrid[-3:3:60j, -3:3:60j, -3:3:60j]
#         t = 0
#         Pf = 0.45+((x*cos(t))*(x*cos(t)) + (y*cos(t))*(y*cos(t))-(z*cos(t))*(z*cos(t)))
#         obj = contour3d(Pf, contours=[0], transparent=False)
#
#     view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
#                      height=250, width=300, show_label=False),
#                 resizable=True )
#
#
# class MayaviQWidget(QtGui.QWidget):
#     def __init__(self, parent=None):
#         QtGui.QWidget.__init__(self, parent)
#         layout = QtGui.QVBoxLayout(self)
#         layout.setContentsMargins(0,0,0,0)
#         layout.setSpacing(0)
#         self.visualization = Visualization()
#
#         self.ui = self.visualization.edit_traits(parent=self,
#                                                  kind='subpanel').control
#         layout.addWidget(self.ui)
#         self.ui.setParent(self)
#
#
# #### PyQt5 GUI ####
# class Ui_MainWindow(object):
#     def setupUi(self, MainWindow):
#
#     ## MAIN WINDOW
#         MainWindow.setObjectName("MainWindow")
#         MainWindow.setGeometry(200,200,1100,700)
#
#     ## CENTRAL WIDGET
#         self.centralwidget = QtWidgets.QWidget(MainWindow)
#         self.centralwidget.setObjectName("centralwidget")
#         MainWindow.setCentralWidget(self.centralwidget)
#
#     ## GRID LAYOUT
#         self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
#         self.gridLayout.setObjectName("gridLayout")
#
#
#     ## BUTTONS
#         self.button_default = QtWidgets.QPushButton(self.centralwidget)
#         self.button_default.setObjectName("button_default")
#         self.gridLayout.addWidget(self.button_default, 0, 0, 1,1)
#
#         self.button_previous_data = QtWidgets.QPushButton(self.centralwidget)
#         self.button_previous_data.setObjectName("button_previous_data")
#         self.gridLayout.addWidget(self.button_previous_data, 1, 1, 1,1)
#
#     ## Mayavi Widget 1
#         container = QtGui.QWidget()
#         mayavi_widget = MayaviQWidget(container)
#         self.gridLayout.addWidget(mayavi_widget, 1, 0,1,1)
#     ## Mayavi Widget 2
#         container1 = QtGui.QWidget()
#         mayavi_widget = MayaviQWidget(container1)
#         self.gridLayout.addWidget(mayavi_widget, 0, 1,1,1)
#
#     ## SET TEXT
#         self.retranslateUi(MainWindow)
#         QtCore.QMetaObject.connectSlotsByName(MainWindow)
#
#     def retranslateUi(self, MainWindow):
#         _translate = QtCore.QCoreApplication.translate
#         MainWindow.setWindowTitle(_translate("MainWindow", "Simulator"))
#         self.button_default.setText(_translate("MainWindow","Default Values"))
#         self.button_previous_data.setText(_translate("MainWindow","Previous Values"))
#
#
# if __name__ == "__main__":
#     import sys
#     app = QtWidgets.QApplication(sys.argv)
#     app.setStyle('Fusion')
#     MainWindow = QtWidgets.QMainWindow()
#
#     ui = Ui_MainWindow()
#     ui.setupUi(MainWindow)
#     MainWindow.show()
#     sys.exit(app.exec_())

# # # load data
# # filename = fr"D:\Dropbox\CEMCodesHub\Cavity800\SimulationData\ea_results\Run3\Generation9.xlsx"
# # filename1 = fr"D:\Dropbox\CEMCodesHub\Cavity800\SimulationData\SLANS\Generation13.xlsx"
# # filename2 = fr'D:\Dropbox\CEMCodesHub\C800MHz\PostprocessingData\Data\ttt.xlsx'
# # df = pd.read_excel(filename, 'Sheet1')
# # df1 = pd.read_excel(filename1, 'Sheet1')
# # df2 = pd.read_excel(filename2, 'Sheet1')
# # v = ['A', 'B', 'a', 'b', 'Ri', 'Req', 'Epk/Eacc', 'Bpk/Eacc', 'R/Q']
# # # v = ['A', 'B', 'a', 'b', 'Ri', 'Req', 'Epk/Eacc_1', 'Bpk/Eacc_1', 'Rsh/Q_1']
# # # v = ['Epk/Eacc', 'Bpk/Eacc', 'R/Q']
# #
# # # for p in v:
# # #     # df[p] = df[p] / df[p].abs().max()
# # #     # print(df['B'].describe())
# # #     df[p].plot.density()
# # sd = np.linspace(30, 80, 1000)
# # dff = pd.DataFrame(sd, columns=['A'])
# # # print(dff)
# # var = 'A'
# # df[var].plot.density()
# # df1[var].plot.density()
# # df2[var].plot.density()
# # dff[var].plot.density()
# # #
# # plt.legend()
# # plt.ylim(0, 0.15)
# # plt.show()

# df = pd.DataFrame(columns=['key', 'A', 'B', 'C'])
# df[['key', 'A', 'B', 'C']] = [['a', 1, 2, 3], ['b', 4, 5, 6], ['c', 7, 8, 9], ['d', 10, 11, 12]]
# df = df.set_index('key')
# print(df)
#
# for index, row in df.iterrows():
#     print(row.tolist())

def stroud(p):
    # Stroud-3 method
    #
    # Input parameters:
    #  p   number of dimensions
    # Output parameters:
    #  nodes   nodes of quadrature rule in [0,1]^p (column-wise)
    #

    nodes = np.zeros((p, 2 * p))
    coeff = np.pi / p
    fac = np.sqrt(2 / 3)

    for i in range(2 * p):
        for r in range(int(np.floor(0.5 * p))):
            k = 2 * r
            nodes[k, i] = fac * np.cos((k+1) * (i+1) * coeff)
            nodes[k + 1, i] = fac * np.sin((k+1) * (i+1) * coeff)

        if 0.5 * p != np.floor(0.5 * p):
            nodes[-1, i] = ((-1) ** (i+1)) / np.sqrt(3)

    # transform nodes from [-1,+1]^p to [0,1]^p
    nodes = 0.5 * nodes + 0.5

    return nodes


def quad_stroud3(rdim, degree):
    # data for Stroud-3 quadrature in [0,1]^k
    # nodes and weights
    nodes = stroud(rdim)
    nodestr = 2. * nodes - 1.
    weights = (1 / (2 * rdim)) * np.ones((2 * rdim, 1))
    print("weights: ", weights)
    print("nodes: ", nodes)

    # evaluation of Legendre polynomials
    bpoly = np.zeros((degree + 1, rdim, 2 * rdim))
    for l in range(rdim):
        for j in range(2 * rdim):
            bpoly[0, l, j] = 1
            bpoly[1, l, j] = nodestr[l, j]
            for i in range(1, degree):
                bpoly[i + 1, l, j] = ((2 * (i+1) - 1) * nodestr[l, j] * bpoly[i, l, j] - i * bpoly[i - 1, l, j]) / (i+1)

    # standardisation of Legendre polynomials
    for i in range(1, degree + 1):
        bpoly[i, :, :] = bpoly[i, :, :] * np.sqrt(2 * (i+1) - 1)

    return nodes, weights, bpoly


def weighted_mean_obj(tab_var, weights):
    rows_sims_no, cols = np.shape(tab_var)
    no_weights, dummy = np.shape(weights)  # z funckji quadr_stroud wekt columnowy

    if rows_sims_no == no_weights:
        expe = np.zeros((cols, 1))
        outvar = np.zeros((cols, 1))
        for i in range(cols):
            expe[i, 0] = np.dot(tab_var[:, i], weights)
            outvar[i, 0] = np.dot(tab_var[:, i]**2, weights)

        stdDev = np.sqrt(outvar - expe**2)
    else:
        expe = 0
        stdDev = 0
        ic('Cols_sims_no!=No_weights')

    return expe, stdDev

#
# def uq():
#     p_true = np.array([1, 2, 3, 4, 5])
#     ic(p_true)
#     rdim = len(p_true)  # How many variabels will be considered as random in our case 5
#     degree = 1
#
#     #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
#     flag_stroud = 1
#     if flag_stroud == 1:
#         nodes, weights, bpoly = quad_stroud3(rdim, degree)
#         nodes = 2. * nodes - 1.
#     elif flag_stroud == 2:
#         nodes, weights, bpoly = quad_stroud3(rdim, degree)  # change to stroud 5 later
#         nodes = 2. * nodes - 1.
#     else:
#         ic('flag_stroud==1 or flag_stroud==2')
#
#     #  mean value of geometrical parameters
#     p_init = np.zeros(np.shape(p_true))
#
#     no_parm, no_sims = np.shape(nodes)
#     delta = 0.05  # or 0.1
#
#     Ttab_val_f = []
#
#     for i in range(no_sims):
#         p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
#         p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
#         p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
#         p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
#         p_init[4] = p_true[4] * (1 + delta * nodes[4, i])
#
#         tab_val_f = [np.random.randint(5), np.random.randint(5)]
#
#         Ttab_val_f.append(tab_val_f)
#
#     ic(np.array(Ttab_val_f).T)
#     v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
#     ic(list(v_expe_fobj.T), list(v_stdDev_fobj.T[0]))
#
#
# p = 5
# degree = 1
# # nodes, weights, bpoly = quad_stroud3(p, degree)
# # ic(nodes, weights, bpoly[:, :, 1])
# uq()
#
# uq_result_dict = {'G0_C0_P': [801.62206,
#                                  3.5669454055770395,
#                                  2.2964732427293555,
#                                  0.07218374831898976,
#                                  5.081215159534722,
#                                  0.05214187332623576,
#                                  96.86471800000001,
#                                  2.179102942247872],
#                      'G0_C1_P': [801.62692,
#                                  4.159920815526951,
#                                  2.036639874786366,
#                                  0.043256932856938476,
#                                  5.069242453781859,
#                                  0.06063893309179642,
#                                  95.755528,
#                                  2.338346393864051],
#                      'G0_C2_P': [801.6367200000001,
#                                  4.818308504389939,
#                                  1.9911435243194302,
#                                  0.03175986811962072,
#                                  5.064680566274672,
#                                  0.06920703401261843,
#                                  94.275848,
#                                  2.4871270745857412],
#                      'G0_C3_P': [801.63716,
#                                  5.59864282593621,
#                                  2.0182945504535574,
#                                  0.02948287780113185,
#                                  5.067588181231035,
#                                  0.07709301994200095,
#                                  92.46944,
#                                  2.6174705920945804],
#                      'G0_C4_P': [801.61972,
#                                  6.618941107281137,
#                                  2.0818413581288855,
#                                  0.035112410132603834,
#                                  5.0838723764953055,
#                                  0.08302605338203956,
#                                  90.347102,
#                                  2.7185604221121293],
#                      'G0_C5_P': [802.76827,
#                                  5.029137168160952,
#                                  2.1266994327525803,
#                                  0.025972610635272203,
#                                  5.160077537544655,
#                                  0.07710495281859536,
#                                  87.333348,
#                                  2.6843850826020694]}
#
# df = pd.DataFrame.from_dict(uq_result_dict, orient='index')
# df.columns = list('abcdefgh')
# ic(df)
#
# df.index.name = 'key'
# df.reset_index(inplace=True)
# ic(df)
# df[['r', 'q']] = np.array([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]]).T
# df[['r', 'q']] = np.array([[4, 5, 3, 4, 5, 6], [6, 7, 3, 4, 5, 6]]).T
# # ic(df.to_numpy().T)
# ic(tuple(df.iloc[0].values))

# def color(row):
#     # if row.isnull().values.any():
#     if row[0] in df['key'].tolist()[0:2]:
#         return ['background-color: #6bbcd1'] * len(row)
#     return [''] * len(row)
#
#
# # df = pd.DataFrame([[7, 5, 6], [1, 2, 3], [4, 5, 6]], index=[list('abc')])
# # Save Styler Object for Later
# styler = df.style
# # Apply Styles (This can be chained or on separate lines)
# # styler.applymap(lambda x: 'background-color : yellow' if x > 1 else '')
# styler.apply(color, axis=1)
# # Export the styler to excel
# styler.to_excel('Output.xlsx')

# conv = []
# for n in np.arange(1):
#     f = []
#     for x in np.linspace(-1, 1, 100):
#         for y in np.linspace(-1, 1, 100):
#             f.append((1-x)**2 + 100*(y - x**2)**2)
#
#     ic(sum(f)/len(f))
#     conv.append(sum(f)/len(f))
#
# plt.plot()
#
# print(2+-3)

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import colors as mcolors
# from matplotlib import cm
#
# ## Figures config
# colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
# LEGEND_SIZE = 15
# TITLE_SIZE = 25
# AXIS_SIZE = 15
#
# from emukit.core.initial_designs import RandomDesign
# from GPy.models import GPRegression
# from emukit.model_wrappers import GPyModelWrapper
# from emukit.sensitivity.monte_carlo import MonteCarloSensitivity
# from emukit.core import ContinuousParameter, ParameterSpace
# from emukit.test_functions.sensitivity import Ishigami ### change this one for the one in the library
# from GPy.models import GPRegression
# from emukit.model_wrappers import GPyModelWrapper
# from emukit.sensitivity.monte_carlo import MonteCarloSensitivity
#
#
# filename = fr'D:\Dropbox\CEMCodesHub\C800MHz\PostprocessingData\Data\GridSimulation_Data.xlsx'
# df = pd.read_excel(filename, 'Sheet1')
#
# # ishigami = Ishigami(a=5, b=0.1)
# # #
# # variable_domain = (-np.pi, np.pi)
# # space = ParameterSpace([ContinuousParameter('x1', variable_domain[0], variable_domain[1]),
# #                         ContinuousParameter('x2', variable_domain[0], variable_domain[1]),
# #                         ContinuousParameter('x3', variable_domain[0], variable_domain[1])])
#
# variable_domain = (-np.pi, np.pi)
# space = ParameterSpace([ContinuousParameter('A', 30, 80),
#                         ContinuousParameter('B', 30, 80),
#                         ContinuousParameter('a2', 10, 60),
#                         ContinuousParameter('b3', 10, 60),
#                         ContinuousParameter('Ri', 60, 85)])
#
# # desing = RandomDesign(space)
# # X = desing.get_samples(500)
# # Y = ishigami.fidelity1(X)[:, None]
#
# ic(df[['A', 'B', 'a2', 'b3', 'Ri']].to_numpy())
# # ic(X)
# ic(np.atleast_2d(df['Epk/Eacc'].to_numpy()).T)
# # ic(Y)
#
# X = df[['A', 'B', 'a2', 'b3', 'Ri']].to_numpy()
# Y = np.atleast_2d(df['Epk/Eacc'].to_numpy()).T
#
# ic("Check 1")
# model_gpy = GPRegression(X, Y)
# ic("Check 2")
# model_emukit = GPyModelWrapper(model_gpy)
# ic("Check 3")
# model_emukit.optimize()
# ic("Check 4")
#
# senstivity_ishigami_gpbased = MonteCarloSensitivity(model = model_emukit, input_domain = space)
# ic("Check 5")
# main_effects_gp, total_effects_gp, _ = senstivity_ishigami_gpbased.compute_effects(num_monte_carlo_points = 10000)
# ic("Check 6")
# main_effects_gp = {ivar: main_effects_gp[ivar][0] for ivar in main_effects_gp}
# ic("Check 7")
#
# d = {#'Sobol True': ishigami.main_effects,
#      'GP Monte Carlo': main_effects_gp}
#
# pd.DataFrame(d).plot(kind='bar', figsize=(12, 5))
# plt.title('First-order Sobol indexes - Ishigami', fontsize=TITLE_SIZE)
# plt.ylabel('% of explained output variance', fontsize=AXIS_SIZE)
# plt.show()


def create_cst_sweep_import_random(p, v, v_int):
    import math
    import random
    # p - polynomial order
    # v - random variables

    # calculate number of random points needed
    n = len(v)
    n_rand = 2 * int((n - 1)*math.factorial(n + p)/(math.factorial(n)*math.factorial(p)))
    ic(p, n, n_rand)

    d = {}
    for i in range(n_rand):
        d[i] = [random.uniform(vi[0], vi[1]) for vi in v_int]

    df = pd.DataFrame.from_dict(d, orient='index', columns=v)
    # ic(df)
    writePath = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\hc_random_vector.txt'

    with open(writePath, 'w') as f:
        dfAsString = df.to_string(header=True, index=False)
        f.write(dfAsString)


var = ['lh1', 'lh3', 'lh4', 'dh3', 'alpha_h', 'ch2', 'r_cyl', 'offset_y']

var_int = [[150, 190], [35, 45], [50, 70], [12, 20], [90, 120], [2, 4], [10, 18], [10, 30]] # big interval
# var_int = [[155, 175], [37.5, 42.5], [55, 65], [15, 18], [100, 115], [2, 4], [12, 16], [15, 25]]  # small interval

# create_cst_sweep_import_random(1, var, var_int)


def create_cst_sweep_import_sobol_sequence(dim, index, columns, l_bounds, u_bounds):
    from scipy.stats import qmc
    sampler = qmc.Sobol(d=dim, scramble=False)
    sample = sampler.random_base2(m=index)
    ic(qmc.discrepancy(sample))
    sample = qmc.scale(sample, l_bounds, u_bounds)
    print(sample)

    df = pd.DataFrame(sample, columns=columns)
    # ic(df)
    writePath = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\dqw_random_vector.txt'

    with open(writePath, 'w') as f:
        dfAsString = df.to_string(header=True, index=False)
        f.write(dfAsString)


columns = ['shaft_y', 'bend_out_sec_prob_length', 'bend_out_cap_gap', 'Window_margin', 'Shift_from_center',
           'Cap_y', 'Cap_thickness', 'Cap_Height', 'Cap_Gap', 'Bend_out_chamber_Length', 'BP_HOM_Penetration']
nominal_values = [215, 60, 5.1, 12, -15, -15.5, 3, 49.6, 15, 30.5, 20]
variation = 0.2
l_bounds, u_bounds = [], []
for nv in nominal_values:
    if nv < 0:
        l_bounds.append(nv*(1 + variation))
        u_bounds.append(nv*(1 - variation))
    else:
        l_bounds.append(nv*(1 - variation))
        u_bounds.append(nv*(1 + variation))
print([list(x) for x in zip(l_bounds, u_bounds)])
# create_cst_sweep_import_sobol_sequence(11, 9, columns, l_bounds, u_bounds)


def plot_settings():
    import matplotlib as mpl
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16

    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['legend.title_fontsize'] = 14

    mpl.rcParams['figure.figsize'] = [5.5, 7]
    mpl.rcParams['figure.dpi'] = 100


plot_settings()


def plot_bar_single():
    import matplotlib.ticker as mticker
    from itertools import combinations, product
    def pairs(*lists):
        for t in combinations(lists, 2):
            for pair in product(*t):
                yield pair


    filename = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\Angle_sweep_2hc_1fpc.xlsx'
    data = pd.read_excel(filename, sheet_name='Sheet1')

    f_multi = data['f'].to_numpy().reshape((56, 10))
    Qext_multi = data['Qext'].to_numpy().reshape((56, 10))
    RQT_multi = data['RQT'].to_numpy().reshape((56, 10))
    ZT_multi = data['ZT'].to_numpy().reshape((56, 10))

    alpha, alpha_fpc = [0, 45, 90, 135, 180, 225, 270, 330], [60, 90, 135, 180, 225, 270, 300]
    angle_pair = []
    for pair in pairs(alpha_fpc, alpha):
        angle_pair.append(pair)
    print(len(angle_pair))
    ic(angle_pair)
    # ic(f_multi)
    fig, ax = plt.subplots(3,sharex=True)
    Qext_max, RQT_max, ZT_max = [], [], []
    for i, f in enumerate(f_multi):
        # ax[0].scatter(f, Qext_multi[i], label=angle_pair[i], marker='o', ec='k')
        # ax[0].scatter(f, RQT_multi[i], label=angle_pair[i], marker='o', ec='k')
        # ax[0].scatter(f, ZT_multi[i], label=angle_pair[i], marker='o', ec='k')

        # save maximum
        Qext_max.append(np.max(Qext_multi[i]))
        RQT_max.append(np.max(RQT_multi[i]))
        ZT_max.append(np.max(ZT_multi[i]))

    N = len(angle_pair)
    ind = np.arange(N)
    width = 0.25


    bars = ax[0].bar(ind, Qext_max, label='$Q_{ext}$', color='lightcoral', edgecolor='k')
    ax[0].plot(ind, Qext_max, lw=2, marker='o', mec='k')
    ax[0].set_ylabel('max($Q_\mathrm{ext})$ $[\cdot]$')
    bars1 = ax[1].bar(ind, RQT_max, label='$(R/Q)_T$', color='tab:blue',  edgecolor='k')
    ax[1].plot(ind, RQT_max, lw=2, marker='o', color='lightcoral', mec='k')
    ax[1].set_ylabel('max($(R/Q)_\mathrm{T}$ $ ~\mathrm{[\Omega]}$')
    bars2 = ax[2].bar(ind, ZT_max, label='$Z_T$', color='black', edgecolor='k')
    ax[2].plot(ind, ZT_max, lw=2, marker='o', mec='white')
    ax[2].set_ylabel('max($Z_\mathrm{T}$) $[\mathrm{[k \Omega/m]}$')

    # ax[0].bar_label(bars, fmt='%.2e')##
    # ax[1].bar_label(bars, fmt='%.2e')##
    # ax[2].bar_label(bars, fmt='%.2e')##

    ticks_loc = ax[0].get_xticks()

    ax[2].xaxis.set_major_locator(mticker.FixedLocator(range(len(angle_pair))))
    ax[2].set_xticklabels(angle_pair, rotation=90)
    ax[2].set_xlabel(r'($\alpha_\mathrm{HC, FPC}, \alpha_\mathrm{HC}$)')


    ax[0].set_yscale('log')
    # ax[1].set_yscale('log')
    ax[2].set_yscale('log')
    # ax.set_ylim(0, 3e6)
    # ax.legend(loc="upper center", ncol=10, bbox_to_anchor=(0.5, 1.65))
    # ax.legend(loc="upper right")


    plt.show()


def plot_bar_double():
    import matplotlib.ticker as mticker
    # create_cst_sweep_import(2, var, var_int)
    from itertools import combinations, product

    def pairs(*lists):
        for t in combinations(lists, 2):
            for pair in product(*t):
                yield pair

    aa = [135, 45]
    filename = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\Angle_sweep_2hc_1fpc.xlsx'
    data_reference = pd.read_excel(filename, sheet_name='ref')
    data_opt = pd.read_excel(filename, sheet_name=f'{aa[0]}_{aa[1]}_opt')
    data1 = pd.read_excel(filename, sheet_name=f'{aa[0]}_{aa[1]}')

    data_list = [data_reference, data1, data_opt]
    fig, ax = plt.subplots(3, sharex=True)

    label = ['(135, 225) (ref)', f'({aa[0]}, {aa[1]})', f'({aa[0]}, {aa[1]})_opt']

    for i, data in enumerate(data_list):
        if i == 0:
            f = data['f'].to_numpy()
            Qext = data['Qext'].to_numpy()
            RQT = data['RQT'].to_numpy()
            ZT = data['ZT'].to_numpy()

            width = 0.25

            ax[0].bar(f, Qext, alpha=0.5, label=label[i])
            ax[0].set_ylabel('$Q_\mathrm{ext}$ $[\cdot]$')

            ax[1].bar(f, RQT, alpha=0.5)
            ax[1].set_ylabel('$(R/Q)_\mathrm{T}$ $ ~\mathrm{[\Omega]}$')

            ax[2].bar(f, ZT, alpha=0.5)
            ax[2].set_ylabel('$Z_\mathrm{T}$ $\mathrm{[k \Omega/m]}$')

            ax[2].set_xlabel(r'$\alpha_\mathrm{HC, FPC} = ' + fr'{aa[0]}' + r'^\circ, \alpha_\mathrm{HC} = ' + fr'{aa[0]}' + '^\circ$')

            ax[0].set_yscale('log')
            ax[2].set_yscale('log')
            # ax.set_ylim(0, 3e6)
            # ax.legend(loc="upper center", ncol=10, bbox_to_anchor=(0.5, 1.65))
            # ax.le
        else:
            f = data['f'].to_numpy()
            Qext = data['Qext'].to_numpy()
            RQT = data['RQT'].to_numpy()
            ZT = data['ZT'].to_numpy()

            width = 0.25

            ax[0].stem(f, Qext, linefmt='k')
            ax[0].scatter(f, Qext, marker='o', linewidth=1, edgecolor='k', zorder=10, label=label[i])
            ax[0].set_ylabel('$Q_\mathrm{ext}$ $[\cdot]$')

            ax[1].stem(f, RQT, linefmt='k')
            ax[1].scatter(f, RQT, marker='o', linewidth=1, edgecolor='k', zorder=10)
            ax[1].set_ylabel('$(R/Q)_\mathrm{T}$ $ ~\mathrm{[\Omega]}$')

            ax[2].stem(f, ZT, linefmt='k')
            ax[2].scatter(f, ZT,  marker='o', linewidth=1, edgecolor='k', zorder=10)
            ax[2].set_ylabel('$Z_\mathrm{T}$ $\mathrm{[k \Omega/m]}$')

            ax[2].set_xlabel(r'$\alpha_\mathrm{HC, FPC} = ' + fr'{aa[0]}' + r'^\circ, \alpha_\mathrm{HC} = ' + fr'{aa[0]}' + '^\circ$')

            ax[0].set_yscale('log')
            ax[2].set_yscale('log')
            # ax.set_ylim(0, 3e6)
            # ax.legend(loc="upper center", ncol=10, bbox_to_anchor=(0.5, 1.65))
            # ax.legend(loc="upper right")

    lines, labels = ax[0].get_legend_handles_labels()

    fig.legend(lines, labels, loc="upper center", title=r'($\alpha_\mathrm{HC, FPC}, \alpha_\mathrm{HC}$)',
               ncol=len(labels), fancybox=True, shadow=False)

    plt.subplots_adjust(left=0.176, right=0.973, top=0.9, bottom=0.112)
    # plt.tight_layout()
    plt.show()


# # plot_bar_double()
#
# import matplotlib
# matplotlib.use('qt5agg')
# import matplotlib.pyplot as plt
#
# L = 300e-3
# b = 10e-3
# Z0 = 377
# c = 299792458
# sigma = 15*1e-3
# sigma_c = 3e3  # conductivity
# df = pd.read_excel(fr"D:\CST Studio\WAKE\compare.xlsx", 'Sheet1')
# f = df['f']*1e9
# ic(f)
# a = (1-1j)*(L/(2*np.pi*b))*np.sqrt(Z0*2*np.pi*f/(2*c*sigma_c))
# a = np.abs(a)
# a.to_csv(fr"D:\CST Studio\WAKE\analytic.txt", header=None, index=None)
#
# plt.plot(f, np.abs(a))
# plt.tight_layout()
# plt.show()
# ############################################
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
# from matplotlib.patches import Rectangle
# from matplotlib.text import Text
# from matplotlib.image import AxesImage
# import numpy as np
# from numpy.random import rand
#
#
# # Fixing random state for reproducibility
# np.random.seed(19680801)
#
# fig, (ax1, ax2) = plt.subplots(2, 1)
# ax1.set_title('click on points, rectangles or text', picker=True)
# ax1.set_ylabel('ylabel', picker=True, bbox=dict(facecolor='red'))
# line, = ax1.plot(rand(100), 'o', picker=True, pickradius=5)
#
# # Pick the rectangle.
# ax2.bar(range(10), rand(10), picker=True)
# for label in ax2.get_xticklabels():  # Make the xtick labels pickable.
#     label.set_picker(True)
#
#
# def onpick1(event):
#     if isinstance(event.artist, Line2D):
#         thisline = event.artist
#         xdata = thisline.get_xdata()
#         ydata = thisline.get_ydata()
#         ind = event.ind
#         print('onpick1 line:', np.column_stack([xdata[ind], ydata[ind]]))
#     elif isinstance(event.artist, Rectangle):
#         patch = event.artist
#         print('onpick1 patch:', patch.get_path())
#     elif isinstance(event.artist, Text):
#         text = event.artist
#         print('onpick1 text:', text.get_text())
#
#
# fig.canvas.mpl_connect('pick_event', onpick1)
# plt.show()
# #################################################


# import quadpy
#
# rdim, degree = 5, 1
# scheme = quadpy.cn.stroud_cn_5_2(5)
# quad_stroud3(rdim, degree)
# print()
# print(scheme.weights)
# print(scheme.points)

#########################################
def func(x):
    return 1/x


# def monte_carlo():
#     # calculate volume
#     V = 1
#     N = 100
#
#     Integral = 0
#     for n in range(1, N, 1):
#         Integral = V*np.sum([func(x)/n for x in np.random.uniform(1, 2, n)])
#
#         plt.scatter(n, abs(0.693417 - Integral)/0.693417)
#
#     plt.yscale('log')
#     plt.show()
#
#
# monte_carlo()


def area_circ(r):
    return np.pi*r**2

area_circ(3)