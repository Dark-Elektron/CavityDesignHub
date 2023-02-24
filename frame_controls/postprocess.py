import ast
import json
import math

# import mayavi
import oapackage
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from icecream import ic
from scipy.interpolate import LinearNDInterpolator, griddata
from termcolor import colored
import scipy.signal as sps
import numpy as np
from frame_controls.postprocess_widgets.pandas_widget import PandasModel
import pandas as pd
from ui_files.postprocess import Ui_Postprocess
from analysis_modules.data_module.abci_data import ABCIData, ABCIDataExtraction
from analysis_modules.data_module.slans_data import *
from utils.file_reader import FileReader
from analysis_modules.plot_module.plotter import Plot
from frame_controls.postprocess_widgets.pandas_widget import PandasTV
from frame_controls.postprocess_widgets.ppplot import PPPlot
from utils.write_geometry import writeCavity
import matplotlib as mpl
from utils.shared_classes import *
from utils.shared_functions import *
# from mayavi import mlab

fr = FileReader()
# abci_data = ABCIData()
abci_data_ex = ABCIDataExtraction()
slans_data_ex = SLANSDataExtraction()

# All_outputs_filtered=(find(All_outputs(:,14)<0.05&All_outputs(:,16)<0.3&All_outputs(:,17)<0.3&All_outputs(:,
# 18)<2.2&All_outputs(:,19)<1&All_outputs(:,20)<3.5&All_outputs(:,29)<2.4&All_outputs(:,30)<6.5 ));%for Ri 150 mm
#
# All_outputs_filtered=(find(All_outputs(:,14)<0.1&All_outputs(:,17)<0.5&All_outputs(:,19)<2&All_outputs(:,
# 20)<1&All_outputs(:,29)<2.4&All_outputs(:,30)<6.7 ));%for Ri 160 mm

file_color = 'red'
DEBUG = True


def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))


class PostprocessControl:
    def __init__(self, parent):
        # pp -> postprocess
        self.w_Postprocess = QWidget()

        self.ui = Ui_Postprocess()
        self.ui.setupUi(self.w_Postprocess)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # initialise pandas object
        self.pandas_object_dict = {}

        # initialise ppplot object
        self.ppplot = PPPlot(self)
        self.pppUI = self.ppplot.ui

        self.ppplot_statistics = PPPlot(self)
        self.pppUI_Stat = self.ppplot_statistics.ui

        self.w_dict = {'Dataset': [self.ui.pb_CEM, self.ui.w_Dataset_From_Simulation],
                       'Combine': [self.ui.pb_CST, self.ui.w_Combine_Dataset]}

        # combine dictionaries list
        self.combine_dict_dict = {}

        # handle dict
        self.OF_slider_dict = {}
        self.datasets = {}
        self.pandas_model_dict = {}
        self.sub_window_dict = {}
        self.x_dict = {}
        self.x_inv_dict = {}
        self.x_global_dict = {}
        self.y_dict = {}
        self.y_inv_dict = {}
        self.y_global_dict = {}
        self.z_dict = {}
        self.z_inv_dict = {}
        self.z_global_dict = {}
        self.filtered_df_dict = {}
        self.filtered_df_inverse_dict = {}
        self.annot_dict = {}

        # plot objects
        self.scatter_plot_object = {}
        self.scatter_plot_object_inv = {}
        self.arrow_patch_list = []

        self.constants()
        self.signals()
        self.initUI()

    def constants(self):
        self.MDI_PLOT_COUNT = 0
        self.MDI_DATAFRAME_COUNT = 0

    def signals(self):
        self.ui.pb_Clear.clicked.connect(self.clear_plots)
        # widget display signals
        self.ui.pb_CEM.clicked.connect(lambda: self.toggle_page('CEM'))
        self.ui.pb_CST.clicked.connect(lambda: self.toggle_page('CST'))

        # load dir
        self.ui.pb_Load_Doc1.clicked.connect(lambda: self.load_dir(self.ui.le_Dir1.text(),
                                                                     self.ui.lw_Doc1,
                                                                     self.ui.lw_Selected_Columns_Doc1,
                                                                     self.column_list_1, 0))
        self.ui.pb_Load_Doc2.clicked.connect(lambda: self.load_dir(self.ui.le_Dir2.text(),
                                                                     self.ui.lw_Doc2,
                                                                     self.ui.lw_Selected_Columns_Doc2,
                                                                     self.column_list_2, 1))

        # add items to combine
        self.column_list_1 = []  # list to hold added columns
        self.column_list_2 = []  # list to hold added columns
        self.ui.pb_Add_Doc1.clicked.connect(
            lambda: self.add_to_list_widget(self.ui.lw_Doc1, self.ui.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ui.pb_Add_Doc2.clicked.connect(
            lambda: self.add_to_list_widget(self.ui.lw_Doc2, self.ui.lw_Selected_Columns_Doc2, self.column_list_2))
        # add all
        self.ui.pb_Add_All_Doc1.clicked.connect(
            lambda: self.add_all(self.ui.lw_Doc1, self.ui.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ui.pb_Add_All_Doc2.clicked.connect(
            lambda: self.add_all(self.ui.lw_Doc2, self.ui.lw_Selected_Columns_Doc2, self.column_list_2))

        # combine from folder parallel run
        self.ui.pb_Combine_Parallel.clicked.connect(
            lambda: self.combine_files_parallel_run(self.ui.sb_Proc_Count.value(),
                                                    self.ui.le_Folder.text(), self.ui.le_Filename_Parallel.text()))

        # remove items
        self.ui.pb_Remove_Doc1.clicked.connect(
            lambda: self.remove_from_list_widget(self.ui.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ui.pb_Remove_Doc2.clicked.connect(
            lambda: self.remove_from_list_widget(self.ui.lw_Selected_Columns_Doc2, self.column_list_2))
        # remove all
        self.ui.pb_Remove_All_Doc1.clicked.connect(
            lambda: self.remove_all(self.ui.lw_Doc1, self.ui.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ui.pb_Remove_All_Doc2.clicked.connect(
            lambda: self.remove_all(self.ui.lw_Doc2, self.ui.lw_Selected_Columns_Doc2, self.column_list_2))

        # combine files
        self.ui.pb_Combine.clicked.connect(lambda: self.combine_files(self.ui.le_Filename.text()))

        # load excel file
        self.ui.pb_Select_File.clicked.connect(lambda: self.load_file())
        self.ui.pb_Select_File.clicked.connect(lambda: self.add_new_dataset())

        # objective function/pareto
        self.ui.pb_Plot_Objective.clicked.connect(lambda: self.plot_objective_function_3D() if self.ui.sc_Toggle_3D.isChecked() else
                                                    self.plot_objective_function())

        self.ui.cb_Pareto.stateChanged.connect(lambda: self.plot_pareto_3D() if self.ui.sc_Toggle_3D.isChecked() else
                                                    self.plot_pareto())

        # mdi signals
        self.ui.pb_MDI_Tiled.clicked.connect(lambda: self.ui.mdiArea.tileSubWindows())
        self.ui.pb_MDI_Cascade.clicked.connect(lambda: self.ui.mdiArea.cascadeSubWindows())

        # add/remove filter
        self.ui.pb_Add_DFFilter.clicked.connect(lambda: self.add_filter())
        self.ui.pb_Remove_Filter.clicked.connect(lambda: self.remove_filter(self.ui.pb_Remove_Filter))

        # enable/disable filter
        self.ui.cb_Filter_Toggle_0.clicked.connect(lambda: self.apply_filter())
        self.ui.hs_Alpha.valueChanged.connect(lambda: self.change_alpha(self.ui.hs_Alpha.value() / 100))

        # enable interactive
        self.ui.cb_Interactive.clicked.connect(lambda: self.toggle_interactive())

        # load shape space
        self.ui.pb_Select_Shape_Space.clicked.connect(
            lambda: self.open_file(self.ui.le_Shape_Space, self.cb_Shape_Space_Keys))

        # select folder
        self.ui.pb_Select_Folder.clicked.connect(lambda: self.open_folder(self.ui.le_SimulationData_Folder))

        # extract data
        self.ui.pb_Extract.clicked.connect(lambda: self.process_folder_data())

        # ABCI/SLANS inputs
        self.ui.cb_Code.currentTextChanged.connect(lambda: self.code_change_control())

        # sequential vs parallel
        self.ui.cb_Run_Mode.currentTextChanged.connect(lambda: self.run_mode_control())

        #
        self.f1_dict = {}
        self.f2_dict = {}
        self.f3_dict = {}
        self.ui.ccb_Objective_Function_F1.currentTextChanged.connect(lambda: self.populate_objective_function(
            self.ui.ccb_Objective_Function_F1.currentText().split(', '), self.ui.ccb_Objective_Function_F1,
            self.ui.tw_Objective_Function_F1, self.f1_dict))

        self.ui.ccb_Objective_Function_F2.currentTextChanged.connect(lambda: self.populate_objective_function(
            self.ui.ccb_Objective_Function_F2.currentText().split(', '), self.ui.ccb_Objective_Function_F2,
            self.ui.tw_Objective_Function_F2, self.f2_dict))

        self.ui.ccb_Objective_Function_F3.currentTextChanged.connect(lambda: self.populate_objective_function(
            self.ui.ccb_Objective_Function_F3.currentText().split(', '), self.ui.ccb_Objective_Function_F3,
            self.ui.tw_Objective_Function_F3, self.f3_dict))

    def initUI(self):

        # splitter
        self.ui.sp_Left_Right_Container.setStretchFactor(1, 3)
        
        # hide code
        self.ui.w_SLANS.hide()

        self.ui.w_Toggle_F3.hide()

        # diable combine button
        self.ui.pb_Combine.setEnabled(False)

        # initial table dataframe
        # self.df = None
        # self.filtered_df = None

        # dataframe filters dict
        self.df_filter_dict = {
            0: [self.ui.w_Filter_0, self.ui.le_LB_0, self.ui.le_RB_0, self.ui.cb_DFFilter_0,
                self.ui.cb_Filter_Toggle_0, self.ui.pb_Remove_Filter]}
        self.filter_count = 1

        # plot variables
        # self.x = []
        # self.y = []
        self.scatter_pareto_object = {}

        # add checkable combobox
        self.cb_Shape_Space_Keys = QCheckableComboBox()
        self.cb_Shape_Space_Keys.setMinimumWidth(75)
        self.ui.gl_Selected_Shapes.addWidget(self.cb_Shape_Space_Keys)

        # shape space initialization
        self._shape_space = {}

        self.code_change_control()

        # sequential vs parallel
        self.run_mode_control()

        # interactive plot
        # self.cid_click = self.ppplot.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cid_pick = self.ppplot.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.ppplot.fig.canvas.mpl_connect("motion_notify_event", self.hover)

        self.ui.w_Constraints.setVisible(self.ui.pb_Constraints.isChecked())
        self.ui.w_Save_Sub_Data.setVisible(self.ui.pb_Save_Sub_Data.isChecked())

    def run_mode_control(self):
        if self.ui.cb_Run_Mode.currentText() == "Sequential":
            self.ui.w_No_Of_Processors.setMinimumWidth(0)
            self.ui.w_No_Of_Processors.setMaximumWidth(0)
        else:
            self.ui.w_No_Of_Processors.setMinimumWidth(0)
            self.ui.w_No_Of_Processors.setMaximumWidth(0)
            animate_width(self.ui.w_No_Of_Processors, 0, 300, True)

    def code_change_control(self):
        if self.ui.cb_Code.currentText() == 'ABCI':
            self.ui.w_Intervals.show()
            self.ui.w_SLANS.hide()

            # set default post processing folder to project directory
            self.ui.le_SimulationData_Folder.setText(fr"{self.main_control.projectDir}\SimulationData\ABCI")
        else:
            self.ui.w_Intervals.hide()
            self.ui.w_SLANS.show()

            # set default post processing folder to project directory
            self.ui.le_SimulationData_Folder.setText(fr"{self.main_control.projectDir}\SimulationData\SLANS")

    def process_ABCI_data(self, shape_space_dir, abci_data_dir, folder_list, request, save_excel=None):
        if not os.path.exists(f'{save_excel}.xlsx'):
            d = fr.json_reader(shape_space_dir)
            # print(d.loc['MC'])
            # shape_var = d.loc['MC']

            l, u = folder_list[0], folder_list[1]

            A, B, a, b, Ri, L, Req = [], [], [], [], [], [], []
            key_list = []
            for key, shape in d.items():
                if l <= key <= u:
                    key_list.append(key)
                    A.append(shape['MC'][0])
                    B.append(shape['MC'][1])
                    a.append(shape['MC'][2])
                    b.append(shape['MC'][3])
                    Ri.append(shape['MC'][4])
                    L.append(shape['MC'][5])
                    Req.append(shape['MC'][6])

            # data arrays initialization
            k_loss_array_transverse = []
            k_loss_array_longitudinal = []
            k_loss_M0 = []
            df_M0D1_array = []
            Zmax_mon_list = []
            Zmax_dip_list = []
            problem_keys = []

            def del_from_list(indx):
                key_list.pop(indx)
                A.pop(indx)
                B.pop(indx)
                a.pop(indx)
                b.pop(indx)
                Ri.pop(indx)
                L.pop(indx)
                Req.pop(indx)

            def delete_problem_keys(p_keys):
                # delete problematic keys
                # remove duplicates
                p_keys = list(dict.fromkeys(p_keys))
                print("After making unique", p_keys)
                for indx in p_keys:
                    del_from_list(indx)

            def calc_k_loss():
                for i in range(l, u + 1):
                    try:
                        print(f"Processing for Cavity {i}")
                        abci_data_long = ABCIData(abci_data_dir, i, 0)
                        abci_data_trans = ABCIData(abci_data_dir, i, 1)

                        # trans
                        x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                        k_loss_trans = abci_data_trans.loss_factor['Transverse']

                        if math.isnan(k_loss_trans):
                            # pop key from lists
                            problem_keys.append(i - l)
                            print_(f"Encountered an exception: Check shape {i}")
                            # problem_keys.append(i-l)
                            continue

                        # long
                        x, y, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                        abci_data_long.get_data('Loss Factor Spectrum Integrated up to F')

                        k_M0 = abci_data_long.y_peaks[0]
                        k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
                        k_loss_HOM = k_loss_long - k_M0

                        # append only after successful run
                        k_loss_M0.append(k_M0)
                        k_loss_array_longitudinal.append(k_loss_HOM)
                        k_loss_array_transverse.append(k_loss_trans)

                    except Exception as e:
                        problem_keys.append(i - l)
                        print_(f"Encountered an exception: Check shape {i} -> {e}")
                        continue

            def calc_df_M0D1():
                for i in range(l, u + 1):
                    try:
                        print(f"Processing for Cavity {i}")
                        abci_data_mon = ABCIData(abci_data_dir, i, 0)
                        abci_data_dip = ABCIData(abci_data_dir, i, 1)

                        abci_data_mon.get_data('Real Part of Longitudinal Impedance')
                        abci_data_dip.get_data('Real Part of Transverse Impedance')

                        df = (abci_data_dip.x_peaks[0] - abci_data_mon.x_peaks[0]) * 1e3  # convert to MHz
                        df_M0D1_array.append(df)
                        # print(df)
                    except Exception as e:
                        print_(f"Encountered an exception: Check shape {i} -> {e}")
                        continue

            def get_Zmax(mask_mon=None, mask_dip=None, plot=False):
                if mask_mon is None:
                    mask_mon = [0.0, 1]
                if mask_dip is None:
                    mask_dip = [0.0, 1]

                for i in range(l, u + 1):
                    try:
                        print(f"Processing for Cavity {i}")
                        abci_data_mon = ABCIData(abci_data_dir, i, 0)
                        abci_data_dip = ABCIData(abci_data_dir, i, 1)

                        # get longitudinal and transverse impedance plot data
                        xr_mon, yr_mon, _ = abci_data_mon.get_data('Real Part of Longitudinal Impedance')
                        xi_mon, yi_mon, _ = abci_data_mon.get_data('Imaginary Part of Longitudinal Impedance')

                        xr_dip, yr_dip, _ = abci_data_dip.get_data('Real Part of Transverse Impedance')
                        xi_dip, yi_dip, _ = abci_data_dip.get_data('Imaginary Part of Transverse Impedance')

                        # calculate magnitude
                        ymag_mon = (yr_mon ** 2 + yi_mon ** 2) ** 0.5
                        ymag_dip = (yr_dip ** 2 + yi_dip ** 2) ** 0.5

                        # get peaks
                        peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
                        xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

                        peaks_dip, _ = sps.find_peaks(ymag_dip, height=0)
                        xp_dip, yp_dip = np.array(xr_dip)[peaks_dip], np.array(ymag_dip)[peaks_dip]

                        # get mask
                        msk_mon = [(mask_mon[0] < x < mask_mon[1]) for x in xp_mon]
                        msk_dip = [(mask_dip[0] < x < mask_dip[1]) for x in xp_dip]

                        # print(max(yp_mon[msk_mon]), max(yp_dip[msk_dip]))

                        if len(yp_mon[msk_mon]) != 0 and len(yp_dip[msk_dip]) != 0:
                            Zmax_mon = max(yp_mon[msk_mon])
                            Zmax_dip = max(yp_dip[msk_dip])

                            Zmax_mon_list.append(Zmax_mon)
                            Zmax_dip_list.append(Zmax_dip)
                        elif len(yp_mon) != 0 and len(yp_dip) != 0:
                            Zmax_mon_list.append(0)
                            Zmax_dip_list.append(0)
                        else:
                            problem_keys.append(i - l)

                    except Exception as e:
                        problem_keys.append(i - l)
                        print(f"Encountered an exception: Check shape {i} -> {e}")

            def all():
                for i in range(l, u + 1):
                    try:
                        print(f"Processing for Cavity {i}")
                        abci_data_long = ABCIData(abci_data_dir, i, 0)
                        abci_data_trans = ABCIData(abci_data_dir, i, 1)

                        # trans
                        x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                        k_loss_trans = abci_data_trans.loss_factor['Transverse']

                        if math.isnan(k_loss_trans):
                            # pop key from lists
                            problem_keys.append(i - l)
                            print_(f"Encountered an exception: Check shape {i}")
                            problem_keys.append(i - l)
                            continue

                        # long
                        x, y, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                        abci_data_long.get_data('Loss Factor Spectrum Integrated up to F')

                        k_M0 = abci_data_long.y_peaks[0]
                        k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
                        k_loss_HOM = k_loss_long - k_M0

                        # dfM0D1
                        df = (abci_data_trans.x_peaks[0] - abci_data_long.x_peaks[0]) * 1e3  # convert to MHz

                        # append only after successful run
                        k_loss_M0.append(k_M0)
                        k_loss_array_longitudinal.append(k_loss_HOM)
                        k_loss_array_transverse.append(k_loss_trans)
                        df_M0D1_array.append(df)

                    except Exception as e:
                        print_(f"Encountered an exception: Check shape {i} -> {e}")
                        continue

            if request == 'k_loss':
                calc_k_loss()
                delete_problem_keys(problem_keys)
                # save excel
                if save_excel:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
                            'k_loss_M0': k_loss_M0, 'k_loss_long': k_loss_array_longitudinal,
                            'k_loss_trans': k_loss_array_transverse}
                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', sheet_name='Sheet_1')

            if request == 'df_M0D1':
                calc_df_M0D1()
                delete_problem_keys(problem_keys)
                # save excel
                if save_excel:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
                            'df_M0D1': df_M0D1_array}
                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', sheet_name='Sheet_1')

            if request == 'Zmax':
                mon_mask, dip_mask = [0.73, 0.77], [0.6, 1]
                get_Zmax(mon_mask, dip_mask)
                problem_keys_sorted = sorted(problem_keys,
                                             reverse=True)  # this is so as to start deleting from the largest key. If this is not done, the index is messed up
                # print(problem_keys_sorted)
                delete_problem_keys(problem_keys_sorted)
                # print(len(Zmax_mon_list), len(Zmax_dip_list), len(A))
                # save excel
                if save_excel:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
                            f'Z_long[max({mon_mask[0]}<f<{mon_mask[1]})]': Zmax_mon_list,
                            f'Z_trans[max({dip_mask[0]}<f<{dip_mask[1]})]': Zmax_dip_list}

                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', sheet_name='Sheet_1')

            if request == 'all':
                all()
                # print(len(k_loss_array_longitudinal), len(k_loss_array_transverse), len(df_M0D1_array), len(A))
                problem_keys_sorted = sorted(problem_keys,
                                             reverse=True)  # this is so as to start deleting from the largest key. If this is not done, the index is messed up
                # print(problem_keys_sorted)
                delete_problem_keys(problem_keys_sorted)
                # print(len(k_loss_array_longitudinal), len(k_loss_array_transverse), len(df_M0D1_array), len(A))

                # save excel
                if save_excel:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
                            'k_loss_M0': k_loss_M0, 'k_loss_long': k_loss_array_longitudinal,
                            'k_loss_trans': k_loss_array_transverse, 'df_M0D1': df_M0D1_array}
                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx')

        else:
            print(f"Hey Chief, seems you've already processed the {save_excel} data for this folder. "
                  "Delete the file to reprocess or rename save_excel argument")

    def process_SLANS_data(self, slans_data_dir, mode, bc, request, save_excel, parallel=False):
        reply = "Yes"
        if os.path.exists(f"{slans_data_dir}\geometric_parameters.json"):
            with open(f"{slans_data_dir}\geometric_parameters.json", 'r') as f:
                d = json.load(f)

    def combine_files(self, save_excel):
        if len(list(self.combine_dict_dict.keys())) == 2:
            f1 = self.combine_dict_dict[0]
            f2 = self.combine_dict_dict[1]

            f3 = f1[self.column_list_1].merge(f2[self.column_list_2], on=self.ui.cb_On.currentText(),
                                              how=self.ui.cb_How.currentText())
            if save_excel.split('.')[-1] != 'xlsx':
                save_excel = f'{save_excel}.xlsx'

            f3.to_excel(save_excel, index=False)

    def combine_files_parallel_run(self, proc_count, folder, save_excel):
        if save_excel != '' and folder != '':
            try:
                df = pd.DataFrame()
                for p in range(0, proc_count):
                    d = fr'{folder}\Proc_{p}.xlsx'
                    df = df.append(pd.read_excel(d), ignore_index=True)

                if save_excel.split('.')[-1] != 'xlsx':
                    save_excel = f'{save_excel}.xlsx'

                df.to_excel(save_excel, index=False)
            except:
                print("Please enter a valid file and/or folder name.")
        else:
            print("Please enter file and/or folder name.")

    def load_dir(self, dirr, lw, selected_columns, l, list_key):
        # remove all items from list and widget if any
        self.remove_all(lw, selected_columns, l)

        d = pd.read_excel(dirr)
        keys = d.keys()

        for k in keys:
            lw.addItem(k)

        self.combine_dict_dict[list_key] = d

    def add_all(self, lw_Doc, selected_columns_widget, l):
        itemsTextList = [str(lw_Doc.item(i).text()) for i in range(lw_Doc.count())]

        for txt in itemsTextList:
            # add to combine list widget
            selected_columns_widget.addItem(txt)

            # add to merge list
            l.append(txt)

        # check for common key
        self.check_for_common_key()

    def remove_all(self, lw, selected_columns, l):
        # remove all items from list and widget if any
        lw.clear()
        l.clear()
        selected_columns.clear()

        # check for common key
        self.check_for_common_key()

    def add_to_list_widget(self, lw_Doc, selected_columns_widget, l):

        # get text from list widget
        txt = lw_Doc.currentItem().text()

        # add to combine list widget
        selected_columns_widget.addItem(txt)

        # add to merge list
        l.append(txt)

        # check for common key
        self.check_for_common_key()

    def remove_from_list_widget(self, selected_columns_widget, l):
        # get current item selected
        item = selected_columns_widget.currentItem()
        txt = item.text()

        # remove from combine list widget
        selected_columns_widget.takeItem(selected_columns_widget.row(item))

        # remove from merge list
        l.remove(txt)

        # check for common key
        self.check_for_common_key()

    def check_for_common_key(self):
        lw1 = self.ui.lw_Selected_Columns_Doc1
        lw2 = self.ui.lw_Selected_Columns_Doc2

        # get item text list from list widgets
        itemsTextList1 = [str(lw1.item(i).text()) for i in range(lw1.count())]
        itemsTextList2 = [str(lw2.item(i).text()) for i in range(lw2.count())]

        # check for intersection in itemTextList
        intersection_list = list(set(itemsTextList1) & set(itemsTextList2))

        # check if intersection_list empty, populate cb_On and disable combine button if yes
        if intersection_list:
            # clear cb_On
            self.ui.cb_On.clear()
            # populate cb_On
            for key in intersection_list:
                self.ui.cb_On.addItem(key)

            # enable combine button
            self.ui.pb_Combine.setEnabled(True)
        else:
            self.ui.pb_Combine.setEnabled(False)

    def show_hide_(self, wid1, wid2):
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()

    def toggle_widgets(self, key):
        for k, wid in self.w_dict.items():
            if k == key:
                wid[1].show()
            else:
                wid[1].hide()

    def toggle_page(self, key):
        if key == 'CEM':
            self.ui.stackedWidget.setCurrentIndex(0)

        if key == 'CST':
            self.ui.stackedWidget.setCurrentIndex(1)

    def load_file(self):
        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Excel Files (*.xlsx)")
        self.ui.le_Pandas_Filename.setText(filename)

        self.dataset_filename = filename
        try:
            if filename not in list(self.datasets.keys()):
                print("loading dict")
                df = pd.read_excel(filename)
                self.datasets[filename] = df
                pandas_model = PandasModel(df)
                self.pandas_model_dict[filename] = pandas_model

                # initialise pandas object
                pandas_object = PandasTV(self)
                pandasUI = pandas_object.pandasUI

                pandasUI.tv_Pandas.setModel(pandas_model)

                if self.MDI_DATAFRAME_COUNT == 0:
                    self.pandas_object_dict[filename] = pandas_object

                    # show in mdi area
                    sub = QMdiSubWindow()
                    sub.setWidget(pandas_object.w_Pandas)
                    sub.setWindowTitle(f'Dataframe: {filename}')
                    self.ui.mdiArea.addSubWindow(sub)
                    sub.show()
                    self.sub_window_dict[filename] = sub
                    # self.MDI_DATAFRAME_COUNT += 1

                # populate first objective functions combo bozed and filter
                for val in self.df_filter_dict.values():
                    self.populate_filter(self.ui.ccb_Objective_Function_F1, df)
                    self.populate_filter(self.ui.ccb_Objective_Function_F2, df)
                    self.populate_filter(self.ui.ccb_Objective_Function_F3, df)
                    self.populate_filter(val[3], df)

                # print("datasets: ", self.datasets)

        except:
            print('Please select file.')

    def add_new_dataset(self):
        # create widget
        w = QWidget()
        gl = QGridLayout()
        w.setLayout(gl)

        # create buttons
        le = QLineEdit(self.dataset_filename)
        le.setMaximumWidth(150)
        gl.addWidget(le, 0, 0, 1, 1)
        cb = QCheckBox()
        cb.clicked.connect(lambda: print("checkbox clicked"))
        gl.addWidget(cb, 0, 1, 1, 1)
        pb = QPushButton('')
        pb.clicked.connect(lambda: print('Cancel clicked'))
        gl.addWidget(pb, 0, 2, 1, 1)

        self.ui.gl_Datasets_Widgets.addWidget(w)

    def process_objective_function_input(self, key, df, exp, norm, weight):
        if df is not None:
            operands = []
            operands_global = []
            operands_inv = []
            for i, e in enumerate(exp):
                # compute global pareto front
                if df is not None:
                    operands_global.append(df[e] * weight[i] / norm[i])

                # compute local pareto front
                if len(self.filtered_df_dict.keys()) != 0:
                    operands.append(self.filtered_df_dict[key][e] * weight[i] / norm[i])
                    operands_inv.append(self.filtered_df_inverse_dict[key][e] * weight[i] / norm[i])
                else:
                    operands.append(df[e] * weight[i] / norm[i])

            var = sum(operands)
            var_inv = sum(operands_inv)
            var_global = sum(operands_global)

            if isinstance(var_inv, int):
                return var.tolist(), [], var_global.tolist()
            else:
                return var.tolist(), var_inv.tolist(), var_global.tolist()
        else:
            return [], [], []

    def update_F_Value_Dict(self, d, ccb, tw):
        f = ccb.currentText().split(', ')
        for i, x in enumerate(f):
            d[x] = {'norm': tw.cellWidget(i, 1).value(), 'weight': tw.cellWidget(i, 2).value()}

    def change_alpha(self, alpha):
        for key in self.datasets.keys():
            self.scatter_plot_object_inv[key][0].set_alpha(alpha)

        self.ppplot.fig.canvas.draw_idle()

    def plot_objective_function(self, alpha=1.0):
        # get operands
        exp_x = self.ui.ccb_Objective_Function_F1.currentText().split(', ')
        exp_y = self.ui.ccb_Objective_Function_F2.currentText().split(', ')

        x_norm, y_norm, x_weight, y_weight = [], [], [], []
        if exp_x != [''] and exp_y != ['']:
            # get norms
            for i in range(len(exp_x)):
                x_norm.append(self.ui.tw_Objective_Function_F1.cellWidget(i, 1).value())
                x_weight.append(self.ui.tw_Objective_Function_F1.cellWidget(i, 2).value())

            for i in range(len(exp_y)):
                y_norm.append(self.ui.tw_Objective_Function_F2.cellWidget(i, 1).value())
                y_weight.append(self.ui.tw_Objective_Function_F2.cellWidget(i, 2).value())

            # check if norm count equal operand count
            if len(exp_x) == len(x_norm) and len(exp_y) == len(y_norm):
                # self.ppplot.ax.cla()
                for key, df in self.datasets.items():
                    x, x_inv, x_global = self.process_objective_function_input(key, df, exp_x, x_norm, x_weight)
                    self.x_dict[key], self.x_inv_dict[key], self.x_global_dict[key] = x, x_inv, x_global
                    y, y_inv, y_global = self.process_objective_function_input(key, df, exp_y, y_norm, y_weight)
                    self.y_inv_dict[key], self.y_inv_dict[key], self.y_global_dict = y, y_inv, y_global

                    # self.scatter_plot_object_inv = self.ppplot.ax.scatter(self.x_inv, self.y_inv, edgecolors='k', s=150,
                    #                                                       picker=True, alpha=alpha, label="Cavity Geometry")
                    # self.scatter_plot_object = self.ppplot.ax.scatter(self.x, self.y, edgecolors='k', s=150,
                    #                                                   picker=True, label="Cavity Geometry (Filtered)")

                    if key in self.scatter_plot_object_inv.keys():
                        self.scatter_plot_object_inv[key][0].set_xdata(x_inv)
                        self.scatter_plot_object_inv[key][0].set_ydata(y_inv)
                    else:
                        self.scatter_plot_object_inv[key] = self.ppplot.ax.plot(x_inv, y_inv, linestyle='None', marker='o', markersize=10.0, markeredgecolor="black",
                                                                          picker=True, alpha=alpha, label=f"{key.split('.')[-1]} Cavity Geometry")

                    if self.scatter_plot_object[key]:
                        self.scatter_plot_object[key][0].set_xdata(x)
                        self.scatter_plot_object[key][0].set_ydata(y)
                    else:
                        if f"{type(self.ppplot.ax)}" == "<class 'matplotlib.axes._subplots.Axes3DSubplot'>":
                            self.ppplot.ax.remove()
                            self.ppplot.ax = self.ppplot.fig.add_subplot()

                        self.scatter_plot_object[key] = self.ppplot.ax.plot(x, y, linestyle='None', marker='o', markersize=10.0, markeredgecolor="black",
                                                                   picker=True, label=f"{key.split('.')[-1]} Cavity Geometry (Filtered)")

                    # add annotation
                    self.annot = self.ppplot.ax.annotate("", xy=(0, 0), xytext=(-20, -10), textcoords="offset points",
                                                         bbox=dict(boxstyle="round", fc="w"),
                                                         arrowprops=dict(arrowstyle="->"))

                    data = x, y, x_inv, y_inv, x_global, y_global
                    # plot pareto
                    self.plot_pareto(key, data)

                # axis decorations
                self.ppplot.ax.set_xlabel("$F_1$")
                self.ppplot.ax.set_ylabel("$F_2$")

                self.ppplot.ax.relim()
                self.ppplot.ax.autoscale()
                self.leg = self.ppplot.ax.legend()
                self.ppplot.fig.canvas.draw_idle()

                # show MDI window. Just one for multiple plots. Maybe change later if need arises
                if self.MDI_PLOT_COUNT == 0:
                    # show in mdi area
                    sub = QMdiSubWindow()
                    sub.setWidget(self.ppplot.w_PPPlot)
                    sub.setWindowTitle('Plot')
                    self.ui.mdiArea.addSubWindow(sub)
                    sub.show()
                    self.MDI_PLOT_COUNT += 1
            else:
                print("Operand and Norm count are not equal. Please check.")
        else:
            pass

    def plot_objective_function_3D(self, alpha=1.0):
        # get operands
        exp_x = self.ui.ccb_Objective_Function_F1.currentText().split(', ')
        exp_y = self.ui.ccb_Objective_Function_F2.currentText().split(', ')
        exp_z = self.ui.ccb_Objective_Function_F3.currentText().split(', ')

        x_norm, y_norm, x_weight, y_weight, z_norm, z_weight = [], [], [], [], [], []
        if exp_x != [''] and exp_y != [''] and exp_z != ['']:
            # get norms
            for i in range(len(exp_x)):
                x_norm.append(self.ui.tw_Objective_Function_F1.cellWidget(i, 1).value())
                x_weight.append(self.ui.tw_Objective_Function_F1.cellWidget(i, 2).value())

            for i in range(len(exp_y)):
                y_norm.append(self.ui.tw_Objective_Function_F2.cellWidget(i, 1).value())
                y_weight.append(self.ui.tw_Objective_Function_F2.cellWidget(i, 2).value())

            for i in range(len(exp_z)):
                z_norm.append(self.ui.tw_Objective_Function_F3.cellWidget(i, 1).value())
                z_weight.append(self.ui.tw_Objective_Function_F3.cellWidget(i, 2).value())

            # check if norm count equal operand count
            if len(exp_x) == len(x_norm) and len(exp_y) == len(y_norm) and len(exp_z) == len(z_norm):
                for key, df in self.datasets.items():
                    # self.ppplot.ax.cla()
                    x, x_inv, x_global = self.process_objective_function_input(key, df, exp_x, x_norm, x_weight)
                    y, y_inv, y_global = self.process_objective_function_input(key, df, exp_y, y_norm, y_weight)
                    z, z_inv, z_global = self.process_objective_function_input(key, df, exp_z, z_norm, z_weight)
                    self.x_dict[key], self.x_inv_dict[key], self.x_global_dict[key] = x, x_inv, x_global
                    self.y_dict[key], self.y_inv_dict[key], self.y_global_dict[key] = y, y_inv, y_global
                    self.z_dict[key], self.z_inv_dict[key], self.z_global_dict[key] = z, z_inv, z_global

                    # self.scatter_plot_object_inv = self.ppplot.ax.scatter(self.x_inv, self.y_inv, edgecolors='k', s=150,
                    #                                                       picker=True, alpha=alpha, label="Cavity Geometry")
                    # self.scatter_plot_object = self.ppplot.ax.scatter(self.x, self.y, edgecolors='k', s=150,
                    #                                                   picker=True, label="Cavity Geometry (Filtered)")
                    # print(self.scatter_plot_object_inv)

                    if key in self.scatter_plot_object_inv.keys():
                        self.scatter_plot_object_inv[key][0].set_xdata(x_inv)
                        self.scatter_plot_object_inv[key][0].set_ydata(y_inv)
                        self.scatter_plot_object_inv[key][0].set_zdata(z_inv)
                    else:
                        print(type(self.ppplot.ax))
                        if f"{type(self.ppplot.ax)}" == "<class 'matplotlib.axes._subplots.AxesSubplot'>":
                            self.ppplot.ax.remove()
                            self.ppplot.ax = self.ppplot.fig.add_subplot(projection='3d')

                        # if isinstance(self.ppplot.ax, mpl.axes._subplots.Axes3DSubplot):

                        self.scatter_plot_object_inv[key] = self.ppplot.ax.plot(x_inv, y_inv, z_inv, linestyle='None',
                                                                           marker='o', markersize=4,
                                                                           markeredgecolor="black",
                                                                           picker=True, alpha=alpha,
                                                                           label=f"{key.split('/')[-1]} Cavity Geometry")

                    if key in self.scatter_plot_object.keys():
                        self.scatter_plot_object[key][0].set_xdata(x)
                        self.scatter_plot_object[key][0].set_ydata(y)
                        self.scatter_plot_object[key][0].set_ydata(z)
                    else:
                        self.scatter_plot_object[key] = self.ppplot.ax.plot(x, y, z, linestyle='None', marker='o',
                                                                       markersize=4, markeredgecolor="black",
                                                                       picker=True, label=f"{key.split('/')[-1]} Cavity Geometry (Filtered)")
                        print("check key:", key, self.scatter_plot_object)

                    # add annotation
                    self.annot_dict[key] = self.ppplot.ax.annotate("", xy=(0, 0), xytext=(-20, -10), textcoords="offset points",
                                                         bbox=dict(boxstyle="round", fc="w"),
                                                         arrowprops=dict(arrowstyle="->"))

                    data = x, y, z, x_inv, y_inv, z_inv, x_global, y_global, z_global
                    # plot pareto
                    self.plot_pareto_3D(key, data)

                # axis decorations
                self.ppplot.ax.set_xlabel("$F_1$")
                self.ppplot.ax.set_ylabel("$F_2$")
                self.ppplot.ax.set_zlabel("$F_3$")

                self.ppplot.ax.relim()
                self.ppplot.ax.autoscale()
                self.ppplot.ax.set_xlim([min(x_global), max(x_global)])
                self.ppplot.ax.set_ylim([min(y_global), max(y_global)])
                self.ppplot.ax.set_zlim([min(z_global), max(z_global)])
                self.leg = self.ppplot.ax.legend()
                self.ppplot.fig.canvas.draw_idle()

                if self.MDI_PLOT_COUNT == 0:
                    # show in mdi area
                    sub = QMdiSubWindow()
                    sub.setWidget(self.ppplot.w_PPPlot)
                    sub.setWindowTitle('Plot')
                    self.ui.mdiArea.addSubWindow(sub)
                    sub.show()
                    self.MDI_PLOT_COUNT += 1

                # self.plot_objective_function_3D_mayavi()
            else:
                print("Operand and Norm count are not equal. Please check.")
        else:
            pass

    # def plot_objective_function_3D_mayavi(self, alpha=1.0):
    #     # get operands
    #     exp_x = self.ui.ccb_Objective_Function_F1.currentText().split(', ')
    #     exp_y = self.ui.ccb_Objective_Function_F2.currentText().split(', ')
    #     exp_z = self.ui.ccb_Objective_Function_F3.currentText().split(', ')
    #
    #     x_norm, y_norm, x_weight, y_weight, z_norm, z_weight = [], [], [], [], [], []
    #     if exp_x != [''] and exp_y != [''] and exp_z != ['']:
    #         # get norms
    #         for i in range(len(exp_x)):
    #             x_norm.append(self.ui.tw_Objective_Function_F1.cellWidget(i, 1).value())
    #             x_weight.append(self.ui.tw_Objective_Function_F1.cellWidget(i, 2).value())
    #
    #         for i in range(len(exp_y)):
    #             y_norm.append(self.ui.tw_Objective_Function_F2.cellWidget(i, 1).value())
    #             y_weight.append(self.ui.tw_Objective_Function_F2.cellWidget(i, 2).value())
    #
    #         for i in range(len(exp_z)):
    #             z_norm.append(self.ui.tw_Objective_Function_F3.cellWidget(i, 1).value())
    #             z_weight.append(self.ui.tw_Objective_Function_F3.cellWidget(i, 2).value())
    #
    #         # check if norm count equal operand count
    #         if len(exp_x) == len(x_norm) and len(exp_y) == len(y_norm) and len(exp_z) == len(z_norm):
    #             for key, df in self.datasets.items():
    #                 # self.ppplot.ax.cla()
    #                 x, x_inv, x_global = self.process_objective_function_input(key, df, exp_x, x_norm, x_weight)
    #                 y, y_inv, y_global = self.process_objective_function_input(key, df, exp_y, y_norm, y_weight)
    #                 z, z_inv, z_global = self.process_objective_function_input(key, df, exp_z, z_norm, z_weight)
    #                 self.x_dict[key], self.x_inv_dict[key], self.x_global_dict[key] = x, x_inv, x_global
    #                 self.y_dict[key], self.y_inv_dict[key], self.y_global_dict[key] = y, y_inv, y_global
    #                 self.z_dict[key], self.z_inv_dict[key], self.z_global_dict[key] = z, z_inv, z_global
    #
    #                 data = x, y, z, x_inv, y_inv, z_inv, x_global, y_global, z_global
    #                 # plot pareto
    #                 self.plot_pareto_3D_mayavi(key, data)
    #                 mlab.show()
    #
    #                 # # axis decorations
    #                 # self.ppplot.ax.set_xlabel("$F_1$")
    #                 # self.ppplot.ax.set_ylabel("$F_2$")
    #                 # self.ppplot.ax.set_zlabel("$F_3$")
    #                 #
    #                 # self.ppplot.ax.relim()
    #                 # self.ppplot.ax.autoscale()
    #                 # self.ppplot.ax.set_xlim([min(self.x_global), max(self.x_global)])
    #                 # self.ppplot.ax.set_ylim([min(self.y_global), max(self.y_global)])
    #                 # self.ppplot.ax.set_zlim([min(self.z_global), max(self.z_global)])
    #                 # self.leg = self.ppplot.ax.legend()
    #                 # self.ppplot.fig.canvas.draw_idle()
    #                 #
    #                 # if self.MDI_PLOT_COUNT == 0:
    #                 #     # show in mdi area
    #                 #     sub = QMdiSubWindow()
    #                 #     sub.setWidget(self.ppplot.w_PPPlot)
    #                 #     sub.setWindowTitle('Plot')
    #                 #     self.ui.mdiArea.addSubWindow(sub)
    #                 #     sub.show()
    #                 #     self.MDI_PLOT_COUNT += 1
    #         else:
    #             print("Operand and Norm count are not equal. Please check.")
    #     else:
    #         pass

    def populate_objective_function(self, f, ccb, tw, d):
        tw.setRowCount(len(f))  # and one row in the table

        for i, x in enumerate(f):
            label = QLabel(x)
            tw.setCellWidget(i, 0, label)

            dsb = QDoubleSpinBox()
            dsb2 = QDoubleSpinBox()
            try:
                dsb.setValue(d[x]['norm'])
                dsb2.setValue(d[x]['weight'])
            except:
                dsb.setValue(1)
                dsb2.setValue(1)

            tw.setCellWidget(i, 1, dsb)
            tw.setCellWidget(i, 2, dsb2)
            dsb.editingFinished.connect(lambda: self.update_F_Value_Dict(d, ccb, tw))

    def plot_pareto(self, key, data):
        x, y, x_inv, y_inv, x_global, y_global = data

        if self.ui.cb_Pareto.checkState() == 2:
            if self.ui.cb_Pareto_Setting.currentText().lower() == "global":
                x_pareto, y_pareto = self.pareto_front(x_global, y_global)
            else:
                x_pareto, y_pareto = self.pareto_front(x, y)

            y_pareto = [y for _, y in sorted(zip(x_pareto, y_pareto))]
            x_pareto.sort()

            if key in self.scatter_pareto_object.keys():
                self.scatter_pareto_object[key][0].set_xdata(x_pareto)
                self.scatter_pareto_object[key][0].set_ydata(y_pareto)
            else:
                self.scatter_pareto_object[key] = self.ppplot.ax.plot(x_pareto, y_pareto, c='r', marker='o', lw=2,
                                                                 markersize=6, mec='k',
                                                                 picker=True, label="Pareto Front")

            if not self.ui.sc_Toggle_3D.isChecked():
                self.ppplot.plt.toggle_ax(True)

                # print pareto shapes
                pareto_list = []
                plotted_pareto_pts = []
                for x, y in zip(x_pareto, y_pareto):
                    if self.ui.cb_Pareto_Setting.currentText().lower() == 'global':
                        xvalues = x_global
                        yvalues = y_global
                    else:
                        xvalues = self.scatter_plot_object[key][0].get_xdata()
                        yvalues = self.scatter_plot_object[key][0].get_ydata()

                    ind = (np.where(xvalues == x) and np.where(yvalues == y))[0][0]
                    plotted_pareto_pts.append((x, y))
                    try:
                        k = self.filtered_df_dict[key].loc[ind][0:8]
                    except (AttributeError, KeyError):
                        k = self.datasets[key].loc[ind][0:8]

                    pareto_list.append(k.tolist())
                subax_count = self.ppplot.plt.n*2-1
                for i in range(subax_count):
                    if i == subax_count-1:
                        id_to_plot = -1
                    else:
                        id_to_plot = int(len(pareto_list)/(subax_count-1)*i)

                    #uncomment to plot cavities
                    # self.plot_cavity(pareto_list[id_to_plot][0], pareto_list[id_to_plot][1:], self.ppplot.plt.plot_list[i])

                    # add arrows
                    # clear old arrows
                    # print(self.ppplot.ax.patches, type(self.ppplot.ax.patches))

                    if i < (subax_count-1)/2:
                        self.add_cross_arrow(self.ppplot.ax, plotted_pareto_pts[id_to_plot], self.ppplot.plt.plot_list[i], (188, 388/2))
                    elif i == (subax_count-1)/2:
                        self.add_cross_arrow(self.ppplot.ax, plotted_pareto_pts[id_to_plot], self.ppplot.plt.plot_list[i], (188, 388))
                    else:
                        self.add_cross_arrow(self.ppplot.ax, plotted_pareto_pts[id_to_plot], self.ppplot.plt.plot_list[i], (0, 388))

                    self.ppplot.fig.tight_layout()

        else:
            for a in self.arrow_patch_list:
                self.ppplot.fig.patches.remove(a)
            self.arrow_patch_list = []

            if self.scatter_pareto_object:
                self.scatter_pareto_object[0].remove()
                self.scatter_pareto_object = None

            self.ppplot.plt.toggle_ax(False)

        self.ppplot.fig.canvas.draw_idle()

    def plot_pareto_3D(self, key, data):
        x, y, z, x_inv, y_inv, z_inv, x_global, y_global, z_global = data
        # if self.ui.cb_Pareto.checkState() == 2:
        if True:
            if self.ui.cb_Pareto_Setting.currentText().lower() == "global":
                x_pareto, y_pareto, z_pareto = self.pareto_front_3D(x_global, y_global, z_global)
            else:
                x_pareto, y_pareto, z_pareto = self.pareto_front_3D(x, y, z)

            y_pareto = [y for _, y in sorted(zip(x_pareto, y_pareto))]
            z_pareto = [z for _, z in sorted(zip(x_pareto, z_pareto))]
            x_pareto.sort()

            # ppppp = pd.DataFrame(columns=["Epk/Eacc", "Bpk/Eacc", "R/Q"])
            # ppppp["Epk/Eacc"] = x_pareto
            # ppppp["Bpk/Eacc"] = y_pareto
            # ppppp["R/Q"] = z_pareto
            # ic(x_pareto)
            # ic(y_pareto)
            # ic(z_pareto)
            # ic(ppppp)

            if key in self.scatter_pareto_object.keys():
                self.scatter_pareto_object[key][0].set_xdata(x_pareto)
                self.scatter_pareto_object[key][0].set_ydata(y_pareto)
                self.scatter_pareto_object[key][0].set_zdata(z_pareto)
            else:

                xi = np.linspace(min(x_pareto), max(x_pareto), 500)
                yi = np.linspace(min(y_pareto), max(y_pareto), 500)
                # VERY IMPORTANT, to tell matplotlib how is your data organized
                zi = griddata((x_pareto, y_pareto), z_pareto, (xi[None, :], yi[:, None]), method='linear')
                xig, yig = np.meshgrid(xi, yi)

                self.scatter_pareto_object[key] = self.ppplot.ax.plot(x_pareto, y_pareto, z_pareto, c='k', marker='o', linestyle='None', markersize=6, mec='k',
                                                  picker=True, label="Pareto Front", zorder=1)
                # surf = self.ppplot.ax.plot_trisurf(x_pareto, y_pareto, z_pareto, cmap='jet', linewidth=0, zorder=0,
                #                                     antialiased=True,)
                surf = self.ppplot.ax.plot_surface(xig, yig, zi, linewidth=0, cmap='jet', antialiased=True,)

            # trial
            import plotly.graph_objects as go
            import plotly.express as px
            df = pd.DataFrame(list(zip(x_pareto, y_pareto, z_pareto)), columns=['Epk', 'Bpk', 'RQ'])
            fig = px.scatter_3d(df, x='Epk', y='Bpk', z='RQ')#, color='species'
            # fig = go.Figure(data=[go.Scatter3d(z=z_pareto, x=x_pareto, y=y_pareto)])
            # fig.update_layout(title='Mt Bruno Elevation', autosize=False,
            #                   width=500, height=500,
            #                   margin=dict(l=65, r=50, b=65, t=90))
            fig.show()

            print("Hereee")

            # print pareto shapes
            pareto_list = []
            plotted_pareto_pts = []
            for x, y, z in zip(x_pareto, y_pareto, z_pareto):
                if self.ui.cb_Pareto_Setting.currentText().lower() == 'global':
                    xvalues = x_global
                    yvalues = y_global
                    zvalues = z_global
                else:
                    xvalues = self.scatter_plot_object[key][0].get_xdata()
                    yvalues = self.scatter_plot_object[key][0].get_ydata()
                    zvalues = self.scatter_plot_object[key][0].get_zdata()

                ind = (np.where(xvalues == x) and np.where(yvalues == y) and np.where(zvalues == z))[0][0]
                plotted_pareto_pts.append((x, y, z))
                # try:
                #     k = self.filtered_df_dict[key].loc[ind][0:8]
                # except (AttributeError, KeyError):
                #     k = df[key].loc[ind][0:8]
                #
                # pareto_list.append(k.tolist())

            # subax_count = self.ppplot.plt.n*2-1
            # for i in range(subax_count):
            #     if i == subax_count-1:
            #         id_to_plot = -1
            #     else:
            #         id_to_plot = int(len(pareto_list)/(subax_count-1)*i)
            #
            #     self.plot_cavity(pareto_list[id_to_plot][0], pareto_list[id_to_plot][1:], self.ppplot.plt.plot_list[i])
            #
            #     # add arrows
            #     # clear old arrows
            #     # print(self.ppplot.ax.patches, type(self.ppplot.ax.patches))
            #
            #     if i < (subax_count-1)/2:
            #         self.add_cross_arrow(self.ppplot.ax, plotted_pareto_pts[id_to_plot], self.ppplot.plt.plot_list[i], (188, 388/2))
            #     elif i == (subax_count-1)/2:
            #         self.add_cross_arrow(self.ppplot.ax, plotted_pareto_pts[id_to_plot], self.ppplot.plt.plot_list[i], (188, 388))
            #     else:
            #         self.add_cross_arrow(self.ppplot.ax, plotted_pareto_pts[id_to_plot], self.ppplot.plt.plot_list[i], (0, 388))

                # self.ppplot.fig.tight_layout()

            print("Hereee 22")
        else:
            for a in self.arrow_patch_list:
                self.ppplot.fig.patches.remove(a)
            self.arrow_patch_list = []

            if self.scatter_pareto_object:
                self.scatter_pareto_object[key][0].remove()
                self.scatter_pareto_object[key] = None

            self.ppplot.plt.toggle_ax(False)

        self.ppplot.fig.canvas.draw_idle()

    # def plot_pareto_3D_mayavi(self, key, data):
    #     x, y, z, x_inv, y_inv, z_inv, x_global, y_global, z_global = data
    #     if self.ui.cb_Pareto_Setting.currentText().lower() == "global":
    #         x_pareto, y_pareto, z_pareto = self.pareto_front_3D(x_global, y_global, z_global)
    #     else:
    #         x_pareto, y_pareto, z_pareto = self.pareto_front_3D(x, y, z)
    #
    #     y_pareto = [y for _, y in sorted(zip(x_pareto, y_pareto))]
    #     z_pareto = [z for _, z in sorted(zip(x_pareto, z_pareto))]
    #     x_pareto.sort()
    #
    #     # self.scatter_pareto_object = self.ppplot.ax.plot(x_pareto, y_pareto, z_pareto, c='k', marker='o', linestyle='None', markersize=6, mec='k',
    #     #                                                  picker=True, label="Pareto Front", zorder=1)
    #     # surf = self.ppplot.ax.plot_trisurf(x_pareto, y_pareto, z_pareto, cmap='jet', linewidth=0, zorder=0)
    #     # s2 = mlab.points3d(x_pareto, y_pareto, z_pareto, color=(1, 0, 0), mode='sphere', extent=[0, 1, 0, 1, 0, 1])
    #
    #     print("heehhwer")
    #     pts = mlab.points3d(x_pareto, y_pareto, z_pareto, z_pareto)
    #     # Triangulate based on X, Y with Delaunay 2D algorithm.
    #     # Save resulting triangulation.
    #
    #     from scipy.spatial import Delaunay
    #     interp = LinearNDInterpolator(list(zip(x_pareto, y_pareto)), z_pareto)
    #     p2d = np.vstack([x_pareto, y_pareto]).T
    #     d2d = Delaunay(p2d)
    #     # Remove the point representation from the plot
    #     pts.remove()
    #
    #     # Draw a surface based on the triangulation
    #     tmesh = mlab.triangular_mesh(x_pareto, y_pareto, z_pareto, d2d.vertices, scalars=z_pareto, colormap='jet', extent=[0, 1, 0, 1, 0, 1])
    #
    #     # surf = mlab.surf(x_pareto, y_pareto, z_pareto, color=(1, 1, 0), extent=[0, 1, 0, 1, 0, 1])
    #
    #     # mlab.outline(extent=(0, 1, 0, 1, 0, 1))
    #     # mlab.axes(extent=(0, 1, 0, 1, 0, 1))
    #     mlab.axes(tmesh, ranges=[min(x_global), max(x_global),
    #                              min(y_global), max(y_global),
    #                              min(z_global), max(z_global)])
    #
    #     # print pareto shapes
    #     pareto_list = []
    #     plotted_pareto_pts = []
    #     for x, y, z in zip(x_pareto, y_pareto, z_pareto):
    #         if self.ui.cb_Pareto_Setting.currentText().lower() == 'global':
    #             xvalues = x_global
    #             yvalues = y_global
    #             zvalues = z_global
    #         else:
    #             xvalues = self.scatter_plot_object[key][0].get_xdata()
    #             yvalues = self.scatter_plot_object[key][0].get_ydata()
    #             zvalues = self.scatter_plot_object[key][0].get_zdata()
    #
    #         ind = (np.where(xvalues == x) and np.where(yvalues == y) and np.where(zvalues == z))[0][0]
    #         plotted_pareto_pts.append((x, y, z))
    #
    #         try:
    #             k = self.filtered_df_dict[key].loc[ind][0:8]
    #         except (AttributeError, KeyError):
    #             k = self.datasets[key].loc[ind][0:8]
    #
    #         pareto_list.append(k.tolist())
    #
    #     print("Hereee 22")

    def pareto_front(self, x, y):
        def reverse_list(l, goal):
            if goal == 'max':
                return l  # to find the pareto maxima
            else:
                return [-x for x in l]  # to find the pareto minima

        datapoints = np.array([reverse_list(x, self.ui.cb_Goal_F1.currentText()), reverse_list(y, self.ui.cb_Goal_F2.currentText())])

        pareto = oapackage.ParetoDoubleLong()

        for ii in range(0, datapoints.shape[1]):
            w = oapackage.doubleVector((datapoints[0, ii], datapoints[1, ii]))
            pareto.addvalue(w, ii)
        # pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs

        optimal_datapoints = datapoints[:, lst]

        return reverse_list(optimal_datapoints[0, :], self.ui.cb_Goal_F1.currentText()), reverse_list(optimal_datapoints[1, :], self.ui.cb_Goal_F2.currentText())

    def pareto_front_3D(self, x, y, z):
        def reverse_list(l, r):
            if r == 'max':
                return l
            else:
                return [-x for x in l]

        datapoints = np.array([reverse_list(x, self.ui.cb_Goal_F1.currentText()),
                               reverse_list(y, self.ui.cb_Goal_F2.currentText()),
                               reverse_list(z, self.ui.cb_Goal_F3.currentText())])

        pareto = oapackage.ParetoDoubleLong()

        for ii in range(0, datapoints.shape[1]):
            w = oapackage.doubleVector((datapoints[0, ii], datapoints[1, ii], datapoints[2, ii]))
            pareto.addvalue(w, ii)
        # pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs

        optimal_datapoints = datapoints[:, lst]

        return reverse_list(optimal_datapoints[0, :], self.ui.cb_Goal_F1.currentText()), \
               reverse_list(optimal_datapoints[1, :], self.ui.cb_Goal_F2.currentText()), \
               reverse_list(optimal_datapoints[2, :], self.ui.cb_Goal_F3.currentText())

    def plot_cavity(self, key, mid_cell, ax):
        # clear axis (improve later
        ax.cla()
        # write geometry
        mid_cell = mid_cell
        lend_cell = mid_cell
        rend_cell = mid_cell
        writeCavity(int(key), mid_cell, lend_cell, rend_cell, beampipe=[0.001, 0.001])
        data = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\C1092V\PostprocessingData\Data\{int(key)}_geom.txt", sep='\s+',
                           header=None)

        ax.plot(data[1] * 1e3, data[0] * 1e3, lw=6, label=f"C{key}")
        ax.legend(loc='lower center')

        x_label = "z [mm]"
        y_label = "r [mm]"
        # plt.xlabel(x_label)
        # plt.ylabel(y_label)
        ax.set_axis_off()
        ax.set_xlim(-94, 94)
        ax.set_ylim(-0.5, 200)
        # ax.relim()

    def add_cross_arrow(self, ax_main, ax_mani_pt, ax_sub, ax_sub_pt):
        # Create the arrow
        # 1. Get transformation operators for axis and figure
        ax0tr = ax_main.transData  # Axis 0 -> Display
        ax1tr = ax_sub.transData  # Axis 1 -> Display
        figtr = self.ppplot.fig.transFigure.inverted()  # Display -> Figure
        # 2. Transform arrow start point from axis 0 to figure coordinates
        ptB = figtr.transform(ax0tr.transform(ax_mani_pt))
        # 3. Transform arrow end point from axis 1 to figure coordinates
        ptE = figtr.transform(ax1tr.transform(ax_sub_pt))
        # 4. Create the patch
        arrow = mpl.patches.FancyArrowPatch(
            ptB, ptE, transform=self.ppplot.fig.transFigure,  # Place arrow in figure coord system
            fc="g", connectionstyle="arc3, rad=0.0", arrowstyle='simple', alpha=0.3,
            mutation_scale=20.
        )
        self.arrow_patch_list.append(arrow)
        # 5. Add patch to list of objects to draw onto the figure
        self.ppplot.fig.patches.append(arrow)

    def populate_filter(self, cb, df):
        # clear cb
        cb.clear()
        for col in df.columns:
            cb.addItem(col)

    def add_filter(self):
        # create filter widgets
        w = QWidget()
        l = QGridLayout(w)

        lb = QLineEdit()
        lb.setFixedWidth(100)
        l.addWidget(lb, 0, 0, 1, 1)
        rb = QLineEdit()
        rb.setFixedWidth(100)
        l.addWidget(rb, 0, 2, 1, 1)
        cb = QComboBox()
        l.addWidget(cb, 0, 1, 1, 1)

        chk = QCheckBox()
        l.addWidget(chk, 0, 4, 1, 1)
        # add signal
        chk.stateChanged.connect(lambda: self.apply_filter())
        chk.setFixedWidth(25)

        pb = QPushButton()
        pb.setMaximumWidth(25)
        pb.setMinimumWidth(25)
        pb.setText('X')
        pb.setStyleSheet('background-color: rgb(255, 89, 67);')
        l.addWidget(pb, 0, 5, 1, 1)
        # add signal
        pb.clicked.connect(lambda: self.remove_filter(pb))

        self.ui.gl_Filter.addWidget(w)

        # add to UI
        self.ui.gl_Filter.addWidget(w)

        # populate combobox
        if self.df is not None:
            self.populate_filter(cb, self.df)

        # add to filter dictionary
        self.df_filter_dict[self.filter_count] = [w, lb, rb, cb, chk, pb]
        self.filter_count += 1

    def remove_filter(self, pb_clicked):
        k = None
        for key, wids in self.df_filter_dict.items():
            # print(key, pb_clicked, self.df_filter_dict)
            if pb_clicked in wids:
                # delete widgets
                for wid in wids:
                    wid.deleteLater()

                # clear widget list
                wids.clear()
                k = key

        if k is not None:
            # remove filter from filter dictionary
            self.df_filter_dict.pop(k)

        self.filter_count -= 1

        # reapply filters
        self.apply_filter()

    def apply_filter(self):

        for key, df in self.datasets.items():
            # check if a dataframe already exists
            if df is not None:
                # copy main data frame
                filtered_df = df.copy(deep=True)
                filtered_df_inverse = df.copy(deep=True)

                # apply filter from to dataframe
                # check if filter dict is not empty
                if len(self.df_filter_dict.keys()) != 0:
                    for k, val in self.df_filter_dict.items():
                        if val[4].checkState() == 2:
                            lb = float(val[1].text())
                            rb = float(val[2].text())
                            k = val[3].currentText()

                            filtered_df = filtered_df[
                                (filtered_df[key] >= lb) & (filtered_df[key] <= rb)]
                            # print(self.filtered_df_inverse)

                            filtered_df = filtered_df.reset_index(drop=True)
                    self.filtered_df_dict[key] = filtered_df

                    filtered_df_inverse = filtered_df_inverse[~filtered_df_inverse.isin(filtered_df)]
                    filtered_df_inverse = filtered_df_inverse.reset_index(drop=True)

                self.filtered_df_inverse_dict[key] = filtered_df_inverse

                self.pandas_model_dict[key] = PandasModel(filtered_df)
                self.pandas_model_dict[key].pandasUI.tv_Pandas.setModel(self.pandas_model_dict[key])

        # update plot if any
        self.plot_objective_function(self.ui.hs_Alpha.value()/100)

    def add_handle(self, i, var, val):
        # combo box
        cbb = QLabel(var)

        # horizontal slider
        dsb = QDoubleSpinBox()
        dsb.setValue(val)

        # slider
        hs = QSlider(Qt.Horizontal)
        hs.setRange(0, 1)
        hs.setValue(0)
        hs.setTickInterval(100)

        # check box
        cb = QCheckBox()

        # add to widget
        self.ui.gl_Handles.addWidget(cbb, i, 0)
        self.ui.gl_Handles.addWidget(dsb, i, 1)
        self.ui.gl_Handles.addWidget(cb, i, 2)
        self.ui.gl_Handles.addWidget(hs, i, 3)

        self.OF_slider_dict[var] = [cbb, dsb, cb, hs]

    def remove_handle(self, handle):
        pass

    def onpick(self, event):
        try:
            ind = event.ind
            # print('onpick scatter:', ind, np.take(self.x, ind), np.take(self.y, ind))
            # print(self.df.loc[ind, 'A':'alpha'])
        except AttributeError:
            pass

    def update_annot(self, ind, event):
        inv = self.ppplot.ax.transData.inverted()
        # pos = self.scatter_plot_object[0].get_offsets()[ind]
        pos = [event.x, event.y]

        for key, df in self.datasets.items():
            self.annot_dict[key].xy = inv.transform(pos)
            try:
                try:
                    text = f"ti: {ind}, {self.filtered_df_dict[key].loc[ind, 'key':'key']}"
                    # text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'alpha']}"
                    # text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'B']}"
                except:
                    text = f"ti: {ind}, {df.loc[ind, 'key':'key']}"
                    # text = f"ti: {ind}, {self.df.loc[ind, 'key':'Req']}"

                self.annot_dict[key].set_text(f'{text}')
                self.annot_dict[key].get_bbox_patch().set_facecolor((167 / 255, 222 / 255, 255 / 255))
                # self.annot.get_bbox_patch().set_alpha(1)
            except KeyError:
                pass

    def update_annot_inv(self, ind, event):
        inv = self.ppplot.ax.transData.inverted()
        # pos = self.scatter_plot_object_inv[0].get_offsets()[ind]
        pos = [event.x, event.y]

        self.annot.xy = inv.transform(pos)
        for key, df in self.datasets.items():
            try:
                text = f"ti: {ind}, {self.filtered_df_inverse_dict[key].loc[ind, 'key':'key']}"
                # text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'alpha']}"
                # text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'B']}"
            except:
                text = f"ti: {ind}, {df.loc[ind, 'key':'key']}"
                # text = f"ti: {ind}, {self.df.loc[ind, 'key':'Req']}"

            self.annot_dict[key].set_text(f'{text}')
            self.annot_dict[key].get_bbox_patch().set_facecolor((167 / 255, 222 / 255, 255 / 255))
            # self.annot.get_bbox_patch().set_alpha(1)

    def hover(self, event):
        for key, df in self.datasets.items():
            vis = self.annot_dict[key].get_visible()
            # print(type(event.inaxes), type(self.ppplot.ax))
            if type(event.inaxes) == type(self.ppplot.ax):  # not so good fix
                cont, ind = self.scatter_plot_object[key][0].contains(event)
                # print(self.scatter_plot_object[0].get_ydata())
                if cont:
                    self.update_annot(ind["ind"][0], event)  # ind returns an array of close points. ind['ind'][0] returns just the first point
                    self.annot_dict[key].set_visible(True)
                    self.ppplot.fig.canvas.draw_idle()
                else:
                    cont, ind = self.scatter_plot_object_inv[key][0].contains(event)
                    if cont:
                        self.update_annot_inv(ind["ind"][0], event)  # ind returns an array of close points. ind['ind'][0] returns just the first point
                        self.annot_dict[key].set_visible(True)
                        self.ppplot.fig.canvas.draw_idle()
                    else:
                        if vis:
                            self.annot_dict[key].set_visible(False)
                            self.ppplot.fig.canvas.draw_idle()

    def toggle_interactive(self):
        if self.ui.cb_Interactive.checkState() == 2:
            self.cid_pick = self.ppplot.fig.canvas.mpl_connect("motion_notify_event", self.hover)
        else:
            self.ppplot.fig.canvas.mpl_disconnect(self.cid_pick)

    def process_folder_data(self):
        folder = self.ui.le_SimulationData_Folder.text()
        filename = fr"{self.main_control.projectDir}\PostprocessingData\Data\{self.ui.le_Save_Filename.text()}"
        proc_count = self.ui.sb_No_Of_Processors.value()

        temp_folder = fr"{self.main_control.projectDir}\PostprocessingData\Data\_temp"
        # create temp folder
        if not os.path.exists(fr"{self.main_control.projectDir}\PostprocessingData\Data"):
            os.mkdir(fr"{self.main_control.projectDir}\PostprocessingData\Data")

        # create temp folder
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)

        # if len(list(self._shape_space.keys())) != 0:
        if self.ui.cb_Run_Mode.currentText() == 'Sequential':
            if self.ui.cb_Code.currentText() == "ABCI":
                mon_interval = text_to_list(self.ui.le_Longitudinal_Intervals.text())
                dip_interval = text_to_list(self.ui.le_Transverse_Intervals.text())

                abci_data_ex.multiple_folders_data(self._shape_space, folder, "all", filename, mon_interval,
                                                   dip_interval)

            else:
                request = self.ui.cb_SLANS_Request.currentText()
                mode = self.ui.sb_SLANS_Mode.value()
                bc = self.ui.cb_BC.currentText()
                bc = bc.replace('m', '3')
                bc = bc.replace('e', '2')

                slans_data_ex.multiple_folders_data(folder, mode, bc, request, filename)

        else:
            if self.ui.cb_Code.currentText() == "ABCI":
                mon_interval = text_to_list(self.ui.le_Longitudinal_Intervals.text())
                dip_interval = text_to_list(self.ui.le_Transverse_Intervals.text())
                abci_data_ex.multiple_folders_data_parallel(self._shape_space, folder, proc_count, 'all', filename,
                                                            temp_folder, mon_interval, dip_interval)
            else:
                request = self.ui.cb_SLANS_Request.currentText()
                mode = self.ui.sb_SLANS_Mode.value()

                bc = self.ui.cb_BC.currentText()
                bc = bc.replace('m', '3')
                bc = bc.replace('e', '2')
                slans_data_ex.multiple_folders_data_parallel(self._shape_space, folder, proc_count, mode, bc, request,
                                                             filename, temp_folder)

        # else:
        #     print("Please select at least one key.")

    def open_folder(self, le):
        project_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
        if project_dir != '':
            le.setText(project_dir)

    def open_file(self, le, cb):
        # clear combobox
        self.cb_Shape_Space_Keys.clear()
        self.cb_Shape_Space_Keys.addItem('All')

        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Json Files (*.json)")
        try:
            le.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            # populate checkboxes with key
            for col in dd.keys():
                cb.addItem(fr'{col}')
                print(f"Added col: {col}")

            self._shape_space = dd

        except Exception as e:
            print('Failed to open file:: ', e)

    def clear_plots(self):
        # clear axis
        self.ppplot.ax.cla()

        # clear plot objects
        self.scatter_plot_object = None
        self.scatter_plot_object_inv = None
        self.scatter_pareto_object = None

        # update plot
        self.ppplot.fig.canvas.draw_idle()

