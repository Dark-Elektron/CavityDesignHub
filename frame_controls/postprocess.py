import ast
import json
import math
import os
import re

import oapackage
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from termcolor import colored
import scipy.signal as sps
import numpy as np
from frame_controls.postprocess_widgets.pandas_widget import PandasModel
import pandas as pd
from ui_files.postprocess import Ui_Postprocess
from modules.data_module.abci_data import ABCIData, ABCIDataExtraction
from modules.data_module.slans_data import *
from utils.file_reader import FileReader
from modules.plot_module.plotter import Plot
from frame_controls.postprocess_widgets.pandas_widget import PandasTV
from frame_controls.postprocess_widgets.ppplot import PPPlot

fr = FileReader()
# abci_data = ABCIData()
abci_data_ex = ABCIDataExtraction()
slans_data_ex = SLANSDataExtraction()


# All_outputs_filtered=(find(All_outputs(:,14)<0.05&All_outputs(:,16)<0.3&All_outputs(:,17)<0.3&All_outputs(:,18)<2.2&All_outputs(:,19)<1&All_outputs(:,20)<3.5&All_outputs(:,29)<2.4&All_outputs(:,30)<6.5 ));%for Ri 150 mm
#
# All_outputs_filtered=(find(All_outputs(:,14)<0.1&All_outputs(:,17)<0.5&All_outputs(:,19)<2&All_outputs(:,20)<1&All_outputs(:,29)<2.4&All_outputs(:,30)<6.7 ));%for Ri 160 mm

file_color = 'red'
DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))


class PostprocessControl:
    def __init__(self, parent):
        # pp -> postprocess
        self.w_Postprocess = QWidget()

        self.ppUI = Ui_Postprocess()
        self.ppUI.setupUi(self.w_Postprocess)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # initialise pandas object
        self.pandas_object = PandasTV(self)
        self.pandasUI = self.pandas_object.pandasUI

        # initialise ppplot object
        self.ppplot = PPPlot(self)
        self.pppUI = self.ppplot.pppUI

        self.w_dict = {'Dataset': [self.ppUI.pb_Dataset_From_Simulation, self.ppUI.w_Dataset_From_Simulation],
                       'Combine': [self.ppUI.pb_Combine_Dataset, self.ppUI.w_Combine_Dataset]}

        # combine dictionaries list
        self.combine_dict_dict = {}

        self.constants()
        self.signals()
        self.initUI()

    def constants(self):
        self.MDI_PLOT_COUNT = 0
        self.MDI_DATAFRAME_COUNT = 0

    def signals(self):

        # widget display signals
        self.ppUI.pb_Dataset_From_Simulation.clicked.connect(lambda: self.toggle_page('Dataset From Simulation'))
        self.ppUI.pb_Combine_Dataset.clicked.connect(lambda: self.toggle_page('Combine Datasets'))
        self.ppUI.pb_Filter_Data.clicked.connect(lambda: self.toggle_page('Filter Data'))

        # load dir
        self.ppUI.pb_Load_Doc1.clicked.connect(lambda: self.load_dir(self.ppUI.le_Dir1.text(),
                                                                              self.ppUI.lw_Doc1,
                                                                              self.ppUI.lw_Selected_Columns_Doc1,
                                                                              self.column_list_1, 0))
        self.ppUI.pb_Load_Doc2.clicked.connect(lambda: self.load_dir(self.ppUI.le_Dir2.text(),
                                                                              self.ppUI.lw_Doc2,
                                                                              self.ppUI.lw_Selected_Columns_Doc2,
                                                                              self.column_list_2, 1))

        # add items to combine
        self.column_list_1 = [] # list to hold added columns
        self.column_list_2 = [] # list to hold added columns
        self.ppUI.pb_Add_Doc1.clicked.connect(lambda: self.add_to_list_widget(self.ppUI.lw_Doc1, self.ppUI.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ppUI.pb_Add_Doc2.clicked.connect(lambda: self.add_to_list_widget(self.ppUI.lw_Doc2, self.ppUI.lw_Selected_Columns_Doc2, self.column_list_2))
        # add all
        self.ppUI.pb_Add_All_Doc1.clicked.connect(lambda: self.add_all(self.ppUI.lw_Doc1, self.ppUI.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ppUI.pb_Add_All_Doc2.clicked.connect(lambda: self.add_all(self.ppUI.lw_Doc2, self.ppUI.lw_Selected_Columns_Doc2, self.column_list_2))

        # combine from folder parallel run
        self.ppUI.pb_Combine_Parallel.clicked.connect(lambda: self.combine_files_parallel_run(self.ppUI.sb_Proc_Count.value(),
                                                               self.ppUI.le_Folder.text(), self.ppUI.le_Filename_Parallel.text()))

        # remove items
        self.ppUI.pb_Remove_Doc1.clicked.connect(lambda: self.remove_from_list_widget(self.ppUI.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ppUI.pb_Remove_Doc2.clicked.connect(lambda: self.remove_from_list_widget(self.ppUI.lw_Selected_Columns_Doc2, self.column_list_2))
        # remove all
        self.ppUI.pb_Remove_All_Doc1.clicked.connect(lambda: self.remove_all(self.ppUI.lw_Doc1, self.ppUI.lw_Selected_Columns_Doc1, self.column_list_1))
        self.ppUI.pb_Remove_All_Doc2.clicked.connect(lambda: self.remove_all(self.ppUI.lw_Doc2, self.ppUI.lw_Selected_Columns_Doc2, self.column_list_2))

        # combine files
        self.ppUI.pb_Combine.clicked.connect(lambda: self.combine_files(self.ppUI.le_Filename.text()))

        # load excel file
        self.ppUI.pb_Select_File.clicked.connect(lambda: self.load_file())

        # objective function/pareto
        self.ppUI.pb_Plot_Objective.clicked.connect(lambda: self.plot_objective_function())
        self.ppUI.cb_Pareto.stateChanged.connect(lambda: self.plot_pareto())

        # mdi signals
        self.ppUI.pb_MDI_Tiled.clicked.connect(lambda: self.ppUI.mdiArea.tileSubWindows())
        self.ppUI.pb_MDI_Cascade.clicked.connect(lambda: self.ppUI.mdiArea.cascadeSubWindows())

        # collapse analysis/right menu
        self.ppUI.pb_Collapse_Analysis_Menu.clicked.connect(lambda: self.main_control.animate_width(self.ppUI.w_Analysis_Menu, 0, 375, True))

        # analysis menu signals
        self.ppUI.pb_Objective_Function.clicked.connect(lambda: self.ppUI.sw_Analysis_Menu.setCurrentIndex(0))
        self.ppUI.pb_Filter.clicked.connect(lambda: self.ppUI.sw_Analysis_Menu.setCurrentIndex(1))

        # add/remove filter
        self.ppUI.pb_Add_DFFilter.clicked.connect(lambda: self.add_filter())
        self.ppUI.pb_Remove_Filter.clicked.connect(lambda: self.remove_filter(self.ppUI.pb_Remove_Filter))

        # enable/disable filter
        self.ppUI.cb_Filter_Toggle_0.clicked.connect(lambda: self.apply_filter())

        # enable interactive
        self.ppUI.cb_Interactive.clicked.connect(lambda: self.toggle_interactive())

        # load shape space
        self.ppUI.pb_Select_Shape_Space.clicked.connect(lambda: self.open_file(self.ppUI.le_Shape_Space, self.cb_Shape_Space_Keys))

        # select folder
        self.ppUI.pb_Select_Folder.clicked.connect(lambda: self.open_folder(self.ppUI.le_SimulationData_Folder))

        # extract data
        self.ppUI.pb_Extract.clicked.connect(lambda: self.process_folder_data())

        # ABCI/SLANS inputs
        self.ppUI.cb_Code.currentTextChanged.connect(lambda: self.code_change_control())

        # sequential vs parallel
        self.ppUI.cb_Run_Mode.currentTextChanged.connect(lambda: self.run_mode_control())

    def initUI(self):
        # hide code
        self.ppUI.w_SLANS.hide()

        # diable combine button
        self.ppUI.pb_Combine.setEnabled(False)

        # initial table dataframe
        self.df = None
        self.filtered_df = None

        # dataframe filters dict
        self.df_filter_dict = {0: [self.ppUI.w_Filter_0, self.ppUI.le_LB_0, self.ppUI.le_RB_0, self.ppUI.cb_DFFilter_0, self.ppUI.cb_Filter_Toggle_0, self.ppUI.pb_Remove_Filter]}
        self.filter_count = 1

        # plot variables
        self.x = []
        self.y = []
        self.x_pareto = []
        self.y_pareto = []

        # add checkable combobox
        self.cb_Shape_Space_Keys = CheckableComboBox()
        self.cb_Shape_Space_Keys.setMinimumWidth(75)
        self.ppUI.gl_Selected_Shapes.addWidget(self.cb_Shape_Space_Keys)

        # shape space initialization
        self._shape_space = {}

        self.code_change_control()

        # sequential vs parallel
        self.run_mode_control()

        # interactive plot
        # self.cid_click = self.ppplot.fig.canvas.mpl_connect('button_press_event', self.onclick)
        # self.cid_pick = self.ppplot.fig.canvas.mpl_connect('pick_event', self.onpick)
        # self.ppplot.fig.canvas.mpl_connect("motion_notify_event", self.hover)

    def run_mode_control(self):
        if self.ppUI.cb_Run_Mode.currentText() == "Sequential":
            self.ppUI.w_No_Of_Processors.setMinimumWidth(0)
            self.ppUI.w_No_Of_Processors.setMaximumWidth(0)
        else:
            self.ppUI.w_No_Of_Processors.setMinimumWidth(0)
            self.ppUI.w_No_Of_Processors.setMaximumWidth(0)
            self.main_control.animate_width(self.ppUI.w_No_Of_Processors, 0, 300, True)

    def code_change_control(self):
        if self.ppUI.cb_Code.currentText() == 'ABCI':
            self.ppUI.w_Intervals.show()
            self.ppUI.w_SLANS.hide()

            # set default post processing folder to project directory
            self.ppUI.le_SimulationData_Folder.setText(fr"{self.main_control.projectDir}\SimulationData\ABCI")
        else:
            self.ppUI.w_Intervals.hide()
            self.ppUI.w_SLANS.show()

            # set default post processing folder to project directory
            self.ppUI.le_SimulationData_Folder.setText(fr"{self.main_control.projectDir}\SimulationData\SLANS")

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
                print(problem_keys_sorted)
                delete_problem_keys(problem_keys_sorted)
                print(len(Zmax_mon_list), len(Zmax_dip_list), len(A))
                # save excel
                if save_excel:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
                            f'Z_long[max({mon_mask[0]}<f<{mon_mask[1]})]': Zmax_mon_list,
                            f'Z_trans[max({dip_mask[0]}<f<{dip_mask[1]})]': Zmax_dip_list}

                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', sheet_name='Sheet_1')

            if request == 'all':
                all()
                print(len(k_loss_array_longitudinal), len(k_loss_array_transverse), len(df_M0D1_array), len(A))
                problem_keys_sorted = sorted(problem_keys,
                                             reverse=True)  # this is so as to start deleting from the largest key. If this is not done, the index is messed up
                print(problem_keys_sorted)
                delete_problem_keys(problem_keys_sorted)
                print(len(k_loss_array_longitudinal), len(k_loss_array_transverse), len(df_M0D1_array), len(A))

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

    def combine_files(self, save_excel):
        if len(list(self.combine_dict_dict.keys())) == 2:
            f1 = self.combine_dict_dict[0]
            f2 = self.combine_dict_dict[1]

            f3 = f1[self.column_list_1].merge(f2[self.column_list_2], on=self.ppUI.cb_On.currentText(), how=self.ppUI.cb_How.currentText())
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
        lw1 = self.ppUI.lw_Selected_Columns_Doc1
        lw2 = self.ppUI.lw_Selected_Columns_Doc2

        # get item text list from list widgets
        itemsTextList1 = [str(lw1.item(i).text()) for i in range(lw1.count())]
        itemsTextList2 = [str(lw2.item(i).text()) for i in range(lw2.count())]

        # check for intersection in itemTextList
        intersection_list = list(set(itemsTextList1) & set(itemsTextList2))

        # check if intersection_list empty, populate cb_On and disable combine button if yes
        if intersection_list:
            # clear cb_On
            self.ppUI.cb_On.clear()
            # populate cb_On
            for key in intersection_list:
                self.ppUI.cb_On.addItem(key)

            # enable combine button
            self.ppUI.pb_Combine.setEnabled(True)
        else:
            self.ppUI.pb_Combine.setEnabled(False)

    def show_hide_(self, wid1, wid2):
        print('here')
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
        if key == 'Dataset From Simulation':
            self.ppUI.stackedWidget.setCurrentIndex(0)

        if key == 'Combine Datasets':
            self.ppUI.stackedWidget.setCurrentIndex(1)

        if key == 'Filter Data':
            self.ppUI.stackedWidget.setCurrentIndex(2)

    def load_file(self):
        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Excel Files (*.xlsx)")
        self.ppUI.le_Pandas_Filename.setText(filename)
        try:
            self.df = pd.read_excel(filename)
            self.pandas_model = PandasModel(self.df)
            self.pandasUI.tv_Pandas.setModel(self.pandas_model)

            if self.MDI_DATAFRAME_COUNT == 0:
                # show in mdi area
                sub = QMdiSubWindow()
                sub.setWidget(self.pandas_object.w_Pandas)
                sub.setWindowTitle('Dataframe')
                self.ppUI.mdiArea.addSubWindow(sub)
                sub.show()
                self.MDI_DATAFRAME_COUNT += 1

            # populate first filter
            for val in self.df_filter_dict.values():
                self.populate_filter(val[3], self.df)

        except:
            print('Please select file.')

    def process_objective_function_input(self, exp, norm):
        if self.df is not None:
            operands = []
            for i, e in enumerate(exp):
                if self.filtered_df is not None:
                    operands.append(self.filtered_df.iloc[:, e]/norm[i])
                else:
                    print(self.df.iloc[:, e])
                    operands.append(self.df.iloc[:, e]/norm[i])

            var = sum(operands)
            # print(var.tolist())

            return var.tolist()

    def plot_objective_function(self):
        # get operands
        exp_x = ast.literal_eval(self.ppUI.le_Objective_X.text())
        exp_y = ast.literal_eval(self.ppUI.le_Objective_Y.text())
        # get norms
        x_norm = ast.literal_eval(self.ppUI.le_X_Norm.text())
        y_norm = ast.literal_eval(self.ppUI.le_Y_Norm.text())

        print(exp_x, x_norm, exp_y, y_norm)
        print('='*50)

        # check if norm count equal operand count
        if len(exp_x) == len(x_norm) and len(exp_y) == len(y_norm):
            self.ppplot.ax.cla()
            self.x = self.process_objective_function_input(exp_x, x_norm)
            self.y = self.process_objective_function_input(exp_y, y_norm)

            self.scatter_plot_object = self.ppplot.ax.scatter(self.x, self.y, facecolors='none', edgecolors='b', picker=True)

            # add annotation
            self.annot = self.ppplot.ax.annotate("", xy=(0,0), xytext=(-20,-10),textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w"),
                                arrowprops=dict(arrowstyle="->"))
            self.annot.set_visible(True)

            self.ppplot.plt.draw()

            if self.MDI_PLOT_COUNT == 0:
                # show in mdi area
                sub = QMdiSubWindow()
                sub.setWidget(self.ppplot.w_PPPlot)
                sub.setWindowTitle('Plot')
                self.ppUI.mdiArea.addSubWindow(sub)
                sub.show()
                self.MDI_PLOT_COUNT += 1
        else:
            print("Operand and Norm count are not equal. Please check.")

    def plot_pareto(self):
        if self.ppUI.cb_Pareto.checkState() == 2:
            self.x_pareto, self.y_pareto = self.pareto_front(self.x, self.y)
            self.scatter_pareto_object = self.ppplot.ax.scatter(self.x_pareto, self.y_pareto, c='r', picker=True)

            self.ppplot.plt.draw()
        else:
            pass

    def pareto_front(self, x, y, reverse='bottom'):
        if reverse == 'top':
            def reverse_list(l):
                return l
        else:
            def reverse_list(l):
                return [-x for x in l]

        datapoints = np.array([reverse_list(x), reverse_list(y)])
        #     print(datapoints)

        pareto = oapackage.ParetoDoubleLong()

        for ii in range(0, datapoints.shape[1]):
            w = oapackage.doubleVector((datapoints[0, ii], datapoints[1, ii]))
            pareto.addvalue(w, ii)

        pareto.show(verbose=1)

        lst = pareto.allindices()  # the indices of the Pareto optimal designs

        optimal_datapoints = datapoints[:, lst]

        return reverse_list(optimal_datapoints[0, :]), reverse_list(optimal_datapoints[1, :])

    def populate_filter(self, cb, df):
        for col in df.columns:
            cb.addItem(col)

    def add_filter(self):
        # create filter widgets
        w = QWidget()
        l = QGridLayout(w)

        lb = QLineEdit()
        l.addWidget(lb, 0, 0, 1, 1)
        rb = QLineEdit()
        l.addWidget(rb, 0, 2, 1, 1)
        cb = QComboBox()
        cb.setMinimumWidth(75)
        cb.setMaximumWidth(75)
        l.addWidget(cb, 0, 1, 1, 1)
        chk = QCheckBox()
        l.addWidget(chk, 0, 3, 1, 1)
        # add signal
        chk.stateChanged.connect(lambda: self.apply_filter())

        pb = QPushButton()
        pb.setMaximumWidth(25)
        pb.setMinimumWidth(25)
        pb.setText('X')
        pb.setStyleSheet('background-color: rgb(255, 89, 67);')
        l.addWidget(pb, 0, 4, 1, 1)
        # add signal
        pb.clicked.connect(lambda: self.remove_filter(pb))

        self.ppUI.gl_Filter.addWidget(w)

        # add to UI
        self.ppUI.gl_Filter.addWidget(w)

        # populate combobox
        if self.df is not None:
            self.populate_filter(cb, self.df)

        # add to filter dictionary
        self.df_filter_dict[self.filter_count] = [w, lb, rb, cb, chk, pb]
        self.filter_count += 1

    def remove_filter(self, pb_clicked):
        k = None
        for key, wids in self.df_filter_dict.items():
            print(key, pb_clicked, self.df_filter_dict)
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
        # check if a dataframe already exists
        if self.df is not None:
            # copy main data frame
            self.filtered_df = self.df.copy(deep=True)

            # apply filter from to dataframe
            # check if filter dict is not empty
            if self.df_filter_dict:
                for key, val in self.df_filter_dict.items():
                    if val[4].checkState() == 2:
                        lb = float(val[1].text())
                        rb = float(val[2].text())
                        key = val[3].currentText()

                        self.filtered_df = self.filtered_df[(self.filtered_df[key] >= lb) & (self.filtered_df[key] <= rb)]
                        self.filtered_df = self.filtered_df.reset_index(drop=True)

            self.pandas_model = PandasModel(self.filtered_df)
            self.pandasUI.tv_Pandas.setModel(self.pandas_model)

        # update plot if any
        self.plot_objective_function()

    def add_handles(self):
        pass

    def onpick(self, event):
        print('here')
        ind = event.ind
        print('onpick scatter:', ind, np.take(self.x, ind), np.take(self.y, ind))
        print(self.df.loc[ind, 'A':'alpha'])

    def update_annot(self, ind):
        pos = self.scatter_plot_object.get_offsets()[ind]
        self.annot.xy = pos
        try:
            text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'key']}"
            # text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'alpha']}"
            # text = f"ti: {ind}, {self.filtered_df.loc[ind, 'key':'B']}"
        except:
            text = f"ti: {ind}, {self.df.loc[ind, 'key':'key']}"
            # text = f"ti: {ind}, {self.df.loc[ind, 'key':'Req']}"

        self.annot.set_text(f'{text}')
        self.annot.get_bbox_patch().set_facecolor((167/255, 222/255, 255/255))
        # self.annot.get_bbox_patch().set_alpha(1)

    def hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ppplot.ax:
            cont, ind = self.scatter_plot_object.contains(event)
            if cont:
                self.update_annot(ind["ind"][0]) # ind returns an array of close points. ind['ind'][0] returns just the first point
                self.annot.set_visible(True)
                self.ppplot.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.ppplot.fig.canvas.draw_idle()

    def toggle_interactive(self):
        if self.ppUI.cb_Interactive.checkState() == 2:
            self.cid_pick = self.ppplot.fig.canvas.mpl_connect("motion_notify_event", self.hover)
        else:
            self.ppplot.fig.canvas.mpl_disconnect(self.cid_pick)

    def process_folder_data(self):
        folder = self.ppUI.le_SimulationData_Folder.text()
        filename = fr"{self.main_control.projectDir}\PostprocessingData\Data\{self.ppUI.le_Save_Filename.text()}"
        proc_count = self.ppUI.sb_No_Of_Processors.value()
        temp_folder = fr"{self.main_control.projectDir}\PostprocessingData\Data\_temp"

        if len(list(self._shape_space.keys())) != 0:
            if self.ppUI.cb_Run_Mode.currentText() == 'Sequential':
                if self.ppUI.cb_Code.currentText() == "ABCI":
                    mon_interval = self.text_to_list(self.ppUI.le_Longitudinal_Intervals.text())
                    dip_interval = self.text_to_list(self.ppUI.le_Transverse_Intervals.text())

                    abci_data_ex.multiple_folders_data(self._shape_space, folder, "all", filename, mon_interval, dip_interval)

                else:
                    request = self.ppUI.cb_SLANS_Request.currentText()
                    mode = self.ppUI.sb_SLANS_Mode.value()
                    bc = self.ppUI.cb_BC.currentText()

                    slans_data_ex.multiple_folders_data(self._shape_space, folder, mode, bc, request, filename)

            else:
                if self.ppUI.cb_Code.currentText() == "ABCI":
                    mon_interval = self.text_to_list(self.ppUI.le_Longitudinal_Intervals.text())
                    dip_interval = self.text_to_list(self.ppUI.le_Transverse_Intervals.text())
                    abci_data_ex.multiple_folders_data_parallel(self._shape_space, folder, proc_count, 'all', filename, temp_folder, mon_interval, dip_interval)
                else:
                    request = self.ppUI.cb_SLANS_Request.currentText()
                    mode = self.ppUI.sb_SLANS_Mode.value()
                    bc = self.ppUI.cb_BC.currentText()

                    slans_data_ex.multiple_folders_data_parallel(self._shape_space, folder, proc_count, mode, bc, request, filename, temp_folder)

        else:
            print("Please select at least one key.")

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

            self._shape_space = dd

        except Exception as e:
            print('Failed to open file:: ', e)

    def text_to_list(self, l):
        if l == '':
            return None
        else:
            l = ast.literal_eval(l)
            if isinstance(l, int) or isinstance(l, float):
                return [l, 2e10]
            else:
                return list(l)

class CheckableComboBox(QComboBox):
    # Subclass Delegate to increase item height
    class Delegate(QStyledItemDelegate):
        def sizeHint(self, option, index):
            size = super().sizeHint(option, index)
            size.setHeight(20)
            return size

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.lineEdit().setReadOnly(True)
        # Make the lineedit the same color as QPushButton
        palette = qApp.palette()
        palette.setBrush(QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # Use custom delegate
        self.setItemDelegate(CheckableComboBox.Delegate())

        # Update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # Hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)
        self.closeOnLineEditClick = False

        # Prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

    def resizeEvent(self, event):
        # Recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, object, event):

        if object == self.lineEdit():
            if event.type() == QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if object == self.view().viewport():
            if event.type() == QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                item = self.model().item(index.row())

                if item.checkState() == Qt.Checked:
                    item.setCheckState(Qt.Unchecked)

                    if item == self.model().item(0):
                        # deselect all items if item check is all
                        for i in range(1, self.model().rowCount()):
                            item = self.model().item(i)
                            item.setCheckState(Qt.Unchecked)
                else:
                    item.setCheckState(Qt.Checked)

                    if item == self.model().item(0):
                        # deselect all items if item check is all
                        for i in range(1, self.model().rowCount()):
                            item = self.model().item(i)
                            item.setCheckState(Qt.Checked)

                return True
        return False

    def showPopup(self):
        super().showPopup()
        # When the popup is displayed, a click on the lineedit should close it
        self.closeOnLineEditClick = True

    def hidePopup(self):
        super().hidePopup()
        # Used to prevent immediate reopening when clicking on the lineEdit
        self.startTimer(100)
        # Refresh the display text when closing
        self.updateText()

    def timerEvent(self, event):
        # After timeout, kill timer, and reenable click on line edit
        self.killTimer(event.timerId())
        self.closeOnLineEditClick = False

    def updateText(self):
        texts = []
        for i in range(1, self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                texts.append(self.model().item(i).text())
        text = ", ".join(texts)
        self.lineEdit().setText(text)

        # # Compute elided text (with "...")
        # metrics = QFontMetrics(self.lineEdit().font())
        # elidedText = metrics.elidedText(text, Qt.ElideRight, self.lineEdit().width())
        # self.lineEdit().setText(elidedText)

    def addItem(self, text, data=None):
        item = QStandardItem()
        item.setText(text)
        if data is None:
            item.setData(text)
        else:
            item.setData(data)
        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
        item.setData(Qt.Unchecked, Qt.CheckStateRole)
        self.model().appendRow(item)

    def addItems(self, texts, datalist=None):
        for i, text in enumerate(texts):
            try:
                data = datalist[i]
            except (TypeError, IndexError):
                data = None
            self.addItem(text, data)

    def currentData(self):
        # Return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                res.append(self.model().item(i).data())
        return res