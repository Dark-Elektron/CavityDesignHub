import math
import os

from PyQt5.QtWidgets import *
from termcolor import colored
import scipy.signal as sps
import numpy as np
import pandas as pd
from ui_files.wg_calculations import Ui_WG_Calc
from modules.data_module.abci_data import ABCIData
from utils.file_reader import FileReader

fr = FileReader()

file_color = 'red'
DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))

c = 299792458 # m / s

class WGCalcControl:
    def __init__(self, parent):
        self.w_WGCalc = QWidget()

        self.wgcalcUI = Ui_WG_Calc()
        self.wgcalcUI.setupUi(self.w_WGCalc)

        # Create parent widget object
        self.wid = parent
        self.misc_control = parent
        self.miscUI = parent.miscUI

        self.signals()
        self.initUI()

    def signals(self):
        # home
        self.wgcalcUI.pb_WGCalc_rMisc.clicked.connect(lambda: self.misc())

        # widget display signals
        # self.wgcalcUI.pb_Dataset_From_Simulation.clicked.connect(lambda: self.toggle_page('Dataset From Simulation'))
        # self.wgcalcUI.pb_Combine_Dataset.clicked.connect(lambda: self.toggle_page('Combine Datasets'))
        # self.wgcalcUI.pb_Filter_Data.clicked.connect(lambda: self.toggle_page('Filter node_editor'))

        # # change code
        # self.wgcalcUI.cb_Run_Mode.currentIndexChanged.connect(lambda: self.show_hide_())

        # frequency calculation signals
        self.wgcalcUI.pb_Calculate_CWG.clicked.connect(lambda: self.calculate_cwg())
        self.wgcalcUI.pb_Calculate_RWG.clicked.connect(lambda: self.calculate_rwg())

    def initUI(self):
        pass

    def calculate_rwg(self):
        try:
            dim_dict = {'mm': 0, 'cm': 1, 'm': 3}
            a = float(self.wgcalcUI.le_A.text())*10**dim_dict[self.wgcalcUI.cb_A_Dim.currentText()]
            b = float(self.wgcalcUI.le_B.text())*10**dim_dict[self.wgcalcUI.cb_B_Dim.currentText()]
            m = float(self.wgcalcUI.le_M.text())
            n = float(self.wgcalcUI.le_N.text())

            f = (c/(2*np.pi))*((m*np.pi/a)**2 + (n*np.pi/b)**2)**0.5
            f = f * 1e-3

            self.wgcalcUI.l_Frequency_RWG.setText(f'{f}')
        except ValueError:
            print("Please enter a valid number.")

    def calculate_cwg(self):
        key = self.wgcalcUI.le_Mode_CWG.text().lower()
        print(key)
        try:
            dim_dict = {'mm': 0, 'cm': 1, 'm': 3}
            dim = dim_dict[self.wgcalcUI.cb_R_Dim.currentText()]
            r = float(self.wgcalcUI.le_R.text())
            r = r*10**dim

            mode_dict = {'te11': 1.841, 'tm01': 2.405, 'te21': 3.054, 'te01': 3.832, 'tm11': 3.832, 'tm21': 5.135,
                         'te12': 5.331, 'tm02': 5.520, 'te22': 6.706, 'te02': 7.016, 'tm12': 7.016, 'tm22': 8.417,
                         'te13': 8.536, 'tm03': 8.654, 'te23': 9.970, 'te03': 10.174, 'tm13': 10.174, 'tm23': 11.620}

            f = mode_dict[key]*c/(2*np.pi*r)
            f = f * 1e-3

            self.wgcalcUI.l_Frequency_CWG.setText(f'{f}')

        except ValueError:
            print("Please enter a valid number.")

    def show_hide_(self, wid1, wid2):
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()

    def toggle_widgets(self, key):
        pass

    def toggle_page(self, key):
        pass

    def misc(self):
        self.misc_control.g_Display.addWidget(self.misc_control.w_Misc)
        self.misc_control.g_Display.removeWidget(self.w_WGCalc)
        self.w_WGCalc.hide()
        self.misc_control.w_Misc.show()
