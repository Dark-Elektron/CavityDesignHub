import math
import os

from PyQt5.QtWidgets import *
from termcolor import colored
import scipy.signal as sps
import numpy as np

import pandas as pd
from ui_files.mode_nomenclature import Ui_Mode_Nom
from modules.data_module.abci_data import ABCIData
from utils.file_reader import FileReader

fr = FileReader()

file_color = 'red'
DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))

class ModeNomControl:
    def __init__(self, parent):
        self.w_ModeNom = QWidget()

        self.modenomUI = Ui_Mode_Nom()
        self.modenomUI.setupUi(self.w_ModeNom)

        # Create parent widget object
        self.wid = parent
        self.misc_control = parent
        self.miscUI = parent.miscUI

        # self.w_dict = {'Dataset': [self.modenomUI.pb_Dataset_From_Simulation, self.modenomUI.w_Dataset_From_Simulation],
        #                'Combine': [self.modenomUI.pb_Combine_Dataset, self.modenomUI.w_Combine_Dataset]}

        self.signals()
        self.initUI()

    def signals(self):
        # home
        self.modenomUI.pb_ModeNom_rMisc.clicked.connect(lambda: self.misc())

    def initUI(self):
        pass

    def misc(self):
        self.misc_control.g_Display.addWidget(self.misc_control.w_Misc)
        self.misc_control.g_Display.removeWidget(self.w_ModeNom)
        self.w_ModeNom.hide()
        self.misc_control.w_Misc.show()
