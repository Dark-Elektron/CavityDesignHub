import math
import os

from PyQt5.QtWidgets import *
from termcolor import colored
import scipy.signal as sps
import numpy as np

import pandas as pd

from frame_controls.misc_widgets.wg_calculations import WGCalcControl
from frame_controls.misc_widgets.mode_nomenclature import ModeNomControl

from ui_files.misc import Ui_Misc
from analysis_modules.data_module.abci_data import ABCIData
from utils.file_reader import FileReader

fr = FileReader()

file_color = 'red'
DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))

class MiscControl:
    def __init__(self, parent):
        self.w_Misc = QWidget()

        self.miscUI = Ui_Misc()
        self.miscUI.setupUi(self.w_Misc)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # g layout
        self.g_Display = self.main_ui.g_Display

        self.frames_dict = {'WGCalc': [self.miscUI.pb_Waveguide_Calculations],
                            'ModeNom': [self.miscUI.pb_Mode_Nomenclature]
                            }
        self.create_frames_ui()
        self.signals()
        self.initUI()

        #

    def signals(self):
        # signal for misc buttons
        for key, button_widget_pair in self.frames_dict.items():
            button_widget_pair[0].clicked.connect(lambda _, b=key: self.change_frame(b))

    def initUI(self):
        pass

    def create_frames_ui(self):
        # frame UIs
        self.wgcalc_widget = WGCalcControl(self)
        self.modenom_widget = ModeNomControl(self)

        self.frames_dict['WGCalc'].append(self.wgcalc_widget.w_WGCalc)
        self.frames_dict['ModeNom'].append(self.modenom_widget.w_ModeNom)

    def change_frame(self, key):
        self.w_Misc.hide()

        # get correct widget
        w = self.frames_dict[key][1]
        self.g_Display.addWidget(w, 0, 1, 1, 1)
        w.show()

    def show_hide_(self, wid1, wid2):
        print('here')
        if wid1.currentText().lower() == 'parallel':
            wid2.show()
        else:
            wid2.hide()
