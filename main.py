"""
Created on 12 December 2022
@author: Sosoho-Abasi Udongwo
"""

import ctypes
import logging
import os
import shutil
import sys
# import pyautogui
from PyQt5 import QtCore
from utils.misc_functions import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from ui_files.main_window import Ui_MainWindow
from ui_files.run_tune import Ui_w_Tune
from ui_files.wakefield import Ui_Wakefield
from utils.file_reader import FileReader
from frame_controls.multipacting_control import MultipactingControl
from frame_controls.misc import MiscControl
from frame_controls.plot_control import PlotControl
from frame_controls.eigenmode import EigenmodeControl
from frame_controls.postprocess import PostprocessControl
from frame_controls.tune_control import TuneControl, OptimizationControl
from frame_controls.wakefield_control import WakefieldControl
from node_editor.node_editor_widget import NodeEditorWidget
import qtvscodestyle as qtvsc

# sys.path.append(r"D:\Dropbox\CavityDesignHub\test_plugin") # to search inside this directory for imports


# pyuic5 -x ui_files/main_window.ui -o ui_files/main_window.py
# pyuic5 -x ui_files/eigenmode.ui -o ui_files/eigenmode.py
# pyuic5 -x ui_files/postprocess.ui -o ui_files/postprocess.py
# pyuic5 -x ui_files/multipacting.ui -o ui_files/multipacting.py
# pyuic5 -x ui_files/run_tune.ui -o ui_files/run_tune.py
# pyuic5 -x ui_files/plot.ui -o ui_files/plot.py
# pyuic5 -x ui_files/wakefield.ui -o ui_files/wakefield.py
# pyuic5 -x ui_files/geometry_input.ui -o ui_files/geometry_input.py
# pyuic5 -x ui_files/misc.ui -o ui_files/misc.py
# pyuic5 -x ui_files/wg_calculations.ui -o ui_files/wg_calculations.py
# pyuic5 -x ui_files/mode_nomenclature.ui -o ui_files/mode_nomenclature.py
# pyuic5 -x ui_files/pandas_table.ui -o ui_files/pandas_table.py
# pyuic5 -x ui_files/pp_plot.ui -o ui_files/pp_plot.py
# pyuic5 -x ui_files/plot_properties.ui -o ui_files/pprops.py
# pyuic5 -x ui_files/edit_annotated_text.ui -o ui_files/edit_annotated_text.py
# pyuic5 -x ui_files/edit_line2D.ui -o ui_files/edit_line2D.py
# pyuic5 -x ui_files/plot_widget_header.ui -o ui_files/plot_widget_header.py
# pyuic5 -x ui_files/abci_plot.ui -o ui_files/abci_plot.py
# pyuic5 -x ui_files/plottypeselector.ui -o ui_files/plottypeselector.py
# pyrcc5 qss/icons.qrc -o icons_rc.py
# git push -f  https://github.com/Dark-Elektron/CavityDesignHub.git master
# sphinx-autobuild ./source ./_build/html
# sphinx-apidoc -o source .. -f


fr = FileReader()

myappid = u'll'  # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

DEBUG = True
AN_DURATION = 250


class MainWindow:
    """
    This does good stuff.

    Here are the details about the good stuff it does.

    Parameters
    ----------


    Returns
    -------
    y : int
        Some other thing
    """
    def __init__(self):
        self.animation = None
        self.wakefield_widget = None
        self.eigenmode_widget = None
        self.misc_widget = None
        self.plot_widget = None
        self.postprocess_widget = None
        self.multipacting_widget = None
        self.tune_widget = None
        self.frames_dict = None
        self.last_frame = None
        self.tray_buttons_dict = None
        self.main_win = QMainWindow()
        self.main_win.closeEvent = self.closeEvent

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_win)

        # set window height
        # width, height = pyautogui.size()
        # self.main_win.setFixedHeight(height)

        self.global_state = 0
        self.last_saved_theme = 'Light VS'

        # self.initUI()
        self.parentDir = os.getcwd()
        self.projectDir = r'D:\Dropbox\CavityDesignHub\SampleProject'

        # add node editor
        # new = NodeEditorWidget()
        # new.addNodes()
        # self.ui.gl_Node_Editor.addWidget(new)

        # initialize logging
        self.logTextBox = QPlainTextEditLogger()
        self.log = logging
        # log options debug(), info(), warning(), error()

        # add to layout
        self.ui.gl_Log.addWidget(self.logTextBox.widget)
        self.logTextBox.setFormatter(self.log.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.log.getLogger().addHandler(self.logTextBox)
        # change formatter to include color
        self.logTextBox.setFormatter(CustomFormatter())

        # You can control the logging level
        self.log.getLogger().setLevel(self.log.DEBUG)
        self.log.info('Session started')

        #
        self.new_open_folder()

        # check if last state should be loaded
        self.load_last_state()

        # Add stylesheet
        # self.stylesheet_filename = 'qss/aqua.qss'
        # self.loadStylesheet(self.stylesheet_filename)
        self.theme_dict = {'Light VS': qtvsc.Theme.LIGHT_VS,
                           'Quiet Light': qtvsc.Theme.QUIET_LIGHT,
                           'Solarized Light': qtvsc.Theme.SOLARIZED_LIGHT,
                           'Abyss': qtvsc.Theme.ABYSS,
                           'Dark VS': qtvsc.Theme.DARK_VS,
                           'Kimbie Dark': qtvsc.Theme.KIMBIE_DARK,
                           'Monokai': qtvsc.Theme.MONOKAI,
                           'Monokai Dimmed': qtvsc.Theme.MONOKAI_DIMMED,
                           'Red': qtvsc.Theme.RED,
                           'Solarized Dark': qtvsc.Theme.SOLARIZED_DARK,
                           'Tomorrow Night Blue': qtvsc.Theme.TOMORROW_NIGHT_BLUE,
                           'Dark High Contrast': qtvsc.Theme.DARK_HIGH_CONTRAST
                           }
        # Theme name             Symbol
        # ___________            ______
        #
        # Light (Visual Studio): LIGHT_VS
        # Quiet Light          : QUIET_LIGHT
        # Solarized Light      : SOLARIZED_LIGHT
        # Abyss                : ABYSS
        # Dark (Visual Studio) : DARK_VS
        # Kimbie Dark          : KIMBIE_DARK
        # Monokai              : MONOKAI
        # Monokai Dimmed       : MONOKAI_DIMMED
        # Red                  : RED
        # Solarized Dark       : SOLARIZED_DARK
        # Tomorrow Night Blue  : TOMORROW_NIGHT_BLUE
        # Dark High Contrast   : DARK_HIGH_CONTRAST
        # stylesheet = qtvsc.load_stylesheet(qtvsc.Theme.DARK_VS)
        stylesheet = qtvsc.load_stylesheet(self.theme_dict[self.last_saved_theme])
        QApplication.instance().setStyleSheet(stylesheet)

        # self.ui.pb_rHome.enterEvent = self.tray_animation
        # self.ui.pb_rHome.leaveEvent = self.tray_animation

        # save event filter to variable to be restored later
        self.default_ef = self.main_win.eventFilter
        self.main_win.eventFilter = self.eventFilter
        self.ui.pb_rBack.installEventFilter(self.main_win)
        self.ui.pb_rHome.installEventFilter(self.main_win)
        self.ui.pb_rTune.installEventFilter(self.main_win)
        self.ui.pb_rEigenmode.installEventFilter(self.main_win)
        self.ui.pb_rWakefield.installEventFilter(self.main_win)
        self.ui.pb_rMultipacting.installEventFilter(self.main_win)
        self.ui.pb_rPlot.installEventFilter(self.main_win)
        self.ui.pb_rPostprocess.installEventFilter(self.main_win)
        self.ui.pb_rMisc.installEventFilter(self.main_win)

        # ui effects
        self.ui_effects()

        # window state before closing
        self.settings = QSettings('MyQtApp', "CavityDesignHub")
        print(self.settings.fileName())
        # restore window state before closing
        try:
            self.readSettings()
        except Exception as e:
            print("Exception occured", e)

    def new_open_folder(self):
        """
        Method to create new project or open old projects.

        """
        # disable all buttons
        # signals to create folder
        self.ui.le_New_Project_Filename.returnPressed.connect(lambda: self.create_project())
        self.ui.pb_Open_Project.clicked.connect(lambda: self.open_project())

        # create new/open project
        # set size to zero
        self.ui.le_New_Project_Filename.setMaximumWidth(0)
        self.ui.pb_New_Project.clicked.connect(lambda:
                                               self.animate_width(self.ui.le_New_Project_Filename, 0, 150, True))

    def initUI(self):
        """
        Initialise class variables and initial GUI state
        """
        # hidden widgets

        # hold last frame
        self.last_frame = None
        if DEBUG: print("Check 2e_i: main_control.py")

        self.frames_dict = {'Home': [self.ui.pb_rHome],
                            'Tune': [self.ui.pb_Tune],
                            'Wakefield': [self.ui.pb_Wakefield_Analysis],
                            'Multipacting': [self.ui.pb_Multipacting_Analysis],
                            'Eigenmode': [self.ui.pb_Eigenmode_Analysis],
                            'Postprocess': [self.ui.pb_Post_Process],
                            'Plot': [self.ui.pb_Plot],
                            'Misc': [self.ui.pb_Misc]
                            }

        self.tray_buttons_dict = {'Home': [self.ui.pb_rHome],
                                  'Tune': [self.ui.pb_rTune],
                                  'Wakefield': [self.ui.pb_rWakefield],
                                  'Multipacting': [self.ui.pb_rMultipacting],
                                  'Eigenmode': [self.ui.pb_rEigenmode],
                                  'Postprocess': [self.ui.pb_rPostprocess],
                                  'Plot': [self.ui.pb_rPlot],
                                  'Misc': [self.ui.pb_rMisc]
                                  }

        if DEBUG: print("Check 2e_ii: main_control.py")
        self.create_frames_ui()
        if DEBUG: print("Check 2e_iii: main_control.py")
        self.signals()
        if DEBUG: print("Check 2e_iii: main_control.py")

        ###########################################

        # Instantiation
        # evolution
        # self.analysis = Analysis(self)

        # house keeping
        # self.hk = HouseKeeping(self)
        if DEBUG: print("Check 2e_iii: main_control.py")

        #############################################
        # set initial session log widget state to collapsed
        self.ui.w_Log.setMinimumWidth(0)
        self.ui.w_Log.setMaximumWidth(0)

    def signals(self):
        """
        All the signals for the PyQt objects can be found here.
        :return:
        """
        # signal for home buttons
        for key, button_widget_pair in self.frames_dict.items():
            button_widget_pair[0].clicked.connect(lambda _, b=key: self.change_frame(b))

        # signal for button panel
        self.ui.pb_rTune.clicked.connect(lambda: self.change_frame('Tune'))
        self.ui.pb_rEigenmode.clicked.connect(lambda: self.change_frame('Eigenmode'))
        self.ui.pb_rWakefield.clicked.connect(lambda: self.change_frame('Wakefield'))
        self.ui.pb_rMultipacting.clicked.connect(lambda: self.change_frame('Multipacting'))
        self.ui.pb_rPostprocess.clicked.connect(lambda: self.change_frame('Postprocess'))
        self.ui.pb_rPlot.clicked.connect(lambda: self.change_frame('Plot'))
        self.ui.pb_rMisc.clicked.connect(lambda: self.change_frame('Misc'))
        self.ui.pb_rHome.clicked.connect(lambda: self.change_frame('Home'))
        self.ui.pb_rBack.clicked.connect(lambda: self.return_to_last_frame())

        # theme signal
        self.ui.pb_Apply_Theme.clicked.connect(lambda: self.change_theme())

        # save state
        self.ui.pb_Save_State.clicked.connect(lambda: self.serialize())

        # # load state file
        # self.ui.pb_Load_State_File.clicked.connect(lambda: self.deserialize())

        # expand/collapse session log
        self.ui.pb_Expand_Collapse_Log.clicked.connect(lambda: self.animate_width(self.ui.w_Log, 0, 375, True, 'min'))

    def create_frames_ui(self):
        """
        Create the GUI frames
        :return:
        """

        if DEBUG: print("Check 2e_i_1: main_control.py")
        # frame UIs
        self.tune_widget = TuneControl(self)
        if DEBUG: print("Check 2e_i_1: main_control.py")
        self.wakefield_widget = WakefieldControl(self)
        if DEBUG: print("Check 2e_i_1: main_control.py")
        self.eigenmode_widget = EigenmodeControl(self)
        if DEBUG: print("Check 2e_i_1: main_control.py")
        self.postprocess_widget = PostprocessControl(self)
        if DEBUG: print("Check 2e_i_1: main_control.py")
        self.misc_widget = MiscControl(self)
        if DEBUG: print("Check 2e_i_1: main_control.py")
        self.plot_widget = PlotControl(self)
        if DEBUG: print("Check 2e_i_1: main_control.py")
        self.multipacting_widget = MultipactingControl(self)
        if DEBUG: print("Check 2e_i_2: main_control.py")

        self.frames_dict['Home'].append(self.ui.sa_Home)
        self.frames_dict['Tune'].append(self.tune_widget.w_Tune)
        self.frames_dict['Wakefield'].append(self.wakefield_widget.w_Wakefield)
        self.frames_dict['Eigenmode'].append(self.eigenmode_widget.w_Eigenmode)
        self.frames_dict['Multipacting'].append(self.multipacting_widget.w_Multipacting)
        self.frames_dict['Postprocess'].append(self.postprocess_widget.w_Postprocess)
        self.frames_dict['Plot'].append(self.plot_widget.w_Plot)
        self.frames_dict['Misc'].append(self.misc_widget.w_Misc)
        if DEBUG: print("Check 2e_i_3: main_control.py")

    def change_frame(self, key):
        """
        Switch between GUI frames
        :param key:
        :return:
        """
        # remove existing widgets
        index = self.ui.g_Display.count() - 1
        # save last frame
        self.last_frame = self.ui.g_Display.itemAt(index).widget()

        while index >= 0:
            widget = self.ui.g_Display.itemAt(index).widget()
            self.ui.g_Display.removeWidget(widget)
            widget.hide()
            index -= 1

        # get correct widget
        w = self.frames_dict[key][1]
        self.ui.g_Display.addWidget(w, 0, 1, 1, 1)
        w.show()

        for k, pb in self.tray_buttons_dict.items():
            if key == k:
                pb[0].setChecked(True)
                pb[0].setMinimumWidth(75)
                pb[0].setMaximumWidth(75)
            else:
                pb[0].setChecked(False)
                pb[0].setMinimumWidth(50)
                pb[0].setMaximumWidth(50)

    def return_to_last_frame(self):
        """
        Is called when return button is clicked on the GUI
        :return:
        """
        if self.last_frame is not None:
            index = self.ui.g_Display.count() - 1
            while index >= 0:
                widget = self.ui.g_Display.itemAt(index).widget()
                self.ui.g_Display.removeWidget(widget)
                widget.hide()
                index -= 1

            # get correct widget
            self.ui.g_Display.addWidget(self.last_frame, 0, 1, 1, 1)
            self.last_frame.show()

    def animate_width(self, widget, min_width, standard, enable, option="max"):
        """
        Animate width for GUI transition effect
        :param widget: PyQt Widget
        :param min_width: int or float
        :param standard: int or float
        :param enable:
        :param option:
        :return:
        """
        if enable:
            # GET WIDTH
            width = widget.width()
            # SET MAX WIDTH
            if width > 0:
                widthCollapsed = min_width
                widget.setMinimumWidth(0)
            else:
                widthCollapsed = standard
                # widget.setMinimumWidth(standard)

            # ANIMATION
            if option == 'max':
                self.animation = QPropertyAnimation(widget, b"maximumWidth")
            else:
                self.animation = QPropertyAnimation(widget, b"minimumWidth")

            self.animation.setDuration(AN_DURATION)
            self.animation.setStartValue(width)
            self.animation.setEndValue(widthCollapsed)
            self.animation.setEasingCurve(QtCore.QEasingCurve.InOutQuart)
            self.animation.start()

    def animate_height(self, widget, min_height, standard, enable, option="max"):
        """
        Animate height for GUI transition effect
        :param widget: PyQt Widget
        :param min_height: int or float
        :param standard: int or float
        :param enable:
        :param option:
        :return:
        """
        if enable:
            # GET WIDTH
            height = widget.height()

            # SET MAX WIDTH
            if height > 0:
                heightCollapsed = min_height
                widget.setMinimumHeight(0)
            else:
                heightCollapsed = standard
                # self.ui.w_Shape_Parameters.setMinimumSize(0, 250)

            # ANIMATION
            if option == 'max':
                self.animation = QPropertyAnimation(widget, b"maximumHeight")
            else:
                self.animation = QPropertyAnimation(widget, b"minimumHeight")
            self.animation.setDuration(AN_DURATION)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.setEasingCurve(QtCore.QEasingCurve.InOutQuart)
            self.animation.start()

    def create_project(self):
        """
        Create project
        :return:
        """
        project_name = self.ui.le_New_Project_Filename.text()

        if project_name != '':
            project_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
            project_dir = self.f2b_slashes(project_dir)

            # check if folder already exist
            e = self.checkIfPathExist(project_dir, project_name)

            if e:
                # self.animate_width(self.ui.le_New_Project_Filename, 0, 150, True)
                self.ui.le_New_Project_Filename.setMinimumWidth(0)
                self.ui.le_New_Project_Filename.setMaximumWidth(0)
                self.ui.l_Project_Name.setText(fr'{project_dir}\{project_name}')

                def make_dirs_from_dict(d, current_dir=project_dir):
                    for key, val in d.items():
                        os.mkdir(os.path.join(current_dir, key))
                        if type(val) == dict:
                            make_dirs_from_dict(val, os.path.join(current_dir, key))

                # create project structure in folders
                project_dir_structure = {f'{project_name}':
                                             {'Cavities': None,
                                              'SimulationData': {
                                                  'SLANS': None,
                                                  'ABCI': None
                                              },
                                              'PostprocessingData': {
                                                  'Plots': None,
                                                  'Data': None,
                                                  'CSTData': None
                                              }
                                              }
                                         }

                make_dirs_from_dict(project_dir_structure)
                self.projectDir = self.f2b_slashes(fr"{project_dir}\{project_name}")

                # only initialize UI after successfully setting folder
                if self.global_state == 0:
                    self.initUI()
                    self.global_state += 1
        else:
            print('Please enter a valid project name')

    def open_project(self, project_dir=None):
        """
        Open project
        :param project_dir:
        :return:
        """
        if DEBUG: print("Check 2a: main_control.py")
        if not project_dir:
            project_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
            self.projectDir = self.f2b_slashes(project_dir)
            print('open project', self.projectDir)

            if os.path.exists(fr'{project_dir}\state_file.json'):
                self.deserialize(fr'{project_dir}\state_file.json')

            # only initialize UI after successfully setting folder and initialise only once
            self.ui.l_Project_Name.setText(self.projectDir)
            if self.global_state == 0:
                self.initUI()
                self.global_state += 1

        elif project_dir != '':
            if DEBUG: print("Check 2b: main_control.py")
            # check if it's a valid project folder
            sub_dirs = [a for a in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, a))]
            compare_dirs = ['Cavities', 'PostprocessingData', 'SimulationData']
            if DEBUG: print("Check 2c: main_control.py")
            if len(set(sub_dirs) & set(sub_dirs)) == len(compare_dirs):
                self.ui.l_Project_Name.setText(project_dir)  # .split('/')[-1]
                self.projectDir = self.f2b_slashes(project_dir)
                if DEBUG: print("Check 2d: main_control.py")

                # only initialize UI after successfully setting folder and initialise only once
                self.ui.l_Project_Name.setText(self.projectDir)
                if self.global_state == 0:
                    if DEBUG: print("Check 2e: main_control.py")
                    self.initUI()
                    self.global_state += 1
                if DEBUG: print("Check 2f: main_control.py")

            else:
                print('Please select a valid project directory')
        else:
            print('Please select a valid project directory')

    def change_theme(self):
        """
        Change GUI theme
        :return:
        :desc:
            Light (Visual Studio): LIGHT_VS
            Quiet Light          : QUIET_LIGHT
            Solarized Light      : SOLARIZED_LIGHT
            Abyss                : ABYSS
            Dark (Visual Studio) : DARK_VS
            Kimbie Dark          : KIMBIE_DARK
            Monokai              : MONOKAI
            Monokai Dimmed       : MONOKAI_DIMMED
            Red                  : RED
            Solarized Dark       : SOLARIZED_DARK
            Tomorrow Night Blue  : TOMORROW_NIGHT_BLUE
            Dark High Contrast   : DARK_HIGH_CONTRAST
            stylesheet = qtvsc.load_stylesheet(qtvsc.Theme.DARK_VS)
        """

        stylesheet = qtvsc.load_stylesheet(self.theme_dict[self.ui.cb_Theme.currentText()])
        QApplication.instance().setStyleSheet(stylesheet)
        self.last_saved_theme = self.ui.cb_Theme.currentText()

        # change plot widget colors
        self.plot_widget.plt.change_background('#586e75')

        # if self.ui.hs_Theme.value() == 0:
        #     self.stylesheet_filename = 'qss/aqua.qss'
        #     self.loadStylesheet(self.stylesheet_filename)
        # else:
        #     self.stylesheet_filename = 'qss/amoled.qss'
        #     self.loadStylesheet(self.stylesheet_filename)

    def serialize(self):
        """
        Save GUI object states
        :return:
        """
        # serialize home
        try:
            # open state file
            state_dict = fr.json_reader('ui_state_files/state_file.json')
        except FileNotFoundError:
            print("state_file.json not found, initializing state_dict = {}")
            state_dict = {}

        # update state file
        # serialize home
        state_dict["Project Directory"] = self.projectDir
        state_dict['Theme'] = self.last_saved_theme

        # serialize tuneUI
        self.tune_widget.serialize(state_dict)

        # serialize eigenmodeUI
        self.eigenmode_widget.serialize(state_dict)

        # serialize wakefieldUI
        self.wakefield_widget.serialize(state_dict)

        # serialize plotUI
        self.plot_widget.serialize(state_dict)

        # dump save state file
        with open('ui_state_files/state_file.json', 'w', encoding='utf-8') as file:
            file.write(json.dumps(state_dict, indent=4, ensure_ascii=False, separators=(',', ': ')))

        # dump save state file in project folder
        with open(f'{self.projectDir}/state_file.json', 'w', encoding='utf-8') as file:
            file.write(json.dumps(state_dict, indent=4, ensure_ascii=False, separators=(',', ': ')))

    def deserialize(self, file):
        """
        Retrieve and update GUI object state from last saved GUI state
        :param file:
        :return:
        """
        if DEBUG: print("Check 1: main_control.py")
        # check if state file exists
        if os.path.exists(file):
            # open state file
            with open(file, 'r') as f:
                state_dict = json.load(f)

            # try:
            self.projectDir = state_dict['Project Directory']
            self.last_saved_theme = state_dict['Theme']
            self.ui.cb_Theme.setCurrentText(self.last_saved_theme)
            # except:
            #     print("Could not deserialize theme and directory!")

            # open project
            if DEBUG: print("Check 2: main_control.py")
            self.open_project(self.projectDir)

            if DEBUG: print("Check 2end: main_control.py")

            # deserialise tuneUI
            # try:
            if DEBUG: print("Check 3: main_control.py")
            self.tune_widget.deserialize(state_dict)
            # except:
            #     print("Could not deserialize tuneUI!")

            # deserialise eigenmodeUI
            # try:
            if DEBUG: print("Check 4: main_control.py")
            self.eigenmode_widget.deserialize(state_dict)
            # except:
            #     print("Could not deserialise eigenmodeUI!")

            # deserialise wakefieldUI
            # try:
            if DEBUG: print("Check 5: main_control.py")
            self.wakefield_widget.deserialize(state_dict)
            # except:
            #     print("Could not deserialise wakefieldUI!")

            # deserialise plotUI
            # try:
            if DEBUG: print("Check 6: main_control.py")
            self.plot_widget.deserialize(state_dict)
            # except:
            #     print("Could not deserialise plotUI!")

    def load_last_state(self):
        # msg = QMessageBox()
        # msg.setWindowTitle("Folder Exist")
        # msg.setText("Do you wish to load the last saved project?")
        # msg.setIcon(QMessageBox.Question)
        # msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        # msg.setDefaultButton(QMessageBox.Yes)
        #
        # msg.buttonClicked.connect(self.button_clicked)
        #
        # x = msg.exec_()
        #
        # if x == msg.Yes:
        #     self.deserialize('ui_state_files/state_file.json')

        try:
            self.deserialize('ui_state_files/state_file.json')
        except:
            print("Could not deserialize last state")

    def checkIfPathExist(self, directory, folder):
        path = f"{directory}/{folder}"
        if os.path.exists(path):
            msg = QMessageBox()
            msg.setWindowTitle("Folder Exist")
            msg.setText("Project already exists. Do you want to overwrite?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            msg.setDefaultButton(QMessageBox.Yes)
            msg.buttonClicked.connect(self.button_clicked)
            x = msg.exec_()

            if x == msg.Yes:
                try:
                    shutil.rmtree(path)
                    return True
                except:
                    return False
            else:
                return False

        else:
            return True

    def eventFilter(self, obj, event):
        if not obj.isChecked():
            if event.type() == QEvent.Enter:  # Enter
                obj.setMaximumWidth(75)
                obj.setMinimumWidth(75)

                return True

            if event.type() == QEvent.Leave:  # Enter
                obj.setMaximumSize(50, 50)
                obj.setMinimumSize(50, 50)
                return True

        return self.default_ef(obj, event)

    def ui_effects(self):
        """
        Control UI effects. Currently turned off.
        :return:
        """
        for push_buttons in self.frames_dict.values():
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setOffset(5)
            shadow_effect.setColor(QColor(0, 0, 0, 77))
            shadow_effect.setBlurRadius(5)

            push_buttons[0].setGraphicsEffect(shadow_effect)

        for push_buttons in self.tray_buttons_dict.values():
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setOffset(2.5)
            shadow_effect.setColor(QColor(0, 0, 0, 77))
            shadow_effect.setBlurRadius(5)

            push_buttons[0].setGraphicsEffect(shadow_effect)

    @staticmethod
    def button_clicked(i):
        return i.text()

    @staticmethod
    def f2b_slashes(path):
        # replaces forward slashes with backward slashes for windows OS
        path = path.replace(r"/", "\\")
        return path

    @staticmethod
    def loadStylesheet(filename):
        file = QFile(filename)
        file.open(QFile.ReadOnly | QFile.Text)
        stylesheet = file.readAll()
        QApplication.instance().setStyleSheet(str(stylesheet, encoding="utf-8"))

    def show(self):
        self.main_win.showMaximized()
        # self.main_win.show()

    def readSettings(self):
        self.main_win.restoreGeometry(self.settings.value("geometry"))
        self.main_win.restoreState(self.settings.value("windowState"))

    def closeEvent(self, event):
        self.settings.setValue("geometry", self.main_win.saveGeometry())
        self.settings.setValue("windowState", self.main_win.saveState(version=0))


class QPlainTextEditLogger(logging.Handler):
    def __init__(self, parent=None):
        super().__init__()
        self.widget = QPlainTextEdit(parent)
        self.widget.setReadOnly(True)

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendHtml(msg)
        # move scrollbar
        scrollbar = self.widget.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())


class CustomFormatter(logging.Formatter):
    FORMATS = {
        logging.ERROR: ("[%(levelname)-8s] %(message)s", QColor("red")),
        logging.DEBUG: ("[%(levelname)-8s] [%(filename)s:%(lineno)d] %(message)s", "green"),
        logging.INFO: ("[%(levelname)-8s] %(message)s", "#0000FF"),
        logging.WARNING: ('%(asctime)s - %(name)s - %(levelname)s - %(message)s', QColor(100, 100, 0))
    }

    def format(self, record):
        last_fmt = self._style._fmt
        opt = CustomFormatter.FORMATS.get(record.levelno)
        if opt:
            fmt, color = opt
            self._style._fmt = "<font color=\"{}\">{}</font>".format(QColor(color).name(), fmt)
        res = logging.Formatter.format(self, record)
        self._style._fmt = last_fmt
        return res


if __name__ == '__main__':
    # # Handle high resolution displays:
    # if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    #     QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    # if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    #     QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)

    main_win = MainWindow()

    main_win.show()
    sys.exit(app.exec_())
