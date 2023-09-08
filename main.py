"""
Created on 12 December 2022
@author: Sosoho-Abasi Udongwo
"""

import ctypes
import logging
import shutil
import sys
from json import JSONDecodeError
from pathlib import Path
from icecream import ic
from utils.misc_functions import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from ui_files.main_window import Ui_MainWindow
from utils.file_reader import FileReader
from frame_controls.multipacting_control import MultipactingControl
from frame_controls.misc import MiscControl
from frame_controls.plot_control import PlotControl
from frame_controls.eigenmode import EigenmodeControl
from frame_controls.postprocess import PostprocessControl
from frame_controls.tune_control import TuneControl
from frame_controls.wakefield_control import WakefieldControl
import qtvscodestyle as qtvsc
from utils.shared_functions import animate_width, f2b_slashes

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
# C:\Users\sosoho\anaconda3\envs\PhD\Scripts\pyrcc5 qss/icons.qrc -o icons_rc.py
# git push -f  https://github.com/Dark-Elektron/CavityDesignHub.git master
# git push https://github.com/Dark-Elektron/CavityDesignHub.git master
# sphinx-autobuild ./source ./_build/html
# sphinx-apidoc -o source .. -f
# make html

fr = FileReader()

myappid = u'll'  # arbitrary string
if os.name == 'nt':
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

DEBUG = True
AN_DURATION = 250


class MainWindow:
    """Main GUI window
    """

    def __init__(self):
        self.tree = None
        self.model = None
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
        self.parentDir = Path(os.getcwd())
        self.projectDir = Path(r'D:/Dropbox/CavityDesignHub/SampleProject')

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

        # # Add stylesheet
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
        self.main_win.setStyleSheet(stylesheet)
        # self.main_win.setStyleSheet("*{font-size: 13px;}")

        # self.ui.pb_rHome.enterEvent = self.tray_animation
        # self.ui.pb_rHome.leaveEvent = self.tray_animation

        # save event filter to variable to be restored later
        # self.default_ef = self.main_win.eventFilter
        # self.main_win.eventFilter = self.eventFilter
        # self.ui.pb_rBack.installEventFilter(self.main_win)
        # self.ui.pb_rHome.installEventFilter(self.main_win)
        # self.ui.pb_rTune.installEventFilter(self.main_win)
        # self.ui.pb_rEigenmode.installEventFilter(self.main_win)
        # self.ui.pb_rWakefield.installEventFilter(self.main_win)
        # self.ui.pb_rMultipacting.installEventFilter(self.main_win)
        # self.ui.pb_rPlot.installEventFilter(self.main_win)
        # self.ui.pb_rPostprocess.installEventFilter(self.main_win)
        # self.ui.pb_rMisc.installEventFilter(self.main_win)

        # ui effects
        try:
            self.ui_effects()
        except AttributeError:
            pass

        # window state before closing
        self.settings = QSettings('MyQtApp', "CavityDesignHub")

        # restore window state before closing
        try:
            self.read_settings()
        except Exception as e:
            print("Exception occurred. ", e)

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
                                               animate_width(self.ui.le_New_Project_Filename, 0, 150, True))

    def initUI(self):
        """
        Initialise class variables and initial GUI state
        """
        # splitter default size
        self.ui.sp_File_Folder_Main.setStretchFactor(1, 4)

        # hold last frame
        self.last_frame = None

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

        self.create_frames_ui()
        self.signals()

        ###########################################

        # Instantiation
        # evolution
        # self.analysis = Analysis(self)

        # house keeping
        # self.hk = HouseKeeping(self)

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
        self.ui.pb_Expand_Collapse_Log.clicked.connect(lambda: animate_width(self.ui.w_Log, 0, 375, True, 'min'))

    def create_frames_ui(self):
        """
        Create the GUI frames
        :return:
        """

        # frame UIs
        self.tune_widget = TuneControl(self)
        self.wakefield_widget = WakefieldControl(self)
        self.eigenmode_widget = EigenmodeControl(self)
        self.postprocess_widget = PostprocessControl(self)
        self.misc_widget = MiscControl(self)
        self.plot_widget = PlotControl(self)
        self.multipacting_widget = MultipactingControl(self)

        self.frames_dict['Home'].append(self.ui.sa_Home)
        self.frames_dict['Tune'].append(self.tune_widget.w_Tune)
        self.frames_dict['Wakefield'].append(self.wakefield_widget.w_Wakefield)
        self.frames_dict['Eigenmode'].append(self.eigenmode_widget.w_Eigenmode)
        self.frames_dict['Multipacting'].append(self.multipacting_widget.w_Multipacting)
        self.frames_dict['Postprocess'].append(self.postprocess_widget.w_Postprocess)
        self.frames_dict['Plot'].append(self.plot_widget.w_Plot)
        self.frames_dict['Misc'].append(self.misc_widget.w_Misc)

    def change_frame(self, key):

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

        # for k, pb in self.tray_buttons_dict.items():
        #     if key == k:
        #         pb[0].setMinimumWidth(75)
        #         pb[0].setMaximumWidth(75)
        #         pb[0].setChecked(True)
        #     else:
        #         pb[0].setMinimumWidth(50)
        #         pb[0].setMaximumWidth(50)
        #         pb[0].setChecked(False)

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

    def create_project(self):
        """
        Create Cavity Design Hub project
        :return:
        """
        project_name = self.ui.le_New_Project_Filename.text()

        if project_name != '':
            project_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
            project_dir = f2b_slashes(project_dir)

            # check if folder already exist
            e = self.checkIfPathExist(project_dir, project_name)

            if e:
                # animate_width(self.ui.le_New_Project_Filename, 0, 150, True)
                self.ui.le_New_Project_Filename.setMinimumWidth(0)
                self.ui.le_New_Project_Filename.setMaximumWidth(0)
                self.ui.l_Project_Name.setText(fr'{project_dir}/{project_name}')

                def make_dirs_from_dict(d, current_dir=project_dir):
                    for key, val in d.items():
                        os.mkdir(f2b_slashes(os.path.join(current_dir, key)))
                        if type(val) == dict:
                            make_dirs_from_dict(val, f2b_slashes(os.path.join(current_dir, key)))

                # create project structure in folders
                project_dir_structure = {
                    f'{project_name}':
                        {
                            'Cavities': None,
                            'OperatingPoints': None,
                            'SimulationData': {
                                'SLANS': None,
                                'NativeEig': None,
                                'ABCI': None,
                                'CavitiesAnalysis': None
                            },
                            'PostprocessingData': {
                                'Plots': None,
                                'Data': None,
                                'CSTData': None
                            },
                            'Reference': None
                        }
                }

                make_dirs_from_dict(project_dir_structure)
                self.projectDir = Path(fr"{project_dir}/{project_name}")

                # only initialize UI after successfully setting folder
                if self.global_state == 0:
                    self.initUI()
                    # add file system tree
                    self.file_system(self.projectDir)
                    self.global_state += 1
        else:
            print('Please enter a valid project name')

    def open_project(self, project_dir=None):
        """

        Parameters
        ----------
        project_dir: str
            Project directory

        Returns
        -------

        """

        if not project_dir:
            project_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
            self.projectDir = Path(project_dir)

            if os.path.exists(self.projectDir / "state_file.json"):
                self.deserialize(self.projectDir / "state_file.json")

            # only initialize UI after successfully setting folder and initialise only once
            self.ui.l_Project_Name.setText(str(self.projectDir))
            if self.global_state == 0:
                print("IOt initialised")
                self.initUI()

                self.global_state += 1

            # add file system tree
            self.file_system(self.projectDir)

        elif project_dir != '':
            try:
                # # check if it's a valid project folder
                # sub_dirs = [a for a in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, a))]
                # compare_dirs = ['Cavities', 'OperatingPoints', 'PostprocessingData', 'SimulationData']
                #
                # if len(set(sub_dirs) & set(sub_dirs)) == len(compare_dirs):
                self.ui.l_Project_Name.setText(project_dir)  # .split('/')[-1]
                self.projectDir = Path(project_dir)

                # only initialize UI after successfully setting folder and initialise only once
                print(self.projectDir)
                self.ui.l_Project_Name.setText(str(self.projectDir))
                if self.global_state == 0:
                    self.initUI()
                    self.global_state += 1

                # add file system tree
                self.file_system(self.projectDir)
            except FileNotFoundError as e:
                print("Invalid project folder not found. ", e)
        else:
            print('Please select a valid project directory')

    def change_theme(self):
        """
        Change GUI theme

        Options available are
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

        Returns
        -------

        """

        stylesheet = qtvsc.load_stylesheet(self.theme_dict[self.ui.cb_Theme.currentText()])
        self.main_win.setStyleSheet(stylesheet)
        self.main_win.setStyleSheet('*{font: 13px "Segoe UI";}')
        self.last_saved_theme = self.ui.cb_Theme.currentText()

        # if self.ui.hs_Theme.value() == 0:
        #     self.stylesheet_filename = 'qss/aqua.qss'
        #     self.loadStylesheet(self.stylesheet_filename)
        # else:
        #     self.stylesheet_filename = 'qss/amoled.qss'
        #     self.loadStylesheet(self.stylesheet_filename)

    def serialize(self):
        """
        Save GUI objects state
        Returns
        -------

        """

        # serialize home
        try:
            # open state file
            with open(Path('ui_state_files/state_file.json'), 'r') as file:
                state_dict = json.load(file)

        except FileNotFoundError:
            print("state_file.json not found, initializing state_dict = {}")
            state_dict = {}

        # update state file○
        # serialize home
        state_dict["Project Directory"] = str(self.projectDir)
        state_dict['Theme'] = self.last_saved_theme

        # serialize tuneUI
        self.tune_widget.serialise(state_dict)

        # serialize eigenmodeUI
        self.eigenmode_widget.serialise(state_dict)

        # serialize wakefieldUI
        self.wakefield_widget.serialise(state_dict)

        # serialize plotUI
        self.plot_widget.serialise(state_dict)

        # dump save state file
        with open(Path('ui_state_files/state_file.json'), 'w', encoding='utf-8') as file:
            file.write(json.dumps(state_dict, indent=4, ensure_ascii=False, separators=(',', ': ')))

        # dump save state file in project folder
        with open(self.projectDir / 'state_file.json', 'w', encoding='utf-8') as file:
            file.write(json.dumps(state_dict, indent=4, ensure_ascii=False, separators=(',', ': ')))

    def deserialize(self, filename):
        """
        Retrieve and update GUI object state from last saved GUI state

        Parameters
        ----------
        filename: str

        Returns
        -------

        """

        # check if state file exists

        filename = Path(filename)
        if os.path.exists(filename):
            # open state file
            try:
                with open(filename, 'r') as f:
                    state_dict = json.load(f)
            except JSONDecodeError:
                print("JSONDecodeError: State file corrupt. Click on save to write new one.")
                state_dict = {}

            # check if state file empty
            if not state_dict:
                # copy default state dict from application directory to project directory
                shutil.copyfile(Path('ui_state_files/_state_file_default.json'), filename)
                with open(filename, 'r') as f:
                    state_dict = json.load(f)

            self.projectDir = state_dict['Project Directory']
            self.last_saved_theme = state_dict['Theme']
            self.ui.cb_Theme.setCurrentText(self.last_saved_theme)

            # open project
            self.open_project(self.projectDir)

            self.tune_widget.deserialise(state_dict)
            self.eigenmode_widget.deserialise(state_dict)
            self.wakefield_widget.deserialise(state_dict)
            self.plot_widget.deserialise(state_dict)

    def load_last_state(self):
        """
        Load GUI last state

        Returns
        -------

        """

        # try:
        self.deserialize(Path('ui_state_files/state_file.json'))
        # except AttributeError as e:
        #     print("Could not deserialize state file: ", e)

    def checkIfPathExist(self, directory, folder):
        path = Path(f"{directory}/{folder}")
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
                # obj.setMaximumWidth(75)
                # obj.setMinimumWidth(75)
                obj.setStyleSheet("background-color: rgb(255, 190, 130); border-radius: 10px; "
                                  "border-style: solid; border-width: 0px; color: white;")
                obj.setIconSize(QSize(75, 50))

                return True

            if event.type() == QEvent.Leave:  # Enter
                # obj.setMaximumSize(50, 50)
                # obj.setMinimumSize(50, 50)
                obj.setStyleSheet("background-color: rgb(255, 170, 127); border-radius: 10px; "
                                  "border-style: solid; border-width: 0px; color: white;")
                obj.setIconSize(QSize(50, 50))
                return True

        return self.default_ef(obj, event)

    def ui_effects(self):
        """
        Control UI effects. Currently turned off.

        Returns
        -------

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
            shadow_effect.setBlurRadius(10)

            push_buttons[0].setStyleSheet("border-radius: 10px; border-style: solid; border-width: 0px; color: white;")
            push_buttons[0].setGraphicsEffect(shadow_effect)

    @staticmethod
    def button_clicked(i):
        return i.text()

    def loadStylesheet(self, filename):
        file = QFile(filename)
        file.open(QFile.ReadOnly | QFile.Text)
        stylesheet = file.readAll()
        self.main_win.setStyleSheet(str(stylesheet, encoding="utf-8"))

    def show(self):
        self.main_win.showMaximized()
        # self.main_win.show()

    def read_settings(self):
        self.main_win.restoreGeometry(self.settings.value("geometry"))
        self.main_win.restoreState(self.settings.value("windowState"))

    def closeEvent(self, event):
        self.settings.setValue("geometry", self.main_win.saveGeometry())
        self.settings.setValue("windowState", self.main_win.saveState(version=0))

    def file_system(self, dir_path):
        if self.model:
            self.model.setRootPath("")
            self.model.setRootPath(str(dir_path))
            self.tree.setRootIndex(self.model.index(str(dir_path)))
        else:
            self.model = QFileSystemModel()
            self.model.setRootPath(str(dir_path))
            self.tree = QTreeView()
            self.tree.setModel(self.model)
            self.tree.setRootIndex(self.model.index(str(dir_path)))
            # self.tree.setColumnWidth(0, 250)
            self.tree.setAlternatingRowColors(True)

            # hide unnecessary details
            for i in range(1, self.tree.model().columnCount()):
                self.tree.header().hideSection(i)

            # set contextMenu
            self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
            self.tree.customContextMenuRequested.connect(self.file_system_context_menu)

            # add functionality to open files
            self.tree.doubleClicked.connect(self.open_file)

            # add to GUI
            self.ui.gl_File_System_View.addWidget(self.tree)

            self.tree.setStyleSheet(
                """
                QTreeView::branch:open:has-children:!has-siblings{image:url(:/icons/icons/PNG/tree_collapse.png);
                    icon-size: 12px 12px;}
                QTreeView::branch:closed:has-children:!has-siblings{image:url(:/icons/icons/PNG/tree_expand.png); 
                   icon-size: 12px 12px;}
                QTreeView::branch:open:has-children{image:url(:/icons/icons/PNG/tree_collapse.png); 
                   icon-size: 12px 12px;}
                QTreeView::branch:closed:has-children{image:url(:/icons/icons/PNG/tree_expand.png);
                    icon-size: 12px 12px;}
                QTreeView::branch:open:{image:url(:/icons/icons/PNG/tree_collapse.png); icon-size: 12px 12px;}
                QTreeView::branch:closed:{image:url(:/icons/icons/PNG/tree_expand.png); icon-size: 12px 12px;}
                """)

    def file_system_context_menu(self):
        menu = QMenu()
        open_ = menu.addAction('Open')
        new = menu.addAction('New')
        delete = menu.addAction('Delete')

        open_.triggered.connect(self.open_file)
        new.triggered.connect(self.new_file)
        delete.triggered.connect(self.delete_file)

        cursor = QCursor()
        menu.exec_(cursor.pos())

    def open_file(self):
        index = self.tree.currentIndex()
        file_path = self.model.filePath(index)

        os.startfile(file_path)

    @staticmethod
    def new_file():
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getSaveFileName(None, "QFileDialog.getSaveFileName()", "",
                                                  "All Files (*);;Text Files (*.txt)", options=options)
        if filename:
            with open(Path(filename), 'w') as file:
                file.write('')

    def delete_file(self):
        index = self.tree.currentIndex()
        file_path = self.model.filePath(index)

        # check if file or folder is selected
        if os.path.isdir(file_path):
            print("Folder")
        else:
            os.remove(file_path)


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
