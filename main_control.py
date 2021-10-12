import ctypes
import logging
import os
import shutil
import sys

from utils.house_keeping import HouseKeeping
from utils.misc_functions import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from modules.analysis_module.analysis_codes import Analysis
from ui_files.main_window import Ui_MainWindow
from ui_files.run_tune import Ui_w_Tune
from ui_files.wakefield import Ui_Wakefield
from utils.file_reader import FileReader
from frame_controls.misc import MiscControl
from frame_controls.plot_control import PlotControl
from frame_controls.eigenmode import EigenmodeControl
from frame_controls.postprocess import PostprocessControl
from frame_controls.tune_control import TuneControl
from frame_controls.wakefield_control import WakefieldControl


# pyuic5 -x ui_files/main_window.ui -o ui_files/main_window.py
# pyuic5 -x ui_files/eigenmode.ui -o ui_files/eigenmode.py
# pyuic5 -x ui_files/postprocess.ui -o ui_files/postprocess.py
# pyuic5 -x ui_files/run_tune.ui -o ui_files/run_tune.py
# pyuic5 -x ui_files/plot.ui -o ui_files/plot.py
# pyuic5 -x ui_files/wakefield.ui -o ui_files/wakefield.py
# pyuic5 -x ui_files/misc.ui -o ui_files/misc.py
# pyuic5 -x ui_files/wg_calculations.ui -o ui_files/wg_calculations.py
# pyuic5 -x ui_files/mode_nomenclature.ui -o ui_files/mode_nomenclature.py
# pyuic5 -x ui_files/pandas_table.ui -o ui_files/pandas_table.py
# pyuic5 -x ui_files/pp_plot.ui -o ui_files/pp_plot.py
# pyrcc5 qss/icons.qrc -o icons_rc.py

fr = FileReader()

myappid = u'll' # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

AN_DURATION = 200

class MainWindow:
    def __init__(self):
        self.main_win = QMainWindow()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_win)

        self.global_state = 0

        # Add stylesheet
        self.stylesheet_filename = 'qss/aqua.qss'
        self.loadStylesheet(self.stylesheet_filename)

        # self.initUI()
        self.parentDir = os.getcwd()
        self.projectDir = r'D:\Dropbox\CEMCodesHub\SampleProject'

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
        self.ui.pb_rPlot.installEventFilter(self.main_win)
        self.ui.pb_rPostprocess.installEventFilter(self.main_win)
        self.ui.pb_rMisc.installEventFilter(self.main_win)

        # ui effects
        self.ui_effects()

    def new_open_folder(self):
        # disable all buttons
        # signals to create folder
        self.ui.le_New_Project_Filename.returnPressed.connect(lambda: self.create_project())
        self.ui.pb_Open_Project.clicked.connect(lambda: self.open_project())

        # create new/open project
        # set size to zero
        self.ui.le_New_Project_Filename.setMaximumWidth(0)
        self.ui.pb_New_Project.clicked.connect(lambda: self.animate_width(self.ui.le_New_Project_Filename, 0, 150, True))

    def initUI(self):
        # hidden widgets

        # hold last frame
        self.last_frame = None

        self.frames_dict = {'Home': [self.ui.pb_rHome],
                            'Tune': [self.ui.pb_Tune],
                            'Wakefield': [self.ui.pb_Wakefield_Analysis],
                            'Eigenmode': [self.ui.pb_Eigenmode_Analysis],
                            'Postprocess': [self.ui.pb_Post_Process],
                            'Plot': [self.ui.pb_Plot],
                            'Misc': [self.ui.pb_Misc]
                            }

        self.tray_buttons_dict = {'Home': [self.ui.pb_rHome],
                                  'Tune': [self.ui.pb_rTune],
                                  'Wakefield': [self.ui.pb_rWakefield],
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
        self.hk = HouseKeeping(self)

        #############################################
        # set initial session log widget state to collapsed
        self.ui.w_Log.setMinimumWidth(0)
        self.ui.w_Log.setMaximumWidth(0)

    def signals(self):
        # signal for home buttons
        for key, button_widget_pair in self.frames_dict.items():
            button_widget_pair[0].clicked.connect(lambda _, b=key: self.change_frame(b))

        # signal for butom panel
        self.ui.pb_rTune.clicked.connect(lambda: self.change_frame('Tune'))
        self.ui.pb_rEigenmode.clicked.connect(lambda: self.change_frame('Eigenmode'))
        self.ui.pb_rWakefield.clicked.connect(lambda: self.change_frame('Wakefield'))
        self.ui.pb_rPostprocess.clicked.connect(lambda: self.change_frame('Postprocess'))
        self.ui.pb_rPlot.clicked.connect(lambda: self.change_frame('Plot'))
        self.ui.pb_rMisc.clicked.connect(lambda: self.change_frame('Misc'))
        self.ui.pb_rHome.clicked.connect(lambda: self.change_frame('Home'))
        self.ui.pb_rBack.clicked.connect(lambda: self.return_to_last_frame())

        # theme signal
        self.ui.hs_Theme.valueChanged.connect(lambda: self.change_theme())

        # save state
        self.ui.pb_Save_State.clicked.connect(lambda: self.serialize())

        # # load state file
        # self.ui.pb_Load_State_File.clicked.connect(lambda: self.deserialize())

        # expand/collapse session log
        self.ui.pb_Expand_Collapse_Log.clicked.connect(lambda: self.animate_width(self.ui.w_Log, 0, 375, True, 'min'))

    def create_frames_ui(self):
        # frame UIs
        self.tune_widget = TuneControl(self)
        self.wakefield_widget = WakefieldControl(self)
        self.eigenmode_widget = EigenmodeControl(self)
        self.postprocess_widget = PostprocessControl(self)
        self.plot_widget = PlotControl(self)
        self.misc_widget = MiscControl(self)

        self.frames_dict['Home'].append(self.ui.sa_Home)
        self.frames_dict['Tune'].append(self.tune_widget.w_Tune)
        self.frames_dict['Wakefield'].append(self.wakefield_widget.w_Wakefield)
        self.frames_dict['Eigenmode'].append(self.eigenmode_widget.w_Eigenmode)
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
        if self.last_frame is not None:
            index = self.ui.g_Display.count() - 1
            while (index >= 0):
                widget = self.ui.g_Display.itemAt(index).widget()
                self.ui.g_Display.removeWidget(widget)
                widget.hide()
                index -= 1

            # get correct widget
            self.ui.g_Display.addWidget(self.last_frame, 0, 1, 1, 1)
            self.last_frame.show()

    def animate_width(self, widget, min_width, standard, enable, option="max"):
        if enable:
            #### GET WIDTH
            width = widget.width()
            #### SET MAX WIDTH
            if width > 0:
                widthCollapsed = min_width
                widget.setMinimumWidth(0)
            else:
                widthCollapsed = standard
                # widget.setMinimumWidth(standard)

            #### ANIMATION
            if option=='max':
                self.animation = QPropertyAnimation(widget, b"maximumWidth")
            else:
                self.animation = QPropertyAnimation(widget, b"minimumWidth")

            self.animation.setDuration(AN_DURATION)
            self.animation.setStartValue(width)
            # print_(widthCollapsed)
            self.animation.setEndValue(widthCollapsed)
            self.animation.start()

    def animate_height(self, widget, min_height, standard, enable, option="max"):
        if enable:
            #### GET WIDTH
            height = widget.height()

            #### SET MAX WIDTH
            if height > 0:
                heightCollapsed = min_height
                widget.setMinimumHeight(0)
            else:
                heightCollapsed = standard
                # self.ui.w_Shape_Parameters.setMinimumSize(0, 250)

            #### ANIMATION
            if option=='max':
                self.animation = QPropertyAnimation(widget, b"maximumHeight")
            else:
                self.animation = QPropertyAnimation(widget, b"minimumHeight")
            self.animation.setDuration(AN_DURATION)
            self.animation.setStartValue(height)
            self.animation.setEndValue(heightCollapsed)
            self.animation.start()

    def create_project(self):
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
        if not project_dir:
            project_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
            self.projectDir = self.f2b_slashes(project_dir)
            print('open project', self.projectDir)

            if os.path.exists(f'{project_dir}\state_file.json'):
                self.deserialize(f'{project_dir}\state_file.json')

            # only initialize UI after successfully setting folder and initialise only once
            self.ui.l_Project_Name.setText(self.projectDir)
            if self.global_state == 0:
                self.initUI()
                self.global_state += 1

        elif project_dir != '':
            # check if it's a valid project folder
            sub_dirs = [a for a in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, a))]
            compare_dirs = ['Cavities', 'PostprocessingData', 'SimulationData']
            if len(set(sub_dirs) & set(sub_dirs)) == len(compare_dirs):
                self.ui.l_Project_Name.setText(project_dir) #.split('/')[-1]
                self.projectDir = self.f2b_slashes(project_dir)

                # only initialize UI after successfully setting folder and initialise only once
                self.ui.l_Project_Name.setText(self.projectDir)
                if self.global_state == 0:
                    self.initUI()
                    self.global_state += 1


            else:
                print('Please select a valid project directory')
        else:
            print('Please select a valid project directory')

    def change_theme(self):
        if self.ui.hs_Theme.value() == 0:
            self.stylesheet_filename = 'qss/aqua.qss'
            self.loadStylesheet(self.stylesheet_filename)
        else:
            self.stylesheet_filename = 'qss/amoled.qss'
            self.loadStylesheet(self.stylesheet_filename)

    def serialize(self):
        # serialize home
        try:
            # open state file
            state_dict = fr.json_reader('ui_state_files/state_file.json')
        except:
            state_dict = {}

        # update state file
        # serialize home
        state_dict["Project Directory"] = self.projectDir

        # serialize plotUI
        self.plot_widget.serialize(state_dict)

        # serialize tunetUI
        self.tune_widget.serialize(state_dict)

        # dump save state file
        with open('ui_state_files/state_file.json', 'w', encoding='utf-8') as file:
            file.write(json.dumps(state_dict, indent=4, ensure_ascii=False, separators=(',', ': ')))

        # dump save state file in project folder
        with open(f'{self.projectDir}/state_file.json', 'w', encoding='utf-8') as file:
            file.write(json.dumps(state_dict, indent=4, ensure_ascii=False, separators=(',', ': ')))

    def deserialize(self, file):
        # check if state file exists
        if os.path.exists(file):
            try:
                # open state file
                with open(file, 'r') as f:
                    state_dict = json.load(f)
    
                self.projectDir = state_dict['Project Directory']
    
                # open project
                self.open_project(self.projectDir)

                # deserialise plotUI
                try:
                    self.plot_widget.deserialize(state_dict)
                except:
                    print("Could not deserialise plot!")

                # deserialise tuneUI
                try:
                    self.tune_widget.deserialize(state_dict)
                except:
                    print("Could not deserialise tune!")

            except:
                print("Corrupt state file!")

    def load_last_state(self):
        msg = QMessageBox()
        msg.setWindowTitle("Folder Exist")
        msg.setText("Do you wish to load the last saved project?")
        msg.setIcon(QMessageBox.Question)
        msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        msg.setDefaultButton(QMessageBox.Yes)

        msg.buttonClicked.connect(self.button_clicked)

        x = msg.exec_()

        if x == msg.Yes:
            self.deserialize('ui_state_files/state_file.json')

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

    def eventFilter(self, object, event):
        if not object.isChecked():
            if event.type() == QEvent.Enter: # Enter
                object.setMaximumWidth(75)
                object.setMinimumWidth(75)

                return True

            if event.type() == QEvent.Leave: # Enter
                object.setMaximumSize(50,50)
                object.setMinimumSize(50, 50)
                return True

        return self.default_ef(object, event)

    def ui_effects(self):
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
        self.main_win.show()


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
        logging.ERROR:   ("[%(levelname)-8s] %(message)s", QColor("red")),
        logging.DEBUG:   ("[%(levelname)-8s] [%(filename)s:%(lineno)d] %(message)s", "green"),
        logging.INFO:    ("[%(levelname)-8s] %(message)s", "#0000FF"),
        logging.WARNING: ('%(asctime)s - %(name)s - %(levelname)s - %(message)s', QColor(100, 100, 0))
    }

    def format( self, record ):
        last_fmt = self._style._fmt
        opt = CustomFormatter.FORMATS.get(record.levelno)
        if opt:
            fmt, color = opt
            self._style._fmt = "<font color=\"{}\">{}</font>".format(QColor(color).name(),fmt)
        res = logging.Formatter.format( self, record )
        self._style._fmt = last_fmt
        return res


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_win = MainWindow()
    main_win.show()

    sys.exit(app.exec_())

