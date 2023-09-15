import json
import os
import psutil
from PyQt5 import QtCore
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_template import FigureCanvas
from matplotlib.figure import Figure
from termcolor import colored
from analysis_modules.data_module.slans_data import SLANSDataExtraction

slans_data_extraction = SLANSDataExtraction()

file_color = 'red'
DEBUG = False


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class QCheckableComboBox(QComboBox):

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
        # self.lineEdit().setReadOnly(False)
        self.lineEdit().setReadOnly(True)
        # Make the lineedit the same color as QPushButton
        palette = qApp.palette()
        palette.setBrush(QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # Use custom delegate
        self.setItemDelegate(QCheckableComboBox.Delegate())

        # Update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # Hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)
        self.closeOnLineEditClick = False

        # Prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

        # combobox text list
        self.texts = []

    def resizeEvent(self, event):
        # Recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, obj, event):
        if obj == self.lineEdit():
            if event.type() == QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if obj == self.view().viewport():
            if event.type() == QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                item = self.model().item(index.row())

                if item.checkState() == Qt.Checked:
                    item.setCheckState(Qt.Unchecked)

                    if item == self.model().item(0):
                        # deselect all items if item check a maximum of 25 elements
                        for i in range(1, 25):
                            try:
                                item = self.model().item(i)
                                item.setCheckState(Qt.Unchecked)
                            except AttributeError:
                                break
                else:
                    item.setCheckState(Qt.Checked)

                    if item == self.model().item(0):
                        # deselect all items if item check a maximum of 25 elements
                        for i in range(1, 25):
                            try:
                                item = self.model().item(i)
                                item.setCheckState(Qt.Checked)
                            except AttributeError:
                                break

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
        self.texts.append(text)
        if data is None:
            item.setData(text)
        else:
            item.setData(data)
        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
        item.setData(Qt.Unchecked, Qt.CheckStateRole)
        self.model().appendRow(item)
        self.updateText()

    def addItems(self, texts, datalist=None):
        for i, text in enumerate(texts):
            try:
                data = datalist[i]
            except (TypeError, IndexError):
                data = None
            self.addItem(text, data)
        self.updateText()

    def currentData(self):
        # Return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                res.append(self.model().item(i).data())
        return res

    def clear(self) -> None:
        super().clear()
        # reset text list
        self.texts = []


class QColorComboBox(QComboBox):
    """ A drop down menu for selecting colors """

    # signal emitted if a color has been selected
    selectedColor = QtCore.pyqtSignal(QColor)

    def __init__(self, parent=None, enableUserDefColors=True):
        """ if the user shall not be able to define colors on its own, then set enableUserDefColors=False """
        # init QComboBox
        super(QColorComboBox, self).__init__(parent)

        # enable the line edit to display the currently selected color
        self.setEditable(True)
        # read only so that there is no blinking cursor or sth editable
        self.lineEdit().setReadOnly(True)

        # text that shall be displayed for the option to pop up the QColorDialog for user defined colors
        self._userDefEntryText = 'Custom'
        # add the option for user defined colors
        if enableUserDefColors:
            self.addItem(self._userDefEntryText)

        self._currentColor = None

        self.activated.connect(self._color_selected)

    # ------------------------------------------------------------------------
    def addColors(self, colors):
        """ Adds colors to the QComboBox """
        for a_color in colors:
            # if input is not a QColor, try to make it one
            if not (isinstance(a_color, QColor)):
                a_color = QColor(a_color)
            # avoid dublicates
            if self.findData(a_color) == -1:
                # add the new color and set the background color of that item
                self.addItem('', userData=a_color)
                self.setItemData(self.count() - 1, QColor(a_color), Qt.BackgroundRole)

    # ------------------------------------------------------------------------
    def addColor(self, color):
        """ Adds the color to the QComboBox """
        self.addColors([color])

    # ------------------------------------------------------------------------
    def setColor(self, color):
        """ Adds the color to the QComboBox and selects it"""
        self.addColor(color)
        self._color_selected(self.findData(color), False)

    # ------------------------------------------------------------------------
    def getCurrentColor(self):
        """ Returns the currently selected QColor
            Returns None if non has been selected yet
        """
        return self._currentColor

    # ------------------------------------------------------------------------
    def _color_selected(self, index, emitSignal=True):
        """ Processes the selection of the QComboBox """
        # if a color is selected, emit the selectedColor signal
        if self.itemText(index) == '':
            self._currentColor = self.itemData(index)
            if emitSignal:
                self.selectedColor.emit(self._currentColor)

        # if the user wants to define a custom color
        elif self.itemText(index) == self._userDefEntryText:
            # get the user defined color
            new_color = QColorDialog.getColor(self._currentColor if self._currentColor else Qt.white)
            if new_color.isValid():
                # add the color to the QComboBox and emit the signal
                self.addColor(new_color)
                self._currentColor = new_color
                if emitSignal:
                    self.selectedColor.emit(self._currentColor)

        # make sure that current color is displayed
        if self._currentColor:
            self.setCurrentIndex(self.findData(self._currentColor))
            self.lineEdit().setStyleSheet("background-color: " + self._currentColor.name())


class ProgressMonitor(QThread):
    sig = QtCore.pyqtSignal(int)

    def __init__(self, frame, projectDir):
        super(QThread, self).__init__()
        self.frame = frame
        self.proc_ids = frame.processes_id
        self.progress_bar = frame.progress_bar
        self.projectDir = projectDir

    def run(self):
        self.progress_monitor()

    def progress_monitor(self):
        while self.frame.show_progress_bar:
            progress = len(self.frame.progress_list)
            self.sig.emit(progress)


class EndRoutine(QThread):
    def __init__(self, frame, projectDir):
        super(QThread, self).__init__()
        self.frame = frame
        self.proc_ids = frame.processes_id
        self.filename = frame.filename
        self.projectDir = projectDir

    def run(self):
        self.end_routine()

    def end_routine(self):
        proc_count = len(self.proc_ids)
        for pid in self.proc_ids:
            try:
                p = psutil.Process(pid)
                while p.is_running():
                    pass

            except psutil.NoSuchProcess:
                pass

        # combine dictionaries
        try:
            combined_dict = self.combine_dict(proc_count, self.filename, self.projectDir)
            self.delete_process_dict(proc_count, self.projectDir)

            print(self.frame.ui.cb_Toggle_Postprocess.checkState())
            if self.frame.ui.cb_Toggle_Postprocess.isChecked():
                print("It's here. if it doesn't end it mean s that thre is a problem")
                # postprocess results
                slans_data_dir = f"{self.projectDir}/SimulationData/SLANS"
                mode = self.frame.ui.sb_Mode.value()
                bc = self.frame.ui.cb_BC.currentText()
                request = self.frame.ui.cb_Request.currentText()
                filename = self.frame.ui.le_Postprocess_Filename.text()
                # filename = self.frame.proof_filename(filename)
                save_excel = f"{self.projectDir}/PostprocessingData/Data/{filename}"
                proc_count = self.frame.ui.sb_No_Of_Processors_Tune.value()
                temp_folder = f"{self.projectDir}/PostprocessingData/Data/_temp"

                # ic(slans_data_dir, proc_count, mode, bc, request, save_excel, temp_folder)
                # print(combined_dict)
                slans_data_extraction.multiple_folders_data_parallel(combined_dict, slans_data_dir, proc_count, mode,
                                                                     bc, request, save_excel, temp_folder)
                print("It's ended")

        except IOError as e:
            self.frame.log.error(f"Some error occurred -> {e}")

        self.frame.cancel()

    @staticmethod
    def combine_dict(proc_count, filename, projectDir):
        # Combining dictionaries
        print_('Combining dictionaries')

        result = {}
        for index in range(proc_count):
            with open(fr'{projectDir}\Cavities\shape_space{index}.json', "r") as infile:
                result.update(json.load(infile))

        # check if extension is included
        if filename:
            if filename.split('.')[-1] != 'json':
                filename = f'{filename}.json'
        else:
            return result

        with open(fr'{projectDir}\Cavities\{filename}', "w") as outfile:
            json.dump(result, outfile, indent=4, separators=(',', ': '))

        print_('Done combining dictionaries')
        return result

    @staticmethod
    def delete_process_dict(proc_count, projectDir):
        for index in range(proc_count):
            os.remove(fr'{projectDir}\Cavities\shape_space{index}.json')


class MonitorConvergence(QThread):
    sig = QtCore.pyqtSignal(list)

    def __init__(self, frame):
        super(QThread, self).__init__()
        self.frame = frame

    def run(self):
        self.monitor_convergence()

    def monitor_convergence(self):
        while not self.frame.tune_ended:
            self.sig.emit(self.frame.convergence_list._getvalue())


class MathTextLabel(QWidget):
    def __init__(self, mathText, parent=None):
        super(QWidget, self).__init__(parent)

        ll = QVBoxLayout(self)
        ll.setContentsMargins(0, 0, 0, 0)

        r, g, b, a = self.palette().base().color().getRgbF()

        self._figure = Figure(edgecolor=(r, g, b), facecolor=(r, g, b))
        self._canvas = FigureCanvas(self._figure)
        ll.addWidget(self._canvas)
        self._figure.clear()
        text = self._figure.suptitle(
            mathText,
            x=0.0,
            y=1.0,
            horizontalalignment='left',
            verticalalignment='top',
            size=QFont().pointSize() * 2
        )
        self._canvas.draw()

        (x0, y0), (x1, y1) = text.get_window_extent().get_points()
        w = x1 - x0
        h = y1 - y0

        self._figure.set_size_inches(w / 80, h / 80)
        self.setFixedSize(w, h)

