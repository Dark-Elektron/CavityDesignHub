from PyQt5 import QtGui
from scipy.special import jn_zeros, jnp_zeros
from ui_files.plot import Ui_Plot
from ui_files.plottypeselector import Ui_PlotTypeSelector
from modules.plot_module.plotter import Plot
from modules.data_module.abci_data import ABCIData
import pandas as pd
from utils.file_reader import FileReader
from utils.shared_classes import *
from utils.shared_functions import *

fr = FileReader()

file_color = 'yellow'
DEBUG = True


def print_(*arg):
    if DEBUG: print(colored(f'\t{arg}', file_color))


class PlotControl:
    def __init__(self, parent):
        print_("Check 1: plot_control.py")
        self.w_Plot = QWidget()

        self.ui = Ui_Plot()
        print_("Check 2: plot_control.py")
        self.ui.setupUi(self.w_Plot)
        print_("Check 3: plot_control.py")

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui
        print_("Check 4: plot_control.py")

        # get logger
        self.log = self.main_control.log

        # Plot objects and variables
        self.plt = Plot(self.ui)
        print_("Check 5: plot_control.py")
        self.ui.gl_Plot_Area.addWidget(self.plt)
        self.fig = self.plt.fig
        self.ax = self.plt.ax
        self.ax_right = self.plt.ax_right
        self.axins = None
        self.indicate_inset = None
        self.leg = None

        print_("Check 6: plot_control.py")

        # class variables
        self.plotID_count = 0

        # class lists
        self.freq_glob = []
        self.baseline_line_objects = []
        print_("Check 7: plot_control.py")

        # class dictionaries
        self.pts_dict = {}  # holds the plot information and plot
        # structure self.pts = {key1: {plot_widget: value, plot: value}}
        self.plot_dict = {}
        self.f_list_sorted = {}
        self.mode_list_sorted = {}
        self.args_dict = {'Code': [],
                          'Folder': [],
                          'Polarization': [],
                          'Id': [],
                          'Request': [],
                          'Toggle': [],
                          'ScaleX': [],
                          'ScaleY': [],
                          'Axis': [],
                          'Type': [],
                          'Filter': []
                          }
        self.args_dict_clone = {}
        self.other_data = {}
        self.ax_obj_dict = {}
        # baseline matplotlib line.Line2D objects
        self.axes_decorations_dict = {"xlabel": "x",
                                      "ylabel": "y",
                                      "ylabel2": "",
                                      "title": "Title",
                                      "legend": []}
        self.other_data_filtered = {}

        self.initUI()
        self.signals()

        print_("Check 3: plot_control.py")

    def createPlotTypeWidget(self):
        # plottypeselector
        w_PTS = QWidget()
        PTS = Ui_PlotTypeSelector()
        PTS.setupUi(w_PTS)

        plotID = self.plotID_count
        # add plottypse selector to dictionary
        self.pts_dict[self.plotID_count] = {"plot inputs": None, "plot data": {}, "plot object": {},
                                            "plot data inputs": None}

        self.init_abci(PTS, plotID)
        self.init_slans(PTS, plotID)
        self.init_other(PTS, plotID)

        # add plottypeselector to layout
        self.ui.gl_Plot_Type_Selector.addWidget(w_PTS, 0, self.plotID_count, 1, 1)
        self.pts_dict[self.plotID_count]["plot inputs"] = PTS

        # add delete signal to delete buttion
        PTS.pb_Remove.clicked.connect(lambda: self.remove_pts(plotID))

        # add plot signal to switch
        PTS.switchControl.stateChanged.connect(lambda: self.make_plot_object(PTS, plotID))

        # increment plotID
        self.plotID_count += 1

    def init_abci(self, PTS, plotID):

        # abci
        PTS.pb_Open_Folder_ABCI.clicked.connect(lambda: self.make_abci_data(PTS))
        request_dict_long = {0: 'Longitudinal Impedance Magnitude',
                             1: r'Cavity Shape Input',
                             2: r'Cavity Shape Used',
                             3: r'Wake Potentials',
                             4: r'Real Part of Longitudinal Impedance',
                             5: r'Imaginary Part of Longitudinal Impedance',
                             6: r'Frequency Spectrum of Loss Factor',
                             7: r'Loss Factor Spectrum Integrated up to F',
                             8: r'Real Part of Long. + Log Impedance',
                             9: r'Imaginary Part of Long. + Log Impedance',
                             10: r'Spectrum of Long. + Log Loss Factor',
                             11: r'Long. + Log Factor Integrated up to F'
                             }

        request_dict_trans = {0: 'Transverse Impedance Magnitude',
                              1: r'Cavity Shape Input',
                              2: r'Cavity Shape Used',
                              3: r'Wake Potentials',
                              4: r'Real Part of Azimuthal Impedance',
                              5: r'Imaginary Part of Azimuthal Impedance',
                              6: r'Real Part of Transverse Impedance',
                              7: r'Imaginary Part of Transverse Impedance',
                              8: r'Real Part of Longitudinal Impedance',
                              9: r'Imaginary Part of Longitudinal Impedance'
                              }

        for req in request_dict_long.values():
            PTS.cb_Request_ABCI.addItem(req)

        PTS.cb_Polarization_ABCI.currentIndexChanged.connect(
            lambda: self.switchRequest(PTS.cb_Request_ABCI, request_dict_long, request_dict_trans,
                                       PTS.cb_Polarization_ABCI))

        type_list = ['Line', 'Scatter', 'Bar']

        for a in type_list:
            PTS.cb_Type_ABCI.addItem(a)

        line_style = ['-', '--', ':', '-.']
        scatter_marker = ['x', 'o', '+', 'P', 'X', 'd', '>', '^', 'v', 'H']
        # fill default
        for a in line_style:
            PTS.cb_Style_ABCI.addItem(a)

        # signal for combobox
        PTS.cb_Type_ABCI.currentIndexChanged.connect(
            lambda: self.populate_combobox_list(PTS.cb_Style_ABCI, line_style)
            if PTS.cb_Type_ABCI.currentText() == "Line"
            else self.populate_combobox_list(PTS.cb_Style_ABCI, scatter_marker))

    def init_slans(self, PTS, plotID):
        # slans
        PTS.pb_Open_Folder_SLANS.clicked.connect(
            lambda: self.open_(PTS.le_Folder_SLANS, PTS.cb_Polarization_SLANS, PTS.ccb_Shape_ID_SLANS, mode='Folder'))

    def init_other(self, PTS, plotID):
        # other
        PTS.pb_Open_File_Other.clicked.connect(lambda: self.open_(PTS.le_Filename, mode='File'))
        PTS.le_Filename.textChanged.connect(
            lambda: self.populate_combobox_new([PTS.ccb_Id_Other, PTS.cb_X, PTS.ccb_Y, PTS.ccb_Filter_Other],
                                               PTS.le_Filename.text(), plotID))

    def make_abci_data(self, PTS):
        self.open_(PTS.le_Folder_ABCI, PTS.cb_Polarization_ABCI, PTS.ccb_Id_ABCI, mode='Folder')

        # abci data is not made on opening because it would be memory consuming to load all simulation results
        # for this reason, the data is made when the plot function is called

    def make_slans_data(self, PTS):
        pass

    def make_other_data(self, PTS):
        pass

    def make_plot_object(self, PTS, plotID):
        if PTS.cb_code.currentText() == "abci":
            self.make_abci_plot_object(plotID)
        elif PTS.cb_code.currentText() == "slans":
            self.make_slans_plot_object(plotID)
        else:
            self.make_other_plot_object(plotID)

    def remove_pts(self, pid):
        for line2D in self.pts_dict[pid]['plot object'].values():
            line2D[0].remove()
        self.pts_dict[pid]["plot inputs"].pb_Remove.clicked.connect(lambda: self.pts_dict.pop(pid))

        self.update_labels()
        self.fig.canvas.draw_idle()

    def make_abci_plot_object(self, pid):
        # check state of switch
        switch = self.pts_dict[pid]["plot inputs"].switchControl
        plottype = self.pts_dict[pid]["plot inputs"].cb_code.currentText()
        if switch.checkState() == 2:

            # check if plot already exist
            if self.pts_dict[pid]['plot data'] != {}:
                # compare plot input with current input
                inputs = self.pts_dict[pid]["plot inputs"]
                ids = [a.strip() for a in inputs.ccb_Id.currentText().split(',')]  # get list
                pol = inputs.cb_Pol.currentText()
                request = inputs.cb_Request.currentText()

                plot_data_inputs = self.pts_dict[pid]["plot data inputs"]

                if plot_data_inputs["ids"] == ids and plot_data_inputs["pol"] == pol and plot_data_inputs[
                    "request"] == request:
                    if self.pts_dict[pid]["plot object"] != {}:
                        for line2D in self.pts_dict[pid]['plot object'].values():
                            line2D[0].set(alpha=1)
                        self.fig.canvas.draw_idle()
                    else:
                        # plot with modified inputs
                        self.plot_impedance_pts(self.pts_dict[pid])
                else:
                    # pop previous entry from plot and dictionary
                    for line2D in self.pts_dict[pid]['plot object'].values():
                        line2D[0].remove()

                    # reset plot data, plot object and plot data inputs
                    self.pts_dict[pid]["plot data"] = {}
                    self.pts_dict[pid]["plot object"] = {}
                    self.pts_dict[pid]["plot data inputs"] = {}

                    self.plot_impedance_pts(self.pts_dict[pid])

            else:
                # self.ax.clear()
                # make plot
                self.plot_impedance_pts(self.pts_dict[pid])

        else:
            for line2D in self.pts_dict[pid]['plot object'].values():
                line2D[0].set(alpha=0)

        self.ax.set_prop_cycle(None)
        cycle = self.ax._get_lines.prop_cycler

        # update colors, loop over all lines and give color in ascending order
        for line2D in self.ax.get_lines():
            line2D.set_color(next(cycle)['color'])

        self.update_labels()
        self.fig.canvas.draw_idle()

    def make_slans_plot_object(self, pid):
        pass

    def make_other_plot_object(self, pid):
        self.plot_other_new(pid)

    def plot_impedance_pts(self, pts_):
        inputs = pts_["plot inputs"]
        ids = [a.strip() for a in inputs.ccb_Id.currentText().split(',')]  # get list
        pol = inputs.cb_Pol.currentText()
        request = inputs.cb_Request.currentText()
        folder = inputs.le_Folder.text()
        axis = inputs.cb_Axis.currentText()

        try:
            scaleX = float(inputs.le_ScaleY.text())
            scaleY = float(inputs.le_ScaleY.text())
        except:
            scaleX = 1
            scaleY = 1

        if folder == '':
            abci_data_dir = fr'{self.main_control.projectDir}/SimulationData/ABCI'
        else:
            abci_data_dir = folder

        # update data input
        pts_["plot data inputs"] = {"ids": ids, "pol": pol, "request": request}

        # plot for multiple ids
        for id_ in ids:
            if pol.lower() == 'long':
                abci_data_long = ABCIData(abci_data_dir, id_, 0)

                if request == 'Longitudinal Impedance Magnitude':
                    xr, yr, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                    xi, yi, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')

                    y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_long.get_data(request)

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # reset global frequency if the impedance magnitude is not requested
                    self.freq_glob = []

                # plot
                ax_selected = self.ax
                if axis.lower() == 'left':
                    ax_selected = self.ax
                    # update plot data
                    pts_["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    pts_["plot object"].update({id_: self.ax.plot(xr, y,
                                                                  label='${' + f'{str(id_).replace("_", ",")}' + '}$' + f' ({abci_data_long.wakelength} m)',
                                                                  linewidth=2)})

                elif axis.lower() == "right":
                    ax_selected = self.ax_right
                    # update plot data
                    pts_["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    pts_["plot object"].update({id_: self.ax_right.plot(xr, y,
                                                                        label='${' + f'{str(id_).replace("_", ",")}' + '}$' + f' ({abci_data_long.wakelength} m)',
                                                                        linewidth=2)})

                # ax_selected.set_xlabel('$f \mathrm{ [GHz]}$')
                # ax_selected.set_ylabel('$Z_{\parallel, \mathrm{HOM}} \mathrm{[k\Omega]}$')
                ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

            else:
                abci_data_trans = ABCIData(abci_data_dir, id_, 1)
                if request == 'Transverse Impedance Magnitude':
                    try:
                        xr, yr, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                        xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')
                    except:
                        xr, yr, _ = abci_data_trans.get_data('Real Part of Azimuthal Impedance')
                        xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Azimuthal Impedance')

                    y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_trans.get_data(request)

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # reset global frequency if the impedance magnitude is not requested
                    self.freq_glob = []

                # plot
                ax_selected = self.ax
                if axis.lower() == 'left':
                    ax_selected = self.ax
                    # update plot data
                    pts_["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    pts_["plot object"].update({id_: self.ax.plot(xr, y,
                                                                  label='${' + f'{str(id_).replace("_", ",")}' + '}$' + f' ({abci_data_trans.wakelength} m)',
                                                                  linewidth=2)})

                elif axis.lower() == 'right':
                    ax_selected = self.ax_right
                    # update plot data
                    pts_["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    pts_["plot object"].update({id_: self.ax_right.plot(xr, y,
                                                                        label='${' + f'{str(id_).replace("_", ",")}' + '}$' + f' ({abci_data_trans.wakelength} m)',
                                                                        linewidth=2)})

                ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

    def signals(self):
        self.ui.pb_Add.clicked.connect(lambda: self.createPlotTypeWidget())

        # plot abci impedance
        self.ui.pb_Plot.clicked.connect(lambda: self.plot())

        # refresh plot
        self.ui.pb_Refresh.clicked.connect(lambda: self.plot())

        # toggle plot menu
        self.ui.pb_Plot_Area_Menu.clicked.connect(
            lambda: self.main_control.animate_width(self.ui.w_Plot_Menu, 0, 440, True))
        self.ui.pb_Plot_Area_Menu.clicked.connect(lambda: self.ui.w_Plot_Area_Buttons.setEnabled(
            True) if self.ui.w_Plot_Menu.maximumWidth() == 0 else self.ui.w_Plot_Area_Buttons.setDisabled(True))

        # signal for plot menu pages
        self.ui.pb_Machine_Parameters.clicked.connect(lambda: self.toggle_page('Machine Parameters'))
        self.ui.pb_Plot_Decorations.clicked.connect(lambda: self.toggle_page('Plot Decorations'))

        # signal for plot argument entrmake widgets rey
        # self.ui.pb_Collapse_Shape_Parameters.clicked.connect(lambda: self.main_control.animate_height(self.ui.w_Plot_Args, 0, 200, True))

        # signal for threshold
        self.ui.cb_Longitudinal_Threshold.stateChanged.connect(lambda: self.calc_limits('monopole'))
        self.ui.cb_Transverse_Threshold.stateChanged.connect(lambda: self.calc_limits('dipole'))

        # signal for inset plot
        self.ui.cb_Show_Hide_Inset.clicked.connect(lambda: self.plot_inset())

        # add plot
        self.ui.pb_Add_Row.clicked.connect(self.add_row)
        # clear plots
        self.ui.pb_Clear.clicked.connect(lambda: self.clear_plots())

        # signals for cutoff
        # self.ui.cb_Check_1.clicked.connect(lambda: self.plot_cutoff(0, self.ui.cb_Check_1))
        # self.ui.cb_Check_2.clicked.connect(lambda: self.plot_cutoff(1, self.ui.cb_Check_2))
        # self.ui.cb_Check_3.clicked.connect(lambda: self.plot_cutoff(2, self.ui.cb_Check_1))

        # update label
        self.ui.pb_Apply_Axis_Labels.clicked.connect(lambda: self.update_labels())

        # plot text
        self.ui.pb_Add_Textt.clicked.connect(lambda: self.plt.add_text(
            text=self.ui.le_Plot_Text_2.text(), box=self.ui.cb_Text_Box.currentText(),
            size=self.ui.sb_Annotation_Text_Size.value(), rotation=self.ui.sb_Rotation.value()))

        # plot axvline
        self.ui.pb_Add_Vline.clicked.connect(lambda: self.plt.add_axvline(self.ui.dsb_Axvline_X.value()))

        # switch axis
        self.ui.cb_Active_Axis.currentTextChanged.connect(lambda: self.switch_axis())

        # frequently used plot labels combobox signal
        self.ui.ccb_Freq_Used_Xlabel.currentTextChanged.connect(lambda: self.ui.le_Xlabel.setText(
            self.ui.ccb_Freq_Used_Xlabel.currentText()))
        self.ui.ccb_Freq_Used_Ylabel.currentTextChanged.connect(lambda: self.ui.le_Ylabel.setText(
            self.ui.ccb_Freq_Used_Ylabel.currentText()))

    def initUI(self):
        self.createPlotTypeWidget()

        print_("Check 5: plot_control.py")
        # add row to table
        self.ui.tableWidget.setRowCount(1)  # and one row in the table
        self.table_control()

        print_("Check 5: plot_control.py")
        # set plot menu initial size to zero and disable plot buttons
        self.ui.w_Plot_Menu.setFixedWidth(0)
        self.ui.w_Plot_Area_Buttons.setDisabled(True)

        print_("Check 5: plot_control.py")
        # tableWidget initializations
        self.ui.tableWidget.mousePressEvent = self.mousePressEvent
        self.ui.tableWidget.setContextMenuPolicy(Qt.CustomContextMenu)
        self.ui.tableWidget.customContextMenuRequested.connect(self.generateMenu)
        self.set_table_size()

        print_("Check 5: plot_control.py")
        # add checkable combobox
        cutoff = QCheckableComboBox()
        cutoff.addItem("All")
        cutoff.setMinimumWidth(150)
        self.ui.gl_Cutoff.addWidget(cutoff, 2, 0, 1, 1)
        self.populate_cutoff_combobox(cutoff)
        print_("Check 5: plot_control.py")

        self.ui.le_Ri_Cutoff_List.editingFinished.connect(lambda: self.populate_cutoff_combobox(cutoff))

        # add signal
        cutoff.currentTextChanged.connect(lambda: self.plot_cutoff(cutoff))

    def plot(self):
        try:
            # args = list(self.args_dict.values())

            # use same color cycler for both axes
            self.ax_right._get_lines.prop_cycler = self.ax._get_lines.prop_cycler
            plot_count = 1
            for key, val in self.plot_dict.items():
                code = self.plot_dict[key]['plot inputs']['Code'].currentText()
                # check plot type
                if code == 'ABCI':
                    self.make_abci_plot(key)
                    # self.plot_impedance(self.args_dict, i, plot_count)
                elif code == 'SLANS':
                    pass
                else:
                    self.make_other_plot(key)

            # toggle axis labels
            self.toggle_axis_labels()

            self.ax.autoscale(True, axis='y')
            self.ax_right.autoscale(True, axis='y')

            # recompute the ax.datalim
            self.ax.relim()
            self.ax_right.relim()

            # show legend
            lines, labels = self.ax.get_legend_handles_labels()
            lines2, labels2 = self.ax_right.get_legend_handles_labels()

            if self.ui.cb_Active_Axis.currentText() == 'Left':
                print_("\tInside here to plot legend on left axis")
                self.leg = self.ax.legend(lines + lines2, labels + labels2, loc='lower left', prop={'size': 18})

                if self.ax_right.get_legend() is not None:
                    self.ax_right.get_legend().remove()
            else:
                self.leg = self.ax_right.legend(lines + lines2, labels + labels2, loc='lower left', prop={'size': 18})
                if self.ax.get_legend() is not None:
                    self.ax.get_legend().remove()

            self.leg.set_zorder(10)
            self.leg.set_draggable(state=True, use_blit=True)

            # plot inset if check box is checked
            self.plot_inset()
            # self.fig.canvas.draw_idle()

            # plot thresholds if threshold is checked
            if self.ui.cb_Longitudinal_Threshold.checkState() == 2:
                self.calc_limits('monopole')
            if self.ui.cb_Transverse_Threshold.checkState() == 2:
                self.calc_limits('dipole')

            print_("it at the very very least got herer")
            # self.switch_axis()

        except Exception as e:
            self.log.error("Please enter a valid argument: Exception: ", e)

    def clear_plots(self):
        self.ax.cla()
        self.ax_right.cla()

        if self.axins is not None:
            self.axins.cla()
            self.axins.remove()
            self.axins = None

        # reset plot dict
        for key, val in self.plot_dict.items():
            self.plot_dict[key]["plot data"] = {}
            self.plot_dict[key]["plot object"] = {}

        # clear annotations
        self.plt.clear()

        # clear cutoff frequency lines
        for key, val in self.ax_obj_dict.items():
            del val
        self.ax_obj_dict = {}

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def switch_axis(self):
        # change order of appearance of axis
        if self.ui.cb_Active_Axis.currentText() == "Left":
            self.ax_right.set_zorder(0)
            self.ax.set_zorder(1)

            # turn at least one frame on, the axis with zorder = 0
            self.ax.set_frame_on(False)
            self.ax_right.set_frame_on(True)
        else:
            self.ax.set_zorder(0)
            self.ax_right.set_zorder(1)

            # turn at least one frame on, the axis with zorder = 0
            self.ax_right.set_frame_on(False)
            self.ax.set_frame_on(True)

        self.update_labels()

    def plot_impedance(self, plot_dict):  # , i, plot_count
        args = plot_dict["plot inputs"]
        ids = [a.strip() for a in args['Id'].currentText().split(',')]  # get list
        pol = args['Polarization'].currentText()
        request = args['Request'][0].currentText()
        folder = args['Folder'][0].text()
        state = args['Toggle'].checkState()
        axis = args['Axis'].currentText()
        type_ = []
        t = args['Type']
        type_.append(t[0].currentText())
        type_.append(t[1].currentText())

        filter_ = []
        f = args['Filter']
        filter_.append(f[0].currentText())
        filter_.append(f[1].text())

        # try:
        scaleX = float(args['ScaleX'].text())
        scaleY = float(args['ScaleY'].text())
        # except:
        #     scaleX = 1
        #     scaleY = 1

        print_("hwere 333wwsdfwerwe")
        if folder == '':
            abci_data_dir = fr'{self.main_control.projectDir}/SimulationData/ABCI'
        else:
            abci_data_dir = folder

        # record data input
        plot_dict["plot data inputs"] = {"Ids": ids, "Polarization": pol, "Request": request, "Folder": folder,
                                         "Axis": axis, "ScaleX": scaleX, "ScaleY": scaleY, "Type": type_,
                                         "Filter": filter_}

        # plot for multiple ids
        for id_ in ids:
            if pol.lower() == 'long':
                abci_data_long = ABCIData(abci_data_dir, id_, 0)

                if request == 'Longitudinal Impedance Magnitude':
                    xr, yr, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                    xi, yi, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')

                    y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_long.get_data(request)

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # reset global frequency if the impedance magnitude is not requested
                    self.freq_glob = []

                # plot
                ax_selected = self.ax
                if axis == 'Left':
                    ax_selected = self.ax

                    # update plot data
                    plot_dict["plot data"].update({id_: {"x": xr, "y": y}})

                    # # process id
                    # self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                    # f' ({abci_data_long.wakelength} m)', linewidth=2)
                    # process id
                    plot_dict["plot object"].update({
                        id_: self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                       f' ({abci_data_long.wakelength} m)',
                                          linewidth=2, picker=True)})
                    # mplcursors.cursor(plot_dict["plot object"][id_])

                elif axis == "Right":
                    ax_selected = self.ax_right
                    # self.ax_right.plot(xr, y,
                    #                    label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                    #                    f' ({abci_data_long.wakelength} m)',
                    #                    linewidth=2)
                    # update plot data
                    plot_dict["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    plot_dict["plot object"].update({
                        id_: self.ax_right.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                             f' ({abci_data_long.wakelength} m)',
                                                linewidth=2, picker=True)})
                    # mplcursors.cursor(plot_dict["plot object"][id_])

                # ax_selected.set_xlabel('$f \mathrm{ [GHz]}$')
                # ax_selected.set_ylabel('$Z_{\parallel, \mathrm{HOM}} \mathrm{[k\Omega]}$')
                ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

            else:
                abci_data_trans = ABCIData(abci_data_dir, id_, 1)
                if request == 'Transverse Impedance Magnitude':
                    try:
                        xr, yr, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                        xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')
                    except:
                        xr, yr, _ = abci_data_trans.get_data('Real Part of Azimuthal Impedance')
                        xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Azimuthal Impedance')

                    y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_trans.get_data(request)

                    # scale x axis]
                    xr = [a * scaleX for a in xr]
                    # scale y axis]
                    y = [a * scaleY for a in y]

                    # reset global frequency if the impedance magnitude is not requested
                    self.freq_glob = []

                # plot
                ax_selected = self.ax
                if axis == 'Left':
                    ax_selected = self.ax
                    # self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                    # f' ({abci_data_trans.wakelength} m)', linewidth=2)
                    # update plot data
                    plot_dict["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    plot_dict["plot object"].update({
                        id_: self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                       f' ({abci_data_trans.wakelength} m)',
                                          linewidth=2, picker=True)})
                    # mplcursors.cursor(plot_dict["plot object"][id_])
                elif axis == 'Right':
                    ax_selected = self.ax_right
                    # self.ax_right.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                    # f' ({abci_data_trans.wakelength} m)', linewidth=2)
                    # update plot data
                    plot_dict["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    plot_dict["plot object"].update({
                        id_: self.ax_right.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                             f' ({abci_data_trans.wakelength} m)',
                                                linewidth=2, picker=True)})
                    # mplcursors.cursor(plot_dict["plot object"][id_])

                ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

            # plot axes labels
            self.update_labels()

            # increment plot count
            # plot_count += 1

    def make_abci_plot(self, pid):
        # check state of switch
        print_("It's here000")
        switch = self.plot_dict[pid]["plot inputs"]["Toggle"]
        plottype = self.plot_dict[pid]["plot inputs"]["Code"].currentText()
        if switch.checkState() == 2:
            print_("It's herel")
            # check if plot already exist
            # print_(self.plot_dict[pid]['plot data'])
            if self.plot_dict[pid]['plot data'] != {}:
                print_("It's here3")
                # compare plot input with current input
                args = self.plot_dict[pid]["plot inputs"]
                ids = [a.strip() for a in args['Id'].currentText().split(',')]  # get list
                pol = args['Polarization'].currentText()
                request = args['Request'][0].currentText()
                folder = args['Folder'][0].text()
                state = args['Toggle'].checkState()
                axis = args['Axis'].currentText()

                type_ = []
                filter_ = []
                t = args['Type']

                type_.append(t[0].currentText())
                type_.append(t[1].currentText())

                f = args['Filter']
                filter_.append(f[0].currentText())
                filter_.append(f[1].text())

                try:
                    scaleX = float(args['ScaleX'].text())
                    scaleY = float(args['ScaleY'].text())
                except:
                    scaleX = 1
                    scaleY = 1

                inputs_compare = {"Ids": ids, "Polarization": pol, "Request": request, "Folder": folder,
                                  "Axis": axis, "ScaleX": scaleX, "ScaleY": scaleY, "Type": type_, "Filter": filter_}

                plot_data_inputs = self.plot_dict[pid]["plot data inputs"]

                if plot_data_inputs == inputs_compare:
                    if self.plot_dict[pid]["plot object"] != {}:
                        for line2D in self.plot_dict[pid]['plot object'].values():
                            line2D[0].set(alpha=1)
                        self.fig.canvas.draw_idle()
                    else:
                        # plot with modified inputs
                        self.plot_impedance(self.plot_dict[pid])
                else:
                    # pop previous entry from plot and dictionary
                    for line2D in self.plot_dict[pid]['plot object'].values():
                        line2D[0].remove()

                    # reset plot data, plot object and plot data inputs
                    self.plot_dict[pid]["plot data"] = {}
                    self.plot_dict[pid]["plot object"] = {}
                    self.plot_dict[pid]["plot data inputs"] = {}
                    self.plot_impedance(self.plot_dict[pid])

            else:
                # self.ax.clear()
                # make plot
                print_("Here to make plot")
                self.plot_impedance(self.plot_dict[pid])
                print_("Here to make plot dfgf")

        else:
            for line2D in self.plot_dict[pid]['plot object'].values():
                line2D[0].set(alpha=0)

        self.ax.set_prop_cycle(None)
        cycle = self.ax._get_lines.prop_cycler

        # update colors, loop over all lines and give color in ascending order
        for line2D in self.ax.get_lines():
            line2D.set_color(next(cycle)['color'])

        self.update_labels()
        self.fig.canvas.draw_idle()

    def make_other_plot(self, pid):
        # check state of switch
        print_("It's here000")
        switch = self.plot_dict[pid]["plot inputs"]["Toggle"]
        plottype = self.plot_dict[pid]["plot inputs"]["Code"].currentText()
        if switch.checkState() == 2:
            print_("It's herel")
            # check if plot already exist
            # print_(self.plot_dict[pid]['plot data'])
            if self.plot_dict[pid]['plot data'] != {}:
                print_("It's here3")
                # compare plot input with current input
                args = self.plot_dict[pid]["plot inputs"]
                ids = [a.strip() for a in args['Id'].currentText().split(',')]  # get list
                filename = args['Folder'][0].text()
                # pol = args['Polarization'].currentText()
                # request = args['Request'][0].currentText()
                print_("Now hwere wsdfad")
                # folder = args['Folder'][0].text()
                state = args['Toggle'].checkState()
                axis = args['Axis'].currentText()

                print_("Now hwere wsdfad")
                requestX = args['Request'][1].currentText()
                requestY = args['Request'][2].currentText().split(', ')

                print_("Now hwere wsdfad")
                type_ = []
                filter_ = []
                t = args['Type']
                type_.append(t[0].currentText())
                type_.append(t[1].currentText())
                print_("Now hwere wsdfad")

                f = args['Filter']
                filter_.append(f[0].currentText())
                filter_.append(f[1].text())

                try:
                    scaleX = float(args['ScaleX'].text())
                    scaleY = float(args['ScaleY'].text())
                except:
                    scaleX = 1
                    scaleY = 1

                print_("Now hwere wsdfad")

                inputs_compare = {"Ids": ids, "RequestX": requestX, "RequestY": requestY, "Folder": filename,
                                  "Axis": axis, "ScaleX": scaleX, "ScaleY": scaleY, "Type": type_, "Filter": filter_}

                plot_data_inputs = self.plot_dict[pid]["plot data inputs"]

                if plot_data_inputs == inputs_compare:
                    print_("It's nowhere")
                    if self.plot_dict[pid]["plot object"] != {}:
                        for line2D in self.plot_dict[pid]['plot object'].values():
                            line2D[0].set(alpha=1)
                        self.fig.canvas.draw_idle()
                    else:
                        # plot with modified inputs
                        self.plot_other(self.plot_dict[pid], pid)
                else:
                    # pop previous entry from plot and dictionary
                    print_("It's her ein alternative")
                    for vals in self.plot_dict[pid]['plot object'].values():
                        for line2D in vals.values():
                            line2D[0].remove()

                    # reset plot data, plot object and plot data inputs
                    self.plot_dict[pid]["plot data"] = {}
                    self.plot_dict[pid]["plot object"] = {}
                    self.plot_dict[pid]["plot data inputs"] = {}
                    print_("Alternative here here")
                    self.plot_other(self.plot_dict[pid], pid)

            else:
                # self.ax.clear()
                # make plot
                print_("Here to make plot")
                self.plot_other(self.plot_dict[pid], pid)
                print_("Here to make plot dfgf")

        else:
            for vals in self.plot_dict[pid]['plot object'].values():
                for line2D in vals.values():
                    line2D[0].set(alpha=0)

        self.ax.set_prop_cycle(None)
        cycle = self.ax._get_lines.prop_cycler

        # update colors, loop over all lines and give color in ascending order
        for line2D in self.ax.get_lines():
            line2D.set_color(next(cycle)['color'])

        self.update_labels()
        self.fig.canvas.draw_idle()

    def plot_inset(self):
        if self.ui.cb_Show_Hide_Inset.checkState() == 2:
            # get lines from axis
            lines = self.ax.get_lines()

            inset_pos = ast.literal_eval(self.ui.le_Inset_Position.text())
            if len(inset_pos) == 4:
                self.axins = self.ax.inset_axes(inset_pos, facecolor='#fafafa')
                # self.axins.axes.xaxis.set_visible(False)
                # self.axins.axes.yaxis.set_visible(False)

            for line in lines:
                if self.axins:
                    if line.get_linestyle() == 'None':
                        self.axins.plot(line.get_xdata(), line.get_ydata(), linestyle='None', marker=line.get_marker(),
                                        markersize=line.get_ms(), c=line.get_color(), mec=line.get_mec())
                    else:
                        self.axins.plot(line.get_xdata(), line.get_ydata(), ls=line.get_linestyle(),
                                        linewidth=line.get_lw(), c=line.get_color())

            # sub region of the original image
            # get values from line edit
            try:
                x0, x1, y0, y1 = ast.literal_eval(self.ui.le_Inset_Window.text())
            except:
                x0, x1, y0, y1 = 0.385, 0.415, 5e-1, 3e1  # default

            self.axins.set_xlim(x0, x1)
            self.axins.set_ylim(y0, y1)

            self.axins.set_xticklabels('')
            self.axins.set_yticklabels('')

            self.axins.set_yscale(self.ui.cb_Y_Scale.currentText())
            self.indicate_inset = self.ax.indicate_inset_zoom(self.axins, edgecolor="black", label='_nolegend_')
        else:
            if self.indicate_inset:
                for x in self.indicate_inset:
                    if isinstance(x, tuple):
                        for y in x:
                            y.remove()
                    else:
                        x.remove()

            if self.axins:
                self.axins.cla()
                self.axins.remove()
                self.axins = None

        self.fig.canvas.draw_idle()

    def plot_other_new(self, pid):
        # print_(self.other_data)
        # print_(i)

        requestX = self.pts_dict[pid]["plot inputs"].cb_X.currentText()
        requestY = self.pts_dict[pid]["plot inputs"].cb_Y.currentText().split(', ')

        filename = self.pts_dict[pid]["plot inputs"].le_Filename.text()
        # state = args['Toggle'][i].checkState()
        axis = self.pts_dict[pid]["plot inputs"].cb_Axis_Other.currentText()
        type_ = self.pts_dict[pid]["plot inputs"].cb_Type_Other.currentText()
        style = self.pts_dict[pid]["plot inputs"].cb_Style.currentText()

        try:
            scaleX = float(self.pts_dict[pid]["plot inputs"].le_ScaleX_Other.text())
            scaleY = float(self.pts_dict[pid]["plot inputs"].cb_ScaleY.text())
        except Exception as e:
            print_(e)
            scaleX = 1
            scaleY = 1

        if filename == '':
            folder = fr'{self.main_control.projectDir}/SimulationData/ABCI'

        # # load file (for now, just excel files)
        # df = fr.excel_reader(filename)
        # sheet_name = list(df.keys())[0]
        #
        # data = df[sheet_name]

        if requestY != [] and requestX != []:
            if filename.split('.')[-1] == 'xlsx':
                print_("xlsx", filename)
                # get sheets
                sheets = [a.strip() for a in self.pts_dict[pid]["plot inputs"].ccb_Id_Other.currentText().split(',')]
                print_("1")
                for sh in sheets:
                    print_("2", sh)

                    print_("3")
                    filter_ = self.pts_dict[pid]["plot inputs"].ccb_Filter_Other.currentText()
                    value = self.pts_dict[pid]["plot inputs"].le_Filter_Value_Other.text()

                    if filter_ == "None" or value == "":
                        print_("3a")
                        self.other_data_filtered = self.pts_dict[pid]["plot data"][sh]
                        print_("4")
                    else:
                        print_("4a")
                        self.other_data_filtered = self.pts_dict[pid]["plot data"][sh][
                            self.pts_dict[pid]["plot data"][sh][filter_] == value]
                        print_("4b")

                    x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                    self.freq_glob = x_data

                    for j in range(len(requestY)):
                        y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                        print_("xlxx: ", x_data, y)

                        if axis == 'Left':
                            if type_ == 'Line':
                                self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                                print_("done plotting xlsx")
                            else:
                                self.ax.plot(x_data, y, linestyle='None', marker=style, markersize=10.0,
                                             markeredgecolor="black", label="Legend")
                            self.ax.set_ylabel('$Y$ []')
                            self.ax.set_xlabel('$X$ []')
                        else:
                            if type_ == 'Line':
                                self.ax_right.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                            else:
                                self.ax_right.plot(x_data, y, linestyle='None', markeredgecolor="black", marker=style,
                                                   markersize=10.0,
                                                   label="Legend")
                            self.ax_right.set_ylabel('$Y$ [dB]')
            else:
                print_("txt", filename)
                # try to filter self.other_data
                try:
                    print_("1a")
                    filter_ = self.pts_dict[pid]["plot inputs"].ccb_Filter_Other.currentText()
                    value = self.pts_dict[pid]["plot inputs"].le_Filter_Value_Other.text()
                    if filter_ == "None" or value == "":
                        print_("4a")
                        self.other_data_filtered = self.pts_dict[pid]["plot data"]
                        print_("5a")
                    else:
                        print_("2a")
                        self.other_data_filtered = self.pts_dict[pid]["plot data"][
                            self.pts_dict[pid]["plot data"][filter_] == value]
                        print_("3a")

                except:
                    print_("6a")
                    self.other_data_filtered = self.pts_dict[pid]["plot data"]
                    print_("7a")

                print_("8a")
                x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                self.freq_glob = x_data
                print_("9a")

                for j in range(len(requestY)):
                    print_("10a")
                    y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                    print_("11a")

                    print_("txt: ", x_data, y)
                    if axis == 'Left':
                        if type_ == 'Line':
                            self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                            print_("Done plotting txt")
                        else:
                            self.ax.plot(x_data, y, linestyle='None', marker=style, markersize=10.0,
                                         markeredgecolor="black",
                                         label="Legend")
                        self.ax.set_ylabel('$Y$ []')
                        self.ax.set_xlabel('$X$ []')
                    else:
                        if type_ == 'Line':
                            self.ax_right.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                        else:
                            self.ax_right.plot(x_data, y, linestyle='None', marker=style, markersize=10.0,
                                               markeredgecolor="black",
                                               label="Legend")
                        self.ax_right.set_ylabel('$Y$ [dB]')
            print_("#" * 50)
        else:
            print_("Please specify columns to plot")

    def plot_other(self, plot_dict, key):
        args = plot_dict
        # print_(self.other_data)
        # print_(i)
        requestX = args["plot inputs"]['Request'][1].currentText()
        requestY = args["plot inputs"]['Request'][2].currentText().split(', ')
        ids = [a.strip() for a in args["plot inputs"]['Id'].currentText().split(',')]  # get list
        filename = args["plot inputs"]['Folder'][0].text()
        state = args["plot inputs"]['Toggle'].checkState()
        axis = args["plot inputs"]['Axis'].currentText()
        type_ = args["plot inputs"]['Type'][0].currentText()
        style = args["plot inputs"]['Type'][1].currentText()
        filter_ = args["plot inputs"]['Filter'][0].currentText()
        value = args["plot inputs"]['Filter'][1].text()

        try:
            scaleX = float(args["plot inputs"]['ScaleX'].text())
            scaleY = float(args["plot inputs"]['ScaleY'].text())
        except Exception as e:
            print_(e)
            scaleX = 1
            scaleY = 1

        if filename == '':
            folder = fr'{self.main_control.projectDir}/SimulationData/ABCI'

        # record data input
        args["plot data inputs"] = {"Ids": ids, "RequestX": requestX, "RequestY": requestY, "Folder": filename,
                                    "Axis": axis, "ScaleX": scaleX, "ScaleY": scaleY, "Type": type_, "Filter": filter_}

        # # load file (for now, just excel files)
        # df = fr.excel_reader(filename)
        # sheet_name = list(df.keys())[0]
        #
        # data = df[sheet_name]

        sheets = [a.strip() for a in args["plot inputs"]['Id'].currentText().split(',')]
        for sh in sheets:
            for id_ in ids:
                args["plot data"][id_] = {}
                args["plot object"][id_] = {}
                if requestY != [] and requestX != []:
                    if filename.split('.')[-1] == 'xlsx':
                        print_("xlsx", filename)
                        # get sheets
                        print_("1")
                        print_("2", sh)
                        if filter_ == "None" or value == "":
                            print_("3a", self.other_data[key][sh])
                            self.other_data_filtered = self.other_data[key][sh]
                            print_("4")
                        else:
                            print_("4a")
                            self.other_data_filtered = self.other_data[key][sh][
                                self.other_data[key][sh][filter_] == value]
                            print_("4b")

                        x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                        self.freq_glob = x_data

                        for j in range(len(requestY)):
                            y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                            print_("xlxx: ", x_data, y)
                            print_(args["plot data"])
                            args["plot data"][id_].update({j: {"x": x_data, "y": y}})

                            if axis == 'Left':
                                if type_ == 'Line':
                                    print_("hsdadf")
                                    args["plot object"][id_].update(
                                        {j: self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)})
                                    print_("done plotting xlsx")
                                else:
                                    print_("y: ", y)
                                    args["plot object"][id_].update({j: self.ax.plot(x_data, y, linestyle='None',
                                                                                     marker=style, markersize=10.0,
                                                                                     markeredgecolor="black",
                                                                                     label=requestY[j], picker=True)})
                                    print_("done plotting xlsx 2")
                                # mplcursors.cursor(args["plot object"][id_][j])
                                self.ax.set_ylabel('$Y$ []')
                                self.ax.set_xlabel('$X$ []')
                            else:
                                if type_ == 'Line':
                                    args["plot object"][id_].update({j: self.ax_right.plot(x_data, y, label=requestY[j],
                                                                                           linewidth=2,
                                                                                           linestyle=style,
                                                                                           picker=True)})
                                else:
                                    args["plot object"][id_].update({j: self.ax_right.plot(x_data, y, linestyle='None',
                                                                                           marker=style,
                                                                                           markeredgecolor="black",
                                                                                           markersize=10.0,
                                                                                           label="Legend",
                                                                                           picker=True)})
                                # mplcursors.cursor(args["plot object"][id_][j])
                                self.ax_right.set_ylabel('$Y$ [dB]')
                else:
                    print_("txt", filename)
                    # try to filter self.other_data
                    try:
                        print_("1a")
                        filter_ = args["plot inputs"]['Filter'][0].currentText()
                        value = args["plot inputs"]['Filter'][1].text()
                        if filter_ == "None" or value == "":
                            print_("4a")
                            self.other_data_filtered = self.other_data[key]
                            print_("5a")
                        else:
                            print_("2a")
                            self.other_data_filtered = self.other_data[key][self.other_data[key][filter_] == value]
                            print_("3a")

                    except:
                        print_("6a")
                        self.other_data_filtered = self.other_data[key]
                        print_("7a")

                    print_("8a")
                    x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                    self.freq_glob = x_data
                    print_("9a")

                    for j in range(len(requestY)):
                        print_("10a")
                        y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                        print_("11a")
                        args["plot data"][id_].update({j: {"x": x_data, "y": y}})

                        print_("txt: ", x_data, y)
                        if axis == 'Left':
                            if type_ == 'Line':
                                args["plot object"][id_].update(
                                    {j: self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)})
                                print_("Done plotting txt")
                            else:
                                args["plot object"][id_].update(
                                    {j: self.ax.plot(x_data, y, linestyle='None', marker=style, markersize=10.0,
                                                     markeredgecolor="black",
                                                     label="Legend", picker=True)})
                            # mplcursors.cursor(args["plot object"][id_][j])
                            self.ax.set_ylabel('$Y$ []')
                            self.ax.set_xlabel('$X$ []')
                        else:
                            if type_ == 'Line':
                                args["plot object"][id_].update(
                                    {j: self.ax_right.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)})
                            else:
                                args["plot object"][id_].update(
                                    {j: self.ax_right.plot(x_data, y, linestyle='None', marker=style, markersize=10.0,
                                                           markeredgecolor="black",
                                                           label="Legend", picker=True)})
                            # mplcursors.cursor(args["plot object"][id_][j])
                            self.ax_right.set_ylabel('$Y$ [dB]')
                print_("#" * 50)
                pass
            else:
                print_("Please specify columns to plot")

    def table_control(self):
        # fill the first line
        self.create_new_row(0, self.ui.tableWidget)

    def create_new_row(self, row_ind, table):

        # Toggle on/off
        w_Toggle_Close = QWidget()
        l_Toggle_Close_Widget = QHBoxLayout()
        w_Toggle_Close.setLayout(l_Toggle_Close_Widget)
        table.setCellWidget(row_ind, 0, w_Toggle_Close)

        cb_toggle = QCheckBox()
        cb_toggle.setCheckState(Qt.Checked)
        cb_toggle.stateChanged.connect(lambda: self.plot())

        pb_Delete_Row = QPushButton()
        pb_Delete_Row.setMinimumSize(QtCore.QSize(20, 20))
        pb_Delete_Row.setMaximumSize(QtCore.QSize(20, 20))
        pb_Delete_Row.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/icons/icons/PNG/stop.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        pb_Delete_Row.setIcon(icon3)
        pb_Delete_Row.setIconSize(QtCore.QSize(20, 20))
        pb_Delete_Row.clicked.connect(lambda: self.remove_row(pb_Delete_Row))

        l_Toggle_Close_Widget.setSpacing(0)
        l_Toggle_Close_Widget.setContentsMargins(0, 0, 0, 0)
        l_Toggle_Close_Widget.addWidget(cb_toggle)
        l_Toggle_Close_Widget.addWidget(pb_Delete_Row)

        # Code
        cb_code = QComboBox()
        cb_code.addItem("ABCI")
        cb_code.addItem("SLANS")
        cb_code.addItem("Other")
        table.setCellWidget(row_ind, 1, cb_code)

        # Folder widget
        w_Folder = QWidget()
        l_Folder_Widget = QHBoxLayout()
        w_Folder.setLayout(l_Folder_Widget)
        le_Folder = QLineEdit()
        le_Folder.setReadOnly(True)
        pb_Open_Folder = QPushButton('...')
        pb_Open_Folder.setMaximumWidth(100)
        pb_Open_Folder.setMinimumWidth(50)
        # signal to change place holder text is 'Other'
        cb_code.currentIndexChanged.connect(lambda: le_Folder.setPlaceholderText(
            'Select file') if cb_code.currentText() == 'Other' else le_Folder.setPlaceholderText('Select Folder'))

        l_Folder_Widget.setSpacing(0)
        l_Folder_Widget.setContentsMargins(0, 0, 0, 0)
        l_Folder_Widget.addWidget(le_Folder)
        l_Folder_Widget.addWidget(pb_Open_Folder)
        table.setCellWidget(row_ind, 2, w_Folder)

        # Polarisation
        cb_pol = QComboBox()
        cb_pol.addItem("Long")
        cb_pol.addItem("Trans")
        table.setCellWidget(row_ind, 3, cb_pol)
        # signal to disable polarisation when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(
            lambda: cb_pol.setDisabled(True) if cb_code.currentText() == 'Other' else cb_pol.setEnabled(True))

        # Id
        ccb_Id = QCheckableComboBox()
        ccb_Id.addItem("All")
        table.setCellWidget(row_ind, 4, ccb_Id)
        # signal to disable id when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(
            lambda: ccb_Id.setDisabled(True) if cb_code.currentText() == 'Other' else ccb_Id.setEnabled(True))

        # push button signal
        pb_Open_Folder.clicked.connect(
            lambda: self.open_(le_Folder, cb_pol, ccb_Id, 'File',
                               start_dir=f"{self.main_control.projectDir}") if cb_code.currentText() == 'Other'
        else self.open_(le_Folder, cb_pol, ccb_Id, 'Folder',
                        start_dir=f"{self.main_control.projectDir}/SimulationData/{cb_code.currentText()}"))
        # le_Folder.textChanged.connect(lambda: )

        cb_pol.currentIndexChanged.connect(
            lambda: self.populate_IDs(le_Folder.text(), cb_pol,
                                      ccb_Id) if cb_code.currentText() == 'Other' else self.populate_IDs(
                le_Folder.text(), cb_pol, ccb_Id))

        # Request
        # Request widget
        w_Request = QWidget()
        l_Request_Widget = QHBoxLayout()
        l_Request_Widget.setSpacing(0)
        l_Request_Widget.setContentsMargins(0, 0, 0, 0)
        w_Request.setLayout(l_Request_Widget)

        request_dict_long = {0: 'Longitudinal Impedance Magnitude',
                             1: r'Cavity Shape Input',
                             2: r'Cavity Shape Used',
                             3: r'Wake Potentials',
                             4: r'Real Part of Longitudinal Impedance',
                             5: r'Imaginary Part of Longitudinal Impedance',
                             6: r'Frequency Spectrum of Loss Factor',
                             7: r'Loss Factor Spectrum Integrated up to F',
                             8: r'Real Part of Long. + Log Impedance',
                             9: r'Imaginary Part of Long. + Log Impedance',
                             10: r'Spectrum of Long. + Log Loss Factor',
                             11: r'Long. + Log Factor Integrated up to F'
                             }

        request_dict_trans = {0: 'Transverse Impedance Magnitude',
                              1: r'Cavity Shape Input',
                              2: r'Cavity Shape Used',
                              3: r'Wake Potentials',
                              4: r'Real Part of Azimuthal Impedance',
                              5: r'Imaginary Part of Azimuthal Impedance',
                              6: r'Real Part of Transverse Impedance',
                              7: r'Imaginary Part of Transverse Impedance',
                              8: r'Real Part of Longitudinal Impedance',
                              9: r'Imaginary Part of Longitudinal Impedance'
                              }

        cb_request = QComboBox()
        for req in request_dict_long.values():
            cb_request.addItem(req)

        # create widget with two comboboxes
        cb_X = QComboBox()
        cb_Y = QCheckableComboBox()

        # place cb_request and le_request in request widget layout
        l_Request_Widget.addWidget(cb_request)
        l_Request_Widget.addWidget(cb_X)
        l_Request_Widget.addWidget(cb_Y)

        # hide le_request
        cb_X.hide()
        cb_Y.hide()

        # add widget to table
        table.setCellWidget(row_ind, 5, w_Request)

        # Filter
        # Filter widget
        w_Filter = QWidget()
        l_Filter_Widget = QHBoxLayout()
        l_Filter_Widget.setSpacing(0)
        l_Filter_Widget.setContentsMargins(0, 0, 0, 0)
        w_Filter.setLayout(l_Filter_Widget)

        # checkable combo box
        ccb_Filter = QCheckableComboBox()
        ccb_Filter.addItem("None")
        ccb_Filter.setMinimumWidth(75)
        l_Filter_Widget.addWidget(ccb_Filter)

        # line edit
        le_Value = QLineEdit()
        l_Filter_Widget.addWidget(le_Value)

        table.setCellWidget(row_ind, 10, w_Filter)

        # add signal to le_Folder after editing to update cb_X and cb_Y
        key = self.plotID_count
        le_Folder.textChanged.connect(lambda: self.populate_combobox([cb_X, cb_Y, ccb_Filter], le_Folder.text(), key))

        # signal to disable polarisation when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(lambda: self.show_hide_request_widgets(cb_code, cb_request, cb_X, cb_Y))

        # add signal to switch between longitudinal and transverse
        cb_pol.currentIndexChanged.connect(
            lambda: self.switchRequest(cb_request, request_dict_long, request_dict_trans, cb_pol))

        # ScaleX
        le_ScaleX = QLineEdit()
        le_ScaleX.setText('1')
        table.setCellWidget(row_ind, 6, le_ScaleX)

        # ScaleY
        le_ScaleY = QLineEdit()
        le_ScaleY.setText('1')
        table.setCellWidget(row_ind, 7, le_ScaleY)
        le_ScaleX.textChanged.connect(lambda: self.validating(le_ScaleX))
        le_ScaleY.textChanged.connect(lambda: self.validating(le_ScaleY))

        # axis
        axis_list = ['Left', 'Right', 'Top', 'Bottom']
        cb_axis = QComboBox()
        for a in axis_list:
            cb_axis.addItem(a)

        table.setCellWidget(row_ind, 8, cb_axis)

        # type
        w_Type = QWidget()
        l_Type_Widget = QHBoxLayout()
        l_Type_Widget.setSpacing(0)
        l_Type_Widget.setContentsMargins(0, 0, 0, 0)
        w_Type.setLayout(l_Type_Widget)
        type_list = ['Line', 'Scatter', 'Bar']

        # check box
        cb_type = QComboBox()
        for a in type_list:
            cb_type.addItem(a)

        l_Type_Widget.addWidget(cb_type)
        # check box
        cb_style = QComboBox()
        l_Type_Widget.addWidget(cb_style)

        line_style = ['-', '--', ':', '-.']
        scatter_marker = ['x', 'o', '+', 'P', 'X', 'd', '>', '^', 'v', 'H']
        # fill default
        for a in line_style:
            cb_style.addItem(a)

        # signal for combobox
        cb_type.currentIndexChanged.connect(
            lambda: self.populate_combobox_list(cb_style, line_style)
            if cb_type.currentText() == "Line"
            else self.populate_combobox_list(cb_style, scatter_marker))

        table.setCellWidget(row_ind, 9, w_Type)

        args_dict = {'Code': cb_code, 'Folder': [le_Folder, pb_Open_Folder, w_Folder, l_Folder_Widget],
                     'Polarization': cb_pol, 'Id': ccb_Id,
                     'Request': [cb_request, cb_X, cb_Y, w_Request, l_Request_Widget], 'Toggle': cb_toggle,
                     'Remove': pb_Delete_Row, 'ScaleX': le_ScaleX, 'ScaleY': le_ScaleY, 'Axis': cb_axis,
                     'Type': [cb_type, cb_style], 'Filter': [ccb_Filter, le_Value]}

        self.plot_dict[key] = {"plot inputs": None, "plot data": {}, "plot object": {},
                               "plot data inputs": None}

        self.plot_dict[key]["plot inputs"] = args_dict
        self.plotID_count += 1

    def generateMenu(self, pos):
        # Get index
        menu = QMenu()
        item1 = menu.addAction("Add Row")
        item2 = menu.addAction("Delete Row")
        # Make the menu display in the normal position
        screenPos = self.ui.tableWidget.mapToGlobal(pos)

        # Click on a menu item to return, making it blocked
        action = menu.exec(screenPos)
        if action == item1:
            self.add_row()
        if action == item2:
            self.remove_row(self.row)
        else:
            return

    def add_row(self):
        # get current number of rows
        n = self.ui.tableWidget.rowCount()

        # add new row
        self.ui.tableWidget.setRowCount(n + 1)  # and one row in the table
        self.create_new_row(n, self.ui.tableWidget)

    def remove_row(self, pb_Remove):

        # remove from table
        n = self.ui.tableWidget.rowCount()
        # print_("\t\t\t\t\trow", row)
        row = 0
        key = 0
        code = "ABCI"
        for i, (k, val) in enumerate(self.plot_dict.items()):
            if val["plot inputs"]['Remove'] == pb_Remove:
                print_("Matched remove button to row")
                row = i
                key = k
                code = val["plot inputs"]['Code'].currentText()

        self.ui.tableWidget.removeRow(row)
        # reset number of rows
        self.ui.tableWidget.setRowCount(n - 1)

        if code == "ABCI":
            for line2D in self.plot_dict[key]['plot object'].values():
                line2D[0].remove()
        else:
            for vals in self.plot_dict[key]['plot object'].values():
                print_(vals)
                for line2D in vals.values():
                    line2D[0].remove()

        # reset plot data, plot object and plot data inputs
        self.plot_dict[key]["plot data"] = {}
        self.plot_dict[key]["plot object"] = {}
        self.plot_dict[key]["plot data inputs"] = {}

        self.plot_dict[key] = {}
        del self.plot_dict[key]
        print_("Here after matched remove button to row")

        self.plot()

    @staticmethod
    def switchRequest(cb_Request, request_dict_long, request_dict_trans, cb_Code):
        # clear combo box
        cb_Request.clear()
        print_('here')
        # add new items
        if cb_Code.currentText().lower() == 'long':
            for req in request_dict_long.values():
                cb_Request.addItem(req)
        else:
            for req in request_dict_trans.values():
                cb_Request.addItem(req)

    @staticmethod
    def switchCode(cb_Code):
        if cb_Code.currentText() == 'ABCI':
            pass
        elif cb_Code.currentText() == 'SLANS':
            pass

    def mousePressEvent(self, event):
        if event.buttons() & Qt.RightButton:
            row = self.ui.tableWidget.currentRow()
            self.row = row

    def calc_limits(self, mode):
        if len(self.freq_glob) > 0:
            if self.ui.cb_Longitudinal_Threshold.checkState() == 2 or self.ui.cb_Transverse_Threshold.checkState() == 2:
                E0 = [45.6, 80, 120, 182.5]  # [GeV] Energy
                nu_s = [0.025, 0.0506, 0.036, 0.087]  # Synchrotron oscillation tune
                I0 = [1280, 135, 29, 10.8]  # [mA] Beam current5.4 * 2
                alpha_c = [1.48, 1.48, 0.73, 0.73]  # [10−5] Momentum compaction factor
                tau_z = [424.6, 78.7, 23.4, 6.8]  # [ms] Longitudinal damping time
                tau_xy = [849.2, 157.4, 46.8, 13.6]  # [ms] Transverse damping time
                frev = [3.07, 3.07, 3.07, 3.07]  # [kHz] Revolution frequency
                beta_xy = 50

                #     Ncav = [52, 52, 136, 584] # Number of cavities per beam
                Ncav = [52, 100, 268, 584]  # 1_2_2_25

                E0 = ast.literal_eval(self.ui.le_E0.text())
                nu_s = ast.literal_eval(self.ui.le_Nu_S.text())
                I0 = ast.literal_eval(self.ui.le_I0.text())
                alpha_c = ast.literal_eval(self.ui.le_Alpha_S.text())
                tau_z = ast.literal_eval(self.ui.le_Tau_Z.text())
                tau_xy = ast.literal_eval(self.ui.le_Tau_XY.text())
                frev = ast.literal_eval(self.ui.le_F_Rev.text())
                beta_xy = ast.literal_eval(self.ui.le_Beta_XY.text())
                Ncav = ast.literal_eval(self.ui.le_N_Cav.text())

                Z_list = []
                if mode == 'monopole':
                    # trim f
                    # f_list = f_list[0:len(f_list) - 100]
                    # f_list = self.freq_glob[0:len(self.freq_glob)]
                    f_list = np.linspace(0, self.freq_glob[-1], num=1000)
                    for i, n in enumerate(Ncav):
                        Z = [(2 * E0[i] * nu_s[i]) * 1e8 / (n * I0[i] * alpha_c[i] * tau_z[i] * a) if a > 1e-8 else 1e5
                             for a in
                             f_list]
                        Z_list.append(Z)

                    self.plot_baselines(f_list, Z_list, mode)

                    return f_list, Z_list

                elif mode == 'dipole':
                    # f_list = self.freq_glob[0:len(self.freq_glob)]
                    f_list = np.linspace(0, self.freq_glob[-1], num=1000)
                    for i, n in enumerate(Ncav):
                        Z = (2 * E0[i]) * 1e9 / (n * I0[i] * beta_xy * tau_xy[i] * frev[i])
                        Z_list.append(Z)

                    self.plot_baselines(f_list, Z_list, mode)

                    return Z_list

            else:
                self.remove_baselines()

    def plot_baselines(self, f_list, Z_list, mode):
        if mode == 'monopole':
            # plot baselines
            text = ['Z', 'W', 'H', 'tt']
            for i, z in enumerate(Z_list):
                aa = self.ax.plot(f_list, z, ls='--', c='gray')
                # ab = self.ax.text(1, z[200], f'{text[i]}')

                # keep record
                self.baseline_line_objects.append(aa[0])
                # self.baseline_line_objects.append(ab)

            self.ax.autoscale(True, axis='y')
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        else:
            # plot baselines
            text = ['Z', 'W', 'H', 'tt']
            for i, z in enumerate(Z_list):
                aa = self.ax.axhline(z, ls='--', c='gray')
                # ab = self.ax.text(0.5, z+1.8, f'{text[i]}', transform = self.ax.transAxes, va='top')
                # coordinate system transformation

                # keep record
                self.baseline_line_objects.append(aa)
                # self.baseline_line_objects.append(ab)

            self.ax.autoscale(True, axis='y')
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

    def remove_baselines(self):
        for line in self.baseline_line_objects:
            if line in self.ax.get_lines():
                line.remove()

            if line in self.ax.findobj():
                line.remove()

        self.ax.autoscale(True, axis='y')
        self.ax.relim()
        self.plt.draw()

    def populate_cutoff_combobox(self, cutoff):
        cutoff.clear()
        cutoff.addItem("All")

        ll = text_to_list(self.ui.le_Ri_Cutoff_List.text())
        if ll is None:
            Ri_list = []
        else:
            Ri_list = text_to_list_(self.ui.le_Ri_Cutoff_List.text())

        c = 299792458

        M = [0, 1, 2, 3]
        N = [1, 2, 3, 4]
        mode_list = {}
        mode_name_list = ['TM', 'TE']
        f_list = {}

        for Ri in Ri_list:
            mode_list[f'{Ri}'] = []
            f_list[f'{Ri}'] = []
            for m in M:
                for n in N:
                    # get jacobian
                    j_mn = jn_zeros(m, n)[n - 1]
                    j_mn_p = jnp_zeros(m, n)[n - 1]
                    J = [j_mn, j_mn_p]
                    # for p in P:
                    for mode_type, j in enumerate(J):
                        # formula
                        # f = c / (2 * np.pi) * ((j / Ri) ** 2 + (p * np.pi / L) ** 2) ** 0.5
                        f = c / (2 * np.pi) * (j / (Ri * 1e-3))

                        # append mode name
                        mode_list[f'{Ri}'].append(f'{mode_name_list[mode_type]}{m}{n}')

                        # append to f list
                        f_list[f'{Ri}'].append(f * 1e-6)

            self.f_list_sorted[f'{Ri}'] = sorted(f_list[f'{Ri}'])
            self.mode_list_sorted[f'{Ri}'] = [x for _, x in sorted(zip(f_list[f'{Ri}'], mode_list[f'{Ri}']))]

        print(self.mode_list_sorted)
        print(self.f_list_sorted)

        if len(Ri_list) != 0:
            cutoff.addItems(self.mode_list_sorted[f"{Ri_list[0]}"])

        color = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']
        step = [0.02, 0.2, 0.4, 0.6, 0.8]
        count = 0

    def plot_cutoff(self, cutoff):
        # selected
        selected = cutoff.currentText().split(", ")
        Ri_list = text_to_list(self.ui.le_Ri_Cutoff_List.text())

        selected_combination = []
        selected_combination_keys = []
        for s in selected:
            if s != '' and s != "All":
                for Ri in Ri_list:
                    selected_combination.append([s, Ri])
                    selected_combination_keys.append(f'{[s, Ri]}')

        for sc in selected_combination:
            if f"{sc}" not in self.ax_obj_dict.keys():
                indx = self.mode_list_sorted[f'{sc[1]}'].index(f"{sc[0]}")
                freq = self.f_list_sorted[f'{sc[1]}'][indx]

                vl = self.ax.axvline(freq, label=f"{sc[0]} cutoff (Ri={sc[1]})", ls='--', c='k')

                ab = self.plt.add_text(r"$f_\mathrm{c," + f"{sc[0]}" + r"} (R_\mathrm{i} = "
                                       + f"{sc[1]}" + r" ~\mathrm{mm}) $",
                                       box="None", xy=(freq, 0.02),
                                       xycoords='data', size=14, rotation=90)

                # update axes object dictionary
                self.ax_obj_dict.update({f"{sc}": [vl, ab]})

        # compare selected to ax_obj_dict
        deleted_lines = []
        for key in self.ax_obj_dict.keys():
            if key not in selected_combination_keys:
                deleted_lines.append(key)

        for k in deleted_lines:
            for obj in self.ax_obj_dict[k]:
                if obj in self.ax.findobj():
                    obj.remove()
                else:
                    # remove object from text dictionary in plotter
                    self.plt.text_dict[f"{id(obj)}"].remove()

            del self.ax_obj_dict[f"{k}"]
        self.fig.canvas.draw()

    @staticmethod
    def calc_cutoff(Ri, mode):
        # calculate frequency from Ri
        p_TM01, p_TE11 = 2.405, 1.841
        c = 299792458  # m/s

        if mode == 'TM01':
            freq = (c * p_TM01) / (2 * np.pi * Ri * 1e9) * 1e3
        else:
            print_("here")
            freq = (c * p_TE11) / (2 * np.pi * Ri * 1e9) * 1e3

        return freq

    def toggle_page(self, key):

        if key == 'Plot Decorations':
            self.ui.sw_Plot_Area_Tools.setCurrentIndex(0)

        if key == 'Machine Parameters':
            self.ui.sw_Plot_Area_Tools.setCurrentIndex(1)

    def toggle_axis_labels(self):
        if len(self.ax.get_lines()) == 0:
            # turn on axis ticks and labels
            self.ax.axes.xaxis.set_visible(False)
            self.ax.axes.yaxis.set_visible(False)
        else:
            # turn on axis ticks and labels
            self.ax.axes.xaxis.set_visible(True)
            self.ax.axes.yaxis.set_visible(True)

        if len(self.ax_right.get_lines()) == 0:
            # turn on axis ticks and labels
            self.ax_right.axes.yaxis.set_visible(False)
            # set colors
            self.ax_right.spines['right'].set_color('#000000')
            self.ax_right.tick_params(axis='y', colors='#000000')
        else:
            colors = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']
            # turn on axis ticks and labels
            self.ax_right.axes.yaxis.set_visible(True)

            # set colors
            self.ax_right.spines['right'].set_color(colors[0])
            self.ax_right.tick_params(axis='y', colors=colors[0])
            self.ax_right.yaxis.label.set_color(colors[0])

    def open_(self, le, cb_pol=None, ccb=None, mode='Folder', row_ind=0, start_dir=''):
        if mode == "Folder":
            data_dir = str(
                QFileDialog.getExistingDirectory(None, "Select Directory", start_dir,
                                                 QFileDialog.ShowDirsOnly))
            if data_dir != '':
                le.setText(data_dir)

                # populate checkable combobox
                self.populate_IDs(data_dir, cb_pol, ccb)
        else:
            filename, _ = QFileDialog.getOpenFileName(None, "Open File", start_dir,
                                                      "Excel Files (*.txt *.xlsx *.json)")
            if filename != '':
                le.setText(filename)

    @staticmethod
    def populate_IDs(dirc, cb_pol, ccb):
        if dirc != "":
            dir_list = os.listdir(dirc)

            # copy selection
            selection = []
            if ccb.currentText():
                selection = ccb.currentText().split(", ")

            # clear checkable check box
            ccb.clear()
            ccb.addItem("All")

            # polarisation
            pol = cb_pol.currentText()
            if pol == "Long":
                i = 0
            else:
                i = 1

            # loop through folder list to check which contained simulation files
            for d in dir_list:
                if os.path.exists(fr"{dirc}\{d}\Cavity_MROT_{i}.pot"):
                    ccb.addItem(d)

            # try to select copied selection
            for txt in selection:
                if txt in dir_list:
                    item = ccb.model().item(dir_list.index(txt) + 1)  # +1 because 'all' takes position 0
                    item.setCheckState(2)

    # def open_file(self):
    #     filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Excel Files (*.xlsx)")
    #
    #     try:
    #         self.ui.le_Filename.setText(filename)
    #         df = pd.read_excel(filename)
    #
    #         # get header
    #         header = df.columns
    #
    #         # create combobox for x checkbox list for y
    #         cb_X_Axis = QComboBox()
    #         cb_Y_Axis = QCheckableComboBox()
    #         cb_Y_Axis_List = []
    #
    #         # populate x and y with checkboxes
    #         for col in header:
    #             cb_X_Axis.addItem(col)
    #             cb_Y_Axis.addItem(col)
    #
    #             cb = QCheckBox(col)
    #             cb_Y_Axis_List.append(cb)
    #
    #         # self.ui.vl_Y_Axis_Checkbox.addWidget(cb_Y_Axis)
    #         #
    #         # # add to layout
    #         # self.ui.vl_X_Axis_Combobox.addWidget(cb_X_Axis)
    #         # self.ui.vl_X_Axis_Combobox.addWidget(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # except Exception as e:
        #     self.log.error(f'Failed to plot:: {e}')

    def set_table_size(self):
        self.ui.tableWidget.setColumnWidth(0, 50)
        self.ui.tableWidget.setColumnWidth(1, 75)
        self.ui.tableWidget.setColumnWidth(2, 200)
        self.ui.tableWidget.setColumnWidth(3, 100)
        self.ui.tableWidget.setColumnWidth(4, 300)
        self.ui.tableWidget.setColumnWidth(5, 200)
        self.ui.tableWidget.setColumnWidth(6, 75)
        self.ui.tableWidget.setColumnWidth(7, 75)
        self.ui.tableWidget.setColumnWidth(8, 75)
        self.ui.tableWidget.setColumnWidth(10, 200)

    @staticmethod
    def validating(le):
        validation_rule = QDoubleValidator(1, 3840, 0)
        if validation_rule.validate(le.text(), 20)[0] == QValidator.Acceptable:
            le.setFocus()
        else:
            le.setText('')

    @staticmethod
    def show_hide_request_widgets(cb_code, cb_request, cb_X, cb_Y):
        if cb_code.currentText() == 'Other':
            cb_request.hide()
            cb_X.show()
            cb_Y.show()
        else:
            cb_request.show()
            cb_X.hide()
            cb_Y.hide()

    def populate_combobox_new(self, cb_list, filename, pid):
        # enable ID combobox
        # PTS.ccb_Id_Other, PTS.cb_X, PTS.ccb_Y, PTS.ccb_Filter_Other], PTS.le_Filename.text())

        if filename.split('.')[-1] == 'xlsx':

            # load file (for now, just excel files)
            df = fr.excel_reader(filename)
            sheet_name = list(df.keys())

            # populate combobox with sheet name
            for sh_name in sheet_name:
                cb_list[0].addItem(sh_name)

            dd = df[sheet_name[0]]
            self.pts_dict[pid]["plot data"] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()

                for a in dd.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1] == 'txt':
            # load file txt
            df = fr.txt_reader(filename, "\t")
            self.pts_dict[pid]["plot data"] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()

                for a in df.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1] == 'json':
            # load file txt
            df = fr.json_reader(filename)
            self.pts_dict[pid]["plot data"] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()
                for a in df.keys():
                    cb.addItem(f"{a}")

    def populate_combobox(self, cb_list, filename, i=0):
        # enable ID combobox
        self.plot_dict[i]["plot inputs"]['Id'].setEnabled(True)
        self.plot_dict[i]["plot inputs"]['Id'].clear()
        self.plot_dict[i]["plot inputs"]['Id'].addItem("All")

        if filename.split('.')[-1] == 'xlsx':

            # load file (for now, just excel files)
            df = fr.excel_reader(filename)
            sheet_name = list(df.keys())

            # populate combobox with sheet name
            for sh_name in sheet_name:
                self.plot_dict[i]["plot inputs"]['Id'].addItem(sh_name)

            dd = df[sheet_name[0]]
            self.other_data[i] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()

                for a in dd.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1] == 'txt':
            # load file txt
            df = fr.txt_reader(filename, "\t")
            self.other_data[i] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()

                for a in df.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1] == 'json':
            # load file txt
            df = fr.json_reader(filename)
            self.other_data[i] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()
                for a in df.keys():
                    cb.addItem(f"{a}")

    @staticmethod
    def populate_combobox_list(cb, ll):
        cb.clear()
        for a in ll:
            cb.addItem(a)

    def update_labels(self):
        # select axes to update
        print_("update labesl", self.ui.cb_Active_Axis.currentText())
        if self.ui.cb_Active_Axis.currentText() == 'Left':
            ax_current = self.ax
            ax_other = self.ax_right
        else:
            ax_current = self.ax_right
            ax_other = self.ax

        xlabel = self.ui.le_Xlabel.text()
        ylabel = self.ui.le_Ylabel.text()
        title = self.ui.le_Title.text()
        xsize = self.ui.sb_XLabel_Size.value()
        ysize = self.ui.sb_YLabel_Size.value()
        title_size = self.ui.sb_Title_Size.value()
        xtick_size = self.ui.sb_XLabel_Tick_Size.value()
        ytick_size = self.ui.sb_YLabel_Tick_Size.value()
        legend_size = self.ui.sb_Legend_Size.value()

        legend_labels = self.ui.le_Legend.text().split("%%")

        # update plot
        ax_current.set_xlabel(xlabel, fontsize=xsize)
        ax_current.set_ylabel(ylabel, fontsize=ysize)
        ax_current.set_title(title, fontsize=title_size)
        ax_current.tick_params(axis='x', labelsize=xtick_size, size=xtick_size)
        ax_current.tick_params(axis='y', labelsize=ytick_size, size=ytick_size)

        # update axes decoration dict
        self.axes_decorations_dict["xlabel"] = xlabel
        self.axes_decorations_dict["ylabel"] = ylabel
        self.axes_decorations_dict["title"] = title

        # update legend
        lines, labels = ax_current.get_legend_handles_labels()
        lines2, labels2 = ax_other.get_legend_handles_labels()
        # handles, labels = ax_other.get_legend_handles_labels()
        # self.leg = self.ax.legend(lines + lines2, labels + labels2, loc='lower left', prop={'size': 18})

        labels = labels + labels2
        try:
            for i in range(len(legend_labels)):
                labels[i] = legend_labels[i] if legend_labels[i] != "" else labels[i]
        except IndexError:
            pass

        legend_other = ax_other.get_legend()
        if legend_other:
            # remove old legend
            legend_other.remove()

        self.leg = ax_current.legend(lines + lines2, labels, fontsize=legend_size)

        self.leg.set_zorder(10)
        self.leg.set_draggable(state=True, use_blit=True)

        self.fig.canvas.draw_idle()

    def serialize(self, state_dict):
        # update state file
        state_dict["Inlet_Position"] = self.ui.le_Inset_Position.text()
        state_dict["Inset_Window"] = self.ui.le_Inset_Window.text()
        state_dict["Show_Hide_Inset"] = self.ui.cb_Show_Hide_Inset.checkState()
        state_dict["Y_Scale"] = self.ui.cb_Y_Scale.currentText()

        state_dict["Show_Wakefield_Parameters"] = self.ui.cb_Show_Wakefield_Parameters.checkState()

        state_dict["Text"] = self.ui.le_Plot_Text.text()
        state_dict["Peaks Threshold"] = self.ui.dsb_Threshold.value()

        state_dict["E0"] = self.ui.le_E0.text()
        state_dict["Nu_S"] = self.ui.le_Nu_S.text()
        state_dict["I0"] = self.ui.le_I0.text()
        state_dict["Alpha_S"] = self.ui.le_Alpha_S.text()
        state_dict["Tau_Z"] = self.ui.le_Tau_Z.text()
        state_dict["Tau_XY"] = self.ui.le_Tau_XY.text()
        state_dict["F_Rev"] = self.ui.le_F_Rev.text()
        state_dict["Beta_XY"] = self.ui.le_Beta_XY.text()
        state_dict["N_Cav"] = self.ui.le_N_Cav.text()
        state_dict["Longitudinal Threshold Checkbox"] = self.ui.cb_Longitudinal_Threshold.checkState()
        state_dict["Transverse Threshold Checkbox"] = self.ui.cb_Transverse_Threshold.checkState()
        state_dict["Threshold Line Color"] = self.ui.cb_Threshold_Line_Color.currentText()
        state_dict["Threshold Linestyle"] = self.ui.cb_Threshold_Linestyle.currentText()

        state_dict["xlabel"] = self.ui.le_Xlabel.text()
        state_dict["ylabel"] = self.ui.le_Ylabel.text()
        state_dict["title"] = self.ui.le_Title.text()
        state_dict["xlabel_size"] = self.ui.sb_XLabel_Size.value()
        state_dict["ylabel_size"] = self.ui.sb_YLabel_Size.value()
        state_dict["title_size"] = self.ui.sb_Title_Size.value()
        state_dict["xlabeltick_size"] = self.ui.sb_XLabel_Tick_Size.value()
        state_dict["ylabeltick_size"] = self.ui.sb_YLabel_Tick_Size.value()
        state_dict["legend_size"] = self.ui.sb_Legend_Size.value()

    def deserialize(self, state_dict):
        self.ui.le_Inset_Position.setText(state_dict["Inlet_Position"])
        self.ui.le_Inset_Window.setText(state_dict["Inset_Window"])

        self.ui.cb_Show_Hide_Inset.setCheckState(state_dict["Show_Hide_Inset"])

        self.ui.cb_Y_Scale.setCurrentText(state_dict["Y_Scale"])

        self.ui.cb_Show_Wakefield_Parameters.setCheckState(state_dict["Show_Wakefield_Parameters"])

        self.ui.le_Plot_Text.setText(state_dict["Text"])
        self.ui.dsb_Threshold.setValue(state_dict["Peaks Threshold"])

        self.ui.le_E0.setText(state_dict["E0"])
        self.ui.le_Nu_S.setText(state_dict["Nu_S"])
        self.ui.le_I0.setText(state_dict["I0"])
        self.ui.le_Alpha_S.setText(state_dict["Alpha_S"])
        self.ui.le_Tau_Z.setText(state_dict["Tau_Z"])
        self.ui.le_Tau_XY.setText(state_dict["Tau_XY"])
        self.ui.le_F_Rev.setText(state_dict["F_Rev"])
        self.ui.le_Beta_XY.setText(state_dict["Beta_XY"])
        self.ui.le_N_Cav.setText(state_dict["N_Cav"])
        self.ui.cb_Longitudinal_Threshold.setCheckState(state_dict["Longitudinal Threshold Checkbox"])
        self.ui.cb_Transverse_Threshold.setCheckState(state_dict["Transverse Threshold Checkbox"])
        self.ui.cb_Threshold_Line_Color.setCurrentText(state_dict["Threshold Line Color"])
        self.ui.cb_Threshold_Linestyle.setCurrentText(state_dict["Threshold Linestyle"])

        self.ui.le_Xlabel.setText(state_dict["xlabel"])
        self.ui.le_Ylabel.setText(state_dict["ylabel"])
        self.ui.le_Title.setText(state_dict["title"])
        self.ui.sb_XLabel_Size.setValue(state_dict["xlabel_size"])
        self.ui.sb_YLabel_Size.setValue(state_dict["ylabel_size"])
        self.ui.sb_Title_Size.setValue(state_dict["title_size"])
        self.ui.sb_XLabel_Tick_Size.setValue(state_dict["xlabeltick_size"])
        self.ui.sb_YLabel_Tick_Size.setValue(state_dict["ylabeltick_size"])
        self.ui.sb_Legend_Size.setValue(state_dict["legend_size"])
