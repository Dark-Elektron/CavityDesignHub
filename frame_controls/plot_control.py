import ast
import json
import os
import time

import numpy as np
# from PyQt5 import Qt
from PyQt5.QtGui import QStandardItemModel, QPalette, QFontMetrics, QStandardItem
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from ui_files.plot import Ui_Plot
from modules.plot_module.plotter import Plot
from modules.data_module.slans_data import SLANSData
from modules.data_module.abci_data import ABCIData
import matplotlib as mpl
import pandas as pd
from utils.file_reader import FileReader

fr = FileReader()

class PlotControl:
    def __init__(self, parent):
        self.w_Plot = QWidget()

        self.plotUI = Ui_Plot()
        self.plotUI.setupUi(self.w_Plot)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # get logger
        self.log = self.main_control.log

        # Plot
        self.plt = Plot(self.plotUI)
        self.plotUI.gl_Plot_Area.addWidget(self.plt)
        self.fig = self.plt.fig
        self.ax = self.plt.ax
        self.ax_right = self.plt.ax_right
        self.axins = None
        self.indicate_inset = None

        # global freq axis
        self.freq_glob = []

        self.initUI()
        self.signals()

    def signals(self):
        # plot abci impedance
        self.plotUI.pb_Plot.clicked.connect(lambda: self.plot())

        # refresh plot
        self.plotUI.pb_Refresh.clicked.connect(lambda: self.plot())

        # toggle plot menu
        self.plotUI.pb_Plot_Area_Menu.clicked.connect(lambda: self.main_control.animate_width(self.plotUI.w_Plot_Menu, 0, 440, True))
        self.plotUI.pb_Plot_Area_Menu.clicked.connect(lambda: self.plotUI.w_Plot_Area_Buttons.setEnabled(True) if self.plotUI.w_Plot_Menu.maximumWidth() == 0 else self.plotUI.w_Plot_Area_Buttons.setDisabled(True))

        # signal for plot menu pages
        self.plotUI.pb_Machine_Parameters.clicked.connect(lambda: self.toggle_page('Machine Parameters'))
        self.plotUI.pb_Plot_Decorations.clicked.connect(lambda: self.toggle_page('Plot Decorations'))

        # signal for plot argument entrmake widgets rey
        # self.plotUI.pb_Collapse_Shape_Parameters.clicked.connect(lambda: self.main_control.animate_height(self.plotUI.w_Plot_Args, 0, 200, True))

        # signal for threshold
        self.plotUI.cb_Longitudinal_Threshold.stateChanged.connect(lambda: self.calc_limits('monopole'))
        self.plotUI.cb_Transverse_Threshold.stateChanged.connect(lambda: self.calc_limits('dipole'))

        # signal for inset plot
        self.plotUI.cb_Show_Hide_Inset.clicked.connect(lambda: self.plot_inset())

        # clear plots
        self.plotUI.pb_Clear.clicked.connect(lambda: self.clear_plots())

    def initUI(self):
        # create argument dictionary
        self.args_dict = {'Code': [],
                          'Folder': [],
                          'Polarization': [],
                          'Id': [],
                          'Request': [],
                          'Toggle': [],
                          'ScaleX': [],
                          'ScaleY': [],
                          'Axis': [],
                          'Type': []
                          }

        # baseline matplotlib line.Line2D objects
        self.baseline_line_objects = []

        # add row to table
        self.plotUI.tableWidget.setRowCount(1)  # and one row in the table
        self.table_control()

        # set plot menu initial size to zero and disable plot buttons
        self.main_control.animate_width(self.plotUI.w_Plot_Menu, 0, 440, True)
        self.plotUI.w_Plot_Area_Buttons.setDisabled(True)

        # tableWidget initializations
        self.plotUI.tableWidget.mousePressEvent = self.mousePressEvent
        self.plotUI.tableWidget.setContextMenuPolicy(Qt.CustomContextMenu)
        self.plotUI.tableWidget.customContextMenuRequested.connect(self.generateMenu)
        self.set_table_size()

    def plot(self):
        try:
            # args = list(self.args_dict.values())
            # print(args)
            self.clear_plots()

            # use same color cycler for both axes
            self.ax_right._get_lines.prop_cycler = self.ax._get_lines.prop_cycler

            plot_count = 1
            for i in range(len(self.args_dict['Code'])):
                code = self.args_dict['Code'][i].currentText()
                # check plot type
                if code == 'ABCI':
                    self.plot_impedance(self.args_dict, i, plot_count)
                elif code == 'SLANS':
                    pass
                else:
                    self.plot_other(self.args_dict, i)

            # toggle axis labels
            self.toggle_axis_labels()

            self.ax.autoscale(True, axis='y')
            self.ax_right.autoscale(True, axis='y')

            # recomputer the ax.datalim
            self.ax.relim()
            self.ax_right.relim()

            # show legend
            lines, labels = self.ax.get_legend_handles_labels()
            lines2, labels2 = self.ax_right.get_legend_handles_labels()

            self.leg = self.ax_right.legend(lines + lines2, labels + labels2, loc='lower left', prop={'size': 18})
            self.leg.set_draggable(state=True, use_blit=True)

            # plot inset if check box is checked
            self.plot_inset()
            # self.fig.canvas.draw_idle()

        except:
            self.log.error("Please enter a valid argument")

    def clear_plots(self):
        self.ax.cla()
        self.ax_right.cla()

        if self.axins is not None:
            self.axins.cla()
            self.axins.remove()
            self.axins = None

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def plot_impedance(self, args, i, plot_count):
        ids = [a.strip() for a in args['Id'][i].currentText().split(',')] # get list

        pol = args['Polarization'][i].currentText()
        request = args['Request'][i][0].currentText()
        folder = args['Folder'][i][0].text()
        state = args['Toggle'][i].checkState()
        axis = args['Axis'][i].currentText()

        try:
            scaleX = float(args['ScaleX'][i].text())
            scaleY = float(args['ScaleY'][i].text())
        except:
            scaleX = 1
            scaleY = 1

        if folder == '':
            abci_data_dir = fr'{self.main_control.projectDir}/SimulationData/ABCI'
        else:
            abci_data_dir = folder

        # plot for multiple ids
        for id in ids:
            if pol.lower() == 'long':
                abci_data_long = ABCIData(abci_data_dir, id, 0)

                if request == 'Longitudinal Impedance Magnitude':
                    xr, yr, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                    xi, yi, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')

                    y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]

                    # scale x axis]
                    xr = [a*scaleX for a in xr]
                    # scale y axis]
                    y = [a*scaleY for a in y]

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_long.get_data(request)

                    # scale x axis]
                    xr = [a*scaleX for a in xr]
                    # scale y axis]
                    y = [a*scaleY for a in y]

                    # reset global frequency if the impedance magnitude is not requested
                    self.freq_glob = []

                # plot
                ax_selected = self.ax
                if axis == 'Left':
                    ax_selected = self.ax
                    # process id
                    self.ax.plot(xr, y, label='${' + f'{str(id).replace("_", ",")}' + '}$'+f' ({abci_data_long.wakelength} m)', linewidth=2)

                elif axis == "Right":
                    ax_selected = self.ax_right
                    self.ax_right.plot(xr, y, label='${' + f'{str(id).replace("_", ",")}' + '}$'+f' ({abci_data_long.wakelength} m)', linewidth=2)

                ax_selected.set_xlabel('$f \mathrm{ [GHz]}$')
                ax_selected.set_ylabel('$Z_{\parallel, \mathrm{HOM}} \mathrm{[k\Omega]}$')
                # ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

            else:
                abci_data_trans = ABCIData(abci_data_dir, id, 1)
                if request == 'Transverse Impedance Magnitude':
                    xr, yr, _ = abci_data_trans.get_data('Real Part of Azimuthal Impedance')
                    xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Azimuthal Impedance')
                    y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]

                    # scale x axis]
                    xr = [a*scaleX for a in xr]
                    # scale y axis]
                    y = [a*scaleY for a in y]

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_trans.get_data(request)

                    # scale x axis]
                    xr = [a*scaleX for a in xr]
                    # scale y axis]
                    y = [a*scaleY for a in y]

                    # reset global frequency if the impedance magnitude is not requested
                    self.freq_glob = []

                # plot
                ax_selected = self.ax
                if axis == 'Left':
                    ax_selected = self.ax
                    self.ax.plot(xr, y, label='${' + f'{str(id).replace("_", ",")}' + '}$'+f' ({abci_data_trans.wakelength} m)', linewidth=2)
                elif axis == 'Right':
                    ax_selected = self.ax_right
                    self.ax_right.plot(xr, y, label='${' + f'{str(id).replace("_", ",")}' + '}$'+f' ({abci_data_trans.wakelength} m)', linewidth=2)

                # plot settings
                ax_selected.set_xlabel('$f \mathrm{[GHz]}$')
                ax_selected.set_ylabel('$Z_{\perp, \mathrm{HOM}} \mathrm{[k\Omega/m]}$')

                ax_selected.set_yscale('log')
                ax_selected.set_ylim(1e-2, 1e2)
                ax_selected.set_xlim(0.25, 2.24)

            # increment plot count
            plot_count += 1

    def plot_inset(self):
        if self.plotUI.cb_Show_Hide_Inset.checkState() == 2:
            # get lines from axis
            lines = self.ax.get_lines()

            inset_pos = ast.literal_eval(self.plotUI.le_Inset_Position.text())
            if len(inset_pos) == 4:
                self.axins = self.ax.inset_axes(inset_pos, facecolor='#fafafa')
                self.axins.axes.xaxis.set_visible(False)
                self.axins.axes.yaxis.set_visible(False)

            for line in lines:
                if self.axins is not None:
                    if line.get_linestyle() == 'None':
                        self.axins.plot(line.get_xdata(), line.get_ydata(), linestyle='None', marker='x', markersize=10.0)
                    else:
                        self.axins.plot(line.get_xdata(), line.get_ydata(), linewidth=2)


            # sub region of the original image
            # get values from line edit
            try:
                x0, x1, y0, y1 = ast.literal_eval(self.plotUI.le_Inset_Window.text())
            except:
                x0, x1, y0, y1 = 0.385, 0.415, 5e-1, 3e1 # default

            self.axins.set_xlim(x0, x1)
            self.axins.set_ylim(y0, y1)

            self.axins.set_xticklabels('')
            self.axins.set_yticklabels('')

            self.axins.set_yscale(self.plotUI.cb_Y_Scale.currentText())
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

    def plot_other(self, args, i):
        # i is the code index
        filename = args['Id'][i].text()

        # # check if filename has correct extenstion
        # if filename.split('.')[-1] != 'xlsx':
        #     filename = fr'{filename}.xlsx'

        requestX = args['Request'][i][1].currentText()
        requestY = args['Request'][i][2].currentText().split(', ')
        print(requestX, requestY)

        filename = args['Folder'][i][0].text()
        state = args['Toggle'][i].checkState()
        axis = args['Axis'][i].currentText()
        type = args['Type'][i].currentText()

        try:
            scaleX = float(args['ScaleX'][i].text())
            scaleY = float(args['ScaleY'][i].text())
        except:
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
            x_data = [a*scaleX for a in self.other_data[requestX].tolist()]
            self.freq_glob = x_data

            for i in range(len(requestY)):
                y = [a*scaleY for a in self.other_data[requestY[i]].tolist()]
                if axis == 'Left':
                    if type == 'Line':
                        self.ax.plot(x_data, y, label=requestY[i], linewidth=2)
                    else:
                        self.ax.plot(x_data, y, linestyle='None', marker='x', markersize=10.0)
                    self.ax.set_ylabel('$S$ [dB]')
                    self.ax.set_xlabel('$f \mathrm{ [GHz]}$')
                else:
                    if type == 'Line':
                        self.ax_right.plot(x_data, self.other_data[requestY[i]].tolist(), label=requestY[i], linewidth=2)
                    else:
                        self.ax_right.plot(x_data, y, linestyle='None', marker='X', markersize=10.0)
                    self.ax_right.set_ylabel('$S$ [dB]')
        else:
            print("Please specify columns to plot")

    def table_control(self):
        # fill the first line
        self.create_new_row(0, self.plotUI.tableWidget)

    def create_new_row(self, row_ind, table):
        # Code
        cb_code = QComboBox()
        cb_code.addItem("ABCI")
        cb_code.addItem("SLANS")
        cb_code.addItem("Other")
        table.setCellWidget(row_ind, 0, cb_code)

        # Folder widget
        w_Folder = QWidget()
        l_Folder_Widget = QHBoxLayout()
        w_Folder.setLayout(l_Folder_Widget)
        le_Folder = QLineEdit()
        le_Folder.setReadOnly(True)
        pb_Open_Folder = QPushButton('...')
        # signal to change place holder text is 'Other'
        cb_code.currentIndexChanged.connect(lambda: le_Folder.setPlaceholderText('Select file') if cb_code.currentText() == 'Other' else le_Folder.setPlaceholderText('Select Folder'))

        l_Folder_Widget.setSpacing(0)
        l_Folder_Widget.setContentsMargins(0, 0, 0, 0)
        l_Folder_Widget.addWidget(le_Folder)
        l_Folder_Widget.addWidget(pb_Open_Folder)
        table.setCellWidget(row_ind, 1, w_Folder)

        # Polarisation
        cb_pol = QComboBox()
        cb_pol.addItem("Long")
        cb_pol.addItem("Trans")
        table.setCellWidget(row_ind, 2, cb_pol)
        # signal to disable polarisation when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(lambda: cb_pol.setDisabled(True) if cb_code.currentText() == 'Other' else cb_pol.setEnabled(True))

        # Id
        ccb_Id = CheckableComboBox()
        ccb_Id.addItem("All")
        table.setCellWidget(row_ind, 3, ccb_Id)
        # signal to disable id when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(lambda: ccb_Id.setDisabled(True) if cb_code.currentText() == 'Other' else ccb_Id.setEnabled(True))

        # push button signal
        pb_Open_Folder.clicked.connect(lambda: self.open_folder(le_Folder, cb_pol, ccb_Id, 'File') if cb_code.currentText() == 'Other' else self.open_folder(le_Folder, cb_pol, ccb_Id, 'Folder'))
        cb_pol.currentIndexChanged.connect(lambda: self.populate_IDs(le_Folder.text(), cb_pol, ccb_Id) if cb_code.currentText() == 'Other' else self.populate_IDs(le_Folder.text(), cb_pol, ccb_Id))

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
        cb_Y = CheckableComboBox()

        # place cb_request and le_request in request widget layout
        l_Request_Widget.addWidget(cb_request)
        l_Request_Widget.addWidget(cb_X)
        l_Request_Widget.addWidget(cb_Y)

        # hide le_request
        cb_X.hide()
        cb_Y.hide()

        # add widget to table
        table.setCellWidget(row_ind, 4, w_Request)

        # add signal to le_Folder after editing to update cb_X and cb_Y
        le_Folder.textChanged.connect(lambda: self.populate_combobox(cb_X, cb_Y, le_Folder.text()))

        # signal to disable polarisation when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(lambda: self.show_hide_request_widgets(cb_code, cb_request, cb_X, cb_Y))

        # add signal to switch between longitudinal and transverse
        cb_pol.currentIndexChanged.connect(lambda: self.switchRequest(cb_request, request_dict_long, request_dict_trans, cb_pol))

        # Toggle on/off
        cb_toggle = QCheckBox()
        table.setCellWidget(row_ind, 5, cb_toggle)

        # ScaleX
        scaleX = QLineEdit()
        scaleX.setText('1')
        table.setCellWidget(row_ind, 6, scaleX)

        # ScaleY
        scaleY = QLineEdit()
        scaleY.setText('1')
        table.setCellWidget(row_ind, 7, scaleY)

        # axis
        axis_list = ['Left', 'Right', 'Top', 'Bottom']
        cb_axis = QComboBox()
        for a in axis_list:
            cb_axis.addItem(a)

        table.setCellWidget(row_ind, 8, cb_axis)

        # type
        type_list = ['Line', 'Scatter', 'Bar']
        cb_type = QComboBox()
        for a in type_list:
            cb_type.addItem(a)

        table.setCellWidget(row_ind, 9, cb_type)

        self.args_dict['Code'].append(cb_code)
        self.args_dict['Folder'].append([le_Folder, pb_Open_Folder, w_Folder, l_Folder_Widget])
        self.args_dict['Polarization'].append(cb_pol)
        self.args_dict['Id'].append(ccb_Id)
        self.args_dict['Request'].append([cb_request, cb_X, cb_Y, w_Request, l_Request_Widget])
        self.args_dict['Toggle'].append(cb_toggle)
        self.args_dict['ScaleX'].append(scaleX)
        self.args_dict['ScaleY'].append(scaleY)
        self.args_dict['Axis'].append(cb_axis)
        self.args_dict['Type'].append(cb_type)

    def generateMenu(self, pos):
        # Get index
        menu = QMenu()
        item1 = menu.addAction("Add Row")
        item2 = menu.addAction("Delete Row")
        # Make the menu display in the normal position
        screenPos = self.plotUI.tableWidget.mapToGlobal(pos)

        # Click on a menu item to return, making it blocked
        action = menu.exec(screenPos)
        if action == item1:
            self.add_row()
        if action == item2:
            self.remove_row(self.row)
        else:
            return

    def add_row(self):
        #get current number of rows
        n = self.plotUI.tableWidget.rowCount()
        # add new row
        self.plotUI.tableWidget.setRowCount(n+1) # and one row in the table
        self.create_new_row(n, self.plotUI.tableWidget)

    def remove_row(self, row):
        # try:
        #     axis = self.args_dict['Axis'][row].currentText()
        #     print(axis)
        #     # select plot
        #     ax_selected = self.ax
        #
        #     if axis == 'Left':
        #         ax_selected = self.ax
        #         axis_index = self.ax_row_plot_dict[row]
        #         print('left: ', row, axis_index, len(self.ax.get_lines()))
        #
        #     else:
        #         ax_selected = self.ax_right
        #         axis_index = self.ax_right_row_plot_dict[row]
        #         print('right: ', row, axis_index, len(self.ax_right.get_lines()))
        #
        #     # remove associated plot
        #     for i in range(len(self.args_dict['Id'][row].text().split(','))):
        #         # convert row index to axis index
        #         ax_selected.get_lines()[axis_index].remove()
        #         try:
        #             self.axins.get_lines()[axis_index].remove()
        #         except:
        #             print("No inset axis.")
        #
        #     # print(ax_selected, len(ax_selected.get_lines()))
        #     # redraw legend
        #     self.ax.legend(loc='lower left')
        #     self.ax.autoscale(True, axis='y')
        #     self.ax.relim()
        #
        #     self.ax_right.legend(loc='lower left')
        #     self.ax_right.autoscale(True, axis='y')
        #     self.ax_right.relim()
        #
        #     # check if plot has no other lines and delete axis labels if so
        #     if len(self.ax.get_lines()) == 0:
        #         self.ax.cla()
        #
        #     if len(self.ax_right.get_lines()) == 0:
        #         self.ax_right.cla()
        #
        #     # toggle axis labels
        #     self.toggle_axis_labels()
        #
        #     self.plt.draw()
        # except:
        #     print("No associated plot(s).")

        # remove from table
        n = self.plotUI.tableWidget.rowCount()

        self.plotUI.tableWidget.removeRow(row)
        # reset number of rows
        self.plotUI.tableWidget.setRowCount(n-1)

        # remove elements from dict
        # delete objects
        del self.args_dict['Code'][row]
        del self.args_dict['Id'][row]
        del self.args_dict['Polarization'][row]

        for i in range(len(self.args_dict['Request'][row])):
            del self.args_dict['Request'][row][-1]
        # delete list afterwards
        del self.args_dict['Request'][row]

        for i in range(len(self.args_dict['Folder'][row])):
            del self.args_dict['Folder'][row][-1]
        # delete list afterwards
        del self.args_dict['Folder'][row]

        del self.args_dict['Toggle'][row]
        del self.args_dict['ScaleX'][row]
        del self.args_dict['ScaleY'][row]
        del self.args_dict['Axis'][row]
        del self.args_dict['Type'][row]

        # you couldn't live with your own failure. Where did that bring you? Back to me.
        # Referring to avoiding just replotting everything when a plot is removed. The long code commented out above
        self.plot()

    def switchRequest(self, cb_Request, request_dict_long, request_dict_trans, cb_Code):
        # clear combo box
        cb_Request.clear()
        print('here')
        # add new items
        if cb_Code.currentText() == 'Long':
            for req in request_dict_long.values():
                cb_Request.addItem(req)
        else:
            for req in request_dict_trans.values():
                cb_Request.addItem(req)

    def switchCode(self, cb_Code):
        if cb_Code.currentText() == 'ABCI':
            pass
        elif cb_Code.currentText() == 'SLANS':
            pass

    def mousePressEvent(self, event):
        if event.buttons() & Qt.RightButton:
            row = self.plotUI.tableWidget.currentRow()
            self.row = row

    def calc_cutoff(self, Ri, mode):
        # calculate frequency from Ri
        p_TM01, p_TE11 = 2.405, 1.841
        c = 299792458  # m/s

        if mode == 'TM01':
            freq = 400.79 * 1e-3
        else:
            freq = (c * p_TE11) / (2 * np.pi * Ri * 1e9) * 1e3

        return freq

    def calc_limits(self, mode):
        if len(self.freq_glob) > 0:
            if self.plotUI.cb_Longitudinal_Threshold.checkState() == 2 or self.plotUI.cb_Transverse_Threshold.checkState() == 2:
                E0 = [45.6, 80, 120, 182.5]  # [GeV] Energy
                nu_s = [0.025, 0.0506, 0.036, 0.087]  # Synchrotron oscillation tune
                I0 = [1390, 147, 29, 10.8]  # [mA] Beam current5.4 * 2
                alpha_c = [1.48, 1.48, 0.73, 0.73]  # [10−5] Momentum compaction factor
                tau_z = [424.6, 78.7, 23.4, 6.8]  # [ms] Longitudinal damping time
                tau_xy = [849.2, 157.4, 46.8, 13.6]  # [ms] Transverse damping time
                frev = [3.07, 3.07, 3.07, 3.07]  # [kHz] Revolution frequency
                beta_xy = 50

                #     Ncav = [52, 52, 136, 584] # Number of cavities per beam
                Ncav = [52, 100, 268, 584]  # 1_2_2_25

                E0 = ast.literal_eval(self.plotUI.le_E0.text())
                nu_s = ast.literal_eval(self.plotUI.le_Nu_S.text())
                I0 = ast.literal_eval(self.plotUI.le_I0.text())
                alpha_c = ast.literal_eval(self.plotUI.le_Alpha_S.text())
                tau_z = ast.literal_eval(self.plotUI.le_Tau_Z.text())
                tau_xy = ast.literal_eval(self.plotUI.le_Tau_XY.text())
                frev = ast.literal_eval(self.plotUI.le_F_Rev.text())
                beta_xy = ast.literal_eval(self.plotUI.le_Beta_XY.text())
                Ncav = ast.literal_eval(self.plotUI.le_N_Cav.text())

                Z_list = []
                if mode == 'monopole':
                    # trim f
                    # f_list = f_list[0:len(f_list) - 100]
                    f_list = self.freq_glob[0:len(self.freq_glob) - 100]
                    for i, n in enumerate(Ncav):
                        Z = [(2 * E0[i] * nu_s[i]) * 1e8 / (n * I0[i] * alpha_c[i] * tau_z[i] * a) if a > 1e-8 else 1e5 for a in
                             f_list]
                        Z_list.append(Z)

                    self.plot_baselines(f_list, Z_list, mode)

                    return f_list, Z_list

                elif mode == 'dipole':
                    f_list = self.freq_glob[0:len(self.freq_glob) - 100]
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
                # ab = self.ax.text(0.5, z+1.8, f'{text[i]}', transform = self.ax.transAxes, va='top')      # coordinate system transformation

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

    def other_plots(self):
        pass

    def toggle_page(self, key):

        if key == 'Plot Decorations':
            self.plotUI.sw_Plot_Area_Tools.setCurrentIndex(0)

        if key == 'Machine Parameters':
            self.plotUI.sw_Plot_Area_Tools.setCurrentIndex(1)

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

    def open_folder(self, le, cb_pol, ccb, mode='Folder'):
        if mode == "Folder":
            data_dir = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
            if data_dir != '':
                le.setText(data_dir)

                # populate checkable combobox
                self.populate_IDs(data_dir, cb_pol, ccb)
        else:
            filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Excel Files (*.xlsx)")
            if filename != '':
                le.setText(filename)

    def populate_IDs(self, dirc, cb_pol, ccb):
        if dirc != "":
            dir_list = os.listdir(dirc)

            # clear checkable check box
            ccb.clear()
            ccb.addItem("All")

            # polarisation
            pol = cb_pol.currentText()
            i = 0
            if pol == "Long":
                i = 0
            else:
                i = 1

            # loop through folder list to check which contained simulation files
            for d in dir_list:
                if os.path.exists(fr"{dirc}\{d}\Cavity_MROT_{i}.pot"):
                    ccb.addItem(d)

    def open_file(self):
        filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Excel Files (*.xlsx)")

        try:
            self.plotUI.le_Filename.setText(filename)
            df = pd.read_excel(filename)

            # get header
            header = df.columns

            # create combobox for x checkbox list for y
            cb_X_Axis = QComboBox()
            cb_Y_Axis = CheckableComboBox()
            cb_Y_Axis_List = []

            # populate x and y with checkboxes
            for col in header:
                cb_X_Axis.addItem(col)
                cb_Y_Axis.addItem(col)

                cb = QCheckBox(col)
                cb_Y_Axis_List.append(cb)

            self.plotUI.vl_Y_Axis_Checkboxes.addWidget(cb_Y_Axis)

            # add to layout
            self.plotUI.vl_X_Axis_Combobox.addWidget(cb_X_Axis)
            # self.plotUI.vl_X_Axis_Combobox.addWidget(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))


        except Exception as e:
            self.log.error(f'Failed to plot:: {e}')

    def set_table_size(self):
        self.plotUI.tableWidget.setColumnWidth(0, 75)
        self.plotUI.tableWidget.setColumnWidth(1, 200)
        self.plotUI.tableWidget.setColumnWidth(2, 100)
        self.plotUI.tableWidget.setColumnWidth(3, 300)
        self.plotUI.tableWidget.setColumnWidth(4, 200)
        self.plotUI.tableWidget.setColumnWidth(5, 75)
        self.plotUI.tableWidget.setColumnWidth(6, 75)
        self.plotUI.tableWidget.setColumnWidth(7, 75)

    def show_hide_request_widgets(self, cb_code, cb_request, cb_X, cb_Y):
        if cb_code.currentText() == 'Other':
            cb_request.hide()
            cb_X.show()
            cb_Y.show()
        else:
            cb_request.show()
            cb_X.hide()
            cb_Y.hide()

    def populate_combobox(self, cb_X, cb_Y, filename):
        # clear combobox item if any
        cb_X.clear()
        cb_Y.clear()

        if filename.split('.')[-1] == 'xlsx':
            # load file (for now, just excel files)
            df = fr.excel_reader(filename)
            sheet_name = list(df.keys())[0]

            self.other_data = df[sheet_name]
            for a in self.other_data.keys():
                cb_X.addItem(a)
                cb_Y.addItem(a)

    def serialize(self, state_dict):
        # update state file
        state_dict["Inlet_Position"] = self.plotUI.le_Inset_Position.text()
        state_dict["Inset_Window"] = self.plotUI.le_Inset_Window.text()
        state_dict["Show_Hide_Inset"] = self.plotUI.cb_Show_Hide_Inset.checkState()
        state_dict["Y_Scale"] = self.plotUI.cb_Y_Scale.currentText()

        state_dict["Show_Wakefield_Parameters"] = self.plotUI.cb_Show_Wakefield_Parameters.checkState()

        state_dict["Text"] = self.plotUI.le_Plot_Text.text()
        state_dict["Peaks Threshold"] = self.plotUI.dsb_Threshold.value()

        state_dict["E0"] = self.plotUI.le_E0.text()
        state_dict["Nu_S"] = self.plotUI.le_Nu_S.text()
        state_dict["I0"] = self.plotUI.le_I0.text()
        state_dict["Alpha_S"] = self.plotUI.le_Alpha_S.text()
        state_dict["Tau_Z"] = self.plotUI.le_Tau_Z.text()
        state_dict["Tau_XY"] = self.plotUI.le_Tau_XY.text()
        state_dict["F_Rev"] = self.plotUI.le_F_Rev.text()
        state_dict["Beta_XY"] = self.plotUI.le_Beta_XY.text()
        state_dict["N_Cav"] = self.plotUI.le_N_Cav.text()
        state_dict["Longitudinal Threshold Checkbox"] = self.plotUI.cb_Longitudinal_Threshold.checkState()
        state_dict["Transverse Threshold Checkbox"] = self.plotUI.cb_Transverse_Threshold.checkState()
        state_dict["Threshold Line Color"] = self.plotUI.cb_Threshold_Line_Color.currentText()
        state_dict["Threshold Linestyle"] = self.plotUI.cb_Threshold_Linestyle.currentText()

    def deserialize(self, state_dict):
        self.plotUI.le_Inset_Position.setText(state_dict["Inlet_Position"])
        self.plotUI.le_Inset_Window.setText( state_dict["Inset_Window"])

        self.plotUI.cb_Show_Hide_Inset.setCheckState(state_dict["Show_Hide_Inset"])

        self.plotUI.cb_Y_Scale.setCurrentText(state_dict["Y_Scale"])

        self.plotUI.cb_Show_Wakefield_Parameters.setCheckState(state_dict["Show_Wakefield_Parameters"])

        self.plotUI.le_Plot_Text.setText(state_dict["Text"])
        self.plotUI.dsb_Threshold.setValue( state_dict["Peaks Threshold"])

        self.plotUI.le_E0.setText(state_dict["E0"])
        self.plotUI.le_Nu_S.setText(state_dict["Nu_S"])
        self.plotUI.le_I0.setText(state_dict["I0"])
        self.plotUI.le_Alpha_S.setText(state_dict["Alpha_S"])
        self.plotUI.le_Tau_Z.setText(state_dict["Tau_Z"])
        self.plotUI.le_Tau_XY.setText(state_dict["Tau_XY"])
        self.plotUI.le_F_Rev.setText(state_dict["F_Rev"])
        self.plotUI.le_Beta_XY.setText(state_dict["Beta_XY"])
        self.plotUI.le_N_Cav.setText(state_dict["N_Cav"])
        self.plotUI.cb_Longitudinal_Threshold.setCheckState(state_dict["Longitudinal Threshold Checkbox"])
        self.plotUI.cb_Transverse_Threshold.setCheckState(state_dict["Transverse Threshold Checkbox"])
        self.plotUI.cb_Threshold_Line_Color.setCurrentText(state_dict["Threshold Line Color"])
        self.plotUI.cb_Threshold_Linestyle.setCurrentText(state_dict["Threshold Linestyle"])


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
