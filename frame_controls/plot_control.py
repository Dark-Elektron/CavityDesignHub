import itertools
import os.path
import matplotlib
import pandas as pd
import scipy
from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from scipy.special import jn_zeros, jnp_zeros
from ui_files.plot import Ui_Plot
from ui_files.plottypeselector import Ui_PlotTypeSelector
from analysis_modules.plot_module.plotter import Plot
from analysis_modules.data_module.abci_data import ABCIData
from utils.shared_classes import *
from utils.shared_functions import *
import scipy.io as spio

fr = FileReader()

file_color = 'yellow'
DEBUG = True

MATPLOTLIB_COLORS = {'Default': plt.rcParams['axes.prop_cycle'].by_key()['color'],
                     'fivethirtyeight': ["#008fd5", "#fc4f30", "#e5ae38", "#6d904f", "#8b8b8b", "#810f7c"],
                     'ggplot': ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8'],
                     'science': ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e'],
                     'seaborn': ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'],
                     'seaborn-bright': ['#003FFF', '#03ED3A', '#E8000B', '#8A2BE2', '#FFC400', '#00D7FF'],
                     'seaborn-colorblind': ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9'],
                     'seaborn-dark-palette': ['#001C7F', '#017517', '#8C0900', '#7600A1', '#B8860B', '#006374'],
                     'seaborn-deep': ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'],
                     'seaborn-muted': ['#4878CF', '#6ACC65', '#D65F5F', '#B47CC7', '#C4AD66', '#77BEDB'],
                     'seaborn-pastel': ['#92C6FF', '#97F0AA', '#FF9F9A', '#D0BBFF', '#FFFEA3', '#B0E0E6'],
                     'Solarize_Light2': ['#268BD2', '#2AA198', '#859900', '#B58900', '#CB4B16', '#DC322F', '#D33682',
                                         '#6C71C4'],
                     'bmh': ['#348ABD', '#A60628', '#7A68A6', '#467821', '#D55E00', '#CC79A7', '#56B4E9', '#009E73',
                             '#F0E442',
                             '#0072B2'],
                     'bright': ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'],
                     'classic': ['b', 'g', 'r', 'c', 'm', 'y', 'k'],
                     'dark_background': ['#8dd3c7', '#feffb3', '#bfbbd9', '#fa8174', '#81b1d2', '#fdb462', '#b3de69',
                                         '#bc82bd',
                                         '#ccebc4', '#ffed6f'],
                     'grayscale': ['0.00', '0.40', '0.60', '0.70'],
                     'high-contrast': ['#004488', '#DDAA33', '#BB5566'],
                     'high-vis': ['#0d49fb', '#e6091c', '#26eb47', '#8936df', '#fec32d', '#25d7fd'],
                     'ieee': ['k', 'r', 'b', 'g'],
                     'light': ['#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00',
                               '#DDDDDD'],
                     'muted': ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933',
                               '#AA4499', '#DDDDDD'],
                     'retro': ['#4165c0', '#e770a2', '#5ac3be', '#696969', '#f79a1e', '#ba7dcd'],
                     'scatter': ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e'],
                     'std-colors': ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e'],
                     'tableau-colorblind10': ['#006BA4', '#FF800E', '#ABABAB', '#595959', '#5F9ED1', '#C85200',
                                              '#898989',
                                              '#A2C8EC', '#FFBC79', '#CFCFCF'],
                     'vibrant': ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']
                     }


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class PlotControl:
    """
    Controls plotting
    """

    def __init__(self, parent):
        self.w_Plot = QWidget()

        self.ui = Ui_Plot()
        self.ui.setupUi(self.w_Plot)

        # Create main window object
        self.win = parent
        self.main_control = parent
        self.main_ui = parent.ui

        # get logger
        self.log = self.main_control.log

        # Plot objects and variables
        self.plt = Plot(self)
        self.ui.gl_Plot_Area.addWidget(self.plt)
        self.fig = self.plt.fig
        self.ax = self.plt.ax
        self.initial_prop_cycler, _ = itertools.tee(self.ax._get_lines.prop_cycler)
        self.ax_right = self.plt.ax_right
        self.axins = None
        self.indicate_inset = None
        self.leg = None
        self.le_Color = QLineEdit()
        self.le_Color.setReadOnly(True)
        self.pb_Color = self.matplotlib_colors()
        self.ui.gl_Color.addWidget(self.le_Color, 0, 0, 1, 1)
        self.ui.gl_Color.addWidget(self.pb_Color, 1, 0, 1, 1)
        self.ui.gl_Custom_Color.addWidget(self.plt.w_Color)

        # class variables
        self.plotID_count = 0

        # class lists
        self.freq_glob = []
        self.baseline_line_objects = []

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
        self.plot_elements_table = {}
        self.operation_points = pd.DataFrame()
        # baseline matplotlib line.Line2D objects
        self.axes_decorations_dict = {"xlabel": "x",
                                      "ylabel": "y",
                                      "ylabel2": "",
                                      "title": "Title",
                                      "legend": []}
        self.other_data_filtered = {}
        self.row = None
        
        # setup QTableWidget
        self.tableWidget = TableWidgetDragDrop(self)
        self.ui.gridLayout_29.addWidget(self.tableWidget, 1, 0, 1, 5)

        self.initUI()
        self.signals()

    def matplotlib_colors(self):
        pb_Color = QPushButton()
        menu = QMenu()
        menu.triggered.connect(lambda x: pb_Color.setStyleSheet(f"background-color: {x.defaultWidget().text()};"))
        menu.triggered.connect(lambda x: pb_Color.setText(x.defaultWidget().text()))
        menu.triggered.connect(lambda x: self.le_Color.setText(x.defaultWidget().text()))
        menu.triggered.connect(lambda x: self.le_Color.setStyleSheet(f"background-color: {x.defaultWidget().text()};"))
        pb_Color.setMenu(menu)

        for k, vals in MATPLOTLIB_COLORS.items():
            sub_menu = menu.addMenu(k)
            for v in vals:
                label = QLabel(str(v))
                action = QWidgetAction(sub_menu)
                label.setStyleSheet(f"background-color: {v};")
                action.setDefaultWidget(label)
                sub_menu.addAction(action)

        return pb_Color

    def createPlotTypeWidget(self):
        # plottypeselector
        w_PTS = QWidget()
        PTS = Ui_PlotTypeSelector()
        PTS.setupUi(w_PTS)

        plotID = self.plotID_count
        # add plottypse selector to dictionary
        self.pts_dict[self.plotID_count] = {"plot inputs": None, "plot data": {}, "plot object": {},
                                            "plot data inputs": None}

        self.init_abci(PTS)
        self.init_slans(PTS)
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

    def init_abci(self, PTS):
        """

        :param PTS:
        :return:
        """
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

    def init_slans(self, PTS):
        """

        :param PTS:
        :return:
        """
        # slans
        PTS.pb_Open_Folder_SLANS.clicked.connect(
            lambda: self.open_(PTS.le_Folder_SLANS, PTS.cb_Polarization_SLANS, PTS.ccb_Shape_ID_SLANS, mode='Folder'))

    def init_other(self, PTS, plotID):
        """

        Parameters
        ----------
        plotID
        PTS

        Returns
        -------

        """
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
        """

        Parameters
        ----------
        PTS
        plotID

        Returns
        -------

        """
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
        # plottype = self.pts_dict[pid]["plot inputs"].cb_code.currentText()
        if switch.checkState() == 2:

            # check if plot already exist
            if self.pts_dict[pid]['plot data'] != {}:
                # compare plot input with current input
                inputs = self.pts_dict[pid]["plot inputs"]
                ids = [a.strip() for a in inputs.ccb_Id.currentText().split(',')]  # get list
                pol = inputs.cb_Pol.currentText()
                request = inputs.cb_Request.currentText()

                plot_data_inputs = self.pts_dict[pid]["plot data inputs"]

                if plot_data_inputs["ids"] == ids and plot_data_inputs["pol"] == pol and \
                        plot_data_inputs["request"] == request:
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
        # cycle, _ = itertools.tee(self.initial_prop_cycler)
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
        except ValueError:
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
                                                                  label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                                        f'({abci_data_long.wakelength} m)',
                                                                  linewidth=2)})

                elif axis.lower() == "right":
                    ax_selected = self.ax_right
                    # update plot data
                    pts_["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    pts_["plot object"].update({id_: self.ax_right.plot(xr, y,
                                                                        label='${' + f'{str(id_).replace("_", ",")}' +
                                                                              '}$' +
                                                                              f' ({abci_data_long.wakelength} m)',
                                                                        linewidth=2)})

                # ax_selected.set_xlabel('$f \mathrm{ [GHz]}$')
                # ax_selected.set_ylabel('$Z_{\parallel, \mathrm{HOM}} \mathrm{[k\Omega]}$')
                ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

            else:
                abci_data_trans = ABCIData(abci_data_dir, id_, 1)
                if request == 'Transverse Impedance Magnitude':
                    # try:
                    xr, yr, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                    xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')
                    # except:
                    #     xr, yr, _ = abci_data_trans.get_data('Real Part of Azimuthal Impedance')
                    #     xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Azimuthal Impedance')

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
                                                                  label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                                        f' ({abci_data_trans.wakelength} m)',
                                                                  linewidth=2)})

                elif axis.lower() == 'right':
                    ax_selected = self.ax_right
                    # update plot data
                    pts_["plot data"].update({id_: {"x": xr, "y": y}})
                    # process id
                    pts_["plot object"].update({id_: self.ax_right.plot(xr, y,
                                                                        label='${' + f'{str(id_).replace("_", ",")}' +
                                                                              '}$' +
                                                                              f' ({abci_data_trans.wakelength} m)',
                                                                        linewidth=2)})

                ax_selected.set_yscale('log')
                ax_selected.set_ylim(min(y), max(y))
                ax_selected.set_xlim(min(xr), max(xr))

    def signals(self):
        self.ui.pb_Reset_Colors.clicked.connect(lambda: self.reset_colors())
        # self.ui.pb_Add.clicked.connect(lambda: self.createPlotTypeWidget())
        self.le_Color.textChanged.connect(lambda: self.plt.update_object_properties({'color': self.le_Color.text()}))
        self.ui.dsb_Line_Width.valueChanged.connect(
            lambda: self.plt.update_object_properties({'lw': self.ui.dsb_Line_Width.value()}))
        self.ui.cb_Line_Style.currentTextChanged.connect(
            lambda: self.plt.update_object_properties({'ls': self.ui.cb_Line_Style.currentText()}))
        self.ui.cb_Marker.currentTextChanged.connect(
            lambda: self.plt.update_object_properties({'marker': self.ui.cb_Marker.currentText(), 'mec': 'k'}))
        self.ui.dsb_Marker_Size.valueChanged.connect(
            lambda: self.plt.update_object_properties({'ms': self.ui.dsb_Marker_Size.value()}))
        self.ui.dsb_Alpha.valueChanged.connect(
            lambda: self.plt.update_object_properties({'alpha': self.ui.dsb_Alpha.value()}))
        self.ui.sb_Layer.valueChanged.connect(
            lambda: self.plt.update_object_properties({'zorder': self.ui.sb_Layer.value()}))

        # plot abci impedance
        self.ui.pb_Plot.clicked.connect(lambda: self.plot())

        # refresh plot
        self.ui.pb_Refresh.clicked.connect(lambda: self.plot())

        # toggle plot menu
        # self.ui.pb_Plot_Area_Menu.clicked.connect(lambda: animate_width(self.ui.w_Plot_Menu, 0, 500, True))
        # self.ui.pb_Plot_Area_Menu.clicked.connect(lambda: self.ui.w_Plot_Area_Buttons.setEnabled(
        #     True) if self.ui.w_Plot_Menu.maximumWidth() == 0 else self.ui.w_Plot_Area_Buttons.setDisabled(True))

        # signal for plot menu pages
        self.ui.pb_Machine_Parameters.clicked.connect(lambda: self.toggle_page('Machine Parameters'))
        self.ui.pb_Plot_Decorations.clicked.connect(lambda: self.toggle_page('Plot Decorations'))
        self.ui.pb_Plot_Element_Format.clicked.connect(lambda: self.toggle_page('Plot Element Format'))

        # signal for plot argument entrmake widgets rey self.ui.pb_Collapse_Shape_Parameters.clicked.connect(lambda:
        # animate_height(self.ui.w_Plot_Args, 0, 200, True))

        # signal for threshold
        self.ui.pb_Load_Machine_Parameters.clicked.connect(
            lambda: open_file(self.ui.le_Machine_Parameters_Filename,
                              start_folder=str(self.main_control.projectDir / "OperatingPoints")))
        self.ui.cb_Longitudinal_Threshold.stateChanged.connect(
            lambda: self.calc_limits('monopole', self.ui.ccb_Operation_Points.currentText()))
        self.ui.cb_Transverse_Threshold.stateChanged.connect(
            lambda: self.calc_limits('dipole', self.ui.ccb_Operation_Points.currentText()))

        self.ui.le_Machine_Parameters_Filename.textChanged.connect(
            lambda: self.load_operating_points(self.ui.le_Machine_Parameters_Filename.text()))

        # signal for inset plot
        self.ui.cb_Show_Hide_Inset.clicked.connect(lambda: self.plot_inset())

        # add plot
        self.ui.pb_Add_Row.clicked.connect(self.tableWidget.add_row)
        # clear plots
        self.ui.pb_Clear.clicked.connect(lambda: self.clear_plots())
        # clear lines
        self.ui.pb_Clear_Lines.clicked.connect(lambda: self.clear_lines())
        # clear texts
        self.ui.pb_Clear_Texts.clicked.connect(lambda: self.clear_texts())

        # signals for cutoff
        # self.ui.cb_Check_1.clicked.connect(lambda: self.plot_cutoff(0, self.ui.cb_Check_1))
        # self.ui.cb_Check_2.clicked.connect(lambda: self.plot_cutoff(1, self.ui.cb_Check_2))
        # self.ui.cb_Check_3.clicked.connect(lambda: self.plot_cutoff(2, self.ui.cb_Check_1))

        # update label
        self.ui.pb_Apply_Axis_Labels.clicked.connect(lambda: self.update_labels())

        # plot text
        # x, y = self.plt.ax.get_xlim()[1]/2, self.plt.ax.get_ylim()[1]/2
        self.ui.pb_Add_Textt.clicked.connect(lambda: self.plt.add_text(
            text=self.ui.le_Plot_Text_2.text(), box=self.ui.cb_Text_Box.currentText(),
            size=self.ui.sb_Annotation_Text_Size.value(), rotation=self.ui.sb_Rotation.value(),
            xycoords='axes fraction'))

        # plot axvline
        self.ui.pb_Add_Vline.clicked.connect(lambda: self.plt.add_axvline(self.ui.dsb_Axvline_X.value()))

        # add rectangular patch
        self.ui.pb_Add_Patch.clicked.connect(lambda: self.plt.add_patch(
            (self.ui.le_Patch_Pos_X.text(), self.ui.le_Patch_Pos_Y.text()), self.ui.le_Patch_Size_A.text(),
            self.ui.le_Patch_Size_B.text(),
            self.ui.cb_Patch_Type.currentText().lower()))

        # switch axis
        self.ui.cb_Active_Axis.currentTextChanged.connect(lambda: self.switch_axis())

        # frequently used plot labels combobox signal
        # self.ui.ccb_Freq_Used_Xlabel.currentTextChanged.connect(lambda: self.le_Xlabel.setText(
        #     self.ui.ccb_Freq_Used_Xlabel.currentText()))
        # self.ui.ccb_Freq_Used_Ylabel.currentTextChanged.connect(lambda: self.le_Ylabel.setText(
        #     self.ui.ccb_Freq_Used_Ylabel.currentText()))

        self.ui.ccb_Operation_Points.currentTextChanged.connect(
            lambda: self.operation_points_selection(self.ui.ccb_Operation_Points.currentText()))

        # change plot mode to draw
        self.ui.cb_Draw.stateChanged.connect(lambda: self.plt.change_plot_mode(self.ui.cb_Draw.checkState()))

    def initUI(self):
        # self.createPlotTypeWidget()
        self.ui.w_Custom_Color.hide()
        # add row to table
        self.tableWidget.setRowCount(1)  # and one row in the table
        self.table_control()

        # set plot menu initial size to zero and disable plot buttons
        # self.ui.w_Plot_Menu.setFixedWidth(0)
        # self.ui.w_Plot_Area_Buttons.setDisabled(True)

        # tableWidget initializations
        self.tableWidget.mousePressEvent = self.mousePressEvent
        self.tableWidget.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableWidget.customContextMenuRequested.connect(self.generateMenu)
        self.set_table_size()

        # add checkable combobox
        cutoff = QCheckableComboBox()
        cutoff.addItem("All")
        cutoff.setMinimumWidth(150)
        self.ui.gl_Cutoff.addWidget(cutoff, 2, 1, 1, 1)
        self.populate_cutoff_combobox(cutoff)

        self.ui.le_Ri_Cutoff_List.editingFinished.connect(lambda: self.populate_cutoff_combobox(cutoff))

        # add signal
        cutoff.currentTextChanged.connect(lambda: self.plot_cutoff(cutoff))

        self.ui.w_Machines_View.hide()

        self.populate_plot_elements_table()
        self.set_latex_text_color()

    def plot(self):
        # try:
        # args = list(self.args_dict.values())

        # use same color cycler for both axes
        # self.ax._get_lines.prop_cycler, _ = itertools.tee(self.initial_prop_cycler)
        # self.ax_right._get_lines.prop_cycler, _ = itertools.tee(self.initial_prop_cycler)
        self.ax_right._get_lines.prop_cycler = self.ax._get_lines.prop_cycler
        plot_count = 1
        for key, val in self.plot_dict.items():
            code = self.plot_dict[key]['plot inputs']['Code'].currentText()
            # check plot type
            if code == 'ABCI':
                try:
                    self.make_abci_plot(key)
                except KeyError:
                    print_("A key error occurred. Please make sure that a directory containing valid ABCI output"
                           "files is loaded.")
                    return

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
            self.leg = self.ax.legend(lines + lines2, labels + labels2, edgecolor='k', facecolor='none', loc='lower left', prop={'size': 18})

            if self.ax_right.get_legend() is not None:
                self.ax_right.get_legend().remove()
        else:
            self.leg = self.ax_right.legend(lines + lines2, labels + labels2, edgecolor='k', facecolor='none', loc='lower left', prop={'size': 18})
            if self.ax.get_legend() is not None:
                self.ax.get_legend().remove()

        self.leg.set_zorder(10)
        self.leg.set_draggable(state=True, use_blit=True)

        # plot inset if check box is checked
        self.plot_inset()

        # plot thresholds if threshold is checked
        if self.ui.cb_Longitudinal_Threshold.checkState() == 2:
            self.calc_limits('monopole')
        if self.ui.cb_Transverse_Threshold.checkState() == 2:
            self.calc_limits('dipole')

        self.update_labels()
        self.fig.canvas.draw_idle()

        # except Exception as e:
        #     self.log.error("Please enter a valid argument: Exception: ", e)

    def clear_lines(self):
        keys = self.plot_dict.keys()
        # reset plot dict
        for key in keys:
            self.plot_dict[key]["plot data"] = {}
            self.plot_dict[key]["plot object"] = {}

        # remove lines from axis
        lines = self.ax.get_lines()
        for line in lines:
            line.remove()

        # reset color cycler
        self.ax.set_prop_cycle(None)
        self.ax_right.set_prop_cycle(None)

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def clear_texts(self):
        # remove texts from axis
        annotations = [child for child in self.ax.get_children() if isinstance(child, matplotlib.text.Annotation)]
        texts_count = len(annotations)
        for i in range(texts_count):
            annotations[i].remove()

        self.plt.text_dict = {}
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def clear_patches(self):
        pass

    def clear_plots(self):

        self.ax.cla()
        self.ax_right.cla()

        if self.axins is not None:
            self.axins.cla()
            self.axins.remove()
            self.axins = None

        keys = self.plot_dict.keys()
        # reset plot dict
        for key in keys:
            self.plot_dict[key]["plot data"] = {}
            self.plot_dict[key]["plot object"] = {}

        # clear annotations
        self.plt.clear()

        keys = self.ax_obj_dict.keys()
        # clear cutoff frequency lines
        for key in keys:
            del self.ax_obj_dict[key]
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
        """

        :param plot_dict:
        :return:
        """
        args = plot_dict["plot inputs"]
        ids = [a.strip() for a in args['Id'].currentText().split(',')]  # get list
        pol = args['Polarization'].currentText()
        request = args['Request'][0].currentText()
        folder = args['Folder'][0].text()
        # state = args['Toggle'].checkState()
        axis = args['Axis'].currentText()
        type_ = []
        t = args['Type']
        type_.append(t[0].currentText())
        type_.append(t[1].currentText())

        filter_ = []
        f = args['Filter']
        filter_.append(f[0].lineEdit().text())
        filter_.append(f[1].text())

        # try:
        scaleX = float(args['ScaleX'].text())
        scaleY = float(args['ScaleY'].text())
        # except:
        #     scaleX = 1
        #     scaleY = 1

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
                    # xr, yr, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                    # xi, yi, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')
                    #
                    # # y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]
                    # y = (np.array(yr) ** 2 + np.array(yi) ** 2) ** 0.5
                    xr, y, _ = abci_data_long.get_data('Longitudinal Impedance Magnitude')

                    # scale x axis]
                    # xr = [a * scaleX for a in xr]
                    xr = np.array(xr)*scaleX
                    # scale y axis]
                    # y = [a * scaleY for a in y]
                    y = np.array(y)*scaleX

                    # set global frequency only if the impedance magnitude is requested
                    self.freq_glob = xr
                else:
                    xr, y, _ = abci_data_long.get_data(request)

                    # to integrate better later. This part of the code is to check the frequency distribution of
                    # HOM power in the cavity
                    if request == 'Frequency Spectrum of Loss Factor':
                        e = 1.602e-19
                        Nb = 2.26e11
                        I0 = 5e-3
                        Pf = np.array(y) * I0 * e * Nb  # k_hom*I0*e*Nb

                        # # Define bin edges and calculate bin centers bin_edges = np.linspace(np.array(xr).min(),
                        # np.array(xr).max(), num=int((np.array(xr).min() + np.array(xr).max())/2)) bin_centers = 0.5
                        # * (bin_edges[1:] + bin_edges[:-1])
                        #
                        # # Bin the power data
                        # power_binned, _ = np.histogram(xr, bins=bin_edges, weights=Pf)
                        #
                        # # Plot the binned data
                        # plt.bar(bin_centers, power_binned, width=np.diff(bin_edges), align='center')
                        plt.bar(xr, Pf)
                        plt.xlabel('Frequency')
                        plt.ylabel('Power (summed in bins)')
                        plt.show()

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

                    # check if peaks is selected
                    if self.ui.cb_Plot_Peaks.checkState() == 2:
                        peaks, _ = scipy.signal.find_peaks(y, height=self.ui.dsb_Threshold.value())
                        plot_dict["plot object"].update({
                            id_: self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                           f' ({abci_data_long.wakelength} m)',
                                              linewidth=2, picker=True, markevery=peaks, marker='o', mec='k')})
                    else:
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
                    if self.ui.cb_Plot_Peaks.checkState() == 2:
                        peaks, _ = scipy.signal.find_peaks(y, height=self.ui.dsb_Threshold.value())
                        plot_dict["plot object"].update({
                            id_: self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                           f' ({abci_data_long.wakelength} m)',
                                              linewidth=2, picker=True, markevery=peaks, marker='o', mec='k')})
                    else:
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
                    # try:
                    # xr, yr, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                    #                     # xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')
                    # # except:
                    # #     xr, yr, _ = abci_data_trans.get_data('Real Part of Azimuthal Impedance')
                    # #     xi, yi, _ = abci_data_trans.get_data('Imaginary Part of Azimuthal Impedance')
                    #
                    # # y = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr, yi)]
                    # y = (np.array(yr) ** 2 + np.array(yi) ** 2) ** 0.5

                    xr, y, _ = abci_data_trans.get_data('Transversal Impedance Magnitude')
                    # scale x axis]
                    # xr = [a * scaleX for a in xr]
                    xr = np.array(xr)*scaleX
                    # scale y axis]
                    # y = [a * scaleY for a in y]
                    y = np.array(y)*scaleX

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
                    if self.ui.cb_Plot_Peaks.checkState() == 2:
                        peaks, _ = scipy.signal.find_peaks(y, height=self.ui.dsb_Threshold.value())
                        plot_dict["plot object"].update({
                            id_: self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                           f' ({abci_data_trans.wakelength} m)',
                                              linewidth=2, picker=True, markevery=peaks, marker='o', mec='k')})
                    else:
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
                    if self.ui.cb_Plot_Peaks.checkState() == 2:
                        peaks, _ = scipy.signal.find_peaks(y, height=self.ui.dsb_Threshold.value())
                        plot_dict["plot object"].update({
                            id_: self.ax.plot(xr, y, label='${' + f'{str(id_).replace("_", ",")}' + '}$' +
                                                           f' ({abci_data_trans.wakelength} m)',
                                              linewidth=2, picker=True, markevery=peaks, marker='o', mec='k')})
                    else:
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
        """

        :param pid:
        :return:
        """
        # check state of switch
        switch = self.plot_dict[pid]["plot inputs"]["Toggle"]
        # plottype = self.plot_dict[pid]["plot inputs"]["Code"].currentText()
        if switch.checkState() == 2:
            # check if plot already exist
            if self.plot_dict[pid]['plot data'] != {}:
                # compare plot input with current input
                args = self.plot_dict[pid]["plot inputs"]
                ids = [a.strip() for a in args['Id'].currentText().split(',')]  # get list
                pol = args['Polarization'].currentText()
                request = args['Request'][0].currentText()
                folder = args['Folder'][0].text()
                # state = args['Toggle'].checkState()
                axis = args['Axis'].currentText()

                type_ = []
                filter_ = []
                t = args['Type']

                type_.append(t[0].currentText())
                type_.append(t[1].currentText())

                f = args['Filter']
                filter_.append(f[0].lineEdit().text())
                filter_.append(f[1].text())

                try:
                    scaleX = float(args['ScaleX'].text())
                    scaleY = float(args['ScaleY'].text())
                except ValueError:
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
                self.plot_impedance(self.plot_dict[pid])

        else:
            for line2D in self.plot_dict[pid]['plot object'].values():
                line2D[0].set(alpha=0)

        self.update_labels()
        self.fig.canvas.draw_idle()

    def make_other_plot(self, pid):
        # check state of switch
        switch = self.plot_dict[pid]["plot inputs"]["Toggle"]
        # plottype = self.plot_dict[pid]["plot inputs"]["Code"].currentText()
        if switch.checkState() == 2:
            # check if plot already exist
            if self.plot_dict[pid]['plot data'] != {}:
                # compare plot input with current input
                args = self.plot_dict[pid]["plot inputs"]
                ids = [a.strip() for a in args['Id'].currentText().split(',')]  # get list
                filename = args['Folder'][0].text()
                # state = args['Toggle'].checkState()
                axis = args['Axis'].currentText()

                requestX = args['Request'][1].currentText()
                requestY = args['Request'][2].currentText().split(', ')

                type_ = []
                filter_ = []
                t = args['Type']
                type_.append(t[0].currentText())
                type_.append(t[1].currentText())

                f = args['Filter']
                filter_.append(f[0].lineEdit().text())
                filter_.append(f[1].text())

                try:
                    scaleX = float(args['ScaleX'].text())
                    scaleY = float(args['ScaleY'].text())
                except ValueError:
                    scaleX = 1
                    scaleY = 1

                inputs_compare = {"Ids": ids, "RequestX": requestX, "RequestY": requestY, "Folder": filename,
                                  "Axis": axis, "ScaleX": scaleX, "ScaleY": scaleY, "Type": type_, "Filter": filter_}

                plot_data_inputs = self.plot_dict[pid]["plot data inputs"]

                if plot_data_inputs == inputs_compare:
                    if self.plot_dict[pid]["plot object"] != {}:
                        for line2D in self.plot_dict[pid]['plot object'].values():
                            line2D[0][0].set(alpha=1)
                        self.fig.canvas.draw_idle()
                    else:
                        # plot with modified inputs
                        self.plot_other(self.plot_dict[pid], pid)
                else:
                    # pop previous entry from plot and dictionary
                    for vals in self.plot_dict[pid]['plot object'].values():
                        for line2D in vals.values():
                            line2D[0].remove()

                    # reset plot data, plot object and plot data inputs
                    self.plot_dict[pid]["plot data"] = {}
                    self.plot_dict[pid]["plot object"] = {}
                    self.plot_dict[pid]["plot data inputs"] = {}

                    self.plot_other(self.plot_dict[pid], pid)

            else:
                self.plot_other(self.plot_dict[pid], pid)

        else:
            for vals in self.plot_dict[pid]['plot object'].values():
                for line2D in vals.values():
                    line2D[0].set(alpha=0)

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
                        if self.ui.cb_Plot_Peaks.checkState() == 2:
                            peaks, _ = scipy.signal.find_peaks(line.get_ydata(), height=self.ui.dsb_Threshold.value())
                            self.axins.plot(line.get_xdata(), line.get_ydata(), ls=line.get_linestyle(),
                                            linewidth=line.get_lw(), c=line.get_color(), markevery=peaks, marker='o',
                                            mec='k')

                        else:
                            self.axins.plot(line.get_xdata(), line.get_ydata(), ls=line.get_linestyle(),
                                            linewidth=line.get_lw(), c=line.get_color())

            # sub region of the original image
            # get values from line edit
            try:
                x0, x1, y0, y1 = ast.literal_eval(self.ui.le_Inset_Window.text())
            except ValueError:
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
        except ValueError:
            scaleX = 1
            scaleY = 1

        # if filename == '':
        #     folder = fr'{self.main_control.projectDir}/SimulationData/ABCI'

        # # load file (for now, just excel files)
        # df = fr.excel_reader(filename)
        # sheet_name = list(df.keys())[0]
        #
        # data = df[sheet_name]

        if requestY != [] and requestX != []:
            if filename.split('.')[-1] == 'xlsx':
                # get sheets
                sheets = [a.strip() for a in self.pts_dict[pid]["plot inputs"].ccb_Id_Other.currentText().split(',')]
                for sh in sheets:
                    filter_ = self.pts_dict[pid]["plot inputs"].ccb_Filter_Other.currentText()
                    value = self.pts_dict[pid]["plot inputs"].le_Filter_Value_Other.text()

                    if filter_ == "None" or value == "":
                        self.other_data_filtered = self.pts_dict[pid]["plot data"][sh]
                    else:
                        self.other_data_filtered = self.pts_dict[pid]["plot data"][sh][
                            self.pts_dict[pid]["plot data"][sh][filter_] == value]

                    x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                    self.freq_glob = x_data

                    for j in range(len(requestY)):
                        y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]

                        if axis == 'Left':
                            if type_ == 'Line':
                                self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                            else:
                                self.ax.plot(x_data, y, linestyle='None', marker=style, markersize=15.0,
                                             markeredgecolor="black", label="Legend")
                            self.ax.set_ylabel('$Y$ []')
                            self.ax.set_xlabel('$X$ []')
                        else:
                            if type_ == 'Line':
                                self.ax_right.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                            else:
                                self.ax_right.plot(x_data, y, linestyle='None', markeredgecolor="black", marker=style,
                                                   markersize=15.0,
                                                   label="Legend")
                            self.ax_right.set_ylabel('$Y$ [dB]')
            else:
                # try to filter self.other_data
                try:
                    filter_ = self.pts_dict[pid]["plot inputs"].ccb_Filter_Other.currentText()
                    value = self.pts_dict[pid]["plot inputs"].le_Filter_Value_Other.text()
                    if filter_ == "None" or value == "":
                        self.other_data_filtered = self.pts_dict[pid]["plot data"]
                    else:
                        self.other_data_filtered = self.pts_dict[pid]["plot data"][
                            self.pts_dict[pid]["plot data"][filter_] == value]
                except Exception as e:
                    print_("plot_control.py:: Exception:: ", e)
                    self.other_data_filtered = self.pts_dict[pid]["plot data"]

                x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                self.freq_glob = x_data

                for j in range(len(requestY)):
                    y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                    if axis == 'Left':
                        if type_ == 'Line':
                            self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                        else:
                            self.ax.plot(x_data, y, linestyle='None', marker=style, markersize=15.0,
                                         markeredgecolor="black",
                                         label="Legend")
                        self.ax.set_ylabel('$Y$ []')
                        self.ax.set_xlabel('$X$ []')
                    else:
                        if type_ == 'Line':
                            self.ax_right.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style)
                        else:
                            self.ax_right.plot(x_data, y, linestyle='None', marker=style, markersize=15.0,
                                               markeredgecolor="black",
                                               label="Legend")
                        self.ax_right.set_ylabel('$Y$ [dB]')
        else:
            print_("Please specify columns to plot")

    def plot_other(self, plot_dict, key):
        args = plot_dict
        requestX = args["plot inputs"]['Request'][1].currentText()
        requestY = args["plot inputs"]['Request'][2].currentText().split(', ')
        ids = [a.strip() for a in args["plot inputs"]['Id'].currentText().split(',')]  # get list
        filename = args["plot inputs"]['Folder'][0].text()
        # state = args["plot inputs"]['Toggle'].checkState()
        axis = args["plot inputs"]['Axis'].currentText()
        type_ = [args["plot inputs"]['Type'][0].currentText(), args["plot inputs"]['Type'][1].currentText()]
        style = args["plot inputs"]['Type'][1].currentText()
        filter_ = [args["plot inputs"]['Filter'][0].lineEdit().text(), args["plot inputs"]['Filter'][1].text()]
        value = args["plot inputs"]['Filter'][1].text()

        try:
            scaleX = float(args["plot inputs"]['ScaleX'].text())
            scaleY = float(args["plot inputs"]['ScaleY'].text())
        except ValueError:
            scaleX = 1
            scaleY = 1

        # if filename == '':
        #     folder = fr'{self.main_control.projectDir}/SimulationData/ABCI'

        # record data input
        args["plot data inputs"] = {"Ids": ids, "RequestX": requestX, "RequestY": requestY, "Folder": filename,
                                    "Axis": axis, "ScaleX": scaleX, "ScaleY": scaleY, "Type": type_, "Filter": filter_}

        for id_ in ids:
            args["plot data"][id_] = {}
            args["plot object"][id_] = {}
            if requestY != [] and requestX != []:
                if filename.split('.')[-1] == 'xlsx':
                    # get sheets
                    if filter_[0] == "None" or value == "":
                        self.other_data_filtered = self.other_data[key][id_]
                    else:
                        self.other_data_filtered = self.other_data[key][id_][
                            self.other_data[key][id_][filter_[0]] == value]

                    x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                    self.freq_glob = x_data
                    for j in range(len(requestY)):
                        y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                        args["plot data"][id_].update({j: {"x": x_data, "y": y}})
                        # ic('In here')
                        if axis == 'Left':
                            if type_[0] == 'Line':
                                args["plot object"][id_].update(
                                    {j: self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style,
                                                                                       picker=True)})
                            else:
                                args["plot object"][id_].update({j: self.ax.plot(x_data, y, linestyle='None',
                                                                                 marker=style, markersize=15.0,
                                                                                 markeredgecolor="black",
                                                                                 label=requestY[j], picker=True)})
                            # mplcursors.cursor(args["plot object"][id_][j])
                            self.ax.set_ylabel('$Y$ []')
                            self.ax.set_xlabel('$X$ []')
                        else:
                            if type_[0] == 'Line':
                                args["plot object"][id_].update({j: self.ax_right.plot(x_data, y, label=requestY[j],
                                                                                       linewidth=2,
                                                                                       linestyle=style,
                                                                                       picker=True)})
                            else:
                                args["plot object"][id_].update({j: self.ax_right.plot(x_data, y, linestyle='None',
                                                                                       marker=style,
                                                                                       markeredgecolor="black",
                                                                                       markersize=15.0,
                                                                                       label="Legend",
                                                                                       picker=True)})
                            # mplcursors.cursor(args["plot object"][id_][j])
                            self.ax_right.set_ylabel('$Y$ [dB]')
                elif filename.split('.')[-1].lower() == 'txt' or filename.split('.')[-1].lower() == 'csv':
                    # get sheets
                    if filter_[0] == "None" or value == "":
                        self.other_data_filtered = self.other_data[key][id_]
                    else:
                        self.other_data_filtered = self.other_data[key][id_][
                            self.other_data[key][id_][filter_[0]] == value]

                    x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                    self.freq_glob = x_data
                    for j in range(len(requestY)):
                        y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                        args["plot data"][id_].update({j: {"x": x_data, "y": y}})

                        if axis == 'Left':
                            if type_[0] == 'Line':
                                args["plot object"][id_].update(
                                    {j: self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style,
                                                     picker=True)})
                            else:
                                args["plot object"][id_].update({j: self.ax.plot(x_data, y, linestyle='None',
                                                                                 marker=style, markersize=15.0,
                                                                                 markeredgecolor="black",
                                                                                 label=requestY[j], picker=True)})
                            # mplcursors.cursor(args["plot object"][id_][j])
                            self.ax.set_ylabel('$Y$ []')
                            self.ax.set_xlabel('$X$ []')
                        else:
                            if type_[0] == 'Line':
                                args["plot object"][id_].update({j: self.ax_right.plot(x_data, y, label=requestY[j],
                                                                                       linewidth=2,
                                                                                       linestyle=style,
                                                                                       picker=True)})
                            else:
                                args["plot object"][id_].update({j: self.ax_right.plot(x_data, y, linestyle='None',
                                                                                       marker=style,
                                                                                       markeredgecolor="black",
                                                                                       markersize=15.0,
                                                                                       label="Legend",
                                                                                       picker=True)})
                            # mplcursors.cursor(args["plot object"][id_][j])
                            self.ax_right.set_ylabel('$Y$ [dB]')
            else:
                # try to filter self.other_data
                filter_ = [args["plot inputs"]['Filter'][0].lineEdit().text(), args["plot inputs"]['Filter'][1].text()]
                value = args["plot inputs"]['Filter'][1].text()
                if filter_[0] == "None" or value == "":
                    self.other_data_filtered = self.other_data[key]
                else:
                    self.other_data_filtered = self.other_data[key][self.other_data[key][filter_[0]] == value]
                x_data = [a * scaleX for a in self.other_data_filtered[requestX].tolist()]
                self.freq_glob = x_data

                for j in range(len(requestY)):
                    y = [a * scaleY for a in self.other_data_filtered[requestY[j]].tolist()]
                    args["plot data"][id_].update({j: {"x": x_data, "y": y}})

                    if axis == 'Left':
                        if type_[0] == 'Line':
                            args["plot object"][id_].update(
                                {j: self.ax.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style,
                                                 picker=True)})
                        else:
                            args["plot object"][id_].update(
                                {j: self.ax.plot(x_data, y, linestyle='None', marker=style, markersize=15.0,
                                                 markeredgecolor="black",
                                                 label="Legend", picker=True)})
                        # mplcursors.cursor(args["plot object"][id_][j])
                        self.ax.set_ylabel('$Y$ []')
                        self.ax.set_xlabel('$X$ []')
                    else:
                        if type_[0] == 'Line':
                            args["plot object"][id_].update(
                                {j: self.ax_right.plot(x_data, y, label=requestY[j], linewidth=2, linestyle=style,
                                                       picker=True)})
                        else:
                            args["plot object"][id_].update(
                                {j: self.ax_right.plot(x_data, y, linestyle='None', marker=style, markersize=15.0,
                                                       markeredgecolor="black",
                                                       label="Legend", picker=True)})
                        # mplcursors.cursor(args["plot object"][id_][j])
                        self.ax_right.set_ylabel('$Y$ [dB]')
            pass

    def reset_colors(self):
        self.ax.set_prop_cycle(None)
        # cycle, _ = itertools.tee(self.initial_prop_cycler)
        cycle = self.ax._get_lines.prop_cycler

        # update colors, loop over all lines and give color in ascending order
        for line2D in self.ax.get_lines():
            line2D.set_color(next(cycle)['color'])

        self.update_labels()
        self.fig.canvas.draw_idle()

    @staticmethod
    def plot_mat(filepath):
        # load mat file
        data = {}
        # files_folder = "D:\Dropbox\multipacting\MPGUI21"
        for f in filepath:
            if ".mat" in f:
                data[f] = spio.loadmat(fr"{filepath}")

        print(data)

    def table_control(self):
        # fill the first line
        self.tableWidget.create_new_row(0, self.tableWidget)

    def generateMenu(self, pos):
        # Get index
        menu = QMenu()
        item1 = menu.addAction("Add Row")
        item2 = menu.addAction("Delete Row")
        # Make the menu display in the normal position
        screenPos = self.tableWidget.mapToGlobal(pos)

        # Click on a menu item to return, making it blocked
        action = menu.exec(screenPos)
        if action == item1:
            self.tableWidget.add_row()
        if action == item2:
            self.tableWidget.remove_row(self.row)
        else:
            return

    def populate_plot_elements_table(self):
        le_Xlabel = QLineEdit("$f ~\mathrm{[MHz]}$")
        le_Xlabel.setObjectName('le_Xlabel')
        self.ui.tw_Plot_Elements.setCellWidget(0, 2, le_Xlabel)
        self.plot_elements_table['le_Xlabel'] = le_Xlabel

        sb_XLabel_Size = QSpinBox()
        sb_XLabel_Size.setObjectName('sb_XLabel_Size')
        self.ui.tw_Plot_Elements.setCellWidget(0, 0, sb_XLabel_Size)
        self.plot_elements_table['sb_XLabel_Size'] = sb_XLabel_Size

        sb_XLabel_Tick_Size = QSpinBox()
        sb_XLabel_Tick_Size.setObjectName('sb_XLabel_Tick_Siz')
        self.ui.tw_Plot_Elements.setCellWidget(0, 1, sb_XLabel_Tick_Size)
        self.plot_elements_table['sb_XLabel_Tick_Size'] = sb_XLabel_Tick_Size

        le_Ylabel = QLineEdit("$Z_{\parallel, \perp} ~[\mathrm{k\Omega}, \mathrm{k\Omega/m}]$")
        le_Ylabel.setObjectName('le_Ylabel')
        self.ui.tw_Plot_Elements.setCellWidget(1, 2, le_Ylabel)
        self.plot_elements_table['le_Ylabel'] = le_Ylabel

        sb_YLabel_Size = QSpinBox()
        sb_YLabel_Size.setObjectName('sb_YLabel_Size')
        self.ui.tw_Plot_Elements.setCellWidget(1, 0, sb_YLabel_Size)
        self.plot_elements_table['sb_YLabel_Size'] = sb_YLabel_Size

        sb_YLabel_Tick_Size = QSpinBox()
        sb_YLabel_Tick_Size.setObjectName('sb_YLabel_Tick_Size')
        self.ui.tw_Plot_Elements.setCellWidget(1, 1, sb_YLabel_Tick_Size)
        self.ui.tw_Plot_Elements.setColumnWidth(0, 150)
        self.plot_elements_table['sb_YLabel_Tick_Size'] = sb_YLabel_Tick_Size

        le_Legend = QLineEdit()
        le_Legend.setObjectName('le_Legend')
        self.ui.tw_Plot_Elements.setCellWidget(2, 2, le_Legend)
        self.plot_elements_table['le_Legend'] = le_Legend

        sb_Legend_Size = QSpinBox()
        sb_Legend_Size.setObjectName('sb_Legend_Size')
        self.ui.tw_Plot_Elements.setCellWidget(2, 0, sb_Legend_Size)
        self.plot_elements_table['sb_Legend_Size'] = sb_Legend_Size

        sb_Legend_Marker_Size = QSpinBox()
        sb_Legend_Marker_Size.setObjectName('sb_Legend_Marker_Size')
        self.ui.tw_Plot_Elements.setCellWidget(2, 1, sb_Legend_Marker_Size)
        self.ui.tw_Plot_Elements.setColumnWidth(1, 150)
        self.plot_elements_table['sb_Legend_Marker_Size'] = sb_Legend_Marker_Size

        le_Title = QLineEdit()
        le_Title.setObjectName('le_Title')
        self.ui.tw_Plot_Elements.setCellWidget(3, 2, le_Title)
        self.plot_elements_table['le_Title'] = le_Title

        sb_Title_Size = QSpinBox()
        sb_Title_Size.setObjectName('sb_Title_Size')
        self.ui.tw_Plot_Elements.setCellWidget(3, 0, sb_Title_Size)
        self.plot_elements_table['sb_Title_Size'] = sb_Title_Size

        # set completers
        self.set_completer()

    @staticmethod
    def switchRequest(cb_Request, request_dict_long, request_dict_trans, cb_Code):
        # clear combo box
        cb_Request.clear()
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
            row = self.tableWidget.currentRow()
            self.row = row

    def load_operating_points(self, filename):
        try:
            self.ui.le_Machine_Parameters_Filename.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            # populate checkboxes with key
            self.ui.ccb_Operation_Points.clear()
            self.ui.ccb_Operation_Points.addItem("All")
            for col in dd.keys():
                self.ui.ccb_Operation_Points.addItem(fr'{col}')

            self.operation_points = pd.DataFrame(dd)

        except Exception as e:
            print('Failed to open file:: ', e)

    def operation_points_selection(self, selection):
        cols = selection.split(', ')
        if not self.operation_points.empty:
            if selection is None or selection == '':
                pass
            else:
                # update machine parameters in View Machine(s)
                self.ui.le_E0.setText(f'{list(self.operation_points.loc["E [GeV]", cols])}')
                self.ui.le_Nu_S.setText(f'{list(self.operation_points.loc["nu_s []", cols])}')
                self.ui.le_I0.setText(f'{list(self.operation_points.loc["I0 [mA]", cols])}')
                self.ui.le_Alpha_P.setText(f'{list(self.operation_points.loc["alpha_p [1e-5]", cols])}')
                self.ui.le_Tau_Z.setText(f'{list(self.operation_points.loc["tau_z [ms]", cols])}')
                self.ui.le_Tau_XY.setText(f'{list(self.operation_points.loc["tau_xy [ms]", cols])}')
                self.ui.le_F_Rev.setText(f'{list(self.operation_points.loc["f_rev [kHz]", cols])}')
                self.ui.le_Beta_XY.setText(f'{list(self.operation_points.loc["beta_xy [m]", cols])}')
                self.ui.le_N_Cav.setText(f'{list(self.operation_points.loc["N_c []", cols])}')

    def calc_limits(self, mode, selection=None):
        if self.operation_points is not None:
            if len(self.freq_glob) > 0:
                if self.ui.cb_Longitudinal_Threshold.checkState() == 2 or \
                        self.ui.cb_Transverse_Threshold.checkState() == 2:
                    E0 = [45.6, 80, 120, 182.5]  # [GeV] Energy
                    nu_s = [0.025, 0.0506, 0.036, 0.087]  # Synchrotron oscillation tune
                    I0 = [1400, 135, 26.7, 5]  # [mA] Beam current5.4 * 2
                    alpha_c = [1.48, 1.48, 0.73, 0.73]  # [10−5] Momentum compaction factor
                    tau_z = [424.6, 78.7, 23.4, 6.8]  # [ms] Longitudinal damping time
                    tau_xy = [849.2, 157.4, 46.8, 13.6]  # [ms] Transverse damping time
                    f_rev = [3.07, 3.07, 3.07, 3.07]  # [kHz] Revolution frequency
                    beta_xy = 50
                    n_cav = [56, 112, 128, 140]  # 1_2_2_25

                    if selection is None or selection == '':
                        if self.operation_points.empty:
                            pass
                        else:
                            E0 = list(self.operation_points.loc["E [GeV]"])
                            nu_s = list(self.operation_points.loc["nu_s []"])
                            I0 = list(self.operation_points.loc["I0 [mA]"])
                            alpha_c = list(self.operation_points.loc["alpha_p [1e-5]"])
                            tau_z = list(self.operation_points.loc["tau_z [ms]"])
                            tau_xy = list(self.operation_points.loc["tau_xy [ms]"])
                            f_rev = list(self.operation_points.loc["f_rev [kHz]"])
                            beta_xy = list(self.operation_points.loc["beta_xy [m]"])
                            n_cav = list(self.operation_points.loc["N_c []"])
                    else:
                        cols = selection.split(', ')
                        E0 = list(self.operation_points.loc["E [GeV]", cols])
                        nu_s = list(self.operation_points.loc["nu_s []", cols])
                        I0 = list(self.operation_points.loc["I0 [mA]", cols])
                        alpha_c = list(self.operation_points.loc["alpha_p [1e-5]", cols])
                        tau_z = list(self.operation_points.loc["tau_z [ms]", cols])
                        tau_xy = list(self.operation_points.loc["tau_xy [ms]", cols])
                        f_rev = list(self.operation_points.loc["f_rev [kHz]", cols])
                        beta_xy = list(self.operation_points.loc["beta_xy [m]", cols])
                        n_cav = list(self.operation_points.loc["N_c []", cols])

                    unit = {'MHz': 1e6,
                            'GHz': 1e9}

                    Z_list, ZT_le = [], []
                    if mode == 'monopole':
                        # trim f
                        # f_list = f_list[0:len(f_list) - 100]
                        # f_list = self.freq_glob[0:len(self.freq_glob)]

                        # f_list = np.linspace(0, self.freq_glob[-1], num=1000)

                        f_list = np.linspace(self.ui.dsb_Frange_Min.value(),
                                             self.ui.dsb_Frange_Max.value(),
                                             num=1000) * unit[self.ui.cb_Frange_Unit.currentText()]
                        Z_le = []
                        try:
                            for i, n in enumerate(n_cav):
                                Z = [(2 * E0[i] * 1e9 * nu_s[i])
                                     / (n * I0[i] * 1e-3 * alpha_c[i] * 1e-5 * tau_z[i] * 1e-3 * f)
                                     if f > 1e-8 else 1e5 for f in f_list]

                                Z_le.append(round((2 * E0[i] * 1e9 * nu_s[i])
                                                  / (n * I0[i] * 1e-3 * alpha_c[i] * 1e-5
                                                     * tau_z[i] * 1e-3) * 1e-9 * 1e-3, 2))
                                Z_list.append(Z)

                            # note that plotted values are normalised to the number of cavities required per module
                            # but the values printed on the line editor are for one cavity
                            self.plot_baselines(f_list / unit[self.ui.cb_Frange_Unit.currentText()],
                                                self.ui.sb_Ncav_per_Mod.value() * np.array(Z_list) * 1e-3,
                                                mode, labels=selection.split(', '))  # to kOhm

                            self.ui.le_Zth_Long.setText(f"[{Z_le}]")
                        except ZeroDivisionError:
                            print("ZeroDivisionError, check input")
                        return f_list, Z_list

                    elif mode == 'dipole':
                        # f_list = self.freq_glob[0:len(self.freq_glob)]
                        # f_list = np.linspace(0, self.freq_glob[-1], num=1000)

                        f_list = np.linspace(self.ui.dsb_Frange_Min.value(),
                                             self.ui.dsb_Frange_Max.value(),
                                             num=1000) * unit[self.ui.cb_Frange_Unit.currentText()]
                        try:
                            for i, n in enumerate(n_cav):
                                ZT = (2 * E0[i]) * 1e9 / (
                                        n * I0[i] * 1e-3 * beta_xy[i] * tau_xy[i] * 1e-3 * f_rev[i] * 1e3)
                                ZT_le.append(round(ZT * 1e-3, 2))

                                # note that plotted values are normalised to the number of cavities required per module
                                # but the values printed on the line editor are for one cavity
                                Z_list.append(ZT)
                            self.ui.le_Zth_Trans.setText(f"[{ZT_le}]")
                            self.plot_baselines(f_list / unit[self.ui.cb_Frange_Unit.currentText()],
                                                self.ui.sb_Ncav_per_Mod.value() * np.array(Z_list) * 1e-3,
                                                mode, labels=selection.split(', '))  # to kOhm

                        except ZeroDivisionError:
                            print("ZeroDivisionError, check input")

                        return Z_list

                else:
                    self.remove_baselines()
        else:
            print("Please load a valid operation point(s) file.")

    def plot_baselines(self, f_list, Z_list, mode, labels):
        if mode == 'monopole':
            # plot baselines
            for i, z in enumerate(Z_list):
                aa = self.ax.plot(f_list, z, ls='--', c='k')

                # transform axes coordinates to data coordinates

                # axis_to_data = self.ax.transLimits
                # pos = axis_to_data.transform((0.05, 0.5))
                # indx = np.argmin(abs(f_list - pos[0]))
                # inv = axis_to_data.inverted()
                # print(pos)
                # print((f_list[indx], z[indx]))
                # x, y = inv.transform((f_list[indx], z[indx]))
                # print(x, y)

                pos = self.plt.axis_data_coords_sys_transform(self.ax, 0.01, 0.5)
                # pos2 = self.axis_data_coords_sys_transform(self.ax, pos[0], pos[1], True)
                indx = np.argmin(abs(f_list - pos[0]))
                x, y = self.plt.axis_data_coords_sys_transform(self.ax, f_list[indx], z[indx], True)

                lab = labels[i].split('_')
                if len(lab) > 1:
                    txt = r"$\mathrm{" + fr"{labels[i].split('_')[0]}_" + "\mathrm{" + fr"{labels[i].split('_')[1]}" + r"}}$"
                else:
                    txt = r"$\mathrm{" + fr"{labels[i]}" + r"}$"

                ab = self.plt.add_text(txt, box="Square", xy=(x, y), xycoords='axes fraction', size=14)
                # ab = self.ax.text(1, z[200], f'{text[i]}')

                # keep record
                self.baseline_line_objects.append(aa[0])
                self.baseline_line_objects.append(ab)
        else:
            # plot baselines
            for i, z in enumerate(Z_list):
                aa = self.ax.axhline(z, ls='--', c='k')

                pos = self.plt.axis_data_coords_sys_transform(self.ax, 0.01, 0.5)
                # pos2 = self.axis_data_coords_sys_transform(self.ax, pos[0], pos[1], True)
                indx = np.argmin(abs(f_list - pos[0]))
                x, y = self.plt.axis_data_coords_sys_transform(self.ax, f_list[indx], z, True)

                lab = labels[i].split('_')
                if len(lab) > 1:
                    txt = r"$\mathrm{" + fr"{labels[i].split('_')[0]}_" + "\mathrm{" + fr"{labels[i].split('_')[1]}" + r"}}$"
                else:
                    txt = r"$\mathrm{" + fr"{labels[i]}" + r"}$"

                ab = self.plt.add_text(txt, box="Square", xy=(x, y), xycoords='axes fraction', size=14)

                # keep record
                self.baseline_line_objects.append(aa)
                self.baseline_line_objects.append(ab)

        self.ax.autoscale(True, axis='y')
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def remove_baselines(self):
        for line in self.baseline_line_objects:
            if line in self.ax.get_lines():
                line.remove()

            if line in self.ax.findobj():
                line.remove()

        self.ax.autoscale(True, axis='y')
        self.ax.relim()
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def populate_cutoff_combobox(self, cutoff):
        cutoff.clear()
        cutoff.addItem("All")

        ll = text_to_list__(self.ui.le_Ri_Cutoff_List.text())
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

        if len(Ri_list) != 0:
            cutoff.addItems(self.mode_list_sorted[f"{Ri_list[0]}"])

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

                vl = self.ax.axvline(freq, ls='--', c='k')  # label=f"{sc[0]} cutoff (Ri={sc[1]})",

                # get y text position from axis position. Position x is not used
                pos = self.plt.axis_data_coords_sys_transform(self.ax, freq, 0.05, inverse=False)

                # ylim = self.ax.get_ylim()
                ab = self.plt.add_text(r"$f_\mathrm{c," + f"{sc[0]}" + r"} (R_\mathrm{i} = "
                                       + f"{sc[1]}" + r" ~\mathrm{mm}) $",
                                       box="None", xy=(freq, pos[1]),
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
        self.fig.canvas.blit()

    @staticmethod
    def calc_cutoff(Ri, mode):
        # calculate frequency from Ri
        p_TM01, p_TE11 = 2.405, 1.841
        c = 299792458  # m/s

        if mode == 'TM01':
            freq = (c * p_TM01) / (2 * np.pi * Ri * 1e9) * 1e3
        else:
            freq = (c * p_TE11) / (2 * np.pi * Ri * 1e9) * 1e3

        return freq

    def toggle_page(self, key):

        if key == 'Plot Decorations':
            self.ui.sw_Plot_Area_Tools.setCurrentIndex(0)

        if key == 'Machine Parameters':
            self.ui.sw_Plot_Area_Tools.setCurrentIndex(1)

        if key == 'Plot Element Format':
            self.ui.sw_Plot_Area_Tools.setCurrentIndex(2)

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

    def open_(self, le, cb_pol=None, ccb=None, mode='Folder', start_dir=''):
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
                                                      "Excel Files (*.txt *.csv *.xlsx *.json)")
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
                if txt in ccb.texts:
                    item = ccb.model().item(ccb.texts.index(txt))
                    item.setCheckState(2)

    def set_table_size(self):
        self.tableWidget.setColumnWidth(0, 50)
        self.tableWidget.setColumnWidth(1, 75)
        self.tableWidget.setColumnWidth(2, 200)
        self.tableWidget.setColumnWidth(3, 100)
        self.tableWidget.setColumnWidth(4, 300)
        self.tableWidget.setColumnWidth(5, 200)
        self.tableWidget.setColumnWidth(6, 75)
        self.tableWidget.setColumnWidth(7, 75)
        self.tableWidget.setColumnWidth(8, 75)
        self.tableWidget.setColumnWidth(10, 100)

    @staticmethod
    def validating(le):
        validation_rule = QDoubleValidator(0, 1e12, 10)
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
            # df = fr.txt_reader(filename, "\t")
            df = pd.read_csv(filename, header=None, sep='\s+')
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
                cb.addItem("All")

                for a in dd.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1].lower() == 'txt':
            self.plot_dict[i]["plot inputs"]['Id'].addItem('TXT')
            self.plot_dict[i]["plot inputs"]['Id'].setCurrentText('TXT')

            # load file txt
            # df = fr.txt_reader(filename, "\t")
            df = pd.read_csv(filename, header=None, sep='\s+')
            # rename columns to string
            df.columns = [f'{cname}' for cname in df.columns]
            self.other_data[i] = {'TXT': df}

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()
                cb.addItem("All")

                for a in df.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1].lower() == 'csv':
            self.plot_dict[i]["plot inputs"]['Id'].addItem('CSV')
            self.plot_dict[i]["plot inputs"]['Id'].setCurrentText('CSV')

            delimeter = ','
            # detect delimeter
            with open(filename) as f:
                firstline = f.readline()
                delim_comma_count = firstline.count(',')
                delim_semi_colon_count = firstline.count(';')

            if delim_comma_count == 0 and delim_semi_colon_count == 0:
                delimeter = '\s+'
            elif delim_semi_colon_count > 0:
                delimeter = ';'

            # load file txt
            # df = fr.txt_reader(filename, "\t")
            df = pd.read_csv(filename, sep=delimeter)
            # rename columns to string
            df.columns = [f'{cname}' for cname in df.columns]
            self.other_data[i] = {'CSV': df}

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()
                cb.addItem("All")

                for a in df.keys():
                    cb.addItem(f"{a}")

        if filename.split('.')[-1] == 'json':
            # load file txt
            # df = fr.json_reader(filename)
            df = pd.read_json(filename)
            self.other_data[i] = df

            for cb in cb_list:
                # clear combobox item if any
                cb.clear()
                cb.addItem("All")
                for a in df.keys():
                    cb.addItem(f"{a}")

    @staticmethod
    def populate_combobox_list(cb, ll):
        cb.clear()
        for a in ll:
            cb.addItem(a)

    def update_labels(self):
        # select axes to update
        if "left" in self.ui.cb_Active_Axis.currentText().lower():
            ax_current = self.ax
            ax_other = self.ax_right
        else:
            ax_current = self.ax_right
            ax_other = self.ax

        xlabel = self.plot_elements_table['le_Xlabel'].text()
        ylabel = self.plot_elements_table['le_Ylabel'].text()
        title = self.plot_elements_table['le_Title'].text()
        xsize = self.plot_elements_table['sb_XLabel_Size'].value()
        ysize = self.plot_elements_table['sb_YLabel_Size'].value()
        title_size = self.plot_elements_table['sb_Title_Size'].value()
        xtick_size = self.plot_elements_table['sb_XLabel_Tick_Size'].value()
        ytick_size = self.plot_elements_table['sb_YLabel_Tick_Size'].value()
        legend_size = self.plot_elements_table['sb_Legend_Size'].value()

        legend_labels = self.plot_elements_table['le_Legend'].text().split("%%")

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

        self.leg = ax_current.legend(lines + lines2, labels, edgecolor='k', facecolor='none', fontsize=legend_size)

        self.leg.set_zorder(10)
        self.leg.set_draggable(state=True, use_blit=True)

        self.fig.canvas.draw_idle()

    @staticmethod
    def get_line_properties(line):
        attr_dict = {'ls': line.get_linestyle(), 'color': line.get_color(), 'lw': line.get_linewidth(),
                     'marker': line.get_marker(),
                     'ms': line.get_markersize(), 'mec': line.get_markeredgecolor(), 'mew': line.get_markeredgewidth(),
                     'mfc': line.get_markerfacecolor()}

        return attr_dict

    @staticmethod
    def get_annotation_properties(text):
        if text.get_bbox_patch():
            box = 'None'
            boxstyle_obj = text.get_bbox_patch().get_boxstyle()
            if isinstance(boxstyle_obj, matplotlib.patches.BoxStyle.Round):
                 box = 'Round'
            elif isinstance(boxstyle_obj, matplotlib.patches.BoxStyle.Square):
                box = 'Square'
            elif isinstance(boxstyle_obj, matplotlib.patches.BoxStyle.Round4):
                box = 'Round4'
            elif isinstance(boxstyle_obj, matplotlib.patches.BoxStyle.Ellipse):
                box = 'Ellipse'

            attr_dict = {'text': text.get_text(), 'color': text.get_color(), 'style': text.get_style(),
                         'fontsize': text.get_fontsize(), 'rotation': text.get_rotation(), 'position': text.get_position(),
                         'fontname': text.get_fontname(),
                         'box': box,
                         'usetex': text.get_usetex()}
        else:
            attr_dict = {'text': text.get_text(), 'color': text.get_color(), 'style': text.get_style(),
                         'fontsize': text.get_fontsize(), 'rotation': text.get_rotation(), 'position': text.get_position(),
                         'fontname': text.get_fontname(),
                         # 'font': text.get_font(), 'bbox_patch': text.get_bbox_patch(),
                         'usetex': text.get_usetex()}
        # print(text.get_bbox_patch())
        # print(text.get_bbox_patch().get_width())
        # print(text.get_bbox_patch().get_boxstyle().stylename)
        # print(text.get_bbox_patch().get_facecolor())
        # print(text.get_bbox_patch().get_edgecolor())
        return attr_dict

    @staticmethod
    def set_line_properties(line, attr_dict):
        line.set(attr_dict)

    def get_stylesheet(self):
        return self.w_Plot.styleSheet()

    def set_latex_text_color(self):
        if self.main_ui.cb_Theme.currentText() not in ['Light VS', 'Solarized Light']:
            color = 'white'
        else:
            color = 'k'

        for i in range(1, 24):
            eval(f"self.ui.latexLabel_{i}.set_color('{color}')")

    def set_completer(self):
        # set XLabel completer
        names = [r"$f ~[\mathrm{MHz}]$"]
        completer = QCompleter(names)
        self.plot_elements_table['le_Xlabel'].setCompleter(completer)

        # set YLabel completer
        names = [r"$Z_{\parallel, \perp} ~[\mathrm{k\Omega}, \mathrm{k\Omega/m}]$",
                 r"$Z_{\parallel} ~[\mathrm{k\Omega}]$",
                 r"$Z_{\perp} ~[\mathrm{k\Omega/m}]$", r"$S~\mathrm{[dB]}$", r"$Q_\mathrm{ext}$"]

        completer = QCompleter(names)
        self.plot_elements_table['le_Ylabel'].setCompleter(completer)

        # set legend completer
        names = [r"Monopole",
                 r"Monopole%%Dipole",
                 r"Monopole%%Dipole%%Quadrupole",
                 r"Monopole%%Dipole%%Quadrupole%%Sextupole",
                 r"Monopole%%Dipole%%Quadrupole%%Sextupole%%Octupole",
                 r"Monopole%%Dipole%%Quadrupole%%Sextupole%%Octupole%%Decapole",
                 r"$Z_{\parallel} ~[\mathrm{k\Omega}]$",
                 r"$Z_{\perp} ~[\mathrm{k\Omega/m}]$",
                 r"$S~\mathrm{[dB]}$",
                 r"$S_\mathrm{21}-\alpha$%%$S_\mathrm{21}-\beta$",
                 r"C3794-2HC-1FPC%%C3794-3HC(1)-1FPC%%C3794-3HC(2)-1FPC%%C3794-4HC(1)-1FPC%%C3794-4HC(2)-1FPC"]

        completer = QCompleter(names)
        self.plot_elements_table['le_Legend'].setCompleter(completer)

    def serialise(self, state_dict):
        serialise(state_dict, self.w_Plot, marker='plot')
        
        # serialise plot objects table widget
        for k, v in self.plot_elements_table.items():
            if v.objectName() != "" and 'qt_' not in v.objectName():
                if isinstance(v, QPushButton):
                    if v.isCheckable():
                        state_dict['plot_' + v.objectName()] = v.isChecked()

                if isinstance(v, QComboBox):
                    state_dict['plot_' + v.objectName()] = v.currentText()

                if isinstance(v, QLineEdit):
                    state_dict['plot_' + v.objectName()] = v.text()

                if isinstance(v, QSpinBox):
                    state_dict['plot_' + v.objectName()] = v.value()

                if isinstance(v, QDoubleSpinBox):
                    state_dict['plot_' + v.objectName()] = v.value()

                if isinstance(v, QCheckBox):
                    state_dict['plot_' + v.objectName()] = v.checkState()

        # serialize table widget
        table_widget_state = {}
        for k, v in self.plot_dict.items():
            table_widget_state.update(
                {k:
                    {
                        "Code": v['plot inputs']["Code"].currentText(),
                        "Folder": v['plot inputs']["Folder"][0].text(),
                        "Polarization": v['plot inputs']["Polarization"].currentText(),
                        "Id": [[v['plot inputs']["Id"].itemText(i) for i in range(v['plot inputs']["Id"].count())],
                               v['plot inputs']["Id"].currentText()],
                        "Request": [v['plot inputs']["Request"][0].currentText(),
                                    v['plot inputs']["Request"][1].currentText(),
                                    v['plot inputs']["Request"][2].currentText()],
                        "Toggle": v['plot inputs']["Toggle"].checkState(),
                        "ScaleX": v['plot inputs']["ScaleX"].text(),
                        "ScaleY": v['plot inputs']["ScaleY"].text(),
                        "Axis": v['plot inputs']["Axis"].currentText(),
                        "Type": [v['plot inputs']["Type"][0].currentText(), v['plot inputs']["Type"][1].currentText()],
                        "Filter": [v['plot inputs']["Filter"][0].lineEdit().text(),
                                   v['plot inputs']["Filter"][1].text()]
                    }}
            )

        state_dict['plot_table_widget'] = table_widget_state

        # print([child for child in ax.get_children() if isinstance(child, matplotlib.text.Annotation)])
        plot_objects_attr = {"lines": [self.get_line_properties(line) for line in self.ax.get_lines()],
                             "collections": [self.get_line_properties(coll) for coll in self.ax.collections],
                             "annotations": [self.get_annotation_properties(text) for text in self.ax.get_children() if isinstance(text, matplotlib.text.Annotation)],
                             }

        state_dict['plot_objects_attr'] = plot_objects_attr

    def deserialise(self, state_dict):
        deserialise(state_dict, self.w_Plot, marker='plot')

        try:
            # deserialise plot element table
            for k, v in self.plot_elements_table.items():
                if v.objectName() != "" and 'qt_' not in v.objectName():
                    if isinstance(v, QPushButton):
                        if v.isCheckable():
                            if not v.isChecked() == state_dict['plot_' + v.objectName()]:
                                v.toggle()

                    if isinstance(v, QComboBox):
                        v.setCurrentText(state_dict['plot_' + v.objectName()])

                        # if isinstance(v, QCheckableComboBox):
                        try:
                            selection = state_dict['plot_' + v.objectName()].split(', ')
                            # try to select copied selection
                            for txt in selection:
                                if txt in v.texts:
                                    item = v.model().item(v.texts.index(txt))
                                    item.setCheckState(2)
                        except AttributeError:
                            pass

                    if isinstance(v, QLineEdit):
                        v.setText(state_dict['plot_' + v.objectName()])

                    if isinstance(v, QSpinBox):
                        v.setValue(int(state_dict['plot_' + v.objectName()]))

                    if isinstance(v, QDoubleSpinBox):
                        v.setValue(float(state_dict['plot_' + v.objectName()]))

                    if isinstance(v, QCheckBox):
                        v.setCheckState(state_dict['plot_' + v.objectName()])

            # serialise plot objects table widget

            # deserialise plots
            for k, v in self.plot_elements_table.items():
                if v.objectName() != "" and 'qt_' not in v.objectName():
                    if isinstance(v, QPushButton):
                        if v.isCheckable():
                            state_dict['plot_' + v.objectName()] = v.isChecked()

                    if isinstance(v, QComboBox):
                        state_dict['plot_' + v.objectName()] = v.currentText()

                    if isinstance(v, QLineEdit):
                        state_dict['plot_' + v.objectName()] = v.text()

                    if isinstance(v, QSpinBox):
                        state_dict['plot_' + v.objectName()] = v.value()

                    if isinstance(v, QDoubleSpinBox):
                        state_dict['plot_' + v.objectName()] = v.value()

                    if isinstance(v, QCheckBox):
                        state_dict['plot_' + v.objectName()] = v.checkState()

            table_widget_state = state_dict['plot_table_widget']
            # print(table_widget_state)
            # self.tableWidget.setRowCount(0)
            for i, (k, v) in enumerate(table_widget_state.items()):
                if i > 0:
                    self.tableWidget.add_row()

                args_dict = self.plot_dict[i]['plot inputs']

                args_dict['Code'].setCurrentText(v["Code"])
                args_dict['Folder'][0].setText(
                    v["Folder"])  # [le_Folder, pb_Open_Folder, w_Folder, l_Folder_Widget]

                args_dict['Polarization'].setCurrentText(v["Polarization"])

                # check checked items
                args_dict['Id'].clear()  # Not optimal
                for r, text in enumerate(v["Id"][0]):
                    # if r != 0:
                    args_dict['Id'].addItem(text)
                    if args_dict['Id'].model().item(r).text() in v["Id"][1].split(','):
                        args_dict['Id'].model().item(r).setCheckState(Qt.Checked)

                args_dict['Request'][0].setCurrentText(v["Request"][0])
                args_dict['Request'][1].setCurrentText(v["Request"][1])

                # check checked items
                for r in range(args_dict['Request'][2].model().rowCount()):
                    if args_dict['Request'][2].model().item(r).text() == v["Request"][2]:
                        args_dict['Request'][2].model().item(r).setCheckState(Qt.Checked)

                args_dict['Toggle'].setCheckState(v["Toggle"])
                args_dict['ScaleX'].setText(v["ScaleX"])
                args_dict['ScaleY'].setText(v["ScaleY"])
                args_dict['Axis'].setCurrentText(v["Axis"])
                args_dict['Type'][0].setCurrentText(v["Type"][0])
                args_dict['Type'][1].setCurrentText(v["Type"][1])

                # check checked items
                for r in range(args_dict['Filter'][0].model().rowCount()):
                    if args_dict['Filter'][0].model().item(r).text() in v['Filter'][0]:
                        args_dict['Filter'][0].model().item(r).setCheckState(Qt.Checked)

                args_dict['Filter'][1].setText(v["Filter"][1])

            self.plot()
            plot_objects_attr = state_dict['plot_objects_attr']
            try:
                for n, line in enumerate(self.ax.get_lines()):
                    line.set(**plot_objects_attr['lines'][n])
            except IndexError:
                print('Index error in deserialising line attributes.')
            try:
                for n, text_attr in enumerate(plot_objects_attr['annotations']):
                    if 'box' in text_attr.keys():
                        self.plt.add_text(text_attr['text'], text_attr['box'], xy=text_attr['position'], xycoords='axes fraction',
                                          xytext=None, textcoords='data', size=text_attr['fontsize'],
                                          rotation=text_attr['rotation'], arrowprops=None)
                    else:
                        self.plt.add_text(text_attr['text'], 'None', xy=text_attr['position'], xycoords='axes fraction',
                                          xytext=None, textcoords='data', size=text_attr['fontsize'],
                                          rotation=text_attr['rotation'], arrowprops=None)
            except IndexError:
                print('Index error in deserialising annotations attributes.')

            for n, coll in enumerate(self.ax.collections):
                coll.set(**plot_objects_attr['collections'][n])

            self.update_labels()
            self.fig.canvas.draw_idle()
        except (KeyError, TypeError) as e:
            print("Could not deserialize plot_control.py: ", e)


class TableWidgetDragDrop(QTableWidget):
    def __init__(self, parent):
        super().__init__()
        self.plot_control = parent

        self.setAcceptDrops(True)
        self.setDragEnabled(True)
        self.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.setDefaultDropAction(QtCore.Qt.TargetMoveAction)
        self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.setObjectName("tableWidget")
        self.setColumnCount(11)
        self.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(6, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(7, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(8, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(9, item)
        item = QtWidgets.QTableWidgetItem()
        self.setHorizontalHeaderItem(10, item)
        self.horizontalHeader().setStretchLastSection(True)

        _translate = QtCore.QCoreApplication.translate

        item = self.horizontalHeaderItem(1)
        item.setText(_translate("Plot", "Code"))
        item = self.horizontalHeaderItem(2)
        item.setText(_translate("Plot", "File/Folder"))
        item = self.horizontalHeaderItem(3)
        item.setText(_translate("Plot", "Polarisation"))
        item = self.horizontalHeaderItem(4)
        item.setText(_translate("Plot", "ID"))
        item = self.horizontalHeaderItem(5)
        item.setText(_translate("Plot", "Request"))
        item = self.horizontalHeaderItem(6)
        item.setText(_translate("Plot", "ScaleX"))
        item = self.horizontalHeaderItem(7)
        item.setText(_translate("Plot", "ScaleY"))
        item = self.horizontalHeaderItem(8)
        item.setText(_translate("Plot", "Axis"))
        item = self.horizontalHeaderItem(9)
        item.setText(_translate("Plot", "Type"))
        item = self.horizontalHeaderItem(10)
        item.setText(_translate("Plot", "Filter"))

    def dragEnterEvent(self, event: QtGui.QDragEnterEvent) -> None:
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event: QtGui.QDragMoveEvent) -> None:
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()
            
    def dropEvent(self, event: QtGui.QDropEvent) -> None:
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    filepath = str(url.toLocalFile())
                    if filepath.split('.')[-1].lower() in ['xlsx', 'txt', 'csv']:
                        self.add_row()
                        key = self.plot_control.plotID_count - 1
                        self.plot_control.plot_dict[key]["plot inputs"]['Code'].setCurrentText('Other')
                        self.plot_control.plot_dict[key]["plot inputs"]['Folder'][0].setText(filepath)

            event.accept()

    def add_row(self):
        # get current number of rows
        n = self.rowCount()

        # add new row
        self.setRowCount(n + 1)  # and one row in the table
        self.create_new_row(n, self)

    def create_new_row(self, row_ind, table):
        """

        Parameters
        ----------
        row_ind
        table

        Returns
        -------

        """

        sizePolicy = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        # Toggle on/off
        w_Toggle_Close = QWidget()
        l_Toggle_Close_Widget = QHBoxLayout()
        w_Toggle_Close.setLayout(l_Toggle_Close_Widget)
        table.setCellWidget(row_ind, 0, w_Toggle_Close)

        cb_toggle = QCheckBox()
        cb_toggle.setCheckState(Qt.Checked)
        cb_toggle.stateChanged.connect(lambda: self.plot_control.plot())

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
        w_Folder.setContentsMargins(0, 0, 0, 0)

        l_Folder_Widget = QHBoxLayout()
        w_Folder.setLayout(l_Folder_Widget)
        le_Folder = QLineEdit()
        le_Folder.setReadOnly(True)
        le_Folder.setSizePolicy(sizePolicy)
        pb_Open_Folder = QPushButton('...')
        pb_Open_Folder.setMaximumWidth(100)
        pb_Open_Folder.setMinimumWidth(50)
        pb_Open_Folder.setSizePolicy(sizePolicy)
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
            lambda: self.plot_control.open_(le_Folder, cb_pol, ccb_Id, 'File',
                               start_dir=f"{self.plot_control.main_control.projectDir}")
            if cb_code.currentText() == 'Other' else
            self.plot_control.open_(le_Folder, cb_pol, ccb_Id, 'Folder',
                       start_dir=f"{self.plot_control.main_control.projectDir}/SimulationData/{cb_code.currentText()}"))
        # le_Folder.textChanged.connect(lambda: )

        cb_pol.currentIndexChanged.connect(
            lambda: self.plot_control.populate_IDs(le_Folder.text(), cb_pol,
                                      ccb_Id) if cb_code.currentText() == 'Other' else self.plot_control.populate_IDs(
                le_Folder.text(), cb_pol, ccb_Id))

        # Request
        # Request widget
        w_Request = QWidget()
        l_Request_Widget = QHBoxLayout()
        l_Request_Widget.setContentsMargins(0, 0, 0, 0)
        l_Request_Widget.setSpacing(0)
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
        cb_request.setSizePolicy(sizePolicy)
        for req in request_dict_long.values():
            cb_request.addItem(req)
        # cb_request.setStyleSheet('background-color: yellow;')

        # create widget with two comboboxes
        cb_X = QComboBox()
        cb_Y = QCheckableComboBox()
        cb_X.setSizePolicy(sizePolicy)
        cb_Y.setSizePolicy(sizePolicy)

        # place cb_request and le_request in request widget layout
        l_Request_Widget.addWidget(cb_request)
        l_Request_Widget.addWidget(cb_X)
        l_Request_Widget.addWidget(cb_Y)

        # hide le_request
        cb_X.hide()
        cb_Y.hide()

        # add widget to table
        # w_Request.setStyleSheet('background-color: blue;')
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
        ccb_Filter.setSizePolicy(sizePolicy)
        ccb_Filter.addItem("All")
        l_Filter_Widget.addWidget(ccb_Filter)

        # line edit
        le_Value = QLineEdit()
        l_Filter_Widget.addWidget(le_Value)
        le_Value.setSizePolicy(sizePolicy)

        table.setCellWidget(row_ind, 10, w_Filter)

        # add signal to le_Folder after editing to update cb_X and cb_Y
        key = self.plot_control.plotID_count
        le_Folder.textChanged.connect(lambda: self.plot_control.populate_combobox([cb_X, cb_Y, ccb_Filter], le_Folder.text(), key))

        # signal to disable polarisation when cb_code item is 'Other'
        cb_code.currentIndexChanged.connect(lambda: self.plot_control.show_hide_request_widgets(cb_code, cb_request, cb_X, cb_Y))

        # add signal to switch between longitudinal and transverse
        cb_pol.currentIndexChanged.connect(
            lambda: self.plot_control.switchRequest(cb_request, request_dict_long, request_dict_trans, cb_pol))

        # ScaleX
        le_ScaleX = QLineEdit()
        le_ScaleX.setText('1')
        table.setCellWidget(row_ind, 6, le_ScaleX)

        # ScaleY
        le_ScaleY = QLineEdit()
        le_ScaleY.setText('1')
        table.setCellWidget(row_ind, 7, le_ScaleY)
        # le_ScaleX.editingFinished.connect(lambda: validating(le_ScaleX, default='1'))
        # le_ScaleY.editingFinished.connect(lambda: validating(le_ScaleY, default='1'))

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
        cb_type.setSizePolicy(sizePolicy)
        for a in type_list:
            cb_type.addItem(a)

        l_Type_Widget.addWidget(cb_type)
        # check box
        cb_style = QComboBox()
        cb_style.setSizePolicy(sizePolicy)
        l_Type_Widget.addWidget(cb_style)

        line_style = ['-', '--', ':', '-.']
        scatter_marker = ['x', 'o', '+', 'P', 'X', 'd', '>', '^', 'v', 'H']

        # fill default
        for a in line_style:
            cb_style.addItem(a)

        # signal for combobox
        cb_type.currentIndexChanged.connect(
            lambda: self.plot_control.populate_combobox_list(cb_style, line_style)
            if cb_type.currentText() == "Line"
            else self.plot_control.populate_combobox_list(cb_style, scatter_marker))

        table.setCellWidget(row_ind, 9, w_Type)

        # pb_Color = QPushButton()
        # menu = QMenu()
        # menu.triggered.connect(lambda x: pb_Color.setStyleSheet(f"background-color: {x.defaultWidget().text()};"))
        # menu.triggered.connect(lambda x: pb_Color.setText(x.defaultWidget().text()))
        # # menu.triggered.connect(lambda x: print(x.defaultWidget().text()))
        # pb_Color.setMenu(menu)
        #
        # for k, vals in MATPLOTLIB_COLORS.items():
        #     sub_menu = menu.addMenu(k)
        #     for v in vals:
        #         l = QLabel(str(v))
        #         action = QWidgetAction(sub_menu)
        #         l.setStyleSheet(f"background-color: {v};")
        #         action.setDefaultWidget(l)
        #         sub_menu.addAction(action)
        #         # action.triggered.connect(lambda: pb_Color.setStyleSheet(f"background-color: {v};"))
        #         # action.triggered.connect(lambda: pb_Color.setText(str(v)))
        #
        # table.setCellWidget(row_ind, 10, pb_Color)

        args_dict = {'Code': cb_code, 'Folder': [le_Folder, pb_Open_Folder, w_Folder, l_Folder_Widget],
                     'Polarization': cb_pol, 'Id': ccb_Id,
                     'Request': [cb_request, cb_X, cb_Y, w_Request, l_Request_Widget], 'Toggle': cb_toggle,
                     'Remove': pb_Delete_Row, 'ScaleX': le_ScaleX, 'ScaleY': le_ScaleY, 'Axis': cb_axis,
                     'Type': [cb_type, cb_style], 'Filter': [ccb_Filter, le_Value]}

        self.plot_control.plot_dict[key] = {"plot inputs": None, "plot data": {}, "plot object": {},
                               "plot data inputs": None}

        self.plot_control.plot_dict[key]["plot inputs"] = args_dict
        self.plot_control.plotID_count += 1

    def remove_row(self, pb_Remove):

        # remove from table
        n = self.rowCount()
        row = 0
        key = 0
        code = "ABCI"
        for i, (k, val) in enumerate(self.plot_control.plot_dict.items()):
            if val["plot inputs"]['Remove'] == pb_Remove:
                row = i
                key = k
                code = val["plot inputs"]['Code'].currentText()

        self.removeRow(row)
        # reset number of rows
        self.setRowCount(n - 1)

        if code == "ABCI":
            for line2D in self.plot_control.plot_dict[key]['plot object'].values():
                line2D[0].remove()
        else:
            for vals in self.plot_control.plot_dict[key]['plot object'].values():
                for line2D in vals.values():
                    line2D[0].remove()

        # reset plot data, plot object and plot data inputs
        self.plot_control.plot_dict[key]["plot data"] = {}
        self.plot_control.plot_dict[key]["plot object"] = {}
        self.plot_control.plot_dict[key]["plot data inputs"] = {}

        self.plot_control.plot_dict[key] = {}
        del self.plot_control.plot_dict[key]

        self.plot_control.plot()


