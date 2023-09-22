from PyQt5.QtCore import qVersion, QSize
from PyQt5.QtGui import *
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from frame_controls.plot_properties import PPropsControl


class MyCustomToolbar(NavigationToolbar):
    def __init__(self, plotCanvas, parent, coordinates=True):

        self.canvas = plotCanvas
        self.ax = self.canvas.ax
        # self.setContentsMargins(0, 0, 0, 0)

        # propertiesl widget
        self.pprops = PPropsControl(self.canvas)

        # create the default toolbar
        NavigationToolbar.__init__(self, plotCanvas, parent, coordinates)

        toolitems = [*NavigationToolbar.toolitems]
        # print(toolitems)

        # check toolitems
        toolitems = [('Home', 'Reset original view', 'home', 'home'),
                     ('Back', 'Back to previous view', 'back', 'back'),
                     ('Forward', 'Forward to next view', 'forward', 'forward'),
                     (None, None, None, None),
                     ('Pan', 'Left button pans, Right button zooms\nx/y fixes axis, CTRL fixes aspect', 'pan', 'pan'),
                     ('Zoom', 'Zoom to rectangle\nx/y fixes axis, CTRL fixes aspect', 'zoom', 'zoom'),
                     ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
                     ('Customize', 'Edit axis, curve and image parameters', 'options', 'edit_parameters'),
                     ('Plot format', 'Edit rcParams, curve and image parameters', 'options', 'rcParams'),
                     (None, None, None, None),
                     ('Linear', 'Linear scale', 'lin', 'linear_scale'),
                     ('Log', 'Log scale', 'log', 'log_scale'),
                     ('Decibel', 'Decibel', 'dB', 'decibel_scale'),
                     (None, None, None, None),
                     ('Grid', 'Toggle Grid', 'grid', 'toggle_grid'),
                     ('Transparency', 'Toggle Transparency', 'transparency', 'toggle_transparency'),
                     (None, None, None, None),
                     ("Size", None, None, None),
                     ('Save', 'Save the figure', 'save', 'save_figure'),
                     ]

        self.canvas.toolbar.setContentsMargins(0, 0, 0, 0)
        self.canvas.toolbar.setIconSize(QSize(20, 20))
        # print(type(self.canvas.toolbar))

        actions = self.findChildren(QAction)
        for a in actions:
            self.removeAction(a)

        self._actions = {}  # mapping of toolitem method names to QActions.

        for text, tooltip_text, image_file, callback in toolitems:
            if text is None:
                self.addSeparator()
            elif text == 'Size':
                # size definitions
                plot_settings_options = ['Presentation', 'Wakefield', 'Square', 'Portrait', 'Landscape', 'Custom']
                self.cb_Save_Plot_Settings = QComboBox()
                for a in plot_settings_options:
                    self.cb_Save_Plot_Settings.addItem(a)
                self.addWidget(self.cb_Save_Plot_Settings)
            elif text == 'Grid':
                self.pb_Grid = QPushButton()
                self.pb_Grid.setCheckable(True)
                self.pb_Grid.setToolTip('Toggle grid')
                self.pb_Grid.clicked.connect(lambda: self.toggle_grid())
                self.addWidget(self.pb_Grid)
            elif text == 'Transparency':
                self.pb_Transparency = QPushButton()
                self.pb_Transparency.setCheckable(True)
                self.pb_Transparency.setToolTip('Toggle transparency')
                self.pb_Transparency.clicked.connect(lambda: self.toggle_transparency())
                self.addWidget(self.pb_Transparency)

            else:
                a = self.addAction(self._icon(image_file + '.png'),
                                   text, getattr(self, callback))
                self._actions[callback] = a

                if callback in ['zoom', 'pan']:
                    a.setCheckable(True)
                if tooltip_text is not None:
                    a.setToolTip(tooltip_text)
        # Add the (x, y) location widget at the right side of the toolbar
        # The stretch factor is 1 which means any resizing of the toolbar
        # will resize this label instead of the buttons.
        if coordinates:
            self.locLabel = QLabel("", self)
            self.locLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.locLabel.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Ignored))
            labelAction = self.addWidget(self.locLabel)
            labelAction.setVisible(True)

    def save_figure(self, *args):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(None, "QFileDialog.getSaveFileName()", "", "PNG Files (*.png)", options=options)
        if fileName:
            if fileName.split('.')[-1] != 'png':
                fileName = f'{fileName}.png'

            try:
                if self.cb_Save_Plot_Settings.currentText() == 'Presentation':
                    # get size of figure before change
                    fig_size = self.canvas.fig.get_size_inches()
                    self.canvas.fig.set_size_inches(8, 6)  # actual size = 8*dpi, 6*dpi; dpi = 100
                elif self.cb_Save_Plot_Settings.currentText() == 'Wakefield':
                    # get size of figure before change
                    fig_size = self.canvas.fig.get_size_inches()
                    self.canvas.fig.set_size_inches(10, 5)  # forward=False
                elif self.cb_Save_Plot_Settings.currentText() == 'Square':
                    # get size of figure before change
                    fig_size = self.canvas.fig.get_size_inches()
                    self.canvas.fig.set_size_inches(8, 8)  # forward=False
                elif self.cb_Save_Plot_Settings.currentText() == 'Landscape':
                    # get size of figure before change
                    fig_size = self.canvas.fig.get_size_inches()
                    self.canvas.fig.set_size_inches(12, 8)  # forward=False
                else:
                    # get size of figure before change
                    fig_size = self.canvas.fig.get_size_inches()

                    # open popup window to get size input
                    size_input_dialog = SizeInputDialog(self)
                    size_input_dialog.show()
                    if size_input_dialog.exec_():
                        width, height = size_input_dialog.getInputs()
                        width, height = float(width)/self.canvas.fig.dpi, float(height)/self.canvas.fig.dpi

                        self.canvas.fig.set_size_inches(width, height)  # , forward=False

                self.canvas.set_axes_elements_color(color='k')

                # change spine, label and ticks to default black color before saving
                self.canvas.fig.savefig(fileName, transparent=True, bbox_inches='tight', pad_inches=0.2)

                # revert to display color
                self.canvas.set_axes_elements_color()

                # restore windows size
                self.canvas.fig.set_size_inches(fig_size[0], fig_size[1])  # actual size = 8*dpi, 6*dpi; dpi = 100
                self.canvas.draw()
                self.canvas.flush_events()

            except Exception as e:
                print("Cannot save file -> ", e)
        else:
            print('Enter a valid filename')

    def _icon(self, name):
        """
        Construct a `.QIcon` from an image file *name*, including the extension
        and relative to Matplotlib's "images" data directory.
        """
        # if qVersion() >= '5.':
        #     name = name.replace('.png', '_large.png')
        pm = QPixmap(f":/icons/icons/PNG/{name}")
        if self.palette().color(self.backgroundRole()).value() < 128:
            icon_color = self.palette().color(self.foregroundRole())
            mask = pm.createMaskFromColor(QColor('black'),
                                          Qt.MaskOutColor)
            pm.fill(icon_color)
            pm.setMask(mask)
        return QIcon(pm)

    def linear_scale(self):
        self.canvas.change_scale('linear')

    def log_scale(self):
        self.canvas.change_scale('log')

    def decibel_scale(self):
        pass

    def toggle_grid(self):
        if self.pb_Grid.isChecked():
            for ax in self.canvas.fig.get_axes():
                ax.grid(True, 'both')
        else:
            for ax in self.canvas.fig.get_axes():
                ax.grid(False)

        self.canvas.draw_idle()
        self.canvas.flush_events()

    def toggle_transparency(self):
        if self.pb_Transparency.isChecked():
            for ax in self.canvas.fig.get_axes():
                ax.set_facecolor('none')
            self.canvas.fig.set_facecolor('none')
        else:
            for ax in self.canvas.fig.get_axes():
                ax.set_facecolor('white')
            self.canvas.fig.set_facecolor('white')
        self.canvas.set_axes_elements_color()

    def rcParams(self):
        # pop up plot properties widget
        # create pop up widget instance
        self.pprops.ppropsUI.cb_Plots.addItems([x.get_label() for x in self.canvas.fig.get_axes()])

        # get and set default parameters
        # legend
        # update legend
        handles, labels = self.canvas.ax.get_legend_handles_labels()
        legend_labels = '%%'.join(labels)
        self.pprops.ppropsUI.le_Legend.setText(legend_labels)

        # axis limits
        self.pprops.ppropsUI.dsb_X0.setValue(0.45*(self.canvas.ax.get_xlim()[1]+self.canvas.ax.get_xlim()[0]))
        self.pprops.ppropsUI.dsb_X1.setValue(0.55*(self.canvas.ax.get_xlim()[1]+self.canvas.ax.get_xlim()[0]))
        self.pprops.ppropsUI.dsb_Y0.setValue(0.45*(self.canvas.ax.get_ylim()[1]+self.canvas.ax.get_ylim()[0]))
        self.pprops.ppropsUI.dsb_Y1.setValue(0.55*(self.canvas.ax.get_ylim()[1]+self.canvas.ax.get_ylim()[0]))

        self.pprops.w_PProps.show()


class SizeInputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.first = QLineEdit(self)
        self.second = QLineEdit(self)

        # add validator for only integers
        self.first.textChanged.connect(lambda: self.validating(self.first))
        self.second.textChanged.connect(lambda: self.validating(self.second))
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)

        layout = QFormLayout(self)
        layout.addRow("Width", self.first)
        layout.addRow("Height", self.second)
        layout.addWidget(buttonBox)

        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

    def validating(self, le):
        validation_rule = QDoubleValidator(1, 3840, 0)
        if validation_rule.validate(le.text(), 20)[0] == QValidator.Acceptable:
            le.setFocus()
        else:
            le.setText('')

    def getInputs(self):
        return self.first.text(), self.second.text()

    def button_press(self):
        if self.sender() == self.btn_ok:
            self.ok = True
        self.close()

    @classmethod
    def isOkay(cls, parent):
        dialog = cls(parent)
        dialog.exec_()
        return dialog.ok

