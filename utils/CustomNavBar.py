from PyQt5.QtCore import qVersion
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
                     (None, None, None, None),
                     ('Plot format', 'Edit rcParams, curve and image parameters', 'options', 'rcParams'),
                     ('Linear', 'Save the figure', 'lin', 'linear_scale'),
                     ('Log', 'Save the figure', 'log', 'log_scale'),
                     ('Decibel', 'Save the figure', 'dB', 'decibel_scale'),
                     ('Grid', 'Toggle Grid', 'dB', 'toggle_grid'),
                     (None, None, None, None),
                     ("Size", None, None, None),
                     ('Save', 'Save the figure', 'save', 'save_figure'),
                     ]

        self.canvas.toolbar.setContentsMargins(0, 0, 0, 0)
        # print(type(self.canvas.toolbar))

        actions = self.findChildren(QAction)
        for a in actions:
            self.removeAction(a)

        self._actions = {} # mapping of toolitem method names to QActions.

        for text, tooltip_text, image_file, callback in toolitems:
            if text is None:
                self.addSeparator()
            elif text == 'Size':
                # size definitions
                plot_settings_options = ['Presentation', 'Square', 'Portrait', 'Landscape', 'Custom']
                self.cb_Save_Plot_Settings = QComboBox()
                for a in plot_settings_options:
                    self.cb_Save_Plot_Settings.addItem(a)
                self.addWidget(self.cb_Save_Plot_Settings)
            elif text == 'Grid':
                self.cb_Grid = QCheckBox()
                self.cb_Grid.clicked.connect(lambda: self.toggle_grid())
                self.addWidget(self.cb_Grid)

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
                    self.canvas.fig.set_size_inches(8, 6) # actual size = 8*dpi, 6*dpi; dpi = 100
                elif self.cb_Save_Plot_Settings.currentText() == 'Square':
                    self.cb_Save_Plot_Settings.set_size_inches(8, 8) #, forward=False
                elif self.cb_Save_Plot_Settings.currentText() == 'Landscape':
                    self.cb_Save_Plot_Settings.set_size_inches(12, 8) #, forward=False
                self.canvas.fig.savefig(fileName, transparent=True, bbox_inches='tight', pad_inches=0.2)

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
        if self.cb_Grid.checkState() == 2:
            self.canvas.fig.get_axes()[0].grid(True, 'both')
        else:
            self.canvas.fig.get_axes()[0].grid(False)
        self.canvas.draw()

    def rcParams(self):
        # pop up plot properties widget
        # create pop up widget instance
        self.pprops.ppropsUI.cb_Plots.addItems([x.get_label() for x in self.canvas.fig.get_axes()])
        self.pprops.w_PProps.show()

