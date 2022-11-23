#  Name   : exampleWidget
#
#          A PyQt5 widget to serve as an example of how to create custom widgets
#          that can be used within Qt designer.
#          This widget embeds Matplotlib canvas (plot space). It contains also
#          defined PyQt5 slots for setting the magnetics IDS parameters and
#          a slot (function) for executing the plot procedure,
#          populating/filling the Matplotlib canvas.
#
#  Author :
#         Dejan Penko
#  E-mail :
#         dejan.penko@lecad.fs.uni-lj.si
#
#****************************************************
#     Copyright(c) 2019- D. Penko

# import module providing system-specific parameters and functions
import sys
# import module providing miscellaneous operating system interfaces
import os
# import module providing log handling
import logging
# import modules providing PyQt5 parameters, functions etc.
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QVBoxLayout, QSizePolicy
from PyQt5.QtCore import pyqtSlot, pyqtSignal
# import module providing matplotlib parameters, functions etc.
import matplotlib
matplotlib.use('Qt5Agg') # Use Qt rendering
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
# import module providing IMAS and IDS-related parameters, functions etc.
# import imas

class exampleWidget(QWidget):
    """A widget for opening magnetics IDS, extracting the flux loop or
    poloidal probe quantities and plotting them to matplotlib figure canvas.
    """

    # Create a custom signal
    idsSet = pyqtSignal(bool)

    def __init__(self, parent=None, ids=None, *args, **kwargs):
        """
        Arguments:
            parent (PyQt5 object) : Qt widget parent (e.g. QMainWindow)
            ids    (IDS object)   : IDS object - optional parameter.
                                    This widget does not require IDS object in
                                    order to work (the default IDS parameters
                                    will be used to open the IDS).
                                    However, if an IDS object is already
                                    available it can be passed to the widget
                                    to be used instead (for example, passing
                                    the IDS object from the IMASViz to this
                                    widget).
        """
        # Run QWidget constructor
        super(QWidget, self).__init__(parent)

        # Check if display is available (display is mandatory, as this is
        # PyQt5 widget)
        self.checkDisplay()

        # Set IDS object
        # Note: if not provided as an argument it will be set to None
        self.ids = ids

        # Set IDS case parameters
        # - Empty dictionary
        self.idsParameters = {}
        # - shot
        self.idsParameters['shot'] = '52344'
        # - run r
        self.idsParameters['run'] = '0'
        # - user
        self.idsParameters['user'] = os.getenv('USER')
        # - device / machine / database name
        self.idsParameters['device'] = 'viztest'
        # - IMAS major version (3.x.y)
        self.idsParameters['IMAS major version'] = '3'
        # - label of the IDS to be used
        self.idsParameters['idsName'] = 'magnetics'

        # Set widget layout
        self.setLayout(QVBoxLayout())
        # Set empty matplotlib canvas
        self.canvas = PlotCanvas(self)
        # Set matplotlib toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        # Add canvas and toolbar to widget layout
        self.layout().addWidget(self.canvas)
        self.layout().addWidget(self.toolbar)

    @pyqtSlot()
    def plotFluxAoS(self):
        """Plot Flux Loop arrays to canvas.
        """

        # IDS check
        if self.ids == None:
            logging.error(' IDS was not set/opened!')
            return
        # Canvas figure check
        if self.canvas.figure != None:
            self.canvas.figure.clear()
        # Plot flux loop Aos
        self.canvas.plotFluxAoS(self.ids)

    @pyqtSlot()
    def plotBPolAoS(self):
        """Plot poloidal field probe arrays to canvas.
        """

        # IDS check
        if self.ids == None:
            logging.error(' IDS was not set/opened!')
            return
        # Canvas figure check
        if self.canvas.figure != None:
            self.canvas.figure.clear()
        # Plot Poloidal field AoS
        self.canvas.plotBPolAoS(self.ids)

    @pyqtSlot()
    def openIDS(self):
        """Open magnetics IDS.
        """
        # Open IDS
        self.ids = imas.ids(int(self.idsParameters['shot']),
                            int(self.idsParameters['run']))
        self.ids.open_env(self.idsParameters['user'],
                          self.idsParameters['device'],
                          self.idsParameters['IMAS major version'])
        # Get magnetics IDS
        self.ids.magnetics.get()

    def setIDS(self, ids):
        self.ids = ids
        # Emit signal indicating that the IDS object is set
        self.idsSet.emit(False)

    def getIDS(self):
        return self.ids

    @pyqtSlot(str)
    def setShot(self, shot):
        self.idsParameters['shot'] = shot

    def getShot(self):
        return self.idsParameters['shot']

    @pyqtSlot(str)
    def setRun(self, run):
        self.idsParameters['run'] = run

    def getRun(self):
        return self.idsParameters['run']

    @pyqtSlot(str)
    def setUser(self, user):
        self.idsParameters['user'] = user

    def getUser(self):
        return self.idsParameters['user']

    @pyqtSlot(str)
    def setDevice(self, device):
        self.idsParameters['device'] = device

    def getDevice(self):
        return self.idsParameters['device']

    @pyqtSlot(str)
    def setIMASmVer(self, ver):
        self.idsParameters['IMAS major version'] = ver

    def getIMASmVer(self):
        return self.idsParameters['IMAS major version']

    @pyqtSlot(str)
    def setIDSname(self, idsName):
        self.idsParameters['idsName'] = idsName

    def getIDSname(self):
        return self.idsParameters['idsName']

    @pyqtSlot()
    def checkDisplay(self):
        try:
            os.environ['DISPLAY']
        except:
            logging.error('No display available!')

class PlotCanvas(FigureCanvas):
    """Matplotlib figure canvas that is to be embedded within the widget.
    FigureCanvas is the area onto which the figure is drawn
    """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        """
        Arguments:
            parent (PyQt5 object) : PyQt5 parent (e.g. QWidget).
            width  (int)          : Canvas width.
            height (int)          : Canvas height.
            dpi    (int)          : Dots per inch.
        """

        # Set figure
        fig = Figure(figsize=(width, height), dpi=dpi)
        # Set canvas (pass figure)
        FigureCanvas.__init__(self, fig)
        # Set canvas parent
        self.setParent(parent)
        # Set canvas size policy
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)

    def plotFluxAoS(self, ids):
        """Plot values found in flux loops AoS.

        Arguments:
            ids (IDS object) : IDS object referring to the IDS from which the
                               data is to be extracted.
        """

        # Add/Update IDS reference to the object (figure canvas)
        self.ids = ids
        # IDS check
        if self.ids == None:
            logging.error('IDS was not set/opened!')
            return
        # Set subplot
        ax = self.figure.add_subplot(111)
        # Extract X-axis values (time)
        time_values = self.ids.magnetics.time
        x = time_values
        # Get the size of AoS (number of arrays)
        num_flux_loop_AoS = len(self.ids.magnetics.flux_loop)
        # For each array extract array values and create a plot
        for i in range(num_flux_loop_AoS):
            # Extract array values
            y = self.ids.magnetics.flux_loop[i].flux.data
            # Set plot (line) defined by X and Y values +
            # set line as full line (-) and add legend label.
            ax.plot(x, y, '-', label='Flux_loop[' + str(i) + ']')
        # Enable grid
        ax.grid()
        # Set axis labels and plot title
        ax.set(xlabel='time [s]', ylabel='Flux Loop values',
               title='Flux loop')
        # Enable legend
        ax.legend()
        # Draw/Show plots
        self.draw()

    def plotBPolAoS(self, ids):
        """Plot poloidal field probe values.

        Arguments:
            ids (IDS object) : IDS object referring to the IDS from which the
                               data is to be extracted.
        """
        # Add/Update IDS reference to the object (figure canvas)
        self.ids = ids
        # IDS check
        if self.ids == None:
            logging.error('IDS was not set/opened!')
            return
        # Set subplot
        ax = self.figure.add_subplot(111)
        # Extract X-axis values (time)
        time_values = self.ids.magnetics.time
        x = time_values
        # Get the size of AoS (number of arrays)
        num_bpol_probe_AoS = len(self.ids.magnetics.bpol_probe)
        # For each array extract array values and create a plot
        for i in range(num_bpol_probe_AoS):
            # Extract array values
            y = self.ids.magnetics.bpol_probe[i].field.data
            # Set plot (line) defined by X and Y values +
            # set line as full line (-) and add legend label.
            ax.plot(x, y, '-', label='bpol_probe[' + str(i) + ']')
        # Enable grid
        ax.grid()
        # Set axis labels and plot title
        ax.set(xlabel='time [s]', ylabel='Poloidal field probe values',
               title='Poloidal field probe')
        # Enable legend
        ax.legend()
        # Draw/Show plots
        self.draw()

if __name__ == '__main__':

    # Set application object
    app = QApplication(sys.argv)
    # Set main PyQt5 window
    mainWindow = QMainWindow()
    # Set window title
    mainWindow.setWindowTitle('Example Widget')
    # Set example widget object
    ew = exampleWidget()
    # Open IDS (magnetics IDS)
    ew.openIDS()
    # Plot Flux Loop arrays
    ew.plotFluxAoS()
    # Plot poloidal field probe arrays (an option other than plotFluxAoS)
    # ew.plotBPolAoS()
    # Set example widget as a central widget of the main window
    mainWindow.setCentralWidget(ew)
    # Show the main window
    mainWindow.show()
    # Keep the application running (until the 'exit application' command is
    # executed
    sys.exit(app.exec_())
