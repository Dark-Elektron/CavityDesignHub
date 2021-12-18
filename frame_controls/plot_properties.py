from PyQt5.QtWidgets import QWidget
from ui_files.pprops import Ui_PProps


class PPropsControl:
    def __init__(self, canvas):
        self.w_PProps = QWidget()

        self.canvas = canvas

        self.ppropsUI = Ui_PProps()
        self.ppropsUI.setupUi(self.w_PProps)

        self.initUI()
        self.signals()

    def initUI(self):
        pass

    def signals(self):
        self.ppropsUI.pb_Apply.clicked.connect(lambda: self.apply())
        self.ppropsUI.pb_Cancel.clicked.connect(lambda: self.close())

    def close(self):
        self.w_PProps.close()

    def apply(self):
        self.canvas.ax.tick_params(axis='x', labelsize=8)
        self.canvas.draw()
