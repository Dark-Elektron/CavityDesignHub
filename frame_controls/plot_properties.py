from PyQt5.QtWidgets import QWidget
from ui_files.pprops import Ui_PProps
import matplotlib as mpl


class PPropsControl:
    def __init__(self, canvas):
        self.w_PProps = QWidget()

        self.canvas = canvas

        self.ppropsUI = Ui_PProps()
        self.ppropsUI.setupUi(self.w_PProps)

        self.initUI()
        self.signals()

        self.axins = None
        self.indicate_inset = None

    def initUI(self):
        self.plot_default_decorations()

    def plot_default_decorations(self):
        # plot settings
        mpl.rcParams['xtick.labelsize'] = 18
        mpl.rcParams['ytick.labelsize'] = 18

        mpl.rcParams['axes.labelsize'] = 18
        mpl.rcParams['axes.titlesize'] = 18
        mpl.rcParams['legend.fontsize'] = 'large'

        mpl.rcParams['figure.figsize'] = [9.8, 6]
        mpl.rcParams['figure.dpi'] = 100

    def signals(self):
        self.ppropsUI.pb_Apply.clicked.connect(lambda: self.apply())
        self.ppropsUI.pb_Cancel.clicked.connect(lambda: self.close())
        self.ppropsUI.cb_Show_Hide_Inset.clicked.connect(lambda: self.plot_inset())

    def plot_inset(self):
        if self.ppropsUI.cb_Show_Hide_Inset.checkState() == 2:
            # get lines from axis
            lines = self.canvas.ax.get_lines()

            # 2HC1FPC (Wakefield)%%3HC1FPC (Wakefield)%%4HC1FPC (Lossy Eigenmode)%%3HC1FPC (Lossy Eigenmode)%%2HC1FPC (Lossy Eigenmode)
            # sub region of the original image
            # get values from line edit
            x0 = self.ppropsUI.dsb_X0.value()
            x1 = self.ppropsUI.dsb_X1.value()
            y0 = self.ppropsUI.dsb_Y0.value()
            y1 = self.ppropsUI.dsb_Y1.value()

            # check if inset bounds is within current plot bound
            if x0 >= self.canvas.ax.get_xlim()[0] and x1 <= self.canvas.ax.get_xlim()[1] and y0 >= self.canvas.ax.get_ylim()[0] and y1 <= self.canvas.ax.get_ylim()[1]:

                inset_pos = [self.ppropsUI.dsb_XPos.value(), self.ppropsUI.dsb_YPos.value(), self.ppropsUI.dsb_XWidth.value(), self.ppropsUI.dsb_YWidth.value()]
                if len(inset_pos) == 4:
                    self.axins = self.canvas.ax.inset_axes(inset_pos, facecolor='#fafafa')
                    # self.axins.axes.xaxis.set_visible(False)
                    # self.axins.axes.yaxis.set_visible(False)

                for line in lines:
                    if self.axins:
                        if line.get_linestyle() == 'None':
                            self.axins.plot(line.get_xdata(), line.get_ydata(), linestyle='None',
                                            marker=line.get_marker(),
                                            markersize=line.get_ms(), c=line.get_color(), mec=line.get_mec(), alpha=line.get_alpha())
                        else:
                            self.axins.plot(line.get_xdata(), line.get_ydata(), ls=line.get_linestyle(), marker=line.get_marker(),
                                            mec=line.get_mec(), linewidth=line.get_lw(), c=line.get_color(), alpha=line.get_alpha())

                self.axins.set_xlim(x0, x1)
                self.axins.set_ylim(y0, y1)

                # self.axins.set_xticklabels('')
                # self.axins.set_yticklabels('')

                self.axins.set_yscale(self.ppropsUI.cb_Y_Scale.currentText())
                self.indicate_inset = self.canvas.ax.indicate_inset_zoom(self.axins, edgecolor="black", label='_nolegend_')
            else:
                print("Axins bounds not in plot bounds")
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

        self.canvas.draw_idle()

    def draw_legend(self):
        legend = self.canvas.ax.get_legend()
        if legend:
            legend.remove()

    def close(self):
        self.w_PProps.close()

    def apply(self):
        self.canvas.ax.tick_params(axis='x', labelsize=self.ppropsUI.sb_XTick.value())
        self.canvas.ax.tick_params(axis='y', labelsize=self.ppropsUI.sb_YTick.value())

        self.canvas.ax.set_xlabel(self.ppropsUI.le_XLabel.text(), fontsize=self.ppropsUI.sb_XLabel.value())
        self.canvas.ax.set_ylabel(self.ppropsUI.le_YLabel.text(), fontsize=self.ppropsUI.sb_YLabel.value())
        self.canvas.ax.set_title(self.ppropsUI.le_Title.text(), fontsize=self.ppropsUI.sb_Title.value())

        # legend

        legend_labels = self.ppropsUI.le_Legend.text().split("%%")
        # update legend
        handles, labels = self.canvas.ax.get_legend_handles_labels()

        try:
            for i in range(len(legend_labels)):
                labels[i] = legend_labels[i] if legend_labels[i] != "" else labels[i]
        except IndexError:
            pass

        self.leg = self.canvas.ax_right.legend(handles, labels, fontsize=self.ppropsUI.sb_Legend.text())
        self.leg.set_zorder(10)
        self.leg.set_draggable(state=True, use_blit=True)
        self.draw_legend()

        self.canvas.fig.tight_layout()
        self.canvas.draw_idle()
