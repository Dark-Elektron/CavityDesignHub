import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from PyQt5.QtWidgets import *
from termcolor import colored
from frame_controls.edit_line2d import EditALine2DDialog
from utils.file_reader import FileReader
from PyQt5 import QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # NavigationToolbar2QT as NavigationToolbar
from utils.CustomNavBar import MyCustomToolbar as NavigationToolbar
from matplotlib.figure import Figure
import warnings
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from analysis_modules.plot_module.matplotlib_annotation_objects import DraggableText, DraggableAxvline, DraggableRectangle
from frame_controls.edit_matplotlibobjects_dialog import EditATextDialog
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.path as mpath

matplotlib.use('Qt5Agg')

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.rcParams['toolbar'] = 'toolmanager'

fr = FileReader()


file_color = 'yellow'
DEBUG = True


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class Plot(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=12, height=9, dpi=100):

        self.event = None
        plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Set2.colors)
        self.event_y = None
        self.event_x = None
        self.eal = None
        self.ax_right = None
        self.eat = None
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.set_tight_layout("True")
        self.fig_size = self.fig.get_size_inches()*self.fig.dpi

        self.parent_ = parent
        self.parentUI = parent.ui

        # new feature
        # add gridspec
        self.plot_list = []
        self.n = 7
        n = self.n
        self.gs = gridspec.GridSpec(1, 1, figure=self.fig)
        self.gs2 = gridspec.GridSpec(n, n, height_ratios=np.ones(n), width_ratios=np.ones(n), figure=self.fig)

        for i in range(n):
            if i < n-1:
                ax = self.fig.add_subplot(self.gs2[i, 0])
                ax2 = self.fig.add_subplot(self.gs2[n-1, n-1-i])
                # turn off axis
                self.plot_list.append(ax)
                self.plot_list.append(ax2)
            else:
                ax = self.fig.add_subplot(self.gs2[n-1, 0])
                self.plot_list.append(ax)

        self.ax = self.fig.add_subplot(self.gs[0])
        self.toggle_ax(False)
        # self.ax = self.fig.add_subplot(111)
        self.ax.set_label("Left axis Default")
        self.ax.set_zorder(1)

        mpl.rcParams['savefig.facecolor'] = 'ffffff'

        super(Plot, self).__init__(self.fig)

        self.widgetUI = parent.ui
        # self.plotUI = self.plot_widget.plotUI
        # self.MROT = self.ui.cb_Polarization_ABCI.currentIndex()

        self.GLOBAL_STATE = 0

        # mouse constants
        self.PRESS = False
        self.RELEASE = False
        self.MOTION = False
        self.PICK = False

        # plot appearance
        # self.aesthetics()

        # handle events
        # self.event_handler()

        # add annotation
        annot = self.ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                 bbox=dict(boxstyle="round", fc="w"),
                                 arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)

        # zoom and pan
        scale = 1.1
        zp = ZoomPan()
        figZoom = zp.zoom_factory(self.ax, base_scale=scale)
        # figPan = zp.pan_factory(self.axes)

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        self.toolbar = NavigationToolbar(self, None)

        # create new axis
        self.create_secondary_axes()

        # # create new widgets
        # self.create_new_widgets()
        self.w_Color = QColorDialog()
        self.w_Color.currentColorChanged.connect(lambda: self.update_object_properties(
            {'color': self.w_Color.currentColor().name()}))

        # redefine save_image function
        self.toolbar.sizeHint()
        self.widgetUI.gl_Plot_Area.addWidget(self.toolbar)

        # handle events
        self.event_handler()

        self.signals()
        self.text_dict = {}
        self.axvline_dict = {}
        self.patch_dict = {}

        self.selected_object = None

        self.picked_points_annotext_id_dict = {}
        self.selected_points_list = []

    def toggle_ax(self, visible):
        """

        :param visible:
        :return:
        """
        n = int((len(self.plot_list)+1)/2)
        if visible:
            self.ax.set_position(self.gs2[0:-2, 2:].get_position(self.fig))
            for i, ax in enumerate(self.plot_list):
                ax.set_visible(visible)

                if i < n:
                    ax.set_position(self.gs2[i, 0].get_position(self.fig))
                else:
                    ax.set_position(self.gs2[n-1, i-(n-1)].get_position(self.fig))
        else:
            for ax in self.plot_list:
                ax.set_visible(False)
            self.ax.set_position(self.gs[0].get_position(self.fig))

    def signals(self):
        pass

    def change_scale(self, arg):
        self.ax.set_yscale(arg)
        self.ax_right.set_yscale(arg)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    @staticmethod
    def plot_settings(mode):
        # plt.style.use(['science', 'no-latex'])
        if mode == 'presentation':
            # plot settings
            mpl.rcParams['xtick.labelsize'] = 20
            mpl.rcParams['ytick.labelsize'] = 20

            mpl.rcParams['axes.labelsize'] = 24
            mpl.rcParams['axes.titlesize'] = 24
            mpl.rcParams['legend.fontsize'] = 'large'

            mpl.rcParams['figure.figsize'] = [11, 10]
            mpl.rcParams['figure.dpi'] = 100

        if mode.lower() == 'square':
            # plot settings
            mpl.rcParams['xtick.labelsize'] = 20
            mpl.rcParams['ytick.labelsize'] = 20

            mpl.rcParams['axes.labelsize'] = 24
            mpl.rcParams['axes.titlesize'] = 24
            mpl.rcParams['legend.fontsize'] = 'large'

            mpl.rcParams['figure.figsize'] = [10, 10]
            mpl.rcParams['figure.dpi'] = 100

    def create_secondary_axes(self):
        # Create multiple axes
        # twin object for two different y-axis on the sample plot
        self.ax_right = self.ax.twinx()
        # turn off axis ticks and labels
        self.ax_right.axes.xaxis.set_visible(False)
        self.ax_right.axes.yaxis.set_visible(False)
        self.ax_right.set_label("Right axis Default")

    def event_handler(self):
        self.cid_button_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_button_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.cid_motion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

        # self.cid_figure_enter = self.fig.canvas.mpl_connect('figure_enter_event', self.enter_figure)
        # self.cid_figure_leave = self.fig.canvas.mpl_connect('figure_leave_event', self.leave_figure)
        # self.cid_axes_enter = self.fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
        # self.cid_axes_leave = self.fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)

    def on_press(self, event):

        # check if zoom or pan mode is active to avoid selection of matplotlib objects
        if self.toolbar.mode.name == 'PAN' or self.toolbar.mode.name == 'ZOOM':
            return 0

        # get event position in fig size
        self.event_x, self.event_y = event.x, event.y

        # set global variable to True
        self.PRESS = True

        # handle identifying objects
        if event.button == 1:  # LEFT
            if event.dblclick:
                if isinstance(self.selected_object, matplotlib.legend.Legend):
                    pass

                if isinstance(self.selected_object, matplotlib.text.Annotation):
                    self.text_dict[f'{id(self.selected_object)}'].connect()
                    pass

                if isinstance(self.selected_object, matplotlib.collections.PathCollection):
                    obj_id = f'{id(self.selected_object)}'
                    if obj_id in self.picked_points_annotext_id_dict.keys():
                        annotext_id = self.picked_points_annotext_id_dict[obj_id]

                        # remove annotated text object from axes
                        self.text_dict[str(annotext_id)].remove()
                        # remove annotated text object from dictionary
                        del self.text_dict[str(annotext_id)]

                        # remove scatter object from axes
                        self.selected_object.remove()
                        # remove scatter object form dictionary
                        del self.picked_points_annotext_id_dict[obj_id]

                        self.draw_idle()
                        self.flush_events()
                        # return
                elif isinstance(self.selected_object, matplotlib.lines.Line2D):
                    ind = self.event.ind
                    x_picked = self.selected_object.get_xdata()[min(ind)]
                    y_picked = self.selected_object.get_ydata()[min(ind)]

                    # add annotated text at picked point
                    # calculate location of annotation from axes extents and event position
                    if "left" in self.parentUI.cb_Active_Axis.currentText().lower():
                        ax = self.ax
                    else:
                        ax = self.ax_right

                    xlim, ylim = ax.get_xlim(), ax.get_ylim()
                    # pos_x, pos_y = (x_picked - xlim[0])/(xlim[1] - xlim[0]), (y_picked - ylim[0])/(ylim[1] - ylim[0])
                    pos_x, pos_y = x_picked, y_picked

                    # check if newly selected points is close to any of the old selected points
                    for i in range(len(self.selected_points_list)):
                        pl = self.selected_points_list[i]
                        if (pl[0] - x_picked) / pl[0] < 0.2 and (pl[1] - y_picked) / pl[1] < 0.2:
                            # delete point from list
                            del self.selected_points_list[i]
                            return

                    # plot circle around the point
                    scatter_obj = ax.scatter(pos_x, pos_y, s=100, marker='o', facecolors="None",
                                             edgecolor='red', picker=True, zorder=5)

                    x_scale, y_scale = ax.get_xscale(), ax.get_yscale()

                    # get x_axis title
                    xlbl = self.ax.xaxis.get_label().get_text()
                    ylbl = self.ax.yaxis.get_label().get_text()

                    if xlim[0] <= pos_x <= xlim[1] / 2:
                        if x_scale == 'linear':
                            text_x = pos_x + (xlim[1] - xlim[0]) / 20
                        else:
                            text_x = pos_x + np.log((xlim[1] - xlim[0]) / 20)
                    else:
                        if x_scale == 'linear':
                            text_size = 30 * len(f"{xlbl}: {round(x_picked, 4)}") / self.fig_size[1]
                            text_x = pos_x - (xlim[1] - xlim[0]) / 20 - (text_size - xlim[0]) / (xlim[1] - xlim[0])
                        else:
                            text_size = 30 * len(f"{xlbl}: {round(x_picked, 4)}") / self.fig_size[1]
                            text_x = pos_x - np.log((xlim[1] - xlim[0]) / 20) - (text_size - xlim[0]) / (
                                        xlim[1] - xlim[0])

                    if ylim[0] <= pos_y <= ylim[1] / 2:
                        if y_scale == 'linear':
                            text_y = pos_y + (ylim[1] - ylim[0]) / 5
                        else:
                            text_y = pos_y
                    else:
                        if y_scale == 'linear':
                            text_y = pos_y - (ylim[1] - ylim[0]) / 5
                        else:
                            text_y = pos_y

                    text = f"{xlbl}: {round(x_picked, 4)}\n {ylbl}: {round(y_picked, 4)}"
                    annotext = self.add_text(text, "Round4", xy=(x_picked, y_picked), xycoords='data',
                                             xytext=(text_x, text_y), textcoords='data',
                                             size=14, rotation=0,
                                             arrowprops=dict(arrowstyle='-', facecolor='black'))

                    self.picked_points_annotext_id_dict[str(id(scatter_obj))] = str(id(annotext))
                    self.selected_points_list.append([x_picked, y_picked])

                    # self.axvline_dict[f"{id(self.selected_object)}"].connect()

                    # xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
                    # x, y = artist.get_xdata(), artist.get_ydata()
                    # ind = event.ind

                self.draw()
            else:
                if event.inaxes is not None:
                    for patch in self.ax.patches:
                        if patch.contains(event)[0]:
                            self.selected_object = patch
                            # check format in which color was given
                            if isinstance(patch.get_facecolor(), list) or isinstance(patch.get_facecolor(), tuple):
                                color = mpl.colors.to_hex(patch.get_facecolor())
                            else:
                                color = patch.get_facecolor()

                            self.color_widgets(color)

                            # patch.set_edgecolor('red')  # Highlight the selected patch
                            # self.fig.canvas.draw_idle()
                # update color and properties
                if isinstance(self.selected_object, matplotlib.lines.Line2D):
                    # set color picker color

                    # check format in which color was given
                    if isinstance(self.selected_object.get_color(), list) or isinstance(self.selected_object.get_color(), tuple):
                        color = mpl.colors.to_hex(self.selected_object.get_color())
                    else:
                        color = self.selected_object.get_color()

                    self.color_widgets(color, self.selected_object)

        if event.button == 2:  # MMB
            pass

        if event.button == 3:  # RIGHT
            # selecting objects right click
            if isinstance(self.selected_object, matplotlib.text.Annotation):
                self.draggableTextContextMenuEvent(event)
            elif isinstance(self.selected_object, matplotlib.lines.Line2D):
                self.draggableAxvlineContextMenuEvent(event)
            elif event.inaxes is not None:
                for patch in self.ax.patches:
                    if patch.contains(event)[0]:
                        self.patchContextMenuEvent(event, patch)
            else:
                self.contextMenuEvent(event)

    def on_release(self, event):
        # check if zoom or pan mode is active to avoid selection of matplotlib objects
        if self.toolbar.mode.name == 'PAN' or self.toolbar.mode.name == 'ZOOM':
            return 0

        # reset selected object to none
        # self.selected_object = None
        self.PRESS = False
        # set press to false

    def on_pick(self, event):

        ########################################################
        # # legend toggle
        # # on the pick event, find the orig line corresponding to the
        # # legend proxy line, and toggle the visibility
        # legline = event.artist

        # origline = self.lined[legline]

        # vis = not origline.get_visible()
        # origline.set_visible(vis)
        # # Change the alpha on the line in the legend so we can see what lines
        # # have been toggled
        # if vis:
        #     legline.set_alpha(1.0)
        # else:
        #     legline.set_alpha(0.2)
        #############################################################

        # check if zoom or pan mode is active to avoid selection of matplotlib objects
        if self.toolbar.mode.name == 'PAN' or self.toolbar.mode.name == 'ZOOM':
            return 0

        if event.artist is not None:
            self.selected_object = event.artist
            self.event = event

    def on_motion(self, event):

        # check if zoom or pan mode is active to avoid selection of matplotlib objects
        if self.toolbar.mode.name == 'PAN' or self.toolbar.mode.name == 'ZOOM':

            return 0

        # vis = self.annot.get_visible()
        # if isinstance(self.ax, event.inaxes):  # not so good fix
        #     cont, ind = self.plot_object.contains(event)
        #     if cont:
        #         self.update_annot(
        #             ind["ind"][0])  # ind returns an array of close points. ind['ind'][0] returns just the first point
        #         self.annot.set_visible(True)
        #         self.fig.canvas.draw_idle()
        #     else:
        #         if vis:
        #             self.annot.set_visible(False)
        #             self.fig.canvas.draw_idle()
        self.MOTION = True

        # if self.PRESS and self.MOTION:

    def enter_axes(self, event):
        event.inaxes.patch.set_facecolor('white')
        event.canvas.draw()

    def leave_axes(self, event):
        event.inaxes.patch.set_facecolor('white')
        event.canvas.draw()

    def enter_figure(self, event):
        fig = event.canvas.figure
        fig.patch.set_facecolor('white')

        # get canvas size
        self.fig_size = fig.get_size_inches()*fig.dpi

        event.canvas.draw()

    def leave_figure(self, event):
        event.canvas.figure.patch.set_facecolor('white')
        event.canvas.draw()

    def update_annot(self, ind):
        x, y = self.line.get_data()
        self.annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
        text = "{}".format(" ".join(list(map(str, ind["ind"]))))
        self.annot.set_text(text)
        self.annot.get_bbox_patch().set_alpha(0.4)

        self.annot.set_text(f'{text}')
        self.annot.get_bbox_patch().set_facecolor((167 / 255, 222 / 255, 255 / 255))
        # self.annot.get_bbox_patch().set_alpha(1)

    def update_object_properties(self, props):
        if self.selected_object is not None:
            if isinstance(self.selected_object, matplotlib.legend.Legend):
                pass
            elif isinstance(self.selected_object, matplotlib.lines.Line2D) or isinstance(self.selected_object, matplotlib.collections.PathCollection):

                if 'color' in props.keys():
                    color = props['color']

                    self.selected_object.set_color(color)
                    self.selected_object.set_markerfacecolor(color)
                    self.parent_.le_Color.setStyleSheet(f'background-color: {color};')
                    self.parent_.pb_Color.setStyleSheet(f"background-color: {color};")
                    self.parent_.le_Color.setText(color)
                    self.parent_.pb_Color.setText(color)
                    # copy the image to the GUI state, but screen might not be changed yet
                if 'alpha' in props.keys():
                    alpha = props['alpha']
                    self.selected_object.set_alpha(alpha)

            elif isinstance(self.selected_object, matplotlib.patches.Rectangle) or isinstance(self.selected_object, matplotlib.patches.Circle):
                if 'color' in props.keys():
                    color = props['color']
                    self.selected_object.set(facecolor=color)
                if 'alpha' in props.keys():
                    alpha = props['alpha']
                    self.selected_object.set_alpha(alpha)

            self.draw_idle()
            self.flush_events()
            self.parent_.update_labels()

    def contextMenuEvent(self, event):
        contextMenu = QMenu(self)
        textAct = contextMenu.addAction("Add Text")
        arrowAct = contextMenu.addAction("Add Arrow")
        arrowAct = contextMenu.addAction("Move axis marker to minimum")
        arrowAct = contextMenu.addAction("Move axis marker to maximum")

        # get figure size and subtract from event.y because event.y counts bottom up.
        action = contextMenu.exec_(
            self.mapToGlobal(QtCore.QPoint(event.x, int(self.fig.get_size_inches()[1] * self.fig.dpi) - event.y)))

        if action == textAct:
            # create pop up widget
            self.eat = EditATextDialog(self, self.selected_object)
            self.eat.show()

    def draggableTextContextMenuEvent(self, event):
        contextMenu = QMenu(self)
        editAct = contextMenu.addAction("Edit")
        deleteAct = contextMenu.addAction("Delete")

        # get figure size and subtract from event.y because event.y counts bottom up
        action = contextMenu.exec_(
            self.mapToGlobal(QtCore.QPoint(event.x, int(self.fig.get_size_inches()[1] * self.fig.dpi) - event.y)))

        if action == deleteAct:
            self.remove_text()

        elif action == editAct:
            # try:
            #     self.text_dict[f"{id(self.selected_object)}"].disconnect()
            # except KeyError:
            #     pass

            # create pop up widget
            self.eat = EditATextDialog(self, self.text_dict[f"{id(self.selected_object)}"])
            self.eat.show()

    def draggableAxvlineContextMenuEvent(self, event):
        contextMenu = QMenu(self)
        editAct = contextMenu.addAction("Edit")
        deleteAct = contextMenu.addAction("Delete")

        # get figure size and subtract from event.y because event.y counts bottom up.
        action = contextMenu.exec_(
            self.mapToGlobal(QtCore.QPoint(event.x, self.fig.get_size_inches()[1] * self.fig.dpi - event.y)))

        if action == deleteAct:
            self.remove_axvline()
            self.axvline_dict.pop(f"{id(self.selected_object)}")
            self.draw()

        elif action == editAct:
            self.axvline_dict[f"{id(self.selected_object)}"].on_release(event)
            # create pop up widget
            self.eal = EditALine2DDialog(self, self.selected_object)
            self.eal.show()

    def patchContextMenuEvent(self, event, patch):
        contextMenu = QMenu(self)
        # editAct = contextMenu.addAction("Edit")
        deleteAct = contextMenu.addAction("Delete")

        # get figure size and subtract from event.y because event.y counts bottom up
        action = contextMenu.exec_(
            self.mapToGlobal(QtCore.QPoint(event.x, int(self.fig.get_size_inches()[1] * self.fig.dpi) - event.y)))

        if action == deleteAct:
            patch.remove()
            self.patch_dict.pop(f"{id(patch)}")
            self.draw_idle()
            self.flush_events()
        #
        # elif action == editAct:
        #     try:
        #         self.patch_dict[f"{id(patch)}"].disconnect()
        #     except KeyError:
        #         pass
        #
        #     # create pop up widget
        #     self.epd = EditPatchDialog(self, self.selected_object)
        #     self.epd.show()

    def add_text(self, text, box, xy=(0.5, 0.5), xycoords='data', xytext=None, textcoords='data',
                 size=14, rotation=0, arrowprops=None):
        """

        Parameters
        ----------
        text
        box
        xy
        xycoords
        xytext
        textcoords
        size
        rotation
        arrowprops

        Returns
        -------

        """
        if text.strip("") == "":
            return

        if "left" in self.parentUI.cb_Active_Axis.currentText().lower():
            ax = self.ax
        else:
            ax = self.ax_right

        # add text
        if xytext:
            bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
            annotext = ax.annotate(text, xy=xy, xycoords=xycoords,
                                   xytext=xytext, textcoords=textcoords, bbox=bbox_props, fontsize=size,
                                   rotation=rotation, arrowprops=arrowprops)
        else:
            if box == "None":
                annotext = ax.annotate(text, xy=xy, xycoords=xycoords, fontsize=size,
                                       rotation=rotation, arrowprops=arrowprops)
            else:
                bbox_props = dict(boxstyle='{}'.format(box), fc='w', ec='k')
                annotext = ax.annotate(text, xy=xy, xycoords=xycoords, bbox=bbox_props, fontsize=size,
                                       rotation=rotation, arrowprops=arrowprops)

        self.draw_idle()
        self.flush_events()

        dt = DraggableText(annotext)
        dt.connect()
        self.text_dict[f'{id(annotext)}'] = dt

        return dt

    def remove_text(self):
        selected_obj_id = f'{id(self.selected_object)}'
        self.text_dict[selected_obj_id].remove()
        # remove key from dictionary
        del self.text_dict[selected_obj_id]

        self.draw_idle()
        self.flush_events()

    def add_axvline(self, x):

        # add vertical line
        axvline = self.ax_right.axvline(x, label=f"{x}", ls='--', c='gray', picker=True)
        self.draw_idle()
        self.flush_events()

        dt = DraggableAxvline(axvline)
        dt.connect()
        self.axvline_dict[f"{id(axvline)}"] = dt

        return axvline

    def remove_axvline(self):
        self.axvline_dict[f"{id(self.selected_object)}"].remove()
        self.draw()

    def add_patch(self, pos, a, b, kind='rectangle', ec='none', fc='red', alpha=0.1):
        # check inputs
        try:
            pos = tuple([eval(n) for n in pos])
            a = eval(a)
            b = eval(b)
            if kind == 'rectangle':
                self.add_rectangle(pos, a, b, ec=ec, fc=fc, alpha=alpha)
            elif kind == 'ellipse':
                self.add_ellipse(pos, a, b, ec=ec, fc=fc, alpha=alpha)
            self.draw_idle()
            self.flush_events()
        except (NameError, SyntaxError) as e:
            print(f"Exception occurred -> {e}")

    def remove_patch(self):
        pass

    def add_rectangle(self, pos, a, b, ec='none', fc='red', alpha=0.3):
        rect = mpatches.Rectangle(pos, a, b, ec=ec, fc=fc, alpha=alpha, lw=2, zorder=10)

        # add to axis
        if "left" in self.parentUI.cb_Active_Axis.currentText().lower():
            ax = self.ax
        else:
            ax = self.ax_right

        ax.add_artist(rect)

        dp = DraggableRectangle(rect)
        dp.connect()

        # add to path dict
        self.patch_dict[f'{id(rect)}'] = dp

    def add_ellipse(self, pos, a, b, ec='b', fc='b', alpha=0.3):
        # add rectangule object
        ellipse = mpatches.Ellipse(pos, a, b, ec=ec, fc=fc, alpha=alpha)

        # add to axis
        if "left" in self.parentUI.cb_Active_Axis.currentText().lower():
            ax = self.ax
        else:
            ax = self.ax_right
        ax.add_artist(ellipse)

        # add to path dict
        self.patch_dict[f'{id(ellipse)}'] = ellipse
        #
        pass

    def remove_rectangle(self):
        pass

    def remove_ellipse(self):
        pass

    def clear(self):
        # remove all text
        for key, val in self.text_dict.items():
            try:
                val.remove()
            except ValueError as e:
                print_("plotter.py: Exception:: ", e)

        # remove all lines
        for key, val in self.axvline_dict.items():
            try:
                val.remove()
            except ValueError as e:
                print_("plotter.py: Exception:: ", e)

        # reset dictionaries
        self.text_dict = {}
        self.axvline_dict = {}

    def change_background(self, color='white'):
        self.ax.set_facecolor(color)
        self.fig.patch.set_facecolor(color)
        self.draw_idle()
        self.flush_events()

    def color_widgets(self, color, obj=None):
        self.w_Color.setCurrentColor(QtGui.QColor(color))
        self.parent_.le_Color.setStyleSheet(f'background-color: {color};')
        self.parent_.pb_Color.setStyleSheet(f"background-color: {color};")
        self.parent_.le_Color.setText(color)
        self.parent_.pb_Color.setText(color)

        alpha = self.selected_object.get_alpha()
        if alpha:
            self.parentUI.dsb_Alpha.setValue(alpha)
        else:
            self.parentUI.dsb_Alpha.setValue(1)

        if obj is not None:
            self.parentUI.dsb_Line_Width.setValue(self.selected_object.get_linewidth())
            self.parentUI.cb_Line_Style.setCurrentText(self.selected_object.get_linestyle())
            self.parentUI.cb_Marker.setCurrentText(self.selected_object.get_marker())
            self.parentUI.db_Marker_Size.setValue(self.selected_object.get_markersize())


class ZoomPan:
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None

    @staticmethod
    def zoom_factory(ax, base_scale=2.):
        def zoom(event):
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()

            inv = ax.transData.inverted()
            pos = [event.x, event.y]

            xdata, ydata = inv.transform(pos)

            if event.button == 'down':
                # deal with zoom in
                scale_factor = 1 / base_scale
            elif event.button == 'up':
                # deal with zoom out
                scale_factor = base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1

            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

            relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])

            ax.set_xlim([xdata - new_width * (1 - relx), xdata + new_width * relx])
            ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * rely])
            ax.figure.canvas.draw()

        fig = ax.get_figure()  # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom

    def pan_factory(self, ax):
        def onPress(event):
            if event.inaxes != ax:
                return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()

            inv = ax.transData.inverted()
            pos = [event.x, event.y]

            xdata, ydata = inv.transform(pos)

            self.press = self.x0, self.y0, xdata, ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press

        def onRelease(event):
            self.press = None
            ax.figure.canvas.draw()

        def onMotion(event):
            if self.press is None:
                return
            if event.inaxes != ax:
                return

            inv = ax.transData.inverted()
            pos = [event.x, event.y]

            xdata, ydata = inv.transform(pos)

            dx = xdata - self.xpress
            dy = ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax.set_xlim(self.cur_xlim)
            ax.set_ylim(self.cur_ylim)

            ax.figure.canvas.draw()

        fig = ax.get_figure()  # get the figure of interest

        # attach the call back
        fig.canvas.mpl_connect('button_press_event', onPress)
        fig.canvas.mpl_connect('button_release_event', onRelease)
        fig.canvas.mpl_connect('motion_notify_event', onMotion)

        # return the function
        return onMotion


if __name__ == '__main__':
    p = Plot()
