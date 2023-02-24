# draggable rectangle with the animation blit techniques; see
# http://www.scipy.org/Cookbook/Matplotlib/Animations
# import numpy as np
# import matplotlib.pyplot as plt

class DraggableRectangle:
    lock = None  # only one can be animated at a time

    def __init__(self, rect):
        self.rect = rect
        self.press = None
        self.background = None

    def connect(self):
        """connect to all the events we need"""
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.rect.axes: return
        if DraggableRectangle.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        print('event contains', self.rect.xy)
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        DraggableRectangle.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        """on motion we will move the rect if the mouse is over us"""
        if DraggableRectangle.lock is not self:
            return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.rect.set_x(x0 + dx)
        self.rect.set_y(y0 + dy)

        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.rect)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """on release we reset the press data"""
        if DraggableRectangle.lock is not self:
            return

        self.press = None
        DraggableRectangle.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)


class DraggableText:
    lock = None  # only one can be animated at a time

    def __init__(self, text):
        self.text = text
        self.press = None
        self.background = None
        self.prev_text = None

        self.text.set_picker(True)

        self.x = 0.5
        self.y = 0.5

        self.cid_press = None
        self.cid_release = None
        self.cid_motion = None

    def connect(self):
        """connect to all the events we need"""
        self.cid_press = self.text.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.text.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.text.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.text.axes:
            return
        if DraggableText.lock is not None:
            return

        contains, attrd = self.text.contains(event)

        if not contains:
            return

        x0, y0 = self.text.xy

        # save location details at press
        self.press = x0, y0, event.xdata, event.ydata
        DraggableText.lock = self

        # draw everything but the selected text and store the pixel buffer
        canvas = self.text.figure.canvas
        axes = self.text.axes
        self.text.set_animated(True)
        canvas.draw_idle()
        canvas.flush_events()
        self.background = canvas.copy_from_bbox(self.text.axes.bbox)

        # now redraw just the Text
        axes.draw_artist(self.text)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        """on motion, we will move the text if the mouse is over us"""
        if DraggableText.lock is not self:
            return
        if event.inaxes != self.text.axes:
            return

        # get text canvas and axes
        fig = self.text.figure
        canvas = self.text.figure.canvas
        axes = self.text.axes

        if self.text.xycoords == 'data':
            # get size of figure
            bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            width, height = bbox.width * fig.dpi, bbox.height * fig.dpi

            # calculate relative position as fraction of figure
            self.x = event.x  # / width
            self.y = event.y  # / height

            # self.x, self.y = fig.transFigure.inversed().transform()
            self.x, self.y = axes.transData.inverted().transform((self.x, self.y))
            self.text.set_position((self.x, self.y))
        elif self.text.xycoords == 'axes fraction':

            # calculate relative position as fraction of figure
            self.x = event.x  # / width
            self.y = event.y  # / height

            self.x, self.y = axes.transAxes.inverted().transform((self.x, self.y))
            self.text.set_position((self.x, self.y))
        elif self.text.xycoords == 'figure':

            # calculate relative position as fraction of figure
            self.x = event.x  # / width
            self.y = event.y  # / height

            self.x, self.y = fig.transFigure.inverted().transform((self.x, self.y))
            self.text.set_position((self.x, self.y))

        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current text
        axes.draw_artist(self.text)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """on release, we reset the press data"""
        if DraggableText.lock is not self:
            return

        # reset press data
        self.press = None
        DraggableText.lock = None

        # turn off the text animation property and reset the background
        self.text.set_animated(False)
        self.background = None

        # redraw the full figure
        self.text.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.text.figure.canvas.mpl_disconnect(self.cid_press)
        self.text.figure.canvas.mpl_disconnect(self.cid_release)
        self.text.figure.canvas.mpl_disconnect(self.cid_motion)

    def remove(self):
        self.text.remove()


class DraggableSquare:
    lock = None  # only one can be animated at a time

    def __init__(self, text):
        self.text = text
        self.press = None
        self.background = None
        self.prev_text = None

        self.text.set_picker(True)

        self.x = 0.5
        self.y = 0.5

        self.cid_press = None
        self.cid_release = None
        self.cid_motion = None

    def connect(self):
        """connect to all the events we need"""
        self.cid_press = self.text.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.text.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.text.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.text.axes:
            return
        if DraggableText.lock is not None:
            return

        contains, attrd = self.text.contains(event)

        if not contains:
            return

        x0, y0 = self.text.xy

        # save location details at press
        self.press = x0, y0, event.xdata, event.ydata
        DraggableText.lock = self

        # draw everything but the selected text and store the pixel buffer
        canvas = self.text.figure.canvas
        axes = self.text.axes
        self.text.set_animated(True)
        canvas.draw_idle()
        canvas.flush_events()
        self.background = canvas.copy_from_bbox(self.text.axes.bbox)

        # now redraw just the Text
        axes.draw_artist(self.text)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        """on motion, we will move the text if the mouse is over us"""
        if DraggableText.lock is not self:
            return
        if event.inaxes != self.text.axes:
            return

        # get text canvas and axes
        fig = self.text.figure
        canvas = self.text.figure.canvas
        axes = self.text.axes

        if self.text.xycoords == 'data':
            # get size of figure
            bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            width, height = bbox.width * fig.dpi, bbox.height * fig.dpi

            # calculate relative position as fraction of figure
            self.x = event.x  # / width
            self.y = event.y  # / height

            # self.x, self.y = fig.transFigure.inversed().transform()
            self.x, self.y = axes.transData.inverted().transform((self.x, self.y))
            self.text.set_position((self.x, self.y))
        elif self.text.xycoords == 'axes fraction':

            # calculate relative position as fraction of figure
            self.x = event.x  # / width
            self.y = event.y  # / height

            self.x, self.y = axes.transAxes.inverted().transform((self.x, self.y))
            self.text.set_position((self.x, self.y))
        elif self.text.xycoords == 'figure':

            # calculate relative position as fraction of figure
            self.x = event.x  # / width
            self.y = event.y  # / height

            self.x, self.y = fig.transFigure.inverted().transform((self.x, self.y))
            self.text.set_position((self.x, self.y))

        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current text
        axes.draw_artist(self.text)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """on release, we reset the press data"""
        if DraggableText.lock is not self:
            return

        # reset press data
        self.press = None
        DraggableText.lock = None

        # turn off the text animation property and reset the background
        self.text.set_animated(False)
        self.background = None

        # redraw the full figure
        self.text.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.text.figure.canvas.mpl_disconnect(self.cid_press)
        self.text.figure.canvas.mpl_disconnect(self.cid_release)
        self.text.figure.canvas.mpl_disconnect(self.cid_motion)

    def remove(self):
        self.text.remove()


class DraggableArrow:
    lock = None  # only one can be animated at a time

    def __init__(self, rect):
        self.cidrelease = None
        self.cidmotion = None
        self.cidpress = None
        self.rect = rect
        self.press = None
        self.background = None

    def connect(self):
        """connect to all the events we need"""
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.rect.axes: return
        if DraggableRectangle.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains:
            return

        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        DraggableRectangle.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        """on motion we will move the rect if the mouse is over us"""
        if DraggableRectangle.lock is not self:
            return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.rect.set_x(x0 + dx)
        self.rect.set_y(y0 + dy)

        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.rect)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """on release we reset the press data"""
        if DraggableRectangle.lock is not self:
            return

        self.press = None
        DraggableRectangle.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)


class DraggableAxvline:
    lock = None  # only one can be animated at a time

    def __init__(self, axvline):
        self.cidmotion = None
        self.cidrelease = None
        self.cidpress = None
        self.axvline = axvline
        self.press = None
        self.background = None

    def connect(self):
        """connect to all the events we need"""
        self.cidpress = self.axvline.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.axvline.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.axvline.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.axvline.axes:
            return
        if DraggableAxvline.lock is not None:
            return
        contains, attrd = self.axvline.contains(event)
        if not contains:
            return

        x0 = self.axvline.get_xdata()
        self.press = x0, event.xdata
        DraggableAxvline.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.axvline.figure.canvas
        axes = self.axvline.axes
        self.axvline.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.axvline.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.axvline)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        """on motion, we will move the axvline if the mouse is over us"""
        if DraggableAxvline.lock is not self:
            return
        if event.inaxes != self.axvline.axes: return
        x0, xpress = self.press
        dx = event.xdata - xpress

        self.axvline.set_xdata([x+dx for x in x0])

        canvas = self.axvline.figure.canvas
        axes = self.axvline.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.axvline)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """on release, we reset the press data"""
        if DraggableAxvline.lock is not self:
            return

        self.press = None
        DraggableAxvline.lock = None

        # turn off the axvline animation property and reset the background
        self.axvline.set_animated(False)
        self.background = None
        # redraw the full figure
        self.axvline.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.axvline.figure.canvas.mpl_disconnect(self.cidpress)
        self.axvline.figure.canvas.mpl_disconnect(self.cidrelease)
        self.axvline.figure.canvas.mpl_disconnect(self.cidmotion)

    def remove(self):
        self.axvline.remove()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# rects = ax.bar(range(10), 20*np.random.rand(10))
# drs = []
#
# for rect in rects:
#     dr = DraggableRectangle(rect)
#     dr.connect()
#     drs.append(dr)
#
# # add text
# bbox_props = dict(boxstyle='round', fc='r', ec='k')
# text = ax.annotate('local max', xy=(3, 17), bbox=bbox_props)
#
# dts = []
# dt = DraggableText(text)
# dt.connect()
# dts.append(dt)
#
# plt.show()
