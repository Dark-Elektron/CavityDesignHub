from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
# from QGraphics.node_graphics_cutline import QCutLine
# from QGraphics.node_graphics_edge import QGraphicsEdge
# from objects.node_edge import Edge, EDGE_TYPE_BEZIER, EDGE_TYPE_DIRECT
from graphics.object_line import Line
# from config import *

MODE_NOOP = 1
MODE_EDGE_DRAG = 2
MODE_EDGE_CUT = 3


EDGE_DRAG_START_THRESHOLD = 10

DEBUG = False
DEBUG_CONTEXT = False


class GraphicsView(QGraphicsView):
    def __init__(self, win, widg, parent=None):
        # self.scene = Scene()
        super().__init__(parent)
        self.win = win

        if widg == 'Wakefield':
            self.ui = win.ui
            self.geom_in = win.ui
        elif widg == "Eigenmode":
            self.ui = win.ui
            self.geom_in = win.geom_in
        else:
            self.ui = win.ui
            self.geom_in = win.ui

        self.scene = win.scene
        self.grScene = win.scene.grScene

        self.initUi()

        # This part ensures that the zoom of the viewport does not affect the scrollbar.
        # I still don't understand the logic but it works.... for now
        self.IN_VIEWPORT = False
        self.wE_original = self.ui.scrollArea.wheelEvent
        self.ui.scrollArea.wheelEvent = self.wE

        self.dragEdge = None # python wanted me to initialize it
        self.last_lmb_click_scene_pos = None # python wanted me to initialize it
        self.previousEdge = None # python wanted me to initialize it
        self.last_start_socket = None # python wanted me to initialize it

        self.setScene(self.grScene)

        self.mode = MODE_NOOP
        self.edge_color = QColor("#001000")
        self.socket_type = 1
        self.editingFlag = False

        self.zoomInFactor = 1.25
        self.zoomClamp = False
        self.zoom = 10
        # set intital scale

        self.zoomStep = 1
        self.zoomRange = [0, 10]

        # line list
        self.line_list = []

        # self.drawCells()
        self.update_signal()

    def initUi(self):
        # Housekeeping QGraphicsView
        self.setRenderHints(QPainter.Antialiasing | QPainter.HighQualityAntialiasing | QPainter.TextAntialiasing | QPainter.SmoothPixmapTransform)
        self.setViewportUpdateMode(QGraphicsView.FullViewportUpdate)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)

        # enable dropping
        self.setAcceptDrops(True)

        # set initial zoom scale
        self.scale(0.2, 0.2)

    def update_signal(self):
        self.geom_in.le_A_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_B_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_a_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_b_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_Ri_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_L_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_Req_i.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))

        self.geom_in.le_A_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_B_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_a_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_b_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_Ri_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_L_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_Req_ol.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))

        self.geom_in.le_A_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_B_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_a_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_b_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_Ri_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_L_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.le_Req_or.editingFinished.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))

        self.geom_in.sb_N_Cells.valueChanged.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))

        # beam pipe signals
        self.geom_in.cb_LBP.clicked.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.cb_RBP.clicked.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))

        # outer cells check box
        self.geom_in.cb_Outer_Cell_L.clicked.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))
        self.geom_in.cb_Outer_Cell_R.clicked.connect(lambda: self.drawCells(color=QColor(0, 0, 0, 255)))

    def drawCells(self, IC=None, OC=None, BP=None, color="#ffffff"):
        line_T = Line(self, line_type='line', section="top", IC=IC, OC=OC, BP=BP, color=color)
        line_B = Line(self, line_type='line', section="bottom", IC=IC, OC=OC, BP=BP, color=color)

        self.line_list.append(line_T)
        self.line_list.append(line_B)

    def removeCells(self):
        for line in self.line_list:
            line.remove()
        self.line_list.clear()

    # override required to enable drops
    def dragMoveEvent(self, event):
        event.acceptProposedAction()

    # def dragEnterEvent(self, event):
    #     if event.mimeData().hasFormat(LISTBOX_MIMETYPE):
    #         event.acceptProposedAction()
    #         # print("Graphics View:: GrView: entered")
    #         # print("Mime:: ", event.mimeData().text())

    def dragLeaveEvent(self, event):
        pass

    # def dropEvent(self, event):
    #     if event.mimeData().hasFormat(LISTBOX_MIMETYPE):
    #
    #         eventData = event.mimeData().data(LISTBOX_MIMETYPE)
    #         dataStream = QDataStream(eventData, QIODevice.ReadOnly)
    #
    #         pixmap = QPixmap()
    #         dataStream >> pixmap
    #         op_code = dataStream.readInt()
    #         node_type = dataStream.readInt()
    #         text = dataStream.readQString()
    #
    #         # print("op_code::", op_code, "text::", text)
    #         event.setDropAction(Qt.MoveAction)
    #         event.accept()
    #
    #         mouse_pos = event.pos()
    #         scene_pos = self.mapToScene(mouse_pos)
    #
    #         if node_type == 1:
    #             if op_code != 9:
    #                 self.addNode(scene_pos.x(), scene_pos.y(), op_code)
    #         else:
    #             self.addOther(scene_pos.x(), scene_pos.y(), node_type)

    def mousePressEvent(self, event):
        if event.button() == Qt.MiddleButton:
            self.middleMouseButtonPress(event)
        elif event.button() == Qt.LeftButton:
            self.leftMouseButtonPress(event)
        # elif event.button() == Qt.RightButton:
        #     self.rightMouseButtonPress(event)

        else:
            super().mousePressEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MiddleButton:
            self.middleMouseButtonRelease(event)
        elif event.button() == Qt.LeftButton:
            self.leftMouseButtonRelease(event)
        elif event.button() == Qt.RightButton:
            self.rightMouseButtonRelease(event)
        else:
            super().mouseReleaseEvent(event)

    def middleMouseButtonPress(self, event):
        releaseEvent = QMouseEvent(QEvent.MouseButtonRelease, event.localPos(), event.screenPos(), Qt.MiddleButton,
                                   Qt.NoButton, event.modifiers())
        super().mouseReleaseEvent(releaseEvent)
        self.setDragMode(QGraphicsView.ScrollHandDrag)

        fakeEvent = QMouseEvent(event.type(), event.localPos(), event.screenPos(), Qt.LeftButton, event.buttons()
                                | Qt.MiddleButton, event.modifiers())
        super().mousePressEvent(fakeEvent)

    def middleMouseButtonRelease(self, event):
        fakeEvent = QMouseEvent(event.type(), event.localPos(), event.screenPos(), Qt.LeftButton, event.buttons()
                                & -Qt.LeftButton, event.modifiers())
        super().mouseReleaseEvent(fakeEvent)
        self.setDragMode(QGraphicsView.NoDrag)

    def leftMouseButtonPress(self, event):
        super().mousePressEvent(event)

        item = self.getItemAtClick(event)

        self.last_lmb_click_scene_pos = self.mapToScene(event.pos())

        if DEBUG: print("LMB Click on", item, self.debug_modifiers(event))
        # print("Item type:", type(item))

        # Logic
        # if hasattr(item, "node") or isinstance(item, QGraphicsEdge) or item is None:
        #     self.setDragMode(QGraphicsView.NoDrag)
        #     if event.modifiers() and Qt.ShiftModifier:
        #         event.ignore()
        #         fakeEvent = QMouseEvent(QEvent.MouseButtonPress, event.localPos(), event.screenPos(),
        #                                 Qt.LeftButton, event.buttons() | Qt.LeftButton,
        #                                 event.modifiers() | Qt.ControlModifier)
        #         super().mousePressEvent(fakeEvent)
        #         return

        # if type(item) is QGraphicsSocket and item.socket in list(item.socket.node.outputs.values()):
        #     print("Node Graphics View:: ", item.socket)
        #     print("Node Graphics View:: ", item.socket.node)
        #     self.setDragMode(QGraphicsView.NoDrag)
        #     if self.mode == MODE_NOOP:
        #         self.mode = MODE_EDGE_DRAG
        #         self.edgeDragStart(item)
        #         return

        if self.mode == MODE_EDGE_DRAG:
            self.setDragMode(QGraphicsView.NoDrag)
            res = self.edgeDragEnd(item)
            if res: return

        QModifiers = QApplication.keyboardModifiers()
        if item is None:
            self.setDragMode(QGraphicsView.RubberBandDrag)
            if (QModifiers & Qt.ShiftModifier) == Qt.ShiftModifier:
                self.setDragMode(QGraphicsView.NoDrag)
                self.mode = MODE_EDGE_CUT
                fakeEvent = QMouseEvent(QEvent.MouseButtonRelease, event.localPos(), event.screenPos(),
                                        Qt.LeftButton, Qt.NoButton, event.modifiers())
                super().mouseReleaseEvent(fakeEvent)
                QApplication.setOverrideCursor(Qt.CrossCursor)
                return

        super().mousePressEvent(event)

    def leftMouseButtonRelease(self, event):
        # Get item clicked
        item = self.getItemAtClick(event)

        # logic
        # if hasattr(item, "node") or isinstance(item, QGraphicsEdge) or isinstance(item, QGraphicsEdge) or item is None:
        #     if event.modifiers() and Qt.ShiftModifier:
        #         event.ignore()
        #         fakeEvent = QMouseEvent(QEvent.MouseButtonRelease, event.localPos(), event.screenPos(),
        #                                 Qt.LeftButton, Qt.NoButton,
        #                                 event.modifiers() | Qt.ControlModifier)
        #         super().mouseReleaseEvent(fakeEvent)
        #         return

        if self.mode == MODE_EDGE_DRAG:
            if self.distanceBetweenCLickAndReleaseIsOff(event):
                res = self.edgeDragEnd(item)
                if res: return

        if self.mode == MODE_EDGE_CUT:
            self.cutIntersectingEdges()
            self.cutline.line_points = []
            self.cutline.update()
            QApplication.setOverrideCursor(Qt.ArrowCursor)
            self.mode = MODE_NOOP

        super().mouseReleaseEvent(event)

    # def rightMouseButtonPress(self, event):
    #     super().mousePressEvent(event)
    #
    #     item = self.getItemAtClick(event)

    def rightMouseButtonRelease(self, event):
        super().mouseReleaseEvent(event)

    def mouseMoveEvent(self, event):
        # update scene position
        mouse_pos = event.pos()
        scene_pos = self.mapToScene(mouse_pos)
        position = "{}/{}".format(scene_pos.x(), scene_pos.y())

        if self.mode == MODE_EDGE_DRAG:
            pos = self.mapToScene(event.pos())
            self.dragEdge.grEdge.setDestination(pos.x(), pos.y())
            self.dragEdge.grEdge.update()

        if self.mode == MODE_EDGE_CUT:
            pos = self.mapToScene(event.pos())
            self.cutline.line_points.append(pos)
            self.cutline.update()

        super().mouseMoveEvent(event)

    def mouseDoubleClickEvent(self, event):
        self.scene.displayDict()
        # item = self.getItemAtClick(event)
        # if isinstance(item, QGraphicsSocket):
        #     print()
        #     # print("Double clicked at ", item.x(), item.y())

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Delete:
            if not self.editingFlag:
                self.deleteSelected()
            else:
                super().keyPressEvent(event)

        elif event.key() == Qt.Key_S and event.modifiers() & Qt.ControlModifier:
            self.grScene.scene.saveToFile("graph.json.txt")

        elif event.key() == Qt.Key_L and event.modifiers() & Qt.ControlModifier:
            self.grScene.scene.loadFromFile("graph.json.txt")
        else:
            super().keyPressEvent(event)

    # def deleteSelected(self):
        # # try:
        # for item in self.grScene.selectedItems():
        #
        #     if isinstance(item, QGraphicsNode):
        #         item.node.remove()
        #     elif isinstance(item, QGraphicsEdge):
        #         item.edge.remove()
        # # except:
        #     print("NODE GRVIEW:: It all starts here")
    def enterEvent(self, event):
        self.IN_VIEWPORT = True
        super().enterEvent(event)

    def leaveEvent(self, event):
        self.IN_VIEWPORT = False
        super().leaveEvent(event)

    def wheelEvent(self, event):
        self.ui.scrollArea.wheelEvent = self.wE

        # calculate our zoom factor
        zoomOutFactor = 1 / self.zoomInFactor

        # calculate zoom
        if event.angleDelta().y() > 0:
            zoomFactor = self.zoomInFactor
            self.zoom += self.zoomStep
        else:
            zoomFactor = zoomOutFactor
            self.zoom -= self.zoomStep

        clamped = False
        if self.zoom < self.zoomRange[0]:
            self.zoom, clamped = self.zoomRange[0], True
        if self.zoom > self.zoomRange[1]:
            self.zoom, clamped = self.zoomRange[1], True

        # set scene scale
        if not clamped or self.zoomClamp is False:
            self.scale(zoomFactor, zoomFactor)

    def wE(self, ev):
        if ev.type() == QEvent.Wheel and self.IN_VIEWPORT:
            ev.ignore()
        else:
            self.ui.scrollArea.wheelEvent = self.wE_original

    def cutIntersectingEdgees(self):
        print()

    def debug_modifiers(self, event):
        out = "MODS: "
        if event.modifiers() and Qt.ShiftModifier: out += "SHIFT "
        elif event.modifiers() and Qt.ControlModifier: out += "CTRL "
        elif event.modifiers() and Qt.AltModifier: out += "ALT "
        return out

    def get_key_modifiers(self):
        QModifiers = QApplication.keyboardModifiers()
        modifiers = []
        if (QModifiers & Qt.ShiftModifier) == Qt.ShiftModifier:
            modifiers.append('shift')
        if (QModifiers & Qt.ControlModifier) == Qt.ControlModifier:
            modifiers.append('control')
        if (QModifiers & Qt.AltModifier) == Qt.AltModifier:
            modifiers.append('alt')
        return modifiers

    def getItemAtClick(self, event):
        """Return object clicked/release the mouse button on"""
        pos = event.pos()
        obj = self.itemAt(pos)
        return obj

    # def edgeDragStart(self, item):
    #     if DEBUG: print("View::edgeDragStart ~ Start dragging edge")
    #     if DEBUG: print("   assigning start socket...")
    #     self.previousEdge = item.socket.edge
    #     self.last_start_socket = item.socket
    #     self.edge_color = item.socket_color
    #     self.socket_type = item.socket_type
    #     self.dragEdge = Edge(self.grScene.scene, item.socket, None, EDGE_TYPE_BEZIER, item.socket_color)
    #     if DEBUG: print("View:edgeDragStart ~ dragEdge:", self.dragEdge)

    # def edgeDragEnd(self, item):
    #     """return true to skip the rest of code"""
    #     self.mode = MODE_NOOP
    #
    #     if type(item) is QGraphicsSocket:
    #         if item.socket != self.last_start_socket:
    #             # check if drop occurred in mode
    #             if item.socket_color == self.edge_color or item.socket_type == 5:
    #                 if DEBUG: print("View::edgeDragEnd ~ previous edge", self.previousEdge)
    #                 if item.socket.hasEdge():
    #                     item.socket.edge.remove()
    #                 if DEBUG: print("   Assigning end socket", item.socket)
    #                 if self.previousEdge is not None:
    #                     self.previousEdge.remove()
    #                 if DEBUG: print('View::edgeDragEnd ~ previous edge removed')
    #                 self.dragEdge.start_socket = self.last_start_socket
    #                 self.dragEdge.end_socket = item.socket
    #                 self.dragEdge.start_socket.setConnectedEdge(self.dragEdge)
    #                 self.dragEdge.end_socket.setConnectedEdge(self.dragEdge)
    #                 if DEBUG: print("View::edgeDragEnd ~ reassigned start and end sockets to drag edge")
    #                 self.dragEdge.updatePositions()
    #
    #                 # change edge type
    #                 if item.socket_type == 5:
    #                     item.socket.socket_type = self.dragEdge.start_socket.socket_type
    #                     end_socket = list(item.socket.node.outputs.values())[0]
    #                     end_socket.socket_type = self.dragEdge.start_socket.socket_type
    #                 return True
    #
    #     if DEBUG: print("End dragging edge")
    #     self.dragEdge.remove()
    #     self.dragEdge = None
    #     if DEBUG: print("View::edgeDragStart ~ about to set socket to previous edge: ", self.previousEdge)
    #     if self.previousEdge is not None:
    #         self.previousEdge.start_socket.edge = self.previousEdge
    #
    #     return False

    def distanceBetweenCLickAndReleaseIsOff(self, event):
        """This measures if we are far enough from last lmb click scene position"""
        new_lmb_release_scene_pos = self.mapToScene(event.pos())
        dist_scene = new_lmb_release_scene_pos - self.last_lmb_click_scene_pos
        return (dist_scene.x() * dist_scene.x() + dist_scene.y() * dist_scene.y()) > EDGE_DRAG_START_THRESHOLD**2

    # def contextMenuEvent(self, event):
    #     try:
    #         item = self.scene.getItemAt(event.pos())
    #
    #         # if DEBUG_CONTEXT: print(item)
    #         if type(item) == QGraphicsProxyWidget:
    #             item = item.widget()
    #
    #         if hasattr(item, 'node'):
    #             self.handleNodeContextMenu(event)
    #         elif hasattr(item, 'socket'):
    #             self.handleSocketContextMenu(event)
    #         elif hasattr(item, 'edge'):
    #             self.handleEdgeContextMenu(event)
    #         else:
    #             self.handleNewNodeContextMenu(event)
    #
    #         return super().contextMenuEvent(event)
    #     except Exception as e:
    #         print("Exception occurred: ", e)
    #
    # def handleNodeContextMenu(self, event):
    #     if DEBUG_CONTEXT: print("CONTEXT:: NODE")
    #
    #     context_menu = QMenu()
    #
    #     inputSocketAct = QMenu("Input Socket")
    #     inputSocket_typeAct = {}
    #     for i, in_socket_type in enumerate(PATHWAYS):
    #         inputSocket_typeAct[i] = inputSocketAct.addAction(in_socket_type)
    #
    #     context_menu.addMenu(inputSocketAct)
    #
    #     outputSocketAct = QMenu("Output Socket")
    #     outputSocket_typeAct = {}
    #     for i, out_socket_type in enumerate(PATHWAYS):
    #         outputSocket_typeAct[i] = outputSocketAct.addAction(out_socket_type)
    #
    #     context_menu.addMenu(outputSocketAct)
    #
    #     action = context_menu.exec_(self.mapToGlobal(event.pos()))
    #
    #     selected = None
    #     item = self.scene.getItemAt(event.pos()).parentItem()
    #     if hasattr(item, 'node'):
    #         selected = item
    #
    #     if DEBUG_CONTEXT: print("Node GraphicsView:: Attempting to add input socket")
    #     for key, itemAct in inputSocket_typeAct.items():
    #         if action == itemAct:
    #             selected.node.addInputSocket(key)
    #
    #     if DEBUG_CONTEXT: print("Node GraphicsView:: Added socket")
    #
    #     if DEBUG_CONTEXT: print("Node GraphicsView:: Attempting to add output socket")
    #     for key, itemAct in outputSocket_typeAct.items():
    #         if action == itemAct:
    #             selected.node.addOutputSocket(key)
    #
    #     if DEBUG_CONTEXT: print("Node GraphicsView:: Added socket")
    #
    # def handleSocketContextMenu(self, event):
    #     if DEBUG_CONTEXT: print("CONTEXT:: NODE")
    #
    #     context_menu = QMenu()
    #
    #     socketAct = context_menu.addAction("Remove Socket")
    #
    #     action = context_menu.exec_(self.mapToGlobal(event.pos()))
    #
    #     selected = None
    #     item = self.scene.getItemAt(event.pos())
    #     if hasattr(item, 'socket'):
    #         selected = item.parentItem()
    #
    #     if DEBUG_CONTEXT: print("Selected:: ", selected)
    #     if selected and action == socketAct:
    #         if DEBUG_CONTEXT: print("Node GraphicsView:: Attempting to remove socket")
    #         selected.node.removeSocket(item.socket)
    #         if DEBUG_CONTEXT: print("Node GraphicsView:: Remove socket")
    #
    # def handleEdgeContextMenu(self, event):
    #     if DEBUG_CONTEXT: print("CONTEXT:: EDGE")
    #     context_menu = QMenu()
    #     bezierAct = context_menu.addAction("Bezier Edge")
    #     directAct = context_menu.addAction("Direct Edge")
    #
    #     action = context_menu.exec_(self.mapToGlobal(event.pos()))
    #
    #     selected = None
    #     item = self.scene.getItemAt(event.pos())
    #
    #     if hasattr(item, 'edge'):
    #         selected = item.edge
    #
    #     if selected and action == bezierAct: selected.edge_type = EDGE_TYPE_BEZIER
    #     if selected and action == directAct: selected.edge_type = EDGE_TYPE_DIRECT
    #
    # def handleNewNodeContextMenu(self, event):
    #     if DEBUG_CONTEXT: print("CONTEXT:: EMPTY SPACE")
    #
    #     context_menu = QMenu()
    #
    #     structureAct = QMenu("Add Structure")
    #     struct_itemsActs = {}
    #     for i, struct_name in enumerate(STRUCTURES):
    #         struct_itemsActs[i] = structureAct.addAction(struct_name)
    #
    #     context_menu.addMenu(structureAct)
    #
    #     substructureAct = context_menu.addAction("Add Substructure")
    #     modeAct = context_menu.addAction("Add Mode")
    #     streamlinesAct = context_menu.addAction("Add Streamlines")
    #     functionalAct = context_menu.addAction("Add Functional")
    #
    #     utilityAct = context_menu.addAction("Utility")
    #
    #     action = context_menu.exec_(self.mapToGlobal(event.pos()))
    #
    #     scene_pos = self.mapToScene(event.pos())
    #
    #     x = scene_pos.x()
    #     y = scene_pos.y()
    #
    #     for key, itemAct in struct_itemsActs.items():
    #         if action == itemAct:
    #             self.addNode(x, y, key+1)
    #
    #     if action == substructureAct:
    #         self.addOther(x, y, SUBSTRUCTURES)
    #     if action == modeAct:
    #         self.addOther(x, y, MODE)
    #     if action == streamlinesAct:
    #         print("this")
    #         self.addOther(x, y, STREAMLINES)
    #     if action == functionalAct:
    #         self.addOther(x, y, FUNCTIONAL_STREAMLINES)
    #     if action == utilityAct:
    #         self.addOther(x, y, UTILITY)
