from graphics.graphics_scene import GraphicsScene

DEBUG = True


class Scene:
    def __init__(self, main_win=None):
        self.main_win = main_win
        super().__init__()
        self.lines = []

        self.scene_width = 6400
        self.scene_height = 6400

        self.grScene = None
        self.cavityShape = None

        self.initUI()

    def initUI(self):
        # DRAW BACKGROUND
        self.grScene = GraphicsScene(self)
        self.grScene.setGrScene(self.scene_width, self.scene_height)

    def getItemAt(self, pos):
        return self.grScene.views()[0].itemAt(pos)

    def addLine(self, line):
        self.lines.append(line)

    def removeLine(self, line):
        self.lines.remove(line)

    def displayItemsOnScene(self):
        for i, line in enumerate(self.lines):
            print("\tLine {}: {}".format(i, line))

    def format(self, obj_list):
        new_obj_list = []
        for obj in obj_list:
            new_obj_list.append("<%s..%s>" % (hex(id(obj))[2:5], hex(id(obj))[-3:]))

        return new_obj_list

    def clear(self):
        while len(self.lines) > 0:
            self.lines[0].remove()
