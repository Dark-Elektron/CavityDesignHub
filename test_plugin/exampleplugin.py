
#!/usr/bin/env python3

"""

An example of a widget as a Qt designer plugin. The widget plots magnetics IDS
data.

"""

from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtDesigner import QPyDesignerCustomWidgetPlugin

from exampleWidget import exampleWidget


class exampleplugin(QPyDesignerCustomWidgetPlugin):
    """Plugin for exampleWidget functionality.
    """
    def __init__(self, parent=None):
        super(exampleplugin, self).__init__(parent)

    def createWidget(self, parent):
        return exampleWidget(parent=parent, ids=None)

    def name(self):
        return "exampleWidget"

    def group(self):
        return "Custom Widgets"

    def icon(self):
        return QIcon(_logo_pixmap)

    def toolTip(self):
        return "Plot magnetics IDS data."

    def whatsThis(self):
        return ""

    def isContainer(self):
        return False

    def domXml(self):
        return '<widget class="exampleWidget" name="exampleWidget">\n</widget>'

    def includeFile(self):
        return "exampleWidget"


# Define the image used for the icon.
# Note: This is IMASviz default pixmap.
_logo_16x16_xpm = [
    "16 16 3 1 ",
    "  c black",
    ". c #0000E3",
    "X c None",
    "XXXXXXXXXXXXXXXX",
    "XXXXXXXXXXXXXXXX",
    "XXXXXXXXXXXXXXXX",
    " XX X  XX  XX  X",
    " XX    XX  XX XX",
    " XX  X X   XXX  ",
    " XX XX X XX X  X",
    "XXXXXXXXXXXXXXXX",
    "................",
    "XXXXXXXXXXXXXXXX",
    "XXX XX X X   XXX",
    "XXX X XX XX XXXX",
    "XXXX  XX X XXXXX",
    "XXXX  XX X   XXX",
    "XXXXXXXXXXXXXXXX",
    "XXXXXXXXXXXXXXXX"]

_logo_pixmap = QPixmap(_logo_16x16_xpm)