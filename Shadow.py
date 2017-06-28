# -*- coding: utf-8 -*-

import sys
from PyQt5.QtWidgets import QApplication
from GUI.mainWindow import Ui_mainWindow

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Ui_mainWindow()
    sys.exit(app.exec_())