import sys

from PyQt5 import QtWidgets
from GUI.mainWindow import Ui_MainWindow

class shadowMainWindow(QtWidgets.QWidget, Ui_MainWindow):
	def __init__(self):
		super(shadowMainWindow,self).__init__()
		self.setupUi(self)

if __name__=="__main__":
	app = QtWidgets.QApplication(sys.argv)
	Shadow = shadowMainWindow()
	Shadow.show()
	sys.exit(app.exec_())
