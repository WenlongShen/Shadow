# -*- coding: utf-8 -*-

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from GUI.annotationWidget import Ui_annoWidget

class Ui_mainWindow(QMainWindow):

	def __init__(self):
		super().__init__()
		self.initUI()

	def initUI(self):

		# mainWindow
		self.setGeometry(200, 200, 800, 600)
		self.setWindowTitle("Shadow")
		self.setWindowIcon(QIcon("Resources/shadow.ico"))
		self.statusBar()

		# menu
		openAction = QAction("Open", self)
		openAction.setShortcut("Ctrl+O")
		openAction.setStatusTip("Open")
		openAction.triggered.connect(self.showOpenDialog)
		saveAction = QAction("Save", self)
		saveAction.setShortcut("Ctrl+S")
		saveAction.setStatusTip("Save")
		saveAction.triggered.connect(self.showSaveDialog)
		exitAction = QAction("Exit", self)
		exitAction.setShortcut("Ctrl+Q")
		exitAction.setStatusTip("Exit")
		exitAction.triggered.connect(qApp.quit)

		annoAction = QAction("Annotation", self)
		annoAction.setStatusTip("Annotation")
		annoAction.triggered.connect(self.showAnnoWidget)

		docAction = QAction("Document", self)
		docAction.setStatusTip("Document")
		aboutAction = QAction("About Shadow", self)
		aboutAction.setStatusTip("About Shadow")
		aboutAction.triggered.connect(self.selectAbout)

		menubar = self.menuBar()
		fileMenu = menubar.addMenu("File")
		fileMenu.addAction(openAction)
		fileMenu.addAction(saveAction)
		fileMenu.addAction(exitAction)
		toolsMenu = menubar.addMenu("Tools")
		toolsMenu.addAction(annoAction)
		helpMenu = menubar.addMenu("Help")
		helpMenu.addAction(docAction)
		helpMenu.addAction(aboutAction)

		#central tab widget
		self.tabWidget = QTabWidget()
		self.setCentralWidget(self.tabWidget)
		self.tabWidget.setTabsClosable(True)

		self.textEdit = QTextEdit()
		self.tabWidget.addTab(self.textEdit, "File")
		self.show()

	def closeEvent(self, event):
		reply = QMessageBox.question(self, "Message", "You sure to quit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
		if reply == QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()

	def showOpenDialog(self):
		fileName = QFileDialog.getOpenFileName(self, "Open file", "./")
		try:
			file = open(fileName[0], "r")
			data = file.read()
		except UnicodeDecodeError:
			data = "Not a valid file..."
		except FileNotFoundError:
			data = "File not found..."
		except IndexError:
			data = ""
		self.textEdit.setText(data)

	def showSaveDialog(self):
		fileName = QFileDialog.getSaveFileName(self, "Save file", "./")
		try:
			file = open(fileName[0], "w")
			data = file.write(self.textEdit.toPlainText())
		except UnicodeDecodeError:
			data = "Not a valid file..."
		except FileNotFoundError:
			data = "File not found..."

	def showAnnoWidget(self):
		annoWidget = Ui_annoWidget()
		self.tabWidget.addTab(annoWidget, "Annotation")
		self.tabWidget.setCurrentIndex(self.tabWidget.currentIndex()+1)

	def selectAbout(self):
		QMessageBox.about(self, "About Shadow", "Copyright @ 2017 Wenlong Shen\n All Right reserved.")
