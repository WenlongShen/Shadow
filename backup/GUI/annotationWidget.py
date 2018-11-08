# -*- coding: utf-8 -*-

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

class Ui_annoWidget(QWidget):

	def __init__(self):
		super().__init__()
		self.initUI()

	def initUI(self):

		annoSiteButton = QPushButton("Sites", self)
		annoSiteButton.setStatusTip("Open sites file")
		annoSiteButton.clicked.connect(self.showAnnoSiteOpenDialog)
		self.annoSiteTextLine = QLineEdit()
		self.annoSiteTextLine.setReadOnly(True)
		
		annoGeneButton = QPushButton("Genes", self)
		annoGeneButton.setStatusTip("Open genes file")
		annoGeneButton.clicked.connect(self.showAnnoGeneOpenDialog)
		self.annoGeneTextLine = QLineEdit()
		self.annoGeneTextLine.setReadOnly(True)
		
		annoResultButton = QPushButton("Annotate", self)
		annoResultButton.clicked.connect(self.annotateProcess)
		self.annoResultTextLine = QLineEdit()
		annoSaveButton = QPushButton("Save", self)
		annoSaveButton.setStatusTip("Save result file")
		annoSaveButton.clicked.connect(self.showSaveDialog)

		self.textBrowser = QTextBrowser()

		gridLayout = QGridLayout()
		gridLayout.setSpacing(10)
		gridLayout.addWidget(annoSiteButton, 1, 1, 1, 1)
		gridLayout.addWidget(self.annoSiteTextLine, 1, 2, 1, 3)
		gridLayout.addWidget(annoGeneButton, 2, 1, 1, 1)
		gridLayout.addWidget(self.annoGeneTextLine, 2, 2, 1, 3)
		gridLayout.addWidget(annoResultButton, 3, 1, 1, 1)
		gridLayout.addWidget(self.annoResultTextLine, 3, 2, 1, 3)
		gridLayout.addWidget(annoSaveButton, 3, 5, 1, 1)
		gridLayout.addWidget(self.textBrowser, 4, 1, 2, 5)

		self.setLayout(gridLayout)


	def showAnnoSiteOpenDialog(self):
		fileName = QFileDialog.getOpenFileName(self, "Open file", "./")
		try:
			file = open(fileName[0], "r")
			data = file.readlines()[0]
		except UnicodeDecodeError:
			data = "Not a valid file..."
		except FileNotFoundError:
			data = "File not found..."
		except IndexError:
			data = ""
		self.textBrowser.setText(data)
		self.annoSiteTextLine.setText(fileName[0])

	def showAnnoGeneOpenDialog(self):
		fileName = QFileDialog.getOpenFileName(self, "Open file", "./")
		try:
			file = open(fileName[0], "r")
			data = file.readlines()[0]
		except UnicodeDecodeError:
			data = "Not a valid file..."
		except FileNotFoundError:
			data = "File not found..."
		except IndexError:
			data = ""
		self.textBrowser.setText(data)
		self.annoGeneTextLine.setText(fileName[0])

	def showSaveDialog(self):
		fileName = QFileDialog.getSaveFileName(self, "Save file", "./")
		try:
			file = open(fileName[0], "w")
			data = file.write(self.textBrowser.toPlainText())
		except UnicodeDecodeError:
			data = "Not a valid file..."
		except FileNotFoundError:
			data = "File not found..."

	def annotateProcess(self):
		return

