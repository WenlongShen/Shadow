# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Main.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(800, 600)

		self.menubar = QtWidgets.QMenuBar(MainWindow)
		self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
		self.menubar.setObjectName("menubar")
		self.menuFile = QtWidgets.QMenu(self.menubar)
		self.menuFile.setObjectName("menuFile")
		self.menuHelp = QtWidgets.QMenu(self.menubar)
		self.menuHelp.setObjectName("menuHelp")

		self.actionOpen = QtWidgets.QAction(MainWindow)
		self.actionOpen.setObjectName("actionOpen")
		self.actionExit = QtWidgets.QAction(MainWindow)
		self.actionExit.setObjectName("actionExit")
		self.actionAbout_Shadow = QtWidgets.QAction(MainWindow)
		self.actionAbout_Shadow.setObjectName("actionAbout_Shadow")

		self.menuFile.addAction(self.actionOpen)
		self.menuFile.addAction(self.actionExit)
		self.menuHelp.addAction(self.actionAbout_Shadow)
		self.menubar.addAction(self.menuFile.menuAction())
		self.menubar.addAction(self.menuHelp.menuAction())

		self.retranslateUi(MainWindow)
		QtCore.QMetaObject.connectSlotsByName(MainWindow)

	def retranslateUi(self, MainWindow):
		_translate = QtCore.QCoreApplication.translate
		MainWindow.setWindowTitle(_translate("MainWindow", "Shadow"))
		self.menuFile.setTitle(_translate("MainWindow", "File"))
		self.menuHelp.setTitle(_translate("MainWindow", "Help"))
		self.actionOpen.setText(_translate("MainWindow", "Open"))
		self.actionExit.setText(_translate("MainWindow", "Exit"))
		self.actionAbout_Shadow.setText(_translate("MainWindow", "About Shadow"))
