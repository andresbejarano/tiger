# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class RegularPolygonDomainForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Regular Polygon Domain")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # sides field
        self.__sides = QtWidgets.QSpinBox(self)
        self.__sides.setMinimum(3)
        self.__sides.setMaximum(1000)
        self.__sides.setSingleStep(1)
        self.__sides.setValue(3)
        formlayout.addRow("Sides:", self.__sides)
        
        # sidelength field
        self.__sidelength = QtWidgets.QDoubleSpinBox(self)
        self.__sidelength.setMinimum(0.0)
        self.__sidelength.setMaximum(1000.0)
        self.__sidelength.setSingleStep(0.1)
        self.__sidelength.setValue(1.0)
        self.__sidelength.setDecimals(3)
        formlayout.addRow("Side length:", self.__sidelength)
        
        # center X
        self.__centerx = QtWidgets.QDoubleSpinBox(self)
        self.__centerx.setMinimum(-1000.0)
        self.__centerx.setMaximum(1000.0)
        self.__centerx.setSingleStep(0.1)
        self.__centerx.setValue(0.0)
        self.__centerx.setDecimals(3)
        formlayout.addRow("Center X:", self.__centerx)
        
        # center Y
        self.__centery = QtWidgets.QDoubleSpinBox(self)
        self.__centery.setMinimum(-1000.0)
        self.__centery.setMaximum(1000.0)
        self.__centery.setSingleStep(0.1)
        self.__centery.setValue(0.0)
        self.__centery.setDecimals(3)
        formlayout.addRow("Center Y:", self.__centery)
        
        # center Z
        self.__centerz = QtWidgets.QDoubleSpinBox(self)
        self.__centerz.setMinimum(-1000.0)
        self.__centerz.setMaximum(1000.0)
        self.__centerz.setSingleStep(0.1)
        self.__centerz.setValue(0.0)
        self.__centerz.setDecimals(3)
        formlayout.addRow("Center Z:", self.__centerz)
        
        # plane field
        self.__plane = QtWidgets.QComboBox(self)
        self.__plane.addItems(['XY', 'XZ', 'YZ'])
        formlayout.addRow("Plane:", self.__plane)
        
        # Define the dialog buttons
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        dialoglayout.addWidget(buttons)
    
    
    def center(self) -> tuple:
        return (float(self.__centerx.value()), float(self.__centery.value()), float(self.__centerz.value()))
    
    
    def plane(self) -> str:
        return self.__plane.currentText()
    
    
    def sides(self) -> float:
        return int(self.__sides.value())
    
    
    def sidelength(self) -> float:
        return float(self.__sidelength.value())
    