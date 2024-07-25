# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class TorusDomainForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Torus Domain")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # major radius field
        self.__majorradius = QtWidgets.QDoubleSpinBox(self)
        self.__majorradius.setMinimum(0.0)
        self.__majorradius.setMaximum(1000.0)
        self.__majorradius.setSingleStep(0.1)
        self.__majorradius.setValue(3.0)
        self.__majorradius.setDecimals(3)
        formlayout.addRow("Major radius:", self.__majorradius)
        
        # minor radius field
        self.__minorradius = QtWidgets.QDoubleSpinBox(self)
        self.__minorradius.setMinimum(0.0)
        self.__minorradius.setMaximum(1000.0)
        self.__minorradius.setSingleStep(0.1)
        self.__minorradius.setValue(1.0)
        self.__minorradius.setDecimals(3)
        formlayout.addRow("Minor radius:", self.__minorradius)
        
        # major radial segments field
        self.__majorradialsegments = QtWidgets.QSpinBox(self)
        self.__majorradialsegments.setMinimum(1)
        self.__majorradialsegments.setMaximum(1000)
        self.__majorradialsegments.setSingleStep(1)
        self.__majorradialsegments.setValue(10)
        formlayout.addRow("Major radial segments:", self.__majorradialsegments)
        
        # minor radial segments field
        self.__minorradialsegments = QtWidgets.QSpinBox(self)
        self.__minorradialsegments.setMinimum(1)
        self.__minorradialsegments.setMaximum(1000)
        self.__minorradialsegments.setSingleStep(1)
        self.__minorradialsegments.setValue(10)
        formlayout.addRow("Minor radial Segments:", self.__minorradialsegments)
        
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
    
    
    def majorradialsegments(self) -> int:
        return int(self.__majorradialsegments.value())
    
    
    def majorradius(self) -> float:
        return float(self.__majorradius.value())
    
    
    def minorradialsegments(self) -> int:
        return int(self.__minorradialsegments.value())
    
    
    def minorradius(self) -> float:
        return float(self.__minorradius.value())
    
    
    def plane(self) -> str:
        return self.__plane.currentText()
    