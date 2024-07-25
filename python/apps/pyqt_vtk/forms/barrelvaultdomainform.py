# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class BarrelVaultDomainForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Barrel Vault Domain")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # length field
        self.__length = QtWidgets.QDoubleSpinBox(self)
        self.__length.setMinimum(0.0)
        self.__length.setMaximum(1000.0)
        self.__length.setSingleStep(0.1)
        self.__length.setValue(1.0)
        self.__length.setDecimals(3)
        formlayout.addRow("Length:", self.__length)
        
        # radius field
        self.__radius = QtWidgets.QDoubleSpinBox(self)
        self.__radius.setMinimum(0.0)
        self.__radius.setMaximum(1000.0)
        self.__radius.setSingleStep(0.1)
        self.__radius.setValue(1.0)
        self.__radius.setDecimals(3)
        formlayout.addRow("Radius:", self.__radius)
        
        # length segments field
        self.__lengthsegments = QtWidgets.QSpinBox(self)
        self.__lengthsegments.setMinimum(1)
        self.__lengthsegments.setMaximum(1000)
        self.__lengthsegments.setSingleStep(1)
        self.__lengthsegments.setValue(3)
        formlayout.addRow("Length Segments:", self.__lengthsegments)
        
        # radial segments field
        self.__radialsegments = QtWidgets.QSpinBox(self)
        self.__radialsegments.setMinimum(1)
        self.__radialsegments.setMaximum(1000)
        self.__radialsegments.setSingleStep(1)
        self.__radialsegments.setValue(3)
        formlayout.addRow("Radial Segments:", self.__radialsegments)
        
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
        
    
    def length(self) -> float:
        return float(self.__length.value())
    
    
    def lengthsegments(self) -> int:
        return int(self.__lengthsegments.value())
    
    
    def plane(self) -> str:
        return self.__plane.currentText()
    
    
    def radialsegments(self) -> int:
        return int(self.__radialsegments.value())
    
    
    def radius(self) -> float:
        return float(self.__radius.value())
    