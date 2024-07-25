# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class ScaleDomainForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Scale Domain")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # factor field
        self.__factor = QtWidgets.QDoubleSpinBox(self)
        self.__factor.setMinimum(0.0)
        self.__factor.setMaximum(1000.0)
        self.__factor.setSingleStep(0.1)
        self.__factor.setValue(1.0)
        self.__factor.setDecimals(3)
        formlayout.addRow("Scale factor:", self.__factor)
        
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
        
        # Define the dialog buttons
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        dialoglayout.addWidget(buttons)
    
    
    def center(self) -> tuple:
        return (float(self.__centerx.value()), float(self.__centery.value()), float(self.__centerz.value()))
        
    
    def factor(self) -> float:
        return float(self.__factor.value())
    