# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class NormalizeDomainForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Normalize Domain")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # radius field
        self.__radius = QtWidgets.QDoubleSpinBox(self)
        self.__radius.setMinimum(0.0)
        self.__radius.setMaximum(1000.0)
        self.__radius.setSingleStep(0.1)
        self.__radius.setValue(1.0)
        self.__radius.setDecimals(3)
        formlayout.addRow("Radius:", self.__radius)
        
        # Define the dialog buttons
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        dialoglayout.addWidget(buttons)
        
    
    def radius(self) -> float:
        return float(self.__radius.value())
    