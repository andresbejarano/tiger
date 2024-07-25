# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class TiltingAngleMethodForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Tilting Angle Method")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # direction field
        self.__direction = QtWidgets.QComboBox(self)
        self.__direction.addItems(['1', '-1'])
        formlayout.addRow("Direction:", self.__direction)
        
        # angle field
        self.__angle = QtWidgets.QDoubleSpinBox(self)
        self.__angle.setMinimum(0.0)
        self.__angle.setMaximum(1000.0)
        self.__angle.setSingleStep(0.1)
        self.__angle.setValue(1.0)
        self.__angle.setDecimals(3)
        formlayout.addRow("Angle:", self.__angle)
        
        # angle unit field
        self.__angleunit = QtWidgets.QComboBox(self)
        self.__angleunit.addItems(['degs', 'rads'])
        formlayout.addRow("Angle Unit:", self.__angleunit)
        
        # Define the dialog buttons
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        dialoglayout.addWidget(buttons)
    
    
    def angle(self) -> float:
        return float(self.__angle.value())
    
    def angleunit(self) -> str:
        return self.__angleunit.currentText()
    
    def direction(self) -> str:
        return self.__direction.currentText()
    