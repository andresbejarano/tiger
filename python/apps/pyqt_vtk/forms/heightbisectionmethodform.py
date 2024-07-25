# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class HeightBisectionMethodForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Height Bisection Method")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # direction field
        self.__direction = QtWidgets.QComboBox(self)
        self.__direction.addItems(['1', '-1'])
        formlayout.addRow("Direction:", self.__direction)
        
        # top height field
        self.__topheight = QtWidgets.QDoubleSpinBox(self)
        self.__topheight.setMinimum(0.0)
        self.__topheight.setMaximum(1000.0)
        self.__topheight.setSingleStep(0.1)
        self.__topheight.setValue(1.0)
        self.__topheight.setDecimals(3)
        formlayout.addRow("Top Height:", self.__topheight)
        
        # bottom height field
        self.__bottomheight = QtWidgets.QDoubleSpinBox(self)
        self.__bottomheight.setMinimum(0.0)
        self.__bottomheight.setMaximum(1000.0)
        self.__bottomheight.setSingleStep(0.1)
        self.__bottomheight.setValue(1.0)
        self.__bottomheight.setDecimals(3)
        formlayout.addRow("Bottom Height:", self.__bottomheight)
        
        # Define the dialog buttons
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        dialoglayout.addWidget(buttons)
    
    
    def bottomheight(self):
        return float(self.__bottomheight.value())
    
    def direction(self) -> str:
        return self.__direction.currentText()
    
    def topheight(self):
        return float(self.__topheight.value())
    