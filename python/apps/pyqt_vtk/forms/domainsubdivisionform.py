# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class DomainSubdivisionForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Domain Subdivision")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # type field
        self.__method = QtWidgets.QComboBox(self)
        self.__method.addItems(['Midpoint', 'Triangulate', 'Uniform'])
        formlayout.addRow("method:", self.__method)
        
        # Define the dialog buttons
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        dialoglayout.addWidget(buttons)
    
    
    def method(self) -> str:
        return self.__method.currentText()
    