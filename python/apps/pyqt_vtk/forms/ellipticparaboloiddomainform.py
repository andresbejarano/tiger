# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets

class EllipticParaboloidDomainForm(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("Elliptic Paraboloid Domain")
        dialoglayout = QtWidgets.QVBoxLayout(self)
        groupbox = QtWidgets.QGroupBox("Parameters")
        dialoglayout.addWidget(groupbox)
        formlayout = QtWidgets.QFormLayout(groupbox)
        
        # a field
        self.__a = QtWidgets.QDoubleSpinBox(self)
        self.__a.setMinimum(-1000.0)
        self.__a.setMaximum(1000.0)
        self.__a.setSingleStep(0.1)
        self.__a.setValue(1.0)
        self.__a.setDecimals(3)
        formlayout.addRow("a:", self.__a)
        
        # b field
        self.__b = QtWidgets.QDoubleSpinBox(self)
        self.__b.setMinimum(-1000.0)
        self.__b.setMaximum(1000.0)
        self.__b.setSingleStep(0.1)
        self.__b.setValue(1.0)
        self.__b.setDecimals(3)
        formlayout.addRow("b:", self.__b)
        
        # width field
        self.__width = QtWidgets.QDoubleSpinBox(self)
        self.__width.setMinimum(0.0)
        self.__width.setMaximum(1000.0)
        self.__width.setSingleStep(0.1)
        self.__width.setValue(1.0)
        self.__width.setDecimals(3)
        formlayout.addRow("Width:", self.__width)
        
        # height field
        self.__height = QtWidgets.QDoubleSpinBox(self)
        self.__height.setMinimum(0.0)
        self.__height.setMaximum(1000.0)
        self.__height.setSingleStep(0.1)
        self.__height.setValue(1.0)
        self.__height.setDecimals(3)
        formlayout.addRow("Height:", self.__height)
        
        # width segments field
        self.__widthsegments = QtWidgets.QSpinBox(self)
        self.__widthsegments.setMinimum(1)
        self.__widthsegments.setMaximum(1000)
        self.__widthsegments.setSingleStep(1)
        self.__widthsegments.setValue(3)
        formlayout.addRow("Width Segments:", self.__widthsegments)
        
        # height segments field
        self.__heightsegments = QtWidgets.QSpinBox(self)
        self.__heightsegments.setMinimum(1)
        self.__heightsegments.setMaximum(1000)
        self.__heightsegments.setSingleStep(1)
        self.__heightsegments.setValue(3)
        formlayout.addRow("Height Segments:", self.__heightsegments)
        
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
    
    
    def a(self):
        return float(self.__a.value())
    
    
    def b(self):
        return float(self.__b.value())
    
    
    def center(self) -> tuple:
        return (float(self.__centerx.value()), float(self.__centery.value()), float(self.__centerz.value()))
        
    
    def height(self) -> float:
        return float(self.__height.value())
    
    
    def heightsegments(self) -> int:
        return int(self.__heightsegments.value())
    
    
    def plane(self) -> str:
        return self.__plane.currentText()
    
    
    def width(self) -> float:
        return float(self.__width.value())
    
    
    def widthsegments(self) -> int:
        return int(self.__widthsegments.value())
    