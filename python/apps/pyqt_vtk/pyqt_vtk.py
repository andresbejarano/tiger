# -*- coding: utf-8 -*-

import sys

from PyQt5 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from toolkit.pvector import PVector
from ds.vf import VF
import geometries.planar as planar
import geometries.platonicsolids as platonic
import geometries.closed as closed
from workspace.workspace import Workspace
from apps.vtk.vtkworkspace import VTKWorkspace
import toolkit.utils as utils


#
def planevectors(label:str):
    assert isinstance(label, str), "label not a string"
    assert label in ['XY', 'XZ', 'YZ'], "Wrong label value"
    if label == 'XY':
        return PVector(1), PVector(y=1)
    elif label == 'XZ':
        return PVector(1), PVector(z=1)
    elif label == 'YZ':
        return PVector(y=1), PVector(z=1)
    else:
        assert False, "Unexpected behavior"



class MainWindow(QtWidgets.QMainWindow):
    
    # Constructor of the class.
    def __init__(self, parent=None):
        
        # Initialize the main Qt window
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setWindowTitle("TIGER - Topological Interlocking GEneratoR")
        self.resize(QtWidgets.QDesktopWidget().availableGeometry(self).size() * 0.85)
        self.statusBar().setSizeGripEnabled(True)
        self.statusBar().show()
        
        # Main widget holds the workspace panel and the VTK render window
        self.centralwidget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.centralwidget)
        self.mainlayout = QtWidgets.QHBoxLayout(self.centralwidget)
        
        # Initialize the menu bar and its elements
        self.__initmenubar()
        
        # The object that handles everything related to TIC generation
        self.workspace = Workspace()
        
        # The object that handles everything related to VTK (except the widget)
        self.vtkworkspace = VTKWorkspace(self.width(), self.height())
        
        # Set up the VTK render window
        self.vtkwidget = QVTKRenderWindowInteractor(self.centralwidget)
        self.vtkwidget.Initialize()
        self.vtkwidget.Start()
        self.vtkwidget.GetRenderWindow().AddRenderer(self.vtkworkspace.renderer)
        self.vtkwidget.GetRenderWindow().Render()
        
        self.mainlayout.addWidget(self.vtkwidget)
        self.show()
    
    
    # Asks the user if they want to load a new geometric domain into the 
    # workspace. If the domain in the workspace is none then its returns True.
    def __confirmnewdomain(self):
        assert self.workspace is not None, "Missing workspace"
        assert hasattr(self.workspace, 'domain'), "workspace has no domain"
        
        # Return True if there is no domain in the workspace (no need of 
        # confirmation)
        if self.workspace.domain is None:
            return True
        
        # Warn the user that a new domain will clear the workspace
        reply = QtWidgets.QMessageBox.question(self, 
                                     'Confirm Action', 
                                     'Setting a new geometric domain clears all existing elements. Do you want to continue?',
                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, 
                                     QtWidgets.QMessageBox.No)
        return reply == QtWidgets.QMessageBox.Yes
    
    # Initializes the submenus and action of the assembly menu
    def __initassemblymenu(self, menubar):
        assemblymenu = menubar.addMenu('Assembly')
        
        # Define the Assembly->Generate menu
        generatemenu = assemblymenu.addMenu('Generate')
        tiltingangleaction = generatemenu.addAction('Tilting Angle Method')
        tiltingangleaction.setStatusTip('Generate the assembly using the Tilting Angle method')
        tiltingangleaction.triggered.connect(self.__ontiltingangleactiontrigger)
        heightbisectionaction = generatemenu.addAction('Height Bisection Method')
        heightbisectionaction.setStatusTip('Generate the assembly using the Height Bisection method')
        heightbisectionaction.triggered.connect(self.__onheightbisectionaction)
    
    
    # Initializes the submenus and actions of the Domain menu
    def __initdomainmenu(self, menubar):
        domainmenu = menubar.addMenu('Domain')
        
        # Define the Domain->Geometries menu
        geometriesmenu = domainmenu.addMenu("Geometries")
        closedgeometriesmenu = geometriesmenu.addMenu("Closed")
        torusaction = closedgeometriesmenu.addAction("Torus")
        torusaction.setStatusTip("Set a torus as the geometric domain")
        torusaction.triggered.connect(self.__ontorusactiontrigger)
        planargeometriesmenu = geometriesmenu.addMenu("Planar")
        barrelvaultaction = planargeometriesmenu.addAction("Barrel Vault")
        barrelvaultaction.setStatusTip("Set a barrel vault as the geometric domain")
        barrelvaultaction.triggered.connect(self.__onbarrelvaultactiontrigger)
        ellipticparaboloidaction = planargeometriesmenu.addAction("Elliptic Paraboloid")
        ellipticparaboloidaction.setStatusTip('Set an elliptic paraboloid as the geometric domain')
        ellipticparaboloidaction.triggered.connect(self.__onellipticparaboloidactiontrigger)
        equilateraltriangleaction = planargeometriesmenu.addAction("Equilateral Triangle")
        equilateraltriangleaction.setStatusTip("Set an equilateral triangle as the geometric domain")
        equilateraltriangleaction.triggered.connect(self.__onequilateraltriangleactiontrigger)
        monkeysaddleaction = planargeometriesmenu.addAction("Monkey Saddle")
        monkeysaddleaction.setStatusTip("Set a monkey saddle as the geometric domain")
        monkeysaddleaction.triggered.connect(self.__onmonkeysaddleactiontrigger)
        regularpolygonaction = planargeometriesmenu.addAction("Regular Polygon")
        regularpolygonaction.setStatusTip('Set a regular polygon as the geometric domain')
        regularpolygonaction.triggered.connect(self.__onregularpolygonactiontrigger)
        saddleaction = planargeometriesmenu.addAction("Saddle")
        saddleaction.setStatusTip("Set a saddle as the geometric domain")
        saddleaction.triggered.connect(self.__onsaddleactiontrigger)
        squareaction = planargeometriesmenu.addAction("Square")
        squareaction.setStatusTip("Set a square as the geometric domain")
        squareaction.triggered.connect(self.__onsquareactiontrigger)
        squaredgridaction = planargeometriesmenu.addAction("Squared Grid")
        squaredgridaction.setStatusTip("Set a squared grid as the geometric domain")
        squaredgridaction.triggered.connect(self.__onsquaredgridactiontrigger)
        waveaction = planargeometriesmenu.addAction("Wave")
        waveaction.setStatusTip("Set a wave as the geometric domain")
        waveaction.triggered.connect(self.__onwaveactiontrigger)
        polyhedramenu = geometriesmenu.addMenu("Polyhedra")
        platonicsolidaction = polyhedramenu.addAction("Platonic Solid")
        platonicsolidaction.setStatusTip('Set a Platonic solid as the geometric domain')
        platonicsolidaction.triggered.connect(self.__onplatonicsolidactiontrigger)
        
        # Define the Domain->Edit menu
        editdomainmenu = domainmenu.addMenu('Edit')
        dualaction = editdomainmenu.addAction('Dual')
        dualaction.setStatusTip('Set the dual as the geometric domain')
        dualaction.triggered.connect(self.__ondualactiontrigger)
        flipaction = editdomainmenu.addAction('Flip')
        flipaction.setStatusTip('Flip the faces of the geometric domain')
        flipaction.triggered.connect(self.__onflipactiontrigger)
        normalizeaction = editdomainmenu.addAction('Normalize')
        normalizeaction.setStatusTip('Normalize the vertices of the geometric domain')
        normalizeaction.triggered.connect(self.__onnormalizeactiontrigger)
        scaleaction = editdomainmenu.addAction('Scale')
        scaleaction.setStatusTip('Scale the geometric domain')
        scaleaction.triggered.connect(self.__onscaleactiontrigger)
        subdivideaction = editdomainmenu.addAction('Subdivide')
        subdivideaction.setStatusTip('Subdivide the faces of the geometric domain')
        subdivideaction.triggered.connect(self.__onsubdivideactiontrigger)
    
    
    # Initializes the submenus and actions of the File menu
    def __initfilemenu(self, menubar):
        filemenu = menubar.addMenu('File')
        exitAction = filemenu.addAction("Exit")
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(QtWidgets.qApp.quit)
    
    
    # Initializes the menu bar
    def __initmenubar(self):
        menubar = self.menuBar()
        self.__initfilemenu(menubar)
        self.__initdomainmenu(menubar)
        self.__initassemblymenu(menubar)
        #TODO self.__initanalysismenu(menubar)
        self.__initviewmenu(menubar)
    
    
    # Initializes the submenus and actions of the View menu
    def __initviewmenu(self, menubar):
        viewmenu = menubar.addMenu("View")
        colorsmenu = viewmenu.addMenu("Colors")
        backgroundcoloraction = colorsmenu.addAction("Background")
        backgroundcoloraction.setStatusTip("Set the background color")
        backgroundcoloraction.triggered.connect(self.__onbackgroundcoloractiontrigger)
        togglemenu = viewmenu.addMenu("Toggle")
        self.viewaxesaction = togglemenu.addAction("Axes")
        self.viewaxesaction.setCheckable(True)
        self.viewaxesaction.setChecked(True)
        self.viewaxesaction.setStatusTip("Toggle the 3D axes viewing")
        self.viewaxesaction.triggered.connect(self.__onviewaxesactiontrigger)
        self.viewxyplaneaction = togglemenu.addAction("Plane")
        self.viewxyplaneaction.setCheckable(True)
        self.viewxyplaneaction.setChecked(True)
        self.viewxyplaneaction.setStatusTip("Toggle the XY plane grid viewing")
        self.viewxyplaneaction.triggered.connect(self.__onviewxyplanetrigger)
    
    
    # Clears the information of the geometric domain, both workspace and 
    # and loads them with the given geometry.
    def __loaddomain(self, vf:VF):
        assert vf is not None, "Missing vf"
        assert isinstance(vf, VF), "vf not a VF"
        assert self.workspace is not None, "Missing workspace"
        assert self.vtkworkspace is not None, "Missing VTK workspace"
        
        # Sets the geometric domain to the workspace. Doing this clears all of
        # the information currently stored in the workspace
        self.workspace.setdomain(vf)
        self.vtkworkspace.loaddomainactors(self.workspace.domain)
        
        # Force a rendering
        self.vtkwidget.GetRenderWindow().Render()
    
    
    # Function called when the user clicks the View->Colors->Background menu action. 
    def __onbackgroundcoloractiontrigger(self):
        color = QtWidgets.QColorDialog.getColor()
        if color.isValid():
            rgb = color.getRgbF()
            self.vtkworkspace.renderer.SetBackground(rgb[0], rgb[1], rgb[2])
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Barrel Vault menu action.
    def __onbarrelvaultactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.barrelvaultdomainform as form
        dialog = form.BarrelVaultDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        length = dialog.length()
        radius = dialog.radius()
        ls = dialog.lengthsegments()
        rs = dialog.radialsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.barrelvault(length, radius, ls, rs, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Edit->Dual menu action.
    def __ondualactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Warn the user in case of no geometric domain in the workspace
        if self.workspace.domain is None:
            QtWidgets.QMessageBox.critical(self, "No Geometric Domain", "There is no geometric domain in the workspace", QtWidgets.QMessageBox().Ok)
            return
        
        vf = self.workspace.domain.dual()
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Elliptic Paraboloid menu action.
    def __onellipticparaboloidactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.ellipticparaboloiddomainform as form
        dialog = form.EllipticParaboloidDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        a = dialog.a()
        b = dialog.b()
        width = dialog.width()
        height = dialog.height()
        ws = dialog.widthsegments()
        hs = dialog.heightsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.ellipticparaboloid(a, b, width, height, ws, hs, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Equilateral Triangle menu action.
    def __onequilateraltriangleactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.equilateraltriangledomainform as form
        dialog = form.EquilateralTriangleDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        length = dialog.length()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.equilateraltriangle(length, C, U, V)
        self.__loaddomain(vf)
    
    
    def __onflipactiontrigger(self):
        QtWidgets.QMessageBox.information(self, "Title", "Flip", QtWidgets.QMessageBox().Ok)
        #TODO
    
    
    def __onheightbisectionaction(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Warn the user in case of no geometric domain in the workspace
        if self.workspace.domain is None:
            QtWidgets.QMessageBox.critical(self, "No Geometric Domain", "There is no geometric domain in the workspace", QtWidgets.QMessageBox().Ok)
            return
        
        if not self.workspace.domain.arefacesevensided():
            QtWidgets.QMessageBox.critical(self, "No Suitable Geometric Domain", "At least one of the faces of the geometric domain does not have an even number of sides", QtWidgets.QMessageBox().Ok)
            return
        
        # Open the dialog for the assembly generation parameters. Exit the 
        # function if the users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.heightbisectionmethodform as form
        dialog = form.HeightBisectionMethodForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        direction = int(dialog.direction())
        topheight = dialog.topheight()
        bottomheight = dialog.bottomheight()
        
        # Run the steps for creating the TIC using the Height Bisection method
        self.workspace.setedgemidpoints()
        self.workspace.setedgedirections(direction)
        self.workspace.setfaceheights(topheight, bottomheight)
        self.workspace.setheightbisectionelements()
        self.workspace.setedgeplanes()
        self.workspace.calculateblocks()
        
        self.vtkworkspace.loadassemblyactors(self.workspace.assembly)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Monkey Saddle menu action.
    def __onmonkeysaddleactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.monkeysaddledomainform as form
        dialog = form.MonkeySaddleDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        width = dialog.width()
        height = dialog.height()
        ws = dialog.widthsegments()
        hs = dialog.heightsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.monkeysaddle(width, height, ws, hs, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Edit->Normalize menu action.
    def __onnormalizeactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Warn the user in case of no geometric domain in the workspace
        if self.workspace.domain is None:
            QtWidgets.QMessageBox.critical(self, "No Geometric Domain", "There is no geometric domain in the workspace", QtWidgets.QMessageBox().Ok)
            return
        
        # Ask the user to confirm the normalization
        question = QtWidgets.QMessageBox.question(self, "Normalize Domain", "Normalizing the domain vertices may result in unexpected results depending on the domain. Do you wish to continue?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
        if question == QtWidgets.QMessageBox.No:
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.normalizedomainform as form
        dialog = form.NormalizeDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        radius = dialog.radius()
        
        self.workspace.domain.normalizevertices(radius)
        self.vtkworkspace.loaddomainactors(self.workspace.domain)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Polyhedra->Platonic Solid menu action.
    def __onplatonicsolidactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.platonicsoliddomainform as form
        dialog = form.PlatonicSolidDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        solid = dialog.solid()
        radius = dialog.radius()
        
        # Get the required geometric domain and load it
        vf = platonic.byname(solid, radius)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Regular Polygon menu action.
    def __onregularpolygonactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.regularpolygondomainform as form
        dialog = form.RegularPolygonDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        sides = dialog.sides()
        sidelength = dialog.sidelength()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.regularpolygon(sides, sidelength, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Saddle menu action.
    def __onsaddleactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.saddledomainform as form
        dialog = form.SaddleDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        width = dialog.width()
        height = dialog.height()
        ws = dialog.widthsegments()
        hs = dialog.heightsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.saddle(width, height, ws, hs, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Edit->Scale menu action.
    def __onscaleactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Warn the user in case of no geometric domain in the workspace
        if self.workspace.domain is None:
            QtWidgets.QMessageBox.critical(self, "No Geometric Domain", "There is no geometric domain in the workspace", QtWidgets.QMessageBox().Ok)
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.scaledomainform as form
        dialog = form.ScaleDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected parameters from the dialog
        s = dialog.factor()
        center = dialog.center()
        C = PVector(center[0], center[1], center[2])
        
        # Edit the domain and load the respective actors
        self.workspace.domain.scale(s, C)
        self.vtkworkspace.loaddomainactors(self.workspace.domain)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Square menu action.
    def __onsquareactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.squaredomainform as form
        dialog = form.SquareDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        length = dialog.length()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.square(length, C, U, V)
        self.__loaddomain(vf)
        
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Squared Grid menu action.
    def __onsquaredgridactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.squaredgriddomainform as form
        dialog = form.SquaredGridDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        width = dialog.width()
        height = dialog.height()
        ws = dialog.widthsegments()
        hs = dialog.heightsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.squaredgrid(width, height, ws, hs, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Domain->Edit->Subdivide menu action.
    def __onsubdivideactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Warn the user in case of no geometric domain in the workspace
        if self.workspace.domain is None:
            QtWidgets.QMessageBox.critical(self, "No Geometric Domain", "There is no geometric domain in the workspace", QtWidgets.QMessageBox().Ok)
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.domainsubdivisionform as form
        dialog = form.DomainSubdivisionForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        method = dialog.method().lower()
        
        if method == "midpoint":
            vf = self.workspace.domain.facemidpointsubdivision()
        elif method == "triangulate":
            vf = self.workspace.domain.facetriangulatesubdivision()
        elif method == 'uniform':
            vf = self.workspace.domain.faceuniformsubdivision()
        else:
            assert False, "Unexpected value"
        
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the 
    # Assembly->Generate->Tilting Angle Method menu action.
    def __ontiltingangleactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Warn the user in case of no geometric domain in the workspace
        if self.workspace.domain is None:
            QtWidgets.QMessageBox.critical(self, "No Geometric Domain", "There is no geometric domain in the workspace", QtWidgets.QMessageBox().Ok)
            return
        
        if not self.workspace.domain.arefacesevensided():
            QtWidgets.QMessageBox.critical(self, "No Suitable Geometric Domain", "At least one of the faces of the geometric domain does not have an even number of sides", QtWidgets.QMessageBox().Ok)
            return
        
        # Open the dialog for the assembly generation parameters. Exit the 
        # function if the users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.tiltinganglemethodform as form
        dialog = form.TiltingAngleMethodForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        direction = int(dialog.direction())
        angle = dialog.angle()
        if dialog.angleunit() == "degs":
            angle = utils.toradians(angle)
        
        # Run the steps for creating the TIC using the Tilting Angle method
        self.workspace.setedgemidpoints()
        self.workspace.setedgedirections(direction)
        self.workspace.setedgevectorsfromnormals()
        self.workspace.setedgerotationangles(angle)
        self.workspace.rotateedgevectors()
        self.workspace.setedgeplanes()
        self.workspace.calculateblocks()
        
        self.vtkworkspace.loadassemblyactors(self.workspace.assembly)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Closed->Torus menu action.
    def __ontorusactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.torusdomainform as form
        dialog = form.TorusDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        R = dialog.majorradius()
        r = dialog.minorradius()
        Rs = dialog.majorradialsegments()
        rs = dialog.minorradialsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = closed.torus(R, r, Rs, rs, C, U, V)
        self.__loaddomain(vf)
    
    
    # Function called when the user clicks the View->Plane menu action.
    def __onviewxyplanetrigger(self):
        assert self.vtkworkspace.xyplaneactor is not None, "Missing XY plane actor"
        viewplane = self.viewxyplaneaction.isChecked()
        self.vtkworkspace.xyplaneactor.SetVisibility(viewplane)
    
    
    # Function called when the user clicks the View->Axes menu action.
    def __onviewaxesactiontrigger(self):
        assert self.vtkworkspace.axesactor is not None, "Missing axes actor"
        viewaxes = self.viewaxesaction.isChecked()
        self.vtkworkspace.axesactor.SetVisibility(viewaxes)
    
    
    # Function called when the user clicks the 
    # Domain->Geometries->Planar->Wave menu action.
    def __onwaveactiontrigger(self):
        assert self.workspace is not None, "Missing workspace"
        
        # Ask the user to confirm (if required) the new goemetric domain in the
        # workspace
        if not self.__confirmnewdomain():
            return
        
        # Open the dialog for the domain parameters. Exit the function if the 
        # users rejects (i.e., cancels or closes) the dialog
        import apps.pyqt_vtk.forms.wavedomainform as form
        dialog = form.WaveDomainForm(self)
        if dialog.exec() == QtWidgets.QDialog.Rejected:
            return
        
        # Get the selected domain parameters from the dialog
        width = dialog.width()
        height = dialog.height()
        ws = dialog.widthsegments()
        hs = dialog.heightsegments()
        center = dialog.center()
        plane = dialog.plane()
        C = PVector(center[0], center[1], center[2])
        U, V = planevectors(plane)
        
        # Get the required geometric domain and load it
        vf = planar.wave(width, height, ws, hs, C, U, V)
        self.__loaddomain(vf)        


def run():
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())


if __name__ == "__main__":
    run()
    