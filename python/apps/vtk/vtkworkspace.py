# -*- coding: utf-8 -*-

from ds.dcel import DCEL
import vtk


# Loads the geometry of a DCEL object into two vtkActor objects (one for 
# triangular faces, and one for the edges). The reason we need two objects is 
# to maintain the logical correspondence between the polygonal nature of the
# geometry faces while drawing the edges. If we use a single object for both
# purposes, the drawn edges will correspond to the triangulation of the faces.
def loadvtkactorswithdcel(dcel:DCEL):
    assert dcel is not None, "Missing dcel"
    assert isinstance(dcel, DCEL), "dcel not a DCEL"
    
    # A dictionary to map a vertex coordinate to its index in the vertices
    # array of the geometry
    vToI = {}
    
    # A set to store the (src, dst) tuples representing the edges of the 
    # geometry.
    edges = set()

    # Define the points of the geometry
    nV = dcel.numvertices()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nV)
    
    for i in range(nV):
        L = dcel.V[i].coords.tolist()
        points.InsertPoint(i, L)
        vToI[str(L)] = i
    
    # Define the triangles representing the faces of the geometry
    nF = dcel.numfaces()
    triangles = vtk.vtkCellArray()
    triangles.SetNumberOfCells(dcel.numtriangles())
    for i in range(nF):
        
        # The start vertex of the incident half edge of the face is the common
        # vertex for all triangles representing the current face.
        v0 = vToI[str(dcel.F[i].halfedge.start.coords.tolist())]
        
        # Walk through the half edges of the current face and define its
        # triangles
        currentHalfedge = dcel.F[i].halfedge.next
        v1 = vToI[str(currentHalfedge.start.coords.tolist())]
        
        # Define the information of the first edge of the current face
        edges.add((min(v0, v1), max(v0, v1)))
        
        while True:
            
            v2 = vToI[str(currentHalfedge.next.start.coords.tolist())]
            
            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, v0)
            triangle.GetPointIds().SetId(1, v1)
            triangle.GetPointIds().SetId(2, v2)
            triangles.InsertNextCell(triangle)
            
            #edges.add((min(v0, v1), max(v0, v1)))
            edges.add((min(v1, v2), max(v1, v2)))
            
            # Move to the next half edge of the face. Break the cycle if
            # necessary
            v1 = v2
            currentHalfedge = currentHalfedge.next
            
            if currentHalfedge == dcel.F[i].halfedge.previous:
                break
        
        # Finish the last edge of the current face
        edges.add((min(v0, v1), max(v0, v1)))
    
    # Define the polydata object that contains the information of the faces of
    # the geometric object to be rendered as an actor
    facesdata = vtk.vtkPolyData()
    facesdata.SetPoints(points)
    facesdata.SetPolys(triangles)
    
    # Define polydatanormals object to calculate the normals of the points
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(facesdata)
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.Update()
    
    # Create a mapper and actor
    facesmapper = vtk.vtkPolyDataMapper()
    facesmapper.SetInputConnection(normals.GetOutputPort())
    
    dcelactor = vtk.vtkActor()
    dcelactor.SetMapper(facesmapper)
    
    # Populate a vtkCellArray object with the information of the edges
    linesegments = vtk.vtkCellArray()
    linesegments.SetNumberOfCells(len(edges))
    for edge in edges:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, edge[0])
        line.GetPointIds().SetId(1, edge[1])
        linesegments.InsertNextCell(line)
    
    edgesdata = vtk.vtkPolyData()
    edgesdata.SetPoints(points)
    edgesdata.SetLines(linesegments)
    
    edgesmapper = vtk.vtkPolyDataMapper()
    edgesmapper.SetInputData(edgesdata)
    
    edgesactor = vtk.vtkActor()
    edgesactor.SetMapper(edgesmapper)
    
    return dcelactor, edgesactor



class VTKWorkspace:
    
    # Constructor of the class.
    def __init__(self, width:int, height:int):
        assert isinstance(width, int), "width is not an integer"
        assert width > 0, "width not positive"
        assert isinstance(height, int), "height not an integer"
        assert height > 0, "height not positive"
        
        # The VTK actors representing the geometric domain
        self.domaindcelactor = None
        self.domainedgesdactor = None
        
        # The actors representing the assembly
        self.assemblyactor = None
        self.assemblyedgesactor = None
        
        # Initialize the renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.1, 0.2, 0.3)
        
        # Initialize the XY plane actor
        xyplane = vtk.vtkPlaneSource()
        xyplane.SetXResolution(30)  # Number of divisions along X
        xyplane.SetYResolution(30)  # Number of divisions along Y
        xyplane.SetOrigin(-3, -3, 0)  # Set the origin of the plane
        xyplane.SetPoint1(3, -3, 0)  # Set the endpoint of X axis on the plane
        xyplane.SetPoint2(-3, 3, 0)  # Set the endpoint of Y axis on the plane
        xyplanemapper = vtk.vtkPolyDataMapper()
        xyplanemapper.SetInputConnection(xyplane.GetOutputPort())
        self.xyplaneactor = vtk.vtkActor()
        self.xyplaneactor.SetMapper(xyplanemapper)
        self.xyplaneactor.GetProperty().SetRepresentationToWireframe()  # Set as wireframe
        self.xyplaneactor.GetProperty().SetColor(0.8, 0.8, 0.8)  # Color the grid
        self.renderer.AddActor(self.xyplaneactor)
        
        # Create the axes actor
        self.axesactor = vtk.vtkAxesActor()
        self.axesactor.SetTotalLength(1, 1, 1)  # Set the length of the axes
        self.axesactor.SetAxisLabels(1)  # Enable axis labels
        self.renderer.AddActor(self.axesactor)
        
        # Set up the camera
        self.camera = vtk.vtkCamera()
        self.camera.SetPosition(8, 8, 8)
        self.camera.SetFocalPoint(0, 0, 0)
        self.camera.SetViewUp(0, 0, 1)
        self.renderer.SetActiveCamera(self.camera)
        
        # Set the SSAO pass
        scenediagonal = (width ** 2 + height ** 2) ** 0.5
        ssaopass = vtk.vtkSSAOPass()
        ssaopass.SetRadius(0.1 * scenediagonal)
        ssaopass.SetBias(0.001 * scenediagonal)
        ssaopass.SetKernelSize(128)
        ssaopass.BlurOff()
        basicpasses = vtk.vtkRenderStepsPass()
        ssaopass.SetDelegatePass(basicpasses)
        self.renderer.SetPass(ssaopass)
    
    
    # Removes all TIC related actos from the scene.
    def clearall(self):
        
        if self.domaindcelactor is not None:
            self.renderer.RemoveActor(self.domaindcelactor)
            self.renderer.RemoveActor(self.domainedgesdactor)
            del self.domaindcelactor
            del self.domainedgesdactor
            self.domaindcelactor = None
            self.domainedgesdactor = None
        
        if self.assemblyactor is not None:
            self.renderer.RemoveActor(self.assemblyactor)
            self.renderer.RemoveActor(self.assemblyedgesactor)
            del self.assemblyactor
            del self.assemblyedgesactor
            self.assemblyactor = None
            self.assemblyedgesactor = None
    
    
    # 
    def loadassemblyactors(self, assembly):
        assert assembly is not None, "missing assembly"
        
        # Remove the assembly actors from the scene
        if self.assemblyactor is not None:
            self.renderer.RemoveActor(self.assemblyactor)
            self.renderer.RemoveActor(self.assemblyedgesactor)
            del self.assemblyactor
            del self.assemblyedgesactor
        
        pointscount = 0
        
        points = vtk.vtkPoints()
        edges = vtk.vtkCellArray()
        triangles = vtk.vtkCellArray()
        
        # Load the new assembly actors
        for f in assembly:
            vertices, faces = assembly[f]
            
            # If the block vertices is None then continue to the next face. 
            # This could happen if the respective edge planes did not intersect
            if faces is None:
                continue
            
            # Traverse through the vertices of the current block and add them 
            # to the points object.
            nvertices = len(vertices)
            for i in range(nvertices):
                points.InsertPoint(pointscount + i, vertices[i].tolist())
            
            # A set to store the tuples representing the edges of the block
            visitededges = set()
            
            # Traverse through the faces of the current block, triangulate them
            # and add them to the triangles object.
            for face in faces:
                
                # Get the number of indices required to define the current 
                # block face
                nindices = len(face)
                
                # Traverse through the face indices, define the triangles, and 
                # add them to the triangles object. All triangles are defined 
                # in terms of the first vertex index of the face.
                for vIdx in range(1, nindices - 1):
                    triangle = vtk.vtkTriangle()
                    triangle.GetPointIds().SetId(0, pointscount + face[0])
                    triangle.GetPointIds().SetId(1, pointscount + face[vIdx])
                    triangle.GetPointIds().SetId(2, pointscount + face[vIdx + 1])
                    triangles.InsertNextCell(triangle)
                
                #
                for vIdx in range(nindices - 1):
                    
                    i = pointscount + face[vIdx]
                    j = pointscount + face[vIdx + 1]
                    edge = (min(i, j), max(i, j))
                    
                    if edge not in visitededges:
                        line = vtk.vtkLine()
                        line.GetPointIds().SetId(0, edge[0])
                        line.GetPointIds().SetId(1, edge[1])
                        edges.InsertNextCell(line)
                        visitededges.add(edge)
            
            pointscount += nvertices
        
        # Define the polydata object that contains the information of the faces of
        # the geometric object to be rendered as an actor
        facesdata = vtk.vtkPolyData()
        facesdata.SetPoints(points)
        facesdata.SetPolys(triangles)
        
        # Define polydatanormals object to calculate the normals of the points
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(facesdata)
        normals.ComputePointNormalsOn()
        normals.ComputeCellNormalsOn()
        normals.Update()
        
        # Create a mapper and actor
        facesmapper = vtk.vtkPolyDataMapper()
        facesmapper.SetInputConnection(normals.GetOutputPort())
        
        self.assemblyactor = vtk.vtkActor()
        self.assemblyactor.SetMapper(facesmapper)
        
        edgesdata = vtk.vtkPolyData()
        edgesdata.SetPoints(points)
        edgesdata.SetLines(edges)
        
        edgesmapper = vtk.vtkPolyDataMapper()
        edgesmapper.SetInputData(edgesdata)
        
        self.assemblyedgesactor = vtk.vtkActor()
        self.assemblyedgesactor.SetMapper(edgesmapper)
        
        self.assemblyactor.GetProperty().BackfaceCullingOn()
        self.assemblyactor.GetProperty().SetColor(1, 1, 1)
        self.assemblyedgesactor.GetProperty().SetColor(0, 0, 0)
        self.assemblyedgesactor.GetProperty().SetLineWidth(2)
        
        # Add the assembly actors to the scene
        self.renderer.AddActor(self.assemblyactor)
        self.renderer.AddActor(self.assemblyedgesactor)
    
    
    # Uses the current geometric domain in the workspace to load the actors 
    # that represent it visually
    def loaddomainactors(self, domain:DCEL):
        assert domain is not None, "Missing domain in the workspace"
        assert isinstance(domain, DCEL), "Workspace domain not a DCEL"
        
        # Remove everything from the scene (with a new geometric domain, any 
        # assembly in the scene becomes irrelevant)
        self.clearall()
        
        # Load the new domain actors
        self.domaindcelactor, self.domainedgesdactor = loadvtkactorswithdcel(domain)
        self.domaindcelactor.GetProperty().BackfaceCullingOn()
        self.domaindcelactor.GetProperty().SetColor(1, 1, 0)
        self.domainedgesdactor.GetProperty().SetColor(0, 0, 0)
        self.domainedgesdactor.GetProperty().SetLineWidth(1)
        
        # Add the domain actors to the scene
        self.renderer.AddActor(self.domaindcelactor)
        self.renderer.AddActor(self.domainedgesdactor)
        